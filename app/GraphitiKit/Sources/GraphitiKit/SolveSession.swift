import SwiftUI

@MainActor
public final class SolveSession: ObservableObject {
    public enum Phase: Sendable { case empty, drawing, searching, locked, noFit }

    public struct Step {
        public let time: Double
        public let tokens: [EquationToken]
        public let params: CurveParams
        public var samples: [Double] = []
        public var latex: String = ""
        public var xs: [Double] = []
        public var tried: Int = 0
        public var terms: Int = 0
        public var others: [[CGPoint]] = []
    }

    @Published public var liveStroke: [CGPoint] = []
    @Published public var showStroke = false
    public var liveCurvePoints: [CGPoint] {
        guard !liveStroke.isEmpty else { return [] }
        let st = step; let n = min(st.xs.count, st.samples.count)
        return (0..<n).map { CGPoint(x: st.xs[$0], y: -st.samples[$0]) }
    }

    @Published public private(set) var phase: Phase
    @Published public private(set) var index: Int
    @Published public private(set) var elapsed: Double
    @Published public private(set) var startedAt: Date?
    @Published public private(set) var dying: Set<String> = []
    @Published public private(set) var flash: Double
    @Published public private(set) var indicate = false
    @Published public private(set) var strokes: Int
    @Published public private(set) var convergedAtStep: Int? = nil

    public var steps: [Step] { liveSteps.isEmpty ? SolveSession.demoSteps : liveSteps }
    public let target = CurveParams(m: 0.5, a: 1, w: 2, q: 0)
    static let demoLatex = "y = \\frac{1}{2}x + \\sin(2x)"
    public var latex: String { liveSteps.isEmpty ? Self.demoLatex : step.latex }
    public var plainText: String {
        guard !liveSteps.isEmpty else { return EquationFormatter.plain(target) }
        var s = step.latex.replacingOccurrences(of: "\\cdot", with: "·").replacingOccurrences(of: "\\pi", with: "π")
        for (pattern, template) in [("\\\\frac\\{([^{}]*)\\}\\{([^{}]*)\\}", "($1)/($2)"), ("\\\\([a-zA-Z]+)", "$1"), ("[{}]", "")] {
            s = s.replacingOccurrences(of: pattern, with: template, options: .regularExpression)
        }
        return s
    }
    var colourKey: String { latex }
    public var colour: Color { EquationColor.color(for: colourKey) }

    public struct Term: Identifiable, Sendable {
        public let id: String
        public let tokenIDs: Set<String>
        public let params: CurveParams
        public let hex: UInt32
        public let name: String
        public let caption: String
        public var colour: Color { Color(hex: hex) }
        public var platformColour: PlatformColor { .graphiti(hex: hex) }
    }

    public let terms: [Term] = [
        Term(id: "tilt", tokenIDs: ["c1", "x1"], params: .init(m: 0.5, a: 0, w: 1, q: 0), hex: 0xF0AC5F,
             name: "The tilt", caption: "up ½ for every 1 across"),
        Term(id: "wave", tokenIDs: ["sin", "c2", "x2", "rp"], params: .init(m: 0, a: 1, w: 2, q: 0), hex: 0x5CD0B3,
             name: "The wave", caption: "height 1, repeats every π"),
    ]

    @Published public var explained: String?
    public var explainedTerm: Term? { terms.first { $0.id == explained } }
    public func term(containing token: String) -> Term? { terms.first { $0.tokenIDs.contains(token) } }

    public func explain(_ id: String?) {
        withAnimation(.easeInOut(duration: 0.25)) { explained = (explained == id) ? nil : id }
    }
    public var platformColour: PlatformColor { EquationColor.platformColor(for: colourKey) }

    public init(step: Int = 0, locked: Bool = false, flash: Double = 0, indicate: Bool = false,
                dying: Set<String> = [], phase: Phase? = nil, elapsed: Double? = nil, explained: String? = nil) {
        self.explained = explained
        let i = min(step, SolveSession.demoSteps.count - 1)
        self.index = i
        self.phase = phase ?? (locked ? .locked : .searching)
        self.flash = flash
        self.indicate = indicate
        self.dying = dying
        self.elapsed = elapsed ?? SolveSession.demoSteps[i].time
        self.strokes = (phase == .empty) ? 0 : 1
    }

    public var locked: Bool { phase == .locked }
    public var isEmpty: Bool { phase == .empty }
    public var step: Step { steps[min(index, steps.count - 1)] }

    public func status(at now: Date = Date()) -> String {
        switch phase {
        case .searching: String(format: "%.1f s", startedAt.map { now.timeIntervalSince($0) } ?? elapsed)
            + (!liveStroke.isEmpty && liveSteps.isEmpty ? " · preparing GPU" : Self.triedText(step.tried))
        case .locked: String(format: "Found in %.1f s", elapsed) + Self.triedText(steps.last?.tried ?? 0)
            + (convergedAtStep.map { " · converged at step \($0)" } ?? "")
        default: ""
        }
    }

    static func triedText(_ n: Int) -> String {
        guard n > 0 else { return "" }
        let v = Double(n)
        let s = v >= 1e9 ? String(format: "%.1fB", v / 1e9) : v >= 1e6 ? String(format: "%.1fM", v / 1e6) : String(format: "%.0fk", v / 1e3)
        return " · \(s) tried"
    }

    public var changed: Set<String> {
        guard index > 0 else { return Set(step.tokens.map(\.id)) }
        let prev = Dictionary(uniqueKeysWithValues: steps[index - 1].tokens.map { ($0.id, $0.latex) })
        return Set(step.tokens.filter { prev[$0.id] != $0.latex }.map(\.id))
    }

    public func undo() {
        guard strokes > 0 else { return }
        withAnimation(.easeOut(duration: 0.25)) { strokes -= 1; if strokes == 0 { phase = .empty } }
    }

    public func new() {
        withAnimation(.easeOut(duration: 0.25)) { strokes = 0; phase = .empty; index = 0; elapsed = 0; startedAt = nil; liveStroke = []; liveSteps = []; engineMap = [] }
    }

    @Published public private(set) var engineMap: [Double] = []

    public private(set) var finalPick: (latex: String, score: Double)? = nil

    public func consume(_ source: BurstSource, map: (() -> [Double])? = nil) async {
        finalPick = nil; convergedAtStep = nil
        var converge: ConvergeStop? = ConvergeStop.enabled ? ConvergeStop() : nil
        withAnimation(.easeOut(duration: 0.15)) {
            index = 0; phase = .searching; flash = 0; dying = []; elapsed = 0; strokes = 1; explained = nil
        }
        liveSteps = []; engineMap = []; showStroke = false
        let start = Date()
        startedAt = start
        var lastShown = Date.distantPast
        let spring = Animation.spring(response: 0.3, dampingFraction: 0.88)
        for await b in source.bursts() {
            if b.isFinal, let c = b.champion { finalPick = (BurstDisplay.latex(c), c.score) }
            guard let champ = b.champion else { if b.isFinal { break } else { continue } }
            if engineMap.isEmpty, let map { engineMap = map() }
            if !b.isFinal, converge != nil, converge!.observe(step: b.step, champion: champ) {
                convergedAtStep = b.step; converge = nil
                EnginePrepare.stopRunning()
            }
            let others = b.pairs.filter { !$0.isChampion }.prefix(8).map { p in
                (0..<min(p.x.count, p.y.count)).map { CGPoint(x: p.x[$0], y: -p.y[$0]) } }
            let step = Step(time: Date().timeIntervalSince(start), tokens: BurstDisplay.tokens(champ),
                            params: .zero, samples: champ.y, latex: BurstDisplay.latex(champ), xs: champ.x,
                            tried: b.isFinal ? (liveSteps.last?.tried ?? 0) : max(b.step, 0) * b.candidates,
                            terms: BurstDisplay.terms(champ), others: Array(others))
            if !b.isFinal, let last = liveSteps.last, last.latex == step.latex {
                liveSteps[liveSteps.count - 1] = step
                continue
            }
            let gap = Date().timeIntervalSince(lastShown)
            if gap < 0.25, lastShown > .distantPast { try? await Task.sleep(for: .seconds(0.25 - gap)) }
            lastShown = Date()
            if let prev = liveSteps.last {
                let changedTree = prev.tokens.map(\.latex) != step.tokens.map(\.latex)
                if changedTree { SolveSounds.shared.play(prev.latex.count > step.latex.count + 8 ? .prune : .term(liveSteps.count)) }
                else { SolveSounds.shared.play(.tune(liveSteps.count)) }
            }
            liveSteps.append(step)
            withAnimation(spring) { index = liveSteps.count - 1 }
            if b.isFinal { break }
        }
        elapsed = Date().timeIntervalSince(start)
        startedAt = nil
        if liveSteps.isEmpty {
            withAnimation(.easeOut(duration: 0.25)) { strokes = 0; phase = .empty; index = 0; elapsed = 0; liveStroke = []; engineMap = [] }
            return
        }
        withAnimation(.easeInOut(duration: 0.6)) { phase = .locked }
        SolveSounds.shared.play(.lock)
        withAnimation(.easeInOut(duration: 0.45)) { flash = 1 }
        try? await Task.sleep(for: .milliseconds(450))
        withAnimation(.spring(response: 0.25, dampingFraction: 0.7)) { indicate = true }
        try? await Task.sleep(for: .milliseconds(260))
        withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) { indicate = false; flash = 0 }
    }

    public func solveEngine(_ points: [CGPoint], budget: Double) async {
        liveStroke = points
        _ = await EngineSolve.run(points: points, budget: budget, session: self)
    }

    public func solveLive(_ source: BurstSource) async {
        await consume(source)
    }

    @Published public private(set) var liveSteps: [Step] = []

    public func play() async {
        withAnimation(.easeOut(duration: 0.15)) {
            index = 0; phase = .searching; flash = 0; dying = []; elapsed = 0; strokes = 1; explained = nil
        }
        try? await Task.sleep(for: .milliseconds(32))
        var start = Date()
        startedAt = start
        var lastShown = Date()
        let spring = Animation.spring(response: 0.3, dampingFraction: 0.88)
        for next in 1..<steps.count {
            let wait = steps[next].time - Date().timeIntervalSince(start)
            if wait > 0 { try? await Task.sleep(for: .seconds(wait)) }
            let gap = Date().timeIntervalSince(lastShown)
            if gap < 0.25 {
                try? await Task.sleep(for: .seconds(0.25 - gap))
                start = start.addingTimeInterval(0.25 - gap)
            }
            lastShown = Date()
            let leaving = Set(steps[next - 1].tokens.map(\.id)).subtracting(steps[next].tokens.map(\.id))
            if !leaving.isEmpty {
                SolveSounds.shared.play(.prune)
                withAnimation(.easeOut(duration: 0.12)) { dying = leaving }
                try? await Task.sleep(for: .milliseconds(140))
            }
            let structureChanged = steps[next].tokens.map(\.id) != steps[next - 1].tokens.map(\.id)
            if leaving.isEmpty { SolveSounds.shared.play(structureChanged ? .term(next) : .tune(next)) }
            withAnimation(spring) { index = next; dying = [] }
        }
        elapsed = startedAt.map { Date().timeIntervalSince($0) } ?? steps[steps.count - 1].time
        startedAt = nil
        withAnimation(.easeInOut(duration: 0.6)) { phase = .locked }
        SolveSounds.shared.play(.lock)
        withAnimation(.easeInOut(duration: 0.45)) { flash = 1 }
        try? await Task.sleep(for: .milliseconds(450))
        withAnimation(.spring(response: 0.25, dampingFraction: 0.7)) { indicate = true }
        try? await Task.sleep(for: .milliseconds(260))
        withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) { indicate = false; flash = 0 }
    }

    public nonisolated static func warmUp() async {
        let colours: [PlatformColor] = [.white, .graphiti(hex: 0xFFFFFF, alpha: 0.28),
                                        EquationColor.platformColor(for: "y = \\frac{1}{2}x + \\sin(2x)")]
        await Task.detached(priority: .utility) {
            for step in demoSteps {
                for t in step.tokens {
                    for size in [Typeface.equationSize, 30, 26, 22] {
                        let s = t.kind == .exponent ? size * 0.68 : size
                        for c in colours { _ = MathRenderCache.shared.render(t.latex, size: s, color: c) }
                    }
                }
            }
        }.value
    }

    public struct Stats: Sendable {
        public let seconds: Double
        public let generations: Int
        public let equationsTried: Int
        public let terms: Int
        public let fit: Double
        public let exactConstants: Int
        public let fitHistory: [(Int, Double)]
        public let properties: [(String, String)]
    }

    static let population = 60

    public var stats: Stats {
        if !liveStroke.isEmpty {
            let fits = liveSteps.enumerated().map { ($0.offset + 1, Self.shapeFit(stroke: liveStroke, curve: Self.points($0.element))) }
            return Stats(seconds: elapsed, generations: liveSteps.count, equationsTried: liveSteps.last?.tried ?? 0,
                         terms: liveSteps.last?.terms ?? 0, fit: fits.last?.1 ?? 0, exactConstants: 0,
                         fitHistory: fits, properties: [])
        }
        return Stats(seconds: elapsed, generations: 31, equationsTried: 31 * SolveSession.population,
              terms: 2, fit: 0.996, exactConstants: 1,
              fitHistory: [(1, 0.42), (4, 0.58), (9, 0.71), (14, 0.83), (19, 0.90), (23, 0.95), (26, 0.97), (31, 0.992), (31, 0.996)],
              properties: [("slope", "\\frac{1}{2}"), ("amplitude", "1"), ("period", "\\pi")])
    }

    static func points(_ s: Step) -> [CGPoint] { (0..<min(s.xs.count, s.samples.count)).map { CGPoint(x: s.xs[$0], y: -s.samples[$0]) } }
    static func shapeFit(stroke: [CGPoint], curve: [CGPoint]) -> Double {
        guard stroke.count > 1, curve.count > 1 else { return 0 }
        let xs = stroke.map(\.x), ys = stroke.map(\.y)
        let diag = hypot((xs.max()! - xs.min()!), (ys.max()! - ys.min()!)); guard diag > 0 else { return 0 }
        var ss = 0.0
        for p in stroke { var best = Double.infinity; for q in curve { best = min(best, Double((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y))) }; ss += best }
        return max(0, 1 - (ss / Double(stroke.count)).squareRoot() / Double(diag))
    }

    static let demoSteps: [Step] = {
        typealias T = EquationToken
        let y = T("y", "y"), eq = T("eq", "=", .op), x1 = T("x1", "x"), p1 = T("p1", "+", .op)
        let sin = T("sin", "\\sin("), x2 = T("x2", "x"), rp = T("rp", ")")
        let p2 = T("p2", "+", .op), c3 = T("c3", "0.1", .constant), x3 = T("x3", "x"), e3 = T("e3", "2", .exponent)
        func c1(_ v: String) -> T { T("c1", v, .constant) }
        func c2(_ v: String) -> T { T("c2", v, .constant) }
        return [
            Step(time: 0.0, tokens: [y, eq, x1], params: .init(m: 1, a: 0, w: 1, q: 0)),
            Step(time: 0.2, tokens: [y, eq, c1("0.4"), x1], params: .init(m: 0.4, a: 0, w: 1, q: 0)),
            Step(time: 0.5, tokens: [y, eq, c1("0.4"), x1, p1, sin, x2, rp], params: .init(m: 0.4, a: 1, w: 1, q: 0)),
            Step(time: 0.8, tokens: [y, eq, c1("0.4"), x1, p1, sin, c2("1.7"), x2, rp], params: .init(m: 0.4, a: 1, w: 1.7, q: 0)),
            Step(time: 1.1, tokens: [y, eq, c1("0.5"), x1, p1, sin, c2("1.9"), x2, rp], params: .init(m: 0.5, a: 1, w: 1.9, q: 0)),
            Step(time: 1.4, tokens: [y, eq, c1("0.5"), x1, p1, sin, c2("2"), x2, rp, p2, c3, x3], params: .init(m: 0.6, a: 1, w: 2, q: 0)),
            Step(time: 1.7, tokens: [y, eq, c1("0.5"), x1, p1, sin, c2("2"), x2, rp, p2, c3, x3, e3], params: .init(m: 0.5, a: 1, w: 2, q: 0.1)),
            Step(time: 2.1, tokens: [y, eq, c1("0.5"), x1, p1, sin, c2("2"), x2, rp], params: .init(m: 0.5, a: 1, w: 2, q: 0)),
            Step(time: 2.4, tokens: [y, eq, c1("\\frac{1}{2}"), x1, p1, sin, c2("2"), x2, rp], params: .init(m: 0.5, a: 1, w: 2, q: 0)),
        ]
    }()
}
