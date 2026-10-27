import SwiftUI
import CoreGraphics
#if canImport(UIKit)
import UIKit
#endif


// MARK: - Palette
private enum P {
    static let bgPage   = Color(red: 0.020, green: 0.027, blue: 0.043)
    static let stage    = Color(red: 0.027, green: 0.035, blue: 0.059)
    static let accent   = Color(red: 0.549, green: 0.420, blue: 1.0)
    static let accentTo = Color(red: 0.388, green: 0.400, blue: 0.945)
    static let cyan     = Color(hue: 0.520, saturation: 0.62, brightness: 0.90)
    static let curve    = Color(red: 0.718, green: 0.580, blue: 0.965)
    static let ink      = Color.white
    static let dim      = Color(red: 0.580, green: 0.639, blue: 0.722)
    static let hairline = Color.white.opacity(0.12)

    static let manimPalette: [Color] = [
        Color(hue: 0.520, saturation: 0.62, brightness: 0.90),
        Color(hue: 0.720, saturation: 0.55, brightness: 0.98),
        Color(hue: 0.360, saturation: 0.58, brightness: 0.82),
        Color(hue: 0.120, saturation: 0.70, brightness: 0.95),
        Color(hue: 0.930, saturation: 0.55, brightness: 0.95),
        Color(hue: 0.590, saturation: 0.60, brightness: 0.95),
        Color(hue: 0.030, saturation: 0.68, brightness: 0.95),
    ]
    static func manimAccentIndex(for pts: [CGPoint]) -> Int {
        guard !pts.isEmpty else { return 1 }
        var h: UInt64 = 1469598103934665603
        for p in pts {
            let qx = Int((p.x / 8).rounded()), qy = Int((p.y / 8).rounded())
            h = (h ^ UInt64(bitPattern: Int64(qx &* 73856093 ^ qy &* 19349663))) &* 1099511628211
        }
        return Int(h % UInt64(manimPalette.count))
    }
    static func manimAccent(for pts: [CGPoint]) -> Color { paletteColor(manimAccentIndex(for: pts)) }
    static func paletteColor(_ i: Int) -> Color {
        let n = manimPalette.count
        return manimPalette[((i % n) + n) % n]
    }
}

// MARK: - Geometry
private struct BBox {
    var minX = 0.0, minY = 0.0, maxX = 1.0, maxY = 1.0
    var cx: Double { (minX + maxX) * 0.5 }
    var cy: Double { (minY + maxY) * 0.5 }
    var span: Double { max(max(maxX - minX, maxY - minY), 1e-6) }
    init() {}
    init(_ pts: [CGPoint]) {
        guard !pts.isEmpty else { return }
        let xs = pts.map { Double($0.x) }, ys = pts.map { Double($0.y) }
        minX = xs.min()!; maxX = xs.max()!; minY = ys.min()!; maxY = ys.max()!
    }
}

private func mapCurve(_ curve: [CGPoint], to box: BBox) -> [CGPoint] {
    guard !curve.isEmpty else { return [] }
    let cb = BBox(curve)
    let s = box.span / cb.span
    return curve.map { CGPoint(x: box.cx + (Double($0.x) - cb.cx) * s,
                               y: box.cy + (Double($0.y) - cb.cy) * s) }
}

private func pathThrough(_ pts: [CGPoint]) -> Path {
    var p = Path()
    guard let f = pts.first else { return p }
    p.move(to: f)
    for pt in pts.dropFirst() { p.addLine(to: pt) }
    return p
}

private func fitPoints(_ curve: [CGPoint], in size: CGSize, pad: Double = 0.82) -> [CGPoint] {
    guard !curve.isEmpty else { return [] }
    let b = BBox(curve)
    let s = Double(min(size.width, size.height)) * pad / b.span
    let w = Double(size.width), h = Double(size.height)
    return curve.map { CGPoint(x: w * 0.5 + (Double($0.x) - b.cx) * s,
                               y: h * 0.5 + (Double($0.y) - b.cy) * s) }
}

// MARK: - Model
enum NoisePredictorCache {
    static let shared: ContinuousNoisePredictor? = {
        #if SWIFT_PACKAGE
        return ContinuousNoisePredictor(bundle: .module)
        #else
        return ContinuousNoisePredictor(bundle: .main)
        #endif
    }()
}

@MainActor
final class SolveModel: ObservableObject {
    @Published var inputStroke: [CGPoint] = []
    @Published var isSolving = false
    @Published var champion: ParametricCandidate?
    @Published var candidates: [ParametricCandidate] = []
    @Published var generation = 0
    @Published var elapsedMs = 0
    @Published var curveTrim: CGFloat = 0
    @Published var locked = false
    @Published var noFit = false
    @Published var placed: [PlacedResult] = []
    @Published var netDiag: String = ""
    @Published private(set) var accent: Color = P.accent
    private(set) var accentIndex: Int = 1

    weak var gallery: GalleryStore?
    private var solveTask: Task<Void, Never>?
    private let noiseEnabled = ProcessInfo.processInfo.environment["SKETCH_CONTINUOUS_NOISE"] != "0"

    var canSolve: Bool { inputStroke.count >= 6 && !isSolving }

    fileprivate func championCurve(in box: BBox) -> [CGPoint] {
        guard let c = champion else { return [] }
        return mapCurve(c.curvePoints, to: box)
    }

    func clear() {
        solveTask?.cancel(); solveTask = nil
        inputStroke = []; champion = nil; candidates = []
        generation = 0; elapsedMs = 0; curveTrim = 0; isSolving = false; locked = false; noFit = false
        placed = []
        accent = P.accent
    }

    func placeCurrentAndReset() {
        guard let c = champion, locked else { return }
        placed.append(PlacedResult(stroke: inputStroke, curvePoints: c.curvePoints, accentIndex: accentIndex))
        solveTask?.cancel(); solveTask = nil
        inputStroke = []; champion = nil; candidates = []
        generation = 0; elapsedMs = 0; curveTrim = 0; isSolving = false; locked = false; noFit = false
        accent = P.accent
    }

    func beginStroke(at p: CGPoint) {
        if isSolving { solveTask?.cancel(); solveTask = nil; isSolving = false }
        champion = nil; candidates = []; curveTrim = 0; locked = false; noFit = false
        inputStroke = [p]
    }
    func extendStroke(to p: CGPoint) { inputStroke.append(p) }

    func solve() {
        let pts = inputStroke
        guard pts.count >= 6, pts.allSatisfy({ $0.x.isFinite && $0.y.isFinite }) else { return }
        accentIndex = P.manimAccentIndex(for: pts)
        accent = P.paletteColor(accentIndex)
        solveTask?.cancel()
        champion = nil; candidates = []; curveTrim = 0; locked = false; noFit = false
        isSolving = true; generation = 0; elapsedMs = 0

        let useNoise = noiseEnabled

        solveTask = Task { [weak self] in
            let field: NativeNoiseField? = useNoise
                ? await Task.detached(priority: .userInitiated) {
                    NoisePredictorCache.shared?.makeNativeNoiseField(points: pts)
                  }.value
                : nil
            let stream = SketchNativeEngine.shared.solveStrokeStream(
                points: pts, population: 500, generations: 60, maxLength: 28, noiseField: field)
            for await burst in stream {
                if Task.isCancelled { break }
                await MainActor.run {
                    guard let self else { return }
                    self.generation = burst.generation
                    self.elapsedMs = burst.elapsedMs
                    self.candidates = burst.candidates.sorted { $0.shapeErr < $1.shapeErr }
                    let isNew = self.champion?.id != burst.champion.id
                    let firstReveal = (self.champion == nil)
                    self.champion = burst.champion
                    if firstReveal {
                        self.curveTrim = 0
                        withAnimation(.timingCurve(0.65, 0, 0.35, 1, duration: 0.7)) { self.curveTrim = 1 }
                    } else if isNew {
                        withAnimation(.easeInOut(duration: 0.28)) { self.curveTrim = 1 }
                    }
                }
            }
            await MainActor.run {
                guard let self else { return }
                self.isSolving = false
                if let c = self.champion, let acc = self.candidates.first {
                    let moved = (c.exprX != acc.exprX) || (c.exprY != acc.exprY)
                    let s2 = (field?.s2 ?? []).sorted()
                    let s2med = s2.isEmpty ? Float(-1) : s2[s2.count / 2]
                    self.netDiag = String(format: "net:%@ cand:%d acc:%dn/%.4f→champ:%dn/%.4f moved:%@ s2med:%.2e",
                                          field == nil ? "NONE" : "ok", self.candidates.count,
                                          acc.nodes, acc.shapeErr, c.nodes, c.shapeErr,
                                          moved ? "YES" : "no", Double(s2med))
                    print("NETDIAG \(self.netDiag)")
                }
                let viable = (self.champion?.r2 ?? 0) >= 0.55
                if viable {
                    withAnimation(.easeOut(duration: 0.5)) { self.curveTrim = 1; self.locked = true }
                    if let c = self.champion {
                        self.gallery?.add(champion: c, stroke: self.inputStroke, accentIndex: self.accentIndex)
                    }
                } else {
                    self.champion = nil; self.candidates = []; self.curveTrim = 0; self.locked = false
                    withAnimation(.easeInOut(duration: 0.35)) { self.noFit = true }
                }
            }
        }
    }

}

struct PlacedResult: Identifiable {
    let id = UUID()
    let stroke: [CGPoint]
    let curvePoints: [CGPoint]
    let accentIndex: Int
    var accent: Color { P.paletteColor(accentIndex) }
}

struct SavedSketch: Codable, Identifiable, Equatable {
    let id: UUID
    let date: Date
    let exprX: String
    let exprY: String
    let r2: Double
    let nodes: Int
    let accentIndex: Int
    let curvePoints: [CGPoint]
    let stroke: [CGPoint]

    var accent: Color { P.paletteColor(accentIndex) }
}

@MainActor
final class GalleryStore: ObservableObject {
    @Published private(set) var sketches: [SavedSketch] = []

    private let fileURL: URL = {
        let dir = (try? FileManager.default.url(for: .applicationSupportDirectory, in: .userDomainMask, appropriateFor: nil, create: true))
            ?? FileManager.default.urls(for: .documentDirectory, in: .userDomainMask)[0]
        return dir.appendingPathComponent("graphart_gallery.json")
    }()

    init() { load() }

    func add(champion c: ParametricCandidate, stroke: [CGPoint], accentIndex: Int) {
        let s = SavedSketch(id: UUID(), date: Date(),
                            exprX: c.displayExprX.isEmpty ? c.exprX : c.displayExprX,
                            exprY: c.displayExprY.isEmpty ? c.exprY : c.displayExprY,
                            r2: c.r2, nodes: c.nodes, accentIndex: accentIndex,
                            curvePoints: c.curvePoints, stroke: stroke)
        sketches.insert(s, at: 0)
        persist()
    }

    func delete(_ s: SavedSketch) {
        sketches.removeAll { $0.id == s.id }
        persist()
    }

    private func load() {
        guard let data = try? Data(contentsOf: fileURL),
              let decoded = try? JSONDecoder().decode([SavedSketch].self, from: data) else { return }
        sketches = decoded
    }

    private func persist() {
        guard let data = try? JSONEncoder().encode(sketches) else { return }
        try? data.write(to: fileURL, options: .atomic)
    }
}

// MARK: - Root
public struct GraphArtNativeRootView: View {
    @StateObject private var model = SolveModel()
    @StateObject private var gallery = GalleryStore()
    @State private var view = 0

    public init() {}

    public var body: some View {
        ZStack {
            P.bgPage.ignoresSafeArea()
            VStack(spacing: 0) {
                HeaderBar(model: model)
                TabView(selection: $view) {
                    GraphView(model: model).tag(0)
                    SymbolicView(model: model).tag(1)
                    GalleryView(store: gallery).tag(2)
                }
                #if os(iOS)
                .tabViewStyle(.page(indexDisplayMode: .never))
                #endif
                Dock(view: $view, count: 3)
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .top)

            if ProcessInfo.processInfo.environment["SKETCH_DEMO_SHOWCARD"] == "1", let c = model.champion {
                ExportCard(e: SketchExport(curvePoints: c.curvePoints, exprX: c.exprX, exprY: c.exprY,
                                           r2: c.r2, nodes: c.nodes, accent: model.accent))
                    .clipShape(RoundedRectangle(cornerRadius: 22))
                    .shadow(color: .black.opacity(0.6), radius: 30)
            }
        }
        .preferredColorScheme(.dark)
        .onAppear {
            model.gallery = gallery
            if let t = ProcessInfo.processInfo.environment["SKETCH_DEMO_TAB"], let n = Int(t) { view = n }
        }
        .task {
            guard ProcessInfo.processInfo.environment["SKETCH_DEMO_SOLVE"] == "1" else { return }
            let noisy = ProcessInfo.processInfo.environment["SKETCH_DEMO_NOISY"] == "1"
            let scribble = ProcessInfo.processInfo.environment["SKETCH_DEMO_SCRIBBLE"] == "1"
            var seed: UInt64 = 0x9E3779B9
            func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
            func gauss() -> Double { let u1 = max(1e-9, rnd()), u2 = rnd(); return (-2*log(u1)).squareRoot() * cos(2*Double.pi*u2) }
            let shape = ProcessInfo.processInfo.environment["SKETCH_DEMO_SHAPE"] ?? (scribble ? "scribble" : "circle")
            var pts: [CGPoint] = []
            switch shape {
            case "scribble":
                var x = 200.0, y = 420.0
                for _ in 0..<120 { x += 26 * (rnd() - 0.5); y += 26 * (rnd() - 0.5); pts.append(CGPoint(x: x, y: y)) }
            case "liss":
                for i in 0..<160 { let t = Double(i)/160*2*Double.pi
                    pts.append(CGPoint(x: 200 + 130*sin(3*t), y: 420 + 130*sin(2*t))) }
            case "spiral":
                for i in 0..<160 { let t = Double(i)/160*4*Double.pi; let r = 20 + 12*t
                    pts.append(CGPoint(x: 200 + r*cos(t), y: 420 + r*sin(t))) }
            default:
                for i in 0..<80 {
                    let a = Double(i) / 80.0 * 2.0 * Double.pi
                    pts.append(CGPoint(x: 200 + 120 * cos(a) + (noisy ? 2.4 * gauss() : 0),
                                       y: 420 + 120 * sin(a) + (noisy ? 2.4 * gauss() : 0)))
                }
            }
            if ProcessInfo.processInfo.environment["SKETCH_DEMO_MULTI"] == "1" {
                var a: [CGPoint] = [], b: [CGPoint] = []
                for i in 0..<80 { let t = Double(i)/80*2*Double.pi
                    a.append(CGPoint(x: 120 + 70*cos(t), y: 300 + 70*sin(t)))
                    b.append(CGPoint(x: 300 + 80*sin(3*t), y: 640 + 80*sin(2*t))) }
                model.placed = [PlacedResult(stroke: a, curvePoints: a, accentIndex: P.manimAccentIndex(for: a)),
                                PlacedResult(stroke: b, curvePoints: b, accentIndex: P.manimAccentIndex(for: b))]
            }
            model.inputStroke = pts
            if ProcessInfo.processInfo.environment["SKETCH_DEMO_NOFIT"] == "1" {
                model.isSolving = false
                withAnimation(.easeInOut(duration: 0.35)) { model.noFit = true }
            } else {
                model.solve()
            }
        }
    }
}

public typealias GraphArtWebViewContainer = GraphArtNativeRootView

private struct BrandGlyph: View {
    var stroke: Color = .white
    var smile: Color? = nil
    var body: some View {
        Canvas { ctx, size in
            let w = size.width, h = size.height
            let lw = max(1.4, w * 0.085)
            func pt(_ x: CGFloat, _ y: CGFloat) -> CGPoint { CGPoint(x: x * w, y: y * h) }
            func head(_ tip: CGPoint, _ dx: CGFloat, _ dy: CGFloat) -> Path {
                let s = w * 0.10
                var p = Path(); p.move(to: CGPoint(x: tip.x - s * dx - s * dy, y: tip.y - s * dy - s * dx))
                p.addLine(to: tip)
                p.addLine(to: CGPoint(x: tip.x - s * dx + s * dy, y: tip.y - s * dy + s * dx))
                return p
            }
            var axes = Path()
            axes.move(to: pt(0.5, 0.10)); axes.addLine(to: pt(0.5, 0.90))
            axes.move(to: pt(0.10, 0.5)); axes.addLine(to: pt(0.90, 0.5))
            ctx.stroke(axes, with: .color(stroke.opacity(0.42)), style: StrokeStyle(lineWidth: lw * 0.8, lineCap: .round))
            let arrowStyle = StrokeStyle(lineWidth: lw * 0.8, lineCap: .round, lineJoin: .round)
            ctx.stroke(head(pt(0.5, 0.10), 0, 1), with: .color(stroke.opacity(0.42)), style: arrowStyle)
            ctx.stroke(head(pt(0.9, 0.5), 1, 0), with: .color(stroke.opacity(0.42)), style: arrowStyle)
            let r = w * 0.085
            for cx in [0.34, 0.66] as [CGFloat] {
                let c = pt(cx, 0.30)
                ctx.stroke(Path(ellipseIn: CGRect(x: c.x - r, y: c.y - r, width: 2 * r, height: 2 * r)),
                           with: .color(stroke), style: StrokeStyle(lineWidth: lw))
            }
            var smilePath = Path(); smilePath.move(to: pt(0.26, 0.64))
            smilePath.addQuadCurve(to: pt(0.74, 0.64), control: pt(0.5, 0.86))
            ctx.stroke(smilePath, with: .color(smile ?? stroke), style: StrokeStyle(lineWidth: lw, lineCap: .round))
        }
        .aspectRatio(1, contentMode: .fit)
    }
}

// MARK: - Header
private struct HeaderBar: View {
    @ObservedObject var model: SolveModel

    var body: some View {
        HStack(spacing: 12) {
            BrandGlyph(stroke: .white, smile: model.accent)
                .frame(width: 30, height: 30)
                .shadow(color: model.accent.opacity(0.35), radius: 6)
            Text("GraphArt").font(.system(size: 20, weight: .bold)).foregroundStyle(P.ink)
                .fixedSize()
            Spacer(minLength: 8)
            Button(action: model.clear) {
                Label("Clear", systemImage: "arrow.counterclockwise")
                    .font(.system(size: 12, weight: .semibold)).foregroundStyle(P.dim)
                    .lineLimit(1).fixedSize()
                    .padding(.horizontal, 13).padding(.vertical, 7)
                    .background(Capsule().fill(.white.opacity(0.07)).overlay(Capsule().stroke(.white.opacity(0.14))))
            }.buttonStyle(.plain)
            Button(action: model.solve) {
                Label(model.isSolving ? "Solving…" : "Solve", systemImage: "sparkles")
                    .font(.system(size: 12.5, weight: .bold)).foregroundStyle(.white)
                    .lineLimit(1).fixedSize()
                    .padding(.horizontal, 16).padding(.vertical, 7)
                    .background(Capsule().fill(LinearGradient(colors: [P.accent, P.accentTo], startPoint: .topLeading, endPoint: .bottomTrailing)))
                    .overlay(Capsule().stroke(.white.opacity(0.35)))
                    .shadow(color: P.accent.opacity(0.5), radius: 10, y: 4)
                    .opacity(model.canSolve ? 1 : 0.45)
            }.buttonStyle(.plain).disabled(!model.canSolve)
        }
        .padding(.horizontal, 20).frame(height: 60)
        .background(.ultraThinMaterial.opacity(0.9))
        .overlay(Rectangle().frame(height: 1).foregroundStyle(.white.opacity(0.08)), alignment: .bottom)
    }
}

private struct GridBackground: View {
    var body: some View {
        Canvas { ctx, size in
            let w = size.width, h = size.height
            let cx = w / 2, cy = h / 2
            let minor: CGFloat = 32, major: CGFloat = 128
            ctx.drawLayer { layer in
                var pMinor = Path()
                var x = cx.truncatingRemainder(dividingBy: minor)
                while x <= w { pMinor.move(to: CGPoint(x: x, y: 0)); pMinor.addLine(to: CGPoint(x: x, y: h)); x += minor }
                var y = cy.truncatingRemainder(dividingBy: minor)
                while y <= h { pMinor.move(to: CGPoint(x: 0, y: y)); pMinor.addLine(to: CGPoint(x: w, y: y)); y += minor }
                layer.stroke(pMinor, with: .color(.white.opacity(0.035)), lineWidth: 0.5)
                var pMajor = Path()
                x = cx.truncatingRemainder(dividingBy: major)
                while x <= w { pMajor.move(to: CGPoint(x: x, y: 0)); pMajor.addLine(to: CGPoint(x: x, y: h)); x += major }
                y = cy.truncatingRemainder(dividingBy: major)
                while y <= h { pMajor.move(to: CGPoint(x: 0, y: y)); pMajor.addLine(to: CGPoint(x: w, y: y)); y += major }
                layer.stroke(pMajor, with: .color(.white.opacity(0.07)), lineWidth: 0.75)
                var axes = Path()
                axes.move(to: CGPoint(x: 0, y: cy)); axes.addLine(to: CGPoint(x: w, y: cy))
                axes.move(to: CGPoint(x: cx, y: 0)); axes.addLine(to: CGPoint(x: cx, y: h))
                layer.stroke(axes, with: .color(.white.opacity(0.34)), lineWidth: 1)
                let a: CGFloat = 6
                var heads = Path()
                heads.move(to: CGPoint(x: w - 2, y: cy)); heads.addLine(to: CGPoint(x: w - 2 - a, y: cy - a * 0.7)); heads.move(to: CGPoint(x: w - 2, y: cy)); heads.addLine(to: CGPoint(x: w - 2 - a, y: cy + a * 0.7))
                heads.move(to: CGPoint(x: cx, y: 2)); heads.addLine(to: CGPoint(x: cx - a * 0.7, y: 2 + a)); heads.move(to: CGPoint(x: cx, y: 2)); heads.addLine(to: CGPoint(x: cx + a * 0.7, y: 2 + a))
                layer.stroke(heads, with: .color(.white.opacity(0.34)), lineWidth: 1)
            }
        }
        .mask(
            RadialGradient(colors: [.white, .white.opacity(0.15), .clear],
                           center: .center, startRadius: 40, endRadius: 460)
        )
        .allowsHitTesting(false)
    }
}

private struct ShakeEffect: GeometryEffect {
    var amount: CGFloat = 9
    var shakes: CGFloat = 3
    var animatableData: CGFloat
    func effectValue(size: CGSize) -> ProjectionTransform {
        let decay = 1 - animatableData
        let dx = amount * decay * sin(animatableData * .pi * 2 * shakes)
        return ProjectionTransform(CGAffineTransform(translationX: dx, y: 0))
    }
}

private struct GraphView: View {
    @ObservedObject var model: SolveModel
    @State private var dragging = false
    @State private var searchPulse = false
    @State private var shakeUnit: CGFloat = 0
    @State private var lockBeat: CGFloat = 0

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .bottom) {
                P.stage
                    .contentShape(Rectangle())
                    .gesture(
                        DragGesture(minimumDistance: 0)
                            .onChanged { v in
                                if !dragging { dragging = true; model.beginStroke(at: v.location) }
                                else { model.extendStroke(to: v.location) }
                            }
                            .onEnded { _ in dragging = false }
                    )
                GridBackground()
                ForEach(model.placed) { p in
                    let pAccent = p.accent
                    let pCurve = mapCurve(p.curvePoints, to: BBox(p.stroke))
                    pathThrough(p.stroke)
                        .stroke(P.ink.opacity(0.22), style: StrokeStyle(lineWidth: 2, lineCap: .round, lineJoin: .round))
                        .allowsHitTesting(false)
                    if pCurve.count > 1 {
                        pathThrough(pCurve).stroke(pAccent.opacity(0.32), style: StrokeStyle(lineWidth: 10, lineCap: .round, lineJoin: .round))
                            .blur(radius: 8).allowsHitTesting(false)
                        pathThrough(pCurve).stroke(pAccent, style: StrokeStyle(lineWidth: 2.6, lineCap: .round, lineJoin: .round))
                            .shadow(color: pAccent.opacity(0.7), radius: 8).allowsHitTesting(false)
                    }
                }
                pathThrough(model.inputStroke)
                    .stroke(P.ink, style: StrokeStyle(lineWidth: 2.5, lineCap: .round, lineJoin: .round))
                    .shadow(color: .white.opacity(0.25), radius: 3)
                    .modifier(ShakeEffect(animatableData: shakeUnit))
                    .allowsHitTesting(false)
                let curve = model.championCurve(in: BBox(model.inputStroke))
                if curve.count > 1 {
                    let trimmed = pathThrough(curve).trim(from: 0, to: model.curveTrim)
                    let solvingBoost = (model.isSolving && searchPulse) ? 0.28 : 0.0
                    Group {
                        if lockBeat > 0.001 {
                            trimmed.stroke(model.accent.opacity(0.55 * lockBeat),
                                           style: StrokeStyle(lineWidth: 18, lineCap: .round, lineJoin: .round))
                                .blur(radius: 18).allowsHitTesting(false)
                        }
                        trimmed.stroke(model.accent.opacity((model.locked ? 0.35 : 0.2) + solvingBoost),
                                       style: StrokeStyle(lineWidth: (model.locked ? 12 : 8) + (model.isSolving && searchPulse ? 6 : 0), lineCap: .round, lineJoin: .round))
                            .blur(radius: model.locked ? 9 : 6).allowsHitTesting(false)
                        trimmed.stroke(model.accent.opacity(0.5),
                                       style: StrokeStyle(lineWidth: model.locked ? 6 : 4.5, lineCap: .round, lineJoin: .round))
                            .blur(radius: 3).allowsHitTesting(false)
                        trimmed.stroke(model.accent, style: StrokeStyle(lineWidth: model.locked ? 3 : 2.4, lineCap: .round, lineJoin: .round))
                            .shadow(color: model.accent.opacity(model.locked ? 0.8 : 0.4), radius: model.locked ? 12 : 6)
                            .allowsHitTesting(false)
                    }
                    .scaleEffect(1 + 0.05 * lockBeat)
                }
                if model.inputStroke.isEmpty {
                    VStack(spacing: 10) {
                        Image(systemName: "scribble.variable").font(.system(size: 32)).foregroundStyle(P.dim.opacity(0.55))
                        Text(model.placed.isEmpty ? "Sketch a curve, then tap Solve" : "Draw another stroke to build the composition")
                            .font(.system(size: 14, weight: .medium)).foregroundStyle(P.dim).multilineTextAlignment(.center)
                    }
                    .allowsHitTesting(false)
                }
                if model.noFit {
                    VStack(spacing: 8) {
                        Image(systemName: "wand.and.stars.inverse").font(.system(size: 26)).foregroundStyle(P.dim.opacity(0.75))
                        Text("No clean law found").font(.system(size: 15, weight: .semibold)).foregroundStyle(P.ink.opacity(0.9))
                        Text("Try a smoother, simpler stroke").font(.system(size: 13)).foregroundStyle(P.dim)
                    }
                    .padding(.horizontal, 20).padding(.vertical, 15)
                    .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 16, style: .continuous))
                    .overlay(RoundedRectangle(cornerRadius: 16, style: .continuous).stroke(P.dim.opacity(0.25), lineWidth: 1))
                    .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .center)
                    .transition(.opacity.combined(with: .scale(scale: 0.95)))
                    .allowsHitTesting(false)
                }
                if model.champion != nil {
                    MathgramHUD(model: model, maxEqWidth: max(220, geo.size.width - 76)).padding(.bottom, 18)
                        .scaleEffect(1 + 0.03 * lockBeat)
                }
            }
            .accessibilityElement()
            .accessibilityLabel("Drawing canvas")
            .accessibilityHint("Draw a curve with one finger, then tap Solve")
            .onChange(of: model.isSolving) { _, solving in
                if solving {
                    withAnimation(.easeInOut(duration: 0.7).repeatForever(autoreverses: true)) { searchPulse = true }
                } else {
                    withAnimation(.easeOut(duration: 0.3)) { searchPulse = false }
                }
            }
            .onChange(of: model.noFit) { _, nf in
                if nf { shakeUnit = 0; withAnimation(.easeInOut(duration: 0.5)) { shakeUnit = 1 } }
                else { shakeUnit = 0 }
            }
            .onChange(of: model.locked) { _, isLocked in
                if isLocked {
                    withAnimation(.spring(response: 0.26, dampingFraction: 0.5)) { lockBeat = 1 }
                    DispatchQueue.main.asyncAfter(deadline: .now() + 0.28) {
                        withAnimation(.easeOut(duration: 0.5)) { lockBeat = 0 }
                    }
                } else {
                    lockBeat = 0
                }
            }
        }
    }
}

private indirect enum MathNode {
    case num(Double)
    case v(String)
    case fn(String, MathNode)
    case sum([(Character, MathNode)])
    case prod([MathNode])
    case frac(MathNode, MathNode)
    case pow(MathNode, MathNode)
    case neg(MathNode)
    case group(MathNode)
}

private enum MathParser {
    static func parse(_ raw: String) -> MathNode? {
        let s = raw.replacingOccurrences(of: "X1", with: "t").replacingOccurrences(of: "x1", with: "t")
        var toks: [String] = []
        let chars = Array(s); var i = 0
        while i < chars.count {
            let c = chars[i]
            if c.isWhitespace { i += 1; continue }
            if c.isNumber || (c == "." && i + 1 < chars.count && chars[i + 1].isNumber) {
                var j = i
                while j < chars.count && (chars[j].isNumber || chars[j] == "." ||
                    ((chars[j] == "e" || chars[j] == "E") && j + 1 < chars.count) ||
                    ((chars[j] == "+" || chars[j] == "-") && j > i && (chars[j - 1] == "e" || chars[j - 1] == "E"))) { j += 1 }
                toks.append(String(chars[i..<j])); i = j
            } else if c.isLetter || c == "_" {
                var j = i
                while j < chars.count && (chars[j].isLetter || chars[j].isNumber || chars[j] == "_") { j += 1 }
                toks.append(String(chars[i..<j])); i = j
            } else { toks.append(String(c)); i += 1 }
        }
        var p = 0
        func peek() -> String? { p < toks.count ? toks[p] : nil }
        func eat() -> String? { defer { p += 1 }; return peek() }
        func parseSum() -> MathNode {
            var terms: [(Character, MathNode)] = [("+", parseProd())]
            while let t = peek(), t == "+" || t == "-" { _ = eat(); terms.append((Character(t), parseProd())) }
            return terms.count == 1 ? terms[0].1 : .sum(terms)
        }
        func parseProd() -> MathNode {
            var factors: [MathNode] = [parsePow()]
            var divs: [MathNode] = []
            while let t = peek(), t == "*" || t == "/" {
                _ = eat()
                if t == "/" { divs.append(parsePow()) } else { factors.append(parsePow()) }
            }
            var node: MathNode = factors.count == 1 ? factors[0] : .prod(factors)
            for d in divs { node = .frac(node, d) }
            return node
        }
        func parsePow() -> MathNode {
            let base = parseAtom()
            if peek() == "^" { _ = eat(); return .pow(base, parsePow()) }
            return base
        }
        func parseAtom() -> MathNode {
            guard let t = peek() else { return .num(0) }
            if t == "-" { _ = eat(); return .neg(parseAtom()) }
            if t == "+" { _ = eat(); return parseAtom() }
            if t == "(" { _ = eat(); let e = parseSum(); if peek() == ")" { _ = eat() }; return .group(e) }
            if let d = Double(t) { _ = eat(); return .num(d) }
            _ = eat()
            if peek() == "(" {
                _ = eat(); let arg = parseSum(); if peek() == ")" { _ = eat() }
                return .fn(t, arg)
            }
            return .v(t)
        }
        guard !toks.isEmpty else { return nil }
        return parseSum()
    }
}

private struct MathExpr: View {
    let raw: String
    var size: CGFloat = 15
    var body: some View {
        if let node = MathParser.parse(raw) { MathRender(node: node, size: size) }
        else { Text(raw).font(.system(size: size, design: .monospaced)) }
    }
}

private struct MathRender: View {
    let node: MathNode
    let size: CGFloat
    private static let fns: Set<String> = ["sin", "cos", "tan", "exp", "log", "ln", "sqrt", "abs", "atan", "asin", "acos", "sinh", "cosh", "tanh"]

    var body: some View {
        switch node {
        case .num(let d):
            Text(fmt(d)).font(.system(size: size, weight: .medium, design: .serif))
        case .v(let name):
            Text(name).font(.system(size: size, weight: .medium, design: .serif)).italic()
        case .group(let inner):
            if case .sum = inner {
                HStack(spacing: 1) { paren("("); MathRender(node: inner, size: size); paren(")") }
            } else {
                MathRender(node: inner, size: size)
            }
        case .neg(let inner):
            HStack(spacing: 1) { Text("−").font(.system(size: size, design: .serif)); MathRender(node: inner, size: size) }
        case .fn(let name, let arg):
            let inner: MathNode = { if case .group(let g) = arg { return g }; return arg }()
            HStack(spacing: 2) {
                Text(name).font(.system(size: size, weight: .medium, design: .serif))
                paren("("); MathRender(node: inner, size: size); paren(")")
            }
        case .pow(let base, let exp):
            HStack(alignment: .firstTextBaseline, spacing: 0) {
                MathRender(node: base, size: size)
                MathRender(node: exp, size: size * 0.72).baselineOffset(size * 0.36)
            }
        case .frac(let n, let d):
            VStack(spacing: 2.5) {
                MathRender(node: n, size: size * 0.92)
                Rectangle().fill(P.ink).frame(height: 1.5).padding(.horizontal, 1)
                MathRender(node: d, size: size * 0.92)
            }.fixedSize()
        case .prod(let factors):
            HStack(spacing: 3) {
                ForEach(Array(factors.enumerated()), id: \.offset) { idx, f in
                    if idx > 0 { Text("·").font(.system(size: size * 0.8, design: .serif)).foregroundStyle(P.dim) }
                    MathRender(node: f, size: size)
                }
            }
        case .sum(let terms):
            HStack(spacing: 4) {
                ForEach(Array(terms.enumerated()), id: \.offset) { idx, t in
                    if idx > 0 { Text(t.0 == "-" ? "−" : "+").font(.system(size: size, design: .serif)) }
                    MathRender(node: t.1, size: size)
                }
            }
        }
    }
    private func paren(_ s: String) -> some View { Text(s).font(.system(size: size, design: .serif)).foregroundStyle(P.dim) }
    private func fmt(_ d: Double) -> String {
        if !d.isFinite { return "∞" }
        if d == d.rounded() && abs(d) < 1e6 { return String(Int(d)) }
        return String(format: "%.3g", d)
    }
}

private struct MathgramHUD: View {
    @ObservedObject var model: SolveModel
    var maxEqWidth: CGFloat = 320
    @State private var expanded = false
    @State private var eqNatural: CGSize = .zero
    #if canImport(UIKit)
    @State private var shareItem: ShareItem?
    #endif

    private func share() {
        #if canImport(UIKit)
        guard let c = model.champion else { return }
        let e = SketchExport(curvePoints: c.curvePoints, exprX: c.exprX, exprY: c.exprY, r2: c.r2, nodes: c.nodes, accent: model.accent)
        if let img = renderExportImage(e) { shareItem = ShareItem(image: img) }
        #endif
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            if let c = model.champion {
                let scale = (eqNatural.width > maxEqWidth && eqNatural.width > 0) ? maxEqWidth / eqNatural.width : 1
                VStack(alignment: .leading, spacing: 8) {
                    mathLine("x(t)", c.exprX, size: 15)
                    mathLine("y(t)", c.exprY, size: 15)
                }
                .fixedSize()
                .background(GeometryReader { g in
                    Color.clear
                        .onAppear { eqNatural = g.size }
                        .onChange(of: g.size) { _, s in eqNatural = s }
                })
                .scaleEffect(scale, anchor: .topLeading)
                .frame(width: eqNatural.width * scale, height: eqNatural.height * scale, alignment: .topLeading)
            }
            HStack(spacing: 8) {
                chip(String(format: "%04.1fs", Double(model.elapsedMs) / 1000.0), icon: model.isSolving ? "hourglass" : "checkmark.seal.fill",
                     tint: model.isSolving ? P.dim : model.accent)
                if let c = model.champion { chip(String(format: "R² %.3f", c.r2), tint: c.r2 > 0.95 ? model.accent : P.dim) }
                if let c = model.champion { chip("\(c.nodes) nodes", tint: model.accent) }
                Spacer()
                if model.locked, model.champion != nil {
                    Button(action: model.placeCurrentAndReset) {
                        Image(systemName: "plus.circle.fill").font(.system(size: 15, weight: .semibold)).foregroundStyle(P.accent)
                    }.buttonStyle(.plain).accessibilityLabel("Add another stroke")
                }
                #if canImport(UIKit)
                if model.locked, model.champion != nil {
                    Button(action: share) {
                        Image(systemName: "square.and.arrow.up").font(.system(size: 13, weight: .semibold)).foregroundStyle(P.accent)
                    }.buttonStyle(.plain).accessibilityLabel("Share result card")
                }
                #endif
                Image(systemName: expanded ? "chevron.down" : "chevron.up").font(.system(size: 11, weight: .bold)).foregroundStyle(P.dim)
            }
            if expanded, let c = model.champion {
                Divider().overlay(P.hairline)
                HStack {
                    metaCol("FAMILY", "PARAMETRIC")
                    metaCol("COMPLEXITY", "\(c.nodes)")
                    metaCol("MDL", String(format: "%.0f", c.dlBits))
                    metaCol("RESIDUAL", String(format: "%.4f", c.shapeErr))
                }
                Text("by lljm").font(.system(size: 11)).foregroundStyle(P.dim)
                if !model.netDiag.isEmpty {
                    Text(model.netDiag).font(.system(size: 9, design: .monospaced)).foregroundStyle(P.dim)
                        .lineLimit(2).fixedSize(horizontal: false, vertical: true)
                }
            }
        }
        .foregroundStyle(P.ink)
        .padding(.horizontal, 20).padding(.vertical, 14)
        .frame(maxWidth: 540)
        .background(RoundedRectangle(cornerRadius: 20).fill(.ultraThinMaterial)
            .overlay(RoundedRectangle(cornerRadius: 20).stroke(.white.opacity(0.13)))
            .shadow(color: P.accent.opacity(0.2), radius: 20, y: 8))
        .padding(.horizontal, 16)
        .onTapGesture { withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) { expanded.toggle() } }
        #if canImport(UIKit)
        .sheet(item: $shareItem) { item in ShareSheet(items: [item.image]) }
        #endif
    }

    private func mathLine(_ lhs: String, _ expr: String, size: CGFloat = 15) -> some View {
        HStack(alignment: .center, spacing: 6) {
            (Text(lhs).font(.system(size: size, weight: .medium, design: .serif)).italic()
             + Text(" =").font(.system(size: size, weight: .medium, design: .serif)))
                .foregroundStyle(P.ink)
            MathExpr(raw: expr, size: size).foregroundStyle(P.ink)
        }
    }

    private func chip(_ text: String, icon: String? = nil, tint: Color) -> some View {
        HStack(spacing: 4) {
            if let icon { Image(systemName: icon).font(.system(size: 10, weight: .bold)) }
            Text(text).font(.system(size: 12, weight: .semibold, design: .monospaced))
        }
        .foregroundStyle(tint)
        .padding(.horizontal, 8).padding(.vertical, 3)
        .background(Capsule().fill(.white.opacity(0.05)).overlay(Capsule().stroke(.white.opacity(0.10))))
    }
    private func metaCol(_ k: String, _ v: String) -> some View {
        VStack(alignment: .leading, spacing: 2) {
            Text(k).font(.system(size: 9, weight: .bold)).tracking(0.5).foregroundStyle(P.dim)
            Text(v).font(.system(size: 12, weight: .semibold, design: .monospaced)).foregroundStyle(P.ink)
        }.frame(maxWidth: .infinity, alignment: .leading)
    }
}

private struct SymbolicView: View {
    @ObservedObject var model: SolveModel
    private let cols = [GridItem(.adaptive(minimum: 100), spacing: 12)]

    var body: some View {
        VStack(spacing: 10) {
            HStack {
                Text("SYMBOLIC EVOLUTION").font(.system(size: 11, weight: .bold)).tracking(1.5).foregroundStyle(P.dim)
                Spacer()
                Text("\(model.candidates.count) · gen \(model.generation)")
                    .font(.system(size: 11, design: .monospaced)).foregroundStyle(P.dim)
            }.padding(.horizontal, 20).padding(.top, 12)

            if model.candidates.isEmpty {
                Spacer()
                VStack(spacing: 10) {
                    Image(systemName: "square.grid.3x3").font(.system(size: 30)).foregroundStyle(P.dim.opacity(0.5))
                    Text("Solve to watch the population evolve").font(.system(size: 13)).foregroundStyle(P.dim)
                }
                Spacer()
            } else {
                ScrollView {
                    LazyVGrid(columns: cols, spacing: 12) {
                        ForEach(Array(model.candidates.prefix(24))) { c in
                            EqNode(c: c, isChamp: c.id == model.champion?.id)
                        }
                    }
                    .padding(.horizontal, 20).padding(.bottom, 14)
                    .animation(.spring(response: 0.4, dampingFraction: 0.82), value: model.candidates.map(\.id))
                }
            }
        }
    }
}

private struct EqNode: View {
    let c: ParametricCandidate
    let isChamp: Bool
    var body: some View {
        VStack(spacing: 6) {
            Canvas { ctx, size in
                ctx.stroke(pathThrough(fitPoints(c.curvePoints, in: size)),
                           with: .color(isChamp ? P.curve : P.cyan.opacity(0.75)),
                           style: StrokeStyle(lineWidth: isChamp ? 2.2 : 1.4, lineCap: .round, lineJoin: .round))
            }.frame(height: 72)
            HStack(spacing: 4) {
                Text(String(format: "R² %.2f", c.r2)).font(.system(size: 9, weight: .bold, design: .monospaced))
                    .foregroundStyle(c.r2 > 0.95 ? P.cyan : P.dim)
                Text("· \(c.nodes)n").font(.system(size: 9, design: .monospaced)).foregroundStyle(P.dim)
            }
        }
        .padding(7)
        .background(RoundedRectangle(cornerRadius: 14).fill(.white.opacity(0.04)))
        .overlay(RoundedRectangle(cornerRadius: 14).stroke(isChamp ? P.accent : .white.opacity(0.08), lineWidth: isChamp ? 1.6 : 1))
        .shadow(color: isChamp ? P.accent.opacity(0.4) : .clear, radius: 8)
        .scaleEffect(isChamp ? 1.03 : 1)
    }
}

// MARK: - Dock
private struct Dock: View {
    @Binding var view: Int
    var count: Int = 2
    private func go(_ delta: Int) {
        let next = min(count - 1, max(0, view + delta))
        guard next != view else { return }
        withAnimation(.spring(response: 0.35, dampingFraction: 0.8)) { view = next }
    }
    private func arrow(_ system: String, _ delta: Int) -> some View {
        let enabled = (delta < 0 && view > 0) || (delta > 0 && view < count - 1)
        return Button { go(delta) } label: {
            Image(systemName: system)
                .font(.system(size: 13, weight: .bold))
                .foregroundStyle(enabled ? P.ink.opacity(0.85) : P.dim.opacity(0.35))
                .frame(width: 30, height: 26)
                .contentShape(Rectangle())
        }
        .buttonStyle(.plain).disabled(!enabled)
        .accessibilityLabel(delta < 0 ? "Previous view" : "Next view")
    }
    var body: some View {
        HStack(spacing: 8) {
            arrow("chevron.left", -1)
            HStack(spacing: 10) {
                ForEach(0..<count, id: \.self) { i in
                    Capsule()
                        .fill(view == i ? P.accent : .white.opacity(0.20))
                        .frame(width: view == i ? 22 : 7, height: 7)
                        .onTapGesture { withAnimation(.spring(response: 0.35, dampingFraction: 0.8)) { view = i } }
                }
            }
            arrow("chevron.right", 1)
        }
        .padding(.vertical, 14)
    }
}

// MARK: - Gallery
private struct GalleryView: View {
    @ObservedObject var store: GalleryStore
    private let cols = [GridItem(.adaptive(minimum: 108), spacing: 14)]
    #if canImport(UIKit)
    @State private var shareItem: ShareItem?
    #endif

    private func share(_ s: SavedSketch) {
        #if canImport(UIKit)
        let e = SketchExport(curvePoints: s.curvePoints, exprX: s.exprX, exprY: s.exprY, r2: s.r2, nodes: s.nodes, accent: s.accent)
        if let img = renderExportImage(e) { shareItem = ShareItem(image: img) }
        #endif
    }

    var body: some View {
        VStack(spacing: 10) {
            HStack {
                Text("GALLERY").font(.system(size: 11, weight: .bold)).tracking(1.5).foregroundStyle(P.dim)
                Spacer()
                Text("\(store.sketches.count) saved").font(.system(size: 11, design: .monospaced)).foregroundStyle(P.dim)
            }.padding(.horizontal, 20).padding(.top, 12)

            if store.sketches.isEmpty {
                Spacer()
                VStack(spacing: 10) {
                    Image(systemName: "square.grid.2x2").font(.system(size: 30)).foregroundStyle(P.dim.opacity(0.5))
                    Text("Solved curves collect here").font(.system(size: 13)).foregroundStyle(P.dim)
                }
                Spacer()
            } else {
                ScrollView {
                    LazyVGrid(columns: cols, spacing: 14) {
                        ForEach(store.sketches) { s in
                            GalleryCard(sketch: s)
                                .contextMenu {
                                    Button { share(s) } label: { Label("Share", systemImage: "square.and.arrow.up") }
                                    Button(role: .destructive) { store.delete(s) } label: { Label("Delete", systemImage: "trash") }
                                }
                        }
                    }
                    .padding(.horizontal, 20).padding(.bottom, 16).padding(.top, 4)
                    .animation(.spring(response: 0.4, dampingFraction: 0.82), value: store.sketches.map(\.id))
                }
            }
        }
        #if canImport(UIKit)
        .sheet(item: $shareItem) { item in ShareSheet(items: [item.image]) }
        #endif
    }
}

private struct GalleryCard: View {
    let sketch: SavedSketch
    var body: some View {
        let accent = sketch.accent
        VStack(spacing: 7) {
            Canvas { ctx, size in
                let pts = fitPoints(sketch.curvePoints, in: size)
                let path = pathThrough(pts)
                ctx.addFilter(.shadow(color: accent.opacity(0.6), radius: 6))
                ctx.stroke(path, with: .color(accent), style: StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
            }
            .frame(height: 92)
            HStack(spacing: 5) {
                Text(String(format: "R² %.2f", sketch.r2)).font(.system(size: 9, weight: .bold, design: .monospaced))
                    .foregroundStyle(sketch.r2 > 0.95 ? P.cyan : P.dim)
                Text("· \(sketch.nodes)n").font(.system(size: 9, design: .monospaced)).foregroundStyle(P.dim)
            }
        }
        .padding(9)
        .background(RoundedRectangle(cornerRadius: 16).fill(.white.opacity(0.04)))
        .overlay(RoundedRectangle(cornerRadius: 16).stroke(accent.opacity(0.55), lineWidth: 1.4))
        .shadow(color: accent.opacity(0.28), radius: 10)
    }
}

struct SketchExport {
    let curvePoints: [CGPoint]
    let exprX: String
    let exprY: String
    let r2: Double
    let nodes: Int
    let accent: Color
}

struct ExportCard: View {
    let e: SketchExport
    var body: some View {
        VStack(spacing: 0) {
            ZStack {
                Rectangle().fill(P.stage)
                Canvas { ctx, size in
                    let step: CGFloat = 40
                    var g = Path()
                    var x: CGFloat = 0; while x < size.width { g.move(to: .init(x: x, y: 0)); g.addLine(to: .init(x: x, y: size.height)); x += step }
                    var y: CGFloat = 0; while y < size.height { g.move(to: .init(x: 0, y: y)); g.addLine(to: .init(x: size.width, y: y)); y += step }
                    ctx.stroke(g, with: .color(.white.opacity(0.05)), lineWidth: 1)
                    let pts = fitPoints(e.curvePoints, in: size, pad: 0.72)
                    let path = pathThrough(pts)
                    ctx.addFilter(.shadow(color: e.accent.opacity(0.7), radius: 18))
                    ctx.stroke(path, with: .color(e.accent), style: StrokeStyle(lineWidth: 5, lineCap: .round, lineJoin: .round))
                    ctx.stroke(path, with: .color(.white.opacity(0.85)), style: StrokeStyle(lineWidth: 1.6, lineCap: .round, lineJoin: .round))
                }
            }
            .frame(height: 300)

            VStack(alignment: .leading, spacing: 8) {
                ViewThatFits(in: .horizontal) {
                    ForEach([18, 16, 14, 12, 10, 9], id: \.self) { s in
                        eqBlock(CGFloat(s))
                    }
                }
                .frame(maxWidth: .infinity, alignment: .leading)
                HStack(spacing: 8) {
                    chip(String(format: "R² %.3f", e.r2), e.accent)
                    chip("\(e.nodes) nodes", e.accent)
                    Spacer()
                    HStack(spacing: 5) {
                        BrandGlyph(stroke: P.ink, smile: e.accent).frame(width: 15, height: 15)
                        Text("GraphArt").font(.system(size: 13, weight: .bold))
                    }.foregroundStyle(P.ink)
                }
            }
            .padding(16)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(P.bgPage)
        }
        .frame(width: 360, height: 440)
        .background(P.bgPage)
        .environment(\.colorScheme, .dark)
    }

    @ViewBuilder private func eqBlock(_ size: CGFloat) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            line("x(t)", e.exprX, size)
            line("y(t)", e.exprY, size)
        }
    }
    private func line(_ lhs: String, _ expr: String, _ size: CGFloat) -> some View {
        HStack(alignment: .center, spacing: 6) {
            (Text(lhs).font(.system(size: size, weight: .medium, design: .serif)).italic()
             + Text(" =").font(.system(size: size, weight: .medium, design: .serif))).foregroundStyle(P.ink)
            MathExpr(raw: expr, size: size).foregroundStyle(P.ink)
        }
    }
    private func chip(_ t: String, _ tint: Color) -> some View {
        Text(t).font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundStyle(tint)
            .padding(.horizontal, 9).padding(.vertical, 4)
            .background(Capsule().fill(.white.opacity(0.06)).overlay(Capsule().stroke(tint.opacity(0.4))))
    }
}

@MainActor func renderExportImage(_ e: SketchExport) -> UIImage? {
    let renderer = ImageRenderer(content: ExportCard(e: e))
    renderer.scale = 3
    #if canImport(UIKit)
    return renderer.uiImage
    #else
    return nil
    #endif
}

#if canImport(UIKit)
struct ShareItem: Identifiable { let id = UUID(); let image: UIImage }

struct ShareSheet: UIViewControllerRepresentable {
    let items: [Any]
    func makeUIViewController(context: Context) -> UIActivityViewController {
        UIActivityViewController(activityItems: items, applicationActivities: nil)
    }
    func updateUIViewController(_ vc: UIActivityViewController, context: Context) {}
}
#endif
