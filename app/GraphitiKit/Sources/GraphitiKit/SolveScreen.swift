import SwiftUI

public struct SolveScreen: View {
    @ObservedObject var session: SolveSession
    var onShare: () -> Void = {}
    var onInfo: () -> Void = {}
    var onStroke: ([CGPoint]) -> Void = { _ in }
    @State private var hint: Double
    @State private var hintOpacity: Double = 1
    @Environment(\.accessibilityReduceMotion) private var reduceMotion
    @Environment(\.scenePhase) private var scenePhase
    private let showsHint: Bool

    private let layout = PlaneLayout.solve
    @State private var strokePts: [CGPoint] = []
    @State private var glide = false
    @Environment(\.accessibilityDifferentiateWithoutColor) private var differentiateWithoutColor
    @Environment(\.openURL) private var openURL

    public init(session: SolveSession, onShare: @escaping () -> Void = {}, onInfo: @escaping () -> Void = {},
                onStroke: @escaping ([CGPoint]) -> Void = { _ in },
                firstLaunchHint: Double? = nil) {
        self.session = session; self.onShare = onShare; self.onInfo = onInfo
        self.onStroke = onStroke
        self.showsHint = firstLaunchHint != nil
        self._hint = State(initialValue: firstLaunchHint ?? 0)
    }

    public var body: some View {
        ZStack {
            Theme.ink.ignoresSafeArea()
            if session.engineMap.count < 4 { ManimPlane(layout: layout).ignoresSafeArea().transition(.opacity) }
            Group {
                if session.isEmpty { if showsHint { exampleStroke } } else if session.liveStroke.isEmpty { demoCurves }
            }
            .drawingGroup()
            .mask(fade)
            .ignoresSafeArea()
            .accessibilityElement()
            .accessibilityLabel(session.isEmpty ? "Empty graph" : "Graph of your curve")
            .accessibilityChartDescriptor(CurveAudioGraph(params: session.step.params, title: "Your curve", points: session.liveCurvePoints))
            VStack(spacing: 14) {
                Spacer()
                if session.isEmpty { emptyEquation } else { statusRow; if awaitingFirstBurst { pendingEquation } else { equation } }
            }
            .padding(.bottom, Space.xxl)
            .padding(.horizontal, Space.l)
        }
        .task { if scenePhase == .active { EnginePrepare.start(); EngineSelfTest.startIfRequested() }; await SolveSession.warmUp() }
        .onChange(of: scenePhase) { _, p in if p == .active { EnginePrepare.start(); EngineSelfTest.startIfRequested() } else if p == .background { EnginePrepare.stopRunning() } }
        .onTapGesture {
            if session.explained != nil { session.explain(nil) }
            else if session.locked && !session.liveStroke.isEmpty { withAnimation(.easeInOut(duration: 0.2)) { session.showStroke.toggle() } }
        }
        .contentShape(Rectangle())
        .overlay(alignment: .topLeading) {
            GeometryReader { geo in
            ZStack(alignment: .topLeading) {
                let pl = placement(in: geo.size)
                if session.engineMap.count == 4 {
                    GaugePlane(map: session.engineMap, fadeStart: layout.fadeStart, scale: pl.s, shift: pl.offset)
                }
                if !session.liveStroke.isEmpty {
                    liveCurves(width: 1 / pl.s).scaleEffect(pl.s, anchor: .topLeading).offset(pl.offset)
                }
                if strokePts.count > 1 {
                    Path { p in p.addLines(strokePts) }
                        .stroke(.white, style: StrokeStyle(lineWidth: 3.5, lineCap: .round, lineJoin: .round))
                }
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .topLeading)
            }
            .allowsHitTesting(false)
        }
        .onChange(of: session.liveStroke) { _, s in
            glide = false
            if !s.isEmpty { withAnimation(.spring(response: 0.55, dampingFraction: 0.86).delay(0.05)) { glide = true } }
        }
        .gesture(DragGesture(minimumDistance: 10).onChanged { v in if strokePts.isEmpty { strokePts = [v.startLocation] }; strokePts.append(v.location) }.onEnded { _ in
            if strokePts.count > 1 { withAnimation(.easeInOut(duration: 0.3)) { session.liveStroke = strokePts } }
            onStroke(strokePts); strokePts = []
        }, including: session.isEmpty ? .all : .subviews)
        #if os(iOS)
        .sensoryFeedback(.impact(weight: .light, intensity: 0.6), trigger: session.step.tokens.map(\.latex))
        .sensoryFeedback(.success, trigger: session.locked)
        #endif
    }

    private var fade: some View {
        GeometryReader { geo in
            LinearGradient(stops: [.init(color: .white, location: 0),
                                   .init(color: .white, location: layout.fadeStart),
                                   .init(color: .clear, location: layout.fadeStart + 0.08)],
                           startPoint: .top, endPoint: .bottom)
                .frame(height: geo.size.height)
        }
    }

    // MARK: Curves

    private func placement(in size: CGSize) -> (s: CGFloat, offset: CGSize) {
        let pts = session.liveStroke
        guard glide, let x0 = pts.map(\.x).min(), let x1 = pts.map(\.x).max(),
              let y0 = pts.map(\.y).min(), let y1 = pts.map(\.y).max() else { return (1, .zero) }
        let top: CGFloat = 120, bottom = size.height * layout.fadeStart - 24
        let bw = max(x1 - x0, 1), bh = max(y1 - y0, 1)
        let s = min(max(min(0.78 * size.width / bw, 0.78 * (bottom - top) / bh), 0.5), 2.5)
        let cx = (x0 + x1) / 2, cy = (y0 + y1) / 2
        return (s, CGSize(width: size.width / 2 - cx * s, height: (top + bottom) / 2 - cy * s))
    }

    private func liveCurves(width w: CGFloat) -> some View {
        let st = session.step
        let champ = Path { p in
            let n = min(st.xs.count, st.samples.count)
            for i in 0..<n {
                let pt = CGPoint(x: st.xs[i], y: -st.samples[i])
                if i == 0 { p.move(to: pt) } else { p.addLine(to: pt) }
            }
        }
        let stroke = Path { p in p.addLines(session.liveStroke) }
        return ZStack {
            stroke.stroke(.white.opacity(session.locked ? (session.showStroke ? 0.34 : 0) : 0.34), style: StrokeStyle(lineWidth: (6) * w, lineCap: .round, lineJoin: .round))
            if !session.locked {
                ForEach(Array(st.others.enumerated()), id: \.offset) { _, pts in
                    Path { p in p.addLines(pts) }.stroke(session.colour.opacity(0.16), lineWidth: (1.6) * w)
                }
            }
            champ.stroke(session.colour.opacity(session.locked ? 1 : 0.8),
                         style: StrokeStyle(lineWidth: (session.locked ? 4 : 3.4) * w, lineCap: .round, lineJoin: .round))
                .shadow(color: session.colour.opacity(session.locked ? 0.55 : 0), radius: 7 * w)
            if session.flash > 0 && session.flash < 1 && !reduceMotion {
                champ.trim(from: max(0, session.flash - 0.2), to: session.flash)
                    .stroke(.white, style: StrokeStyle(lineWidth: (5) * w, lineCap: .round))
                    .shadow(color: .white.opacity(0.9), radius: 6 * w)
            }
        }
    }

    private var demoCurves: some View {
        let p = session.step.params
        let progress = Double(session.index) / Double(session.steps.count - 1)
        return ZStack {
            chalk
            if !session.locked {
                ForEach(0..<8, id: \.self) { j in
                    let off = [0.42, -0.36, 0.28, -0.22, 0.5, -0.46, 0.16, -0.12][j] * (1 - progress)
                    CurveShape(params: CurveParams(m: p.m + off * 0.25, a: p.a + off * 0.9, w: p.w + off * 0.35, q: p.q + off * 0.04),
                               layout: layout)
                        .stroke(session.colour.opacity(0.08 + 0.22 * (1 - progress)), lineWidth: 1.6)
                }
            }
            CurveShape(params: p, layout: layout)
                .stroke(session.colour.opacity(session.locked ? 1 : 0.75),
                        style: StrokeStyle(lineWidth: session.locked ? 4 : 3.4, lineCap: .round, lineJoin: .round))
                .shadow(color: session.colour.opacity(session.locked ? 0.55 : 0), radius: 7)
                .opacity(session.explainedTerm == nil ? 1 : 0.3)
            if let term = session.explainedTerm {
                CurveShape(params: term.params, layout: layout)
                    .stroke(session.colour, style: StrokeStyle(lineWidth: 3.5, lineCap: .round,
                                                            dash: differentiateWithoutColor ? [8, 6] : []))
                    .transition(.opacity)
            }
            if session.flash > 0 && session.flash < 1 && !reduceMotion {
                CurveShape(params: p, layout: layout)
                    .trim(from: max(0, session.flash - 0.2), to: session.flash)
                    .stroke(.white, style: StrokeStyle(lineWidth: 5, lineCap: .round))
                    .shadow(color: .white.opacity(0.9), radius: 6)
            }
        }
    }

    private var chalk: some View {
        let shape = CurveShape(params: session.target, wobble: 0.05, layout: layout)
        let o = session.locked ? 0.0 : 1.0
        return ZStack {
            shape.stroke(.white.opacity(0.16 * o), style: StrokeStyle(lineWidth: 9, lineCap: .round, lineJoin: .round))
            shape.stroke(.white.opacity(0.34 * o), style: StrokeStyle(lineWidth: 3.2, lineCap: .round, lineJoin: .round, dash: [5, 1.5, 2, 1.2]))
        }
    }

    private var exampleStroke: some View {
        GeometryReader { geo in
            let f = PlaneFrame(size: geo.size, layout: layout)
            let demo = CurveParams(m: 0.35, a: 0.9, w: 1.6, q: 0)
            let r = f.xRange
            let x = r.lowerBound + 0.6 + (r.upperBound - r.lowerBound - 1.2) * hint
            ZStack {
                CurveShape(params: demo, layout: layout)
                    .trim(from: 0.07, to: 0.07 + 0.86 * hint)
                    .stroke(.white.opacity(0.4), style: StrokeStyle(lineWidth: 3, lineCap: .round, dash: [2, 7]))
                Circle().fill(.white.opacity(0.2)).frame(width: 40, height: 40).position(f.point(x, demo.y(x)))
                Circle().fill(.white.opacity(0.85)).frame(width: 14, height: 14).position(f.point(x, demo.y(x)))
            }
            .opacity(hintOpacity)
        }
        .task {
            guard showsHint, hint == 0 else { return }
            for _ in 0..<2 {
                hint = 0; hintOpacity = 1
                try? await Task.sleep(for: .milliseconds(600))
                let start = Date()
                while hint < 1 {
                    let t = min(1, Date().timeIntervalSince(start) / 1.6)
                    hint = t < 0.5 ? 2 * t * t : 1 - pow(-2 * t + 2, 2) / 2
                    try? await Task.sleep(for: .milliseconds(16))
                }
                try? await Task.sleep(for: .milliseconds(700))
                withAnimation(.easeOut(duration: 0.4)) { hintOpacity = 0 }
                try? await Task.sleep(for: .milliseconds(900))
            }
        }
    }


    private var equation: some View {
        let toks = session.step.tokens
        let cut = toks.firstIndex { $0.id == "sep" }
        return VStack(spacing: 4) {
            if let cut {
                evolving(28, Array(toks[..<cut])).fixedSize().scaleToFitWidth().frame(height: 40)
                evolving(28, Array(toks[(cut + 1)...])).fixedSize().scaleToFitWidth().frame(height: 40)
            } else {
                evolving(Typeface.equationSize, toks).fixedSize().scaleToFitWidth().frame(height: 64)
            }
        }
        .scaleEffect(session.indicate && !reduceMotion ? 1.08 : 1)
        .frame(maxWidth: .infinity)
        .frame(height: 84)
        #if os(iOS)
        .contextMenu {
            Button("Copy Equation", systemImage: "doc.on.doc") { UIPasteboard.general.string = session.plainText }
            Button("Copy LaTeX", systemImage: "chevron.left.forwardslash.chevron.right") { UIPasteboard.general.string = session.latex }
            Button("Open in Desmos", systemImage: "function") {
                UIPasteboard.general.string = session.plainText
                openURL(URL(string: "https://www.desmos.com/calculator")!)
            }
        }
        #endif
        .sensoryFeedback(.selection, trigger: session.explained)
    }

    private func evolving(_ size: CGFloat, _ tokens: [EquationToken]) -> some View {
        EvolvingEquation(tokens: tokens, size: size,
                         highlight: session.indicate ? Set(session.step.tokens.map(\.id)) : (session.locked ? [] : session.changed),
                         tint: session.platformColour, dying: session.dying,
                         colours: explainColours,
                         onTap: session.locked && session.liveStroke.isEmpty ? { id in session.explain(session.term(containing: id)?.id) } : nil)
    }

    private var explainColours: [String: PlatformColor] {
        guard let term = session.explainedTerm else { return [:] }
        var c: [String: PlatformColor] = [:]
        for t in session.step.tokens {
            c[t.id] = term.tokenIDs.contains(t.id) ? session.platformColour : .graphiti(hex: 0xFFFFFF, alpha: 0.4)
        }
        return c
    }

    private var emptyEquation: some View {
        VStack(spacing: 14) {
            Text("Draw any curve").font(Typeface.status).foregroundStyle(.secondary)
                .frame(height: 48)
            MathText("y = \\,?", size: Typeface.equationSize, color: .graphiti(hex: 0xFFFFFF, alpha: 0.35))
                .frame(height: 64)
        }
    }

    private var awaitingFirstBurst: Bool { session.liveSteps.isEmpty && session.liveStroke.count > 1 && !session.locked }

    private var pendingEquation: some View {
        MathText("y = \\,?", size: Typeface.equationSize, color: .graphiti(hex: 0xFFFFFF, alpha: 0.35))
            .frame(height: 64).frame(maxWidth: .infinity).frame(height: 84)
    }

    private var statusRow: some View {
        HStack {
            IconButton("info", "About this equation", action: onInfo).opacity(session.locked ? 1 : 0).accessibilityHidden(!session.locked)
            Spacer()
            if let term = session.explainedTerm {
                (Text(term.name).foregroundStyle(session.colour) + Text(" · " + term.caption).foregroundStyle(.secondary))
                    .font(Typeface.status)
                    .lineLimit(1).minimumScaleFactor(0.8)
                    .transition(.opacity)
            } else {
                TimelineView(.animation(minimumInterval: 0.1, paused: session.startedAt == nil)) { ctx in
                    Text(session.status(at: ctx.date))
                        .font(Typeface.status)
                        .foregroundStyle(.secondary)
                        .accessibilityAddTraits(.updatesFrequently)
                }
            }
            Spacer()
            IconButton("square.and.arrow.up", "Share math card", action: onShare).opacity(session.locked ? 1 : 0).accessibilityHidden(!session.locked)
        }
        .frame(height: 48)
        .allowsHitTesting(session.locked)
    }
}

extension View {
    func scaleToFitWidth(alignment: Alignment = .center) -> some View { modifier(ScaleToFitWidth(alignment: alignment)) }
}

private struct ScaleToFitWidth: ViewModifier {
    var alignment: Alignment = .center
    @State private var natural: CGFloat = 0
    func body(content: Content) -> some View {
        GeometryReader { geo in
            content
                .background(GeometryReader { g in Color.clear.onAppear { natural = g.size.width }
                    .onChange(of: g.size.width) { _, w in natural = w } })
                .scaleEffect(min(1, geo.size.width / max(natural, 1)), anchor: alignment == .leading ? .leading : .center)
                .frame(width: geo.size.width, height: geo.size.height, alignment: alignment)
                .opacity(natural > 0 ? 1 : 0)
        }
    }
}
