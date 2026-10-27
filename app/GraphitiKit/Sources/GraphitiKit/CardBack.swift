import SwiftUI

public struct MathCardBack: View {
    let stats: SolveSession.Stats
    let latex: String
    let author: String
    var date = Date()

    public init(stats: SolveSession.Stats, latex: String, author: String) {
        self.stats = stats; self.latex = latex; self.author = author
    }

    private var colour: Color { EquationColor.color(for: latex) }

    public var body: some View {
        VStack(alignment: .leading, spacing: Space.m) {
            FitMath(latex, size: 18, color: EquationColor.platformColor(for: latex))


            Grid(horizontalSpacing: Space.s, verticalSpacing: Space.s) {
                GridRow {
                    tile(String(format: "%.1f s", stats.seconds), "solve time")
                    tile(stats.equationsTried.formatted(), "equations tried")
                    tile("\(stats.generations)", "generations")
                }
                GridRow {
                    tile("\(stats.terms)", "terms")
                    tile(String(format: "%.1f%%", stats.fit * 100), "fit")
                    if stats.exactConstants > 0 { tile("\(stats.exactConstants)", "exact constant") } else { Color.clear }
                }
            }

            HStack(spacing: Space.l) {
                ForEach(stats.properties, id: \.0) { name, value in
                    HStack(alignment: .center, spacing: 5) {
                        Text(name).font(.caption).foregroundStyle(.secondary)
                        MathText(value, size: 15)
                    }
                }
            }
            .frame(maxWidth: .infinity)

            Spacer(minLength: 0)
            HStack {
                Text("@\(author.isEmpty ? "me" : author)").foregroundStyle(.secondary)
                Spacer()
                Text("Graphiti").fontWeight(.semibold).foregroundStyle(.white)
            }
            .font(Typeface.meta)
        }
        .padding(.horizontal, Space.s + Radius.cardPadding)
        .padding(.top, Space.l).padding(.bottom, Radius.cardPadding + Space.xs)
        .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .top)
        .cardSurface(radius: Radius.heroCard)
        .accessibilityElement(children: .combine)
    }

    private var sparkline: some View {
        GeometryReader { geo in
            let pts = stats.fitHistory
            let maxGen = Double(pts.last?.0 ?? 1)
            let score = { (f: Double) -> Double in -log10(max(1e-4, 1 - f)) }
            let lo = score(pts.first?.1 ?? 0.4), hi = score(pts.last?.1 ?? 0.99)
            let point = { (g: Int, f: Double) -> CGPoint in
                CGPoint(x: geo.size.width * Double(g - 1) / max(1, maxGen - 1),
                        y: geo.size.height * (1 - (score(f) - lo) / max(1e-6, hi - lo)))
            }
            ZStack(alignment: .topLeading) {
                Path { p in
                    for (i, (g, f)) in pts.enumerated() {
                        if i == 0 { p.move(to: point(g, f)) } else { p.addLine(to: point(g, f)) }
                    }
                }
                .stroke(colour, style: StrokeStyle(lineWidth: 2, lineCap: .round, lineJoin: .round))
                if let last = pts.last {
                    Circle().fill(colour).frame(width: 7, height: 7).position(point(last.0, last.1))
                }
                Text("fit by generation").font(.caption2).foregroundStyle(.tertiary)
                    .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .bottomTrailing)
            }
        }
        .padding(Space.s)
        .background(Theme.ink, in: ConcentricRectangle(corners: .concentric, isUniform: true))
    }

    private func tile(_ value: String, _ label: String) -> some View {
        VStack(spacing: 2) {
            Text(value).font(.headline).monospacedDigit().lineLimit(1).minimumScaleFactor(0.7)
            Text(label).font(.caption2).foregroundStyle(.secondary).lineLimit(1).minimumScaleFactor(0.8)
        }
        .frame(maxWidth: .infinity).padding(.vertical, Space.s)
        .background(Theme.surfaceRaised, in: .rect(cornerRadius: 12))
    }
}

public struct FlipCard<Front: View, Back: View>: View {
    let front: Front
    let back: Back
    @Binding var flipped: Bool
    @State private var peek = 0.0
    @Environment(\.accessibilityReduceMotion) private var reduceMotion

    public init(flipped: Binding<Bool>, @ViewBuilder front: () -> Front, @ViewBuilder back: () -> Back) {
        _flipped = flipped; self.front = front(); self.back = back()
    }

    public var body: some View {
        front
            .opacity(flipped ? 0 : 1)
            .overlay {
                back
                    .rotation3DEffect(.degrees(reduceMotion ? 0 : 180), axis: (0, 1, 0))
                    .opacity(flipped ? 1 : 0)
            }
            .rotation3DEffect(.degrees((flipped && !reduceMotion ? 180 : 0) + peek), axis: (0, 1, 0), perspective: 0.6)
            .contentShape(Rectangle())
            .onTapGesture { withAnimation(.spring(response: 0.55, dampingFraction: 0.82)) { flipped.toggle() } }
            .sensoryFeedback(.impact(weight: .light), trigger: flipped)
            .accessibilityAction(named: flipped ? "Show the curve" : "Show stats") { flipped.toggle() }
            .task {
                guard !reduceMotion, !flipped else { return }
                try? await Task.sleep(for: .milliseconds(900))
                withAnimation(.spring(response: 0.35, dampingFraction: 0.6)) { peek = 16 }
                try? await Task.sleep(for: .milliseconds(260))
                withAnimation(.spring(response: 0.45, dampingFraction: 0.7)) { peek = 0 }
            }
    }
}
