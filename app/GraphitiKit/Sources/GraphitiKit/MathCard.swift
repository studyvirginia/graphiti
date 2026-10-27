import SwiftUI

public struct MathCard: View {
    var params: CurveParams
    var latex: String
    var author: String
    var points: [CGPoint] = []
    var stroke: [CGPoint] = []
    var date = Date()
    @State private var showStroke = false

    public init(params: CurveParams, latex: String, author: String, points: [CGPoint] = [], stroke: [CGPoint] = []) {
        self.params = params; self.latex = latex; self.author = author; self.points = points; self.stroke = stroke
    }

    public var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            PlotView(radius: Radius.heroPlot) {
                if points.count > 1 {
                    let box = points + stroke
                    if showStroke && stroke.count > 1 {
                        FittedPolyline(points: stroke, frame: box).stroke(.white.opacity(0.35), style: StrokeStyle(lineWidth: 6, lineCap: .round, lineJoin: .round))
                            .transition(.opacity)
                    }
                    FittedPolyline(points: points, frame: box).stroke(EquationColor.color(for: latex), style: StrokeStyle(lineWidth: 3.5, lineCap: .round, lineJoin: .round))
                } else {
                    CurveShape(params: params, layout: .card).stroke(EquationColor.color(for: latex), style: StrokeStyle(lineWidth: 3.5, lineCap: .round))
                }
            }
            .frame(height: 226)
            .contentShape(Rectangle())
            .onTapGesture { if stroke.count > 1 { withAnimation(.easeInOut(duration: 0.2)) { showStroke.toggle() } } }
            .accessibilityAction(named: "Show your stroke") { showStroke.toggle() }

            FitMath(latex, size: Typeface.cardEquationSize).padding(.top, Space.l).padding(.horizontal, Space.s)
            HStack {
                Text("@\(author.isEmpty ? "me" : author)").foregroundStyle(Theme.secondary)
                Spacer()
                Text("Graphiti").fontWeight(.semibold).foregroundStyle(.white)
            }
            .font(Typeface.meta)
            .padding(.horizontal, Space.s).padding(.top, Space.m).padding(.bottom, Space.xs)
        }
        .padding(Radius.cardPadding)
        .cardSurface(radius: Radius.heroCard)
        .shadow(color: .black.opacity(0.6), radius: 35, y: 30)
    }
}

public struct SharePopup<Background: View>: View {
    let background: Background
    let card: MathCard
    let back: MathCardBack?
    var onClose: () -> Void
    var onSave: () -> Void
    var onPost: () -> Void
    var posted: Bool
    @State private var flipped: Bool
    @State private var saved = false

    public init(card: MathCard, back: MathCardBack? = nil, startFlipped: Bool = false,
                onClose: @escaping () -> Void = {}, onSave: @escaping () -> Void = {}, onPost: @escaping () -> Void = {},
                posted: Bool = false, @ViewBuilder background: () -> Background) {
        self.card = card; self.back = back; self.onClose = onClose; self.background = background()
        self.onSave = onSave; self.onPost = onPost; self.posted = posted
        _flipped = State(initialValue: startFlipped)
    }

    public var body: some View {
        ZStack {
            background.blur(radius: 26).overlay(Color.black.opacity(0.4))
            VStack(spacing: 0) {
                HStack { Spacer(); IconButton("xmark", "Close", action: onClose) }
                    .padding(.horizontal, Space.l)
                Spacer(minLength: Space.l)
                Group {
                    if let back {
                        FlipCard(flipped: $flipped) { card.shadow(radius: 0) } back: { back }
                            .shadow(color: .black.opacity(0.6), radius: 35, y: 30)
                    } else {
                        card
                    }
                }
                .padding(.horizontal, 28)
                HStack(spacing: 18) {
                    destination("message.fill", "Messages")
                    destination("camera", "Instagram")
                    destination("function", "Desmos")
                    destination(saved ? "checkmark" : "square.and.arrow.down", saved ? "Saved" : "Save") { onSave(); saved = true }
                    destination("ellipsis", "More")
                }
                .padding(.top, 36)
                PrimaryButton(posted ? "Posted to Discover" : "Post to Discover", systemImage: posted ? "checkmark" : nil, action: onPost)
                    .disabled(posted)
                    .padding(.horizontal, 28).padding(.top, 28)
                Spacer(minLength: Space.l)
            }
        }
    }

    private func destination(_ symbol: String, _ label: String, action: @escaping () -> Void = {}) -> some View {
        VStack(spacing: Space.s) {
            IconButton(symbol, label, size: 52, action: action)
            Text(label).font(Typeface.label).foregroundStyle(Theme.secondary)
        }
    }
}


struct FittedPolyline: Shape {
    var points: [CGPoint]
    var frame: [CGPoint]? = nil
    func path(in rect: CGRect) -> Path {
        let xs = (frame ?? points).map(\.x), ys = (frame ?? points).map(\.y)
        guard let x0 = xs.min(), let x1 = xs.max(), let y0 = ys.min(), let y1 = ys.max() else { return Path() }
        let inset = rect.insetBy(dx: rect.width * 0.1, dy: rect.height * 0.12)
        let s = min(inset.width / max(x1 - x0, 1e-6), inset.height / max(y1 - y0, 1e-6))
        let ox = inset.midX - (x0 + x1) / 2 * s, oy = inset.midY - (y0 + y1) / 2 * s
        var p = Path()
        for (i, q) in points.enumerated() {
            let pt = CGPoint(x: q.x * s + ox, y: q.y * s + oy)
            if i == 0 { p.move(to: pt) } else { p.addLine(to: pt) }
        }
        return p
    }
}
