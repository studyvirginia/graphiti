import SwiftUI

public struct PlaneLayout: Equatable, Sendable {
    public var originY: CGFloat
    public var unitsAcross: CGFloat
    public var fadeStart: CGFloat
    public init(originY: CGFloat, unitsAcross: CGFloat, fadeStart: CGFloat) {
        self.originY = originY; self.unitsAcross = unitsAcross; self.fadeStart = fadeStart
    }
    public static let standard = PlaneLayout(originY: 0.42, unitsAcross: 8.6, fadeStart: 0.64)
    public static let solve = PlaneLayout(originY: 0.40, unitsAcross: 10, fadeStart: 0.6)
    public static let card = PlaneLayout(originY: 0.5, unitsAcross: 8.6, fadeStart: 1)
}

public struct PlaneFrame: Equatable {
    public var size: CGSize
    public var unit: CGFloat
    public var origin: CGPoint

    public init(size: CGSize, layout: PlaneLayout = .standard) {
        self.size = size
        self.unit = size.width / layout.unitsAcross
        self.origin = CGPoint(x: size.width / 2, y: size.height * layout.originY)
    }

    public func point(_ x: Double, _ y: Double) -> CGPoint {
        CGPoint(x: origin.x + CGFloat(x) * unit, y: origin.y - CGFloat(y) * unit)
    }
    public var xRange: ClosedRange<Double> {
        Double(-origin.x / unit)...Double((size.width - origin.x) / unit)
    }
}

public struct ManimPlane: View {
    var fadeStart: CGFloat
    var layout: PlaneLayout
    public init(fadeStart: CGFloat? = nil, layout: PlaneLayout = .standard) {
        self.fadeStart = fadeStart ?? layout.fadeStart; self.layout = layout
    }

    public var body: some View {
        GeometryReader { geo in
            let f = PlaneFrame(size: geo.size, layout: layout)
            ZStack(alignment: .topLeading) {
                Canvas { ctx, size in
                    drawGrid(ctx, f)
                    drawAxes(ctx, f)
                }
                if geo.size.height > 160 { labels(f) }
            }
            .mask(LinearGradient(stops: [.init(color: .white, location: 0),
                                         .init(color: .white, location: fadeStart),
                                         .init(color: .clear, location: fadeStart + 0.08)],
                                 startPoint: .top, endPoint: .bottom))
        }
        .accessibilityHidden(true)
    }

    @ViewBuilder private func labels(_ f: PlaneFrame) -> some View {
                ForEach([-4, -2, 2, 4], id: \.self) { i in
                    tickLabel("\(i)").position(x: f.point(Double(i), 0).x, y: f.origin.y + 16)
                }
                ForEach([-2, 2], id: \.self) { j in
                    tickLabel("\(j)").position(x: f.origin.x - 14, y: f.point(0, Double(j)).y)
                }
    }

    private func tickLabel(_ s: String) -> some View {
        MathText(s, size: 13, color: .graphiti(hex: 0xFFFFFF, alpha: 0.8))
            .padding(.horizontal, 2)
            .background(Theme.ink.opacity(0.85), in: .rect(cornerRadius: 3))
    }

    private func drawGrid(_ ctx: GraphicsContext, _ f: PlaneFrame) {
        let half = f.unit / 2
        var major = Path(), minor = Path()
        var i = -40
        while i <= 40 {
            let x = f.origin.x + CGFloat(i) * half, y = f.origin.y + CGFloat(i) * half
            if i != 0 {
                var p = Path()
                if x >= 0 && x <= f.size.width { p.move(to: CGPoint(x: x, y: 0)); p.addLine(to: CGPoint(x: x, y: f.size.height)) }
                if y >= 0 && y <= f.size.height { p.move(to: CGPoint(x: 0, y: y)); p.addLine(to: CGPoint(x: f.size.width, y: y)) }
                if i % 2 == 0 { major.addPath(p) } else { minor.addPath(p) }
            }
            i += 1
        }
        ctx.stroke(minor, with: .color(.white.opacity(0.06)), lineWidth: 0.7)
        ctx.stroke(major, with: .color(.white.opacity(0.15)), lineWidth: 1)
    }

    private func drawAxes(_ ctx: GraphicsContext, _ f: PlaneFrame) {
        let top: CGFloat = f.size.height > 500 ? 116 : 10
        var axes = Path()
        axes.move(to: CGPoint(x: 0, y: f.origin.y)); axes.addLine(to: CGPoint(x: f.size.width - 12, y: f.origin.y))
        axes.move(to: CGPoint(x: f.origin.x, y: f.size.height)); axes.addLine(to: CGPoint(x: f.origin.x, y: top))
        for i in -6...6 where i != 0 {
            let p = f.point(Double(i), 0)
            if p.x > 14 && p.x < f.size.width - 20 { axes.move(to: CGPoint(x: p.x, y: p.y - 5)); axes.addLine(to: CGPoint(x: p.x, y: p.y + 5)) }
            let q = f.point(0, Double(i))
            if q.y > top + 14 && q.y < f.size.height - 20 { axes.move(to: CGPoint(x: q.x - 5, y: q.y)); axes.addLine(to: CGPoint(x: q.x + 5, y: q.y)) }
        }
        ctx.stroke(axes, with: .color(.white), lineWidth: 2)

        var tips = Path()
        let w = f.size.width
        tips.move(to: CGPoint(x: w - 2, y: f.origin.y)); tips.addLine(to: CGPoint(x: w - 15, y: f.origin.y - 6)); tips.addLine(to: CGPoint(x: w - 15, y: f.origin.y + 6)); tips.closeSubpath()
        tips.move(to: CGPoint(x: f.origin.x, y: top - 8)); tips.addLine(to: CGPoint(x: f.origin.x - 6, y: top + 5)); tips.addLine(to: CGPoint(x: f.origin.x + 6, y: top + 5)); tips.closeSubpath()
        ctx.fill(tips, with: .color(.white))
    }
}

struct GaugePlane: View {
    let map: [Double]
    let fadeStart: CGFloat
    var scale: CGFloat = 1
    var shift: CGSize = .zero

    var body: some View {
        GeometryReader { geo in
            let ox = map[1] * scale + shift.width, oy = -map[3] * scale + shift.height, ux = map[0] * scale, uy = map[2] * scale
            let sx = Self.step(ux), sy = Self.step(uy)
            ZStack(alignment: .topLeading) {
                Canvas { ctx, size in
                    var grid = Path(), axes = Path()
                    for k in -60...60 where k != 0 {
                        let x = ox + Double(k) * sx * ux, y = oy - Double(k) * sy * uy
                        if x > 0 && x < size.width { grid.move(to: CGPoint(x: x, y: 0)); grid.addLine(to: CGPoint(x: x, y: size.height)) }
                        if y > 0 && y < size.height { grid.move(to: CGPoint(x: 0, y: y)); grid.addLine(to: CGPoint(x: size.width, y: y)) }
                    }
                    ctx.stroke(grid, with: .color(.white.opacity(0.12)), lineWidth: 1)
                    axes.move(to: CGPoint(x: 0, y: oy)); axes.addLine(to: CGPoint(x: size.width, y: oy))
                    axes.move(to: CGPoint(x: ox, y: 0)); axes.addLine(to: CGPoint(x: ox, y: size.height))
                    ctx.stroke(axes, with: .color(.white.opacity(0.9)), lineWidth: 1.6)
                }
                ForEach([-2, -1, 1, 2], id: \.self) { k in
                    label(Double(k) * sx).position(x: min(max(ox + Double(k) * sx * ux, 16), geo.size.width - 16),
                                                   y: min(max(oy + 14, 10), geo.size.height - 10))
                    label(Double(k) * sy).position(x: min(max(ox - 18, 16), geo.size.width - 16),
                                                   y: min(max(oy - Double(k) * sy * uy, 10), geo.size.height - 10))
                }
            }
            .frame(width: geo.size.width, height: geo.size.height)
            .mask(LinearGradient(stops: [.init(color: .white, location: 0), .init(color: .white, location: fadeStart),
                                         .init(color: .clear, location: fadeStart + 0.08)], startPoint: .top, endPoint: .bottom))
        }
        .accessibilityHidden(true)
    }

    private func label(_ v: Double) -> some View {
        MathText(Self.fmt(v), size: 12, color: .graphiti(hex: 0xFFFFFF, alpha: 0.75))
            .padding(.horizontal, 2).background(Theme.ink.opacity(0.85), in: .rect(cornerRadius: 3))
    }

    static func step(_ unitPoints: Double) -> Double {
        let raw = 64 / max(abs(unitPoints), 1e-9), e = pow(10, floor(log10(raw))), f = raw / e
        return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * e
    }
    static func fmt(_ v: Double) -> String {
        if abs(v - v.rounded()) < 1e-9 { return String(Int(v.rounded())).replacingOccurrences(of: "-", with: "-") }
        var s = String(format: "%.3f", v); while s.hasSuffix("0") { s.removeLast() }
        return s
    }
}
