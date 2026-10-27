import SwiftUI

public struct CurveParams: VectorArithmetic, Sendable {
    public var m: Double, a: Double, w: Double, q: Double

    public init(m: Double, a: Double, w: Double, q: Double) {
        self.m = m; self.a = a; self.w = w; self.q = q
    }
    public func y(_ x: Double) -> Double { m * x + a * sin(w * x) + q * x * x }

    public static var zero: CurveParams { .init(m: 0, a: 0, w: 0, q: 0) }
    public static func + (l: CurveParams, r: CurveParams) -> CurveParams { .init(m: l.m + r.m, a: l.a + r.a, w: l.w + r.w, q: l.q + r.q) }
    public static func - (l: CurveParams, r: CurveParams) -> CurveParams { .init(m: l.m - r.m, a: l.a - r.a, w: l.w - r.w, q: l.q - r.q) }
    public mutating func scale(by s: Double) { m *= s; a *= s; w *= s; q *= s }
    public var magnitudeSquared: Double { m * m + a * a + w * w + q * q }
}

public struct CurveShape: Shape {
    public var params: CurveParams
    var wobble: Double = 0
    var layout: PlaneLayout
    var samples = 360

    public init(params: CurveParams, wobble: Double = 0, layout: PlaneLayout = .standard) {
        self.params = params; self.wobble = wobble; self.layout = layout
    }

    public var animatableData: CurveParams {
        get { params }
        set { params = newValue }
    }

    public func path(in rect: CGRect) -> Path {
        let f = PlaneFrame(size: rect.size, layout: layout)
        let r = f.xRange
        var p = Path()
        for i in 0...samples {
            let x = r.lowerBound + (r.upperBound - r.lowerBound) * Double(i) / Double(samples)
            var y = params.y(x)
            if wobble > 0 { y += wobble * (sin(x * 9.3) * 0.6 + sin(x * 3.1 + 1) * 0.8) }
            let pt = f.point(x, y)
            if i == 0 { p.move(to: pt) } else { p.addLine(to: pt) }
        }
        return p
    }
}
