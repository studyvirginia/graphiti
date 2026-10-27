import Foundation

public enum CurveEditor {
    public static var weights = CurveParams(m: 1.0, a: 1.4, w: 0.25, q: 0.6)

    public static func drag(_ p: CurveParams, at x: Double, by dy: Double) -> CurveParams {
        let h = 1e-4
        let base = p.y(x)
        let j = CurveParams(
            m: (CurveParams(m: p.m + h, a: p.a, w: p.w, q: p.q).y(x) - base) / h,
            a: (CurveParams(m: p.m, a: p.a + h, w: p.w, q: p.q).y(x) - base) / h,
            w: (CurveParams(m: p.m, a: p.a, w: p.w + h, q: p.q).y(x) - base) / h,
            q: (CurveParams(m: p.m, a: p.a, w: p.w, q: p.q + h).y(x) - base) / h)
        let wj = CurveParams(m: weights.m * j.m * active(p.m), a: weights.a * j.a * active(p.a),
                             w: weights.w * j.w * active(p.w), q: weights.q * j.q * active(p.q))
        let denom = wj.m * j.m + wj.a * j.a + wj.w * j.w + wj.q * j.q
        guard abs(denom) > 1e-9 else { return p }
        let s = dy / denom
        return CurveParams(m: p.m + wj.m * s, a: p.a + wj.a * s, w: p.w + wj.w * s, q: p.q + wj.q * s)
    }

    public static func drag(_ p: CurveParams, at x: Double, to targetY: Double, iterations: Int = 8) -> CurveParams {
        var q = p
        for _ in 0..<iterations {
            let miss = targetY - q.y(x)
            if abs(miss) < 1e-6 { break }
            q = drag(q, at: x, by: miss)
        }
        return q
    }

    private static func active(_ v: Double) -> Double { abs(v) > 1e-9 ? 1 : 0 }
}
