import Foundation

public enum ExactSnap {
    public struct Match: Equatable, Sendable {
        public let latex: String
        public let value: Double
    }

    public static var tolerance = 0.002

    public static func snap(_ v: Double) -> Match? {
        guard v.isFinite, abs(v) > 1e-9 else { return nil }
        let sign = v < 0 ? "-" : ""
        let x = abs(v)
        var best: (Match, Double)?
        func consider(_ latex: String, _ value: Double, complexity: Double) {
            let err = abs(x - value) / max(value, 1e-9)
            guard err <= tolerance else { return }
            let score = complexity + err * 1000
            if best == nil || score < best!.1 { best = (Match(latex: sign + latex, value: v < 0 ? -value : value), score) }
        }
        for d in [1, 2, 3, 4, 6] {
            let n = (x * Double(d)).rounded()
            guard n >= 1, n <= 60 else { continue }
            if gcd(Int(n), d) != 1 { continue }
            let latex = d == 1 ? "\(Int(n))" : "\\frac{\(Int(n))}{\(d)}"
            consider(latex, n / Double(d), complexity: Double(d) + (d == 1 ? 0 : 1))
        }
        let named: [(String, Double)] = [("\\pi", .pi), ("\\sqrt{2}", 2.0.squareRoot()), ("\\sqrt{3}", 3.0.squareRoot()), ("e", M_E)]
        for (sym, c) in named {
            for d in 1...6 {
                for n in 1...6 where gcd(n, d) == 1 {
                    let latex: String
                    switch (n, d) {
                    case (1, 1): latex = sym
                    case (_, 1): latex = "\(n)\(sym)"
                    case (1, _): latex = "\\frac{\(sym)}{\(d)}"
                    default: latex = "\\frac{\(n)\(sym)}{\(d)}"
                    }
                    consider(latex, Double(n) * c / Double(d), complexity: 2 + Double(n) + Double(d))
                }
            }
        }
        if let b = best?.0, b.latex.contains("\\") || abs(b.value - v) > 1e-9 { return b }
        return nil
    }

    private static func gcd(_ a: Int, _ b: Int) -> Int { b == 0 ? a : gcd(b, a % b) }
}
