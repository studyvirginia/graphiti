import Foundation

public struct CardPayload: Equatable, Sendable {
    public var params: CurveParams
    public var author: String
    public var version: Int = 1

    public init(params: CurveParams, author: String) {
        self.params = params; self.author = author
    }

    public var url: URL {
        var c = URLComponents()
        c.scheme = "graphiti"; c.host = "card"
        c.queryItems = [
            .init(name: "v", value: "\(version)"),
            .init(name: "m", value: Self.fmt(params.m)), .init(name: "a", value: Self.fmt(params.a)),
            .init(name: "w", value: Self.fmt(params.w)), .init(name: "q", value: Self.fmt(params.q)),
            .init(name: "by", value: author),
        ]
        return c.url!
    }

    public init?(url: URL) {
        guard let c = URLComponents(url: url, resolvingAgainstBaseURL: false), c.host == "card" else { return nil }
        let q = Dictionary(uniqueKeysWithValues: (c.queryItems ?? []).map { ($0.name, $0.value ?? "") })
        guard let m = Double(q["m"] ?? ""), let a = Double(q["a"] ?? ""),
              let w = Double(q["w"] ?? ""), let qq = Double(q["q"] ?? "") else { return nil }
        self.params = CurveParams(m: m, a: a, w: w, q: qq)
        self.author = q["by"] ?? "someone"
        self.version = Int(q["v"] ?? "1") ?? 1
    }

    public var latex: String { EquationFormatter.latex(params) }

    private static func fmt(_ v: Double) -> String { String(format: "%.4g", v) }
}

public enum EquationFormatter {
    public static func latex(_ p: CurveParams) -> String {
        var terms: [String] = []
        if abs(p.m) > 1e-6 { terms.append(coef(p.m) + "x") }
        if abs(p.a) > 1e-6 && abs(p.w) > 1e-6 { terms.append(coef(p.a) + "\\sin(" + coef(p.w) + "x)") }
        if abs(p.q) > 1e-6 { terms.append(coef(p.q) + "x^{2}") }
        guard !terms.isEmpty else { return "y = 0" }
        var s = "y = " + terms[0]
        for t in terms.dropFirst() {
            s += t.hasPrefix("-") ? " - " + t.dropFirst() : " + " + t
        }
        return s
    }

    static func coef(_ v: Double) -> String {
        if abs(v - 1) < 1e-6 { return "" }
        if abs(v + 1) < 1e-6 { return "-" }
        if let s = ExactSnap.snap(v) { return s.latex.contains("\\") ? "{\(s.latex)}" : s.latex }
        return String(format: "%.2g", v)
    }
}
