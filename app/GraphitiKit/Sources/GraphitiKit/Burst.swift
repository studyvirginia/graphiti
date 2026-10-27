import Foundation

public struct Burst: Sendable, Equatable {
    public struct Pair: Sendable, Equatable {
        public var x: [Double]
        public var y: [Double]
        public var xTokens: [BurstToken]
        public var yTokens: [BurstToken]
        public var score: Double
        public var complexity: Int
        public var isChampion: Bool
        public var fit: Fit? = nil
    }

    public struct Fit: Sendable, Equatable {
        public var kind: Int
        public var overflow: Bool
        public var rootCount: [Int]
        public var thetaCount: [Int]
        public var isDriver: [Bool]
        public var warp: [Double]
        public var theta: [[Double]]
        public var coefficients: [[Double]]
    }

    public var seq: Int
    public var step: Int
    public var wallMs: Double
    public var candidates: Int
    public var isFinal: Bool
    public var tLo: Double
    public var tHi: Double
    public var pairs: [Pair]

    public var champion: Pair? { pairs.first { $0.isChampion } }

    public init(seq: Int, step: Int, wallMs: Double, candidates: Int, isFinal: Bool,
                tLo: Double, tHi: Double, pairs: [Pair]) {
        self.seq = seq; self.step = step; self.wallMs = wallMs; self.candidates = candidates
        self.isFinal = isFinal; self.tLo = tLo; self.tHi = tHi; self.pairs = pairs
    }
}

public struct BurstToken: Sendable, Equatable {
    public enum Op: UInt8, Sendable { case constant = 1, t = 2, add = 3, sub = 4, mul = 5, div = 6,
        sin = 7, cos = 8, exp = 9, log = 10, sq = 11, cube = 12, sqrt = 13, cbrt = 14, abs = 15 }
    public var op: Op
    public var value: Double
    public var constID: Int?
}

private indirect enum TexExpr {
    case num(Double)
    case variable(String)
    case add(TexExpr, TexExpr)
    case sub(TexExpr, TexExpr)
    case mul(TexExpr, TexExpr)
    case div(TexExpr, TexExpr)
    case call(String, TexExpr)
    case absCall(String, TexExpr)
    case power(TexExpr, Int)
    case exp(TexExpr)
    case sqrt(TexExpr, Bool)
    case abs(TexExpr)
}

private enum TexKind { case sum, product, atom }

public extension BurstToken {
    static func latex(_ tokens: [BurstToken], variable: String = "x") -> String {
        var stack: [TexExpr] = []
        for t in tokens {
            switch t.op {
            case .constant: stack.append(.num(t.value))
            case .t: stack.append(.variable(variable))
            case .add, .sub, .mul, .div:
                guard stack.count >= 2 else { return "?" }
                let r = stack.removeLast(), l = stack.removeLast()
                switch t.op {
                case .add: stack.append(.add(l, r))
                case .sub: stack.append(.sub(l, r))
                case .mul: stack.append(.mul(l, r))
                default: stack.append(.div(l, r))
                }
            case .sin, .cos, .exp, .log, .sq, .cube, .sqrt, .cbrt, .abs:
                guard let a = stack.popLast() else { return "?" }
                switch t.op {
                case .sin: stack.append(.call("\\sin", a))
                case .cos: stack.append(.call("\\cos", a))
                case .log: stack.append(.absCall("\\ln", a))
                case .exp: stack.append(.exp(a))
                case .sq: stack.append(.power(a, 2))
                case .cube: stack.append(.power(a, 3))
                case .sqrt: stack.append(.sqrt(a, false))
                case .abs: stack.append(.abs(a))
                default: stack.append(.sqrt(a, true))
                }
            }
        }
        guard stack.count == 1 else { return "?" }
        return texRender(stack[0]).0
    }

    private static func texSplit(_ e: TexExpr) -> (neg: Bool, mag: TexExpr) {
        switch e {
        case .num(let v) where v < 0: return (true, .num(-v))
        case .mul(let a, let b):
            let (na, ma) = texSplit(a), (nb, mb) = texSplit(b)
            return (na != nb, .mul(ma, mb))
        case .div(let a, let b):
            let (na, ma) = texSplit(a), (nb, mb) = texSplit(b)
            return (na != nb, .div(ma, mb))
        default: return (false, e)
        }
    }

    private static func texIsNum(_ e: TexExpr, _ value: Double) -> Bool {
        if case .num(let v) = e { return v == value }
        return false
    }

    private static func texIsNumber(_ e: TexExpr) -> Bool {
        if case .num = e { return true }
        return false
    }

    private static func texRender(_ e: TexExpr) -> (String, TexKind) {
        let (neg, mag) = texSplit(e)
        let (body, kind) = texBody(mag)
        guard neg else { return (body, kind) }
        return ("-" + (kind == .sum ? "(\(body))" : body), .product)
    }

    private static func texBody(_ e: TexExpr) -> (String, TexKind) {
        switch e {
        case .num(let v): return (formatNumber(v), .atom)
        case .variable(let name): return (name, .atom)
        case .add(let a, let b):
            if texIsNum(a, 0) { return texRender(b) }
            if texIsNum(b, 0) { return texRender(a) }
            let left = texRender(a).0
            let (nb, mb) = texSplit(b)
            if nb {
                let (rm, km) = texBody(mb)
                return (left + " - " + (km == .sum ? "(\(rm))" : rm), .sum)
            }
            return (left + " + " + texRender(b).0, .sum)
        case .sub(let a, let b):
            if texIsNum(a, 0) { return texRender(.mul(.num(-1), b)) }
            if texIsNum(b, 0) { return texRender(a) }
            let left = texRender(a).0
            let (nb, mb) = texSplit(b)
            if nb { return (left + " + " + texBody(mb).0, .sum) }
            let (rb, kb) = texRender(b)
            return (left + " - " + (kb == .sum ? "(\(rb))" : rb), .sum)
        case .mul(let a, let b):
            if texIsNum(a, 1) { return texBody(b) }
            if texIsNum(b, 1) { return texBody(a) }
            if case .num = b, !texIsNumber(a) { return texBody(.mul(b, a)) }
            var (ra, ka) = texBody(a)
            var (rb, kb) = texBody(b)
            if ka == .sum { ra = "(\(ra))" }
            if kb == .sum { rb = "(\(rb))" }
            let leftNumber: Bool = { if case .num = a { return true } else { return false } }()
            let rightDigit = rb.first.map { $0.isNumber || $0 == "." } ?? false
            let implicit = (leftNumber && !rightDigit) || rb.hasPrefix("\\") || rb.hasPrefix("(")
            return (ra + (implicit ? "" : " \\cdot ") + rb, .product)
        case .div(let a, let b):
            return ("\\frac{\(texRender(a).0)}{\(texRender(b).0)}", .atom)
        case .call(let name, let a):
            return ("\(name)(\(texRender(a).0))", .atom)
        case .absCall(let name, let a):
            return ("\(name)\\left|\(texRender(a).0)\\right|", .atom)
        case .exp(let a):
            return ("e^{\(texRender(a).0)}", .product)
        case .sqrt(let a, let cube):
            return (cube ? "\\sqrt[3]{\(texRender(a).0)}" : "\\sqrt{\\left|\(texRender(a).0)\\right|}", .atom)
        case .abs(let a):
            return ("\\left|\(texRender(a).0)\\right|", .atom)
        case .power(let a, let n):
            let (ra, ka) = texRender(a)
            return ((ka == .atom ? ra : "(\(ra))") + "^{\(n)}", .product)
        }
    }

    static func formatNumber(_ v: Double) -> String {
        let r = (v * 10).rounded() / 10
        if r == 0, v != 0 { return String(format: "%.1g", v) }
        if abs(r - r.rounded()) < 1e-9 { return String(Int(r.rounded())) }
        return String(format: "%.1f", r)
    }
}

public enum BurstDisplay {
    public static func isSingleValued(_ b: Burst.Pair) -> Bool {
        b.xTokens.count == 1 && b.xTokens.first?.op == .t
    }

    public static func isSingleValuedY(_ b: Burst.Pair) -> Bool {
        b.yTokens.count == 1 && b.yTokens.first?.op == .t && !isSingleValued(b)
    }

    public static func tokens(_ b: Burst.Pair) -> [EquationToken] {
        if isSingleValuedY(b) {
            return [EquationToken("x", "x"), EquationToken("eq", "=", .op), EquationToken("xs", BurstToken.latex(b.xTokens, variable: "y"))]
        }
        if isSingleValued(b) {
            return [EquationToken("y", "y"), EquationToken("eq", "=", .op), EquationToken("ys", BurstToken.latex(b.yTokens, variable: "x"))]
        }
        return [EquationToken("xl", "x(t)"), EquationToken("eqx", "=", .op), EquationToken("xs", BurstToken.latex(b.xTokens, variable: "t")),
                EquationToken("sep", ",\\;", .op),
                EquationToken("yl", "y(t)"), EquationToken("eqy", "=", .op), EquationToken("ys", BurstToken.latex(b.yTokens, variable: "t"))]
    }

    public static func latex(_ b: Burst.Pair) -> String {
        if isSingleValuedY(b) { return "x = " + BurstToken.latex(b.xTokens, variable: "y") }
        if isSingleValued(b) { return "y = " + BurstToken.latex(b.yTokens, variable: "x") }
        return "x(t) = " + BurstToken.latex(b.xTokens, variable: "t") + ",\\quad y(t) = " + BurstToken.latex(b.yTokens, variable: "t")
    }

    public static func terms(_ b: Burst.Pair) -> Int {
        let trees = isSingleValued(b) ? [b.yTokens] : isSingleValuedY(b) ? [b.xTokens] : [b.xTokens, b.yTokens]
        return trees.reduce(0) { n, t in n + 1 + t.filter { $0.op == .add || $0.op == .sub }.count }
    }
}
