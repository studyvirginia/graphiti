import Foundation

public struct ConvergeStop {
    public static let defaultsKey = "convergeStop"
    public static let patience = 3
    public static let epsNats = 0.5
    public static let minSteps = 4
    public static let nResiduals = 96.0

    public static var enabled: Bool { UserDefaults.standard.bool(forKey: defaultsKey) }

    private var best = Double.infinity
    private var structure: Int?
    private var quiet = 0
    private var seen = 0
    public private(set) var stoppedAtStep: Int?

    public init() {}

    static func structureID(_ p: Burst.Pair) -> Int {
        var h = Hasher()
        for t in p.xTokens { h.combine(t.op.rawValue) }
        h.combine(0xFF as UInt8)
        for t in p.yTokens { h.combine(t.op.rawValue) }
        return h.finalize()
    }

    public mutating func observe(step: Int, champion p: Burst.Pair) -> Bool {
        seen += 1
        let sid = Self.structureID(p)
        let s = p.score
        var gain = 0.0
        if s.isFinite && s > 0 {
            if best.isFinite && best > 0 { gain = (Self.nResiduals / 2) * log(best / s) } else if !best.isFinite { gain = .infinity }
            best = min(best, s)
        }
        let quietStep = gain < Self.epsNats && structure == sid
        structure = sid
        quiet = quietStep ? quiet + 1 : 0
        if seen >= Self.minSteps && quiet >= Self.patience { stoppedAtStep = step; return true }
        return false
    }
}
