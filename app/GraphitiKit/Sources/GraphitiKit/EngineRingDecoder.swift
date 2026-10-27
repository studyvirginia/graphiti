import Foundation

public enum EngineRingDecoder {
    public struct Layout: Sendable, Equatable {
        public var burst = 45008, pairs = 40, pair = 4996, tree = 1796
        public var treeOp = 1, treeCid = 258, treeC = 772, maxLen = 256
        public var score = 3592, cx = 3596, champion = 3598, fitted = 3599, xs = 3600, ys = 4080, fit = 4560, fitSize = 436
        public init() {}
        public static let standard = Layout()
        public init(engine o: [Int32]) {
            precondition(o.count >= 16)
            burst = Int(o[0]); pairs = Int(o[1]); pair = Int(o[2]); tree = Int(o[3]); treeOp = Int(o[4]); treeCid = Int(o[5]); treeC = Int(o[6])
            score = Int(o[7]); cx = Int(o[8]); champion = Int(o[9]); fitted = Int(o[10]); xs = Int(o[11]); ys = Int(o[12]); fit = Int(o[13]); fitSize = Int(o[15])
            maxLen = treeCid - treeOp
            if (maxLen & 1) == 1 { maxLen -= 1 }
            _ = o[14]
        }
    }
    public static let standard = Layout()

    public static func decode(_ raw: UnsafeRawPointer, layout L: Layout = .standard) -> Burst {
        let seq = Int(raw.load(fromByteOffset: 0, as: UInt32.self))
        let step = Int(raw.load(fromByteOffset: 4, as: Int32.self))
        let nPairs = Int(raw.load(fromByteOffset: 8, as: Int32.self))
        let finalPick = raw.load(fromByteOffset: 12, as: Int32.self) != 0
        let wall = raw.load(fromByteOffset: 16, as: Double.self)
        let cand = Int(raw.load(fromByteOffset: 24, as: Int64.self))
        let tLo = Double(raw.load(fromByteOffset: 32, as: Float.self))
        let tHi = Double(raw.load(fromByteOffset: 36, as: Float.self))
        var pairs: [Burst.Pair] = []
        for p in 0..<min(nPairs, 9) {
            let base = L.pairs + p * L.pair
            let x = decodeTree(raw, base, L)
            let y = decodeTree(raw, base + L.tree, L)
            let score = Double(raw.load(fromByteOffset: base + L.score, as: Float.self))
            let cx = Int(raw.load(fromByteOffset: base + L.cx, as: UInt16.self))
            let champ = raw.load(fromByteOffset: base + L.champion, as: UInt8.self) != 0
            let fittedKind = Int(raw.load(fromByteOffset: base + L.fitted, as: UInt8.self))
            let xs = (0..<120).map { Double(raw.load(fromByteOffset: base + L.xs + $0 * 4, as: Float.self)) }
            let ys = (0..<120).map { Double(raw.load(fromByteOffset: base + L.ys + $0 * 4, as: Float.self)) }
            let fit = fittedKind != 0 ? decodeFit(raw, base + L.fit) : nil
            pairs.append(Burst.Pair(x: xs, y: ys, xTokens: x, yTokens: y, score: score, complexity: cx, isChampion: champ, fit: fit))
        }
        return Burst(seq: seq, step: step, wallMs: wall, candidates: cand, isFinal: finalPick,
                     tLo: tLo, tHi: tHi, pairs: pairs)
    }

    static func decodeTree(_ raw: UnsafeRawPointer, _ base: Int, _ L: Layout) -> [BurstToken] {
        let len = min(Int(raw.load(fromByteOffset: base, as: UInt8.self)), L.maxLen)
        var out: [BurstToken] = []
        for k in 0..<len {
            let op = raw.load(fromByteOffset: base + L.treeOp + k, as: UInt8.self)
            let cid = Int(raw.load(fromByteOffset: base + L.treeCid + k * 2, as: Int16.self))
            let c = Double(raw.load(fromByteOffset: base + L.treeC + k * 4, as: Float.self))
            guard let o = BurstToken.Op(rawValue: op) else { continue }
            out.append(BurstToken(op: o, value: c, constID: o == .constant && cid >= 0 ? cid : nil))
        }
        return out
    }

    static func decodeFit(_ raw: UnsafeRawPointer, _ b: Int) -> Burst.Fit {
        func f(_ off: Int) -> Double { Double(raw.load(fromByteOffset: b + off, as: Float.self)) }
        func u8(_ off: Int) -> Int { Int(raw.load(fromByteOffset: b + off, as: UInt8.self)) }
        let q = [u8(2), u8(3)], nth = [u8(4), u8(5)]
        return Burst.Fit(kind: u8(0), overflow: u8(1) != 0, rootCount: q, thetaCount: nth, isDriver: [u8(6) != 0, u8(7) != 0],
                         warp: (0..<10).map { f(12 + $0 * 4) },
                         theta: (0..<2).map { s in (0..<min(nth[s], 32)).map { f(52 + (s * 32 + $0) * 4) } },
                         coefficients: (0..<2).map { s in (0..<min(q[s], 16)).map { f(308 + (s * 16 + $0) * 4) } })
    }
}
