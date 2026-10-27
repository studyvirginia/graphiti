import Foundation
import CoreGraphics


public struct NativeNoiseField: Sendable {
    public let s2: [Float]
}

public struct ParametricCandidate: Identifiable, Sendable, Equatable {
    public let id = UUID()
    public let nodes: Int
    public let shapeErr: Double
    public let r2: Double
    public let dlBits: Double
    public let exprX: String
    public let exprY: String
    public let displayExprX: String
    public let displayExprY: String
    public let curvePoints: [CGPoint]
    public var rank: Int = 0

    public static func == (lhs: ParametricCandidate, rhs: ParametricCandidate) -> Bool { lhs.id == rhs.id }
}

public struct EvolutionBurst: Sendable {
    public let generation: Int
    public let elapsedMs: Int
    public let champion: ParametricCandidate
    public let candidates: [ParametricCandidate]
    public let isFinal: Bool
}


public final class SketchNativeEngine: @unchecked Sendable {
    public static let shared = SketchNativeEngine()
    public var isAvailable: Bool { true }

    public func solveStrokeStream(
        points: [CGPoint],
        population: Int = 500,
        generations: Int = 60,
        maxLength: Int = 28,
        noiseField: NativeNoiseField? = nil
    ) -> AsyncStream<EvolutionBurst> {
        AsyncStream<EvolutionBurst>(bufferingPolicy: .unbounded) { continuation in
            guard points.count >= 4 else { continuation.finish(); return }
            Task.detached {
                let steps = 8
                let start = Date()
                for g in 1...steps {
                    if Task.isCancelled { break }
                    let t = Double(g) / Double(steps)
                    let curve = Self.smooth(points, passes: g * 3)
                    let err = 0.30 * (1 - t) + 0.01
                    let champ = Self.candidate(curve, err: err, nodes: 14 - g, g: g)
                    let decoys = (0..<4).map { i in
                        Self.candidate(Self.smooth(points, passes: max(1, g * 3 - 2 * (i + 1))),
                                       err: err + 0.05 * Double(i + 1), nodes: 14 - g + 3 * (i + 1), g: g)
                    }
                    continuation.yield(EvolutionBurst(
                        generation: g * (generations / steps),
                        elapsedMs: Int(Date().timeIntervalSince(start) * 1000),
                        champion: champ, candidates: [champ] + decoys, isFinal: g == steps))
                    try? await Task.sleep(nanoseconds: 350_000_000)
                }
                continuation.finish()
            }
        }
    }

    private static func smooth(_ pts: [CGPoint], passes: Int) -> [CGPoint] {
        var p = pts
        for _ in 0..<max(passes, 0) where p.count > 2 {
            var q = p
            for i in 1..<(p.count - 1) {
                q[i] = CGPoint(x: (p[i-1].x + 2 * p[i].x + p[i+1].x) / 4,
                               y: (p[i-1].y + 2 * p[i].y + p[i+1].y) / 4)
            }
            p = q
        }
        return p
    }

    private static func candidate(_ curve: [CGPoint], err: Double, nodes: Int, g: Int) -> ParametricCandidate {
        ParametricCandidate(
            nodes: max(nodes, 3), shapeErr: err, r2: min(0.999, 1 - err * 0.4), dlBits: 100 - Double(g),
            exprX: "x(t) = demo", exprY: "y(t) = demo",
            displayExprX: "x(t) = a·cos(t) + b·cos(3t)", displayExprY: "y(t) = c·sin(t) + d·sin(2t)",
            curvePoints: curve)
    }
}


public final class ContinuousNoisePredictor: @unchecked Sendable {
    public init?(bundle: Bundle) { }
    public func makeNativeNoiseField(points: [CGPoint]) -> NativeNoiseField? { nil }
}
