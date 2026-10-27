import XCTest
@testable import GraphitiKit

final class BurstStreamTests: XCTestCase {
    private static func constant(_ id: Int, _ v: Double) -> BurstToken {
        BurstToken(op: .constant, value: v, constID: id)
    }

    private static func pair(xt: [BurstToken], yt: [BurstToken], champion: Bool,
                             score: Double, complexity: Int) -> Burst.Pair {
        let xs = (0..<120).map { Double($0) / 119 }
        let ys = xs.map { $0 * 1.4 }
        return Burst.Pair(x: xs, y: ys, xTokens: xt, yTokens: yt, score: score,
                          complexity: complexity, isChampion: champion)
    }

    private static func recordedLog() -> [Burst] {
        let t = BurstToken(op: .t, value: 0, constID: nil)
        let x1 = [BurstToken(op: .constant, value: 1, constID: 0), t, BurstToken(op: .mul, value: 0, constID: nil)]
        let champ = pair(xt: [t], yt: [constant(1, 0.5), t, BurstToken(op: .mul, value: 0, constID: nil)],
                         champion: true, score: 0.02, complexity: 3)
        let rival = pair(xt: x1, yt: [constant(2, 2), BurstToken(op: .sin, value: 0, constID: nil)],
                         champion: false, score: 0.4, complexity: 5)
        return [
            Burst(seq: 0, step: 1, wallMs: 3, candidates: 96, isFinal: false, tLo: 0, tHi: 1, pairs: [champ]),
            Burst(seq: 1, step: 2, wallMs: 9, candidates: 192, isFinal: false, tLo: 0, tHi: 1, pairs: [champ, rival]),
            Burst(seq: 2, step: 3, wallMs: 15, candidates: 288, isFinal: true, tLo: 0, tHi: 1, pairs: [champ]),
        ]
    }

    @MainActor
    func testConsumeBuildsStepsFromRecordedLog() async {
        let session = SolveSession()
        await session.consume(RecordedBurstSource(log: Self.recordedLog()))
        XCTAssertEqual(session.liveSteps.count, 2, "one step per non-final burst with a champion")
        XCTAssertEqual(session.steps.count, 2)
        let first = session.liveSteps[0]
        XCTAssertEqual(first.samples.count, 120)
        XCTAssertEqual(first.latex, "y = 0.5 \\cdot x")
        XCTAssertEqual(session.phase, .locked)
    }

    @MainActor
    func testFinalBurstEndsStreamWithoutStep() async {
        let log = Self.recordedLog()
        let session = SolveSession()
        await session.consume(RecordedBurstSource(log: log))
        XCTAssertFalse(session.liveSteps.contains { $0.samples.isEmpty })
        XCTAssertTrue(log.last?.isFinal == true)
    }

    func testNonSingleValuedChampionShowsPairLatex() {
        let rival = Self.recordedLog()[1].pairs[1]
        XCTAssertFalse(BurstDisplay.isSingleValued(rival))
        XCTAssertTrue(BurstDisplay.latex(rival).hasPrefix("x(t) = "))
        XCTAssertTrue(BurstDisplay.latex(rival).contains(",\\quad y(t) = "))
    }

    func testTokensKeyConstantsByStableID() {
        let champ = Self.recordedLog()[0].pairs[0]
        let ids = BurstDisplay.tokens(champ).map(\.id)
        XCTAssertTrue(ids.contains("c1"))
    }
}
