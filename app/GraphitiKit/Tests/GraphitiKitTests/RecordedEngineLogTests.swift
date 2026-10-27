import XCTest
@testable import GraphitiKit

final class RecordedEngineLogTests: XCTestCase {
    private static let slot = 11848

    private static func writeSlot(_ raw: UnsafeMutableRawPointer, seq: UInt32, step: Int32, pairs: Int32, final: Int32) {
        raw.storeBytes(of: seq, toByteOffset: 0, as: UInt32.self)
        raw.storeBytes(of: step, toByteOffset: 4, as: Int32.self)
        raw.storeBytes(of: pairs, toByteOffset: 8, as: Int32.self)
        raw.storeBytes(of: final, toByteOffset: 12, as: Int32.self)
        raw.storeBytes(of: Double(step) * 10, toByteOffset: 16, as: Double.self)
        raw.storeBytes(of: Int64(96 * Int(step)), toByteOffset: 24, as: Int64.self)
        raw.storeBytes(of: Float(0), toByteOffset: 32, as: Float.self)
        raw.storeBytes(of: Float(1), toByteOffset: 36, as: Float.self)
        let p = 40
        raw.storeBytes(of: UInt8(1), toByteOffset: p, as: UInt8.self)
        raw.storeBytes(of: UInt8(2), toByteOffset: p + 1, as: UInt8.self)
        raw.storeBytes(of: Int16(-1), toByteOffset: p + 26, as: Int16.self)
        let y = p + 172
        raw.storeBytes(of: UInt8(3), toByteOffset: y, as: UInt8.self)
        raw.storeBytes(of: UInt8(1), toByteOffset: y + 1, as: UInt8.self)
        raw.storeBytes(of: UInt8(2), toByteOffset: y + 2, as: UInt8.self)
        raw.storeBytes(of: UInt8(5), toByteOffset: y + 3, as: UInt8.self)
        raw.storeBytes(of: Int16(0), toByteOffset: y + 26, as: Int16.self)
        raw.storeBytes(of: Float(Double(step) * 0.1), toByteOffset: y + 76, as: Float.self)
        let q = p + 344
        raw.storeBytes(of: Float(0.02), toByteOffset: q, as: Float.self)
        raw.storeBytes(of: UInt16(3), toByteOffset: q + 4, as: UInt16.self)
        raw.storeBytes(of: UInt8(1), toByteOffset: q + 6, as: UInt8.self)
        for i in 0..<120 {
            raw.storeBytes(of: Float(Double(i) / 119), toByteOffset: q + 8 + i * 4, as: Float.self)
            raw.storeBytes(of: Float(Double(i) / 119 * Double(step) * 0.1), toByteOffset: q + 8 + 480 + i * 4, as: Float.self)
        }
    }

    @MainActor
    func testRecordedEngineLogDrivesSolveLive() async {
        var log = [UInt8](repeating: 0, count: Self.slot * 4)
        log.withUnsafeMutableBytes { buf in
            for k in 0..<3 {
                Self.writeSlot(buf.baseAddress! + k * Self.slot, seq: UInt32(k), step: Int32(k + 1), pairs: 1, final: 0)
            }
            Self.writeSlot(buf.baseAddress! + 3 * Self.slot, seq: 3, step: -1, pairs: 1, final: 1)
        }
        let bursts: [Burst] = (0..<4).map { k in
            log.withUnsafeBytes { EngineRingDecoder.decode($0.baseAddress! + k * Self.slot) }
        }
        XCTAssertTrue(bursts[3].isFinal)
        XCTAssertFalse(bursts[0].isFinal)
        let session = SolveSession()
        let t0 = Date()
        await session.solveLive(RecordedBurstSource(log: bursts))
        let took = Date().timeIntervalSince(t0)
        XCTAssertEqual(session.liveSteps.count, 3, "one step per non-final champion burst")
        XCTAssertEqual(session.liveSteps.first?.latex, "y = 0.1 \\cdot x")
        XCTAssertGreaterThanOrEqual(took, 0.5, "250 ms minimum on screen between champions")
        XCTAssertEqual(session.phase, .locked)
    }
}
