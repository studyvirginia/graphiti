import XCTest
@testable import GraphitiKit

final class EngineRingDecoderTests: XCTestCase {
    func testDecodesCLayoutSlot() {
        var buf = [UInt8](repeating: 0, count: 45008)
        buf.withUnsafeMutableBytes { raw in
            raw.storeBytes(of: UInt32(7), toByteOffset: 0, as: UInt32.self)
            raw.storeBytes(of: Int32(3), toByteOffset: 4, as: Int32.self)
            raw.storeBytes(of: Int32(1), toByteOffset: 8, as: Int32.self)
            raw.storeBytes(of: Int32(0), toByteOffset: 12, as: Int32.self)
            raw.storeBytes(of: 12.5, toByteOffset: 16, as: Double.self)
            raw.storeBytes(of: Int64(96), toByteOffset: 24, as: Int64.self)
            raw.storeBytes(of: Float(0), toByteOffset: 32, as: Float.self)
            raw.storeBytes(of: Float(1), toByteOffset: 36, as: Float.self)
            let p = 40
            raw.storeBytes(of: UInt8(1), toByteOffset: p, as: UInt8.self)
            raw.storeBytes(of: UInt8(2), toByteOffset: p + 1, as: UInt8.self)
            raw.storeBytes(of: Int16(-1), toByteOffset: p + 258, as: Int16.self)
            let y = p + 1796
            raw.storeBytes(of: UInt8(3), toByteOffset: y, as: UInt8.self)
            raw.storeBytes(of: UInt8(1), toByteOffset: y + 1, as: UInt8.self)
            raw.storeBytes(of: UInt8(2), toByteOffset: y + 2, as: UInt8.self)
            raw.storeBytes(of: UInt8(5), toByteOffset: y + 3, as: UInt8.self)
            raw.storeBytes(of: Int16(0), toByteOffset: y + 258, as: Int16.self)
            raw.storeBytes(of: Float(0.5), toByteOffset: y + 772, as: Float.self)
            let q = p + 3592
            raw.storeBytes(of: Float(0.02), toByteOffset: q, as: Float.self)
            raw.storeBytes(of: UInt16(3), toByteOffset: q + 4, as: UInt16.self)
            raw.storeBytes(of: UInt8(1), toByteOffset: q + 6, as: UInt8.self)
            for i in 0..<120 {
                raw.storeBytes(of: Float(Double(i) / 119), toByteOffset: q + 8 + i * 4, as: Float.self)
                raw.storeBytes(of: Float(Double(i) / 119 * 0.5), toByteOffset: q + 8 + 480 + i * 4, as: Float.self)
            }
        }
        let b = buf.withUnsafeBytes { EngineRingDecoder.decode($0.baseAddress!) }
        XCTAssertEqual(b.seq, 7)
        XCTAssertEqual(b.step, 3)
        XCTAssertFalse(b.isFinal)
        XCTAssertEqual(b.pairs.count, 1)
        XCTAssertEqual(b.champion?.xTokens.map(\.op), [.t])
        XCTAssertEqual(b.champion?.yTokens.map(\.op), [.constant, .t, .mul])
        XCTAssertEqual(b.champion?.yTokens.first?.constID, 0)
        XCTAssertEqual(b.champion?.yTokens.first?.value ?? 0, 0.5, accuracy: 1e-6)
        XCTAssertEqual(b.champion?.complexity, 3)
        XCTAssertEqual(b.champion?.y.count, 120)
        XCTAssertEqual(b.champion?.y.last ?? 0, 0.5, accuracy: 1e-6)
    }

    func testBudgetClampsAndDefaults() {
        XCTAssertEqual(SearchBudget().seconds, 8)
        XCTAssertEqual(SearchBudget(seconds: 0).seconds, 2)
        XCTAssertEqual(SearchBudget(seconds: 999).seconds, 60)
    }
}
