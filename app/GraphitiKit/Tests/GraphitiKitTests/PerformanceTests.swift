import XCTest
import SwiftUI
@testable import GraphitiKit

final class PerformanceTests: XCTestCase {
    func testCurvePathsPerFrame() {
        let rect = CGRect(x: 0, y: 0, width: 393, height: 852)
        let p = CurveParams(m: 0.5, a: 1, w: 2, q: 0.1)
        let n = 1000
        let start = Date()
        for i in 0..<n {
            let t = Double(i) / Double(n)
            for _ in 0..<10 { _ = CurveShape(params: p + CurveParams(m: t, a: t, w: t, q: 0), layout: .solve).path(in: rect) }
        }
        let ms = Date().timeIntervalSince(start) * 1000 / Double(n)
        print("curve paths per frame: \(String(format: "%.3f", ms)) ms")
    }

    func testTypesetting() {
        let cache = MathRenderCache.shared
        let n = 200
        var start = Date()
        for i in 0..<n { _ = cache.render("\\sin(\(i)x)", size: 34, color: .white) }
        let miss = Date().timeIntervalSince(start) * 1000 / Double(n)
        start = Date()
        for i in 0..<n { _ = cache.render("\\sin(\(i)x)", size: 34, color: .white) }
        let hit = Date().timeIntervalSince(start) * 1000 / Double(n)
        print("typeset miss: \(String(format: "%.3f", miss)) ms, hit: \(String(format: "%.4f", hit)) ms")
    }
}
