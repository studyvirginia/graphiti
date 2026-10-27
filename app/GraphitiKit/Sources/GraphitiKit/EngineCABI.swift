import Foundation
#if os(iOS)
import UIKit
#endif

@_silgen_name("sketch_ui_ring_create") func sketch_ui_ring_create() -> UnsafeMutableRawPointer?
@_silgen_name("sketch_ui_ring_destroy") func sketch_ui_ring_destroy(_ ring: UnsafeMutableRawPointer?)
@_silgen_name("sketch_ui_pop") func sketch_ui_pop(_ ring: UnsafeMutableRawPointer?, _ out: UnsafeMutableRawPointer?) -> Int32
@_silgen_name("sketch_ui_stop") func sketch_ui_stop()
@_silgen_name("sketch_ui_map") func sketch_ui_map(_ out4: UnsafeMutablePointer<Double>)
@_silgen_name("sketch_ui_burst_bytes") func sketch_ui_burst_bytes() -> UInt
@_silgen_name("sketch_ui_layout") func sketch_ui_layout(_ out16: UnsafeMutablePointer<Int32>)
@_silgen_name("graphiti_stroke_log_path") func graphiti_stroke_log_path() -> UnsafePointer<CChar>?
@_silgen_name("graphiti_solve") func graphiti_solve(_ xy: UnsafePointer<Float>?, _ n: Int32, _ budget: Double,
                                                    _ ring: UnsafeMutableRawPointer?) -> Int32

public enum EngineSolve {
    public static func run(points: [CGPoint], budget: Double, session: SolveSession) async -> Int32 {
        var xy: [Float] = []
        xy.reserveCapacity(points.count * 2)
        for p in points { xy.append(Float(p.x)); xy.append(Float(-p.y)) }
        guard let ring = sketch_ui_ring_create() else { return -1 }
        defer { sketch_ui_ring_destroy(ring) }
        var lay = [Int32](repeating: 0, count: 16); lay.withUnsafeMutableBufferPointer { sketch_ui_layout($0.baseAddress!) }
        let bridge = EngineBridge(handle: .init(raw: ring), slotBytes: Int(sketch_ui_burst_bytes()), layout: .init(engine: lay)) { r, out in
            sketch_ui_pop(r, out) == 1
        }
        let rcBox = RCBox()
        let t0 = Date()
        Diag.log("solve", "start pts=\(points.count) budget=\(budget)s mem=\(Diag.availableMB())MB")
        let worker = Thread {
            let rc = xy.withUnsafeBufferPointer { graphiti_solve($0.baseAddress, Int32(points.count), budget, ring) }
            rcBox.set(rc)
        }
        worker.stackSize = 16 << 20
        worker.start()
        await session.consume(bridge.source()) {
            var m = [Double](repeating: 0, count: 4); m.withUnsafeMutableBufferPointer { sketch_ui_map($0.baseAddress!) }; return m
        }
        let pick = await MainActor.run { session.finalPick }
        while !rcBox.done { try? await Task.sleep(for: .milliseconds(50)) }
        let rc = rcBox.value
        print("[graphiti] solve rc=\(rc) pts=\(points.count) budget=\(budget)"); fflush(stdout)
        let ms = Int(Date().timeIntervalSince(t0) * 1000)
        if rc != 0 { Diag.log("solve", "engine rc=\(rc) after \(ms) ms (nonzero: see engine log)", level: .error) }
        Diag.log("solve", "done rc=\(rc) ms=\(ms) final=\(pick.map { "score=\($0.score) latex=\($0.latex)" } ?? "none") mem=\(Diag.availableMB())MB")
        appendFinalPick(pick)
        return rc
    }

    private static func appendFinalPick(_ pick: (latex: String, score: Double)?) {
        guard let p = graphiti_stroke_log_path() else { return }
        let path = String(cString: p)
        let tex = String((pick?.latex ?? "").replacingOccurrences(of: "\n", with: " ").prefix(300))
        let line = pick.map { "# FINAL score=\(String(format: "%.6g", $0.score)) latex=\(tex)\n" } ?? "# FINAL none\n"
        guard let h = FileHandle(forWritingAtPath: path) else { return }
        h.seekToEndOfFile(); h.write(Data(line.utf8)); h.closeFile()
    }
}

final class RCBox: @unchecked Sendable {
    private let lock = NSLock()
    private var _value: Int32 = 0
    private var _done = false
    func set(_ v: Int32) { lock.lock(); _value = v; _done = true; lock.unlock() }
    var value: Int32 { lock.lock(); defer { lock.unlock() }; return _value }
    var done: Bool { lock.lock(); defer { lock.unlock() }; return _done }
}

@_silgen_name("graphiti_diag_init") func graphiti_diag_init(_ logdir: UnsafePointer<CChar>?, _ session: UnsafePointer<CChar>?, _ verbose: Int32)
@_silgen_name("graphiti_prepare") func graphiti_prepare() -> Int32
@_silgen_name("graphiti_selftest") func graphiti_selftest() -> Int32

public enum EngineSelfTest {
    nonisolated(unsafe) private static var started = false
    @MainActor public static func startIfRequested() {
        guard !started, CommandLine.arguments.contains("-selftest") else { return }
        started = true
        let t = Thread {
            _ = graphiti_prepare()
            let rc = graphiti_selftest()
            print("[graphiti] selftest done rc=\(rc)"); fflush(stdout)
        }
        t.stackSize = 16 << 20
        t.qualityOfService = .userInitiated
        t.start()
    }
}

public enum EnginePrepare {
    public static func stopRunning() { sketch_ui_stop() }
    nonisolated(unsafe) private static var started = false
    @MainActor public static func start() {
        guard !started else { return }
        #if os(iOS)
        guard UIApplication.shared.applicationState == .active else { return }
        #endif
        started = true
        let t = Thread {
            let rc = graphiti_prepare()
            Diag.log("prepare", "rc=\(rc)", level: rc == 0 ? .default : .error)
            if rc != 0 { DispatchQueue.main.async { started = false } }
        }
        t.stackSize = 16 << 20
        t.qualityOfService = .userInitiated
        t.start()
    }
}
