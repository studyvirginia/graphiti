import Foundation

public final class EngineBridge: @unchecked Sendable {
    public struct Handle: @unchecked Sendable { public let raw: UnsafeMutableRawPointer }
    public let handle: Handle
    public let slotBytes: Int
    public let layout: EngineRingDecoder.Layout

    public typealias Pop = (_ ring: UnsafeMutableRawPointer, _ out: UnsafeMutableRawPointer) -> Bool

    private let pop: Pop

    public init(handle: Handle, slotBytes: Int, layout: EngineRingDecoder.Layout = .standard, pop: @escaping Pop) {
        self.handle = handle; self.slotBytes = slotBytes; self.layout = layout; self.pop = pop
    }

    public func source() -> LiveBurstSource {
        let h = handle, size = slotBytes, pop = self.pop, layout = self.layout
        return LiveBurstSource(pop: {
            let buf = UnsafeMutableRawPointer.allocate(byteCount: size, alignment: 8)
            defer { buf.deallocate() }
            guard pop(h.raw, buf) else { return nil }
            return EngineRingDecoder.decode(UnsafeRawPointer(buf), layout: layout)
        })
    }
}
