import Foundation

public protocol BurstSource: Sendable {
    func bursts() -> AsyncStream<Burst>
}

public final class LiveBurstSource: BurstSource, @unchecked Sendable {
    public typealias Pop = @Sendable () -> Burst?
    private let pop: Pop
    private let interval: Duration

    public init(pop: @escaping Pop, interval: Duration = .milliseconds(16)) {
        self.pop = pop; self.interval = interval
    }

    public func bursts() -> AsyncStream<Burst> {
        AsyncStream(bufferingPolicy: .bufferingNewest(64)) { continuation in
            let task = Task.detached(priority: .userInitiated) { [pop, interval] in
                while !Task.isCancelled {
                    var drained = false
                    while let b = pop() {
                        drained = true
                        continuation.yield(b)
                        if b.isFinal { continuation.finish(); return }
                    }
                    if !drained { try? await Task.sleep(for: interval) }
                }
                continuation.finish()
            }
            continuation.onTermination = { _ in task.cancel() }
        }
    }
}

public struct RecordedBurstSource: BurstSource {
    public let log: [Burst]
    public let spacing: Duration

    public init(log: [Burst], spacing: Duration = .zero) { self.log = log; self.spacing = spacing }

    public func bursts() -> AsyncStream<Burst> {
        AsyncStream { continuation in
            let task = Task {
                for b in log {
                    if Task.isCancelled { break }
                    if spacing > .zero { try? await Task.sleep(for: spacing) }
                    continuation.yield(b)
                    if b.isFinal { break }
                }
                continuation.finish()
            }
            continuation.onTermination = { _ in task.cancel() }
        }
    }
}
