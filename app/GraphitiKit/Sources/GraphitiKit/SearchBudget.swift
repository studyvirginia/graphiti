import Foundation

public struct SearchBudget: Equatable, Sendable {
    public static let defaultSeconds: Double = 8
    public static let range: ClosedRange<Double> = 2...60

    public var seconds: Double
    public var generations: Int?

    public init(seconds: Double = SearchBudget.defaultSeconds, generations: Int? = nil) {
        self.seconds = min(max(seconds, Self.range.lowerBound), Self.range.upperBound)
        self.generations = generations
    }
}
