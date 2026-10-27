import SwiftUI

public struct GridCrossingDetector {
    public enum Event: Equatable, Sendable { case gridLine, axis }

    let frame: PlaneFrame
    let spacing: Double
    private var last: (gx: Int, gy: Int)?

    public init(frame: PlaneFrame, spacing: Double = 0.5) {
        self.frame = frame; self.spacing = spacing
    }

    public mutating func add(_ point: CGPoint) -> Event? {
        let x = Double((point.x - frame.origin.x) / frame.unit)
        let y = Double((frame.origin.y - point.y) / frame.unit)
        let cell = (gx: Int((x / spacing).rounded(.down)), gy: Int((y / spacing).rounded(.down)))
        defer { last = cell }
        guard let prev = last, prev != cell else { return nil }
        let crossedXAxis = (prev.gy < 0) != (cell.gy < 0)
        let crossedYAxis = (prev.gx < 0) != (cell.gx < 0)
        return crossedXAxis || crossedYAxis ? .axis : .gridLine
    }
}

public struct DrawingCanvas: View {
    var onLift: ([CGPoint]) -> Void
    @State private var stroke: [CGPoint] = []
    @State private var detector: GridCrossingDetector?
    @State private var tick = 0
    @State private var axisTick = 0

    public init(onLift: @escaping ([CGPoint]) -> Void) { self.onLift = onLift }

    public var body: some View {
        GeometryReader { geo in
            Path { p in
                guard let first = stroke.first else { return }
                p.move(to: first); stroke.dropFirst().forEach { p.addLine(to: $0) }
            }
            .stroke(.white.opacity(0.55), style: StrokeStyle(lineWidth: 4, lineCap: .round, lineJoin: .round))
            .contentShape(Rectangle())
            .gesture(DragGesture(minimumDistance: 0)
                .onChanged { v in
                    if stroke.isEmpty { detector = GridCrossingDetector(frame: PlaneFrame(size: geo.size)) }
                    stroke.append(v.location)
                    switch detector?.add(v.location) {
                    case .gridLine: tick += 1
                    case .axis: axisTick += 1
                    case nil: break
                    }
                }
                .onEnded { _ in onLift(stroke); stroke = [] })
        }
        #if os(iOS)
        .sensoryFeedback(.impact(flexibility: .soft, intensity: 0.35), trigger: tick)
        .sensoryFeedback(.impact(flexibility: .rigid, intensity: 0.8), trigger: axisTick)
        #endif
    }
}
