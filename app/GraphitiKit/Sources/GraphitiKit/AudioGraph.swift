import SwiftUI
import Accessibility

public struct CurveAudioGraph: AXChartDescriptorRepresentable {
    let params: CurveParams
    let title: String
    let range: ClosedRange<Double>
    let points: [CGPoint]

    public init(params: CurveParams, title: String, range: ClosedRange<Double> = -4...4, points: [CGPoint] = []) {
        self.params = params; self.title = title; self.range = range; self.points = points
    }

    public func makeChartDescriptor() -> AXChartDescriptor {
        let n = 120
        let curve: [(Double, Double)] = points.count > 1
            ? points.map { (Double($0.x), -Double($0.y)) }
            : (0...n).map { i in let x = range.lowerBound + (range.upperBound - range.lowerBound) * Double(i) / Double(n); return (x, params.y(x)) }
        let xs = curve.map(\.0), ys = curve.map(\.1)
        let xRange = (xs.min() ?? range.lowerBound)...(xs.max() ?? range.upperBound)
        let data = zip(xs, ys).map { AXDataPoint(x: $0, y: $1) }
        let lo = ys.min() ?? -1, hi = ys.max() ?? 1
        let xAxis = AXNumericDataAxisDescriptor(title: "x", range: xRange, gridlinePositions: []) { String(format: "%.1f", $0) }
        let yAxis = AXNumericDataAxisDescriptor(title: "y", range: lo...hi, gridlinePositions: []) { String(format: "%.1f", $0) }
        let series = AXDataSeriesDescriptor(name: title, isContinuous: true, dataPoints: data)
        return AXChartDescriptor(title: title,
                                 summary: "Curve from x \(Int(range.lowerBound)) to \(Int(range.upperBound)), y between \(String(format: "%.1f", lo)) and \(String(format: "%.1f", hi)).",
                                 xAxis: xAxis, yAxis: yAxis, additionalAxes: [], series: [series])
    }
}
