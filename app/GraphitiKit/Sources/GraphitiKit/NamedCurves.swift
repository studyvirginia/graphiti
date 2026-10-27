import SwiftUI

public struct NamedCurve: Equatable, Sendable {
    public let name: String
    public let history: String
    let pattern: String

    public static let catalog: [NamedCurve] = [
        .init(name: "Cardioid", history: "Named by Johann Castillon in 1741, for its heart shape.",
              pattern: #"^r=[\d.\\a-z{}]*\(1[-+]\\(cos|sin)\\theta\)$"#),
        .init(name: "Lemniscate of Bernoulli", history: "Described by Jakob Bernoulli in 1694, a figure eight.",
              pattern: #"^r\^\{?2\}?=[\d.\\a-z{}]*\\cos\(2\\theta\)$"#),
        .init(name: "Archimedean spiral", history: "Studied by Archimedes in On Spirals, around 225 BC.",
              pattern: #"^r=[\d.\\a-z{}]*\\theta$"#),
        .init(name: "Catenary", history: "The shape of a hanging chain. Named by Christiaan Huygens in 1690.",
              pattern: #"^y=[\d.\\a-z{}]*\\cosh\(.*\)$"#),
        .init(name: "Witch of Agnesi", history: "Studied by Maria Gaetana Agnesi in her 1748 textbook.",
              pattern: #"^y=\\frac\{[\d.]+\}\{x\^\{?2\}?\+[\d.]+\}$"#),
    ]

    public static func identify(_ latex: String) -> NamedCurve? {
        let s = latex.replacingOccurrences(of: " ", with: "")
        return catalog.first { s.range(of: $0.pattern, options: .regularExpression) != nil }
    }
}

public struct NamedCurvePlate: View {
    let curve: NamedCurve
    public init(_ curve: NamedCurve) { self.curve = curve }

    public var body: some View {
        HStack(spacing: 12) {
            Image(systemName: "seal.fill").font(.system(size: 22)).foregroundStyle(Theme.goldC)
            VStack(alignment: .leading, spacing: 2) {
                Text("You drew a \(curve.name.lowercased())").font(.subheadline.weight(.semibold)).foregroundStyle(.white)
                Text(curve.history).font(.caption).foregroundStyle(Theme.secondary).fixedSize(horizontal: false, vertical: true)
            }
            Spacer(minLength: 0)
        }
        .padding(14)
        .cardSurface(radius: Radius.notice)
        .accessibilityElement(children: .combine)
    }
}

public struct PolarCurveShape: Shape {
    var r: @Sendable (Double) -> Double
    var turns: Double = 1
    var layout: PlaneLayout = .standard
    public init(turns: Double = 1, layout: PlaneLayout = .standard, r: @escaping @Sendable (Double) -> Double) {
        self.r = r; self.turns = turns; self.layout = layout
    }

    public func path(in rect: CGRect) -> Path {
        let f = PlaneFrame(size: rect.size, layout: layout)
        var p = Path()
        let n = 480
        for i in 0...n {
            let t = 2 * .pi * turns * Double(i) / Double(n)
            let pt = f.point(r(t) * cos(t), r(t) * sin(t))
            if i == 0 { p.move(to: pt) } else { p.addLine(to: pt) }
        }
        return p
    }
}
