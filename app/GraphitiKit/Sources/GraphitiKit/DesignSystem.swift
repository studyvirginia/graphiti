import SwiftUI


public enum Space {
    public static let xs: CGFloat = 4
    public static let s: CGFloat = 8
    public static let m: CGFloat = 12
    public static let l: CGFloat = 16
    public static let xl: CGFloat = 24
    public static let xxl: CGFloat = 32
}

public enum Radius {
    public static let card: CGFloat = 26
    public static let heroCard: CGFloat = 32
    public static let cardPadding: CGFloat = 8
    public static var plot: CGFloat { card - cardPadding }
    public static var heroPlot: CGFloat { heroCard - cardPadding }
    public static let notice: CGFloat = 20
}

public enum Typeface {
    public static let section = Font.title2.weight(.bold)
    public static let headline = Font.headline
    public static let status = Font.subheadline.weight(.medium).monospacedDigit()
    public static let body = Font.body
    public static let meta = Font.footnote
    public static let hint = Font.footnote
    public static let label = Font.caption.weight(.medium)

    public static let equationSize: CGFloat = 34
    public static let cardEquationSize: CGFloat = 22
    public static let thumbEquationSize: CGFloat = 16
}

// MARK: - Surfaces

public extension View {
    func cardSurface(radius: CGFloat = Radius.card) -> some View {
        self.background(Theme.surface, in: .rect(cornerRadius: radius))
            .containerShape(.rect(cornerRadius: radius))
    }
}

public struct PlotView<Curve: View>: View {
    let curve: Curve
    public init(radius: CGFloat = Radius.plot, @ViewBuilder curve: () -> Curve) { self.curve = curve() }
    public var body: some View {
        ZStack { Theme.ink; ManimPlane(layout: .card); curve }
            .clipShape(ConcentricRectangle(corners: .concentric, isUniform: true))
    }
}


public extension EnvironmentValues {
    @Entry var screenAccent: Color = Theme.accent
}


public extension View {
    func primaryButton() -> some View {
        self.buttonStyle(.glassProminent).controlSize(.large).buttonBorderShape(.capsule)
    }
    func secondaryButton() -> some View {
        self.buttonStyle(.glass).buttonBorderShape(.capsule).foregroundStyle(.primary)
    }
}

public struct PrimaryButton: View {
    @Environment(\.screenAccent) private var accent
    let title: String
    var systemImage: String?
    var action: () -> Void
    public init(_ title: String, systemImage: String? = nil, action: @escaping () -> Void = {}) {
        self.title = title; self.systemImage = systemImage; self.action = action
    }
    public var body: some View {
        Button(action: action) {
            Group {
                if let systemImage { Label(title, systemImage: systemImage) } else { Text(title) }
            }
            .font(.headline)
            .foregroundStyle(.black)
            .frame(maxWidth: .infinity)
        }
        .primaryButton()
        .tint(accent)
    }
}

public struct IconButton: View {
    let systemName: String
    let title: String
    let size: CGFloat
    var action: () -> Void
    public init(_ systemName: String, _ title: String = "", size: CGFloat = 44, action: @escaping () -> Void = {}) {
        self.systemName = systemName; self.title = title; self.size = size; self.action = action
    }
    public var body: some View {
        Button(action: action) {
            Image(systemName: systemName)
                .font(.system(size: size * 0.39, weight: .medium))
                .frame(width: size, height: size)
                .contentShape(.circle)
        }
        .buttonStyle(.plain)
        .glassEffect(.regular.interactive(), in: .circle)
        .foregroundStyle(.primary)
        .accessibilityLabel(title)
    }
}

public struct FilterChip: View {
    @Environment(\.screenAccent) private var accent
    let text: String
    var selected: Bool
    public init(_ text: String, selected: Bool = false) { self.text = text; self.selected = selected }
    public var body: some View {
        if selected {
            Button {} label: { Text(text).foregroundStyle(.black) }
                .buttonStyle(.glassProminent).buttonBorderShape(.capsule).tint(accent)
        } else {
            Button(text) {}.buttonStyle(.glass).buttonBorderShape(.capsule).foregroundStyle(.primary)
        }
    }
}

public struct Tag: View {
    let text: String
    public init(_ text: String) { self.text = text }
    public var body: some View {
        Text(text)
            .font(Typeface.label).foregroundStyle(.secondary)
            .padding(.horizontal, 8).padding(.vertical, 3)
            .background(.quaternary, in: .capsule)
    }
}

public struct SectionHeader: View {
    let title: String
    var chevron: Bool
    public init(_ title: String, chevron: Bool = true) { self.title = title; self.chevron = chevron }
    public var body: some View {
        HStack(spacing: 4) {
            Text(title).font(Typeface.section)
            if chevron { Image(systemName: "chevron.forward").font(.title3.weight(.semibold)).foregroundStyle(.secondary) }
        }
        .accessibilityAddTraits(.isHeader)
    }
}
