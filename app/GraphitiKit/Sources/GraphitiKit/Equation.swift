import SwiftUI

public struct EquationToken: Identifiable, Hashable, Sendable {
    public enum Kind: Sendable { case symbol, op, constant, exponent }
    public let id: String
    public let latex: String
    public let kind: Kind

    public init(_ id: String, _ latex: String, _ kind: Kind = .symbol) {
        self.id = id; self.latex = latex; self.kind = kind
    }
}

public extension Array where Element == EquationToken {
    var latex: String {
        map { t in
            switch t.kind {
            case .exponent: return "^{\(t.latex)}"
            case .op: return " \(t.latex) "
            default: return t.latex
            }
        }.joined()
    }
}

public struct EvolvingEquation: View {
    let tokens: [EquationToken]
    let size: CGFloat
    let highlight: Set<String>
    let tint: PlatformColor
    let dying: Set<String>
    let colours: [String: PlatformColor]
    let onTap: ((String) -> Void)?
    @Environment(\.accessibilityReduceMotion) private var reduceMotion

    public init(tokens: [EquationToken], size: CGFloat = 30, highlight: Set<String> = [],
                tint: PlatformColor = .white, dying: Set<String> = [],
                colours: [String: PlatformColor] = [:], onTap: ((String) -> Void)? = nil) {
        self.tokens = tokens; self.size = size; self.highlight = highlight; self.tint = tint; self.dying = dying
        self.colours = colours; self.onTap = onTap
    }

    public var body: some View {
        HStack(alignment: .mathBaseline, spacing: 0) {
            ForEach(tokens) { t in
                token(t)
                    .contentShape(Rectangle())
                    .onTapGesture { onTap?(t.id) }
                    .transition(reduceMotion ? .opacity : .asymmetric(
                        insertion: .scale(scale: 0.7, anchor: .bottomLeading)
                            .combined(with: .opacity)
                            .combined(with: .offset(y: t.kind == .exponent ? 10 : 6)),
                        removal: .opacity.combined(with: .offset(y: 10)).combined(with: .scale(scale: 0.85))))
            }
        }
        .fixedSize()
        .accessibilityElement(children: .ignore)
        .accessibilityLabel(spoken)
    }

    private func color(for t: EquationToken) -> PlatformColor {
        if dying.contains(t.id) { return .graphiti(hex: 0xFFFFFF, alpha: 0.28) }
        if let c = colours[t.id] { return c }
        return highlight.contains(t.id) ? tint : .white
    }

    @ViewBuilder private func token(_ t: EquationToken) -> some View {
        switch t.kind {
        case .op:
            MathText(t.latex, size: size, color: color(for: t)).padding(.horizontal, size * 0.2)
        case .exponent:
            MathText(t.latex, size: size * 0.68, color: color(for: t))
                .alignmentGuide(.mathBaseline) { d in d[.mathBaseline] + size * 0.42 }
        case .constant:
            ZStack {
                MathText(t.latex, size: size, color: color(for: t))
                    .id(t.latex)
                    .transition(reduceMotion ? .opacity : .push(from: .bottom))
            }
            .clipped()
            .alignmentGuide(.mathBaseline) { d in d[.mathBaseline] }
        case .symbol:
            MathText(t.latex, size: size, color: color(for: t))
        }
    }

    private var spoken: String {
        tokens.map { t in
            switch t.latex {
            case "=": return "equals"
            case "+": return "plus"
            case "-": return "minus"
            case "\\sin(": return "sine of"
            case ")": return ""
            case "\\frac{1}{2}": return "one half"
            default: return t.kind == .exponent ? "to the power \(t.latex)" : t.latex
            }
        }.filter { !$0.isEmpty }.joined(separator: " ")
    }
}
