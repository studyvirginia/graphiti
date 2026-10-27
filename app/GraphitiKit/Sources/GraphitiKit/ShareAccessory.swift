import SwiftUI

public struct ShareAccessory: View {
    let params: CurveParams
    let latex: String
    var onShare: () -> Void
    public init(params: CurveParams, latex: String, onShare: @escaping () -> Void = {}) {
        self.params = params; self.latex = latex; self.onShare = onShare
    }

    public var body: some View {
        HStack(spacing: Space.m) {
            ZStack {
                Theme.ink
                CurveShape(params: params, layout: .card).stroke(EquationColor.color(for: latex), lineWidth: 1.6)
            }
            .frame(width: 34, height: 34)
            .clipShape(.rect(cornerRadius: 8))
            VStack(alignment: .leading, spacing: 1) {
                Text("Math card ready").font(.subheadline.weight(.semibold))
                Text("Saved to You").font(.caption).foregroundStyle(.secondary)
            }
            Spacer(minLength: 0)
            Button("Share", systemImage: "square.and.arrow.up", action: onShare)
                .labelStyle(.iconOnly)
                .font(.system(size: 17, weight: .semibold))
        }
        .padding(.leading, 8).padding(.trailing, 16)
        .contentShape(Rectangle())
        .onTapGesture(perform: onShare)
        .accessibilityElement(children: .combine)
        .accessibilityLabel("Math card ready. Share")
        .accessibilityAddTraits(.isButton)
    }
}

public struct ShareSheet: View {
    let card: MathCard
    public init(card: MathCard) { self.card = card }
    public var body: some View {
        VStack(spacing: Space.xl) {
            card.padding(.horizontal, Space.l)
            HStack(spacing: Space.l) {
                destination("message.fill", "Messages")
                destination("camera", "Instagram")
                destination("music.note", "TikTok")
                destination("function", "Desmos")
                destination("ellipsis", "More")
            }
            VStack(spacing: Space.m) {
                PrimaryButton("Post to Discover")
                HStack(spacing: Space.m) {
                    Button {} label: { Label("Save Image", systemImage: "photo").frame(maxWidth: .infinity) }.secondaryButton().controlSize(.large)
                    Button {} label: { Label("Save Video", systemImage: "video").frame(maxWidth: .infinity) }.secondaryButton().controlSize(.large)
                }
            }
            .padding(.horizontal, Space.l)
        }
        .padding(.top, Space.xl)
    }
    private func destination(_ symbol: String, _ label: String) -> some View {
        VStack(spacing: Space.s) {
            Button(label, systemImage: symbol) {}
                .labelStyle(.iconOnly).font(.system(size: 20, weight: .medium))
                .frame(width: 34, height: 34)
                .buttonStyle(.glass).buttonBorderShape(.circle).controlSize(.extraLarge)
            Text(label).font(Typeface.label).foregroundStyle(.secondary)
        }
    }
}
