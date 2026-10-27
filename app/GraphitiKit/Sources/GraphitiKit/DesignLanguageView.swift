import SwiftUI

public struct DesignLanguageView: View {
    public init() {}
    public var body: some View {
        NavigationStack {
            ScrollView {
                VStack(alignment: .leading, spacing: Space.xl) {
                    group("Type: system text styles; Latin Modern for maths") {
                        Text("Curve of the Week").font(Typeface.section)
                        Text("Found in 2.4 s").font(Typeface.status).foregroundStyle(.secondary)
                        Text("@physicsjay · 3.3k likes").font(Typeface.meta).foregroundStyle(.secondary)
                        HStack(alignment: .mathBaseline, spacing: Space.l) {
                            MathText("y = \\frac{1}{2}x", size: Typeface.equationSize)
                            MathText("y = \\frac{1}{2}x", size: Typeface.cardEquationSize)
                            MathText("y = \\frac{1}{2}x", size: Typeface.thumbEquationSize)
                        }
                    }
                    group("Buttons: Liquid Glass, prominent once per screen") {
                        PrimaryButton("Open in Solve", systemImage: "function")
                        HStack(spacing: Space.m) {
                            Button("Open") {}.secondaryButton()
                            IconButton("info", "Info"); IconButton("square.and.arrow.up", "Share"); IconButton("heart", "Like")
                        }
                        HStack(spacing: Space.s) { FilterChip("All", selected: true); FilterChip("Periodic"); Tag("Odd"); Tag("Polar") }
                    }
                    group("Corners: outer radius only; insides are concentric") {
                        HStack(alignment: .top, spacing: Space.m) {
                            PlotView { CurveShape(params: CurveParams(m: 0.5, a: 1, w: 2, q: 0)).stroke(EquationColor.color(for: "a"), lineWidth: 2.5) }
                                .frame(width: 170, height: 100).padding(Radius.cardPadding).cardSurface(radius: Radius.heroCard)
                            PlotView { CurveShape(params: CurveParams(m: 0, a: 1.5, w: 1.2, q: 0)).stroke(EquationColor.color(for: "b"), lineWidth: 2.5) }
                                .frame(width: 120, height: 100).padding(Radius.cardPadding).cardSurface()
                        }
                    }
                    group("Colour: one accent per screen; blue while searching, then the equation's own hue") {
                        HStack(spacing: Space.s) {
                            swatch(Theme.ink); swatch(Theme.surface); swatch(Theme.surfaceRaised)
                            swatch(Theme.blueC)
                            ForEach(["p", "q", "r"], id: \.self) { swatch(EquationColor.color(for: $0)) }
                        }
                    }
                }
                .padding(.horizontal, Space.l)
            }
            .navigationTitle("Design Language")
        }
        .preferredColorScheme(.dark)
        .tint(Theme.accent)
    }

    private func group<C: View>(_ title: String, @ViewBuilder _ c: () -> C) -> some View {
        VStack(alignment: .leading, spacing: Space.m) {
            Text(title).font(Typeface.label).foregroundStyle(.secondary).textCase(.uppercase)
            c()
        }
    }
    private func swatch(_ c: Color) -> some View {
        RoundedRectangle(cornerRadius: 10).fill(c).frame(width: 34, height: 34)
            .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(Theme.hairline, lineWidth: 0.5))
    }
}
