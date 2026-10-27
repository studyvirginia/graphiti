import SwiftUI

public struct LiveCardEditor: View {
    @State var payload: CardPayload
    var onSend: (CardPayload) -> Void

    public init(payload: CardPayload, onSend: @escaping (CardPayload) -> Void = { _ in }) {
        _payload = State(initialValue: payload); self.onSend = onSend
    }

    public var body: some View {
        VStack(spacing: 18) {
            MathCard(params: payload.params, latex: payload.latex, author: payload.author.isEmpty ? "" : "@" + payload.author)
                .shadow(radius: 0)
            VStack(spacing: 14) {
                if abs(payload.params.a) > 1e-9 {
                    constantSlider("Height", value: $payload.params.a, range: 0...3)
                }
                if abs(payload.params.w) > 1e-9 {
                    constantSlider("Wiggle", value: $payload.params.w, range: 0.5...5)
                }
                if abs(payload.params.m) > 1e-9 {
                    constantSlider("Tilt", value: $payload.params.m, range: -2...2)
                }
            }
            .padding(.horizontal, 6)
            PrimaryButton("Send your version") { onSend(payload) }
        }
        .padding(20)
        .background(Theme.ink)
    }

    private func constantSlider(_ label: String, value: Binding<Double>, range: ClosedRange<Double>) -> some View {
        HStack(spacing: 12) {
            Text(label).font(Typeface.status).foregroundStyle(Theme.secondary).frame(width: 52, alignment: .leading)
            ConstantSlider(value: value, range: range)
            Text(String(format: "%.2f", value.wrappedValue)).font(.subheadline.weight(.semibold)).monospacedDigit()
                .foregroundStyle(.white).frame(width: 40, alignment: .trailing)
        }
    }
}

public struct ConstantSlider: View {
    @Binding var value: Double
    let range: ClosedRange<Double>
    public init(value: Binding<Double>, range: ClosedRange<Double>) { _value = value; self.range = range }

    public var body: some View {
        GeometryReader { geo in
            let w = geo.size.width
            let t = (value - range.lowerBound) / (range.upperBound - range.lowerBound)
            ZStack(alignment: .leading) {
                Capsule().fill(.white.opacity(0.18)).frame(height: 4)
                Capsule().fill(Theme.blueC).frame(width: max(4, w * t), height: 4)
                Circle().fill(.white).frame(width: 26, height: 26)
                    .shadow(color: .black.opacity(0.4), radius: 2, y: 1)
                    .offset(x: w * t - 13)
            }
            .frame(maxHeight: .infinity)
            .contentShape(Rectangle())
            .gesture(DragGesture(minimumDistance: 0).onChanged { g in
                let f = min(1, max(0, g.location.x / w))
                value = range.lowerBound + f * (range.upperBound - range.lowerBound)
            })
        }
        .frame(height: 30)
        .accessibilityRepresentation { Slider(value: $value, in: range) }
    }
}

extension EquationFormatter {
    public static func plain(_ p: CurveParams) -> String {
        var s = latex(p).replacingOccurrences(of: "{", with: "").replacingOccurrences(of: "}", with: "")
        let vulgar = ["12": "½", "13": "⅓", "23": "⅔", "14": "¼", "34": "¾", "15": "⅕", "16": "⅙", "18": "⅛"]
        while let r = s.range(of: #"\\frac\d\d"#, options: .regularExpression) {
            let nm = String(s[r].suffix(2))
            s.replaceSubrange(r, with: vulgar[nm] ?? "\(nm.first!)/\(nm.last!)")
        }
        return s
            .replacingOccurrences(of: "\\pi", with: "π")
            .replacingOccurrences(of: "\\sin", with: "sin")
            .replacingOccurrences(of: "x^2", with: "x²")
    }
}

#if os(iOS) && canImport(Messages)
import Messages

open class GraphitiMessagesViewController: MSMessagesAppViewController {
    private var host: UIHostingController<AnyView>?

    open override func willBecomeActive(with conversation: MSConversation) {
        super.willBecomeActive(with: conversation)
        let incoming = conversation.selectedMessage?.url.flatMap(CardPayload.init(url:))
        let payload = incoming ?? CardPayload(params: CurveParams(m: 0.5, a: 1, w: 2, q: 0), author: "")
        show(LiveCardEditor(payload: payload) { [weak self] edited in self?.send(edited, in: conversation) })
    }

    private func show<V: View>(_ view: V) {
        host?.view.removeFromSuperview(); host?.removeFromParent()
        let h = UIHostingController(rootView: AnyView(view.preferredColorScheme(.dark)))
        addChild(h); h.view.frame = self.view.bounds
        h.view.autoresizingMask = [.flexibleWidth, .flexibleHeight]
        self.view.addSubview(h.view); h.didMove(toParent: self)
        host = h
    }

    @MainActor private func send(_ payload: CardPayload, in conversation: MSConversation) {
        let message = MSMessage(session: conversation.selectedMessage?.session ?? MSSession())
        let layout = MSMessageTemplateLayout()
        let renderer = ImageRenderer(content: MathCard(params: payload.params, latex: payload.latex, author: payload.author.isEmpty ? "" : "@" + payload.author)
            .frame(width: 337).environment(\.colorScheme, .dark))
        renderer.scale = 3
        layout.image = renderer.uiImage
        layout.caption = EquationFormatter.plain(payload.params)
        layout.subcaption = "Tap to change it"
        message.layout = layout
        message.url = payload.url
        conversation.insert(message)
        requestPresentationStyle(.compact)
    }
}
#endif
