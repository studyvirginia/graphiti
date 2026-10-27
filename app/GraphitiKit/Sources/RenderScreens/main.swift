import SwiftUI
import AppKit
import GraphitiKit


let out = URL(fileURLWithPath: CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "renders")
try? FileManager.default.createDirectory(at: out, withIntermediateDirectories: true)

if let i = CommandLine.arguments.firstIndex(of: "--appstore"), i + 1 < CommandLine.arguments.count {
    MainActor.assumeIsolated { renderAppStoreShots(from: URL(fileURLWithPath: CommandLine.arguments[i + 1]), to: out) }
    exit(0)
}

struct StatusBar: View {
    var body: some View {
        HStack {
            Text("9:41").font(.system(size: 17, weight: .semibold))
            Spacer()
            Image(systemName: "cellularbars"); Image(systemName: "battery.100")
        }
        .foregroundStyle(.white)
        .padding(.horizontal, 34).frame(height: 54)
        .overlay(Capsule().fill(.black).frame(width: 124, height: 36))
    }
}

struct MockTabBar: View {
    var selected = 0
    var body: some View {
        HStack {
            item("function", "Solve", 0); item("safari", "Discover", 1); item("person.crop.circle", "You", 2)
        }
        .padding(4)
        .background(Capsule().fill(Color(white: 0.16).opacity(0.88)))
        .overlay(Capsule().strokeBorder(LinearGradient(colors: [.white.opacity(0.22), .white.opacity(0.04)], startPoint: .top, endPoint: .bottom), lineWidth: 0.7))
        .padding(.horizontal, 22)
    }
    func item(_ symbol: String, _ label: String, _ i: Int) -> some View {
        VStack(spacing: 3) {
            Image(systemName: symbol).font(.system(size: 19, weight: .medium))
            Text(label).font(.system(size: 10, weight: i == selected ? .semibold : .medium))
        }
        .foregroundStyle(i == selected ? Theme.accent : .white)
        .frame(maxWidth: .infinity).frame(height: 54)
        .background(Capsule().fill(.white.opacity(i == selected ? 0.1 : 0)))
    }
}

struct MockTopBar: View {
    var undo = true, new = true
    var body: some View {
        HStack {
            button("arrow.uturn.backward", enabled: undo)
            Spacer()
            button("square.and.pencil", enabled: new)
        }
        .padding(.horizontal, 16).padding(.top, 58)
    }
    func button(_ symbol: String, enabled: Bool) -> some View {
        Image(systemName: symbol).font(.system(size: 17, weight: .medium))
            .foregroundStyle(.white.opacity(enabled ? 1 : 0.3))
            .frame(width: 44, height: 44)
            .background(Circle().fill(Color(white: 0.16).opacity(0.88)))
            .overlay(Circle().strokeBorder(LinearGradient(colors: [.white.opacity(0.22), .white.opacity(0.04)], startPoint: .top, endPoint: .bottom), lineWidth: 0.7))
    }
}

@MainActor
func phone<V: View>(_ content: V, tabBar: Bool = true, tab: Int = 0, topBar: MockTopBar? = nil, accessory: Bool = false) -> some View {
    ZStack(alignment: .top) {
        content
            .safeAreaPadding(.top, topBar == nil ? 54 : 106)
            .safeAreaPadding(.bottom, tabBar ? 92 : 34)
            .frame(width: 393, height: 852, alignment: .top).clipped()
        if let topBar { topBar }
        StatusBar()
        if tabBar {
            VStack(spacing: 8) {
                Spacer()
                if accessory {
                    ShareAccessory(params: CurveParams(m: 0.5, a: 1, w: 2, q: 0), latex: "y = \\frac{1}{2}x + \\sin(2x)")
                        .frame(height: 48)
                        .background(Capsule().fill(Color(white: 0.16).opacity(0.88)))
                        .overlay(Capsule().strokeBorder(LinearGradient(colors: [.white.opacity(0.22), .white.opacity(0.04)], startPoint: .top, endPoint: .bottom), lineWidth: 0.7))
                        .padding(.horizontal, 22)
                }
                MockTabBar(selected: tab).padding(.bottom, 26)
            }
        }
        VStack { Spacer(); Capsule().fill(.white).frame(width: 139, height: 5).padding(.bottom, 8) }
    }
    .frame(width: 393, height: 852)
    .background(.black)
    .environment(\.colorScheme, .dark)
}

@MainActor
func save<V: View>(_ view: V, _ name: String) {
    let r = ImageRenderer(content: view)
    r.scale = 3
    guard let cg = r.cgImage,
          let png = NSBitmapImageRep(cgImage: cg).representation(using: .png, properties: [:]) else {
        print("failed: \(name)"); return
    }
    try? png.write(to: out.appendingPathComponent("\(name).png"))
    print("wrote \(name).png")
}

MainActor.assumeIsolated {
    let bar = MockTopBar()
    save(phone(SolveScreen(session: SolveSession(phase: .empty)), topBar: MockTopBar(undo: false, new: false)), "swift-empty")
    save(phone(SolveScreen(session: SolveSession(phase: .empty), firstLaunchHint: 0.62), topBar: MockTopBar(undo: false, new: false)), "swift-empty-first-launch")
    for s in [0, 2, 5, 6] {
        save(phone(SolveScreen(session: SolveSession(step: s)), topBar: bar), "swift-evolving-\(s)")
    }
    save(phone(SolveScreen(session: SolveSession(step: 6, dying: ["p2", "c3", "x3", "e3"], elapsed: 1.9)), topBar: bar), "swift-prune")
    save(phone(SolveScreen(session: SolveSession(step: 8, locked: true, flash: 0.55)), topBar: bar), "swift-lock-flash")
    save(phone(SolveScreen(session: SolveSession(step: 8, locked: true, indicate: true)), topBar: bar), "swift-lock-indicate")
    save(phone(SolveScreen(session: SolveSession(step: 8, locked: true)), topBar: bar), "swift-locked")
    let card = MathCard(params: CurveParams(m: 0.5, a: 1, w: 2, q: 0), latex: "y = \\frac{1}{2}x + \\sin(2x)", author: "@lincoln")
    save(phone(SharePopup(card: card) { SolveScreen(session: SolveSession(step: 8, locked: true)) }, tabBar: false), "swift-share")
}


struct MessagesMock: View {
    let payload: CardPayload
    var editing: Bool
    var body: some View {
        VStack(spacing: 0) {
            VStack(spacing: 4) {
                Circle().fill(Color(white: 0.35)).frame(width: 44, height: 44)
                    .overlay(Text("A").font(.headline).foregroundStyle(.white))
                Text("Ana").font(.caption).foregroundStyle(.white)
            }
            .padding(.top, 60).frame(maxWidth: .infinity).padding(.bottom, 10)
            .background(Color(white: 0.09))
            VStack(alignment: .leading, spacing: 10) {
                bubble("did u see the math thing", mine: false)
                bubble("draw anything and it finds the equation", mine: true)
                VStack(alignment: .leading, spacing: 0) {
                    MathCard(params: payload.params, latex: payload.latex, author: "@" + payload.author)
                        .shadow(radius: 0).scaleEffect(0.72, anchor: .topLeading)
                        .frame(width: 245, height: 270, alignment: .topLeading)
                    VStack(alignment: .leading, spacing: 2) {
                        Text(EquationFormatter.plain(payload.params)).font(.subheadline.weight(.semibold)).foregroundStyle(.white)
                        Text("Tap to change it").font(.caption).foregroundStyle(.white.opacity(0.6))
                    }
                    .padding(10).frame(width: 245, alignment: .leading).background(Color(white: 0.17))
                }
                .clipShape(RoundedRectangle(cornerRadius: 18, style: .continuous))
                .frame(maxWidth: .infinity, alignment: .trailing)
                if !editing { bubble("ok drag the wiggle and send it back", mine: true) }
            }
            .padding(14)
            Spacer()
        }
        .overlay(alignment: .bottom) {
            if editing {
                VStack(spacing: 0) {
                    Capsule().fill(.white.opacity(0.3)).frame(width: 36, height: 5).padding(.top, 8)
                    LiveCardEditor(payload: CardPayload(params: CurveParams(m: 0.5, a: 1.6, w: 3, q: 0), author: "ana"))
                }
                .background(Theme.ink)
                .clipShape(UnevenRoundedRectangle(topLeadingRadius: 30, topTrailingRadius: 30))
            }
        }
    }
    func bubble(_ t: String, mine: Bool) -> some View {
        Text(t).font(.body).foregroundStyle(.white)
            .padding(.horizontal, 14).padding(.vertical, 9)
            .background(RoundedRectangle(cornerRadius: 18, style: .continuous).fill(mine ? Color(red: 0.04, green: 0.52, blue: 1) : Color(white: 0.17)))
            .frame(maxWidth: .infinity, alignment: mine ? .trailing : .leading)
    }
}

struct DragProof: View {
    var body: some View {
        let p = CurveParams(m: 0.5, a: 1, w: 2, q: 0)
        let x = 0.785
        let q = CurveEditor.drag(p, at: x, to: p.y(x) + 0.9)
        return GeometryReader { geo in
            let f = PlaneFrame(size: geo.size)
            ZStack {
                Theme.ink
                ManimPlane()
                CurveShape(params: p).stroke(Theme.blueC.opacity(0.28), style: StrokeStyle(lineWidth: 2.5, lineCap: .round, dash: [6, 6]))
                CurveShape(params: q).stroke(Theme.blueC, style: StrokeStyle(lineWidth: 4, lineCap: .round))
                Circle().fill(.white.opacity(0.18)).frame(width: 46, height: 46).position(f.point(x, q.y(x)))
                Circle().fill(.white).frame(width: 12, height: 12).position(f.point(x, q.y(x)))
                VStack(spacing: 12) {
                    Spacer()
                    Text("Dragging the peak").font(.subheadline.weight(.medium)).foregroundStyle(Theme.secondary)
                    HStack(alignment: .mathBaseline, spacing: 0) {
                        MathText("y = \\frac{1}{2}x", size: 30)
                        MathText("+", size: 30).padding(.horizontal, 6)
                        MathText(String(format: "%.2f", q.a), size: 30, color: .graphiti(hex: 0x58C4DD))
                        MathText("\\sin(", size: 30)
                        MathText(String(format: "%.2f", q.w), size: 30, color: .graphiti(hex: 0x58C4DD))
                        MathText("x)", size: 30)
                    }
                    Text("Height moves most; the wiggle adjusts a little").font(.caption).foregroundStyle(.white.opacity(0.35))
                }
                .padding(.bottom, 118)
            }
        }
    }
}

struct SnapProof: View {
    var body: some View {
        ZStack {
            Theme.ink
            ManimPlane(fadeStart: 0.52)
            CurveShape(params: CurveParams(m: 0, a: 1.5, w: .pi / 2, q: 0)).stroke(Theme.blueC, style: StrokeStyle(lineWidth: 4, lineCap: .round))
            VStack(spacing: 14) {
                Spacer()
                Text("Generation 31").font(.subheadline.weight(.medium)).foregroundStyle(Theme.secondary)
                MathText("y = 1.4998\\sin(1.5709x)", size: 26, color: .graphiti(hex: 0xFFFFFF, alpha: 0.45))
                Image(systemName: "arrow.down").foregroundStyle(Theme.secondary)
                MathText(EquationFormatter.latex(CurveParams(m: 0, a: 1.4998, w: 1.5709, q: 0)), size: 32)
                Text("Found in 2.1 s · snapped to exact").font(.caption).foregroundStyle(.white.opacity(0.35))
            }
            .padding(.bottom, 118)
        }
    }
}

struct NamedProof: View {
    var body: some View {
        let latex = "r = 1.6(1 - \\cos\\theta)"
        return ZStack {
            Theme.ink
            ManimPlane()
            PolarCurveShape { 1.6 * (1 - cos($0)) }.stroke(Theme.tealC, style: StrokeStyle(lineWidth: 4, lineCap: .round))
            VStack(spacing: 14) {
                Spacer()
                if let named = NamedCurve.identify(latex) { NamedCurvePlate(named).padding(.horizontal, 20) }
                Text("Found in 1.8 s").font(.subheadline.weight(.medium)).foregroundStyle(Theme.secondary)
                MathText(latex, size: 30)
            }
            .padding(.bottom, 118)
        }
    }
}

MainActor.assumeIsolated {
    let payload = CardPayload(params: CurveParams(m: 0.5, a: 1, w: 2, q: 0), author: "lincoln")
    save(phone(MessagesMock(payload: payload, editing: false), tabBar: false), "poc-imessage-thread")
    save(phone(MessagesMock(payload: payload, editing: true), tabBar: false), "poc-imessage-edit")
    save(phone(DragProof()), "poc-drag-curve")
    save(phone(SnapProof()), "poc-snap-exact")
    save(phone(NamedProof()), "poc-named-curve")
}


struct BlendShape: Shape {
    var from: CurveParams, to: CurveParams, t: Double
    func path(in rect: CGRect) -> Path {
        let f = PlaneFrame(size: rect.size)
        var p = Path()
        for i in 0...360 {
            let x = f.xRange.lowerBound + (f.xRange.upperBound - f.xRange.lowerBound) * Double(i) / 360
            let pt = f.point(x, (1 - t) * from.y(x) + t * to.y(x))
            if i == 0 { p.move(to: pt) } else { p.addLine(to: pt) }
        }
        return p
    }
}

struct MorphFrame: View {
    let from: CurveParams, to: CurveParams, t: Double, blend: Bool, title: String
    var body: some View {
        let e = t * t * (3 - 2 * t)
        var mid = from; var d = to - from; d.scale(by: e); mid += d
        return ZStack(alignment: .top) {
            Theme.ink
            ManimPlane(fadeStart: 1)
            CurveShape(params: from).stroke(.white.opacity(0.18), style: StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
            CurveShape(params: to).stroke(.white.opacity(0.18), style: StrokeStyle(lineWidth: 1.5, dash: [4, 4]))
            let style = StrokeStyle(lineWidth: 3.5, lineCap: .round)
            if blend { BlendShape(from: from, to: to, t: e).stroke(Theme.blueC, style: style) }
            else { CurveShape(params: mid).stroke(Theme.blueC, style: style) }
            Text(title).font(.system(size: 13, weight: .semibold)).foregroundStyle(.white.opacity(0.8)).padding(.top, 10)
        }
        .frame(width: 300, height: 220)
    }
}

if CommandLine.arguments.contains("--morph") { MainActor.assumeIsolated {
    let dir = out.appendingPathComponent("morph")
    try? FileManager.default.createDirectory(at: dir, withIntermediateDirectories: true)
    let cases: [(String, CurveParams, CurveParams)] = [
        ("tune", CurveParams(m: 0.4, a: 1, w: 1, q: 0), CurveParams(m: 0.4, a: 1, w: 1.7, q: 0)),
        ("newterm", CurveParams(m: 0.4, a: 0, w: 0, q: 0), CurveParams(m: 0.4, a: 1, w: 1, q: 0)),
        ("prune", CurveParams(m: 0.5, a: 1, w: 2, q: 0.1), CurveParams(m: 0.5, a: 1, w: 2, q: 0)),
    ]
    for (name, a, b) in cases {
        for blend in [false, true] {
            for k in 0...24 {
                let t = Double(k) / 24
                let r = ImageRenderer(content: MorphFrame(from: a, to: b, t: t, blend: blend,
                                                          title: blend ? "Point blend (Manim)" : "Parameter morph (now)"))
                r.scale = 2
                if let cg = r.cgImage, let png = NSBitmapImageRep(cgImage: cg).representation(using: .png, properties: [:]) {
                    try? png.write(to: dir.appendingPathComponent("\(name)-\(blend ? "blend" : "param")-\(String(format: "%02d", k)).png"))
                }
            }
        }
    }
    print("wrote morph frames")
} }


MainActor.assumeIsolated {
    let cards: [(CurveParams, String)] = [
        (CurveParams(m: 0.5, a: 1, w: 2, q: 0), "y = \\frac{1}{2}x + \\sin(2x)"),
        (CurveParams(m: 0, a: 1.5, w: .pi / 2, q: 0), "y = \\frac{3}{2}\\sin(\\frac{\\pi}{2}x)"),
        (CurveParams(m: 0, a: 0.8, w: 3, q: 0.12), "y = 0.12x^2 + 0.8\\sin(3x)"),
        (CurveParams(m: -0.3, a: 1.2, w: 1, q: 0), "y = -0.3x + 1.2\\sin(x)"),
        (CurveParams(m: 0.2, a: 0, w: 1, q: 0.2), "y = 0.2x^2 + 0.2x"),
        (CurveParams(m: 0, a: 2, w: 0.8, q: 0), "y = 2\\sin(0.8x)"),
        (CurveParams(m: 1, a: 0.5, w: 4, q: 0), "y = x + \\frac{1}{2}\\sin(4x)"),
        (CurveParams(m: 0, a: 1, w: 2, q: -0.1), "y = \\sin(2x) - 0.1x^2"),
    ]
    let grid = VStack(spacing: 14) {
        ForEach(0..<4, id: \.self) { row in
            HStack(spacing: 14) {
                ForEach(0..<2, id: \.self) { col in
                    let c = cards[row * 2 + col]
                    MathCard(params: c.0, latex: c.1, author: "@lincoln")
                        .frame(width: 340).scaleEffect(0.5).frame(width: 170, height: 182)
                }
            }
        }
    }
    .padding(20).background(Color.black)
    save(grid.environment(\.colorScheme, .dark), "equation-colours")
}

