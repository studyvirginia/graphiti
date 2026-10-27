import SwiftUI

public struct DiscoverItem: Identifiable, Sendable {
    public enum Curve: Sendable {
        case function(CurveParams)
        case polar(@Sendable (Double) -> Double)
        case parametric(@Sendable (Double) -> (Double, Double))
        case points([CGPoint])
    }
    public let id: String
    public let latex: String
    public let curve: Curve
    public let name: String?
    public let author: String
    var handle: String? { author.isEmpty ? nil : "@" + author }
    public let likes: Int
    public let tags: [String]
    public let stroke: [CGPoint]

    public init(_ id: String, _ latex: String, _ curve: Curve, name: String? = nil, author: String, likes: Int, tags: [String],
                stroke: [CGPoint] = []) {
        self.id = id; self.latex = latex; self.curve = curve; self.name = name
        self.author = author; self.likes = likes; self.tags = tags; self.stroke = stroke
    }

    public var colour: Color { EquationColor.color(for: latex) }
    public var likesText: String { likes >= 1000 ? String(format: "%.1fk", Double(likes) / 1000) : "\(likes)" }

    public static func tags(for p: CurveParams) -> [String] {
        var t: [String] = []
        if p.a != 0 && p.m == 0 && p.q == 0 { t.append("Periodic") }
        if p.a == 0 { t.append("Polynomial") }
        if p.q == 0 { t.append("Odd") } else if p.m == 0 && p.a == 0 { t.append("Even") }
        if p.a != 0 && (p.m != 0 || p.q != 0) { t.append("Wave on a trend") }
        return t
    }

    public static let samples: [DiscoverItem] = [
        DiscoverItem("lissajous", "x = \\sin 3t,\\ y = \\sin 2t", .parametric { t in (2.6 * sin(3 * t), 2.2 * sin(2 * t)) },
                     name: "Lissajous", author: "physicsjay", likes: 3300, tags: ["Parametric", "Periodic", "Closed", "3:2 frequency"]),
        DiscoverItem("wave", "y = \\frac{1}{2}x + \\sin(2x)", .function(CurveParams(m: 0.5, a: 1, w: 2, q: 0)),
                     author: "lincoln", likes: 1200, tags: tags(for: CurveParams(m: 0.5, a: 1, w: 2, q: 0))),
        DiscoverItem("cardioid", "r = 1.4(1 - \\sin\\theta)", .polar { 1.4 * (1 - sin($0)) },
                     name: "Cardioid", author: "mrs.ortiz", likes: 842, tags: ["Polar", "Closed", "Symmetric"]),
        DiscoverItem("spiral", "r = 0.28\\theta", .polar { 0.28 * $0 },
                     name: "Archimedean spiral", author: "noor", likes: 517, tags: ["Polar", "Unbounded"]),
        DiscoverItem("bump", "y = 0.2x^2", .function(CurveParams(m: 0, a: 0, w: 1, q: 0.2)),
                     author: "ana.k", likes: 391, tags: ["Polynomial", "Even"]),
        DiscoverItem("chirp", "y = 2\\sin(0.8x)", .function(CurveParams(m: 0, a: 2, w: 0.8, q: 0)),
                     author: "dev", likes: 268, tags: tags(for: CurveParams(m: 0, a: 2, w: 0.8, q: 0))),
    ]
}

struct ItemCurve: View {
    let item: DiscoverItem
    var lineWidth: CGFloat = 3
    var showStroke = false
    var body: some View {
        let style = StrokeStyle(lineWidth: lineWidth, lineCap: .round, lineJoin: .round)
        switch item.curve {
        case .function(let p):
            CurveShape(params: p, layout: .card).stroke(item.colour, style: style)
        case .polar(let r):
            PolarCurveShape(turns: item.id == "spiral" ? 3 : 1, layout: .card, r: r).stroke(item.colour, style: style)
        case .parametric(let xy):
            ParametricCurveShape(xy: xy).stroke(item.colour, style: style)
        case .points(let pts):
            ZStack {
                if showStroke && item.stroke.count > 1 {
                    FittedPolyline(points: item.stroke, frame: pts + item.stroke).stroke(.white.opacity(0.35), style: StrokeStyle(lineWidth: lineWidth * 1.7, lineCap: .round, lineJoin: .round))
                }
                FittedPolyline(points: pts, frame: pts + item.stroke).stroke(item.colour, style: style)
            }
        }
    }
}

public struct ParametricCurveShape: Shape {
    var xy: @Sendable (Double) -> (Double, Double)
    var layout: PlaneLayout = .card
    public init(layout: PlaneLayout = .card, xy: @escaping @Sendable (Double) -> (Double, Double)) { self.xy = xy; self.layout = layout }
    public func path(in rect: CGRect) -> Path {
        let f = PlaneFrame(size: rect.size, layout: layout)
        var p = Path()
        for i in 0...600 {
            let (x, y) = xy(2 * .pi * Double(i) / 600)
            if i == 0 { p.move(to: f.point(x, y)) } else { p.addLine(to: f.point(x, y)) }
        }
        return p
    }
}


public struct DiscoverScreen: View {
    let sampleItems: [DiscoverItem]
    @ObservedObject private var library = Library.shared
    var items: [DiscoverItem] { [sampleItems[0]] + library.posted.map(\.item) + sampleItems.dropFirst() }
    var selectedTag = "All"
    var scrolls: Bool
    @State private var path: [String]
    @Namespace private var zoom
    @State private var query = ""
    public init(items: [DiscoverItem] = DiscoverItem.samples, scrolls: Bool = true, initialPath: [String] = []) {
        self.sampleItems = items; self.scrolls = scrolls; _path = State(initialValue: initialPath)
    }

    public var body: some View {
        NavigationStack(path: $path) {
            Group {
                if scrolls { ScrollView { feed } } else { feed.fixedSize(horizontal: false, vertical: true).frame(maxHeight: .infinity, alignment: .top).clipped() }
            }
            .navigationTitle("Discover")
            .searchable(text: $query, prompt: "Curves, tags, people")
            .navigationDestination(for: String.self) { id in
                if let item = items.first(where: { $0.id == id }) {
                    DiscoverDetail(item: item)
                        #if os(iOS)
                        .navigationTransition(.zoom(sourceID: id, in: zoom))
                        #endif
                }
            }
        }
    }

    private var chips: some View {
        HStack(spacing: Space.s) {
            ForEach(["All", "Periodic", "Polar", "Closed", "Symmetric", "Polynomial"], id: \.self) {
                FilterChip($0, selected: $0 == selectedTag)
            }
        }
    }

    private var feed: some View {
        VStack(alignment: .leading, spacing: Space.xl) {
            ScrollView(.horizontal, showsIndicators: false) { chips.padding(.horizontal, Space.l) }
                .scrollClipDisabled()
                .padding(.horizontal, -Space.l)

            VStack(alignment: .leading, spacing: Space.m) {
                SectionHeader("Curve of the Week")
                NavigationLink(value: items[0].id) { HeroCard(item: items[0]) }.buttonStyle(.plain)
                    .matchedTransitionSource(id: items[0].id, in: zoom)
            }

            VStack(alignment: .leading, spacing: Space.m) {
                SectionHeader("New")
                LazyVGrid(columns: [GridItem(.flexible(), spacing: Space.m), GridItem(.flexible(), spacing: Space.m)], spacing: Space.m) {
                    ForEach(items.dropFirst()) { item in
                        NavigationLink(value: item.id) { ThumbCard(item: item) }.buttonStyle(.plain)
                            .matchedTransitionSource(id: item.id, in: zoom)
                            .scrollTransition(.animated) { card, phase in
                                card.opacity(phase.isIdentity ? 1 : 0.7).scaleEffect(phase.isIdentity ? 1 : 0.97)
                            }
                    }
                }
            }
        }
        .padding(.horizontal, Space.l)
        .padding(.bottom, Space.xxl)
    }
}

struct HeroCard: View {
    let item: DiscoverItem
    var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            PlotView { ItemCurve(item: item, lineWidth: 3.5) }
                .frame(height: 210)
            HStack(alignment: .center) {
                VStack(alignment: .leading, spacing: Space.xs) {
                    FitMath(item.latex, size: Typeface.cardEquationSize, alignment: .leading)
                    Text([item.name, item.handle].compactMap { $0 }.joined(separator: " · "))
                        .font(Typeface.meta).foregroundStyle(.secondary)
                }
                Spacer()
                Button("Remix", systemImage: "arrow.trianglehead.branch") {}.secondaryButton()
            }
            .padding(.horizontal, Space.s).padding(.top, Space.m).padding(.bottom, Space.xs)
        }
        .padding(Radius.cardPadding)
        .cardSurface(radius: Radius.heroCard)
    }
}

struct ThumbCard: View {
    let item: DiscoverItem
    var body: some View {
        VStack(alignment: .leading, spacing: Space.s) {
            PlotView { ItemCurve(item: item, lineWidth: 2.5) }
                .frame(height: 118)
            FitMath(item.latex, size: Typeface.thumbEquationSize, alignment: .leading)
                .frame(minHeight: 28, alignment: .leading)
                .padding(.horizontal, Space.xs)
            HStack {
                Text(item.handle ?? "")
                Spacer()
                Label(item.likesText, systemImage: "heart")
            }
            .font(Typeface.meta).foregroundStyle(.secondary)
            .padding(.horizontal, Space.xs).padding(.bottom, Space.xs)
        }
        .padding(Radius.cardPadding)
        .cardSurface()
    }
}


public struct DiscoverDetail: View {
    let item: DiscoverItem
    var scrolls: Bool
    @State private var showStroke = false
    public init(item: DiscoverItem, scrolls: Bool = true) { self.item = item; self.scrolls = scrolls }

    public var body: some View {
        Group {
            if scrolls { ScrollView { page } } else { page.fixedSize(horizontal: false, vertical: true).frame(maxHeight: .infinity, alignment: .top).clipped() }
        }
        .navigationTitle(item.name ?? "")
        #if os(iOS)
        .navigationBarTitleDisplayMode(.inline)
        #endif
        .toolbar {
            ToolbarItemGroup(placement: .primaryAction) {
                Button("Like", systemImage: "heart") {}.tint(.primary)
                Button("Share", systemImage: "square.and.arrow.up") {}.tint(.primary)
            }
        }
    }

    private var page: some View {
        VStack(alignment: .leading, spacing: Space.xl) {
            VStack(alignment: .leading, spacing: 0) {
                PlotView { ItemCurve(item: item, lineWidth: 3.5, showStroke: showStroke) }
                    .frame(height: 260)
                    .contentShape(Rectangle())
                    .onTapGesture { if item.stroke.count > 1 { withAnimation(.easeInOut(duration: 0.2)) { showStroke.toggle() } } }
                VStack(alignment: .leading, spacing: Space.xs) {
                    FitMath(item.latex, size: 26)
                    HStack {
                        Text([item.handle.map { "drawn by \($0)" }, "\(item.likesText) likes"].compactMap { $0 }.joined(separator: " · "))
                        Spacer()
                        Text("Graphiti").fontWeight(.semibold).foregroundStyle(.primary)
                    }
                    .font(Typeface.meta).foregroundStyle(.secondary)
                }
                .padding(.horizontal, Space.s).padding(.top, Space.m).padding(.bottom, Space.xs)
            }
            .padding(Radius.cardPadding)
            .cardSurface(radius: Radius.heroCard)

            FlowLayout(spacing: Space.s) { ForEach(item.tags, id: \.self) { Tag($0) } }

            PrimaryButton("Open in Solve", systemImage: "function")

            VStack(alignment: .leading, spacing: Space.m) {
                HStack(alignment: .firstTextBaseline, spacing: Space.s) {
                    SectionHeader("Comments", chevron: false)
                    Text("14").font(Typeface.status).foregroundStyle(.secondary)
                }
                VStack(spacing: 0) {
                    comment("ana.k", "The 3:2 ratio is why it crosses itself four times.")
                    Divider().padding(.leading, 56)
                    comment("mrs.ortiz", "Using this in period 3 tomorrow.")
                    Divider().padding(.leading, 56)
                    comment("noor", "Changed 3 to 5 and got a flower")
                }
                .background(Theme.surface, in: .rect(cornerRadius: Radius.notice))
            }
        }
        .padding(.horizontal, Space.l)
        .padding(.bottom, Space.xxl)
    }

    private func comment(_ who: String, _ text: String) -> some View {
        HStack(alignment: .top, spacing: Space.m) {
            Image(systemName: "person.crop.circle.fill")
                .font(.system(size: 32)).foregroundStyle(.secondary)
            VStack(alignment: .leading, spacing: 2) {
                Text("@" + who).font(Typeface.meta.weight(.semibold))
                Text(text).font(Typeface.body)
            }
            Spacer(minLength: 0)
        }
        .padding(Space.m)
    }
}

struct FlowLayout: Layout {
    var spacing: CGFloat
    func sizeThatFits(proposal: ProposedViewSize, subviews: Subviews, cache: inout ()) -> CGSize {
        let rows = arrange(width: proposal.width ?? .infinity, subviews)
        return CGSize(width: proposal.width ?? rows.width, height: rows.height)
    }
    func placeSubviews(in bounds: CGRect, proposal: ProposedViewSize, subviews: Subviews, cache: inout ()) {
        let rows = arrange(width: bounds.width, subviews)
        for (i, p) in rows.points.enumerated() {
            subviews[i].place(at: CGPoint(x: bounds.minX + p.x, y: bounds.minY + p.y), proposal: .unspecified)
        }
    }
    private func arrange(width: CGFloat, _ subviews: Subviews) -> (points: [CGPoint], width: CGFloat, height: CGFloat) {
        var x: CGFloat = 0, y: CGFloat = 0, rowH: CGFloat = 0, maxW: CGFloat = 0
        var pts: [CGPoint] = []
        for s in subviews {
            let d = s.sizeThatFits(.unspecified)
            if x > 0 && x + d.width > width { x = 0; y += rowH + spacing; rowH = 0 }
            pts.append(CGPoint(x: x, y: y))
            x += d.width + spacing; rowH = max(rowH, d.height); maxW = max(maxW, x)
        }
        return (pts, maxW, y + rowH)
    }
}
