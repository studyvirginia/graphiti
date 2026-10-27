import SwiftUI

public struct GraphitiRoot: View {
    public enum Demo: String, Sendable { case live, firstLaunch, empty, locked, explain, shareBack, discover, detail, share, language, info, howItWorks, you }

    @StateObject private var session: SolveSession
    @State private var tab: Int
    @State private var sharing: Bool
    @State private var showingInfo: Bool
    @State private var showingHowItWorks: Bool
    @AppStorage("hasSeenHowItWorks") private var hasSeenHowItWorks = false
    @AppStorage("searchSeconds") private var searchSeconds = SearchBudget.defaultSeconds
    @AppStorage("displayName") private var displayName = ""
    @ObservedObject private var library = Library.shared
    @State private var savedID: String?
    private let demo: Demo

    public init(demo: Demo = .live) {
        self.demo = demo
        let s: SolveSession
        switch demo {
        case .locked, .share, .shareBack, .info: s = SolveSession(step: 8, locked: true)
        case .explain: s = SolveSession(step: 8, locked: true, explained: "wave")
        default: s = SolveSession(phase: .empty)
        }
        _session = StateObject(wrappedValue: s)
        _tab = State(initialValue: [.discover: 1, .detail: 1, .you: 2][demo] ?? 0)
        _sharing = State(initialValue: demo == .share || demo == .shareBack)
        _showingInfo = State(initialValue: demo == .info)
        _showingHowItWorks = State(initialValue: demo == .howItWorks)
    }

    private var card: MathCard {
        MathCard(params: session.target, latex: session.latex, author: displayName, points: session.liveCurvePoints, stroke: session.liveStroke)
    }

    public var body: some View {
        if demo == .language { DesignLanguageView() } else {
            tabs
        }
    }

    private var tabs: some View {
        TabView(selection: $tab) {
            Tab("Solve", systemImage: "function", value: 0) {
                NavigationStack {
                    SolveScreen(session: session,
                                onShare: { show(true) },
                                onInfo: { showingInfo = true },
                                onStroke: { pts in Task {
                                    hasSeenHowItWorks = true
                                    if demo == .live, pts.count >= 4 {
                                        await session.solveEngine(pts, budget: SearchBudget(seconds: searchSeconds).seconds)
                                    } else {
                                        await session.play()
                                    }
                                    try? await Task.sleep(for: .milliseconds(650))
                                    if session.locked {
                                        savedID = demo == .live ? saveSolve() : nil
                                        show(true)
                                    }
                                } },
                                firstLaunchHint: (demo == .live && !hasSeenHowItWorks) || demo == .firstLaunch ? 0 : nil)
                        .toolbar { solveToolbar }
                        .tint(.primary)
                }
            }
            Tab("Discover", systemImage: "safari", value: 1) {
                DiscoverScreen(initialPath: demo == .detail ? ["lissajous"] : [])
            }
            Tab("You", systemImage: "person.crop.circle", value: 2) { YouScreen() }
        }
        #if os(iOS)
        .tabBarMinimizeBehavior(.onScrollDown)
        #endif
        .overlay {
            if sharing {
                SharePopup(card: card,
                           back: MathCardBack(stats: session.stats, latex: session.latex, author: displayName),
                           startFlipped: demo == .shareBack,
                           onClose: { show(false) },
                           onSave: { savedID = saveSolve() },
                           onPost: { if let id = savedID ?? saveSolve() { savedID = id; library.post(id) } },
                           posted: library.isPosted(savedID)) { Color.clear }
                    .background(.ultraThinMaterial)
                    .contentShape(Rectangle())
                    .onTapGesture { show(false) }
                    .transition(.opacity.combined(with: .scale(scale: 1.04)))
                    .zIndex(1)
            }
        }
        .sheet(isPresented: $showingInfo) {
            EquationDetailsSheet(params: session.target, latex: session.latex, seconds: session.elapsed,
                                 live: session.liveStroke.isEmpty ? nil : (session.stats, session.liveCurvePoints, session.liveStroke))
                .presentationDetents([.medium, .large])
        }
        .sheet(isPresented: $showingHowItWorks) {
            HowItWorksSheet().presentationDetents([.medium])
        }
        .tint(accent)
        .environment(\.screenAccent, accent)
        .preferredColorScheme(.dark)
    }

    private var accent: Color { tab == 0 ? session.colour : Theme.accent }

    private func saveSolve() -> String? {
        guard !session.liveStroke.isEmpty else { return nil }
        let st = session.stats
        return library.save(latex: session.latex, points: session.liveCurvePoints, stroke: session.liveStroke, seconds: st.seconds,
                            tried: st.equationsTried, author: displayName)
    }

    private func show(_ on: Bool) {
        withAnimation(.spring(response: 0.42, dampingFraction: 0.86)) { sharing = on }
    }

    @ToolbarContentBuilder private var solveToolbar: some ToolbarContent {
        #if os(iOS)
        if !session.isEmpty {
            ToolbarItem(placement: .topBarLeading) {
                IconButton("arrow.uturn.backward", "Undo") { session.undo() }
            }
            .sharedBackgroundVisibility(.hidden)
            ToolbarItem(placement: .topBarTrailing) {
                IconButton("square.and.pencil", "New") { session.new() }
            }
            .sharedBackgroundVisibility(.hidden)
        }
        #else
        ToolbarItem(placement: .navigation) {
            Button("Undo", systemImage: "arrow.uturn.backward") { session.undo() }
        }
        ToolbarItem(placement: .primaryAction) {
            Button("New", systemImage: "square.and.pencil") { session.new() }
        }
        #endif
    }
}
