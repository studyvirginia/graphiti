import SwiftUI


public struct EquationDetailsSheet: View {
    let params: CurveParams
    let latex: String
    let seconds: Double
    let live: (stats: SolveSession.Stats, points: [CGPoint], stroke: [CGPoint])?
    @State private var explanation: String?
    @State private var explainError: String?
    @State private var showStroke = false
    @Environment(\.dismiss) private var dismiss

    public init(params: CurveParams, latex: String, seconds: Double, live: (SolveSession.Stats, [CGPoint], [CGPoint])? = nil) {
        self.params = params; self.latex = latex; self.seconds = seconds
        self.live = live.map { (stats: $0.0, points: $0.1, stroke: $0.2) }
    }

    public var body: some View {
        NavigationStack {
            ScrollView {
                VStack(alignment: .leading, spacing: Space.xl) {
                    VStack(spacing: Space.m) {
                        PlotView {
                            if let live {
                                let box = live.points + live.stroke
                                if showStroke && live.stroke.count > 1 {
                                    FittedPolyline(points: live.stroke, frame: box).stroke(.white.opacity(0.35), style: StrokeStyle(lineWidth: 6, lineCap: .round, lineJoin: .round))
                                }
                                FittedPolyline(points: live.points, frame: box).stroke(EquationColor.color(for: latex), lineWidth: 3)
                            }
                            else { CurveShape(params: params, layout: .card).stroke(EquationColor.color(for: latex), lineWidth: 3) }
                        }
                            .frame(height: 170)
                            .contentShape(Rectangle())
                            .onTapGesture { if let live, live.stroke.count > 1 { withAnimation(.easeInOut(duration: 0.2)) { showStroke.toggle() } } }
                        FitMath(latex, size: 28)
                    }
                    .padding(Radius.cardPadding)
                    .cardSurface(radius: Radius.heroCard)

                    if let live { liveSections(live.stats) } else {
                    VStack(alignment: .leading, spacing: Space.m) {
                        SectionHeader("What each part does", chevron: false)
                        VStack(spacing: 0) {
                            term("\\frac{1}{2}x", EquationColor.color(for: latex), EquationColor.platformColor(for: latex),
                                 "The tilt", "The curve climbs half a unit for every unit it moves right.")
                            Divider().padding(.leading, Space.m)
                            term("\\sin(2x)", EquationColor.color(for: latex), EquationColor.platformColor(for: latex),
                                 "The wave", "Rises 1 above the tilt and 1 below, and repeats every π ≈ 3.14 units.")
                        }
                        .background(Theme.surface, in: .rect(cornerRadius: Radius.notice))
                    }

                    VStack(alignment: .leading, spacing: Space.m) {
                        SectionHeader("How it was found", chevron: false)
                        HStack(spacing: Space.m) {
                            stat(String(format: "%.1f s", seconds), "to find")
                            stat("31", "generations")
                            stat("99.6%", "fit")
                        }
                        Text("Candidates competed to match your stroke. Each generation kept the best and tried small changes: a new term, a different constant, a term removed. Then ½ was recognised as exact.")
                            .font(Typeface.meta).foregroundStyle(.secondary)
                    }
                    }

                    if live == nil {
                        Button {} label: { Label("Hear the Curve", systemImage: "waveform").frame(maxWidth: .infinity) }
                            .secondaryButton().controlSize(.large)
                    }
                }
                .padding(.horizontal, Space.l)
                .padding(.bottom, Space.xl)
            }
            .navigationTitle("This Equation")
            #if os(iOS)
            .navigationBarTitleDisplayMode(.inline)
            #endif
            .toolbar {
                ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } }
            }
        }
    }

    @ViewBuilder private func liveSections(_ st: SolveSession.Stats) -> some View {
        VStack(alignment: .leading, spacing: Space.m) {
            SectionHeader("What it means", chevron: false)
            Group {
                if let explanation { Text(explanation) }
                else if let explainError { Text(explainError).foregroundStyle(.secondary) }
                else { HStack(spacing: Space.s) { ProgressView(); Text("Asking Claude…").foregroundStyle(.secondary) } }
            }
            .font(Typeface.meta)
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(Space.m)
            .background(Theme.surface, in: .rect(cornerRadius: Radius.notice))
            .task(id: latex) {
                explanation = nil; explainError = nil
                do { explanation = try await EquationExplainer.explain(latex: latex, stats: st) }
                catch { explainError = error.localizedDescription }
            }
        }
        VStack(alignment: .leading, spacing: Space.m) {
            SectionHeader("How it was found", chevron: false)
            HStack(spacing: Space.m) {
                stat(String(format: "%.1f s", st.seconds), "to find")
                stat(String(SolveSession.triedText(st.equationsTried).dropFirst(3).dropLast(6)), "tried")
                stat(String(format: "%.0f%%", st.fit * 100), "fit")
            }
            Text("\(st.generations) champions shown, \(st.terms) terms. Candidate equations competed on your stroke, scored by how likely your drawing is under each, and the simplest one within the noise of the best won.")
                .font(Typeface.meta).foregroundStyle(.secondary)
        }
    }

    private func term(_ tex: String, _ colour: Color, _ platform: PlatformColor, _ name: String, _ detail: String) -> some View {
        HStack(alignment: .center, spacing: Space.m) {
            MathText(tex, size: 22, color: platform)
                .fixedSize()
                .frame(minWidth: 72, alignment: .leading)
            VStack(alignment: .leading, spacing: 2) {
                Text(name).font(Typeface.headline)
                Text(detail).font(Typeface.meta).foregroundStyle(.secondary)
            }
            Spacer(minLength: 0)
        }
        .padding(Space.m)
        .accessibilityElement(children: .combine)
    }

    private func stat(_ value: String, _ label: String) -> some View {
        VStack(spacing: 2) {
            Text(value).font(.title3.weight(.semibold)).monospacedDigit()
            Text(label).font(Typeface.label).foregroundStyle(.secondary)
        }
        .frame(maxWidth: .infinity).padding(.vertical, Space.m)
        .background(Theme.surface, in: .rect(cornerRadius: Radius.notice))
    }
}


public struct HowItWorksSheet: View {
    @Environment(\.dismiss) private var dismiss
    @State private var page = 0
    public init() {}

    private let pages: [(String, String, String)] = [
        ("scribble", "Draw One Stroke", "Any curve: a wave, a hill, a spiral. One finger, one stroke."),
        ("arrow.triangle.branch", "Watch It Evolve", "Equations compete to match it. The best one changes term by term until it locks."),
        ("hand.tap", "Tap the Math", "Tap a term to see what it does. Hold the equation to copy it."),
    ]

    public var body: some View {
        VStack(spacing: Space.xl) {
            TabView(selection: $page) {
                ForEach(pages.indices, id: \.self) { i in
                    VStack(spacing: Space.l) {
                        Image(systemName: pages[i].0)
                            .font(.system(size: 64, weight: .light))
                            .foregroundStyle(Theme.accent)
                            .symbolEffect(.bounce, value: page == i)
                            .frame(height: 90)
                        Text(pages[i].1).font(.title.weight(.bold))
                        Text(pages[i].2).font(Typeface.body).foregroundStyle(.secondary)
                            .multilineTextAlignment(.center).padding(.horizontal, Space.xxl)
                    }
                    .tag(i)
                }
            }
            #if os(iOS)
            .tabViewStyle(.page(indexDisplayMode: .always))
            #endif
            PrimaryButton(page == pages.count - 1 ? "Start Drawing" : "Next") {
                if page < pages.count - 1 { withAnimation { page += 1 } } else { dismiss() }
            }
            .padding(.horizontal, Space.l)
        }
        .padding(.top, Space.xxl).padding(.bottom, Space.l)
    }
}


public struct YouScreen: View {
    @ObservedObject private var library = Library.shared
    @State private var showingSettings = false
    private var items: [DiscoverItem] { library.entries.map(\.item) }
    public init() {}

    public var body: some View {
        NavigationStack {
            ScrollView {
                VStack(alignment: .leading, spacing: Space.xl) {
                    HStack(spacing: Space.m) {
                        stat("\(library.entries.count)", "Solved")
                        stat("\(library.posted.count)", "Posted")
                        stat("0", "Likes")
                    }
                    VStack(alignment: .leading, spacing: Space.m) {
                        SectionHeader("Saved")
                        if items.isEmpty {
                            Text("Curves you solve are saved here automatically.").font(Typeface.meta).foregroundStyle(.secondary)
                        }
                        LazyVGrid(columns: [GridItem(.flexible(), spacing: Space.m), GridItem(.flexible(), spacing: Space.m)], spacing: Space.m) {
                            ForEach(items) { ThumbCard(item: $0) }
                        }
                    }
                }
                .padding(.horizontal, Space.l)
                .padding(.bottom, Space.xxl)
            }
            .navigationTitle("You")
            .toolbar {
                ToolbarItem(placement: .primaryAction) {
                    Button("Settings", systemImage: "gearshape") { showingSettings = true }.tint(.primary)
                }
            }
            .sheet(isPresented: $showingSettings) { SettingsSheet().presentationDetents([.medium]) }
        }
    }

    private func stat(_ value: String, _ label: String) -> some View {
        VStack(spacing: 2) {
            Text(value).font(.title2.weight(.bold)).monospacedDigit()
            Text(label).font(Typeface.label).foregroundStyle(.secondary)
        }
        .frame(maxWidth: .infinity).padding(.vertical, Space.m)
        .background(Theme.surface, in: .rect(cornerRadius: Radius.notice))
    }
}

// MARK: - Settings

public struct SettingsSheet: View {
    @AppStorage("soundsOn") private var soundsOn = true
    @AppStorage("displayName") private var name = ""
    @AppStorage("searchSeconds") private var searchSeconds = SearchBudget.defaultSeconds
    @AppStorage("searchGenerations") private var searchGenerations = 0
    @AppStorage(ConvergeStop.defaultsKey) private var convergeStop = false
    @AppStorage("noiseModelMode") private var noiseModelMode = 0
    @State private var showingHowItWorks = false
    @State private var apiKey = APIKeyStore.key ?? ""
    @Environment(\.dismiss) private var dismiss
    public init() {}

    public var body: some View {
        NavigationStack {
            Form {
                Section {
                    TextField("Name", text: $name)
                } header: { Text("Shown on your math cards") }
                Section {
                    Stepper("Search time: \(Int(searchSeconds)) s", value: $searchSeconds,
                            in: SearchBudget.range.lowerBound...SearchBudget.range.upperBound, step: 1)
                    Stepper("Generations (advanced): \(searchGenerations == 0 ? "off" : "\(searchGenerations)")",
                            value: $searchGenerations, in: 0...500, step: 5)
                    Toggle("Stop when converged", isOn: $convergeStop)
                    Picker("Noise model", selection: $noiseModelMode) {
                        Text("As shipped").tag(0); Text("Recalibrated (0.1×)").tag(1); Text("Off").tag(2)
                    }
                } header: { Text("Search") } footer: { Text("The engine stops after the step in which the search time is reached, then shows its final pick. A generation count (if set) overrides the time. Stop when converged ends the search early once the best equation has stopped improving for \(ConvergeStop.patience) steps (after at least \(ConvergeStop.minSteps)). Noise model: how much hand shake the scoring expects. Recalibrated applies on the next drawing; Off (plain fit error) applies after you close and reopen the app.") }
                Section {
                    SecureField("Anthropic API key", text: $apiKey)
                        .textContentType(.password).autocorrectionDisabled()
                        #if os(iOS)
                        .textInputAutocapitalization(.never)
                        #endif
                        .onSubmit { APIKeyStore.set(apiKey) }
                        .onChange(of: apiKey) { _, v in APIKeyStore.set(v) }
                } header: { Text("Explanations") } footer: { Text("Your own key, kept in this iPhone's Keychain. This Equation sends the equation and its stats to Claude Haiku for a plain-words explanation.") }
                Section {
                    Toggle("Sounds", isOn: $soundsOn)
                } footer: { Text("Quiet sounds for new terms, pruning and the lock. They follow the silent switch; haptics stay on.") }
                Section {
                    Button("How It Works") { showingHowItWorks = true }
                    NavigationLink("Diagnostics") { DiagnosticsView() }
                }
            }
            .navigationTitle("Settings")
            #if os(iOS)
            .navigationBarTitleDisplayMode(.inline)
            #endif
            .toolbar { ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } } }
            .sheet(isPresented: $showingHowItWorks) { HowItWorksSheet().presentationDetents([.medium]) }
        }
        .tint(Theme.accent)
    }
}
