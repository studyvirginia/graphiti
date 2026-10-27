import Foundation
import SwiftUI
#if os(iOS)
import MetricKit
#endif

public final class MetricsSink: NSObject, MXMetricManagerSubscriber {
    public static let shared = MetricsSink()
    public func start() { MXMetricManager.shared.add(self) }
    public func didReceive(_ payloads: [MXMetricPayload]) {
        for p in payloads { save(p.jsonRepresentation(), kind: "metric") }
    }
    public func didReceive(_ payloads: [MXDiagnosticPayload]) {
        for p in payloads { save(p.jsonRepresentation(), kind: "diagnostic") }
    }
    private func save(_ data: Data, kind: String) {
        let stamp = Int(Date().timeIntervalSince1970)
        let name = "\(kind == "diagnostic" ? "metric-diag" : "metric")-\(stamp).json"
        let url = Diag.dir.appendingPathComponent(name)
        try? data.write(to: url, options: .atomic)
        Diag.log("metrickit", "saved \(name) (\(data.count) bytes)", level: .error)
    }
}

public struct DiagnosticsView: View {
    @AppStorage("diagVerboseEngine") private var verbose = true
    @State private var files: [(url: URL, bytes: Int)] = []
    @State private var confirmClear = false

    public init() {}

    public var body: some View {
        List {
            Section {
                Toggle("Verbose engine logging", isOn: $verbose)
            } footer: { Text("Writes the engine's stdout and stderr into engine-<session>.log. Takes effect at next launch.") }
            Section("Log files") {
                if files.isEmpty { Text("No log files yet.").foregroundStyle(.secondary) }
                ForEach(files, id: \.url) { f in
                    NavigationLink {
                        LogTailView(url: f.url)
                    } label: {
                        HStack {
                            Text(f.url.lastPathComponent).font(.system(.footnote, design: .monospaced))
                            Spacer()
                            Text(ByteCountFormatter.string(fromByteCount: Int64(f.bytes), countStyle: .file))
                                .foregroundStyle(.secondary).font(.footnote)
                        }
                    }
                }
            }
            Section {
                ShareLink(items: files.map(\.url)) { Label("Share all log files", systemImage: "square.and.arrow.up") }
                    .disabled(files.isEmpty)
                Button("Clear logs", role: .destructive) { confirmClear = true }
            } footer: { Text("Logs are also visible in the Files app under GraphArt.") }
        }
        .navigationTitle("Diagnostics")
        .onAppear(perform: reload)
        .refreshable { reload() }
        .confirmationDialog("Delete all log files except this session's?", isPresented: $confirmClear, titleVisibility: .visible) {
            Button("Clear", role: .destructive) { Diag.clear(); reload() }
        }
    }

    private func reload() { files = Diag.files() }
}

struct LogTailView: View {
    let url: URL
    @State private var text = ""

    var body: some View {
        ScrollView {
            Text(text)
                .font(.system(.caption2, design: .monospaced))
                .textSelection(.enabled)
                .frame(maxWidth: .infinity, alignment: .leading)
                .padding()
        }
        .navigationTitle(url.lastPathComponent)
        .toolbar { ShareLink(item: url) }
        .onAppear { text = Diag.tail(url, lines: 500) }
    }
}
