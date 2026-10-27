import Foundation
import os
#if os(iOS)
import UIKit
#endif

public enum Diag {
    static let osLog = Logger(subsystem: "com.graphart.app", category: "diag")
    public static let session: String = {
        let f = DateFormatter(); f.dateFormat = "yyyyMMdd-HHmmss"; f.locale = Locale(identifier: "en_US_POSIX")
        return f.string(from: Date())
    }()
    public static var verboseEngine: Bool { UserDefaults.standard.object(forKey: "diagVerboseEngine") as? Bool ?? true }

    public static var dir: URL {
        let d = FileManager.default.urls(for: .documentDirectory, in: .userDomainMask)[0].appendingPathComponent("logs", isDirectory: true)
        try? FileManager.default.createDirectory(at: d, withIntermediateDirectories: true)
        return d
    }
    public static var appLogURL: URL { dir.appendingPathComponent("app-\(session).log") }
    private static let queue = DispatchQueue(label: "graphart.diag")
    private static let stamp: DateFormatter = {
        let f = DateFormatter(); f.dateFormat = "HH:mm:ss.SSS"; f.locale = Locale(identifier: "en_US_POSIX"); return f
    }()

    public static func log(_ category: String, _ message: String, level: OSLogType = .default) {
        osLog.log(level: level, "[\(category, privacy: .public)] \(message, privacy: .public)")
        let line = "\(stamp.string(from: Date())) [\(category)] \(message)\n"
        queue.async { append(line, to: appLogURL) }
    }

    static func append(_ line: String, to url: URL) {
        guard let data = line.data(using: .utf8) else { return }
        if !FileManager.default.fileExists(atPath: url.path) { FileManager.default.createFile(atPath: url.path, contents: nil) }
        guard let h = try? FileHandle(forWritingTo: url) else { return }
        defer { try? h.close() }
        _ = try? h.seekToEnd()
        try? h.write(contentsOf: data)
    }

    public static func tail(_ url: URL, lines count: Int = 500) -> String {
        guard let s = try? String(contentsOf: url, encoding: .utf8) else { return "(unreadable)" }
        let all = s.split(separator: "\n", omittingEmptySubsequences: false)
        return all.suffix(count).joined(separator: "\n")
    }

    public static func files() -> [(url: URL, bytes: Int)] {
        let urls = (try? FileManager.default.contentsOfDirectory(at: dir, includingPropertiesForKeys: [.contentModificationDateKey, .fileSizeKey])) ?? []
        return urls.filter { $0.pathExtension == "log" || $0.pathExtension == "json" }
            .map { u -> (url: URL, date: Date, bytes: Int) in
                let v = try? u.resourceValues(forKeys: [.contentModificationDateKey, .fileSizeKey])
                return (u, v?.contentModificationDate ?? .distantPast, v?.fileSize ?? 0)
            }
            .sorted { $0.date > $1.date }
            .map { ($0.url, $0.bytes) }
    }

    public static func clear() {
        queue.sync {
            for u in (try? FileManager.default.contentsOfDirectory(at: dir, includingPropertiesForKeys: nil)) ?? []
                where u.lastPathComponent != "app-\(session).log" && u.lastPathComponent != "engine-\(session).log" {
                try? FileManager.default.removeItem(at: u)
            }
        }
    }

    public static func lastSessionCrash() -> String? {
        let urls = (try? FileManager.default.contentsOfDirectory(at: dir, includingPropertiesForKeys: nil)) ?? []
        return urls.map(\.lastPathComponent).filter { $0.hasPrefix("crash-") || $0.hasPrefix("metric-") }.sorted().last
    }

    public static func start() {
        #if os(iOS)
        let dev = UIDevice.current
        let bundle = Bundle.main
        let v = bundle.object(forInfoDictionaryKey: "CFBundleShortVersionString") as? String ?? "?"
        log("launch", "GraphArt \(v) session \(session) on \(dev.model) iOS \(dev.systemVersion) verbose=\(verboseEngine)")
        #else
        log("launch", "session \(session)")
        #endif
        if let c = lastSessionCrash() { log("launch", "LAST SESSION CRASHED: \(c)", level: .error) }
    }

    nonisolated(unsafe) private static var started = false
    public static func startLogging() {
        guard !started else { return }
        started = true
        dir.path.withCString { d in session.withCString { s in graphiti_diag_init(d, s, verboseEngine ? 1 : 0) } }
        start()
        observe()
        #if os(iOS)
        MetricsSink.shared.start()
        #endif
    }

    private static var observers: [NSObjectProtocol] = []
    private static func observe() {
        #if os(iOS)
        let nc = NotificationCenter.default
        observers.append(nc.addObserver(forName: UIApplication.didReceiveMemoryWarningNotification, object: nil, queue: nil) { _ in
            log("memory", "MEMORY WARNING available=\(availableMB()) MB", level: .error)
        })
        observers.append(nc.addObserver(forName: UIApplication.didBecomeActiveNotification, object: nil, queue: nil) { _ in
            log("scene", "active (available=\(availableMB()) MB)")
        })
        observers.append(nc.addObserver(forName: UIApplication.willResignActiveNotification, object: nil, queue: nil) { _ in
            log("scene", "resign active (available=\(availableMB()) MB)")
        })
        observers.append(nc.addObserver(forName: UIApplication.didEnterBackgroundNotification, object: nil, queue: nil) { _ in
            log("scene", "background (available=\(availableMB()) MB)")
        })
        #endif
    }

    public static func availableMB() -> Int {
        #if os(iOS)
        return Int(os_proc_available_memory() / (1 << 20))
        #else
        return -1
        #endif
    }
}
