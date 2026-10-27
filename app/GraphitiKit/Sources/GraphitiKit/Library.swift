import SwiftUI

@MainActor
public final class Library: ObservableObject {
    public static let shared = Library()

    public struct Entry: Codable, Identifiable, Sendable, Equatable {
        public var id: String
        public var latex: String
        public var xs: [Double]
        public var ys: [Double]
        public var date: Date
        public var seconds: Double
        public var tried: Int
        public var author: String
        public var posted: Bool
        public var sx: [Double]? = nil
        public var sy: [Double]? = nil
        var points: [CGPoint] { zip(xs, ys).map { CGPoint(x: $0, y: $1) } }
        var stroke: [CGPoint] { zip(sx ?? [], sy ?? []).map { CGPoint(x: $0, y: $1) } }
    }

    @Published public private(set) var entries: [Entry] = []
    private let url: URL = FileManager.default.urls(for: .documentDirectory, in: .userDomainMask)[0]
        .appendingPathComponent("library.json")

    private init() {
        if let d = try? Data(contentsOf: url) {
            do { entries = try JSONDecoder().decode([Entry].self, from: d) }
            catch { Diag.log("library", "JSON decode failed for library.json: \(error)", level: .error) }
        }
    }

    @discardableResult
    public func save(latex: String, points: [CGPoint], stroke: [CGPoint] = [], seconds: Double, tried: Int, author: String) -> String? {
        guard points.count > 1 else { return nil }
        let xs = points.map { Double($0.x) }, ys = points.map { Double($0.y) }
        if let e = entries.first(where: { $0.latex == latex && $0.xs == xs }) { return e.id }
        let e = Entry(id: UUID().uuidString, latex: latex, xs: xs, ys: ys, date: Date(), seconds: seconds,
                      tried: tried, author: author, posted: false,
                      sx: stroke.map { Double($0.x) }, sy: stroke.map { Double($0.y) })
        entries.insert(e, at: 0)
        write()
        return e.id
    }

    public func post(_ id: String) {
        guard let i = entries.firstIndex(where: { $0.id == id }) else { return }
        entries[i].posted = true
        write()
    }

    public func isPosted(_ id: String?) -> Bool { entries.first { $0.id == id }?.posted ?? false }

    public var posted: [Entry] { entries.filter(\.posted) }

    private func write() {
        do {
            try JSONEncoder().encode(entries).write(to: url, options: .atomic)
        } catch {
            Diag.log("library", "save failed: \(error)", level: .error)
        }
    }
}

extension Library.Entry {
    var item: DiscoverItem {
        DiscoverItem(id, latex, .points(points), author: author == "you" ? "" : author, likes: 0, tags: [], stroke: stroke)
    }
}
