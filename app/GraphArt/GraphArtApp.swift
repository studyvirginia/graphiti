import SwiftUI
import GraphitiKit

@main
struct GraphArtApp: App {
    private let demo = UserDefaults.standard.string(forKey: "demo").flatMap(GraphitiRoot.Demo.init(rawValue:))

    init() {
        Diag.startLogging()
    }

    var body: some Scene {
        WindowGroup {
            GraphitiRoot(demo: demo ?? .live)
        }
        #if os(macOS)
        .windowStyle(.hiddenTitleBar)
        .defaultSize(width: 1280, height: 850)
        #endif
    }
}
