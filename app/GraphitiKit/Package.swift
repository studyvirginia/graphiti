import PackageDescription

let package = Package(
    name: "GraphitiKit",
    platforms: [.iOS("26.1"), .macOS(.v26)],
    products: [
        .library(name: "GraphitiKit", targets: ["GraphitiKit"]),
        .executable(name: "render-screens", targets: ["RenderScreens"]),
    ],
    dependencies: [
        .package(url: "https://github.com/mgriebling/SwiftMath.git", from: "1.7.3"),
    ],
    targets: [
        .target(name: "GraphitiKit", dependencies: ["SwiftMath"]),
        .executableTarget(name: "RenderScreens", dependencies: ["GraphitiKit"]),
        .testTarget(name: "GraphitiKitTests", dependencies: ["GraphitiKit"]),
    ],
    swiftLanguageModes: [.v5]
)
