import SwiftUI

#if os(iOS)
import UIKit
public typealias PlatformColor = UIColor
public typealias PlatformImage = UIImage
#else
import AppKit
public typealias PlatformColor = NSColor
public typealias PlatformImage = NSImage
#endif

public enum Theme {
    public static let ink = Color.black
    public static let blueC = Color(hex: 0x58C4DD)
    public static let blueE = Color(hex: 0x1C758A)
    public static let goldC = Color(hex: 0xF0AC5F)
    public static let redC = Color(hex: 0xFC6255)
    public static let tealC = Color(hex: 0x5CD0B3)
    public static let chalk = Color.white
    public static let secondary = Color.white.opacity(0.6)
    public static let tertiary = Color.white.opacity(0.32)
    public static let fill = Color(white: 0.46).opacity(0.24)
    #if os(iOS)
    public static let surface = Color(uiColor: .secondarySystemBackground)
    public static let surfaceRaised = Color(uiColor: .tertiarySystemBackground)
    public static let hairline = Color(uiColor: .separator)
    #else
    public static let surface = Color(hex: 0x1C1C1E)
    public static let surfaceRaised = Color(hex: 0x2C2C2E)
    public static let hairline = Color.white.opacity(0.15)
    #endif
    public static let accent = blueC

}

extension Color {
    init(hex: UInt32, opacity: Double = 1) {
        self.init(.sRGB,
                  red: Double((hex >> 16) & 0xFF) / 255,
                  green: Double((hex >> 8) & 0xFF) / 255,
                  blue: Double(hex & 0xFF) / 255,
                  opacity: opacity)
    }
}

public extension PlatformColor {
    static func graphiti(hex: UInt32, alpha: CGFloat = 1) -> PlatformColor {
        PlatformColor(red: CGFloat((hex >> 16) & 0xFF) / 255,
                      green: CGFloat((hex >> 8) & 0xFF) / 255,
                      blue: CGFloat(hex & 0xFF) / 255,
                      alpha: alpha)
    }
}

extension Image {
    init(platformImage: PlatformImage) {
        #if os(iOS)
        self.init(uiImage: platformImage)
        #else
        self.init(nsImage: platformImage)
        #endif
    }
}
