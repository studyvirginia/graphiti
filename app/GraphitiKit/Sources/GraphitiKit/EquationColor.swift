import SwiftUI

public enum EquationColor {
    public static var lightness = 0.80
    public static var chroma = 0.115

    static let pinnedHue = UserDefaults.standard.string(forKey: "hue").flatMap(Double.init)

    public static func hue(for latex: String) -> Double {
        if let pinnedHue { return pinnedHue }
        var h: UInt64 = 0xcbf29ce484222325
        for b in canonical(latex).utf8 { h ^= UInt64(b); h = h &* 0x100000001b3 }
        return Double(h % 3600) / 10
    }

    public static func color(for latex: String) -> Color {
        let (r, g, b) = rgb(hue: hue(for: latex))
        return Color(.sRGB, red: r, green: g, blue: b, opacity: 1)
    }

    public static func color(hue: Double) -> Color {
        let (r, g, b) = rgb(hue: hue)
        return Color(.sRGB, red: r, green: g, blue: b, opacity: 1)
    }

    public static func platformColor(for latex: String) -> PlatformColor {
        let (r, g, b) = rgb(hue: hue(for: latex))
        return PlatformColor(red: r, green: g, blue: b, alpha: 1)
    }

    static func canonical(_ latex: String) -> String {
        latex.filter { !$0.isWhitespace && $0 != "{" && $0 != "}" }
    }

    static func rgb(hue: Double) -> (Double, Double, Double) {
        let hr = hue * .pi / 180
        let L = lightness, a = chroma * cos(hr), b = chroma * sin(hr)
        let l_ = L + 0.3963377774 * a + 0.2158037573 * b
        let m_ = L - 0.1055613458 * a - 0.0638541728 * b
        let s_ = L - 0.0894841775 * a - 1.2914855480 * b
        let l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_
        let lin = (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
                   -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
                   -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)
        func enc(_ x: Double) -> Double {
            let c = min(max(x, 0), 1)
            return c <= 0.0031308 ? 12.92 * c : 1.055 * pow(c, 1 / 2.4) - 0.055
        }
        return (enc(lin.0), enc(lin.1), enc(lin.2))
    }
}
