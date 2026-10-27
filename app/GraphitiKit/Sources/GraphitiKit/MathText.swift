import SwiftUI
import SwiftMath

extension VerticalAlignment {
    private enum MathBaseline: AlignmentID {
        static func defaultValue(in d: ViewDimensions) -> CGFloat { d[VerticalAlignment.center] }
    }
    public static let mathBaseline = VerticalAlignment(MathBaseline.self)
}

public struct MathText: View {
    let latex: String
    let size: CGFloat
    let color: PlatformColor

    public init(_ latex: String, size: CGFloat, color: PlatformColor = .white) {
        self.latex = latex
        self.size = size
        self.color = color
    }

    public var body: some View {
        if let r = MathRenderCache.shared.render(latex, size: size, color: color) {
            Image(platformImage: r.image)
                .alignmentGuide(.mathBaseline) { _ in r.baselineFromTop }
        } else {
            Text(latex).font(.system(size: size, design: .serif)).italic()
        }
    }
}

struct MathRender {
    let image: PlatformImage
    let baselineFromTop: CGFloat
}

final class MathRenderCache: @unchecked Sendable {
    static let shared = MathRenderCache()
    private var cache: [String: MathRender] = [:]
    private let lock = NSLock()

    func render(_ latex: String, size: CGFloat, color: PlatformColor) -> MathRender? {
        let key = "\(latex)|\(size)|\(color.description)"
        lock.lock(); defer { lock.unlock() }
        if let hit = cache[key] { return hit }
        if cache.count > 600 { cache.removeAll(keepingCapacity: true) }

        let pad = size * 0.3
        var math = MathImage(latex: latex, fontSize: size, textColor: color, labelMode: .text, textAlignment: .left)
        math.font = .latinModernFont
        math.contentInsets = MTEdgeInsets(top: pad, left: 0, bottom: pad, right: 0)
        let (error, image, info) = math.asImage()
        guard error == nil, let image, let info else { return nil }

        let h = max(info.ascent + info.descent, size / 2)
        let available = image.size.height - 2 * pad
        let baselineFromBottom = (available - h) / 2 + info.descent + pad
        let r = MathRender(image: image, baselineFromTop: image.size.height - baselineFromBottom)
        cache[key] = r
        return r
    }
}

public struct FitMath: View {
    let latex: String
    let size: CGFloat
    var color: PlatformColor = .white
    var alignment: Alignment = .center
    public init(_ latex: String, size: CGFloat, color: PlatformColor = .white, alignment: Alignment = .center) {
        self.latex = latex; self.size = size; self.color = color; self.alignment = alignment
    }
    static func lines(_ latex: String) -> [String] {
        for sep in [",\\quad ", ",\\;"] where latex.contains(sep) { return latex.components(separatedBy: sep) }
        return [latex]
    }
    public var body: some View {
        let ls = Self.lines(latex)
        VStack(alignment: alignment.horizontal, spacing: 2) {
            ForEach(Array(ls.enumerated()), id: \.offset) { _, l in
                MathText(l, size: ls.count > 1 ? size * 0.85 : size, color: color).fixedSize().scaleToFitWidth(alignment: alignment)
                    .frame(height: (ls.count > 1 ? size * 0.85 : size) * 1.5)
            }
        }
    }
}
