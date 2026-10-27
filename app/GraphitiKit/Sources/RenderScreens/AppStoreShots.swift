import SwiftUI
import AppKit
import GraphitiKit


struct AppStoreShot: View {
    let headline: String
    let sub: String
    let image: NSImage
    var accent: Color

    var body: some View {
        ZStack(alignment: .top) {
            Color.black
            Canvas { ctx, size in
                var p = Path()
                stride(from: 0.0, through: size.width, by: 43).forEach { p.move(to: CGPoint(x: $0, y: 0)); p.addLine(to: CGPoint(x: $0, y: size.height)) }
                stride(from: 0.0, through: size.height, by: 43).forEach { p.move(to: CGPoint(x: 0, y: $0)); p.addLine(to: CGPoint(x: size.width, y: $0)) }
                ctx.stroke(p, with: .color(.white.opacity(0.06)), lineWidth: 1)
            }
            VStack(spacing: 0) {
                VStack(spacing: 10) {
                    Text(headline)
                        .font(.system(size: 42, weight: .bold))
                        .multilineTextAlignment(.center)
                        .foregroundStyle(.white)
                        .fixedSize(horizontal: false, vertical: true)
                    Text(sub)
                        .font(.system(size: 21, weight: .regular))
                        .multilineTextAlignment(.center)
                        .foregroundStyle(.white.opacity(0.62))
                        .fixedSize(horizontal: false, vertical: true)
                }
                .padding(.horizontal, 30)
                .frame(width: 430, height: 236)
                Image(nsImage: image)
                    .resizable()
                    .aspectRatio(contentMode: .fit)
                    .frame(width: 340)
                    .clipShape(RoundedRectangle(cornerRadius: 46, style: .continuous))
                    .overlay(RoundedRectangle(cornerRadius: 46, style: .continuous).strokeBorder(.white.opacity(0.22), lineWidth: 2))
                    .shadow(color: accent.opacity(0.18), radius: 40)
            }
            .padding(.top, 30)
        }
        .frame(width: 430, height: 932)
        .environment(\.colorScheme, .dark)
    }
}

@MainActor
func renderAppStoreShots(from dir: URL, to out: URL) {
    let hue = { (d: Double) in EquationColor.color(hue: d) }
    let shots: [(String, String, String, Color)] = [
        ("locked", "Draw a curve.\nGet its equation.", "One stroke on the plane. The exact maths behind it.", hue(230)),
        ("searching", "Watch the maths evolve.", "Equations compete to fit your stroke, term by term.", hue(150)),
        ("prune", "It keeps what works.", "Wrong terms are pruned. 0.5 snaps to ½.", hue(70)),
        ("explain", "Tap a term to see\nwhat it does.", "The wave, the tilt: each part drawn on its own.", hue(20)),
        ("share", "Every curve\nbecomes a card.", "Send it, save it, or share the equation.", hue(330)),
        ("shareBack", "Turn it over\nfor the stats.", "1,860 equations tried in 2.4 seconds.", hue(285)),
    ]
    let dest = out.appendingPathComponent("appstore")
    try? FileManager.default.createDirectory(at: dest, withIntermediateDirectories: true)
    for (i, (name, headline, sub, accent)) in shots.enumerated() {
        guard let img = NSImage(contentsOf: dir.appendingPathComponent("\(name).png")) else { print("missing \(name)"); continue }
        let r = ImageRenderer(content: AppStoreShot(headline: headline, sub: sub, image: img, accent: accent))
        r.scale = 3
        guard let cg = r.cgImage,
              let png = NSBitmapImageRep(cgImage: cg).representation(using: .png, properties: [:]) else { continue }
        try? png.write(to: dest.appendingPathComponent(String(format: "%02d-%@.png", i + 1, name)))
        print("wrote appstore \(i + 1) \(name) \(cg.width)x\(cg.height)")
    }
}
