import XCTest
@testable import GraphitiKit

final class ExactSnapTests: XCTestCase {
    func testSnapsNearExactValues() {
        XCTAssertEqual(ExactSnap.snap(0.4999)?.latex, "\\frac{1}{2}")
        XCTAssertEqual(ExactSnap.snap(3.1416)?.latex, "\\pi")
        XCTAssertEqual(ExactSnap.snap(1.4142)?.latex, "\\sqrt{2}")
        XCTAssertEqual(ExactSnap.snap(1.5708)?.latex, "\\frac{\\pi}{2}")
        XCTAssertEqual(ExactSnap.snap(-0.3333)?.latex, "-\\frac{1}{3}")
    }

    func testLeavesOrdinaryNumbersAlone() {
        XCTAssertNil(ExactSnap.snap(0.37))
        XCTAssertNil(ExactSnap.snap(2.0))
        XCTAssertNil(ExactSnap.snap(1.47))
        XCTAssertNil(ExactSnap.snap(1.6))
    }

    func testFormatterUsesSnaps() {
        XCTAssertEqual(EquationFormatter.latex(CurveParams(m: 0.5, a: 1, w: 2, q: 0)), "y = {\\frac{1}{2}}x + \\sin(2x)")
        XCTAssertEqual(EquationFormatter.latex(CurveParams(m: 0, a: 1, w: 3.14159, q: 0)), "y = \\sin({\\pi}x)")
    }
}

final class CardPayloadTests: XCTestCase {
    func testRoundTripsThroughMessageURL() {
        let card = CardPayload(params: CurveParams(m: 0.5, a: 1.25, w: 2, q: 0), author: "lincoln")
        let back = CardPayload(url: card.url)
        XCTAssertEqual(back, card)
        XCTAssertEqual(card.url.scheme, "graphiti")
    }

    func testPlainCaption() {
        XCTAssertEqual(EquationFormatter.plain(CurveParams(m: 0.5, a: 1, w: 2, q: 0.25)), "y = ½x + sin(2x) + ¼x²")
    }

    func testRejectsForeignURLs() {
        XCTAssertNil(CardPayload(url: URL(string: "https://example.com/card?m=1")!))
    }
}

final class CurveEditorTests: XCTestCase {
    func testDragLandsUnderTheFinger() {
        let p = CurveParams(m: 0.5, a: 1, w: 2, q: 0)
        let x = 0.785
        let target = p.y(x) + 0.8
        let q = CurveEditor.drag(p, at: x, to: target)
        XCTAssertEqual(q.y(x), target, accuracy: 1e-4)
    }

    func testDragNeverInventsTerms() {
        let p = CurveParams(m: 0.5, a: 1, w: 2, q: 0)
        let q = CurveEditor.drag(p, at: 1, to: 3)
        XCTAssertEqual(q.q, 0)
    }

    func testDraggingAPeakMostlyChangesHeight() {
        let p = CurveParams(m: 0.5, a: 1, w: 2, q: 0)
        let q = CurveEditor.drag(p, at: 0.785, by: 0.5)
        XCTAssertGreaterThan(abs(q.a - p.a), abs(q.w - p.w))
    }
}

final class GridHapticsTests: XCTestCase {
    func testAxisCrossingIsStrongerThanGridLine() {
        let f = PlaneFrame(size: CGSize(width: 393, height: 852))
        var d = GridCrossingDetector(frame: f)
        XCTAssertNil(d.add(f.point(0.2, 0.2)))
        XCTAssertEqual(d.add(f.point(0.7, 0.2)), .gridLine)
        XCTAssertEqual(d.add(f.point(0.7, -0.2)), .axis)
        XCTAssertNil(d.add(f.point(0.72, -0.21)))
    }
}

final class NamedCurveTests: XCTestCase {
    func testRecognisesFamousCurves() {
        XCTAssertEqual(NamedCurve.identify("r = 1.6(1 - \\cos\\theta)")?.name, "Cardioid")
        XCTAssertEqual(NamedCurve.identify("r = 0.28\\theta")?.name, "Archimedean spiral")
        XCTAssertEqual(NamedCurve.identify("y = 2\\cosh(x/2)")?.name, "Catenary")
        XCTAssertNil(NamedCurve.identify("y = {\\frac{1}{2}}x + \\sin(2x)"))
    }
}

final class EquationColorTests: XCTestCase {
    func testSameEquationSameColour() {
        XCTAssertEqual(EquationColor.hue(for: "y = {\\frac{1}{2}}x + \\sin(2x)"),
                       EquationColor.hue(for: "y=\\frac12x+\\sin(2x)".replacingOccurrences(of: "12", with: "{1}{2}")))
        XCTAssertNotEqual(EquationColor.hue(for: "y = \\sin(2x)"), EquationColor.hue(for: "y = \\sin(3x)"))
    }

    func testEveryHueStaysInGamut() {
        for h in stride(from: 0.0, to: 360, by: 5) {
            let (r, g, b) = EquationColor.rgb(hue: h)
            for c in [r, g, b] { XCTAssert((0...1).contains(c)) }
        }
    }
}
