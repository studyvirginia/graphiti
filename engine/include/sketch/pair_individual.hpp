#ifndef SKETCH_PAIR_INDIVIDUAL_HPP
#define SKETCH_PAIR_INDIVIDUAL_HPP

#include <operon/operon/core/tree.hpp>
#include <operon/operon/core/types.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <vector>

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif

#include "grid_sdf.hpp"
#include "point_segment_kernel.hpp"

namespace Sketch {

struct AmxData {
    bool ok{false};
    std::vector<float> pMat;
    std::vector<float> pSq;
    std::size_t n{0};
};

inline void buildAmxData(std::vector<Operon::Scalar> const& dataX,
                         std::vector<Operon::Scalar> const& dataY,
                         AmxData& out) {
    auto n = dataX.size();
    if (n < 2) { out.ok = false; return; }
    out.n = n;
    out.pMat.resize(2 * n);
    out.pSq.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        float const x = static_cast<float>(dataX[i]);
        float const y = static_cast<float>(dataY[i]);
        out.pMat[2 * i + 0] = x;
        out.pMat[2 * i + 1] = y;
        out.pSq[i] = x * x + y * y;
    }
    out.ok = true;
}

struct PairIndividual {
    Operon::Tree x;
    Operon::Tree y;
    Operon::Scalar shapeErr{std::numeric_limits<Operon::Scalar>::max()};
    Operon::Scalar complexity{0};
    Operon::Scalar constBits{0};
    int rank{0};
    Operon::Scalar crowding{0};
    int age{0};
    std::vector<std::size_t> warmSegmentIndex;
    std::vector<double> warmSegmentT;
    std::size_t warmCompletedSweeps{0};
    std::vector<Operon::Scalar> regionErrCache;
    int sf2LvX{-1};
    int sf2LvY{-1};
    Operon::Scalar sf2ParentErr{0};
    bool sf2Child{false};
    Operon::Scalar errX{0};
    Operon::Scalar errY{0};
    Operon::Scalar chamferErr{0};
    std::vector<Operon::Scalar> perPointResid_;
    std::vector<Operon::Scalar> perPointResidSigned_;
    std::vector<Vec2> perPointNormal_;

    PairIndividual() = default;
    PairIndividual(PairIndividual const&) = default;
    PairIndividual(PairIndividual&&) noexcept = default;
    PairIndividual& operator=(PairIndividual const&) = default;
    PairIndividual& operator=(PairIndividual&&) noexcept = default;

    [[nodiscard]] auto totalLength() const -> std::size_t {
        return x.Length() + y.Length();
    }
};

inline auto structurallyValid(Operon::Tree const& t) -> bool {
    for (auto const& n : t.Nodes()) {
        if ((n.IsConstant() || n.IsVariable()) && n.Arity > 0) { return false; }
    }
    return true;
}

inline auto perTreeCreditMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto constBitsOf(Operon::Tree const& t) -> Operon::Scalar {
    Operon::Scalar bits = 0;
    for (auto const& n : t.Nodes()) {
        if (!n.IsConstant()) { continue; }
        auto const v = static_cast<double>(n.Value);
        if (!std::isfinite(v)) { continue; }
        auto const b = std::round(std::log2(1.0 + std::abs(v)));
        if (std::isfinite(b)) { bits += static_cast<Operon::Scalar>(b); }
    }
    return bits;
}

inline auto constMagnitudeComplexityMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto eliasGammaLength(Operon::Scalar c) -> Operon::Scalar {
    auto const v = static_cast<double>(c);
    if (!std::isfinite(v) || v < 1.0) { return Operon::Scalar{1}; }
    double const L = 2.0 * std::floor(std::log2(v)) + 1.0;
    return static_cast<Operon::Scalar>(L);
}

inline auto priorMinimalComplexityMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto effectiveComplexity(PairIndividual const& p) -> Operon::Scalar {
    Operon::Scalar const base = priorMinimalComplexityMode()
        ? eliasGammaLength(p.complexity) : p.complexity;
    if (!constMagnitudeComplexityMode()) { return base; }
    return base + p.constBits;
}

inline auto errorClampValue() -> Operon::Scalar& {
    static thread_local Operon::Scalar clamp = 0;
    return clamp;
}

inline auto clampedErr(Operon::Scalar e) -> Operon::Scalar {
    auto const c = errorClampValue();
    if (c <= Operon::Scalar{0}) { return e; }
    return (e < c) ? c : e;
}

inline auto floorOnChamferMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto cisFitnessMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto arcLengthGuardRatio() -> Operon::Scalar& {
    static thread_local Operon::Scalar ratio = 0;
    return ratio;
}

inline auto fastNonDominatedSortMode() -> bool& {
    static thread_local bool mode = true;
    return mode;
}

inline auto useExactNlogNSortMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto clampedShapeErr(PairIndividual const& p) -> Operon::Scalar {
    if (!floorOnChamferMode()) { return clampedErr(p.shapeErr); }
    auto const c = errorClampValue();
    if (c <= Operon::Scalar{0}) { return p.shapeErr; }
    return (p.chamferErr <= c) ? c : p.shapeErr;
}

inline auto dominates(PairIndividual const& a, PairIndividual const& b) -> bool {
    auto const ac = effectiveComplexity(a), bc = effectiveComplexity(b);
    if (perTreeCreditMode()) {
        auto axe = clampedErr(a.errX), aye = clampedErr(a.errY);
        auto bxe = clampedErr(b.errX), bye = clampedErr(b.errY);
        bool le = axe <= bxe && aye <= bye && ac <= bc;
        bool lt = axe <  bxe || aye <  bye || ac <  bc;
        return le && lt;
    }
    auto ae = clampedShapeErr(a), be = clampedShapeErr(b);
    bool le = ae <= be && ac <= bc;
    bool lt = ae <  be || ac <  bc;
    return le && lt;
}

inline auto crowdedLess(PairIndividual const& a, PairIndividual const& b) -> bool {
    if (a.rank != b.rank) { return a.rank < b.rank; }
    return a.crowding > b.crowding;
}

inline void normalizeInPlace(std::vector<Operon::Scalar>& x,
                             std::vector<Operon::Scalar>& y) {
    auto [xlo, xhi] = std::minmax_element(x.begin(), x.end());
    auto [ylo, yhi] = std::minmax_element(y.begin(), y.end());
    auto cx = (*xhi + *xlo) / 2, cy = (*yhi + *ylo) / 2;
    auto sp = std::max({*xhi - *xlo, *yhi - *ylo, Operon::Scalar{1e-9}});
    for (auto& v : x) { v = (v - cx) / sp; }
    for (auto& v : y) { v = (v - cy) / sp; }
}

inline void applyRotationGauge(std::vector<Operon::Scalar>& x,
                                std::vector<Operon::Scalar>& y,
                                double& appliedAngle) {
    appliedAngle = 0.0;
    auto const n = x.size();
    if (n < 2 || y.size() != n) { return; }

    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= static_cast<double>(n);
    my /= static_cast<double>(n);

    double cxx = 0.0, cyy = 0.0, cxy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const dx = static_cast<double>(x[i]) - mx;
        double const dy = static_cast<double>(y[i]) - my;
        cxx += dx * dx; cyy += dy * dy; cxy += dx * dy;
    }
    cxx /= static_cast<double>(n); cyy /= static_cast<double>(n); cxy /= static_cast<double>(n);

    double const tr = cxx + cyy;
    double const diff = cxx - cyy;
    double const disc = std::sqrt(diff * diff / 4.0 + cxy * cxy);
    double const lambda1 = tr / 2.0 + disc;
    double const lambda2 = tr / 2.0 - disc;

    double const denom = lambda1 + lambda2;
    if (!(denom > 1e-12) || (lambda1 - lambda2) < 1e-6 * denom) {
        return;
    }

    double const theta = 0.5 * std::atan2(2.0 * cxy, diff);
    double const ct = std::cos(-theta), st = std::sin(-theta);

    for (std::size_t i = 0; i < n; ++i) {
        double const dx = static_cast<double>(x[i]) - mx;
        double const dy = static_cast<double>(y[i]) - my;
        x[i] = static_cast<Operon::Scalar>(mx + dx * ct - dy * st);
        y[i] = static_cast<Operon::Scalar>(my + dx * st + dy * ct);
    }
    appliedAngle = -theta;
}

struct NoiseModel {
    std::vector<Operon::Scalar> s2;
    std::vector<Operon::Scalar> nu;
    std::vector<Operon::Scalar> logVar;
    [[nodiscard]] auto empty() const -> bool { return s2.empty(); }
    [[nodiscard]] auto size() const -> std::size_t { return s2.size(); }
};

inline auto noiseIndexAtArc(std::size_t i, std::size_t N, std::size_t K) -> std::size_t {
    if (K == 0 || N <= 1) { return 0; }
    return std::min(i * K / (N - 1), K - 1);
}

inline auto noiseRegionLUT(std::size_t N, std::size_t K) -> std::vector<std::size_t> const& {
    thread_local std::vector<std::size_t> lut;
    thread_local std::size_t cachedN = static_cast<std::size_t>(-1);
    thread_local std::size_t cachedK = static_cast<std::size_t>(-1);
    if (N != cachedN || K != cachedK) {
        lut.resize(N);
        for (std::size_t i = 0; i < N; ++i) { lut[i] = noiseIndexAtArc(i, N, K); }
        cachedN = N;
        cachedK = K;
    }
    return lut;
}

inline auto pointClampMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

inline auto pointLoss(Operon::Scalar r2, Operon::Scalar s2, Operon::Scalar nu)
    -> Operon::Scalar
{
    if (!(s2 > Operon::Scalar{0}) || !std::isfinite(static_cast<double>(s2))) {
        s2 = Operon::Scalar{1};
    }
    if (nu <= Operon::Scalar{0}) {
        auto val = r2 / s2;
        if (pointClampMode()) {
            return std::max(val, Operon::Scalar{1.0});
        }
        return val;
    }
    auto val = Operon::Scalar{0.5} * (nu + Operon::Scalar{2})
         * std::log1p(r2 / (nu * s2));
    if (pointClampMode()) {
        return std::max(val, Operon::Scalar{1.0});
    }
    return val;
}

inline auto pointSegDist2(Operon::Scalar px, Operon::Scalar py,
                          Operon::Scalar ax, Operon::Scalar ay,
                          Operon::Scalar bx, Operon::Scalar by) -> Operon::Scalar {
    auto dx = bx - ax, dy = by - ay;
    auto len2 = dx*dx + dy*dy;
    Operon::Scalar t = 0;
    if (len2 > Operon::Scalar{1e-12}) {
        t = ((px-ax)*dx + (py-ay)*dy) / len2;
        t = std::clamp(t, Operon::Scalar{0}, Operon::Scalar{1});
    }
    auto qx = ax + t*dx, qy = ay + t*dy;
    auto ex = px-qx, ey = py-qy;
    return ex*ex + ey*ey;
}

inline auto pointToPolyline2(Operon::Scalar px, Operon::Scalar py,
                             std::vector<Operon::Scalar> const& X,
                             std::vector<Operon::Scalar> const& Y) -> Operon::Scalar {
    Operon::Scalar best = std::numeric_limits<Operon::Scalar>::max();
    for (std::size_t j = 0; j + 1 < X.size(); ++j) {
        auto d2 = pointSegDist2(px, py, X[j], Y[j], X[j+1], Y[j+1]);
        if (d2 < best) { best = d2; }
    }
    return best;
}

#ifdef SKETCH_COUNT_OPS
struct NnCounters {
    unsigned long long items = 0;
    unsigned long long boxes = 0;
    unsigned long long queries = 0;
};
inline auto nnCounters() -> NnCounters& {
    static thread_local NnCounters c;
    return c;
}
#define SKETCH_COUNT_ITEM(n) (::Sketch::nnCounters().items += (n))
#define SKETCH_COUNT_BOX(n)  (::Sketch::nnCounters().boxes += (n))
#define SKETCH_COUNT_QUERY() (++::Sketch::nnCounters().queries)
#else
#define SKETCH_COUNT_ITEM(n) ((void)0)
#define SKETCH_COUNT_BOX(n)  ((void)0)
#define SKETCH_COUNT_QUERY() ((void)0)
#endif

struct CurveBlocks {
    bool ok = false;
    int nblocks = 0;
    int blockSize = 0;
    int nitems = 0;
    std::vector<Operon::Scalar> box;
    int nsuper = 0;
    int superSize = 0;
    std::vector<Operon::Scalar> sbox;
};

inline constexpr int kBlockMinN = 24;
inline constexpr int kBlockSize = 16;
#ifndef SKETCH_SUPER_SIZE
#define SKETCH_SUPER_SIZE 4
#endif
inline constexpr int kSuperSize = SKETCH_SUPER_SIZE;
inline constexpr int kSuperMinBlocks = 8;

inline void buildBlocks(std::vector<Operon::Scalar> const& X, std::vector<Operon::Scalar> const& Y,
                        bool segments, CurveBlocks& b)
{
    b.ok = false;
    auto const n = X.size();
    if (Y.size() != n) { return; }
    auto const nitems = static_cast<int>(segments ? (n >= 1 ? n - 1 : 0) : n);
    if (nitems < kBlockMinN) { return; }
    b.nitems = nitems;
    b.blockSize = kBlockSize;
    b.nblocks = (nitems + b.blockSize - 1) / b.blockSize;
    b.box.resize(4 * static_cast<std::size_t>(b.nblocks));
    for (int k = 0; k < b.nblocks; ++k) {
        int const lo = k * b.blockSize;
        int const hi = std::min(lo + b.blockSize, nitems);
        int const vlo = lo;
        int const vhi = segments ? hi : hi - 1;
        auto xlo = X[static_cast<std::size_t>(vlo)], xhi = xlo;
        auto ylo = Y[static_cast<std::size_t>(vlo)], yhi = ylo;
        for (int v = vlo; v <= vhi; ++v) {
            auto const x = X[static_cast<std::size_t>(v)], y = Y[static_cast<std::size_t>(v)];
            if (!std::isfinite(x) || !std::isfinite(y)) { return; }
            xlo = std::min(xlo, x); xhi = std::max(xhi, x);
            ylo = std::min(ylo, y); yhi = std::max(yhi, y);
        }
        auto* o = &b.box[4 * static_cast<std::size_t>(k)];
        o[0] = xlo; o[1] = ylo; o[2] = xhi; o[3] = yhi;
    }
    b.nsuper = 0; b.superSize = 0; b.sbox.clear();
    if (b.nblocks >= kSuperMinBlocks) {
        b.superSize = kSuperSize;
        b.nsuper = (b.nblocks + b.superSize - 1) / b.superSize;
        b.sbox.resize(4 * static_cast<std::size_t>(b.nsuper));
        for (int s = 0; s < b.nsuper; ++s) {
            int const klo = s * b.superSize;
            int const khi = std::min(klo + b.superSize, b.nblocks);
            auto const* f = &b.box[4 * static_cast<std::size_t>(klo)];
            auto xlo = f[0], ylo = f[1], xhi = f[2], yhi = f[3];
            for (int k = klo + 1; k < khi; ++k) {
                auto const* o = &b.box[4 * static_cast<std::size_t>(k)];
                xlo = std::min(xlo, o[0]); ylo = std::min(ylo, o[1]);
                xhi = std::max(xhi, o[2]); yhi = std::max(yhi, o[3]);
            }
            auto* t = &b.sbox[4 * static_cast<std::size_t>(s)];
            t[0] = xlo; t[1] = ylo; t[2] = xhi; t[3] = yhi;
        }
    }
    b.ok = true;
}

inline auto pointBoxDist2(double px, double py, Operon::Scalar const* o) -> double
{
    double const oxlo = static_cast<double>(o[0]);
    double const oylo = static_cast<double>(o[1]);
    double const oxhi = static_cast<double>(o[2]);
    double const oyhi = static_cast<double>(o[3]);
    double const dx = std::max(std::max(oxlo - px, px - oxhi), 0.0);
    double const dy = std::max(std::max(oylo - py, py - oyhi), 0.0);
    return dx * dx + dy * dy;
}

inline auto pointBoxDist2(Operon::Scalar px, Operon::Scalar py, Operon::Scalar const* o) -> double
{
    return pointBoxDist2(static_cast<double>(px), static_cast<double>(py), o);
}

inline auto blockExcluded(double lb2, Operon::Scalar best) -> bool {
    return lb2 * (1.0 - 1e-6) > static_cast<double>(best);
}

inline auto pointToPolylineBlocks2(Operon::Scalar px, Operon::Scalar py,
                                   std::vector<Operon::Scalar> const& X,
                                   std::vector<Operon::Scalar> const& Y,
                                   CurveBlocks const& b,
                                   int* hint = nullptr) -> Operon::Scalar
{
    if (!b.ok || !std::isfinite(px) || !std::isfinite(py)) { return pointToPolyline2(px, py, X, Y); }
    SKETCH_COUNT_QUERY();
    Operon::Scalar best = std::numeric_limits<Operon::Scalar>::max();
    int bestBlock = -1;
    auto scan = [&](int k) {
        int const lo = k * b.blockSize;
        int const hi = std::min(lo + b.blockSize, b.nitems);
        SKETCH_COUNT_ITEM(hi - lo);
        for (int j = lo; j < hi; ++j) {
            auto const s = static_cast<std::size_t>(j);
            auto const d2 = pointSegDist2(px, py, X[s], Y[s], X[s + 1], Y[s + 1]);
            if (d2 < best) { best = d2; bestBlock = k; }
        }
    };
    int const hintK = (hint != nullptr) ? *hint : -1;
    if (hintK >= 0 && hintK < b.nblocks) { scan(hintK); }
    double const dpx = static_cast<double>(px);
    double const dpy = static_cast<double>(py);
    if (b.nsuper > 0) {
        for (int sb = 0; sb < b.nsuper; ++sb) {
            if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dpx, dpy, &b.sbox[4 * static_cast<std::size_t>(sb)]), best)) { continue; }
            int const klo = sb * b.superSize;
            int const khi = std::min(klo + b.superSize, b.nblocks);
            for (int k = klo; k < khi; ++k) {
                if (k == hintK) { continue; }
                if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dpx, dpy, &b.box[4 * static_cast<std::size_t>(k)]), best)) { continue; }
                scan(k);
            }
        }
    } else {
        for (int k = 0; k < b.nblocks; ++k) {
            if (k == hintK) { continue; }
            if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dpx, dpy, &b.box[4 * static_cast<std::size_t>(k)]), best)) { continue; }
            scan(k);
        }
    }
    if (hint != nullptr && bestBlock >= 0) { *hint = bestBlock; }
    return best;
}


inline auto blockExcludedD(double lb2, double best) -> bool {
    return lb2 * (1.0 - 1e-6) > best;
}
inline auto pointBoxDist2D(double px, double py, Operon::Scalar const* o) -> double {
    return pointBoxDist2(px, py, o);
}

struct SegArgmin {
    double d2 = std::numeric_limits<double>::max();
    std::size_t j = 0;
    double tau = 0.0;
};

inline auto nearestSegmentD(double px, double py,
                            std::vector<Operon::Scalar> const& X,
                            std::vector<Operon::Scalar> const& Y,
                            CurveBlocks const& b, int* hint = nullptr) -> SegArgmin
{
    auto const M = X.size();
    SegArgmin r;
    std::size_t found = M;
    auto one = [&](std::size_t j, double& d2, double& t) {
        double const ax = static_cast<double>(X[j]), ay = static_cast<double>(Y[j]);
        double const bx = static_cast<double>(X[j + 1]), by = static_cast<double>(Y[j + 1]);
        double const ex = bx - ax, ey = by - ay;
        double const len2 = ex * ex + ey * ey;
        t = 0.0;
        if (len2 > 1e-12) { t = std::clamp(((px - ax) * ex + (py - ay) * ey) / len2, 0.0, 1.0); }
        double const zx = ax + t * ex, zy = ay + t * ey;
        d2 = (px - zx) * (px - zx) + (py - zy) * (py - zy);
    };
    auto take = [&](std::size_t j) {
        double d2 = 0, t = 0;
        one(j, d2, t);
        if (d2 < r.d2) { r.d2 = d2; r.j = j; r.tau = t; found = j; }
        else if (found < M && d2 == r.d2 && j < found) { r.j = j; r.tau = t; found = j; }
    };
    if (!b.ok || !std::isfinite(px) || !std::isfinite(py)) {
        for (std::size_t j = 0; j + 1 < M; ++j) { take(j); }
        if (found == M) { r.j = 0; r.tau = 0.0; }
        return r;
    }
    int bestBlock = -1;
    auto scan = [&](int k) {
        int const lo = k * b.blockSize;
        int const hi = std::min(lo + b.blockSize, b.nitems);
        SKETCH_COUNT_ITEM(hi - lo);
        for (int j = lo; j < hi; ++j) {
            auto const before = found;
            take(static_cast<std::size_t>(j));
            if (found != before) { bestBlock = k; }
        }
    };
    SKETCH_COUNT_QUERY();
    int const hintK = (hint != nullptr) ? *hint : -1;
    if (hintK >= 0 && hintK < b.nblocks) { scan(hintK); }
    if (b.nsuper > 0) {
        for (int sb = 0; sb < b.nsuper; ++sb) {
            if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.sbox[4 * static_cast<std::size_t>(sb)]), r.d2)) { continue; }
            int const klo = sb * b.superSize;
            int const khi = std::min(klo + b.superSize, b.nblocks);
            for (int k = klo; k < khi; ++k) {
                if (k == hintK) { continue; }
                if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.box[4 * static_cast<std::size_t>(k)]), r.d2)) { continue; }
                scan(k);
            }
        }
    } else {
        for (int k = 0; k < b.nblocks; ++k) {
            if (k == hintK) { continue; }
            if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.box[4 * static_cast<std::size_t>(k)]), r.d2)) { continue; }
            scan(k);
        }
    }
    if (found == M) { r.j = 0; r.tau = 0.0; }
    if (hint != nullptr && bestBlock >= 0) { *hint = bestBlock; }
    return r;
}

struct PtArgmin {
    double d2 = std::numeric_limits<double>::max();
    std::size_t i = 0;
};

inline auto nearestPointD(double px, double py,
                          std::vector<Operon::Scalar> const& X,
                          std::vector<Operon::Scalar> const& Y,
                          CurveBlocks const& b, int* hint = nullptr) -> PtArgmin
{
    auto const N = X.size();
    PtArgmin r;
    std::size_t found = N;
    auto take = [&](std::size_t i) {
        double const ddx = px - static_cast<double>(X[i]);
        double const ddy = py - static_cast<double>(Y[i]);
        double const d2 = ddx * ddx + ddy * ddy;
        if (d2 < r.d2) { r.d2 = d2; r.i = i; found = i; }
        else if (found < N && d2 == r.d2 && i < found) { r.i = i; found = i; }
    };
    if (!b.ok || !std::isfinite(px) || !std::isfinite(py)) {
        for (std::size_t i = 0; i < N; ++i) { take(i); }
        if (found == N) { r.i = 0; }
        return r;
    }
    int bestBlock = -1;
    auto scan = [&](int k) {
        int const lo = k * b.blockSize;
        int const hi = std::min(lo + b.blockSize, b.nitems);
        SKETCH_COUNT_ITEM(hi - lo);
        for (int i = lo; i < hi; ++i) {
            auto const before = found;
            take(static_cast<std::size_t>(i));
            if (found != before) { bestBlock = k; }
        }
    };
    SKETCH_COUNT_QUERY();
    int const hintK = (hint != nullptr) ? *hint : -1;
    if (hintK >= 0 && hintK < b.nblocks) { scan(hintK); }
    if (b.nsuper > 0) {
        for (int sb = 0; sb < b.nsuper; ++sb) {
            if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.sbox[4 * static_cast<std::size_t>(sb)]), r.d2)) { continue; }
            int const klo = sb * b.superSize;
            int const khi = std::min(klo + b.superSize, b.nblocks);
            for (int k = klo; k < khi; ++k) {
                if (k == hintK) { continue; }
                if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.box[4 * static_cast<std::size_t>(k)]), r.d2)) { continue; }
                scan(k);
            }
        }
    } else {
        for (int k = 0; k < b.nblocks; ++k) {
            if (k == hintK) { continue; }
            if (SKETCH_COUNT_BOX(1), blockExcludedD(pointBoxDist2D(px, py, &b.box[4 * static_cast<std::size_t>(k)]), r.d2)) { continue; }
            scan(k);
        }
    }
    if (found == N) { r.i = 0; }
    if (hint != nullptr && bestBlock >= 0) { *hint = bestBlock; }
    return r;
}

inline auto shapeDistanceBrute(std::vector<Operon::Scalar> const& gridX,
                          std::vector<Operon::Scalar> const& gridY,
                          std::vector<Operon::Scalar> const& dataX,
                          std::vector<Operon::Scalar> const& dataY,
                          NoiseModel const& noise = {},
                          std::vector<Operon::Scalar>* regionOut = nullptr,
                          int  = 4) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
    auto const M = gridX.size();
    auto const N = dataX.size();
    if (M < 2 || N < 2) { return BIG; }
    for (std::size_t j = 0; j < M; ++j) {
        if (!std::isfinite(gridX[j]) || !std::isfinite(gridY[j])) { return BIG; }
    }

    std::vector<Operon::Scalar> cx(gridX), cy(gridY);
    normalizeInPlace(cx, cy);

    if (arcLengthGuardRatio() > Operon::Scalar{0}) {
        Operon::Scalar candLen = 0;
        for (std::size_t j = 0; j + 1 < M; ++j) {
            candLen += std::hypot(cx[j + 1] - cx[j], cy[j + 1] - cy[j]);
        }
        Operon::Scalar dataLen = 0;
        for (std::size_t i = 0; i + 1 < N; ++i) {
            dataLen += std::hypot(dataX[i + 1] - dataX[i], dataY[i + 1] - dataY[i]);
        }
        if (dataLen > Operon::Scalar{1e-9} && candLen > arcLengthGuardRatio() * dataLen) { return BIG; }
    }

    auto const K = noise.size();
    auto const& rlut = noiseRegionLUT(N, K);
    auto regionOf = [&rlut](std::size_t i) -> std::size_t { return rlut[i]; };

    Operon::Scalar accData = 0;
    for (std::size_t i = 0; i < N; ++i) {
        auto r2 = pointToPolyline2(dataX[i], dataY[i], cx, cy);
        if (K == 0) { accData += r2; }
        else { auto k = regionOf(i); accData += pointLoss(r2, noise.s2[k], noise.nu[k]); }
    }
    accData /= static_cast<Operon::Scalar>(N);

    Operon::Scalar accCurve = 0;
    for (std::size_t j = 0; j < M; ++j) {
        Operon::Scalar bestd = BIG; std::size_t bestI = 0;
        for (std::size_t i = 0; i < N; ++i) {
            auto ddx = cx[j] - dataX[i], ddy = cy[j] - dataY[i];
            auto d2 = ddx*ddx + ddy*ddy;
            if (d2 < bestd) { bestd = d2; bestI = i; }
        }
        if (K == 0) { accCurve += bestd; }
        else { auto k = regionOf(bestI); accCurve += pointLoss(bestd, noise.s2[k], noise.nu[k]); }
    }
    accCurve /= static_cast<Operon::Scalar>(M);

    Operon::Scalar best = Operon::Scalar{0.5} * (accData + accCurve);

    if (regionOut != nullptr) {
        auto const R = (K == 0) ? std::size_t{1} : K;
        regionOut->assign(R, Operon::Scalar{0});
        std::vector<std::size_t> cnt(R, 0);
        for (std::size_t i = 0; i < N; ++i) {
            auto r2 = pointToPolyline2(dataX[i], dataY[i], cx, cy);
            auto k = std::min(i * R / N, R - 1);
            (*regionOut)[k] += (K == 0) ? r2 : pointLoss(r2, noise.s2[k], noise.nu[k]);
            ++cnt[k];
        }
        for (std::size_t k = 0; k < R; ++k) {
            if (cnt[k] > 0) { (*regionOut)[k] /= static_cast<Operon::Scalar>(cnt[k]); }
        }
    }
    return best;
}

inline auto shapeDistance(std::vector<Operon::Scalar> const& gridX,
                          std::vector<Operon::Scalar> const& gridY,
                          std::vector<Operon::Scalar> const& dataX,
                          std::vector<Operon::Scalar> const& dataY,
                          NoiseModel const& noise = {},
                          std::vector<Operon::Scalar>* regionOut = nullptr,
                          int  = 4,
                          CurveBlocks const* dataBlocks = nullptr,
                          GridSdf const* gridSdf = nullptr,
                          AmxData const* amxData = nullptr) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
    auto const M = gridX.size();
    auto const N = dataX.size();
    if (M < 2 || N < 2) { return BIG; }

    thread_local std::vector<Operon::Scalar> cx, cy;
    cx.resize(M);
    cy.resize(M);
    {
        Operon::Scalar gx0 = gridX[0], gy0 = gridY[0];
        if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return BIG; }
        Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
        cx[0] = gx0;
        cy[0] = gy0;
        for (std::size_t j = 1; j < M; ++j) {
            Operon::Scalar const gx = gridX[j], gy = gridY[j];
            if (!std::isfinite(gx) || !std::isfinite(gy)) { return BIG; }
            cx[j] = gx;
            cy[j] = gy;
            if (gx < xlo) { xlo = gx; } else if (gx > xhi) { xhi = gx; }
            if (gy < ylo) { ylo = gy; } else if (gy > yhi) { yhi = gy; }
        }
        Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
        Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
        for (std::size_t j = 0; j < M; ++j) {
            cx[j] = (cx[j] - cxc) / sp;
            cy[j] = (cy[j] - cyc) / sp;
        }
    }

    if (arcLengthGuardRatio() > Operon::Scalar{0}) {
        Operon::Scalar candLen = 0;
        for (std::size_t j = 0; j + 1 < M; ++j) {
            candLen += std::hypot(cx[j + 1] - cx[j], cy[j + 1] - cy[j]);
        }
        Operon::Scalar dataLen = 0;
        for (std::size_t i = 0; i + 1 < N; ++i) {
            dataLen += std::hypot(dataX[i + 1] - dataX[i], dataY[i + 1] - dataY[i]);
        }
        if (dataLen > Operon::Scalar{1e-9} && candLen > arcLengthGuardRatio() * dataLen) { return BIG; }
    }

    auto const K = noise.size();
    auto const& rlut = noiseRegionLUT(N, K);
    auto regionOf = [&rlut](std::size_t i) -> std::size_t { return rlut[i]; };

    thread_local CurveBlocks segB, ptB;
    buildBlocks(cx, cy, true, segB);
    if (dataBlocks == nullptr) { buildBlocks(dataX, dataY, false, ptB); }
    CurveBlocks const& ptRef = (dataBlocks != nullptr) ? *dataBlocks : ptB;

    Operon::Scalar accData = 0;
    int segHint = -1;
    for (std::size_t i = 0; i < N; ++i) {
        auto r2 = pointToPolylineBlocks2(dataX[i], dataY[i], cx, cy, segB, &segHint);
        if (K == 0) { accData += r2; }
        else { auto k = regionOf(i); accData += pointLoss(r2, noise.s2[k], noise.nu[k]); }
    }
    accData /= static_cast<Operon::Scalar>(N);

    Operon::Scalar accCurve = 0;
    if (gridSdf != nullptr && gridSdf->ok) {
        for (std::size_t j = 0; j < M; ++j) {
            uint8_t reg = 0;
            Operon::Scalar const d2 = queryGridSdf(*gridSdf, cx[j], cy[j], reg);
            if (K == 0) { accCurve += d2; }
            else { accCurve += pointLoss(d2, noise.s2[reg], noise.nu[reg]); }
        }
#if defined(__APPLE__)
    } else if (amxData != nullptr && amxData->ok) {
        thread_local std::vector<float> qMat, qSq, dMat;
        qMat.resize(2 * M);
        qSq.resize(M);
        dMat.resize(N * M);

        for (std::size_t j = 0; j < M; ++j) {
            float const x = static_cast<float>(cx[j]);
            float const y = static_cast<float>(cy[j]);
            qMat[2 * j + 0] = x;
            qMat[2 * j + 1] = y;
            qSq[j] = x * x + y * y;
        }

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    static_cast<int>(N), static_cast<int>(M), 2,
                    -2.0f, amxData->pMat.data(), 2,
                    qMat.data(), 2,
                    0.0f, dMat.data(), static_cast<int>(M));

        constexpr float kNearFloorThresh = 4e-4f;

        for (std::size_t j = 0; j < M; ++j) {
            float bestd = BIG;
            std::size_t bestI = 0;
            float const q2 = qSq[j];
            for (std::size_t i = 0; i < N; ++i) {
                float const d2 = dMat[i * M + j] + amxData->pSq[i] + q2;
                if (d2 < bestd) {
                    bestd = d2;
                    bestI = i;
                }
            }
            if (bestd <= kNearFloorThresh) {
                float const dx = static_cast<float>(cx[j]) - static_cast<float>(dataX[bestI]);
                float const dy = static_cast<float>(cy[j]) - static_cast<float>(dataY[bestI]);
                bestd = dx * dx + dy * dy;
            }
            if (K == 0) { accCurve += static_cast<Operon::Scalar>(bestd); }
            else {
                auto k = regionOf(bestI);
                accCurve += pointLoss(static_cast<Operon::Scalar>(bestd), noise.s2[k], noise.nu[k]);
            }
        }
#endif
    } else {
        int ptHint = -1;
        for (std::size_t j = 0; j < M; ++j) {
            Operon::Scalar bestd = BIG; std::size_t bestI = 0;
            if (ptRef.ok) {
                std::size_t found = N;
                int foundBlock = -1;
                SKETCH_COUNT_QUERY();
                auto scan = [&](int k) {
                    int const lo = k * ptRef.blockSize;
                    int const hi = std::min(lo + ptRef.blockSize, ptRef.nitems);
                    SKETCH_COUNT_ITEM(hi - lo);
                    for (int ii = lo; ii < hi; ++ii) {
                        auto const i = static_cast<std::size_t>(ii);
                        auto ddx = cx[j] - dataX[i], ddy = cy[j] - dataY[i];
                        auto d2 = ddx*ddx + ddy*ddy;
                        if (d2 < bestd) { bestd = d2; found = i; foundBlock = k; }
                        else if (found < N && d2 == bestd && i < found) { found = i; foundBlock = k; }
                    }
                };
                if (ptHint >= 0 && ptHint < ptRef.nblocks) { scan(ptHint); }
                double const dcxj = static_cast<double>(cx[j]);
                double const dcyj = static_cast<double>(cy[j]);
                if (ptRef.nsuper > 0) {
                    for (int sb = 0; sb < ptRef.nsuper; ++sb) {
                        if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dcxj, dcyj, &ptRef.sbox[4 * static_cast<std::size_t>(sb)]), bestd)) { continue; }
                        int const klo = sb * ptRef.superSize;
                        int const khi = std::min(klo + ptRef.superSize, ptRef.nblocks);
                        for (int k = klo; k < khi; ++k) {
                            if (k == ptHint) { continue; }
                            if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dcxj, dcyj, &ptRef.box[4 * static_cast<std::size_t>(k)]), bestd)) { continue; }
                            scan(k);
                        }
                    }
                } else {
                    for (int k = 0; k < ptRef.nblocks; ++k) {
                        if (k == ptHint) { continue; }
                        if (SKETCH_COUNT_BOX(1), blockExcluded(pointBoxDist2(dcxj, dcyj, &ptRef.box[4 * static_cast<std::size_t>(k)]), bestd)) { continue; }
                        scan(k);
                    }
                }
                bestI = (found < N) ? found : std::size_t{0};
                if (foundBlock >= 0) { ptHint = foundBlock; }
            } else {
                for (std::size_t i = 0; i < N; ++i) {
                    auto ddx = cx[j] - dataX[i], ddy = cy[j] - dataY[i];
                    auto d2 = ddx*ddx + ddy*ddy;
                    if (d2 < bestd) { bestd = d2; bestI = i; }
                }
            }
            if (K == 0) { accCurve += bestd; }
            else { auto k = regionOf(bestI); accCurve += pointLoss(bestd, noise.s2[k], noise.nu[k]); }
        }
    }
    accCurve /= static_cast<Operon::Scalar>(M);

    Operon::Scalar best = Operon::Scalar{0.5} * (accData + accCurve);

    if (regionOut != nullptr) {
        auto const R = (K == 0) ? std::size_t{1} : K;
        regionOut->assign(R, Operon::Scalar{0});
        std::vector<std::size_t> cnt(R, 0);
        for (std::size_t i = 0; i < N; ++i) {
            auto r2 = pointToPolylineBlocks2(dataX[i], dataY[i], cx, cy, segB);
            auto k = std::min(i * R / N, R - 1);
            (*regionOut)[k] += (K == 0) ? r2 : pointLoss(r2, noise.s2[k], noise.nu[k]);
            ++cnt[k];
        }
        for (std::size_t k = 0; k < R; ++k) {
            if (cnt[k] > 0) { (*regionOut)[k] /= static_cast<Operon::Scalar>(cnt[k]); }
        }
    }
    return best;
}

inline void arcResample(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                        int n, std::vector<Operon::Scalar>& ox, std::vector<Operon::Scalar>& oy)
{
    auto const M = x.size();
    ox.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    oy.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    if (M < 2) { if (M == 1) { std::fill(ox.begin(), ox.end(), x[0]); std::fill(oy.begin(), oy.end(), y[0]); } return; }
    thread_local std::vector<Operon::Scalar> d;
    d.assign(M, Operon::Scalar{0});
    for (std::size_t i = 1; i < M; ++i) {
        auto dx = x[i] - x[i-1], dy = y[i] - y[i-1];
        d[i] = d[i-1] + std::sqrt(dx*dx + dy*dy);
    }
    auto L = d[M-1];
    if (L < Operon::Scalar{1e-12}) { std::fill(ox.begin(), ox.end(), x[0]); std::fill(oy.begin(), oy.end(), y[0]); return; }
    std::size_t seg = 0;
    for (int k = 0; k < n; ++k) {
        auto s = L * static_cast<Operon::Scalar>(k) / static_cast<Operon::Scalar>(n - 1);
        while (seg + 1 < M - 1 && d[seg + 1] < s) { ++seg; }
        auto d0 = d[seg], d1 = d[seg + 1];
        auto tt = (d1 > d0) ? (s - d0) / (d1 - d0) : Operon::Scalar{0};
        ox[static_cast<std::size_t>(k)] = x[seg] + tt * (x[seg + 1] - x[seg]);
        oy[static_cast<std::size_t>(k)] = y[seg] + tt * (y[seg + 1] - y[seg]);
    }
}

inline auto polylineTurning(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                            bool closed) -> double
{
    constexpr double TWO_PI = 6.283185307179586;
    constexpr double PI = 3.141592653589793;
    std::size_t const M = std::min(x.size(), y.size());
    if (M < 3) { return 0.0; }
    double total = 0.0, prev = 0.0;
    bool have = false;
    std::size_t const nseg = closed ? M : M - 1;
    for (std::size_t i = 0; i < nseg; ++i) {
        std::size_t const j = (i + 1) % M;
        double const dx = static_cast<double>(x[j]) - static_cast<double>(x[i]);
        double const dy = static_cast<double>(y[j]) - static_cast<double>(y[i]);
        if (dx == 0.0 && dy == 0.0) { continue; }
        double const ang = std::atan2(dy, dx);
        if (have) {
            double d = ang - prev;
            while (d > PI) { d -= TWO_PI; }
            while (d < -PI) { d += TWO_PI; }
            total += std::abs(d);
        }
        prev = ang; have = true;
    }
    if (closed && have) {
        std::size_t i = 0;
        while (i < M) {
            std::size_t const j = (i + 1) % M;
            double const dx = static_cast<double>(x[j]) - static_cast<double>(x[i]);
            double const dy = static_cast<double>(y[j]) - static_cast<double>(y[i]);
            if (dx != 0.0 || dy != 0.0) {
                double d = std::atan2(dy, dx) - prev;
                while (d > PI) { d -= TWO_PI; }
                while (d < -PI) { d += TWO_PI; }
                total += std::abs(d);
                break;
            }
            ++i;
        }
    }
    return total / TWO_PI;
}

inline auto fourierLowPass(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                           bool closed, int K, int n,
                           std::vector<Operon::Scalar>& sx, std::vector<Operon::Scalar>& sy) -> bool
{
    constexpr double TWO_PI = 6.283185307179586;
    sx.clear(); sy.clear();
    if (x.size() < 3 || y.size() < 3 || n < 8 || K < 1) { return false; }
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) { return false; }
    }
    std::vector<Operon::Scalar> rx, ry;
    arcResample(x, y, n, rx, ry);
    std::size_t const P = closed ? static_cast<std::size_t>(n) : static_cast<std::size_t>(2 * n);
    std::vector<double> zr(P), zi(P);
    for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j) {
        zr[j] = static_cast<double>(rx[j]); zi[j] = static_cast<double>(ry[j]);
    }
    if (!closed) {
        for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j) {
            zr[static_cast<std::size_t>(n) + j] = static_cast<double>(rx[static_cast<std::size_t>(n) - 1 - j]);
            zi[static_cast<std::size_t>(n) + j] = static_cast<double>(ry[static_cast<std::size_t>(n) - 1 - j]);
        }
    }
    int const KK = closed ? K : 2 * K;
    std::vector<double> cr(static_cast<std::size_t>(2 * KK + 1)), ci(static_cast<std::size_t>(2 * KK + 1));
    double const Pd = static_cast<double>(P);
    for (int k = -KK; k <= KK; ++k) {
        double sr = 0.0, si = 0.0;
        for (std::size_t j = 0; j < P; ++j) {
            double const ang = -TWO_PI * static_cast<double>(k) * static_cast<double>(j) / Pd;
            double const c = std::cos(ang), s = std::sin(ang);
            sr += zr[j] * c - zi[j] * s;
            si += zr[j] * s + zi[j] * c;
        }
        cr[static_cast<std::size_t>(k + KK)] = sr / Pd;
        ci[static_cast<std::size_t>(k + KK)] = si / Pd;
    }
    sx.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    sy.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j) {
        double vr = 0.0, vi = 0.0;
        for (int k = -KK; k <= KK; ++k) {
            double const ang = TWO_PI * static_cast<double>(k) * static_cast<double>(j) / Pd;
            double const c = std::cos(ang), s = std::sin(ang);
            double const ar = cr[static_cast<std::size_t>(k + KK)], ai = ci[static_cast<std::size_t>(k + KK)];
            vr += ar * c - ai * s;
            vi += ar * s + ai * c;
        }
        sx[j] = static_cast<Operon::Scalar>(vr); sy[j] = static_cast<Operon::Scalar>(vi);
    }
    return true;
}

inline auto gaussianLowPass(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                            bool closed, double width, int n,
                            std::vector<Operon::Scalar>& sx, std::vector<Operon::Scalar>& sy) -> bool
{
    sx.clear(); sy.clear();
    if (x.size() < 3 || y.size() < 3 || n < 8 || !(width > 0.0) || !std::isfinite(width)) { return false; }
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) { return false; }
    }
    std::vector<Operon::Scalar> rx, ry;
    arcResample(x, y, n, rx, ry);
    int const half = std::max(1, static_cast<int>(std::ceil(4.0 * width)));
    std::vector<double> w(static_cast<std::size_t>(2 * half + 1));
    for (int t = -half; t <= half; ++t) {
        double const u = static_cast<double>(t) / width;
        w[static_cast<std::size_t>(t + half)] = std::exp(-0.5 * u * u);
    }
    sx.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    sy.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    for (int j = 0; j < n; ++j) {
        double ax = 0.0, ay = 0.0, aw = 0.0;
        for (int t = -half; t <= half; ++t) {
            int k = j + t;
            if (closed) { k = ((k % n) + n) % n; }
            else if (k < 0 || k >= n) { continue; }
            double const wt = w[static_cast<std::size_t>(t + half)];
            ax += wt * static_cast<double>(rx[static_cast<std::size_t>(k)]);
            ay += wt * static_cast<double>(ry[static_cast<std::size_t>(k)]);
            aw += wt;
        }
        if (!(aw > 0.0)) { return false; }
        sx[static_cast<std::size_t>(j)] = static_cast<Operon::Scalar>(ax / aw);
        sy[static_cast<std::size_t>(j)] = static_cast<Operon::Scalar>(ay / aw);
    }
    return true;
}

inline auto gaussianWidthForK(int K, int n) -> double
{
    if (K < 1 || n < 8) { return 0.0; }
    return 0.1325037 * static_cast<double>(n) / static_cast<double>(K);
}

inline auto gaussianNoiseGain(double sigma) -> double
{
    if (!(sigma > 0.0) || !std::isfinite(sigma)) { return 1.0; }
    int const half = std::max(1, static_cast<int>(std::ceil(4.0 * sigma)));
    double s1 = 0.0, s2 = 0.0;
    for (int t = -half; t <= half; ++t) {
        double const u = static_cast<double>(t) / sigma;
        double const w = std::exp(-0.5 * u * u);
        s1 += w; s2 += w * w;
    }
    if (!(s1 > 0.0)) { return 1.0; }
    return s2 / (s1 * s1);
}

inline auto gaussianWidthForKEqualNoise(int K, int n) -> double
{
    if (K < 1 || n < 8) { return 0.0; }
    double const target = static_cast<double>(2 * K + 1) / static_cast<double>(n);
    if (!(target < 1.0)) { return 0.0; }
    thread_local int lastK = 0, lastN = 0;
    thread_local double lastW = 0.0;
    if (K == lastK && n == lastN) { return lastW; }
    constexpr double kSqrtPi = 1.7724538509055159;
    double const cont = static_cast<double>(n)
                      / (2.0 * kSqrtPi * static_cast<double>(2 * K + 1));
    double lo = 0.1, hi = 1024.0;
    if (cont > 0.2) {
        double const l = 0.5 * cont, h = 2.0 * cont;
        if (gaussianNoiseGain(l) > target && gaussianNoiseGain(h) < target) { lo = l; hi = h; }
    }
    for (int it = 0; it < 200; ++it) {
        double const mid = 0.5 * (lo + hi);
        if (gaussianNoiseGain(mid) > target) { lo = mid; } else { hi = mid; }
    }
    lastK = K; lastN = n; lastW = 0.5 * (lo + hi);
    return lastW;
}

inline auto savgolWeightRows(int window, int order, std::vector<double>& rows) -> bool
{
    rows.clear();
    if (window < 3 || (window % 2) == 0 || order < 0 || order >= window) { return false; }
    int const m = (window - 1) / 2;
    int const P = order + 1;
    std::vector<double> M(static_cast<std::size_t>(P * P), 0.0);
    std::vector<double> Minv(static_cast<std::size_t>(P * P), 0.0);
    for (int a = 0; a < P; ++a) {
        Minv[static_cast<std::size_t>(a * P + a)] = 1.0;
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int t = -m; t <= m; ++t) {
                double v = 1.0;
                for (int e = 0; e < a + b; ++e) { v *= static_cast<double>(t); }
                s += v;
            }
            M[static_cast<std::size_t>(a * P + b)] = s;
        }
    }
    for (int c = 0; c < P; ++c) {
        int piv = c;
        for (int r = c + 1; r < P; ++r) {
            if (std::abs(M[static_cast<std::size_t>(r * P + c)])
                > std::abs(M[static_cast<std::size_t>(piv * P + c)])) { piv = r; }
        }
        double const d = M[static_cast<std::size_t>(piv * P + c)];
        if (!(std::abs(d) > 0.0) || !std::isfinite(d)) { return false; }
        if (piv != c) {
            for (int k = 0; k < P; ++k) {
                std::swap(M[static_cast<std::size_t>(c * P + k)], M[static_cast<std::size_t>(piv * P + k)]);
                std::swap(Minv[static_cast<std::size_t>(c * P + k)], Minv[static_cast<std::size_t>(piv * P + k)]);
            }
        }
        double const inv = 1.0 / M[static_cast<std::size_t>(c * P + c)];
        for (int k = 0; k < P; ++k) {
            M[static_cast<std::size_t>(c * P + k)] *= inv;
            Minv[static_cast<std::size_t>(c * P + k)] *= inv;
        }
        for (int r = 0; r < P; ++r) {
            if (r == c) { continue; }
            double const f = M[static_cast<std::size_t>(r * P + c)];
            if (f == 0.0) { continue; }
            for (int k = 0; k < P; ++k) {
                M[static_cast<std::size_t>(r * P + k)] -= f * M[static_cast<std::size_t>(c * P + k)];
                Minv[static_cast<std::size_t>(r * P + k)] -= f * Minv[static_cast<std::size_t>(c * P + k)];
            }
        }
    }
    std::vector<double> tp(static_cast<std::size_t>(window * P));
    for (int i = 0; i < window; ++i) {
        double v = 1.0;
        for (int b = 0; b < P; ++b) {
            tp[static_cast<std::size_t>(i * P + b)] = v;
            v *= static_cast<double>(i - m);
        }
    }
    rows.assign(static_cast<std::size_t>(window * window), 0.0);
    for (int r = 0; r < window; ++r) {
        double const d = static_cast<double>(r - m);
        std::vector<double> g(static_cast<std::size_t>(P), 0.0);
        double da = 1.0;
        for (int a = 0; a < P; ++a) {
            for (int b = 0; b < P; ++b) {
                g[static_cast<std::size_t>(b)] += da * Minv[static_cast<std::size_t>(a * P + b)];
            }
            da *= d;
        }
        for (int i = 0; i < window; ++i) {
            double s = 0.0;
            for (int b = 0; b < P; ++b) {
                s += g[static_cast<std::size_t>(b)] * tp[static_cast<std::size_t>(i * P + b)];
            }
            rows[static_cast<std::size_t>(r * window + i)] = s;
        }
    }
    return true;
}

inline auto savgolNoiseGain(int window, int order) -> double
{
    std::vector<double> rows;
    if (!savgolWeightRows(window, order, rows)) { return 1.0; }
    int const m = (window - 1) / 2;
    double s2 = 0.0;
    for (int i = 0; i < window; ++i) {
        double const w = rows[static_cast<std::size_t>(m * window + i)];
        s2 += w * w;
    }
    return s2;
}

inline auto savgolWindowForKEqualNoise(int K, int n, int order = 3) -> int
{
    if (K < 1 || n < 8 || order < 1) { return 0; }
    double const target = static_cast<double>(2 * K + 1) / static_cast<double>(n);
    if (!(target < 1.0)) { return 0; }
    double const cont = 9.0 * static_cast<double>(n) / (4.0 * static_cast<double>(2 * K + 1));
    thread_local int lastK = 0, lastN = 0, lastOrd = 0, lastW = 0;
    if (K == lastK && n == lastN && order == lastOrd) { return lastW; }
    int best = 0;
    double bestErr = 0.0;
    int const hi = std::min(n, static_cast<int>(cont) * 2 + 9);
    for (int W = 3; W <= hi; ++W) {
        if ((W % 2) == 0 || W <= order) { continue; }
        double const e = std::abs(savgolNoiseGain(W, order) - target);
        if (best == 0 || e < bestErr) { best = W; bestErr = e; }
    }
    lastK = K; lastN = n; lastOrd = order; lastW = best;
    return best;
}

inline auto savgolLowPass(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                          bool closed, int window, int order, int n,
                          std::vector<Operon::Scalar>& sx, std::vector<Operon::Scalar>& sy) -> bool
{
    sx.clear(); sy.clear();
    if (x.size() < 3 || y.size() < 3 || n < 8 || window < 3 || (window % 2) == 0
        || order < 0 || order >= window || n < window) { return false; }
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) { return false; }
    }
    thread_local int memoW = 0, memoO = -1;
    thread_local std::vector<double> memoRows;
    if (window != memoW || order != memoO) {
        if (!savgolWeightRows(window, order, memoRows)) { memoW = 0; memoO = -1; return false; }
        memoW = window; memoO = order;
    }
    std::vector<double> const& rows = memoRows;
    std::vector<Operon::Scalar> rx, ry;
    arcResample(x, y, n, rx, ry);
    int const m = (window - 1) / 2;
    sx.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    sy.assign(static_cast<std::size_t>(n), Operon::Scalar{0});
    for (int j = 0; j < n; ++j) {
        int base = j - m;
        int row = m;
        if (!closed) {
            if (j < m) { base = 0; row = j; }
            else if (j > n - 1 - m) { base = n - window; row = j - (n - window) ; }
        }
        double ax = 0.0, ay = 0.0;
        for (int i = 0; i < window; ++i) {
            int k = base + i;
            if (closed) { k = ((k % n) + n) % n; }
            double const wt = rows[static_cast<std::size_t>(row * window + i)];
            ax += wt * static_cast<double>(rx[static_cast<std::size_t>(k)]);
            ay += wt * static_cast<double>(ry[static_cast<std::size_t>(k)]);
        }
        sx[static_cast<std::size_t>(j)] = static_cast<Operon::Scalar>(ax);
        sy[static_cast<std::size_t>(j)] = static_cast<Operon::Scalar>(ay);
    }
    return true;
}

constexpr int kSavgolOrder = 3;

inline auto lowPassSelect(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                          bool closed, int K, int n, int smoother,
                          std::vector<Operon::Scalar>& sx, std::vector<Operon::Scalar>& sy) -> bool
{
    if (smoother == 1) { return gaussianLowPass(x, y, closed, gaussianWidthForK(K, n), n, sx, sy); }
    if (smoother == 2) { return gaussianLowPass(x, y, closed, gaussianWidthForKEqualNoise(K, n), n, sx, sy); }
    if (smoother == 3) {
        return savgolLowPass(x, y, closed, savgolWindowForKEqualNoise(K, n, kSavgolOrder),
                             kSavgolOrder, n, sx, sy);
    }
    return fourierLowPass(x, y, closed, K, n, sx, sy);
}

inline auto turningFourier(std::vector<Operon::Scalar> const& x, std::vector<Operon::Scalar> const& y,
                           bool closed, int K = 8, int n = 240, int smoother = 0) -> double
{
    std::vector<Operon::Scalar> sx, sy;
    if (!lowPassSelect(x, y, closed, K, n, smoother, sx, sy)) { return 0.0; }
    return polylineTurning(sx, sy, closed);
}

inline auto polylineTurningSigned(std::vector<Operon::Scalar> const& x,
                                  std::vector<Operon::Scalar> const& y, bool closed) -> double
{
    constexpr double TWO_PI = 6.283185307179586;
    constexpr double PI = 3.141592653589793;
    std::size_t const M = std::min(x.size(), y.size());
    if (M < 3) { return 0.0; }
    double total = 0.0, prev = 0.0;
    bool have = false;
    std::size_t const nseg = closed ? M : M - 1;
    for (std::size_t i = 0; i < nseg; ++i) {
        std::size_t const j = (i + 1) % M;
        double const dx = static_cast<double>(x[j]) - static_cast<double>(x[i]);
        double const dy = static_cast<double>(y[j]) - static_cast<double>(y[i]);
        if (dx == 0.0 && dy == 0.0) { continue; }
        double const ang = std::atan2(dy, dx);
        if (have) {
            double d = ang - prev;
            while (d > PI) { d -= TWO_PI; }
            while (d < -PI) { d += TWO_PI; }
            total += d;
        }
        prev = ang; have = true;
    }
    if (closed && have) {
        std::size_t i = 0;
        while (i < M) {
            std::size_t const j = (i + 1) % M;
            double const dx = static_cast<double>(x[j]) - static_cast<double>(x[i]);
            double const dy = static_cast<double>(y[j]) - static_cast<double>(y[i]);
            if (dx != 0.0 || dy != 0.0) {
                double d = std::atan2(dy, dx) - prev;
                while (d > PI) { d -= TWO_PI; }
                while (d < -PI) { d += TWO_PI; }
                total += d;
                break;
            }
            ++i;
        }
    }
    return total / TWO_PI;
}

inline auto turningFourierSigned(std::vector<Operon::Scalar> const& x,
                                 std::vector<Operon::Scalar> const& y,
                                 bool closed, int K = 8, int n = 240, int smoother = 0) -> double
{
    std::vector<Operon::Scalar> sx, sy;
    if (!lowPassSelect(x, y, closed, K, n, smoother, sx, sy)) { return 0.0; }
    return polylineTurningSigned(sx, sy, closed);
}

inline auto bandedDTW(std::vector<Operon::Scalar> const& px, std::vector<Operon::Scalar> const& py,
                      std::vector<Operon::Scalar> const& qx, std::vector<Operon::Scalar> const& qy,
                      int w) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max() / 4;
    int const n = static_cast<int>(px.size());
    int const m = static_cast<int>(qx.size());
    thread_local std::vector<Operon::Scalar> prev, cur;
    prev.assign(static_cast<std::size_t>(m + 1), BIG);
    cur.assign(static_cast<std::size_t>(m + 1), BIG);
    prev[0] = 0;
    for (int i = 1; i <= n; ++i) {
        std::fill(cur.begin(), cur.end(), BIG);
        int jlo = std::max(1, i - w), jhi = std::min(m, i + w);
        auto pix = px[static_cast<std::size_t>(i-1)], piy = py[static_cast<std::size_t>(i-1)];
        for (int j = jlo; j <= jhi; ++j) {
            auto dx = qx[static_cast<std::size_t>(j-1)] - pix, dy = qy[static_cast<std::size_t>(j-1)] - piy;
            auto c = dx*dx + dy*dy;
            auto best = prev[static_cast<std::size_t>(j)];
            best = std::min(best, cur[static_cast<std::size_t>(j-1)]);
            best = std::min(best, prev[static_cast<std::size_t>(j-1)]);
            cur[static_cast<std::size_t>(j)] = c + best;
        }
        std::swap(prev, cur);
    }
    return prev[static_cast<std::size_t>(m)] / static_cast<Operon::Scalar>(n + m);
}

inline auto bestShift(std::vector<Operon::Scalar> const& px, std::vector<Operon::Scalar> const& py,
                      std::vector<Operon::Scalar> const& qx, std::vector<Operon::Scalar> const& qy) -> int
{
    int const n = static_cast<int>(px.size());
    Operon::Scalar best = -std::numeric_limits<Operon::Scalar>::max(); int bs = 0;
    thread_local std::vector<Operon::Scalar> qx2, qy2;
    qx2.resize(static_cast<std::size_t>(2 * n));
    qy2.resize(static_cast<std::size_t>(2 * n));
    std::copy(qx.begin(), qx.end(), qx2.begin());
    std::copy(qx.begin(), qx.end(), qx2.begin() + n);
    std::copy(qy.begin(), qy.end(), qy2.begin());
    std::copy(qy.begin(), qy.end(), qy2.begin() + n);

    for (int s = 0; s < n; ++s) {
        Operon::Scalar dot = 0;
        float const* ppx = px.data();
        float const* ppy = py.data();
        float const* pqx = qx2.data() + s;
        float const* pqy = qy2.data() + s;
#if defined(__ARM_NEON) || defined(__aarch64__)
        float32x4_t vdot = vdupq_n_f32(0.0f);
        for (int k = 0; k < n; k += 4) {
            float32x4_t vx1 = vld1q_f32(ppx + k);
            float32x4_t vy1 = vld1q_f32(ppy + k);
            float32x4_t vx2 = vld1q_f32(pqx + k);
            float32x4_t vy2 = vld1q_f32(pqy + k);
            vdot = vfmaq_f32(vdot, vx1, vx2);
            vdot = vfmaq_f32(vdot, vy1, vy2);
        }
        dot = vaddvq_f32(vdot);
#else
        for (int k = 0; k < n; ++k) {
            dot += ppx[k] * pqx[k] + ppy[k] * pqy[k];
        }
#endif
        if (dot > best) { best = dot; bs = s; }
    }
    return bs;
}

inline auto arcCorrespondenceDistance(std::vector<Operon::Scalar> const& px, std::vector<Operon::Scalar> const& py,
                                      std::vector<Operon::Scalar> const& qx, std::vector<Operon::Scalar> const& qy)
    -> Operon::Scalar
{
    auto const n = px.size();
    if (n == 0 || qx.size() != n) { return std::numeric_limits<Operon::Scalar>::max(); }
    Operon::Scalar acc = 0;
    for (std::size_t k = 0; k < n; ++k) {
        auto const dx = px[k] - qx[k], dy = py[k] - qy[k];
        acc += dx * dx + dy * dy;
    }
    return acc / static_cast<Operon::Scalar>(n);
}

inline auto orderDistanceArc(std::vector<Operon::Scalar> const& gridX, std::vector<Operon::Scalar> const& gridY,
                             std::vector<Operon::Scalar> const& dataX, std::vector<Operon::Scalar> const& dataY,
                             bool closed, int n = 64,
                             std::vector<Operon::Scalar> const* precomputedDataArcX = nullptr,
                             std::vector<Operon::Scalar> const* precomputedDataArcY = nullptr) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
    if (gridX.size() < 2 || dataX.size() < 2) { return BIG; }
    for (auto v : gridX) { if (!std::isfinite(v)) { return BIG; } }
    for (auto v : gridY) { if (!std::isfinite(v)) { return BIG; } }
    thread_local std::vector<Operon::Scalar> cx, cy;
    cx.assign(gridX.begin(), gridX.end());
    cy.assign(gridY.begin(), gridY.end());
    normalizeInPlace(cx, cy);
    thread_local std::vector<Operon::Scalar> pxScratch, pyScratch, qx, qy;
    std::vector<Operon::Scalar> const* pxPtr = nullptr;
    std::vector<Operon::Scalar> const* pyPtr = nullptr;
    if (precomputedDataArcX != nullptr && precomputedDataArcY != nullptr &&
        precomputedDataArcX->size() == static_cast<std::size_t>(n) &&
        precomputedDataArcY->size() == static_cast<std::size_t>(n)) {
        pxPtr = precomputedDataArcX;
        pyPtr = precomputedDataArcY;
    } else {
        arcResample(dataX, dataY, n, pxScratch, pyScratch);
        pxPtr = &pxScratch;
        pyPtr = &pyScratch;
    }
    auto const& px = *pxPtr;
    auto const& py = *pyPtr;
    arcResample(cx, cy, n, qx, qy);
    Operon::Scalar best = BIG;
    thread_local std::vector<Operon::Scalar> rx, ry;
    thread_local std::vector<Operon::Scalar> sx, sy;
    for (int dir = 0; dir < 2; ++dir) {
        rx.assign(qx.begin(), qx.end());
        ry.assign(qy.begin(), qy.end());
        if (dir == 1) { std::reverse(rx.begin(), rx.end()); std::reverse(ry.begin(), ry.end()); }
        if (closed) {
            int bs = bestShift(px, py, rx, ry);
            for (int off = bs - 1; off <= bs + 1; ++off) {
                int r = ((off % n) + n) % n;
                sx.resize(static_cast<std::size_t>(n));
                sy.resize(static_cast<std::size_t>(n));
                for (int k = 0; k < n; ++k) { sx[static_cast<std::size_t>(k)] = rx[static_cast<std::size_t>((k + r) % n)];
                                              sy[static_cast<std::size_t>(k)] = ry[static_cast<std::size_t>((k + r) % n)]; }
                best = std::min(best, arcCorrespondenceDistance(px, py, sx, sy));
            }
        } else {
            best = std::min(best, arcCorrespondenceDistance(px, py, rx, ry));
        }
    }
    return best;
}

inline auto orderDistanceBanded(std::vector<Operon::Scalar> const& gridX, std::vector<Operon::Scalar> const& gridY,
                                std::vector<Operon::Scalar> const& dataX, std::vector<Operon::Scalar> const& dataY,
                                bool closed, int band, int n = 64,
                                std::vector<Operon::Scalar> const* precomputedDataArcX = nullptr,
                                std::vector<Operon::Scalar> const* precomputedDataArcY = nullptr) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
    if (gridX.size() < 2 || dataX.size() < 2) { return BIG; }
    for (auto v : gridX) { if (!std::isfinite(v)) { return BIG; } }
    for (auto v : gridY) { if (!std::isfinite(v)) { return BIG; } }
    thread_local std::vector<Operon::Scalar> cx, cy;
    cx.assign(gridX.begin(), gridX.end());
    cy.assign(gridY.begin(), gridY.end());
    normalizeInPlace(cx, cy);
    thread_local std::vector<Operon::Scalar> pxScratch, pyScratch, qx, qy;
    std::vector<Operon::Scalar> const* pxPtr = nullptr;
    std::vector<Operon::Scalar> const* pyPtr = nullptr;
    if (precomputedDataArcX != nullptr && precomputedDataArcY != nullptr &&
        precomputedDataArcX->size() == static_cast<std::size_t>(n) &&
        precomputedDataArcY->size() == static_cast<std::size_t>(n)) {
        pxPtr = precomputedDataArcX;
        pyPtr = precomputedDataArcY;
    } else {
        arcResample(dataX, dataY, n, pxScratch, pyScratch);
        pxPtr = &pxScratch;
        pyPtr = &pyScratch;
    }
    auto const& px = *pxPtr;
    auto const& py = *pyPtr;
    arcResample(cx, cy, n, qx, qy);
    int const w = std::max(1, std::min(band, n));
    Operon::Scalar best = BIG;
    thread_local std::vector<Operon::Scalar> rx, ry;
    thread_local std::vector<Operon::Scalar> sx, sy;
    for (int dir = 0; dir < 2; ++dir) {
        rx.assign(qx.begin(), qx.end());
        ry.assign(qy.begin(), qy.end());
        if (dir == 1) { std::reverse(rx.begin(), rx.end()); std::reverse(ry.begin(), ry.end()); }
        if (closed) {
            int bs = bestShift(px, py, rx, ry);
            for (int off = bs - 1; off <= bs + 1; ++off) {
                int r = ((off % n) + n) % n;
                sx.resize(static_cast<std::size_t>(n));
                sy.resize(static_cast<std::size_t>(n));
                for (int k = 0; k < n; ++k) { sx[static_cast<std::size_t>(k)] = rx[static_cast<std::size_t>((k + r) % n)];
                                              sy[static_cast<std::size_t>(k)] = ry[static_cast<std::size_t>((k + r) % n)]; }
                best = std::min(best, bandedDTW(px, py, sx, sy, w));
            }
        } else {
            best = std::min(best, bandedDTW(px, py, rx, ry, w));
        }
    }
    return best;
}

inline auto orderDistance(std::vector<Operon::Scalar> const& gridX, std::vector<Operon::Scalar> const& gridY,
                          std::vector<Operon::Scalar> const& dataX, std::vector<Operon::Scalar> const& dataY,
                          bool closed, int n = 64,
                          std::vector<Operon::Scalar> const* precomputedDataArcX = nullptr,
                          std::vector<Operon::Scalar> const* precomputedDataArcY = nullptr) -> Operon::Scalar
{
    constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
    if (gridX.size() < 2 || dataX.size() < 2) { return BIG; }
    for (auto v : gridX) { if (!std::isfinite(v)) { return BIG; } }
    for (auto v : gridY) { if (!std::isfinite(v)) { return BIG; } }
    thread_local std::vector<Operon::Scalar> cx, cy;
    cx.assign(gridX.begin(), gridX.end());
    cy.assign(gridY.begin(), gridY.end());
    normalizeInPlace(cx, cy);
    thread_local std::vector<Operon::Scalar> pxScratch, pyScratch, qx, qy;
    std::vector<Operon::Scalar> const* pxPtr = nullptr;
    std::vector<Operon::Scalar> const* pyPtr = nullptr;
    if (precomputedDataArcX != nullptr && precomputedDataArcY != nullptr &&
        precomputedDataArcX->size() == static_cast<std::size_t>(n) &&
        precomputedDataArcY->size() == static_cast<std::size_t>(n)) {
        pxPtr = precomputedDataArcX;
        pyPtr = precomputedDataArcY;
    } else {
        arcResample(dataX, dataY, n, pxScratch, pyScratch);
        pxPtr = &pxScratch;
        pyPtr = &pyScratch;
    }
    auto const& px = *pxPtr;
    auto const& py = *pyPtr;
    arcResample(cx, cy, n, qx, qy);
    int w = n;
    Operon::Scalar best = BIG;
    thread_local std::vector<Operon::Scalar> rx, ry;
    thread_local std::vector<Operon::Scalar> sx, sy;
    for (int dir = 0; dir < 2; ++dir) {
        rx.assign(qx.begin(), qx.end());
        ry.assign(qy.begin(), qy.end());
        if (dir == 1) { std::reverse(rx.begin(), rx.end()); std::reverse(ry.begin(), ry.end()); }
        if (closed) {
            int bs = bestShift(px, py, rx, ry);
            for (int off = bs - 1; off <= bs + 1; ++off) {
                int r = ((off % n) + n) % n;
                sx.resize(static_cast<std::size_t>(n));
                sy.resize(static_cast<std::size_t>(n));
                for (int k = 0; k < n; ++k) { sx[static_cast<std::size_t>(k)] = rx[static_cast<std::size_t>((k + r) % n)];
                                              sy[static_cast<std::size_t>(k)] = ry[static_cast<std::size_t>((k + r) % n)]; }
                best = std::min(best, bandedDTW(px, py, sx, sy, w));
            }
        } else {
            best = std::min(best, bandedDTW(px, py, rx, ry, w));
        }
    }
    return best;
}

}
#endif
