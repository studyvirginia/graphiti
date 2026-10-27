#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <vector>

#include "sketch/net_band_select.hpp"

namespace Sketch {

inline constexpr int kCisExactModeDefault = 4;
inline auto cisExactModeRef() -> std::atomic<int>&
{
    static std::atomic<int> v{[] {
        char const* e = std::getenv("SKETCH_CIS_EXACT");
        return (e != nullptr) ? std::atoi(e) : kCisExactModeDefault;
    }()};
    return v;
}
inline auto cisExactModeDefault() -> int { return cisExactModeRef().load(std::memory_order_relaxed); }
inline void setCisExactMode(int mode)
{
    if (mode < 0) {
        char const* e = std::getenv("SKETCH_CIS_EXACT");
        mode = (e != nullptr) ? std::atoi(e) : kCisExactModeDefault;
    }
    cisExactModeRef().store(mode, std::memory_order_relaxed);
}

struct CisScoringConfig {
    bool useCisScoring{false};
    bool scoreDcCell{true};
    bool useSavitzkyGolayBoundary{false};
    bool useScaleInvariantFitness{false};
    int exactLikelihood{cisExactModeDefault()};
    int residualPrecision{0};
};

struct CisCell {
    int level{0};
    std::size_t index{0};
    double energy{0.0};
    double predictedEnergy{0.0};
    double predictedLogVar{0.0};
    bool hasEpi{false};
    double dof{1.0};
    double nll{0.0};
    double recordDelta{0.0};
    bool intervalCell{false};
};

struct CisScoringResult {
    double fitness{0.0};
    double shapeNll{0.0};
    double scaleNll{0.0};
    double sHat{0.0};
    std::vector<CisCell> cells;
    bool ok{true};
};

struct CisPerCellEnergyInput {
    std::vector<double> perCellEnergy;
    bool present{false};
    std::vector<std::vector<double>> memberEnergy;
};

inline void applyDirectCellEnergyStub(std::vector<CisCell>& cells, CisPerCellEnergyInput const& direct)
{
    if (!direct.present || direct.perCellEnergy.size() != cells.size()) { return; }
    for (std::size_t i = 0; i < cells.size(); ++i) {
        cells[i].predictedEnergy = std::max(direct.perCellEnergy[i], 0.0);
    }
}

inline constexpr int kCisGrid48Default = 7;
inline auto cisGrid48ModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_CIS_GRID48");
        return (e != nullptr) ? std::atoi(e) : kCisGrid48Default;
    }();
    return v;
}

inline constexpr int kB3Grid48GpuDefault = 1;
inline auto b3Grid48GpuFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_B3_GRID48_GPU");
        return (e != nullptr) ? std::atoi(e) : kB3Grid48GpuDefault;
    }();
    return v;
}

struct CisGrid48Counters {
    std::atomic<long long> badStroke{0};
    std::atomic<long long> missingModel{0};
    std::atomic<long long> cellMismatch{0};
    std::atomic<long long> residualFail{0};
    std::atomic<long long> bit1NoBit0{0};
    std::atomic<long long> ensembleSkipped{0};
    std::atomic<long long> absFallback{0};
    std::atomic<long long> law1FallbackToLaw3{0};
};
inline auto cisGrid48Counters() -> CisGrid48Counters&
{
    static CisGrid48Counters c;
    return c;
}

inline constexpr int kCisEnsembleDefault = 1;
inline auto cisEnsembleModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_CIS_ENSEMBLE");
        return (e != nullptr) ? std::atoi(e) : kCisEnsembleDefault;
    }();
    return v;
}

inline constexpr int kNmlAbsDefault = 0;
inline auto nmlAbsModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_NML_ABS");
        return (e != nullptr) ? std::atoi(e) : kNmlAbsDefault;
    }();
    return v;
}

namespace detail {

constexpr double kCisPi = 3.14159265358979323846;
constexpr double kCisEnergyEps = 1e-300;

inline void cisHaarLevel(std::vector<double> const& x, std::vector<double>& approx, std::vector<double>& det)
{
    std::size_t const n = x.size();
    std::size_t const pairs = n / 2;
    approx.assign(pairs + (n % 2), 0.0);
    det.assign(pairs, 0.0);
    constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
    for (std::size_t i = 0; i < pairs; ++i) {
        double const a = x[2 * i];
        double const b = x[2 * i + 1];
        approx[i] = (a + b) * kInvSqrt2;
        det[i] = (a - b) * kInvSqrt2;
    }
    if (n % 2 == 1) {
        approx[pairs] = x[n - 1];
    }
}

inline std::vector<double> cisHaarLevelInverse(std::vector<double> const& approx, std::vector<double> const& det, std::size_t n)
{
    std::size_t const pairs = n / 2;
    std::vector<double> x(n, 0.0);
    constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
    for (std::size_t i = 0; i < pairs; ++i) {
        double const s = approx[i];
        double const d = det[i];
        x[2 * i] = (s + d) * kInvSqrt2;
        x[2 * i + 1] = (s - d) * kInvSqrt2;
    }
    if (n % 2 == 1) {
        x[n - 1] = approx[pairs];
    }
    return x;
}

inline double cisTrigamma(double x)
{
    double result = 0.0;
    while (x < 6.0) {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    double const xi = 1.0 / x;
    double const xi2 = xi * xi;
    result += xi + 0.5 * xi2
        + xi2 * xi * (1.0 / 6.0 - xi2 * (1.0 / 30.0 - xi2 * (1.0 / 42.0 - xi2 * (1.0 / 30.0))));
    return result;
}

inline int cisNearestSupportedK(std::size_t count)
{
    static constexpr int kSupported[] = {1, 2, 4, 8, 16};
    int best = kSupported[0];
    long long bestDiff = std::numeric_limits<long long>::max();
    for (int k : kSupported) {
        long long const diff = std::llabs(static_cast<long long>(count) - static_cast<long long>(k));
        if (diff < bestDiff) {
            bestDiff = diff;
            best = k;
        }
    }
    return best;
}

constexpr int kSgPolyOrder = 2;
constexpr int kSgWindowSize = kSgPolyOrder + 1;

inline auto cisSolveSmallLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b, bool& ok) -> std::vector<double>
{
    std::size_t const n = b.size();
    ok = true;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        double best = std::fabs(A[col][col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
        }
        if (best < 1e-14) { ok = false; return std::vector<double>(n, 0.0); }
        if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }
        for (std::size_t r = col + 1; r < n; ++r) {
            double const f = A[r][col] / A[col][col];
            for (std::size_t c = col; c < n; ++c) { A[r][c] -= f * A[col][c]; }
            b[r] -= f * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (std::size_t ri = n; ri-- > 0;) {
        double s = b[ri];
        for (std::size_t c = ri + 1; c < n; ++c) { s -= A[ri][c] * x[c]; }
        x[ri] = s / A[ri][ri];
    }
    return x;
}

inline auto cisSgOneSidedCoeffs(int windowSize, int order, double evalPos, bool& ok) -> std::vector<double>
{
    int const m = order + 1;
    ok = (windowSize >= m);
    if (!ok) { return std::vector<double>(static_cast<std::size_t>(std::max(windowSize, 0)), 0.0); }
    std::vector<std::vector<double>> V(static_cast<std::size_t>(windowSize), std::vector<double>(static_cast<std::size_t>(m)));
    for (int i = 0; i < windowSize; ++i) {
        double p = 1.0;
        for (int j = 0; j < m; ++j) { V[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = p; p *= static_cast<double>(i); }
    }
    std::vector<std::vector<double>> A(static_cast<std::size_t>(m), std::vector<double>(static_cast<std::size_t>(m), 0.0));
    for (int a = 0; a < m; ++a) {
        for (int b = 0; b < m; ++b) {
            double s = 0.0;
            for (int i = 0; i < windowSize; ++i) { s += V[static_cast<std::size_t>(i)][static_cast<std::size_t>(a)] * V[static_cast<std::size_t>(i)][static_cast<std::size_t>(b)]; }
            A[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = s;
        }
    }
    std::vector<double> e(static_cast<std::size_t>(m));
    { double p = 1.0; for (int j = 0; j < m; ++j) { e[static_cast<std::size_t>(j)] = p; p *= evalPos; } }
    bool solveOk = true;
    std::vector<double> const coefInBasis = cisSolveSmallLinearSystem(A, e, solveOk);
    ok = solveOk;
    std::vector<double> c(static_cast<std::size_t>(windowSize), 0.0);
    if (!ok) { return c; }
    for (int i = 0; i < windowSize; ++i) {
        double s = 0.0;
        for (int j = 0; j < m; ++j) { s += V[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * coefInBasis[static_cast<std::size_t>(j)]; }
        c[static_cast<std::size_t>(i)] = s;
    }
    return c;
}

inline void cisHaarLevelSgBoundary(std::vector<double> const& x, std::vector<double>& approx, std::vector<double>& det)
{
    cisHaarLevel(x, approx, det);
    std::size_t const n = x.size();
    if (n % 2 != 1) { return; }
    if (n < static_cast<std::size_t>(kSgWindowSize) + 1) { return; }

    std::size_t const boundaryIdx = n - 1;
    std::size_t const windowStart = boundaryIdx - static_cast<std::size_t>(kSgWindowSize);
    std::vector<double> const window(x.begin() + static_cast<std::ptrdiff_t>(windowStart), x.begin() + static_cast<std::ptrdiff_t>(boundaryIdx));
    bool ok = true;
    auto const coeffs = cisSgOneSidedCoeffs(kSgWindowSize, kSgPolyOrder, static_cast<double>(kSgWindowSize), ok);
    if (!ok) { return; }
    double pred = 0.0;
    for (int i = 0; i < kSgWindowSize; ++i) { pred += coeffs[static_cast<std::size_t>(i)] * window[static_cast<std::size_t>(i)]; }
    if (!std::isfinite(pred)) { return; }
    double const boundaryDetail = x[boundaryIdx] - pred;
    if (!std::isfinite(boundaryDetail)) { return; }

    approx.back() = pred;
    det.push_back(boundaryDetail);
}

}

struct HaarDecomposition {
    std::vector<std::vector<double>> details;
    std::vector<double> approxFinal;
    std::vector<std::size_t> lengthsAtLevel;
};

inline auto haarForward(std::vector<double> const& residual) -> HaarDecomposition
{
    HaarDecomposition out;
    std::vector<double> cur = residual;
    while (cur.size() >= 2) {
        out.lengthsAtLevel.push_back(cur.size());
        std::vector<double> approx;
        std::vector<double> det;
        detail::cisHaarLevel(cur, approx, det);
        out.details.push_back(std::move(det));
        cur = std::move(approx);
    }
    out.approxFinal = cur;
    return out;
}

inline auto haarForwardSgBoundary(std::vector<double> const& residual) -> HaarDecomposition
{
    HaarDecomposition out;
    std::vector<double> cur = residual;
    bool first = true;
    while (cur.size() >= 2) {
        out.lengthsAtLevel.push_back(cur.size());
        std::vector<double> approx;
        std::vector<double> det;
        if (first) {
            detail::cisHaarLevelSgBoundary(cur, approx, det);
            first = false;
        } else {
            detail::cisHaarLevel(cur, approx, det);
        }
        out.details.push_back(std::move(det));
        cur = std::move(approx);
    }
    out.approxFinal = cur;
    return out;
}

inline auto haarInverse(HaarDecomposition const& dec) -> std::vector<double>
{
    std::vector<double> cur = dec.approxFinal;
    for (std::size_t li = dec.details.size(); li-- > 0;) {
        cur = detail::cisHaarLevelInverse(cur, dec.details[li], dec.lengthsAtLevel[li]);
    }
    return cur;
}

inline auto haarFlattenCoeffs(HaarDecomposition const& dec) -> std::vector<double>
{
    std::vector<double> flat;
    std::size_t total = dec.approxFinal.size();
    for (auto const& d : dec.details) { total += d.size(); }
    flat.reserve(total);
    for (auto const& d : dec.details) { flat.insert(flat.end(), d.begin(), d.end()); }
    flat.insert(flat.end(), dec.approxFinal.begin(), dec.approxFinal.end());
    return flat;
}

inline auto cisLognormalNllOfRatio(double ratio, double dof, double extraVar = 0.0) -> double
{
    double const dofSafe = std::max(dof, 1e-6);
    double var = detail::cisTrigamma(dofSafe / 2.0) + std::max(extraVar, 0.0);
    var = std::max(var, 1e-12);
    double const L = std::log(std::max(ratio, detail::kCisEnergyEps));
    return 0.5 * (L * L) / var + 0.5 * std::log(2.0 * detail::kCisPi * var);
}

inline auto cisExactCellNllData(double y, double dof, double logSp) -> double
{
    return 0.5 * y + 0.5 * dof * (std::log(2.0 * detail::kCisPi) + logSp);
}
inline auto cisExactCellNllLog(double y, double L, double dof) -> double
{
    if (L == -std::numeric_limits<double>::infinity()) { return std::numeric_limits<double>::infinity(); }
    return 0.5 * y - 0.5 * dof * L + 0.5 * dof * std::log(2.0) + std::lgamma(0.5 * dof);
}

namespace detail {
inline auto cisUlpF32Half(double r) -> double
{
    float const f = std::abs(static_cast<float>(r));
    if (!std::isfinite(f)) { return std::numeric_limits<double>::infinity(); }
    float const nx = std::nextafter(f, std::numeric_limits<float>::infinity());
    return 0.5 * (static_cast<double>(nx) - static_cast<double>(f));
}

inline auto cisUlpF16Half(double mag) -> double
{
    if (!std::isfinite(mag)) { return std::numeric_limits<double>::infinity(); }
    double const m = std::abs(mag);
    if (m == 0.0) { return 0.5 * std::ldexp(1.0, -24); }
    int ex = 0;
    std::frexp(m, &ex);
    int const exp2 = ex - 1;
    double const ulp = (exp2 >= -14) ? std::ldexp(1.0, exp2 - 10) : std::ldexp(1.0, -24);
    return 0.5 * ulp;
}

inline auto cisFp16Mode1FallbackHalfWidth() -> double { return cisUlpF16Half(1.0); }
inline auto cisLogRegLowerGamma(double a, double x) -> double
{
    constexpr double kInf = std::numeric_limits<double>::infinity();
    if (!(x > 0.0)) { return -kInf; }
    if (x == kInf) { return 0.0; }
    if (a == 0.5) {
        double const z = std::sqrt(x);
        return (x < 1.0) ? std::log(std::erf(z)) : std::log1p(-std::erfc(z));
    }
    constexpr double kEps = std::numeric_limits<double>::epsilon();
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 0; n < 100000; ++n) {
            ap += 1.0; del *= x / ap; sum += del;
            if (std::abs(del) < std::abs(sum) * kEps) { break; }
        }
        return a * std::log(x) - x - std::lgamma(a) + std::log(sum);
    }
    double const tiny = std::numeric_limits<double>::min() / kEps;
    double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
    for (int i = 1; i < 100000; ++i) {
        double const an = -i * (i - a);
        b += 2.0;
        d = an * d + b; if (std::abs(d) < tiny) { d = tiny; }
        c = b + an / c; if (std::abs(c) < tiny) { c = tiny; }
        d = 1.0 / d;
        double const del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < kEps) { break; }
    }
    double const logQ = a * std::log(x) - x - std::lgamma(a) + std::log(h);
    return std::log1p(-std::exp(logQ));
}
}

inline auto cisRecordingHalfWidths(std::vector<double> const& residual, int precision = 0) -> HaarDecomposition
{
    std::vector<double> w(residual.size());
    if (precision == 2 || precision == 3) {
        for (std::size_t i = 0; i < residual.size(); ++i) { w[i] = detail::cisUlpF16Half(residual[i]); }
    } else if (precision == 1) {
        double const fallback = detail::cisFp16Mode1FallbackHalfWidth();
        for (std::size_t i = 0; i < residual.size(); ++i) { w[i] = fallback; }
    } else {
        for (std::size_t i = 0; i < residual.size(); ++i) { w[i] = detail::cisUlpF32Half(residual[i]); }
    }
    HaarDecomposition out;
    std::vector<double> cur = std::move(w);
    constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
    while (cur.size() >= 2) {
        out.lengthsAtLevel.push_back(cur.size());
        std::size_t const n = cur.size(), pairs = n / 2;
        std::vector<double> approx(pairs + (n % 2), 0.0), det(pairs, 0.0);
        for (std::size_t i = 0; i < pairs; ++i) {
            double const v = (cur[2 * i] + cur[2 * i + 1]) * kInvSqrt2;
            approx[i] = v; det[i] = v;
        }
        if (n % 2 == 1) { approx[pairs] = cur[n - 1]; }
        out.details.push_back(std::move(det));
        cur = std::move(approx);
    }
    out.approxFinal = cur;
    return out;
}

inline auto buildCisCells(
    HaarDecomposition const& dec,
    std::vector<double> const& sigma2,
    std::vector<double> const& logVar,
    bool scoreDcCell = true) -> std::vector<CisCell>
{
    std::vector<CisCell> cells;
    bool const haveLv = !logVar.empty() && logVar.size() == sigma2.size();

    auto pushLevel = [&](std::vector<double> const& coeffs, int level) {
        std::size_t const n = coeffs.size();
        if (n == 0) { return; }
        int const K = detail::cisNearestSupportedK(n);
        auto const s2K = Sketch::detail::aggToK(sigma2, K);
        std::vector<double> lvK;
        if (haveLv) { lvK = Sketch::detail::aggToK(logVar, K); }
        std::size_t const kCount = s2K.size();
        for (std::size_t i = 0; i < n; ++i) {
            CisCell cell;
            cell.level = level;
            cell.index = i;
            cell.energy = coeffs[i] * coeffs[i];
            std::size_t k = 0;
            if (kCount > 0) {
                k = std::min<std::size_t>((i * static_cast<std::size_t>(K)) / std::max<std::size_t>(n, 1), kCount - 1);
                cell.predictedEnergy = std::max(s2K[k], 0.0);
            }
            if (haveLv && k < lvK.size()) {
                cell.predictedLogVar = lvK[k];
                cell.hasEpi = true;
            }
            cell.dof = 1.0;
            cells.push_back(cell);
        }
    };

    for (std::size_t li = 0; li < dec.details.size(); ++li) {
        pushLevel(dec.details[li], static_cast<int>(li + 1));
    }
    if (scoreDcCell) {
        pushLevel(dec.approxFinal, static_cast<int>(dec.details.size() + 1));
    }

    return cells;
}

inline auto cisEpiExtraVar(CisCell const& cell) -> double
{
    if (!cell.hasEpi) { return 0.0; }
    return std::exp(std::clamp(cell.predictedLogVar, -8.0, 4.0));
}

inline void computeCisNllExact(std::vector<CisCell>& cells, int mode, double& shapeNll, double& scaleNll,
                               double& sHat, bool& ok)
{
    constexpr double kInf = std::numeric_limits<double>::infinity();
    ok = true;
    shapeNll = 0.0;
    scaleNll = 0.0;
    sHat = 0.0;
    std::size_t const n = cells.size();
    if (n == 0) { ok = false; return; }
    double S = 0.0, D = 0.0;
    for (auto const& c : cells) {
        double const p = c.predictedEnergy;
        double const e = c.energy;
        if (!(p > 0.0) || !std::isfinite(p) || !std::isfinite(e) || e < 0.0 || !(c.dof > 0.0)) { ok = false; return; }
        S += e / p;
        D += c.dof;
    }
    if (!std::isfinite(S)) { ok = false; return; }
    sHat = S / D;
    if (S == 0.0) {
        if (mode == 3) {
            for (auto& c : cells) { c.intervalCell = true; c.nll = 0.0; }
            shapeNll = 0.0;
            scaleNll = kInf;
            return;
        }
        for (auto& c : cells) { c.nll = (mode == 1) ? -kInf : kInf; }
        shapeNll = (mode == 1) ? -kInf : kInf;
        scaleNll = kInf;
        return;
    }
    double const logS = std::log(S) - std::log(D);
    for (auto& c : cells) {
        double const p = c.predictedEnergy;
        double const y = (c.energy / p) / sHat;
        double nll = 0.0;
        if (mode == 1) {
            nll = cisExactCellNllData(y, c.dof, logS + std::log(p));
        } else if (mode == 3 && c.energy < c.recordDelta) {
            c.intervalCell = true;
            nll = -detail::cisLogRegLowerGamma(0.5 * c.dof, (c.recordDelta / p) / sHat * 0.5);
        } else {
            double const L = (c.energy > 0.0) ? (std::log(c.energy) - std::log(p) - logS)
                                              : -kInf;
            nll = cisExactCellNllLog(y, L, c.dof);
        }
        c.nll = nll;
        shapeNll += nll;
    }
    scaleNll = 0.5 * D * (sHat - 1.0 - logS);
    if (std::isnan(shapeNll) || std::isnan(scaleNll)) { ok = false; }
}

inline double cisLseNegMean(std::vector<double> const& x)
{
    std::size_t const M = x.size();
    double mn = x[0];
    for (double v : x) { mn = std::min(mn, v); }
    if (mn == -std::numeric_limits<double>::infinity()) { return -std::numeric_limits<double>::infinity(); }
    if (!std::isfinite(mn)) { return mn; }
    double sum = 0.0;
    for (double v : x) {
        double const d = v - mn;
        sum += std::isfinite(d) ? std::exp(-d) : 0.0;
    }
    if (!(sum > 0.0)) { return std::numeric_limits<double>::infinity(); }
    return mn - std::log(sum) + std::log(static_cast<double>(M));
}

inline void computeCisNllEnsemble(std::vector<CisCell>& cells, int mode,
                                  std::vector<std::vector<double>> const& memberEnergy,
                                  double& shapeNll, double& scaleNll, double& sHat, bool& ok)
{
    ok = true;
    shapeNll = 0.0;
    scaleNll = 0.0;
    sHat = 0.0;
    std::size_t const n = cells.size();
    std::size_t const M = memberEnergy.size();
    if (n == 0 || M == 0) { ok = false; return; }
    for (auto const& me : memberEnergy) {
        if (me.size() != n) { ok = false; return; }
    }

    std::vector<double> savedPredicted(n);
    for (std::size_t i = 0; i < n; ++i) { savedPredicted[i] = cells[i].predictedEnergy; }
    std::vector<char> savedInterval(n);
    for (std::size_t i = 0; i < n; ++i) { savedInterval[i] = cells[i].intervalCell ? 1 : 0; }

    std::vector<double> shapeM(M), fitM(M), sHatM(M);
    std::vector<std::vector<double>> nllPerCellM(M, std::vector<double>(n));
    for (std::size_t m = 0; m < M; ++m) {
        for (std::size_t i = 0; i < n; ++i) { cells[i].predictedEnergy = std::max(memberEnergy[m][i], 0.0); }
        double sN = 0.0, scN = 0.0, sH = 0.0;
        bool mOk = true;
        computeCisNllExact(cells, mode, sN, scN, sH, mOk);
        if (!mOk) { ok = false; }
        shapeM[m] = sN;
        fitM[m] = sN + scN;
        sHatM[m] = sH;
        for (std::size_t i = 0; i < n; ++i) { nllPerCellM[m][i] = cells[i].nll; }
    }

    for (std::size_t i = 0; i < n; ++i) {
        cells[i].predictedEnergy = savedPredicted[i];
        cells[i].intervalCell = savedInterval[i] != 0;
    }
    if (!ok) { return; }

    shapeNll = cisLseNegMean(shapeM);
    double const total = cisLseNegMean(fitM);
    scaleNll = total - shapeNll;

    double mn = shapeM[0];
    for (double v : shapeM) { mn = std::min(mn, v); }
    std::vector<double> w(M, 0.0);
    double wsum = 0.0;
    for (std::size_t m = 0; m < M; ++m) {
        double const d = shapeM[m] - mn;
        w[m] = std::isfinite(d) ? std::exp(-d) : 0.0;
        wsum += w[m];
    }
    if (wsum > 0.0) {
        for (std::size_t m = 0; m < M; ++m) { w[m] /= wsum; }
    } else {
        for (std::size_t m = 0; m < M; ++m) { w[m] = 1.0 / static_cast<double>(M); }
    }
    for (std::size_t m = 0; m < M; ++m) { sHat += w[m] * sHatM[m]; }
    for (std::size_t i = 0; i < n; ++i) {
        double nll = 0.0;
        for (std::size_t m = 0; m < M; ++m) { nll += w[m] * nllPerCellM[m][i]; }
        cells[i].nll = nll;
    }
    if (std::isnan(shapeNll) || std::isnan(scaleNll)) { ok = false; }
}

inline auto cisJointCoarseFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("CIS_JOINT_COARSE");
        return (e != nullptr && e[0] == '1') ? 1 : 0;
    }();
    return v;
}

inline void computeCisNllEnsembleJointCoarse(std::vector<CisCell> const& cells,
                                             std::vector<std::vector<double>> const& memberEnergy,
                                             double& fitness, bool& ok)
{
    constexpr double kInf = std::numeric_limits<double>::infinity();
    ok = true;
    fitness = 0.0;
    std::size_t const n = cells.size();
    std::size_t const M = memberEnergy.size();
    if (n == 0 || M == 0) { ok = false; return; }
    for (auto const& me : memberEnergy) {
        if (me.size() != n) { ok = false; return; }
    }
    int maxLevel = 0;
    for (auto const& c : cells) { maxLevel = std::max(maxLevel, c.level); }
    std::vector<std::vector<std::size_t>> blocks(static_cast<std::size_t>(maxLevel) + 1);
    for (std::size_t c = 0; c < n; ++c) { blocks[static_cast<std::size_t>(cells[c].level)].push_back(c); }

    std::vector<double> totalM(M, 0.0);
    for (std::size_t m = 0; m < M; ++m) {
        double nllM = 0.0;
        for (auto const& idx : blocks) {
            if (idx.empty()) { continue; }
            double D = 0.0, S = 0.0, sumLogP = 0.0;
            for (std::size_t c : idx) {
                double const p = std::max(memberEnergy[m][c], 0.0);
                double const e = cells[c].energy;
                double const d = cells[c].dof;
                if (!(p > 0.0) || !std::isfinite(p) || !std::isfinite(e) || e < 0.0 || !(d > 0.0)) {
                    ok = false;
                    return;
                }
                D += d;
                S += e / p;
                sumLogP += 0.5 * d * std::log(p);
            }
            if (S == 0.0) { nllM = -kInf; break; }
            nllM += 0.5 * D * std::log(2.0 * detail::kCisPi) + sumLogP
                  - std::lgamma(0.5 * D) + 0.5 * D * std::log(0.5 * S);
        }
        totalM[m] = nllM;
    }
    fitness = cisLseNegMean(totalM);
    if (std::isnan(fitness)) { ok = false; }
}

inline auto cisAbsDataNllEnsemble(std::vector<CisCell> const& cells,
                                  std::vector<std::vector<double>> const& memberEnergy,
                                  std::vector<double>& memberNll,
                                  std::vector<double>& memberWeight,
                                  bool& ok) -> double
{
    constexpr double kInf = std::numeric_limits<double>::infinity();
    ok = true;
    std::size_t const n = cells.size();
    std::size_t const M = memberEnergy.size();
    memberNll.assign(M, 0.0);
    memberWeight.assign(M, 0.0);
    if (n == 0 || M == 0) { ok = false; return kInf; }
    for (auto const& me : memberEnergy) {
        if (me.size() != n) { ok = false; return kInf; }
    }
    for (std::size_t m = 0; m < M; ++m) {
        double sum = 0.0;
        for (std::size_t c = 0; c < n; ++c) {
            double const p = std::max(memberEnergy[m][c], 0.0);
            double const e = cells[c].energy;
            if (!(p > 0.0) || !std::isfinite(p) || !std::isfinite(e) || e < 0.0 || !(cells[c].dof > 0.0)) {
                ok = false;
                sum = kInf;
                break;
            }
            sum += cisExactCellNllData(e / p, cells[c].dof, std::log(p));
        }
        memberNll[m] = sum;
    }
    if (!ok) { return kInf; }

    double const total = cisLseNegMean(memberNll);

    double mn = memberNll[0];
    for (double v : memberNll) { mn = std::min(mn, v); }
    double wsum = 0.0;
    for (std::size_t m = 0; m < M; ++m) {
        double const d = memberNll[m] - mn;
        memberWeight[m] = std::isfinite(d) ? std::exp(-d) : 0.0;
        wsum += memberWeight[m];
    }
    if (wsum > 0.0) {
        for (auto& w : memberWeight) { w /= wsum; }
    } else {
        for (auto& w : memberWeight) { w = 1.0 / static_cast<double>(M); }
    }
    if (std::isnan(total)) { ok = false; }
    return total;
}

inline void computeCisNll(std::vector<CisCell>& cells, double& shapeNll, double& scaleNll, bool& ok)
{
    ok = true;
    std::size_t const n = cells.size();
    if (n == 0) {
        shapeNll = 0.0;
        scaleNll = 0.0;
        ok = false;
        return;
    }

    std::vector<double> ratio(n);
    double sumRatio = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const pred = std::max(cells[i].predictedEnergy, detail::kCisEnergyEps);
        double const e = std::max(cells[i].energy, 0.0);
        double r = (e + detail::kCisEnergyEps) / (pred + detail::kCisEnergyEps);
        if (!std::isfinite(r) || r <= 0.0) {
            r = 1.0;
            ok = false;
        }
        ratio[i] = r;
        sumRatio += r;
    }
    double const meanRatio = sumRatio / static_cast<double>(n);
    double const meanRatioSafe = (std::isfinite(meanRatio) && meanRatio > detail::kCisEnergyEps) ? meanRatio : 1.0;

    double sumExtraVar = 0.0;
    int epiCount = 0;
    for (auto const& c : cells) {
        if (c.hasEpi) {
            sumExtraVar += cisEpiExtraVar(c);
            ++epiCount;
        }
    }
    double const avgExtraVar = (epiCount > 0) ? (sumExtraVar / epiCount) : 0.0;
    scaleNll = cisLognormalNllOfRatio(meanRatioSafe, static_cast<double>(n), avgExtraVar);

    shapeNll = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const rNorm = ratio[i] / meanRatioSafe;
        double const extraVar = cisEpiExtraVar(cells[i]);
        double const cellNll = cisLognormalNllOfRatio(rNorm, 1.0, extraVar);
        cells[i].nll = cellNll;
        shapeNll += cellNll;
        if (!std::isfinite(cellNll)) { ok = false; }
    }

    if (!std::isfinite(shapeNll) || !std::isfinite(scaleNll)) { ok = false; }
}

inline auto evaluateCisScoring(
    std::vector<double> const& residual,
    std::vector<double> const& sigma2,
    std::vector<double> const& logVar = {},
    CisScoringConfig const& cfg = {},
    CisPerCellEnergyInput const& directEnergy = {}) -> CisScoringResult
{
    CisScoringResult out;
    if (residual.size() < 2 || sigma2.empty()) {
        out.ok = false;
        return out;
    }

    auto const dec = cfg.useSavitzkyGolayBoundary ? haarForwardSgBoundary(residual) : haarForward(residual);
    out.cells = buildCisCells(dec, sigma2, logVar, cfg.scoreDcCell);
    applyDirectCellEnergyStub(out.cells, directEnergy);

    bool ok = true;
    int const dispatchMode = (cfg.exactLikelihood == 4) ? 3 : cfg.exactLikelihood;
    if (cfg.exactLikelihood == 4) {
        cisGrid48Counters().law1FallbackToLaw3.fetch_add(1, std::memory_order_relaxed);
    }
    if (dispatchMode == 3) {
        auto const wdec = cisRecordingHalfWidths(residual, cfg.residualPrecision);
        std::size_t k = 0;
        auto attach = [&](std::vector<double> const& ws) {
            for (double const w : ws) { if (k < out.cells.size()) { out.cells[k].recordDelta = w * w; } ++k; }
        };
        for (auto const& d : wdec.details) { attach(d); }
        if (cfg.scoreDcCell) { attach(wdec.approxFinal); }
    }
    bool const ensembleOn = cisEnsembleModeFromEnv() != 0;
    bool ensembleUsed = false;
    if (dispatchMode >= 1 && dispatchMode <= 3) {
        if (ensembleOn && !directEnergy.memberEnergy.empty()) {
            bool sizeOk = true;
            for (auto const& me : directEnergy.memberEnergy) {
                if (me.size() != out.cells.size()) { sizeOk = false; break; }
            }
            if (sizeOk) {
                computeCisNllEnsemble(out.cells, dispatchMode, directEnergy.memberEnergy,
                                      out.shapeNll, out.scaleNll, out.sHat, ok);
                ensembleUsed = true;
                if (cisJointCoarseFromEnv() != 0) {
                    bool jOk = true;
                    double fJ = 0.0;
                    computeCisNllEnsembleJointCoarse(out.cells, directEnergy.memberEnergy, fJ, jOk);
                    if (jOk) { out.shapeNll = fJ; out.scaleNll = 0.0; ok = true; }
                    else { ok = false; }
                }
            }
        }
        if (!ensembleUsed) {
            if (ensembleOn) { cisGrid48Counters().ensembleSkipped.fetch_add(1, std::memory_order_relaxed); }
            computeCisNllExact(out.cells, dispatchMode, out.shapeNll, out.scaleNll, out.sHat, ok);
        }
    } else {
        if (ensembleOn) { cisGrid48Counters().ensembleSkipped.fetch_add(1, std::memory_order_relaxed); }
        computeCisNll(out.cells, out.shapeNll, out.scaleNll, ok);
    }

    out.fitness = cfg.useScaleInvariantFitness ? out.shapeNll : (out.shapeNll + out.scaleNll);
    out.ok = ok && std::isfinite(out.fitness);
    return out;
}

}
