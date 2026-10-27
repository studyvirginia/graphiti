#ifndef SKETCH_SEGMENT_HPP
#define SKETCH_SEGMENT_HPP

#ifndef SKETCH_SEGMENT_HPP_NO_INCLUDES
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "pair_individual.hpp"
#endif

namespace Sketch {

struct Corner {
    std::size_t index{0};
    double score{0};
    double angle{0};
    int scales{0};
    int finestK{0};
};

inline constexpr int kSegmentLadder[3] = {8, 12, 16};
inline constexpr int kSegmentBootstrap = 64;
inline constexpr double kSegmentLevel = 0.95;
inline constexpr int kSegmentResample = 240;
inline constexpr unsigned kSegmentSeed = 20260908u;

namespace segdetail {

inline auto unwrappedTangentAngle(std::vector<Operon::Scalar> const& x,
                                  std::vector<Operon::Scalar> const& y,
                                  bool closed, double* loopTotal = nullptr) -> std::vector<double>
{
    constexpr double TWO_PI = 6.283185307179586;
    constexpr double PI = 3.141592653589793;
    std::size_t const M = std::min(x.size(), y.size());
    std::vector<double> th;
    if (loopTotal != nullptr) { *loopTotal = 0.0; }
    if (M < 3) { return th; }
    std::size_t const nseg = closed ? M : M - 1;
    th.assign(nseg, 0.0);
    double acc = 0.0, prev = 0.0, first = 0.0;
    bool have = false;
    for (std::size_t i = 0; i < nseg; ++i) {
        std::size_t const j = (i + 1) % M;
        double const dx = static_cast<double>(x[j]) - static_cast<double>(x[i]);
        double const dy = static_cast<double>(y[j]) - static_cast<double>(y[i]);
        if (dx == 0.0 && dy == 0.0) { th[i] = acc; continue; }
        double const ang = std::atan2(dy, dx);
        if (!have) { acc = ang; prev = ang; first = ang; have = true; }
        else {
            double d = ang - prev;
            while (d > PI) { d -= TWO_PI; }
            while (d < -PI) { d += TWO_PI; }
            acc += d;
            prev = ang;
        }
        th[i] = acc;
    }
    if (closed && have && loopTotal != nullptr) {
        double d = first - prev;
        while (d > PI) { d -= TWO_PI; }
        while (d < -PI) { d += TWO_PI; }
        *loopTotal = (acc + d) - th.front();
    }
    return th;
}

inline auto totalArcLength(std::vector<Operon::Scalar> const& x,
                           std::vector<Operon::Scalar> const& y, bool closed) -> double
{
    std::size_t const M = std::min(x.size(), y.size());
    if (M < 2) { return 0.0; }
    double L = 0.0;
    for (std::size_t i = 1; i < M; ++i) {
        L += std::hypot(static_cast<double>(x[i]) - static_cast<double>(x[i - 1]),
                        static_cast<double>(y[i]) - static_cast<double>(y[i - 1]));
    }
    if (closed) {
        L += std::hypot(static_cast<double>(x[0]) - static_cast<double>(x[M - 1]),
                        static_cast<double>(y[0]) - static_cast<double>(y[M - 1]));
    }
    return L;
}

inline auto cumulativeArcLength(std::vector<Operon::Scalar> const& x,
                                std::vector<Operon::Scalar> const& y) -> std::vector<double>
{
    std::size_t const M = std::min(x.size(), y.size());
    std::vector<double> d(M, 0.0);
    for (std::size_t i = 1; i < M; ++i) {
        d[i] = d[i - 1] + std::hypot(static_cast<double>(x[i]) - static_cast<double>(x[i - 1]),
                                     static_cast<double>(y[i]) - static_cast<double>(y[i - 1]));
    }
    return d;
}

inline auto windowFor(int K, int n) -> int
{
    return std::max(2, static_cast<int>(std::lround(static_cast<double>(n) / (4.0 * K))));
}

inline auto turningProfile(std::vector<Operon::Scalar> const& sx,
                           std::vector<Operon::Scalar> const& sy,
                           bool closed, int w) -> std::vector<double>
{
    double loop = 0.0;
    auto const th = unwrappedTangentAngle(sx, sy, closed, &loop);
    auto const T = static_cast<long>(th.size());
    std::vector<double> S(static_cast<std::size_t>(std::max<long>(T, 0)),
                          -std::numeric_limits<double>::infinity());
    if (T < 8L * w + 1) { return S; }

    auto at = [&](long j) -> double {
        if (!closed) {
            long const k = std::max<long>(0, std::min<long>(T - 1, j));
            return th[static_cast<std::size_t>(k)];
        }
        long q = j / T;
        long r = j % T;
        if (r < 0) { r += T; q -= 1; }
        return th[static_cast<std::size_t>(r)] + static_cast<double>(q) * loop;
    };

    long const lo = closed ? 0 : 3 * w;
    long const hi = closed ? T : T - 3 * w;
    for (long j = lo; j < hi; ++j) {
        double const centre = std::abs(at(j + w) - at(j - w));
        double const right = std::abs(at(j + 3 * w) - at(j + w));
        double const left = std::abs(at(j - w) - at(j - 3 * w));
        S[static_cast<std::size_t>(j)] = centre - 0.5 * (right + left);
    }
    return S;
}

}

inline auto detectCorners(std::vector<Operon::Scalar> const& x,
                          std::vector<Operon::Scalar> const& y,
                          bool closed, double sigma, int nullMode = 0,
                          int smoother = 0) -> std::vector<Corner>
{
    std::vector<Corner> out;
    std::size_t const N = std::min(x.size(), y.size());
    if (N < 16 || !(sigma > 0.0) || !std::isfinite(sigma)) { return out; }
    for (std::size_t i = 0; i < N; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) { return out; }
    }

    std::vector<Operon::Scalar> gx(x.begin(), x.begin() + static_cast<long>(N));
    std::vector<Operon::Scalar> gy(y.begin(), y.begin() + static_cast<long>(N));
    normalizeInPlace(gx, gy);

    double const L = segdetail::totalArcLength(gx, gy, closed);
    if (!(L > 0.0) || !std::isfinite(L)) { return out; }

    int const n = kSegmentResample;

    std::vector<Operon::Scalar> nx(static_cast<std::size_t>(N)), ny(static_cast<std::size_t>(N));
    if (nullMode == 0) {
        double const Tturn = turningFourierSigned(gx, gy, closed, kSegmentLadder[0], n, smoother);
        double const kappa = 6.283185307179586 * Tturn / L;
        double const ds = L / static_cast<double>(N);
        double px = 0.0, py = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            nx[i] = static_cast<Operon::Scalar>(px);
            ny[i] = static_cast<Operon::Scalar>(py);
            double const th = kappa * (static_cast<double>(i) + 0.5) * ds;
            px += ds * std::cos(th);
            py += ds * std::sin(th);
        }
        normalizeInPlace(nx, ny);
    }
    bool const haveNull = true;

    std::mt19937 rng(kSegmentSeed);
    std::normal_distribution<double> gauss(0.0, sigma);

    struct ScaleHit { double pos; double s; double ratio; };
    std::vector<std::vector<ScaleHit>> perScale;
    std::vector<int> ladderK;

    auto statistic = [&](std::vector<Operon::Scalar> const& px,
                         std::vector<Operon::Scalar> const& py,
                         int K, int Kc, int w) -> std::vector<double> {
        std::vector<Operon::Scalar> sx, sy;
        if (!lowPassSelect(px, py, closed, K, n, smoother, sx, sy)) { return {}; }
        auto S = segdetail::turningProfile(sx, sy, closed, w);
        if (nullMode == 0) { return S; }
        std::vector<Operon::Scalar> cx, cy;
        if (!lowPassSelect(px, py, closed, Kc, n, smoother, cx, cy)) { return {}; }
        auto const C = segdetail::turningProfile(cx, cy, closed, w);
        if (C.size() != S.size()) { return {}; }
        for (std::size_t i = 0; i < S.size(); ++i) {
            if (std::isfinite(S[i]) && std::isfinite(C[i])) { S[i] -= C[i]; }
            else { S[i] = -std::numeric_limits<double>::infinity(); }
        }
        return S;
    };

    for (int K : kSegmentLadder) {
        int const w = segdetail::windowFor(K, n);
        int const Kc = std::max(2, K / 2);

        if (nullMode != 0) {
            std::vector<Operon::Scalar> cx, cy;
            if (!lowPassSelect(gx, gy, closed, Kc, n, smoother, cx, cy)) { continue; }
            nx.assign(cx.begin(), cx.end());
            ny.assign(cy.begin(), cy.end());
        }

        auto const S = statistic(gx, gy, K, Kc, w);
        auto const T = static_cast<long>(S.size());
        if (T < 8L * w + 1) { continue; }

        double tau = std::numeric_limits<double>::infinity();
        std::vector<double> tauPos;
        {
            std::vector<double> maxima;
            maxima.reserve(kSegmentBootstrap);
            std::vector<double> sum(S.size(), 0.0), sumsq(S.size(), 0.0);
            std::vector<int> cnt(S.size(), 0);
            std::vector<Operon::Scalar> rx(nx.size()), ry(ny.size());
            for (int b = 0; b < kSegmentBootstrap; ++b) {
                for (std::size_t i = 0; i < nx.size(); ++i) {
                    rx[i] = static_cast<Operon::Scalar>(static_cast<double>(nx[i]) + gauss(rng));
                    ry[i] = static_cast<Operon::Scalar>(static_cast<double>(ny[i]) + gauss(rng));
                }
                auto const Sb = statistic(rx, ry, K, Kc, w);
                double m = -std::numeric_limits<double>::infinity();
                for (double v : Sb) { if (std::isfinite(v) && v > m) { m = v; } }
                if (std::isfinite(m)) { maxima.push_back(m); }
                if (Sb.size() == S.size()) {
                    for (std::size_t i = 0; i < Sb.size(); ++i) {
                        if (!std::isfinite(Sb[i])) { continue; }
                        sum[i] += Sb[i]; sumsq[i] += Sb[i] * Sb[i]; cnt[i] += 1;
                    }
                }
            }
            if (maxima.size() >= 8) {
                std::sort(maxima.begin(), maxima.end());
                auto q = static_cast<std::size_t>(kSegmentLevel * (static_cast<double>(maxima.size()) - 1.0));
                tau = maxima[std::min(q, maxima.size() - 1)];
            }
            if (nullMode == 2) {
                double const alpha = (1.0 - kSegmentLevel) / (2.0 * static_cast<double>(K));
                double zlo = 0.0, zhi = 10.0;
                for (int it = 0; it < 80; ++it) {
                    double const zm = 0.5 * (zlo + zhi);
                    double const tail = 0.5 * std::erfc(zm / std::sqrt(2.0));
                    if (tail > alpha) { zlo = zm; } else { zhi = zm; }
                }
                double const z = 0.5 * (zlo + zhi);
                tauPos.assign(S.size(), std::numeric_limits<double>::infinity());
                for (std::size_t i = 0; i < S.size(); ++i) {
                    if (cnt[i] < 8) { continue; }
                    double const mu = sum[i] / cnt[i];
                    double const var = std::max(0.0, sumsq[i] / cnt[i] - mu * mu);
                    tauPos[i] = mu + z * std::sqrt(var);
                }
                tau = 1.0;
            }
        }
        if (!std::isfinite(tau) || !(tau > 0.0)) { continue; }

        std::vector<ScaleHit> hits;
        for (long j = 0; j < T; ++j) {
            double const v = S[static_cast<std::size_t>(j)];
            double const thr = tauPos.empty() ? tau : tauPos[static_cast<std::size_t>(j)];
            if (!std::isfinite(v) || !std::isfinite(thr) || !(thr > 0.0) || !(v > thr)) { continue; }
            bool best = true;
            for (long d = -w; d <= w && best; ++d) {
                if (d == 0) { continue; }
                long k = j + d;
                if (closed) { k = ((k % T) + T) % T; }
                else if (k < 0 || k >= T) { continue; }
                double const u = S[static_cast<std::size_t>(k)];
                if (!std::isfinite(u)) { continue; }
                if (u > v || (u == v && k < j)) { best = false; }
            }
            if (!best) { continue; }
            hits.push_back(ScaleHit{static_cast<double>(j) / static_cast<double>(T), v, v / thr});
        }
        perScale.push_back(std::move(hits));
        ladderK.push_back(K);
    }

    if (perScale.size() < 2) { return out; }

    struct Merged { double pos; double ratio; double s; int scales; int finestK; };
    std::vector<Merged> merged;
    for (std::size_t a = 0; a + 1 < perScale.size(); ++a) {
        double const tol = static_cast<double>(segdetail::windowFor(ladderK[a], n))
                         / static_cast<double>(n);
        for (auto const& hf : perScale[a + 1]) {
            for (auto const& hc : perScale[a]) {
                double dp = std::abs(hf.pos - hc.pos);
                if (closed) { dp = std::min(dp, 1.0 - dp); }
                if (dp > tol) { continue; }
                Merged m{hf.pos, std::min(hf.ratio, hc.ratio), hf.s, 2, ladderK[a + 1]};
                bool absorbed = false;
                for (auto& e : merged) {
                    double d2 = std::abs(e.pos - m.pos);
                    if (closed) { d2 = std::min(d2, 1.0 - d2); }
                    if (d2 <= tol) {
                        absorbed = true;
                        e.scales += 1;
                        if (m.finestK > e.finestK) { e.pos = m.pos; e.s = m.s; e.finestK = m.finestK; }
                        e.ratio = std::max(e.ratio, m.ratio);
                        break;
                    }
                }
                if (!absorbed) { merged.push_back(m); }
                break;
            }
        }
    }

    auto const cum = segdetail::cumulativeArcLength(gx, gy);
    for (auto const& m : merged) {
        double const target = m.pos * L;
        std::size_t idx = 0;
        double bestd = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < cum.size(); ++i) {
            double const d = std::abs(cum[i] - target);
            if (d < bestd) { bestd = d; idx = i; }
        }
        out.push_back(Corner{idx, m.ratio, m.s, m.scales, m.finestK});
    }

    std::sort(out.begin(), out.end(),
              [](Corner const& a, Corner const& b) { return a.index < b.index; });
    std::vector<Corner> dedup;
    std::size_t const minSep = std::max<std::size_t>(2, N / 60);
    for (auto const& c : out) {
        if (!dedup.empty() && c.index - dedup.back().index < minSep) {
            if (c.score > dedup.back().score) { dedup.back() = c; }
            continue;
        }
        dedup.push_back(c);
    }
    return dedup;
}

}
#endif
