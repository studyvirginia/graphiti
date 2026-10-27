#pragma once


#include "sketch/point_segment_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <vector>

namespace Sketch {

constexpr int kOneResidualModeDefault = 35;
inline auto oneResidualModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_ONE_RESIDUAL");
        return (e != nullptr) ? std::atoi(e) : kOneResidualModeDefault;
    }();
    return v;
}

constexpr int kOneResidualUnsignedBit = 32;
inline auto oneResidualUnsigned(int mode) -> bool
{
    return (mode & 1) != 0 && (mode & kOneResidualUnsignedBit) != 0;
}

inline auto regEvalModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_REG_EVAL");
        return (e != nullptr) ? std::atoi(e) : 0;
    }();
    return v;
}

inline auto useParallelBreedFromEnv() -> bool
{
    static bool const v = [] {
        char const* e = std::getenv("SKETCH_PARALLEL_BREED");
        return (e != nullptr) ? (std::atoi(e) != 0) : false;
    }();
    return v;
}

inline auto useGpuReverseHalfFromEnv() -> bool
{
    static bool const v = [] {
        char const* e = std::getenv("SKETCH_GPU_REVERSE_HALF");
        return (e != nullptr) ? (std::atoi(e) != 0) : true;
    }();
    return v;
}

inline auto useVarProLinearFromEnv() -> bool
{
    static bool const v = [] {
        char const* e = std::getenv("SKETCH_VARPRO_LINEAR");
        return (e != nullptr) && (std::atoi(e) != 0);
    }();
    return v;
}

template <class T>
struct PointSegOneResult {
    std::vector<T> residual;
    std::vector<T> tau;
    std::vector<T> normalX;
    std::vector<T> normalY;
    T fitness{0};
    bool ok{true};
};

template <class T, class P>
inline auto oneResidualMeanSegLen(P const* poly, std::size_t n) -> T
{
    constexpr float kEpsF = std::numeric_limits<float>::epsilon();
    T const lenSqEps = static_cast<T>(kEpsF) * T{100};
    double sum = 0.0;
    std::size_t cnt = 0;
    for (std::size_t j = 0; j + 1 < n; ++j) {
        T const dx = static_cast<T>(poly[j + 1].x) - static_cast<T>(poly[j].x);
        T const dy = static_cast<T>(poly[j + 1].y) - static_cast<T>(poly[j].y);
        T const lenSq = dx * dx + dy * dy;
        if (lenSq > lenSqEps) { sum += std::sqrt(static_cast<double>(lenSq)); ++cnt; }
    }
    return cnt > 0 ? static_cast<T>(sum / static_cast<double>(cnt)) : T{0};
}

template <class T, class P, class Q>
inline void evaluatePointSegmentOneInto(P const* strokePts, std::size_t N,
                                        Q const* cand, std::size_t Nc,
                                        bool wantTau, bool wantNormal,
                                        PointSegOneResult<T>& out, bool unsignedMag = false)
{
    out.ok = true;
    out.fitness = T{0};
    out.residual.assign(N, T{0});
    if (wantTau) { out.tau.assign(N, T{0}); } else { out.tau.clear(); }
    if (wantNormal) { out.normalX.assign(N, T{0}); out.normalY.assign(N, T{0}); }
    else { out.normalX.clear(); out.normalY.clear(); }
    if (N == 0 || Nc < 2) {
        out.ok = false;
        out.fitness = std::numeric_limits<T>::max();
        return;
    }
    constexpr float kEpsF = 1.1920929e-7f;
    T const kEps = static_cast<T>(kEpsF);
    T const lenSqEps = kEps * T{100};
    T const softMinScale = std::max(oneResidualMeanSegLen<T>(cand, Nc), kEps);
    T const kMax = std::numeric_limits<T>::max();

    T fit{0};
    for (std::size_t i = 0; i < N; ++i) {
        T const px = static_cast<T>(strokePts[i].x);
        T const py = static_cast<T>(strokePts[i].y);
        T mRun = kMax;
        T sumExp{0}, blended{0}, blendedTau{0}, blendedNx{0}, blendedNy{0};
        T blendedAbs{0};
        T hardMin = kMax;
        for (std::size_t j = 0; j + 1 < Nc; ++j) {
            T const ax = static_cast<T>(cand[j].x), ay = static_cast<T>(cand[j].y);
            T const vx = static_cast<T>(cand[j + 1].x) - ax;
            T const vy = static_cast<T>(cand[j + 1].y) - ay;
            T const lenSq = vx * vx + vy * vy;
            T t{0};
            if (lenSq > lenSqEps) {
                T const traw = ((px - ax) * vx + (py - ay) * vy) / lenSq;
                t = std::min(std::max(traw, T{0}), T{1});
            }
            T const pax = px - ax, pay = py - ay;
            T const paSq = pax * pax + pay * pay;
            T const paDotV = pax * vx + pay * vy;
            T const d2raw = paSq - T{2} * t * paDotV + t * t * lenSq;
            T const d2 = std::max(d2raw, T{0});
            T const diffX = pax - t * vx, diffY = pay - t * vy;
            hardMin = std::min(hardMin, d2);
            T const dVal = std::sqrt(std::max(d2, T{0}));
            T nX{0}, nY{1};
            if (lenSq > lenSqEps) {
                T const invLen = T{1} / std::sqrt(lenSq);
                T const tx = vx * invLen, ty = vy * invLen;
                nX = -ty; nY = tx;
            }
            T const signedVal = diffX * nX + diffY * nY;
            T const tauVal = static_cast<T>(j) + t;
            if (dVal < mRun) {
                T const rescale = std::exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                blendedNx *= rescale;
                blendedNy *= rescale;
                if (unsignedMag) { blendedAbs *= rescale; }
                mRun = dVal;
            }
            T const e = std::exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
            if (unsignedMag) {
                blendedAbs += e * std::abs(signedVal);
                if (wantNormal) {
                    T const sg = (signedVal < T{0}) ? T{-1} : T{1};
                    blendedNx += e * (sg * nX); blendedNy += e * (sg * nY);
                }
            } else if (wantNormal) { blendedNx += e * nX; blendedNy += e * nY; }
        }
        T sgnI{1};
        if (unsignedMag) {
            T const u = static_cast<T>(std::numeric_limits<T>::epsilon()) / T{2};
            T const nRound = T{2} * static_cast<T>(Nc);
            T const gammaN = (nRound * u) / (T{1} - nRound * u);
            sgnI = (std::abs(blended) > gammaN * blendedAbs && blended < T{0}) ? T{-1} : T{1};
            blended = sgnI * blendedAbs;
            if (wantNormal) { blendedNx = sgnI * blendedNx; blendedNy = sgnI * blendedNy; }
        }
        T const r = (sumExp > T{0}) ? (blended / sumExp) : T{0};
        if (!std::isfinite(r)) { out.ok = false; out.residual[i] = T{0}; }
        else { out.residual[i] = r; }
        if (wantTau) {
            T const tv = (sumExp > T{0}) ? (blendedTau / sumExp) : T{0};
            if (!std::isfinite(tv)) { out.ok = false; } else { out.tau[i] = tv; }
        }
        if (wantNormal) {
            T const nx = (sumExp > T{0}) ? (blendedNx / sumExp) : T{0};
            T const ny = (sumExp > T{0}) ? (blendedNy / sumExp) : T{0};
            if (!std::isfinite(nx) || !std::isfinite(ny)) { out.ok = false; }
            else { out.normalX[i] = nx; out.normalY[i] = ny; }
        }
        fit += hardMin;
    }
    out.fitness = fit;
    if (!std::isfinite(out.fitness)) { out.fitness = kMax; out.ok = false; }
}

struct ReverseSegHit {
    double d2{std::numeric_limits<double>::max()};
    std::size_t nearestVertex{0};
    double zx{0}, zy{0};
};

struct OneResP2d {
    double x{0}, y{0};
};

template <class P>
inline auto pointToPolylineSegOne(double qx, double qy, P const* poly, std::size_t n) -> ReverseSegHit
{
    ReverseSegHit h;
    if (n == 0) { return h; }
    if (n == 1) {
        double const dx = qx - static_cast<double>(poly[0].x), dy = qy - static_cast<double>(poly[0].y);
        h.d2 = dx * dx + dy * dy;
        h.zx = static_cast<double>(poly[0].x); h.zy = static_cast<double>(poly[0].y);
        return h;
    }
    double const lenSqEps = static_cast<double>(1.1920929e-7f) * 100.0;
    for (std::size_t j = 0; j + 1 < n; ++j) {
        double const ax = static_cast<double>(poly[j].x), ay = static_cast<double>(poly[j].y);
        double const vx = static_cast<double>(poly[j + 1].x) - ax;
        double const vy = static_cast<double>(poly[j + 1].y) - ay;
        double const lenSq = vx * vx + vy * vy;
        double t = 0.0;
        if (lenSq > lenSqEps) {
            double const traw = ((qx - ax) * vx + (qy - ay) * vy) / lenSq;
            t = std::min(std::max(traw, 0.0), 1.0);
        }
        double const pax = qx - ax, pay = qy - ay;
        double const paSq = pax * pax + pay * pay;
        double const paDotV = pax * vx + pay * vy;
        double const d2 = std::max(paSq - 2.0 * t * paDotV + t * t * lenSq, 0.0);
        if (d2 < h.d2) {
            h.d2 = d2;
            h.nearestVertex = (t < 0.5) ? j : j + 1;
            h.zx = ax + t * vx;
            h.zy = ay + t * vy;
        }
    }
    return h;
}

inline auto evaluatePointSegmentOne(std::vector<Vec2> const& strokePoints,
                                    std::vector<Vec2> const& candidatePolyline,
                                    bool computeTau, bool computeNormal,
                                    bool unsignedMag = false) -> PointSegmentResult
{
    thread_local PointSegOneResult<Operon::Scalar> r;
    evaluatePointSegmentOneInto<Operon::Scalar>(strokePoints.data(), strokePoints.size(),
                                                candidatePolyline.data(), candidatePolyline.size(),
                                                computeTau, computeNormal, r, unsignedMag);
    PointSegmentResult out;
    out.ok = r.ok;
    out.fitness = r.fitness;
    out.residual = r.residual;
    if (computeTau) { out.tau = r.tau; }
    if (computeNormal) {
        out.normalAtClosest.resize(r.normalX.size());
        for (std::size_t i = 0; i < r.normalX.size(); ++i) { out.normalAtClosest[i] = Vec2{r.normalX[i], r.normalY[i]}; }
    }
    return out;
}

}
