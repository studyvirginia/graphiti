#pragma once


#include "sketch/topology_detect.hpp"

#include <cmath>
#include <cstddef>

namespace Sketch {

struct ParametrizationGauge {
    int  n = 0;
    bool closed = false;
    bool reversed = false;
    int  shift = 0;
    bool landmarkTied = false;
};

namespace detail {

inline void paramCentroid(float const* x, float const* y, int n, double& cx, double& cy)
{
    double sx = 0.0;
    double sy = 0.0;
    for (int i = 0; i < n; ++i) {
        sx += static_cast<double>(x[i]);
        sy += static_cast<double>(y[i]);
    }
    cx = sx / static_cast<double>(n);
    cy = sy / static_cast<double>(n);
}

inline double signedArea(float const* x, float const* y, int n)
{
    double a = 0.0;
    for (int i = 0; i < n; ++i) {
        int const j = (i + 1) % n;
        a += static_cast<double>(x[i]) * static_cast<double>(y[j])
           - static_cast<double>(x[j]) * static_cast<double>(y[i]);
    }
    return 0.5 * a;
}

inline bool lexSmaller(float const* x, float const* y, int n, double cx, double cy,
                       bool revA, int shiftA, bool revB, int shiftB)
{
    for (int k = 0; k < n; ++k) {
        int const ia = revA ? (n - 1 - ((k + shiftA) % n)) : ((k + shiftA) % n);
        int const ib = revB ? (n - 1 - ((k + shiftB) % n)) : ((k + shiftB) % n);
        double const ax = static_cast<double>(x[ia]) - cx;
        double const bx = static_cast<double>(x[ib]) - cx;
        if (ax < bx) { return true; }
        if (ax > bx) { return false; }
        double const ay = static_cast<double>(y[ia]) - cy;
        double const by = static_cast<double>(y[ib]) - cy;
        if (ay < by) { return true; }
        if (ay > by) { return false; }
    }
    return false;
}

inline double turningAngleAt(float const* x, float const* y, int i)
{
    double const e1x = static_cast<double>(x[i]) - static_cast<double>(x[i - 1]);
    double const e1y = static_cast<double>(y[i]) - static_cast<double>(y[i - 1]);
    double const e2x = static_cast<double>(x[i + 1]) - static_cast<double>(x[i]);
    double const e2y = static_cast<double>(y[i + 1]) - static_cast<double>(y[i]);
    double const cross = e1x * e2y - e1y * e2x;
    double const dot = e1x * e2x + e1y * e2y;
    return std::atan2(cross, dot);
}

inline double totalTurningAngle(float const* x, float const* y, int n)
{
    double theta = 0.0;
    for (int i = 1; i <= n - 2; ++i) { theta += turningAngleAt(x, y, i); }
    return theta;
}

}

inline auto canonicalizeParametrization(float const* x, float const* y, int n, int* perm,
                                        float gapFraction = kClosedGapFraction) -> ParametrizationGauge
{
    ParametrizationGauge g;
    g.n = n;
    if (x == nullptr || y == nullptr || perm == nullptr || n < 2) {
        for (int k = 0; k < (n > 0 ? n : 0); ++k) { perm[k] = k; }
        return g;
    }

    g.closed = isClosedStroke(x, y, n, gapFraction);

    double cx = 0.0;
    double cy = 0.0;
    detail::paramCentroid(x, y, n, cx, cy);

    double sumR2 = 0.0;
    double maxR2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double const rx = static_cast<double>(x[i]) - cx;
        double const ry = static_cast<double>(y[i]) - cy;
        double const r2 = rx * rx + ry * ry;
        sumR2 += r2;
        if (r2 > maxR2) { maxR2 = r2; }
    }
    double const rmsR = std::sqrt(sumR2 / static_cast<double>(n));
    double const scaleRef = rmsR > 0.0 ? rmsR : 1.0;
    double const epsR2 = 1e-9 * (scaleRef * scaleRef) + 1e-30;

    if (g.closed) {
        double const area = detail::signedArea(x, y, n);
        double const areaEps = 1e-9 * (scaleRef * scaleRef);
        if (area < -areaEps) {
            g.reversed = true;
        } else if (area <= areaEps) {
            g.reversed = detail::lexSmaller(x, y, n, cx, cy, true, 0, false, 0);
        } else {
            g.reversed = false;
        }
    } else if (n < 3) {
        g.reversed = false;
    } else {
        double const theta = detail::totalTurningAngle(x, y, n);
        double const epsTheta = 1e-9;
        if (theta > epsTheta) {
            g.reversed = false;
        } else if (theta < -epsTheta) {
            g.reversed = true;
        } else {
            bool tied = true;
            bool reverse = false;
            for (int i = 1; i <= n - 2; ++i) {
                double const kFwd = detail::turningAngleAt(x, y, i);
                double const kRev = -detail::turningAngleAt(x, y, n - 1 - i);
                if (kFwd < kRev - epsTheta) { tied = false; reverse = false; break; }
                if (kFwd > kRev + epsTheta) { tied = false; reverse = true; break; }
            }
            g.reversed = tied ? false : reverse;
        }
    }

    g.shift = 0;
    if (g.closed) {
        auto inputIndex = [&](int k) -> int {
            return g.reversed ? (n - 1 - (k % n)) : (k % n);
        };

        double best = 0.0;
        for (int k = 0; k < n; ++k) {
            int const i = inputIndex(k);
            double const rx = static_cast<double>(x[i]) - cx;
            double const ry = static_cast<double>(y[i]) - cy;
            double const r2 = rx * rx + ry * ry;
            if (r2 > best) { best = r2; }
        }

        int chosen = -1;
        int tieCount = 0;
        for (int k = 0; k < n; ++k) {
            int const i = inputIndex(k);
            double const rx = static_cast<double>(x[i]) - cx;
            double const ry = static_cast<double>(y[i]) - cy;
            double const r2 = rx * rx + ry * ry;
            if (r2 >= best - epsR2) {
                ++tieCount;
                if (chosen < 0
                    || detail::lexSmaller(x, y, n, cx, cy,
                                          g.reversed, k, g.reversed, chosen)) {
                    chosen = k;
                }
            }
        }
        g.shift = (chosen >= 0) ? chosen : 0;
        g.landmarkTied = (tieCount > 1);
    }

    for (int k = 0; k < n; ++k) {
        int const pos = (k + g.shift) % n;
        perm[k] = g.reversed ? (n - 1 - pos) : pos;
    }
    return g;
}

inline void applyParametrizationPermutation(float const* src, int const* perm, int n, float* dst)
{
    if (src == nullptr || perm == nullptr || dst == nullptr) { return; }
    for (int k = 0; k < n; ++k) { dst[k] = src[perm[k]]; }
}

inline void undoParametrizationPermutation(float const* canonical, int const* perm, int n, float* dst)
{
    if (canonical == nullptr || perm == nullptr || dst == nullptr) { return; }
    for (int k = 0; k < n; ++k) { dst[perm[k]] = canonical[k]; }
}

}
