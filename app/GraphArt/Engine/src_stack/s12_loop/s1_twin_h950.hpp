#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace S1 {

inline double pointSegmentSoftMin(double px, double py, const double* cx, const double* cy, std::size_t n, double softScale)
{
    const double kEps = 1.1920929e-7;
    const double lenSqEps = kEps * 100.0;
    const double kMax = std::numeric_limits<double>::max();
    double mRun = kMax, sumExp = 0, blendedAbs = 0, blended = 0;
    for (std::size_t j = 0; j + 1 < n; ++j) {
        const double ax = cx[j], ay = cy[j];
        const double vx = cx[j + 1] - ax, vy = cy[j + 1] - ay;
        const double lenSq = vx * vx + vy * vy;
        double t = 0;
        if (lenSq > lenSqEps) t = std::min(std::max(((px - ax) * vx + (py - ay) * vy) / lenSq, 0.0), 1.0);
        const double pax = px - ax, pay = py - ay;
        const double d2 = std::max(pax * pax + pay * pay - 2.0 * t * (pax * vx + pay * vy) + t * t * lenSq, 0.0);
        const double dX = pax - t * vx, dY = pay - t * vy;
        double nX = 0, nY = 1;
        if (lenSq > lenSqEps) { const double inv = 1.0 / std::sqrt(lenSq); nX = -vy * inv; nY = vx * inv; }
        const double sv = dX * nX + dY * nY;
        const double dVal = std::sqrt(d2);
        if (dVal < mRun) {
            const double rs = std::exp(-(mRun - dVal) / softScale);
            sumExp *= rs; blended *= rs; blendedAbs *= rs; mRun = dVal;
        }
        const double e = std::exp(-(dVal - mRun) / softScale);
        sumExp += e; blended += e * sv; blendedAbs += e * std::abs(sv);
    }
    const double u = 5.9604645e-8, gN = (94.0 * u) / (1.0 - 94.0 * u);
    const double sgn = (std::abs(blended) > gN * blendedAbs && blended < 0.0) ? -1.0 : 1.0;
    return sumExp > 0 ? sgn * blendedAbs / sumExp : 0.0;
}

inline double meanSegLen(const double* cx, const double* cy, std::size_t n)
{
    const double lenSqEps = 1.1920929e-7 * 100.0;
    double sum = 0; std::size_t cnt = 0;
    for (std::size_t j = 0; j + 1 < n; ++j) {
        const double dx = cx[j + 1] - cx[j], dy = cy[j + 1] - cy[j];
        const double l2 = dx * dx + dy * dy;
        if (l2 > lenSqEps) { sum += std::sqrt(l2); ++cnt; }
    }
    return cnt > 0 ? sum / (double)cnt : 0.0;
}

inline void s1Residual(const double* sx, const double* sy, const double* fx, const double* fy, std::size_t N, double* out)
{
    const double sc = std::max(meanSegLen(fx, fy, N), 1.1920929e-7);
    for (std::size_t i = 0; i < N; ++i) out[i] = pointSegmentSoftMin(sx[i], sy[i], fx, fy, N, sc);
}

}
