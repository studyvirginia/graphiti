#pragma once


#include <cmath>
#include <cstddef>

namespace Sketch {

inline auto strokeClosedness(float const* x, float const* y, int n) -> float
{
    if (x == nullptr || y == nullptr || n < 2) { return 0.0f; }

    double total = 0.0;
    for (int i = 1; i < n; ++i) {
        double const dx = static_cast<double>(x[i]) - static_cast<double>(x[i - 1]);
        double const dy = static_cast<double>(y[i]) - static_cast<double>(y[i - 1]);
        double const seg = std::sqrt(dx * dx + dy * dy);
        if (std::isfinite(seg)) { total += seg; }
    }

    double const gx = static_cast<double>(x[n - 1]) - static_cast<double>(x[0]);
    double const gy = static_cast<double>(y[n - 1]) - static_cast<double>(y[0]);
    double const gap = std::sqrt(gx * gx + gy * gy);

    constexpr double kEps = 1e-12;
    if (!(total > kEps) || !std::isfinite(gap)) { return 0.0f; }

    double closedness = 1.0 - gap / total;
    if (closedness < 0.0) { closedness = 0.0; }
    if (closedness > 1.0) { closedness = 1.0; }
    return static_cast<float>(closedness);
}

inline constexpr float kClosedGapFraction = 0.05f;

inline auto isClosedStroke(float const* x, float const* y, int n, float gapFraction = kClosedGapFraction) -> bool
{
    return strokeClosedness(x, y, n) > (1.0f - gapFraction);
}

}
