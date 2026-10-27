#pragma once


#include <cmath>
#include <limits>
#include <vector>

#include <operon/operon/core/types.hpp>

#include "sketch/grid_sdf.hpp"
#include "sketch/point_segment_kernel.hpp"

namespace Sketch {

struct SymmetricFitnessConfig {
    bool useSymmetricFitness{false};
};

struct SymmetricFitnessResult {
    Operon::Scalar forward{0};
    Operon::Scalar reverse{0};
    Operon::Scalar symmetric{0};
    bool usedReverse{false};
    bool ok{true};
};

inline auto evaluateSymmetricFitness(
    std::vector<Vec2> const& strokePoints,
    std::vector<Vec2> const& candidatePolyline,
    GridSdf const& strokeGrid,
    SymmetricFitnessConfig const&  = {}) -> SymmetricFitnessResult
{
    SymmetricFitnessResult out;

    PointSegmentResult const fwd = evaluatePointSegment(strokePoints, candidatePolyline);
    out.ok = fwd.ok;
    out.forward = fwd.fitness;

    if (!std::isfinite(out.forward)) {
        out.forward = std::numeric_limits<Operon::Scalar>::max();
        out.ok = false;
    }

    Operon::Scalar reverseSum{0};
    bool reverseOk = strokeGrid.ok && !candidatePolyline.empty();
    if (reverseOk) {
        for (auto const& v : candidatePolyline) {
            if (!std::isfinite(v.x) || !std::isfinite(v.y)) {
                reverseOk = false;
                out.ok = false;
                continue;
            }
            Operon::Scalar const d2 = queryGridSdf(strokeGrid, v.x, v.y);
            if (!std::isfinite(d2) || d2 < 0) {
                out.ok = false;
                continue;
            }
            reverseSum += d2;
        }
    }

    if (reverseOk && std::isfinite(reverseSum)) {
        out.reverse = reverseSum;
        out.usedReverse = true;
    } else {
        out.reverse = 0;
        out.usedReverse = false;
    }

    out.symmetric = out.forward + out.reverse;
    if (!std::isfinite(out.symmetric)) {
        out.symmetric = std::numeric_limits<Operon::Scalar>::max();
        out.ok = false;
    }

    return out;
}

}
