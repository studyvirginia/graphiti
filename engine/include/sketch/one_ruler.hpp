#pragma once

#include "sketch/cis_scoring.hpp"
#include "sketch/point_segment_kernel.hpp"
#include "sketch/one_residual.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Sketch {

inline auto oneRulerSiCisShapeNll(std::vector<Operon::Scalar> const& curveX,
                                  std::vector<Operon::Scalar> const& curveY,
                                  std::vector<Operon::Scalar> const& dataX,
                                  std::vector<Operon::Scalar> const& dataY,
                                  std::vector<Operon::Scalar> const& pointS2,
                                  int exactMode = -1) -> double
{
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::size_t const M = std::min(curveX.size(), curveY.size());
    std::size_t const N = std::min(dataX.size(), dataY.size());
    if (M < 2 || N < 2) { return kInf; }
    std::vector<Vec2> strokePts(N), candPts(M);
    for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX[i], dataY[i]}; }
    for (std::size_t j = 0; j < M; ++j) {
        if (!std::isfinite(curveX[j]) || !std::isfinite(curveY[j])) { return kInf; }
        candPts[j] = Vec2{curveX[j], curveY[j]};
    }
    auto const psRes = ((oneResidualModeFromEnv() & 1) != 0)
        ? evaluatePointSegmentOne(strokePts, candPts, false, false,
                                  oneResidualUnsigned(oneResidualModeFromEnv()))
        : evaluatePointSegment(strokePts, candPts);
    if (!psRes.ok) { return kInf; }
    std::vector<double> residualD(psRes.residual.size());
    for (std::size_t i = 0; i < psRes.residual.size(); ++i) { residualD[i] = static_cast<double>(psRes.residual[i]); }
    std::vector<double> sigma2D;
    if (!pointS2.empty()) {
        sigma2D.assign(pointS2.begin(), pointS2.end());
    } else {
        sigma2D.assign(1, 1.0);
    }
    CisScoringConfig cfg{};
    cfg.useScaleInvariantFitness = true;
    if (exactMode >= 0) { cfg.exactLikelihood = exactMode; }
    auto const cisRes = evaluateCisScoring(residualD, sigma2D, {}, cfg);
    if (!cisRes.ok || !std::isfinite(cisRes.fitness)) { return kInf; }
    return cisRes.fitness;
}

}
