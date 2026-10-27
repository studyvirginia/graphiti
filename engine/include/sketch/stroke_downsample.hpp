#pragma once


#include "sketch/point_segment_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Sketch {

struct DownsampleConfig {
    bool useDownsample{false};
    double relTol{std::sqrt(static_cast<double>(std::numeric_limits<Operon::Scalar>::epsilon()))};
};

struct DownsampleResult {
    std::vector<std::size_t> indices;
    std::size_t requestedN{0};
    bool usedAllPoints{false};
    bool ok{true};
};

struct SelfValidatingResult {
    std::vector<std::size_t> indices;
    std::size_t n{0};
    bool converged{false};
    bool ok{true};
};

namespace ds_detail {

inline auto cumulativeArcLength(std::vector<Vec2> const& pts) -> std::vector<Operon::Scalar>
{
    std::vector<Operon::Scalar> s(pts.size(), Operon::Scalar{0});
    for (std::size_t i = 1; i < pts.size(); ++i) {
        Operon::Scalar const dx = pts[i].x - pts[i - 1].x;
        Operon::Scalar const dy = pts[i].y - pts[i - 1].y;
        Operon::Scalar seg2 = dx * dx + dy * dy;
        if (!std::isfinite(seg2) || seg2 < 0) { seg2 = 0; }
        Operon::Scalar seg = std::sqrt(seg2);
        if (!std::isfinite(seg)) { seg = 0; }
        s[i] = s[i - 1] + seg;
    }
    return s;
}

inline auto uniformFloorIndices(std::vector<Operon::Scalar> const& arcLen, std::size_t floorCount)
    -> std::vector<std::size_t>
{
    std::size_t const numPts = arcLen.size();
    std::vector<std::size_t> out;
    if (numPts == 0 || floorCount == 0) { return out; }
    floorCount = std::min(floorCount, numPts);

    if (floorCount == 1) {
        out.push_back(numPts / 2);
        return out;
    }
    if (numPts <= floorCount) {
        out.reserve(numPts);
        for (std::size_t i = 0; i < numPts; ++i) { out.push_back(i); }
        return out;
    }

    Operon::Scalar const total = arcLen.back();
    out.reserve(floorCount);
    constexpr auto kEps = std::numeric_limits<Operon::Scalar>::epsilon();
    if (total <= kEps * Operon::Scalar{100}) {
        for (std::size_t k = 0; k < floorCount; ++k) {
            std::size_t const idx = (numPts - 1) * k / (floorCount - 1);
            out.push_back(idx);
        }
    } else {
        for (std::size_t k = 0; k < floorCount; ++k) {
            Operon::Scalar const target = total * static_cast<Operon::Scalar>(k)
                / static_cast<Operon::Scalar>(floorCount - 1);
            auto it = std::lower_bound(arcLen.begin(), arcLen.end(), target);
            std::size_t idx = static_cast<std::size_t>(std::distance(arcLen.begin(), it));
            if (idx >= numPts) { idx = numPts - 1; }
            out.push_back(idx);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

inline void farthestPointFill(std::vector<Vec2> const& pts, std::vector<std::size_t>& selected, std::size_t target)
{
    std::size_t const numPts = pts.size();
    if (selected.empty() && numPts > 0) { selected.push_back(0); }

    std::vector<char> isSelected(numPts, 0);
    for (auto idx : selected) {
        if (idx < numPts) { isSelected[idx] = 1; }
    }
    std::vector<Operon::Scalar> minDist(numPts, std::numeric_limits<Operon::Scalar>::max());
    auto updateDist = [&](std::size_t s) {
        Vec2 const& sp = pts[s];
        for (std::size_t i = 0; i < numPts; ++i) {
            if (isSelected[i]) { continue; }
            Operon::Scalar const dx = pts[i].x - sp.x;
            Operon::Scalar const dy = pts[i].y - sp.y;
            Operon::Scalar d2 = dx * dx + dy * dy;
            if (!std::isfinite(d2)) { d2 = 0; }
            if (d2 < minDist[i]) { minDist[i] = d2; }
        }
    };
    for (auto idx : selected) {
        if (idx < numPts) { updateDist(idx); }
    }

    while (selected.size() < target && selected.size() < numPts) {
        std::size_t best = numPts;
        Operon::Scalar bestD = Operon::Scalar{-1};
        for (std::size_t i = 0; i < numPts; ++i) {
            if (isSelected[i]) { continue; }
            if (minDist[i] > bestD) {
                bestD = minDist[i];
                best = i;
            }
        }
        if (best == numPts) { break; }
        isSelected[best] = 1;
        selected.push_back(best);
        updateDist(best);
    }
}

}

inline auto selectDownsampleIndices(std::vector<Vec2> const& strokePoints, std::size_t n) -> DownsampleResult
{
    DownsampleResult out;
    std::size_t const numPts = strokePoints.size();
    out.requestedN = n;
    if (numPts == 0) {
        out.ok = false;
        return out;
    }
    if (n >= numPts) {
        out.indices.resize(numPts);
        for (std::size_t i = 0; i < numPts; ++i) { out.indices[i] = i; }
        out.usedAllPoints = true;
        return out;
    }
    if (n < 2) { n = std::min<std::size_t>(2, numPts); }

    auto const arcLen = ds_detail::cumulativeArcLength(strokePoints);

    std::size_t const naturalFloor = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(numPts)))));
    std::size_t const halfBudget = n / 2 + (n % 2);
    std::size_t const floorCount = std::min({ naturalFloor, n, std::max<std::size_t>(1, halfBudget) });

    auto floorIdx = ds_detail::uniformFloorIndices(arcLen, floorCount);
    if (floorIdx.size() > n) { floorIdx.resize(n); }

    std::vector<std::size_t> selected = std::move(floorIdx);
    if (selected.empty()) { selected.push_back(0); }
    ds_detail::farthestPointFill(strokePoints, selected, n);

    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
    out.indices = std::move(selected);
    return out;
}

inline auto subsetMeanFitness(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::size_t> const& indices,
    std::vector<Vec2> const& candidatePolyline) -> Operon::Scalar
{
    std::vector<Vec2> subset;
    subset.reserve(indices.size());
    for (auto idx : indices) {
        if (idx < strokePoints.size()) { subset.push_back(strokePoints[idx]); }
    }
    if (subset.empty()) { return Operon::Scalar{0}; }
    auto const res = evaluatePointSegment(subset, candidatePolyline);
    return res.fitness / static_cast<Operon::Scalar>(subset.size());
}

inline auto downsampleSelfValidating(
    std::vector<Vec2> const& strokePoints,
    std::vector<Vec2> const& candidatePolyline,
    DownsampleConfig const& cfg = {}) -> SelfValidatingResult
{
    SelfValidatingResult out;
    std::size_t const numPts = strokePoints.size();
    if (numPts == 0) {
        out.ok = false;
        return out;
    }
    if (numPts <= 3) {
        auto const sel = selectDownsampleIndices(strokePoints, numPts);
        out.indices = sel.indices;
        out.n = out.indices.size();
        out.converged = true;
        return out;
    }

    std::size_t n = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(numPts)))));
    auto prevSel = selectDownsampleIndices(strokePoints, n);
    Operon::Scalar prevFit = subsetMeanFitness(strokePoints, prevSel.indices, candidatePolyline);
    if (!std::isfinite(prevFit)) { out.ok = false; }

    constexpr auto kEps = std::numeric_limits<Operon::Scalar>::epsilon();
    while (true) {
        if (n >= numPts) {
            out.indices = prevSel.indices;
            out.n = out.indices.size();
            out.converged = false;
            return out;
        }
        std::size_t const nextN = std::min(numPts, n * 2);
        auto sel = selectDownsampleIndices(strokePoints, nextN);
        Operon::Scalar const fit = subsetMeanFitness(strokePoints, sel.indices, candidatePolyline);
        if (!std::isfinite(fit)) { out.ok = false; }

        Operon::Scalar const denom = std::max(std::abs(prevFit), kEps * Operon::Scalar{100});
        Operon::Scalar const relChange = std::abs(fit - prevFit) / denom;

        n = nextN;
        prevSel = sel;
        prevFit = fit;

        if (relChange < static_cast<Operon::Scalar>(cfg.relTol)) {
            out.indices = sel.indices;
            out.n = sel.indices.size();
            out.converged = true;
            return out;
        }
    }
}

}
