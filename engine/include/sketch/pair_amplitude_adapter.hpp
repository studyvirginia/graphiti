#pragma once

#include "bracket_newton_runtime.hpp"

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>
#include <optional>

namespace Sketch {

enum class PairAmplitudeAxis { x, y };

struct PairAmplitudeNodeSlot {
    PairAmplitudeAxis axis{PairAmplitudeAxis::x};
    std::size_t nodeIndex{0};
    BracketNewtonParameterKind kind{BracketNewtonParameterKind::other};
};

template <typename Pair, typename EvaluateFn>
auto buildPerturbedAmplitudeContext(Pair const& candidate,
                                    std::vector<PairAmplitudeNodeSlot> const& slots,
                                    std::vector<double> sampleT, std::vector<double> targetX,
                                    std::vector<double> targetY, EvaluateFn&& evaluate,
                                    double delta = 1.0e-5)
    -> std::optional<BracketNewtonRuntimeContext> {
    if (slots.empty() || sampleT.empty() || sampleT.size() != targetX.size() ||
        sampleT.size() != targetY.size() || !std::isfinite(delta) || delta == 0.0) {
        return std::nullopt;
    }
    Pair baseline = candidate;
    auto baseValues = evaluate(baseline);
    if (baseValues.size() != sampleT.size()) { return std::nullopt; }
    for (double value : baseValues) {
        if (!std::isfinite(value)) { return std::nullopt; }
    }
    BracketNewtonRuntimeContext context{std::move(sampleT), std::move(targetX), std::move(targetY),
                                        0, {}, {}};
    context.candidateX = baseValues;
    context.linearBasisX.assign(context.sampleT.size() * slots.size(), 0.0);
    for (std::size_t j = 0; j < slots.size(); ++j) {
        if (slots[j].kind != BracketNewtonParameterKind::amplitude) { return std::nullopt; }
        Pair trial = candidate;
        auto& nodes = slots[j].axis == PairAmplitudeAxis::x ? trial.x.Nodes() : trial.y.Nodes();
        if (slots[j].nodeIndex >= nodes.size() || !nodes[slots[j].nodeIndex].IsConstant()) {
            return std::nullopt;
        }
        context.nodeSlotCount = std::max(context.nodeSlotCount, slots[j].nodeIndex + 1U);
        nodes[slots[j].nodeIndex].Value += delta;
        auto perturbed = evaluate(trial);
        if (perturbed.size() != context.sampleT.size()) { return std::nullopt; }
        for (std::size_t i = 0; i < perturbed.size(); ++i) {
            if (!std::isfinite(perturbed[i])) { return std::nullopt; }
            context.linearBasisX[i * slots.size() + j] = (perturbed[i] - baseValues[i]) / delta;
        }
        context.parameterMap.nodeSlots.push_back(slots[j].nodeIndex);
        context.parameterMap.kinds.push_back(slots[j].kind);
    }
    return context;
}

template <typename Pair, typename ScoreFn>
auto applyMappedAmplitudeNodes(Pair& candidate, std::vector<PairAmplitudeNodeSlot> const& slots,
                               std::vector<double> const& proposed, ScoreFn&& score) -> bool {
    if (slots.empty() || slots.size() != proposed.size()) { return false; }
    for (double value : proposed) {
        if (!std::isfinite(value) || std::abs(value) > 1.0e6) { return false; }
    }
    double const before = score(candidate);
    if (!std::isfinite(before)) { return false; }
    Pair trial = candidate;
    for (std::size_t i = 0; i < slots.size(); ++i) {
        auto const& slot = slots[i];
        if (slot.kind != BracketNewtonParameterKind::amplitude) { return false; }
        auto& nodes = slot.axis == PairAmplitudeAxis::x ? trial.x.Nodes() : trial.y.Nodes();
        if (slot.nodeIndex >= nodes.size() || !nodes[slot.nodeIndex].IsConstant()) { return false; }
        nodes[slot.nodeIndex].Value = proposed[i];
    }
    double const after = score(trial);
    if (!std::isfinite(after) || after > before) { return false; }
    candidate = std::move(trial);
    return true;
}

}
