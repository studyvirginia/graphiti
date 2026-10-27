#pragma once


#include <cmath>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>
#include "fp16_gram_runtime.hpp"
#include "uma_pingpong_runtime.hpp"

namespace Sketch {

enum class WarmStartStatus {
    NotRequested,
    StateAccepted,
    BlockedMissingCorrespondenceState,
};

struct BracketNewtonConfig {
    std::size_t population{0};
    std::size_t samples{0};
    std::size_t parameters{0};
    std::size_t sweeps{1};
    double ridgeFraction{1e-3};
    bool useFp16Storage{false};
    bool warmStart{false};
    std::size_t warmSegmentRadius{2};
    bool useNoiseWeightedNewton{false};
    bool useGeneralPathCorrespondence{false};
};

struct BracketNewtonFootpointState {
    std::vector<std::size_t> segmentIndex;
    std::vector<double> segmentT;
    std::size_t completedSweeps{0};
};

struct BracketNewtonPoint { double x{0.0}; double y{0.0}; };

inline auto recomputeFootpointsCpu(std::vector<BracketNewtonPoint> const& target,
                                   std::vector<BracketNewtonPoint> const& curve,
                                   std::size_t sweeps,
                                   std::size_t completedBefore = 0)
    -> std::optional<BracketNewtonFootpointState> {
    if (target.empty() || curve.size() < 2U || completedBefore > sweeps) return std::nullopt;
    BracketNewtonFootpointState state;
    state.segmentIndex.resize(target.size());
    state.segmentT.resize(target.size());
    for (std::size_t i = 0; i < target.size(); ++i) {
        if (!std::isfinite(target[i].x) || !std::isfinite(target[i].y)) return std::nullopt;
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j + 1U < curve.size(); ++j) {
            auto const ax = curve[j].x, ay = curve[j].y;
            auto const vx = curve[j + 1U].x - ax, vy = curve[j + 1U].y - ay;
            auto const vv = vx * vx + vy * vy;
            double t = vv > 0.0 ? ((target[i].x - ax) * vx + (target[i].y - ay) * vy) / vv : 0.0;
            t = std::clamp(t, 0.0, 1.0);
            auto const dx = target[i].x - (ax + t * vx);
            auto const dy = target[i].y - (ay + t * vy);
            auto const d2 = dx * dx + dy * dy;
            if (std::isfinite(d2) && d2 < best) { best = d2; state.segmentIndex[i] = j; state.segmentT[i] = t; }
        }
        if (!std::isfinite(best)) return std::nullopt;
    }
    state.completedSweeps = std::min(sweeps, completedBefore + 1U);
    return state;
}

struct WarmFootpointRecompute {
    std::optional<BracketNewtonFootpointState> state;
    std::size_t segmentEvaluations{0};
    bool usedWarmStart{false};
};

inline auto recomputeFootpointsCpuWarm(std::vector<BracketNewtonPoint> const& target,
                                       std::vector<BracketNewtonPoint> const& curve,
                                       std::size_t sweeps,
                                       BracketNewtonFootpointState const* prior,
                                       std::size_t radius = 2U) -> WarmFootpointRecompute {
    WarmFootpointRecompute out;
    if (target.empty() || curve.size() < 2U) return out;
    auto const segments = curve.size() - 1U;
    bool const warm = prior != nullptr && prior->segmentIndex.size() == target.size() &&
        prior->segmentT.size() == target.size() && prior->completedSweeps <= sweeps &&
        std::all_of(prior->segmentIndex.begin(), prior->segmentIndex.end(),
                    [segments](std::size_t index) { return index < segments; });
    if (!warm) {
        out.state = recomputeFootpointsCpu(target, curve, sweeps);
        out.segmentEvaluations = target.size() * segments;
        return out;
    }
    BracketNewtonFootpointState state;
    state.segmentIndex.resize(target.size()); state.segmentT.resize(target.size());
    for (std::size_t i = 0; i < target.size(); ++i) {
        auto const center = prior->segmentIndex[i];
        auto const begin = center > radius ? center - radius : 0U;
        auto const end = std::min(segments, center + radius + 1U);
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = begin; j < end; ++j) {
            ++out.segmentEvaluations;
            auto const ax = curve[j].x, ay = curve[j].y;
            auto const vx = curve[j + 1U].x - ax, vy = curve[j + 1U].y - ay;
            auto const vv = vx * vx + vy * vy;
            double t = vv > 0.0 ? ((target[i].x - ax) * vx + (target[i].y - ay) * vy) / vv : 0.0;
            t = std::clamp(t, 0.0, 1.0);
            auto const dx = target[i].x - (ax + t * vx), dy = target[i].y - (ay + t * vy);
            auto const d2 = dx * dx + dy * dy;
            if (std::isfinite(d2) && d2 < best) {
                best = d2; state.segmentIndex[i] = j; state.segmentT[i] = t;
            }
        }
        if (!std::isfinite(best)) return {};
    }
    state.completedSweeps = std::min(sweeps, prior->completedSweeps + 1U);
    out.state = std::move(state); out.usedWarmStart = true;
    return out;
}

inline auto validFootpointState(BracketNewtonFootpointState const& state,
                                std::size_t samples, std::size_t sweeps) -> bool {
    if (state.segmentIndex.size() != samples || state.segmentT.size() != samples ||
        state.completedSweeps > sweeps) return false;
    for (double t : state.segmentT) if (!std::isfinite(t) || t < 0.0 || t > 1.0) return false;
    return true;
}

inline auto generalPathCorrespondence(BracketNewtonConfig const& cfg,
                                      std::vector<BracketNewtonPoint> const& target,
                                      std::vector<BracketNewtonPoint> const& curve,
                                      std::size_t completedBefore = 0)
    -> std::optional<BracketNewtonFootpointState> {
    return recomputeFootpointsCpu(target, curve, cfg.sweeps, completedBefore);
}

struct BracketNewtonGeneration {
    std::vector<double> residual;
    std::vector<double> gram;
    std::vector<double> rhs;
    std::vector<double> update;
    std::vector<double> parentPacingState;
    std::vector<double> fp16Basis;
    std::vector<double> fp16Residual;
    std::optional<std::size_t> parentIndex;
    BracketNewtonFootpointState footpoints;
    std::vector<BracketNewtonPoint> targetPoints;
    std::vector<BracketNewtonPoint> candidateCurveSamples;
};

struct BracketNewtonParentPacingMap {
    std::vector<std::size_t> parentByOffspring;
    std::vector<std::vector<double>> pacingStateByParent;
    std::vector<BracketNewtonFootpointState> footpointStateByParent;
};

inline auto validParentPacingMap(BracketNewtonParentPacingMap const& map,
                                 std::size_t offspringCount, std::size_t parentCount) -> bool {
    if (map.parentByOffspring.size() != offspringCount ||
        map.pacingStateByParent.size() != parentCount) return false;
    if (!map.footpointStateByParent.empty() && map.footpointStateByParent.size() != parentCount)
        return false;
    for (auto parent : map.parentByOffspring) if (parent >= parentCount) return false;
    for (auto const& state : map.pacingStateByParent)
        for (double value : state) if (!std::isfinite(value)) return false;
    return true;
}

struct BracketNewtonResult {
    bool ok{false};
    bool usedFallback{false};
    double ridge{0.0};
    std::vector<double> update;
    std::vector<double> nextPacingState;
    WarmStartStatus warmStartStatus{WarmStartStatus::NotRequested};
    BracketNewtonFootpointState nextFootpoints;
    std::size_t completedSweeps{0};
    std::size_t correspondenceEvaluations{0};
};

enum class BracketNewtonParameterKind { amplitude, frequency, phase, other };

struct BracketNewtonParameterMap {
    std::vector<std::size_t> nodeSlots;
    std::vector<BracketNewtonParameterKind> kinds;
};

struct BracketNewtonRuntimeContext {
    std::vector<double> sampleT;
    std::vector<double> targetX;
    std::vector<double> targetY;
    std::size_t nodeSlotCount{0};
    BracketNewtonParameterMap parameterMap;
    std::vector<double> linearBasisX;
    std::vector<double> candidateX;
    std::vector<double> candidateY;
};

using BracketNewtonJacobianBuilder = std::function<bool(
    BracketNewtonRuntimeContext const&, BracketNewtonGeneration&)>;

inline auto validParameterMap(BracketNewtonParameterMap const& map, std::size_t parameters,
                             std::size_t nodeSlots) -> bool {
    if (parameters == 0U || map.nodeSlots.size() != parameters || map.kinds.size() != parameters ||
        nodeSlots == 0U) { return false; }
    std::vector<bool> seen(nodeSlots, false);
    for (std::size_t i = 0; i < parameters; ++i) {
        if (map.nodeSlots[i] >= nodeSlots || seen[map.nodeSlots[i]]) { return false; }
        seen[map.nodeSlots[i]] = true;
    }
    return true;
}

inline auto validRuntimeContext(BracketNewtonRuntimeContext const& context,
                                std::size_t parameters) -> bool {
    if (context.sampleT.empty() || context.sampleT.size() != context.targetX.size() ||
        context.sampleT.size() != context.targetY.size()) { return false; }
    for (double value : context.sampleT) {
        if (!std::isfinite(value)) { return false; }
    }
    for (double value : context.targetX) {
        if (!std::isfinite(value)) { return false; }
    }
    for (double value : context.targetY) {
        if (!std::isfinite(value)) { return false; }
    }
    if ((!context.candidateX.empty() || !context.candidateY.empty()) &&
        (context.candidateX.size() != context.sampleT.size() ||
         context.candidateY.size() != context.sampleT.size())) return false;
    return validParameterMap(context.parameterMap, parameters, context.nodeSlotCount);
}

inline auto newtonPointWeights(BracketNewtonConfig const& cfg, std::size_t samples,
                               std::vector<double> const* pointSigma2) -> std::vector<double> {
    std::vector<double> weights(samples, 1.0);
    if (!cfg.useNoiseWeightedNewton || pointSigma2 == nullptr || pointSigma2->size() != samples) {
        return weights;
    }
    for (std::size_t i = 0; i < samples; ++i) {
        double const s2 = (*pointSigma2)[i];
        if (!std::isfinite(s2) || s2 <= 0.0) { return std::vector<double>(samples, 1.0); }
        weights[i] = 1.0 / s2;
    }
    return weights;
}

inline auto buildLinearAmplitudeGeneration(BracketNewtonRuntimeContext const& context,
                                            BracketNewtonConfig const& cfg,
                                            std::vector<double> const* pointSigma2 = nullptr)
    -> std::optional<BracketNewtonGeneration> {
    if (cfg.population != 1U || cfg.samples != context.sampleT.size() ||
        !validRuntimeContext(context, cfg.parameters) ||
        context.linearBasisX.size() != cfg.samples * cfg.parameters) {
        return std::nullopt;
    }
    for (auto kind : context.parameterMap.kinds) {
        if (kind != BracketNewtonParameterKind::amplitude) { return std::nullopt; }
    }
    BracketNewtonGeneration generation;
    generation.residual = context.targetX;
    generation.gram.assign(cfg.parameters * cfg.parameters, 0.0);
    generation.rhs.assign(cfg.parameters, 0.0);
    auto const weights = newtonPointWeights(cfg, cfg.samples, pointSigma2);
    for (std::size_t i = 0; i < cfg.samples; ++i) {
        double const wi = weights[i];
        for (std::size_t a = 0; a < cfg.parameters; ++a) {
            double const ba = context.linearBasisX[i * cfg.parameters + a];
            if (!std::isfinite(ba)) { return std::nullopt; }
            generation.rhs[a] += ba * context.targetX[i] * wi;
            for (std::size_t b = 0; b < cfg.parameters; ++b) {
                double const bb = context.linearBasisX[i * cfg.parameters + b];
                generation.gram[a * cfg.parameters + b] += ba * bb * wi;
            }
        }
    }
    if (context.candidateX.size() == cfg.samples && context.candidateY.size() == cfg.samples) {
        generation.targetPoints.resize(cfg.samples);
        generation.candidateCurveSamples.resize(cfg.samples);
        for (std::size_t i = 0; i < cfg.samples; ++i) {
            generation.targetPoints[i] = {context.targetX[i], context.targetY[i]};
            generation.candidateCurveSamples[i] = {context.candidateX[i], context.candidateY[i]};
        }
    }
    return generation;
}

inline auto buildMappedGaussNewtonGeneration(BracketNewtonRuntimeContext const& context,
                                             BracketNewtonConfig const& cfg,
                                             std::vector<double> const& jacobianX,
                                             std::vector<double> const* pointSigma2 = nullptr)
    -> std::optional<BracketNewtonGeneration> {
    if (cfg.population != 1U || cfg.samples != context.sampleT.size() ||
        !validRuntimeContext(context, cfg.parameters) ||
        context.candidateX.size() != cfg.samples || context.candidateY.size() != cfg.samples ||
        jacobianX.size() != cfg.samples * cfg.parameters) return std::nullopt;
    BracketNewtonGeneration generation;
    generation.residual.resize(cfg.samples);
    generation.gram.assign(cfg.parameters * cfg.parameters, 0.0);
    generation.rhs.assign(cfg.parameters, 0.0);
    generation.targetPoints.resize(cfg.samples);
    generation.candidateCurveSamples.resize(cfg.samples);
    auto const weights = newtonPointWeights(cfg, cfg.samples, pointSigma2);
    for (std::size_t i = 0; i < cfg.samples; ++i) {
        double const wi = weights[i];
        auto const residual = context.targetX[i] - context.candidateX[i];
        if (!std::isfinite(residual)) return std::nullopt;
        generation.residual[i] = residual;
        generation.targetPoints[i] = {context.targetX[i], context.targetY[i]};
        generation.candidateCurveSamples[i] = {context.candidateX[i], context.candidateY[i]};
        for (std::size_t a = 0; a < cfg.parameters; ++a) {
            auto const ja = jacobianX[i * cfg.parameters + a];
            if (!std::isfinite(ja)) return std::nullopt;
            generation.rhs[a] += ja * residual * wi;
            for (std::size_t b = 0; b < cfg.parameters; ++b) {
                auto const jb = jacobianX[i * cfg.parameters + b];
                if (!std::isfinite(jb)) return std::nullopt;
                generation.gram[a * cfg.parameters + b] += ja * jb * wi;
            }
        }
    }
    return generation;
}

inline auto linearAmplitudeSse(BracketNewtonRuntimeContext const& context,
                              std::vector<double> const& amplitudes) -> double {
    if (context.sampleT.size() != context.targetX.size() ||
        context.linearBasisX.size() != context.sampleT.size() * amplitudes.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double sse = 0.0;
    for (std::size_t i = 0; i < context.sampleT.size(); ++i) {
        double prediction = 0.0;
        for (std::size_t j = 0; j < amplitudes.size(); ++j) {
            prediction += context.linearBasisX[i * amplitudes.size() + j] * amplitudes[j];
        }
        double const error = context.targetX[i] - prediction;
        sse += error * error;
    }
    return std::isfinite(sse) ? sse : std::numeric_limits<double>::infinity();
}

inline auto applyMappedAmplitudeUpdate(std::vector<double>& coefficients,
                                       std::vector<double> const& proposed,
                                       BracketNewtonRuntimeContext const& context) -> bool {
    if (coefficients.size() != proposed.size() || proposed.empty() ||
        !validRuntimeContext(context, proposed.size()) ||
        context.linearBasisX.size() != context.sampleT.size() * proposed.size()) {
        return false;
    }
    for (auto kind : context.parameterMap.kinds) {
        if (kind != BracketNewtonParameterKind::amplitude) { return false; }
    }
    for (double value : proposed) {
        if (!std::isfinite(value) || std::abs(value) > 1.0e6) { return false; }
    }
    double const before = linearAmplitudeSse(context, coefficients);
    double const after = linearAmplitudeSse(context, proposed);
    if (!std::isfinite(before) || !std::isfinite(after) || after > before) { return false; }
    coefficients = proposed;
    return true;
}

template <typename Score>
inline auto applyMappedNonlinearDelta(std::vector<double>& coefficients,
                                      std::vector<double> const& delta,
                                      BracketNewtonRuntimeContext const& context,
                                      Score&& score,
                                      double maxAbsCoefficient = 1.0e6) -> bool {
    if (coefficients.size() != delta.size() || delta.empty() ||
        !validRuntimeContext(context, delta.size()) ||
        !(maxAbsCoefficient > 0.0) || !std::isfinite(maxAbsCoefficient)) return false;

    std::vector<double> trial = coefficients;
    for (std::size_t i = 0; i < trial.size(); ++i) {
        if (!std::isfinite(trial[i]) || !std::isfinite(delta[i])) return false;
        trial[i] += delta[i];
        if (!std::isfinite(trial[i]) || std::abs(trial[i]) > maxAbsCoefficient) return false;
    }

    double const before = std::invoke(score, coefficients);
    double const after = std::invoke(score, trial);
    if (!std::isfinite(before) || !std::isfinite(after) || after > before) return false;
    coefficients = std::move(trial);
    return true;
}

inline auto strictTraceRidge(std::vector<double> const& gram, std::size_t parameters,
                             double fraction = 1e-3) -> double {
    if (parameters == 0U || gram.size() < parameters * parameters ||
        !(fraction >= 0.0) || !std::isfinite(fraction)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double trace = 0.0;
    for (std::size_t i = 0; i < parameters; ++i) {
        double const diagonal = gram[i * parameters + i];
        if (!std::isfinite(diagonal) || diagonal < 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        trace += diagonal;
    }
    return std::isfinite(trace) ? fraction * trace
                                : std::numeric_limits<double>::quiet_NaN();
}

inline auto validGeneration(BracketNewtonConfig const& cfg,
                            BracketNewtonGeneration const& generation) -> bool {
    if (cfg.population == 0U || cfg.samples == 0U || cfg.parameters == 0U || cfg.sweeps == 0U ||
        generation.residual.size() != cfg.population * cfg.samples ||
        generation.gram.size() != cfg.parameters * cfg.parameters ||
        generation.rhs.size() != cfg.parameters) {
        return false;
    }
    if (cfg.warmStart && generation.parentPacingState.size() != cfg.population) return false;
    if (cfg.warmStart && !generation.parentIndex.has_value()) return false;
    if (cfg.useFp16Storage && (!generation.fp16Basis.empty() || !generation.fp16Residual.empty())) {
        if (generation.fp16Basis.size() != cfg.samples * cfg.parameters ||
            generation.fp16Residual.size() != cfg.samples) return false;
        Fp16GramConfig fp16Config{cfg.samples, cfg.parameters, true, true};
        Fp16GramBatch fp16Batch{generation.fp16Basis, generation.fp16Residual};
        if (!validFp16Gram(fp16Config, fp16Batch)) return false;
    }
    if (!generation.footpoints.segmentIndex.empty() || !generation.footpoints.segmentT.empty()) {
        if (!validFootpointState(generation.footpoints, cfg.samples, cfg.sweeps)) return false;
    }
    if (!cfg.warmStart && !generation.parentPacingState.empty()) return false;
    for (double value : generation.parentPacingState) {
        if (!std::isfinite(value)) return false;
    }
    for (double value : generation.residual) {
        if (!std::isfinite(value)) { return false; }
    }
    return std::isfinite(strictTraceRidge(generation.gram, cfg.parameters, cfg.ridgeFraction));
}

inline auto solveRidgeCpu(std::vector<double> const& gram, std::vector<double> const& rhs,
                          std::size_t parameters, double fraction = 1e-3)
    -> std::vector<double> {
    double const ridge = strictTraceRidge(gram, parameters, fraction);
    if (!std::isfinite(ridge) || rhs.size() != parameters ||
        gram.size() != parameters * parameters) { return {}; }
    std::vector<double> a(gram);
    std::vector<double> b(rhs);
    for (std::size_t i = 0; i < parameters; ++i) {
        a[i * parameters + i] += ridge;
    }
    for (std::size_t col = 0; col < parameters; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < parameters; ++row) {
            if (std::abs(a[row * parameters + col]) > std::abs(a[pivot * parameters + col])) {
                pivot = row;
            }
        }
        double const scale = std::abs(a[pivot * parameters + col]);
        if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::epsilon()) { return {}; }
        if (pivot != col) {
            for (std::size_t j = col; j < parameters; ++j) {
                std::swap(a[col * parameters + j], a[pivot * parameters + j]);
            }
            std::swap(b[col], b[pivot]);
        }
        for (std::size_t row = col + 1; row < parameters; ++row) {
            double const factor = a[row * parameters + col] / a[col * parameters + col];
            for (std::size_t j = col + 1; j < parameters; ++j) {
                a[row * parameters + j] -= factor * a[col * parameters + j];
            }
            b[row] -= factor * b[col];
        }
    }
    std::vector<double> x(parameters, 0.0);
    for (std::size_t i = parameters; i-- > 0;) {
        double value = b[i];
        for (std::size_t j = i + 1; j < parameters; ++j) {
            value -= a[i * parameters + j] * x[j];
        }
        x[i] = value / a[i * parameters + i];
        if (!std::isfinite(x[i])) { return {}; }
    }
    return x;
}

class BracketNewtonRuntime {
public:
    virtual ~BracketNewtonRuntime() = default;
    virtual auto submitGeneration(BracketNewtonConfig const& cfg,
                                  BracketNewtonGeneration const& generation)
        -> BracketNewtonResult = 0;
};

static auto submitGenerationWithFallback(
    BracketNewtonRuntime& runtime, BracketNewtonConfig const& cfg,
    std::vector<double> residual, std::vector<double> gram, std::vector<double> rhs,
    std::vector<double> parentPacingState = {},
    std::optional<std::size_t> parentIndex = {},
    BracketNewtonFootpointState footpoints = {},
    std::vector<BracketNewtonPoint> targetPoints = {},
    std::vector<BracketNewtonPoint> candidateCurveSamples = {}) -> BracketNewtonResult;

static auto submitGenerationWithFallback(
    BracketNewtonRuntime& runtime, BracketNewtonConfig const& cfg,
    BracketNewtonGeneration generation) -> BracketNewtonResult;

struct UmaBracketNewtonResult {
    bool ownershipOk{false};
    BracketNewtonResult transaction{};
};

class UmaBracketNewtonRunner {
public:
    explicit UmaBracketNewtonRunner(BracketNewtonRuntime& runtime) : runtime_(runtime) {}

    auto submit(std::size_t slot, BracketNewtonConfig const& cfg,
                std::vector<double> residual, std::vector<double> gram,
                std::vector<double> rhs, std::vector<double> parentPacingState = {},
                std::optional<std::size_t> parentIndex = {},
                BracketNewtonFootpointState footpoints = {},
                std::vector<BracketNewtonPoint> targetPoints = {},
                std::vector<BracketNewtonPoint> candidateCurveSamples = {})
        -> UmaBracketNewtonResult {
        if (!submitAsync(slot, cfg, std::move(residual), std::move(gram), std::move(rhs),
                         std::move(parentPacingState), parentIndex, std::move(footpoints),
                         std::move(targetPoints), std::move(candidateCurveSamples))) {
            UmaBracketNewtonResult out;
            out.transaction.usedFallback = true;
            return out;
        }
        return complete(slot);
    }

    auto submitAsync(std::size_t slot, BracketNewtonConfig const& cfg,
                     std::vector<double> residual, std::vector<double> gram,
                     std::vector<double> rhs, std::vector<double> parentPacingState = {},
                     std::optional<std::size_t> parentIndex = {},
                     BracketNewtonFootpointState footpoints = {},
                     std::vector<BracketNewtonPoint> targetPoints = {},
                     std::vector<BracketNewtonPoint> candidateCurveSamples = {}) -> bool {
        if (slot >= UmaPingPongContract::slotCount || !ownership_.acquireWritable(slot) ||
            !ownership_.markSubmitted(slot)) return false;
        pending_[slot] = Pending{cfg, std::move(residual), std::move(gram), std::move(rhs),
                                 std::move(parentPacingState), parentIndex, std::move(footpoints),
                                 std::move(targetPoints), std::move(candidateCurveSamples)};
        return true;
    }

    auto complete(std::size_t slot) -> UmaBracketNewtonResult {
        if (slot >= UmaPingPongContract::slotCount ||
            ownership_.state(slot) != UmaSlotState::submitted || !pending_[slot].has_value()) {
            UmaBracketNewtonResult out;
            out.transaction.usedFallback = true;
            return out;
        }
        auto pending = std::move(*pending_[slot]);
        pending_[slot].reset();
        auto result = submitGenerationWithFallback(runtime_, pending.cfg,
            std::move(pending.residual), std::move(pending.gram), std::move(pending.rhs),
            std::move(pending.parentPacingState), pending.parentIndex, pending.footpoints,
            pending.targetPoints, pending.candidateCurveSamples);
        return finish(slot, std::move(result));
    }

    auto completeFromBackend(std::size_t slot, BracketNewtonResult backendResult)
        -> UmaBracketNewtonResult {
        if (slot >= UmaPingPongContract::slotCount ||
            ownership_.state(slot) != UmaSlotState::submitted || !pending_[slot].has_value()) {
            UmaBracketNewtonResult out;
            out.transaction.usedFallback = true;
            return out;
        }
        auto pending = std::move(*pending_[slot]);
        pending_[slot].reset();
        bool finiteUpdate = backendResult.update.size() == pending.cfg.parameters;
        for (double value : backendResult.update) finiteUpdate = finiteUpdate && std::isfinite(value);
        if (!backendResult.ok || backendResult.usedFallback || !finiteUpdate) {
            backendResult.update = solveRidgeCpu(pending.gram, pending.rhs, pending.cfg.parameters,
                                                 pending.cfg.ridgeFraction);
            backendResult.ok = !backendResult.update.empty();
            backendResult.usedFallback = true;
            backendResult.ridge = strictTraceRidge(pending.gram, pending.cfg.parameters,
                                                   pending.cfg.ridgeFraction);
        }
        if (pending.cfg.warmStart && backendResult.warmStartStatus == WarmStartStatus::NotRequested)
            backendResult.warmStartStatus = WarmStartStatus::BlockedMissingCorrespondenceState;
        if (backendResult.nextPacingState.empty() && pending.cfg.warmStart)
            backendResult.nextPacingState = pending.parentPacingState;
        if (backendResult.nextFootpoints.segmentIndex.empty() &&
            !pending.footpoints.segmentIndex.empty()) {
            backendResult.nextFootpoints = pending.footpoints;
        }
        return finish(slot, std::move(backendResult));
    }

    auto slotState(std::size_t slot) const -> UmaSlotState { return ownership_.state(slot); }
    auto canPrepareNext(std::size_t currentSlot, std::size_t nextSlot) const -> bool {
        return ownership_.canPrepareNext(currentSlot, nextSlot);
    }

private:
    auto finish(std::size_t slot, BracketNewtonResult result) -> UmaBracketNewtonResult {
        UmaBracketNewtonResult out;
        out.transaction = std::move(result);
        out.ownershipOk = ownership_.markCompleted(slot);
        return out;
    }

    struct Pending {
        BracketNewtonConfig cfg;
        std::vector<double> residual;
        std::vector<double> gram;
        std::vector<double> rhs;
        std::vector<double> parentPacingState;
        std::optional<std::size_t> parentIndex;
        BracketNewtonFootpointState footpoints;
        std::vector<BracketNewtonPoint> targetPoints;
        std::vector<BracketNewtonPoint> candidateCurveSamples;
    };
    BracketNewtonRuntime& runtime_;
    UmaPingPongContract ownership_;
    std::array<std::optional<Pending>, UmaPingPongContract::slotCount> pending_{};
};

class MetalBracketNewtonRuntime final : public BracketNewtonRuntime {
public:
    MetalBracketNewtonRuntime();
    ~MetalBracketNewtonRuntime() override;
    MetalBracketNewtonRuntime(MetalBracketNewtonRuntime const&) = delete;
    auto operator=(MetalBracketNewtonRuntime const&) -> MetalBracketNewtonRuntime& = delete;
    auto available() const -> bool;
    auto submitGeneration(BracketNewtonConfig const& cfg,
                          BracketNewtonGeneration const& generation)
        -> BracketNewtonResult override;
#if defined(__APPLE__)
    auto submitGenerationAsync(BracketNewtonConfig const& cfg,
                               BracketNewtonGeneration const& generation,
                               std::function<void(BracketNewtonResult)> completion) -> bool;
#endif

private:
    void* impl_{nullptr};
};

static auto submitGenerationWithFallback(
    BracketNewtonRuntime& runtime, BracketNewtonConfig const& cfg,
    BracketNewtonGeneration generation) -> BracketNewtonResult {
    BracketNewtonResult result;
    result.ridge = strictTraceRidge(generation.gram, cfg.parameters, cfg.ridgeFraction);
    if (!validGeneration(cfg, generation)) {
        result.usedFallback = true;
        return result;
    }

    std::vector<double> fallbackGram = generation.gram;
    std::vector<double> fallbackRhs = generation.rhs;
    if (cfg.useFp16Storage && !generation.fp16Basis.empty()) {
        Fp16GramConfig fp16Config{cfg.samples, cfg.parameters, true, true};
        Fp16GramBatch fp16Batch{generation.fp16Basis, generation.fp16Residual};
        auto fp16Accumulated = accumulateFp16GramCpu(fp16Config, fp16Batch);
        if (!fp16Accumulated.ok) {
            result.usedFallback = true;
            return result;
        }
        fallbackGram = std::move(fp16Accumulated.gram);
        fallbackRhs = std::move(fp16Accumulated.rhs);
    }

    result = runtime.submitGeneration(cfg, generation);
    if (cfg.warmStart && result.warmStartStatus == WarmStartStatus::NotRequested) {
        result.warmStartStatus = WarmStartStatus::BlockedMissingCorrespondenceState;
    }
    if (result.nextFootpoints.segmentIndex.empty() && !generation.footpoints.segmentIndex.empty())
        result.nextFootpoints = generation.footpoints;
    if (!generation.targetPoints.empty() && !generation.candidateCurveSamples.empty()) {
        auto recomputed = cfg.useGeneralPathCorrespondence
            ? generalPathCorrespondence(cfg, generation.targetPoints,
                                        generation.candidateCurveSamples,
                                        generation.footpoints.completedSweeps)
            : recomputeFootpointsCpu(generation.targetPoints,
                                     generation.candidateCurveSamples,
                                     cfg.sweeps,
                                     generation.footpoints.completedSweeps);
        if (recomputed.has_value()) result.nextFootpoints = *recomputed;
    }
    if (!result.nextFootpoints.segmentIndex.empty())
        result.completedSweeps = result.nextFootpoints.completedSweeps;
    if (result.nextPacingState.empty() && cfg.warmStart) {
        result.nextPacingState = generation.parentPacingState;
    }
    if (result.ok && !result.usedFallback && result.update.size() == cfg.parameters) { return result; }
    result.update = solveRidgeCpu(fallbackGram, fallbackRhs, cfg.parameters, cfg.ridgeFraction);
    result.nextPacingState = generation.parentPacingState;
    if (cfg.warmStart) {
        result.warmStartStatus = WarmStartStatus::BlockedMissingCorrespondenceState;
    }
    result.ridge = strictTraceRidge(fallbackGram, cfg.parameters, cfg.ridgeFraction);
    result.ok = !result.update.empty();
    result.usedFallback = true;
    return result;
}

static auto submitGenerationWithFallback(
    BracketNewtonRuntime& runtime, BracketNewtonConfig const& cfg,
    std::vector<double> residual, std::vector<double> gram, std::vector<double> rhs,
    std::vector<double> parentPacingState, std::optional<std::size_t> parentIndex,
    BracketNewtonFootpointState footpoints,
    std::vector<BracketNewtonPoint> targetPoints,
    std::vector<BracketNewtonPoint> candidateCurveSamples)
    -> BracketNewtonResult {
    BracketNewtonGeneration generation;
    generation.residual = std::move(residual);
    generation.gram = std::move(gram);
    generation.rhs = std::move(rhs);
    generation.parentPacingState = std::move(parentPacingState);
    generation.parentIndex = parentIndex;
    generation.footpoints = std::move(footpoints);
    generation.targetPoints = std::move(targetPoints);
    generation.candidateCurveSamples = std::move(candidateCurveSamples);
    return submitGenerationWithFallback(runtime, cfg, std::move(generation));
}

}
