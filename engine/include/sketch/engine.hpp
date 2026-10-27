#ifndef SKETCH_RESEARCH_ENGINE_HPP
#define SKETCH_RESEARCH_ENGINE_HPP


#include "sketch/multitree_gp.hpp"
#include "sketch/portfolio.hpp"
#include "sketch/selection_core.hpp"
#include "sketch/pair_individual.hpp"
#include "sketch/sf2_guard.hpp"

#include <cmath>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

namespace Sketch::research {

struct ExpressionCandidate {
    int nodes{0};
    float shapeErr{0.0F};
    float dlBits{0.0F};
    float r2{0.0F};
    std::string exprX;
    std::string exprY;
    std::vector<float> curveX;
    std::vector<float> curveY;
};

struct EngineResult {
    int status{0};
    std::vector<ExpressionCandidate> front;
    std::size_t championIndex{0};
    int championRule{0};
    int generationsUsed{0};
};

namespace detail {

inline void applyGpuPointSegEnv(Sketch::MultitreeConfig& c)
{
    char const* v = std::getenv("SKETCH_GPU_POINTSEG");
    if (v != nullptr && std::string_view(v) == "0") { c.useBatchedGpuPointSeg = false; }
}

inline void applyGradMutEnv(Sketch::MultitreeConfig& c)
{
    if (char const* v = std::getenv("SKETCH_GRADMUT_MODE")) { c.gradientMutationMode = std::atoi(v); }
    if (std::getenv("SKETCH_GRADMUT_INRUN") != nullptr) { c.gradientMutationInRunErr = true; }
    if (std::getenv("SKETCH_G6E") != nullptr) { c.gradientMutationG6e = true; }
}

inline int applySolveThreadsEnv(int defaultThreads)
{
    char const* v = std::getenv("SKETCH_SOLVE_THREADS");
    if (v == nullptr) { return defaultThreads; }
    int const n = std::atoi(v);
    return n > 0 ? n : defaultThreads;
}

}

template <typename DTable>
EngineResult solveStrokeStub(
    DTable& dtable,
    std::vector<Operon::Scalar> const& dx,
    std::vector<Operon::Scalar> const& dy,
    Sketch::NoiseModel const& nativeNoise,
    double,
    int pop,
    int gens,
    int max_length,
    float floorClamp = 0.0F,
    float floorStop = 0.0F,
    int floorSimplifyGens = 0)
{
    EngineResult out;

    Sketch::MultitreeConfig base;
    base.populationSize = pop > 0 ? static_cast<std::size_t>(pop) : 500;
    base.generations = gens > 0 ? static_cast<std::size_t>(gens) : 60;
    int const sf2 = Sketch::sf2Mask();
    if (Sketch::sf2On(sf2, 15)) {
        constexpr double kBudgetS = 5.0;
        constexpr double kMeasuredWallAt60S = 4.8975;
        base.generations = static_cast<std::size_t>(std::floor(60.0 * kBudgetS / kMeasuredWallAt60S));
    }
    base.maxLength = max_length > 0 ? static_cast<std::size_t>(max_length) : 28;
    base.optFraction = 0.15;
    base.gridSize = 200;
    base.tournamentSize = 7;
    base.crossTreeProb = 0.3;
    base.orderWeight = 1.0;
    base.useAccelerateTranscendentals = true;
    base.useGridSdf = true;
    base.annealLowN = 64;
    base.annealResolution = 0.75;
    base.annealLowM = 64;
    detail::applyGpuPointSegEnv(base);
    detail::applyGradMutEnv(base);

    std::vector<Sketch::MultitreeConfig> members;
    auto c0 = base; c0.seed = 42;
    members.push_back(c0);
    auto c1 = base; c1.seed = 43; c1.crossTreeProb = 0.6;
    members.push_back(c1);
    auto c2 = base; c2.seed = 44; c2.perTreeCredit = true;
    members.push_back(c2);
    auto c3 = base; c3.seed = 45;
    c3.floorClamp = nativeNoise.empty() ? 0.0004f : static_cast<Operon::Scalar>(floorClamp);
    members.push_back(c3);
    auto c4 = base; c4.seed = 46;
    members.push_back(c4);
    auto c5 = base; c5.seed = 47; c5.perTreeCredit = true;
    c5.floorClamp = nativeNoise.empty() ? 0.0004f : static_cast<Operon::Scalar>(floorClamp);
    members.push_back(c5);

    if (!nativeNoise.empty()) {
        auto const fStop = static_cast<Operon::Scalar>(floorStop);
        auto const fSimplifyGens = floorSimplifyGens > 0 ? static_cast<std::size_t>(floorSimplifyGens) : 0U;
        for (auto& member : members) {
            member.floorStop = fStop;
            member.floorSimplifyGens = fSimplifyGens;
        }
    }

    bool const sf2Hw = Sketch::sf2On(sf2, 14);
    auto portfolio_results = Sketch::runPortfolio<DTable>(dx, dy, members, nativeNoise,
detail::applySolveThreadsEnv(sf2Hw ? 0 : 6),
true,
sf2Hw ? base.generations / 4 : 10,
sf2Hw ? base.generations / 2 : 25);

    std::vector<Sketch::FrontEntry> pooled;
    std::vector<int> pooledSrc;
    for (std::size_t mi = 0; mi < portfolio_results.size(); ++mi) {
        auto& pr = portfolio_results[mi];
        if (!pr.ok) continue;
        for (auto& fe : pr.result.front) {
            pooled.push_back(std::move(fe));
            pooledSrc.push_back(static_cast<int>(mi));
        }
    }

    constexpr int kOneRulerSelectDefault = 2;
    Sketch::SelectionOptions selOpt;
    char const* orEnv = std::getenv("SKETCH_ONE_RULER_SELECT");
    selOpt.rulerMode = (orEnv == nullptr) ? kOneRulerSelectDefault
                     : (std::string_view(orEnv) == "0") ? 0
                     : (std::string_view(orEnv) == "2") ? 2 : 1;
    char const* k1LogEnv = std::getenv("SKETCH_K1_MEMBER_LOG");
    selOpt.memberLog = (k1LogEnv != nullptr) && (std::string_view(k1LogEnv) == "1");
    char const* nfsEnv = std::getenv("SKETCH_NOISE_FLOOR_SELECT");
    selOpt.noiseFloorSelect = (nfsEnv == nullptr) || (std::string_view(nfsEnv) != "0");
    selOpt.adequacyBand = base.championAdequacyBand;

    auto sel = Sketch::selectFromPool<DTable>(dtable, std::move(pooled), pooledSrc, members, dx, dy,
                                              nativeNoise, selOpt);
    out.status = sel.status;
    if (sel.status != 0) { return out; }

    auto const& skyline = sel.skyline;
    out.front.reserve(skyline.size());
    for (std::size_t i = 0; i < skyline.size(); ++i) {
        auto const& e = skyline[i];
        ExpressionCandidate c;
        c.nodes = e.nodes;
        c.shapeErr = sel.stats[i].shapeErr;
        c.r2 = sel.stats[i].r2;
        c.dlBits = sel.stats[i].dlBits;
        c.exprX = e.exprX;
        c.exprY = e.exprY;
        c.curveX = e.curveX;
        c.curveY = e.curveY;
        out.front.push_back(std::move(c));
    }
    out.championIndex = sel.championIndex;
    out.championRule = sel.championRule;
    out.generationsUsed = static_cast<int>(base.generations);
    return out;
}

}

#endif
