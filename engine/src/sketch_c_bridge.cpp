#include "sketch/sketch_c_bridge.h"
#include "sketch/multitree_gp.hpp"
#include "sketch/portfolio.hpp"
#include "sketch/selection_core.hpp"
#include "sketch/sigma_v2_native.hpp"
#include "sketch/engine.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <atomic>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>

static void applyGpuPointSegEnv(Sketch::MultitreeConfig& c)
{
    char const* v = std::getenv("SKETCH_GPU_POINTSEG");
    if (v != nullptr && std::string_view(v) == "0") { c.useBatchedGpuPointSeg = false; }
}

static int applySolveThreadsEnv(int defaultThreads)
{
    char const* v = std::getenv("SKETCH_SOLVE_THREADS");
    if (v == nullptr) { return defaultThreads; }
    int const n = std::atoi(v);
    return n > 0 ? n : defaultThreads;
}

static void applyGradMutEnv(Sketch::MultitreeConfig& c)
{
    if (char const* v = std::getenv("SKETCH_GRADMUT_MODE")) { c.gradientMutationMode = std::atoi(v); }
    if (std::getenv("SKETCH_GRADMUT_INRUN") != nullptr) { c.gradientMutationInRunErr = true; }
    if (std::getenv("SKETCH_G6E") != nullptr) { c.gradientMutationG6e = true; }
}

static bool researchEngineOn()
{
    char const* v = std::getenv("SKETCH_RESEARCH_ENGINE");
    return v != nullptr && std::string_view(v) == "1";
}

static bool resampleAndNormalize(
    const float* in_x, const float* in_y, int n_in, int target_n,
    std::vector<Operon::Scalar>& out_x, std::vector<Operon::Scalar>& out_y)
{
    return Sketch::resampleAndNormalizeStroke(in_x, in_y, n_in, target_n, out_x, out_y);
}

extern "C" int sketch_solve_stroke_c(
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int max_length,
    SketchBridgeCallback callback,
    void* user_data)
{
    return sketch_solve_stroke_with_noise_c(
        pts_x, pts_y, n_pts, pop, gens, max_length, nullptr, callback, user_data);
}

extern "C" int sketch_solve_stroke_with_noise_c(
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int max_length,
    const SketchNoiseField* noise_field,
    SketchBridgeCallback callback,
    void* user_data)
{
    using DTable = Operon::DispatchTable<Operon::Scalar>;
    DTable dtable;

    std::vector<Operon::Scalar> dx, dy;
    if (!resampleAndNormalize(pts_x, pts_y, n_pts, 120, dx, dy)) {
        return -1;
    }

    auto fillStats = [&dx, &dy](SketchBridgeCandidate& c, Sketch::FrontEntry const& e) {
        auto const st = Sketch::candidateStats(e, dx, dy);
        c.shape_err = st.shapeErr;
        c.r2 = st.r2;
        c.dl_bits = st.dlBits;
    };

    Sketch::NoiseModel nativeNoise;
    constexpr int kResampledPoints = 120;
    if (noise_field != nullptr) {
        if (noise_field->abi_version != SKETCH_NOISE_FIELD_ABI_VERSION ||
            noise_field->point_s2 == nullptr ||
            noise_field->point_count < 1 || noise_field->point_count > kResampledPoints) {
            return -3;
        }
        nativeNoise.s2.reserve(static_cast<std::size_t>(noise_field->point_count));
        nativeNoise.nu.reserve(static_cast<std::size_t>(noise_field->point_count));
        bool const haveLogVar = noise_field->point_log_var != nullptr;
        if (haveLogVar) { nativeNoise.logVar.reserve(static_cast<std::size_t>(noise_field->point_count)); }
        for (int i = 0; i < noise_field->point_count; ++i) {
            auto const s2 = noise_field->point_s2[i];
            auto const nu = noise_field->point_nu != nullptr ? noise_field->point_nu[i] : 0.0F;
            if (!std::isfinite(s2) || s2 <= 0.0F || !std::isfinite(nu)) { return -3; }
            nativeNoise.s2.push_back(static_cast<Operon::Scalar>(s2));
            nativeNoise.nu.push_back(static_cast<Operon::Scalar>(nu));
            if (haveLogVar) {
                auto const lv = noise_field->point_log_var[i];
                if (!std::isfinite(lv)) { return -3; }
                nativeNoise.logVar.push_back(static_cast<Operon::Scalar>(lv));
            }
        }
    }

    auto t0 = std::chrono::steady_clock::now();

    auto emitFront = [&](int gensLabel, std::vector<Sketch::FrontEntry> const& front) {
        if (!callback || front.empty()) { return; }
        auto now = std::chrono::steady_clock::now();
        int elapsed_ms = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count());
        std::vector<SketchBridgeCandidate> cands;
        cands.reserve(front.size());
        for (const auto& e : front) {
            SketchBridgeCandidate c;
            c.nodes = e.nodes;
            fillStats(c, e);
            c.expr_x = e.exprX.c_str();
            c.expr_y = e.exprY.c_str();
            c.curve_count = static_cast<int>(e.curveX.size());
            c.curve_x = e.curveX.data();
            c.curve_y = e.curveY.data();
            cands.push_back(c);
        }
        std::sort(cands.begin(), cands.end(), [](const auto& a, const auto& b) {
            if (std::abs(a.shape_err - b.shape_err) > 1e-7f) { return a.shape_err < b.shape_err; }
            return a.nodes < b.nodes;
        });
        callback(gensLabel, elapsed_ms, &cands.front(), cands.data(), static_cast<int>(cands.size()), user_data);
    };

    {
        Sketch::MultitreeConfig scoutCfg;
        scoutCfg.populationSize = 40;
        scoutCfg.generations = 12;
        scoutCfg.maxLength = 20;
        scoutCfg.gridSize = 120;
        scoutCfg.tournamentSize = 7;
        scoutCfg.crossTreeProb = 0.3;
        scoutCfg.orderWeight = 1.0;
        scoutCfg.useAccelerateTranscendentals = true;
        scoutCfg.useGridSdf = true;
        applyGpuPointSegEnv(scoutCfg);

        Sketch::MultitreeGP<DTable> scoutGp(dtable, dx, dy, scoutCfg);
        auto scoutRes = scoutGp.Run();
        if (callback && !scoutRes.front.empty()) {
            auto now = std::chrono::steady_clock::now();
            int elapsed_ms = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count());

            std::vector<SketchBridgeCandidate> cands;
            for (const auto& e : scoutRes.front) {
                SketchBridgeCandidate c;
                c.nodes = e.nodes;
                fillStats(c, e);
                c.expr_x = e.exprX.c_str();
                c.expr_y = e.exprY.c_str();
                c.curve_count = static_cast<int>(e.curveX.size());
                c.curve_x = e.curveX.data();
                c.curve_y = e.curveY.data();
                cands.push_back(c);
            }
            std::sort(cands.begin(), cands.end(), [](const auto& a, const auto& b) {
                if (std::abs(a.shape_err - b.shape_err) > 1e-7f) {
                    return a.shape_err < b.shape_err;
                }
                return a.nodes < b.nodes;
            });
            callback(12, elapsed_ms, &cands.front(), cands.data(), static_cast<int>(cands.size()), user_data);
        }
    }

    if (callback) {
        Sketch::MultitreeConfig rc;
        rc.populationSize = 200;
        rc.generations = 30;
        rc.maxLength = max_length > 0 ? static_cast<std::size_t>(max_length) : 28;
        rc.gridSize = 160;
        rc.tournamentSize = 7;
        rc.crossTreeProb = 0.3;
        rc.orderWeight = 1.0;
        rc.optFraction = 0.15;
        rc.useAccelerateTranscendentals = true;
        rc.useGridSdf = true;
        applyGpuPointSegEnv(rc);
        Sketch::MultitreeGP<DTable> rg(dtable, dx, dy, rc);
        auto rr = rg.Run();
        emitFront(30, rr.front);
    }

    std::vector<SketchBridgeCandidate> candidates;
    std::size_t championIndex = 0;
    std::size_t generationsForCallback = 0;
    int t13ChampionRule = 0;
    bool t13InsertedOffSkyline = false;
    std::size_t t13NNmlAbsOk = 0;
    std::size_t t13NNmlAbsMissing = 0;
    std::vector<Sketch::FrontEntry> skylineOwner;
    std::vector<Sketch::research::ExpressionCandidate> researchFrontOwner;

    if (researchEngineOn()) {
        auto res = Sketch::research::solveStrokeStub<DTable>(
            dtable, dx, dy, nativeNoise, 5.0, pop, gens, max_length,
            noise_field != nullptr ? noise_field->floor_clamp : 0.0F,
            noise_field != nullptr ? noise_field->floor_stop : 0.0F,
            noise_field != nullptr ? noise_field->floor_simplify_gens : 0);
        if (res.status != 0) {
            return res.status;
        }
        researchFrontOwner = std::move(res.front);
        candidates.reserve(researchFrontOwner.size());
        for (auto const& e : researchFrontOwner) {
            SketchBridgeCandidate c;
            c.nodes = e.nodes;
            c.shape_err = e.shapeErr;
            c.r2 = e.r2;
            c.dl_bits = e.dlBits;
            c.expr_x = e.exprX.c_str();
            c.expr_y = e.exprY.c_str();
            c.curve_count = static_cast<int>(e.curveX.size());
            c.curve_x = e.curveX.data();
            c.curve_y = e.curveY.data();
            candidates.push_back(c);
        }
        championIndex = res.championIndex;
        generationsForCallback = static_cast<std::size_t>(res.generationsUsed);
        t13ChampionRule = res.championRule;
    } else {
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
    applyGpuPointSegEnv(base);
    applyGradMutEnv(base);

    std::vector<Sketch::MultitreeConfig> members;
    auto c0 = base; c0.seed = 42;
    members.push_back(c0);
    auto c1 = base; c1.seed = 43; c1.crossTreeProb = 0.6;
    members.push_back(c1);
    auto c2 = base; c2.seed = 44; c2.perTreeCredit = true;
    members.push_back(c2);
    auto c3 = base; c3.seed = 45;
    c3.floorClamp = nativeNoise.empty()
        ? 0.0004f
        : static_cast<Operon::Scalar>(noise_field->floor_clamp);
    members.push_back(c3);
    auto c4 = base; c4.seed = 46;
    members.push_back(c4);
    auto c5 = base; c5.seed = 47; c5.perTreeCredit = true;
    c5.floorClamp = nativeNoise.empty()
        ? 0.0004f
        : static_cast<Operon::Scalar>(noise_field->floor_clamp);
    members.push_back(c5);

    if (!nativeNoise.empty()) {
        auto const floorStop = static_cast<Operon::Scalar>(noise_field->floor_stop);
        auto const simplifyGens = noise_field->floor_simplify_gens > 0
            ? static_cast<std::size_t>(noise_field->floor_simplify_gens) : 0U;
        for (auto& member : members) {
            member.floorStop = floorStop;
            member.floorSimplifyGens = simplifyGens;
        }
    }

    bool const sf2Hw = Sketch::sf2On(sf2, 14);
    auto portfolio_results = Sketch::runPortfolio<DTable>(dx, dy, members, nativeNoise,
                                                         applySolveThreadsEnv(sf2Hw ? 0 : 6),
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
    if (sel.status != 0) {
        return sel.status;
    }
    skylineOwner = std::move(sel.skyline);

    candidates.reserve(skylineOwner.size());
    for (std::size_t i = 0; i < skylineOwner.size(); ++i) {
        auto const& e = skylineOwner[i];
        SketchBridgeCandidate c;
        c.nodes = e.nodes;
        c.shape_err = sel.stats[i].shapeErr;
        c.r2 = sel.stats[i].r2;
        c.dl_bits = sel.stats[i].dlBits;
        c.expr_x = e.exprX.c_str();
        c.expr_y = e.exprY.c_str();
        c.curve_count = static_cast<int>(e.curveX.size());
        c.curve_x = e.curveX.data();
        c.curve_y = e.curveY.data();
        candidates.push_back(c);
    }
    championIndex = sel.championIndex;
    generationsForCallback = base.generations;
    t13ChampionRule = sel.championRule;
    t13InsertedOffSkyline = sel.championInsertedOffSkyline;
    t13NNmlAbsOk = sel.nNmlAbsOk;
    t13NNmlAbsMissing = sel.nNmlAbsMissing;
    }

    auto now = std::chrono::steady_clock::now();
    int elapsed_ms = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count());

    const SketchBridgeCandidate* champ = &candidates[championIndex];

    if (std::getenv("SKETCH_T13_DIAG") != nullptr) {
        std::fprintf(stderr,
                     "[t13-diag] championRule=%d insertedOffSkyline=%d nNmlAbsOk=%zu nNmlAbsMissing=%zu "
                     "champNodes=%d champDlBits=%.6f champExprX=%s champExprY=%s\n",
                     t13ChampionRule, t13InsertedOffSkyline ? 1 : 0, t13NNmlAbsOk, t13NNmlAbsMissing,
                     champ->nodes, static_cast<double>(champ->dl_bits), champ->expr_x, champ->expr_y);
    }

    if (callback) {
        callback(static_cast<int>(generationsForCallback), elapsed_ms, champ, candidates.data(), static_cast<int>(candidates.size()), user_data);
    }

    return 0;
}

extern "C" int sketch_gpu_pointseg_counters_c(long long* out4)
{
    auto& c = Sketch::pointSegGpuCounters();
    if (out4 != nullptr) {
        out4[0] = c.gpuBatches.load(std::memory_order_relaxed);
        out4[1] = c.gpuCandidates.load(std::memory_order_relaxed);
        out4[2] = c.cpuInGpuBatch.load(std::memory_order_relaxed);
        out4[3] = c.fallbackBatches.load(std::memory_order_relaxed);
    }
    return Sketch::kPointSegGpuCompiledIn ? 1 : 0;
}

extern "C" void sketch_gpu_pointseg_counters_reset_c(void)
{
    auto& c = Sketch::pointSegGpuCounters();
    c.gpuBatches.store(0);
    c.gpuCandidates.store(0);
    c.cpuInGpuBatch.store(0);
    c.fallbackBatches.store(0);
}

extern "C" void sketch_op_diag_print_c(void)
{
    Sketch::printOpDiagTable(stderr);
}

extern "C" void sketch_nsafety_counters_print_c(void)
{
    auto& c = Sketch::nSafetyCounters();
    std::fprintf(stderr, "[op1-nsafety] nonFiniteMutatedConstRejected=%lld nonFiniteLmConstRejected=%lld\n",
                 c.nonFiniteMutatedConst.load(std::memory_order_relaxed),
                 c.nonFiniteLmConst.load(std::memory_order_relaxed));
}

extern "C" void sketch_op2_diag_print_c(void)
{
    Sketch::printOp2DiagTable(stderr);
}

extern "C" void sketch_e10_diag_print_c(void)
{
    auto& g48 = Sketch::cisGrid48Counters();
    std::fprintf(stderr,
        "[e10-cisgrid48] badStroke=%lld missingModel=%lld cellMismatch=%lld residualFail=%lld "
        "bit1NoBit0=%lld ensembleSkipped=%lld absFallback=%lld law1FallbackToLaw3=%lld\n",
        g48.badStroke.load(std::memory_order_relaxed), g48.missingModel.load(std::memory_order_relaxed),
        g48.cellMismatch.load(std::memory_order_relaxed), g48.residualFail.load(std::memory_order_relaxed),
        g48.bit1NoBit0.load(std::memory_order_relaxed), g48.ensembleSkipped.load(std::memory_order_relaxed),
        g48.absFallback.load(std::memory_order_relaxed), g48.law1FallbackToLaw3.load(std::memory_order_relaxed));
    auto& gauge = Sketch::op2GaugeCounters();
    long long gaugeTotal = 0;
    for (auto const& bucket : gauge.nullDimHist) { gaugeTotal += bucket.load(std::memory_order_relaxed); }
    std::fprintf(stderr, "[e10-op2gauge] solvesWithGaugeProject=%lld maxAbsConst=%.6g\n",
                 gaugeTotal, Sketch::op2MaxAbsConst());
    auto& ps = Sketch::pointSegGpuCounters();
    std::fprintf(stderr,
        "[e10-pointseg] gpuBatches=%lld gpuCandidates=%lld cpuInGpuBatch=%lld fallbackBatches=%lld\n",
        ps.gpuBatches.load(std::memory_order_relaxed), ps.gpuCandidates.load(std::memory_order_relaxed),
        ps.cpuInGpuBatch.load(std::memory_order_relaxed), ps.fallbackBatches.load(std::memory_order_relaxed));
    auto& ns = Sketch::nSafetyCounters();
    std::fprintf(stderr, "[e10-nsafety] nonFiniteMutatedConstRejected=%lld nonFiniteLmConstRejected=%lld\n",
                 ns.nonFiniteMutatedConst.load(std::memory_order_relaxed), ns.nonFiniteLmConst.load(std::memory_order_relaxed));
    std::fprintf(stderr,
        "[e10-masks] SKETCH_SF2=0x%x SKETCH_OP2=0x%x SKETCH_CIS_GRID48=%d SKETCH_CIS_EXACT=%d "
        "SKETCH_CIS_ENSEMBLE=%d SKETCH_NML_ABS=%d SKETCH_T13=%d SKETCH_B3_GRID48_GPU=%d SKETCH_ONE_RESIDUAL=%d\n",
        Sketch::sf2Mask(), Sketch::op2Mask(), Sketch::cisGrid48ModeFromEnv(), Sketch::cisExactModeDefault(),
        Sketch::cisEnsembleModeFromEnv(), Sketch::nmlAbsModeFromEnv(), Sketch::sketchT13Mask(),
        Sketch::b3Grid48GpuFromEnv(), Sketch::oneResidualModeFromEnv());
}

namespace { std::atomic<int> gDiagOneResidualMode{-1}; std::atomic<long long> gDiagFpSeed{-1}; }
extern "C" void sketch_diag_set_one_residual_c(int mode) { gDiagOneResidualMode.store(mode); }
extern "C" void sketch_diag_set_fp_seed_c(long long seed) { gDiagFpSeed.store(seed); }

namespace {
auto diagFingerprintFit(const float* pts_x, const float* pts_y, int n_pts, int pop, int gens,
                        int n_threads, int gpu_mode) -> Sketch::MultitreeResult
{
    using DTable = Operon::DispatchTable<Operon::Scalar>;
    DTable dtable;
    std::vector<Operon::Scalar> x(pts_x, pts_x + n_pts), y(pts_y, pts_y + n_pts);
    Sketch::MultitreeConfig c;
    c.seed = 424242ULL;
    c.populationSize = static_cast<std::size_t>(pop);
    c.generations = static_cast<std::size_t>(gens);
    c.useCisFitness = true;
    c.useSiCisFitness = true;
    c.useGridSdf = true;
    c.useAccelerateTranscendentals = true;
    c.nThreads = static_cast<std::size_t>(n_threads);
    if (gpu_mode == 0) { c.useBatchedGpuPointSeg = false; }
    if (gpu_mode == 1) { c.useBatchedGpuPointSeg = true; }
    if (char const* pv = std::getenv("SKETCH_PROFILE_TIMING"); pv != nullptr && std::atoi(pv) != 0) { c.profileTiming = true; }
    if (int const om = gDiagOneResidualMode.load(); om >= 0) { c.oneResidualMode = om; }
    if (long long const sd = gDiagFpSeed.load(); sd >= 0) { c.seed = static_cast<std::uint64_t>(sd); }
    applyGradMutEnv(c);
    Sketch::MultitreeGP<DTable> gp(dtable, x, y, c);
    return gp.Run();
}
}

extern "C" int sketch_diag_multitree_fingerprint_c(
    const char* id, const float* pts_x, const float* pts_y, int n_pts, int pop, int gens,
    int n_threads, int gpu_mode, char* out, int out_cap)
{
    if (pts_x == nullptr || pts_y == nullptr || n_pts < 2 || pop < 1 || gens < 1 || n_threads < 0) { return -1; }
    auto r = diagFingerprintFit(pts_x, pts_y, n_pts, pop, gens, n_threads, gpu_mode);
    std::string text;
    auto appendf = [&text](char const* fmt, auto... args) {
        int const k = std::snprintf(nullptr, 0, fmt, args...);
        if (k <= 0) { return; }
        std::size_t const at = text.size();
        text.resize(at + static_cast<std::size_t>(k) + 1);
        std::snprintf(text.data() + at, static_cast<std::size_t>(k) + 1, fmt, args...);
        text.resize(at + static_cast<std::size_t>(k));
    };
    appendf("FP %s shapeErr=%.9g front=%zu descents=%zu iters=%zu losscalls=%zu x=%s y=%s\n",
                  id != nullptr ? id : "-", static_cast<double>(r.shapeError), r.front.size(),
                  static_cast<std::size_t>(r.nConstDescents), static_cast<std::size_t>(r.nConstIters),
                  static_cast<std::size_t>(r.nConstLossCalls), r.exprX.c_str(), r.exprY.c_str());
    std::fprintf(stderr, "[op2-lm] %s lmStopUnresolvableDamped=%zu lmStopRealizedSmall=%zu lmStopFixedPoint=%zu\n",
                 id != nullptr ? id : "-", static_cast<std::size_t>(r.nLmStopUnresolvableDamped),
                 static_cast<std::size_t>(r.nLmStopRealizedSmall), static_cast<std::size_t>(r.nLmStopFixedPoint));
    for (auto const& fe : r.front) {
        appendf("   front %.9g %s | %s\n", static_cast<double>(fe.shapeErr), fe.exprX.c_str(), fe.exprY.c_str());
    }
    if (char const* pv = std::getenv("SKETCH_PROFILE_TIMING"); pv != nullptr && std::atoi(pv) != 0) {
        std::fprintf(stderr, "[prof] %s gen=%.1f eval=%.1f lm=%.1f lmBuild=%.1f lmLoss=%.1f score=%.1f scoreDist=%.1f scoreCis=%.1f "
                             "cisPointSeg=%.1f cisNll=%.1f cisRev=%.1f cisGrid48=%.1f gpuPrepass=%.1f breed=%.1f select=%.1f parallelWall=%.1f ms\n",
                     id != nullptr ? id : "-", r.profGenMs, r.profEvalMs, r.profLmMs, r.profLmBuildMs, r.profLmLossMs, r.profScoreMs,
                     r.profScoreDistMs, r.profScoreCisMs, r.profCisPointSegMs, r.profCisNllMs, r.profCisReverseMs, r.profCisGrid48Ms,
                     r.profGpuPrepassMs, r.profBreedMs, r.profSelectMs, r.profParallelWallMs);
    }
    if (out != nullptr && out_cap > 0) {
        std::size_t const k = std::min(text.size(), static_cast<std::size_t>(out_cap - 1));
        std::memcpy(out, text.data(), k);
        out[k] = '\0';
    }
    return static_cast<int>(text.size());
}

extern "C" void sketch_diag_set_cis_exact_c(int mode) { Sketch::setCisExactMode(mode); }
extern "C" void sketch_diag_set_sf2_c(int mask) { Sketch::setSf2Mask(mask); }

extern "C" int sketch_diag_multitree_rulers_c(
    const char* id, const float* pts_x, const float* pts_y, int n_pts, int pop, int gens,
    int n_threads, char* out, int out_cap)
{
    if (pts_x == nullptr || pts_y == nullptr || n_pts < 2 || pop < 1 || gens < 1 || n_threads < 0) { return -1; }
    auto r = diagFingerprintFit(pts_x, pts_y, n_pts, pop, gens, n_threads, -1);
    std::vector<Operon::Scalar> dx(pts_x, pts_x + n_pts), dy(pts_y, pts_y + n_pts);
    Sketch::normalizeInPlace(dx, dy);
    constexpr double kInf = std::numeric_limits<double>::infinity();
    auto rmsFR = [&](Sketch::FrontEntry const& e, double& f, double& rv) {
        f = kInf; rv = kInf;
        std::size_t const M = std::min(e.curveX.size(), e.curveY.size());
        if (M < 2) { return; }
        std::vector<Sketch::Vec2> sp(dx.size()), cp(M);
        for (std::size_t i = 0; i < dx.size(); ++i) { sp[i] = Sketch::Vec2{dx[i], dy[i]}; }
        for (std::size_t j = 0; j < M; ++j) { cp[j] = Sketch::Vec2{e.curveX[j], e.curveY[j]}; }
        auto const a = Sketch::evaluatePointSegmentOne(sp, cp, false, false);
        auto const b = Sketch::evaluatePointSegmentOne(cp, sp, false, false);
        if (a.ok) { f = std::sqrt(static_cast<double>(a.fitness) / static_cast<double>(sp.size())); }
        if (b.ok) { rv = std::sqrt(static_cast<double>(b.fitness) / static_cast<double>(M)); }
    };
    double ch[5] = {kInf, kInf, kInf, kInf, kInf}, mn[5] = {kInf, kInf, kInf, kInf, kInf};
    std::vector<Operon::Scalar> none;
    double bestOwn = kInf;
    for (auto const& e : r.front) {
        double v[5];
        for (int m = 0; m < 3; ++m) { v[m] = Sketch::oneRulerSiCisShapeNll(e.curveX, e.curveY, dx, dy, none, m); }
        rmsFR(e, v[3], v[4]);
        for (int m = 0; m < 5; ++m) { if (std::isfinite(v[m]) && v[m] < mn[m]) { mn[m] = v[m]; } }
        if (static_cast<double>(e.shapeErr) < bestOwn) { bestOwn = static_cast<double>(e.shapeErr); for (int m = 0; m < 5; ++m) { ch[m] = v[m]; } }
    }
    char line[640];
    std::snprintf(line, sizeof line, "RUL %s champ %.9g %.9g %.9g %.9g %.9g | min %.9g %.9g %.9g %.9g %.9g | front %zu shapeErr %.9g\n",
                  id != nullptr ? id : "-", ch[0], ch[1], ch[2], ch[3], ch[4], mn[0], mn[1], mn[2], mn[3], mn[4], r.front.size(),
                  static_cast<double>(r.shapeError));
    std::string const text(line);
    if (out != nullptr && out_cap > 0) {
        std::size_t const k = std::min(text.size(), static_cast<std::size_t>(out_cap - 1));
        std::memcpy(out, text.data(), k);
        out[k] = '\0';
    }
    return static_cast<int>(text.size());
}

namespace {
struct SigmaV2Cache {
    std::mutex m;
    bool attempted = false;
    bool ok = false;
    sketch::SigmaV2Model model;
};
SigmaV2Cache& sigmaV2Cache()
{
    static SigmaV2Cache cache;
    return cache;
}
}

extern "C" int sketch_sigma_v2_cells_c(
    const float* x,
    const float* y,
    int n,
    double* lnE_out,
    double* sd_out,
    int cap)
{
    if (x == nullptr || y == nullptr || n < 2 || lnE_out == nullptr || sd_out == nullptr || cap <= 0) {
        return -3;
    }
    auto& cache = sigmaV2Cache();
    {
        std::lock_guard<std::mutex> lock(cache.m);
        if (!cache.attempted) {
            cache.attempted = true;
            char const* envPath = std::getenv("SKETCH_SIGMA_V2_BIN");
            std::string const path = (envPath != nullptr && *envPath != '\0') ? std::string(envPath) : std::string("models/sigma_v2_zi5.bin");
            try {
                cache.model = sketch::loadSigmaV2Model(path);
                cache.ok = true;
            } catch (std::exception const& e) {
                cache.ok = false;
                std::fprintf(stderr, "[sigma-net] FAILED to load ZI5 model '%s': %s -- served sigma is DISABLED (set SKETCH_SIGMA_V2_BIN)\n", path.c_str(), e.what());
            }
        }
    }
    if (!cache.ok) { return -1; }
    if (static_cast<int>(cache.model.NC) > cap) { return -2; }
    auto const out = sketch::sigmaV2Infer(cache.model, x, y, n);
    for (std::size_t i = 0; i < out.lnE.size(); ++i) {
        lnE_out[i] = out.lnE[i];
        sd_out[i] = out.sd[i];
    }
    return static_cast<int>(out.lnE.size());
}
