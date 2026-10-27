#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>
#include <stdexcept>
#include <string>

#include <operon/operon/core/dispatch.hpp>
#include "sketch/multitree_gp.hpp"
#include "sketch/portfolio.hpp"
#include "sketch/segment.hpp"
#include "sketch/selection_core.hpp"
#include "sketch/symmetry_detector.hpp"

namespace nb = nanobind;

NB_MODULE(sketch_native, m) {
    m.doc() = "Graphiti native search engine";

    using DTable = Operon::DispatchTable<Operon::Scalar>;

    nb::class_<Sketch::MultitreeConfig>(m, "MultitreeConfig")
        .def(nb::init<>())
        .def_rw("population_size", &Sketch::MultitreeConfig::populationSize)
        .def_rw("generations", &Sketch::MultitreeConfig::generations)
        .def_rw("max_length", &Sketch::MultitreeConfig::maxLength)
        .def_rw("max_depth", &Sketch::MultitreeConfig::maxDepth)
        .def_rw("grid_size", &Sketch::MultitreeConfig::gridSize)
        .def_rw("opt_fraction", &Sketch::MultitreeConfig::optFraction)
        .def_rw("lm_targeted", &Sketch::MultitreeConfig::lmTargeted)
        .def_rw("coeff_spread", &Sketch::MultitreeConfig::coeffSpread)
        .def_rw("seed", &Sketch::MultitreeConfig::seed)
        .def_rw("verbose", &Sketch::MultitreeConfig::verbose)
        .def_rw("floor_stop", &Sketch::MultitreeConfig::floorStop)
        .def_rw("alps_layers", &Sketch::MultitreeConfig::alpsLayers)
        .def_rw("age_gap", &Sketch::MultitreeConfig::ageGap)
        .def_rw("order_weight", &Sketch::MultitreeConfig::orderWeight)
        .def_rw("whiten_order_term", &Sketch::MultitreeConfig::whitenOrderTerm)
        .def_rw("use_lm", &Sketch::MultitreeConfig::useLm)
        .def_rw("const_optimizer", &Sketch::MultitreeConfig::constOptimizer)
        .def_rw("lm_block_nn", &Sketch::MultitreeConfig::lmBlockNn)
        .def_rw("anneal_resolution", &Sketch::MultitreeConfig::annealResolution)
        .def_rw("anneal_low_n", &Sketch::MultitreeConfig::annealLowN)
        .def_rw("anneal_stage2_resolution", &Sketch::MultitreeConfig::annealStage2Resolution)
        .def_rw("anneal_mid_n", &Sketch::MultitreeConfig::annealMidN)
        .def_rw("anneal_low_m", &Sketch::MultitreeConfig::annealLowM)
        .def_rw("anneal_mid_m", &Sketch::MultitreeConfig::annealMidM)
        .def_rw("use_lexicase", &Sketch::MultitreeConfig::useLexicase)
        .def_rw("semantic_crossover", &Sketch::MultitreeConfig::semanticCrossover)
        .def_rw("tournament_size", &Sketch::MultitreeConfig::tournamentSize)
        .def_rw("per_tree_credit", &Sketch::MultitreeConfig::perTreeCredit)
        .def_rw("floor_clamp", &Sketch::MultitreeConfig::floorClamp)
        .def_rw("floor_simplify_gens", &Sketch::MultitreeConfig::floorSimplifyGens)
        .def_rw("floor_on_chamfer", &Sketch::MultitreeConfig::floorOnChamfer)
        .def_rw("arc_length_guard_ratio", &Sketch::MultitreeConfig::arcLengthGuardRatio)
        .def_rw("closed_gap", &Sketch::MultitreeConfig::closedGap)
        .def_rw("max_fit_seconds", &Sketch::MultitreeConfig::maxFitSeconds)
        .def_rw("cross_tree_prob", &Sketch::MultitreeConfig::crossTreeProb)
        .def_rw("const_multi_start", &Sketch::MultitreeConfig::constMultiStart)
        .def_rw("const_multi_start_data_seed", &Sketch::MultitreeConfig::constMultiStartDataSeed)
        .def_rw("polish_full_objective", &Sketch::MultitreeConfig::polishFullObjective)
        .def_rw("guard_const_opt", &Sketch::MultitreeConfig::guardConstOpt)
        .def_rw("use_archive", &Sketch::MultitreeConfig::useArchive)
        .def_rw("use_complexity_hall_of_fame", &Sketch::MultitreeConfig::useComplexityHallOfFame)
        .def_rw("curvature_weight", &Sketch::MultitreeConfig::curvatureWeight)
        .def_rw("curvature_multiplicative", &Sketch::MultitreeConfig::curvatureMultiplicative)
        .def_rw("curvature_floor_turns", &Sketch::MultitreeConfig::curvatureFloorTurns)
        .def_rw("curvature_floor_rel", &Sketch::MultitreeConfig::curvatureFloorRel)
        .def_rw("curvature_harmonics", &Sketch::MultitreeConfig::curvatureHarmonics)
        .def_rw("n_threads", &Sketch::MultitreeConfig::nThreads)
        .def_rw("use_grid_sdf", &Sketch::MultitreeConfig::useGridSdf)
        .def_rw("grid_sdf_resolution", &Sketch::MultitreeConfig::gridSdfResolution)
        .def_rw("use_accelerate_transcendentals", &Sketch::MultitreeConfig::useAccelerateTranscendentals)
        .def_rw("use_bytecode_evaluator", &Sketch::MultitreeConfig::useBytecodeEvaluator)
        .def_rw("point_clamp", &Sketch::MultitreeConfig::pointClamp)
        .def_rw("in_loop_const_iter", &Sketch::MultitreeConfig::inLoopConstIter)
        .def_rw("polish_const_iter", &Sketch::MultitreeConfig::polishConstIter)
        .def_rw("use_amx_chamfer", &Sketch::MultitreeConfig::useAmxChamfer)
        .def_rw("const_magnitude_complexity", &Sketch::MultitreeConfig::constMagnitudeComplexity)
        .def_rw("champion_adequacy_band", &Sketch::MultitreeConfig::championAdequacyBand)
        .def_rw("direct_finalist_eval", &Sketch::MultitreeConfig::directFinalistEval)
        .def_rw("use_cis_fitness", &Sketch::MultitreeConfig::useCisFitness)
        .def_rw("use_si_cis_fitness", &Sketch::MultitreeConfig::useSiCisFitness)
        .def_rw("use_global_robust_nll_fitness", &Sketch::MultitreeConfig::useGlobalRobustNllFitness)
        .def_rw("use_symmetric_fitness", &Sketch::MultitreeConfig::useSymmetricFitness)
        .def_rw("robust_continuous_weight", &Sketch::MultitreeConfig::robustContinuousWeight);

    m.def("shape_err_objective", [](Sketch::MultitreeConfig const& c) -> std::string {
        std::string s;
        if (c.useCisFitness) {
            s = c.useGlobalRobustNllFitness ? "robust_t_nll" : (c.useSiCisFitness ? "si_cis_nll" : "cis_nll");
            s += "+revcov";
            if (c.robustContinuousWeight > 0.0) { s += "+robust_cont"; }
        } else {
            s = c.useSymmetricFitness ? "symmetric" : "pointloss";
        }
        if (c.orderWeight > Operon::Scalar{0}) { s += "+order"; }
        if (c.curvatureWeight > Operon::Scalar{0}) { s += "+curv"; }
        return s;
    }, nb::arg("config"),
    "Tag naming what FrontEntry.shape_err contains under this MultitreeConfig. Config-only; no computation.");

    nb::class_<Sketch::FrontEntry>(m, "FrontEntry")
        .def_ro("expr_x", &Sketch::FrontEntry::exprX)
        .def_ro("expr_y", &Sketch::FrontEntry::exprY)
        .def_ro("shape_err", &Sketch::FrontEntry::shapeErr)
        .def_ro("nodes", &Sketch::FrontEntry::nodes)
        .def_ro("region_err", &Sketch::FrontEntry::regionErr)
        .def_ro("curve_x", &Sketch::FrontEntry::curveX)
        .def_ro("curve_y", &Sketch::FrontEntry::curveY)
        .def_ro("order_err", &Sketch::FrontEntry::orderErr)
        .def_ro("err_x", &Sketch::FrontEntry::errX)
        .def_ro("err_y", &Sketch::FrontEntry::errY)
        .def_ro("turning", &Sketch::FrontEntry::turning)
        .def_ro("nml_abs_total", &Sketch::FrontEntry::nmlAbsTotal);

    nb::class_<Sketch::MultitreeResult>(m, "MultitreeResult")
        .def_ro("front", &Sketch::MultitreeResult::front)
        .def_ro("hall_of_fame", &Sketch::MultitreeResult::hallOfFame)
        .def_ro("shape_error", &Sketch::MultitreeResult::shapeError)
        .def_ro("generations", &Sketch::MultitreeResult::generations)
        .def_ro("truncated", &Sketch::MultitreeResult::truncated)
        .def_ro("n_eval_errors", &Sketch::MultitreeResult::nEvalErrors)
        .def_ro("n_const_descents", &Sketch::MultitreeResult::nConstDescents)
        .def_ro("n_const_loss_calls", &Sketch::MultitreeResult::nConstLossCalls)
        .def_ro("n_const_iters", &Sketch::MultitreeResult::nConstIters)
        .def_ro("n_nonfinite_fitness", &Sketch::MultitreeResult::nNonFiniteFitness)
        .def_ro("simplify_gens_ran", &Sketch::MultitreeResult::simplifyGensRan)
        .def_ro("anneal_switch_gen", &Sketch::MultitreeResult::annealSwitchGen)
        .def_ro("data_turning", &Sketch::MultitreeResult::dataTurning)
        .def_ro("precompute_time_ms", &Sketch::MultitreeResult::precomputeTimeMs);

    m.def("fit_stroke_multitree", [](
        nb::ndarray<float, nb::ndim<1>> data_x,
        nb::ndarray<float, nb::ndim<1>> data_y,
        Sketch::MultitreeConfig config,
        std::optional<nb::ndarray<float, nb::ndim<1>>> region_s2,
        std::optional<nb::ndarray<float, nb::ndim<1>>> region_nu
    ) -> Sketch::MultitreeResult {
        auto n = data_x.shape(0);
        if (data_y.shape(0) != n) {
            throw std::invalid_argument(
                "fit_stroke_multitree: data_x and data_y must have equal length (got "
                + std::to_string(n) + " and " + std::to_string(data_y.shape(0)) + ")");
        }
        constexpr std::int64_t kMinPoints = 8;
        if (n < kMinPoints) {
            throw std::invalid_argument(
                "fit_stroke_multitree: need at least " + std::to_string(kMinPoints)
                + " points, got " + std::to_string(n)
                + " (below this the shape distance and the 16-region residual are "
                  "interpolation, not a fit; an empty or 1-point stroke is not a stroke)");
        }
        auto const* xp = data_x.data();
        auto const* yp = data_y.data();
        for (std::int64_t i = 0; i < n; ++i) {
            if (!std::isfinite(xp[i]) || !std::isfinite(yp[i])) {
                throw std::invalid_argument(
                    "fit_stroke_multitree: data_x/data_y must be finite; found a non-finite "
                    "value at index " + std::to_string(i)
                    + " (one NaN or inf poisons the whole stroke through normalizeInPlace)");
            }
        }
        {
            auto xlo = xp[0], xhi = xp[0], ylo = yp[0], yhi = yp[0];
            for (std::int64_t i = 1; i < n; ++i) {
                xlo = std::min(xlo, xp[i]); xhi = std::max(xhi, xp[i]);
                ylo = std::min(ylo, yp[i]); yhi = std::max(yhi, yp[i]);
            }
            if (xhi == xlo && yhi == ylo) {
                throw std::invalid_argument(
                    "fit_stroke_multitree: degenerate stroke --- zero span in BOTH x and y (all "
                    + std::to_string(n) + " points identical); there is no shape to fit");
            }
        }
        std::vector<Operon::Scalar> dx(data_x.data(), data_x.data() + n);
        std::vector<Operon::Scalar> dy(data_y.data(), data_y.data() + n);
        Sketch::NoiseModel noise;
        if (region_s2.has_value() != region_nu.has_value()) {
            throw std::invalid_argument(
                "fit_stroke_multitree: region_s2 and region_nu must be supplied together "
                "(pass both for a noise model, neither for uniform MSE)");
        }
        if (region_s2.has_value() && region_nu.has_value()) {
            auto k = region_s2->shape(0);
            if (k != region_nu->shape(0)) {
                throw std::invalid_argument(
                    "fit_stroke_multitree: region_s2 and region_nu must have equal length (got "
                    + std::to_string(k) + " and " + std::to_string(region_nu->shape(0)) + ")");
            }
            if (k < 1 || k > n) {
                throw std::invalid_argument(
                    "fit_stroke_multitree: region_s2/region_nu length must be in [1, N] (got K="
                    + std::to_string(k) + ", N=" + std::to_string(n) + ")");
            }
            for (std::int64_t j = 0; j < k; ++j) {
                auto const s2v = region_s2->data()[j];
                auto const nuv = region_nu->data()[j];
                if (!std::isfinite(s2v) || s2v <= 0.0F) {
                    throw std::invalid_argument(
                        "fit_stroke_multitree: region_s2 must be finite and > 0 (it is a "
                        "squared scale and every residual is divided by it); bad value at "
                        "index " + std::to_string(j));
                }
                if (!std::isfinite(nuv)) {
                    throw std::invalid_argument(
                        "fit_stroke_multitree: region_nu must be finite (nu <= 0 selects the "
                        "Gaussian limit); bad value at index " + std::to_string(j));
                }
            }
            noise.s2.assign(region_s2->data(), region_s2->data() + k);
            noise.nu.assign(region_nu->data(), region_nu->data() + k);
        }
        if (config.gridSize < 2) {
            throw std::invalid_argument(
                "fit_stroke_multitree: grid_size must be >= 2 (got "
                + std::to_string(config.gridSize)
                + "; a 0- or 1-row grid makes normalizeInPlace deref an empty min/max range)");
        }
        if (config.populationSize < 1) {
            throw std::invalid_argument(
                "fit_stroke_multitree: population_size must be >= 1 (got "
                + std::to_string(config.populationSize)
                + "; an empty population makes the tournament index out of bounds)");
        }
        if (config.generations < 1) {
            throw std::invalid_argument(
                "fit_stroke_multitree: generations must be >= 1 (got "
                + std::to_string(config.generations) + ")");
        }
        DTable dtable;
        Sketch::MultitreeGP<DTable> gp(dtable, dx, dy, config, noise);
        nb::gil_scoped_release release;
        return gp.Run();
    },
    nb::arg("data_x"), nb::arg("data_y"),
    nb::arg("config") = Sketch::MultitreeConfig{},
    nb::arg("region_s2") = nb::none(),
    nb::arg("region_nu") = nb::none(),
    "Fit a 2D stroke: two-tree shape-matching GP, returns the Pareto front. "
    "Optional region_s2/region_nu supply the sigma-net per-region predictive "
    "density (probability-aware fitness); omit both for uniform MSE.");

    m.def("fit_portfolio", [](
        nb::ndarray<float, nb::ndim<1>> data_x,
        nb::ndarray<float, nb::ndim<1>> data_y,
        std::vector<Sketch::MultitreeConfig> members,
        int n_threads,
        std::optional<nb::ndarray<float, nb::ndim<1>>> region_s2,
        std::optional<nb::ndarray<float, nb::ndim<1>>> region_nu,
        bool cpp_halving,
        std::size_t scout_gens,
        std::size_t mid_gens,
        int n_survivors
    ) -> nb::dict {
        auto n = data_x.shape(0);
        if (data_y.shape(0) != n) {
            throw std::invalid_argument(
                "fit_portfolio: data_x and data_y must have equal length (got "
                + std::to_string(n) + " and " + std::to_string(data_y.shape(0)) + ")");
        }
        constexpr std::int64_t kMinPoints = 8;
        if (n < kMinPoints) {
            throw std::invalid_argument(
                "fit_portfolio: need at least " + std::to_string(kMinPoints)
                + " points, got " + std::to_string(n));
        }
        auto const* xp = data_x.data();
        auto const* yp = data_y.data();
        for (std::int64_t i = 0; i < n; ++i) {
            if (!std::isfinite(xp[i]) || !std::isfinite(yp[i])) {
                throw std::invalid_argument(
                    "fit_portfolio: data_x/data_y must be finite; non-finite value at index "
                    + std::to_string(i));
            }
        }
        {
            auto xlo = xp[0], xhi = xp[0], ylo = yp[0], yhi = yp[0];
            for (std::int64_t i = 1; i < n; ++i) {
                xlo = std::min(xlo, xp[i]); xhi = std::max(xhi, xp[i]);
                ylo = std::min(ylo, yp[i]); yhi = std::max(yhi, yp[i]);
            }
            if (xhi == xlo && yhi == ylo) {
                throw std::invalid_argument(
                    "fit_portfolio: degenerate stroke --- zero span in BOTH x and y");
            }
        }
        if (members.empty()) {
            throw std::invalid_argument("fit_portfolio: members must not be empty");
        }
        for (std::size_t i = 0; i < members.size(); ++i) {
            auto const& c = members[i];
            if (c.gridSize < 2 || c.populationSize < 1 || c.generations < 1) {
                throw std::invalid_argument(
                    "fit_portfolio: member " + std::to_string(i)
                    + " has grid_size < 2, population_size < 1 or generations < 1");
            }
        }
        std::vector<Operon::Scalar> dx(data_x.data(), data_x.data() + n);
        std::vector<Operon::Scalar> dy(data_y.data(), data_y.data() + n);

        Sketch::NoiseModel noise;
        if (region_s2.has_value() != region_nu.has_value()) {
            throw std::invalid_argument(
                "fit_portfolio: region_s2 and region_nu must be supplied together");
        }
        if (region_s2.has_value() && region_nu.has_value()) {
            auto k = region_s2->shape(0);
            if (k != region_nu->shape(0) || k < 1 || k > n) {
                throw std::invalid_argument(
                    "fit_portfolio: region_s2/region_nu must have equal length in [1, N]");
            }
            noise.s2.assign(region_s2->data(), region_s2->data() + k);
            noise.nu.assign(region_nu->data(), region_nu->data() + k);
        }

        std::vector<Sketch::PortfolioMemberResult> res;
        auto const t0 = std::chrono::steady_clock::now();
        {
            nb::gil_scoped_release release;
            res = Sketch::runPortfolio<DTable>(dx, dy, members, noise, n_threads,
                                               cpp_halving, scout_gens, mid_gens, n_survivors);
        }
        auto const wall = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        nb::list per_member;
        for (std::size_t i = 0; i < res.size(); ++i) {
            nb::dict d;
            d["index"] = static_cast<int>(i);
            d["seed"] = members[i].seed;
            d["population_size"] = members[i].populationSize;
            d["generations_requested"] = members[i].generations;
            d["per_tree_credit"] = members[i].perTreeCredit;
            d["ok"] = res[i].ok;
            d["error"] = res[i].error;
            d["seconds"] = res[i].seconds;
            d["front"] = nb::cast(res[i].result.front);
            d["shape_error"] = res[i].result.shapeError;
            d["generations"] = res[i].result.generations;
            per_member.append(d);
        }
        nb::dict out;
        out["members"] = per_member;
        out["seconds"] = wall;
        out["n_threads"] = n_threads;
        out["cpp_halving"] = cpp_halving;
        out["scout_gens"] = scout_gens;
        out["mid_gens"] = mid_gens;
        out["n_ok"] = static_cast<int>(
            std::count_if(res.begin(), res.end(),
                          [](Sketch::PortfolioMemberResult const& r) { return r.ok; }));
        return out;
    },
    nb::arg("data_x"), nb::arg("data_y"), nb::arg("members"),
    nb::arg("n_threads") = 1,
    nb::arg("region_s2") = nb::none(),
    nb::arg("region_nu") = nb::none(),
    nb::arg("cpp_halving") = false,
    nb::arg("scout_gens") = 10,
    nb::arg("mid_gens") = 25,
    nb::arg("n_survivors") = 0,
    "Run N parametric-engine restarts in one process on a thread pool "
    "with optional in-process multi-wave successive halving (cpp_halving). "
    "`members` is a list of MultitreeConfig, each with its own seed/pop/generations. Returns "
    "per-member fronts and timings plus the portfolio wall time; MERGING IS THE CALLER'S "
    "responsibility. Byte-identical to running the same members serially.");

    m.def("turning_estimate", [](
        nb::ndarray<float, nb::ndim<1>> data_x,
        nb::ndarray<float, nb::ndim<1>> data_y,
        bool closed,
        int K,
        int n,
        int smoother
    ) -> nb::dict {
        auto len = data_x.shape(0);
        if (data_y.shape(0) != len) {
            throw std::invalid_argument(
                "turning_estimate: data_x and data_y must have equal length");
        }
        std::vector<Operon::Scalar> dx(data_x.data(), data_x.data() + len);
        std::vector<Operon::Scalar> dy(data_y.data(), data_y.data() + len);
        nb::dict out;
        out["turning_abs"] = Sketch::turningFourier(dx, dy, closed, K, n, smoother);
        out["turning_signed"] = Sketch::turningFourierSigned(dx, dy, closed, K, n, smoother);
        std::vector<Operon::Scalar> sx, sy;
        bool const ok = Sketch::lowPassSelect(dx, dy, closed, K, n, smoother, sx, sy);
        out["ok"] = ok;
        out["gauss_width"] = Sketch::gaussianWidthForK(K, n);
        out["gauss_width_eqnoise"] = Sketch::gaussianWidthForKEqualNoise(K, n);
        out["width_used"] = (smoother == 1) ? Sketch::gaussianWidthForK(K, n)
                          : (smoother == 2) ? Sketch::gaussianWidthForKEqualNoise(K, n) : 0.0;
        out["savgol_window"] = Sketch::savgolWindowForKEqualNoise(K, n, Sketch::kSavgolOrder);
        out["savgol_order"] = Sketch::kSavgolOrder;
        out["noise_gain_fourier"] = static_cast<double>(2 * K + 1) / static_cast<double>(n);
        out["noise_gain_used"] = (smoother == 1) ? Sketch::gaussianNoiseGain(Sketch::gaussianWidthForK(K, n))
                               : (smoother == 2) ? Sketch::gaussianNoiseGain(Sketch::gaussianWidthForKEqualNoise(K, n))
                               : (smoother == 3) ? Sketch::savgolNoiseGain(
                                     Sketch::savgolWindowForKEqualNoise(K, n, Sketch::kSavgolOrder),
                                     Sketch::kSavgolOrder)
                               : static_cast<double>(2 * K + 1) / static_cast<double>(n);
        if (ok) {
            nb::list lx, ly;
            for (std::size_t i = 0; i < sx.size(); ++i) {
                lx.append(static_cast<double>(sx[i]));
                ly.append(static_cast<double>(sy[i]));
            }
            out["smooth_x"] = lx;
            out["smooth_y"] = ly;
        }
        return out;
    },
    nb::arg("data_x"), nb::arg("data_y"), nb::arg("closed") = false, nb::arg("K") = 8,
    nb::arg("n") = 240, nb::arg("smoother") = 0,
    "Turning number (rectified and signed, in TURNS) of a stroke under the engine's low-pass, plus "
    "the smoothed polyline itself. Pure: fits nothing.");

    m.def("detect_corners", [](
        nb::ndarray<float, nb::ndim<1>> data_x,
        nb::ndarray<float, nb::ndim<1>> data_y,
        bool closed,
        double sigma,
        int null_mode,
        int smoother
    ) -> nb::list {
        auto n = data_x.shape(0);
        if (data_y.shape(0) != n) {
            throw std::invalid_argument(
                "detect_corners: data_x and data_y must have equal length");
        }
        std::vector<Operon::Scalar> dx(data_x.data(), data_x.data() + n);
        std::vector<Operon::Scalar> dy(data_y.data(), data_y.data() + n);
        auto const corners = Sketch::detectCorners(dx, dy, closed, sigma, null_mode, smoother);
        nb::list out;
        for (auto const& c : corners) {
            nb::dict d;
            d["index"] = static_cast<int>(c.index);
            d["score"] = c.score;
            d["angle"] = c.angle;
            d["scales"] = c.scales;
            d["finest_k"] = c.finestK;
            out.append(d);
        }
        return out;
    },
    nb::arg("data_x"), nb::arg("data_y"), nb::arg("closed") = false, nb::arg("sigma") = 0.02,
    nb::arg("null_mode") = 0, nb::arg("smoother") = 0,
    "Corner cuts for segmentation: a sigma-aware turning spike that persists "
    "across >= 2 adjacent harmonic scales. Returns a list of dicts with index/score/angle/scales/"
    "finest_k. Pure: fits nothing. null_mode 0 (default) = the constant-curvature-arc null; "
    "1 / 2 = a smooth null (the stroke's own K/2 low-pass) with a family-wise / position-wise band.");

    nb::class_<Sketch::RotationCandidate>(m, "RotationCandidate")
        .def_ro("order", &Sketch::RotationCandidate::order)
        .def_ro("residual", &Sketch::RotationCandidate::residual);

    nb::class_<Sketch::ReflectionAxis>(m, "ReflectionAxis")
        .def_ro("angle", &Sketch::ReflectionAxis::angle)
        .def_ro("residual", &Sketch::ReflectionAxis::residual);

    nb::class_<Sketch::SymmetryReport>(m, "SymmetryReport")
        .def_ro("cx", &Sketch::SymmetryReport::cx)
        .def_ro("cy", &Sketch::SymmetryReport::cy)
        .def_ro("mean_spacing", &Sketch::SymmetryReport::meanSpacing)
        .def_ro("rotational_order", &Sketch::SymmetryReport::rotationalOrder)
        .def_ro("rotational_residual", &Sketch::SymmetryReport::rotationalResidual)
        .def_ro("rotation_candidates", &Sketch::SymmetryReport::rotationCandidates)
        .def_ro("reflection_axes", &Sketch::SymmetryReport::reflectionAxes);

    nb::class_<Sketch::SymmetryDetectOptions>(m, "SymmetryDetectOptions")
        .def(nb::init<>())
        .def_rw("max_order", &Sketch::SymmetryDetectOptions::maxOrder)
        .def_rw("tol_gaps", &Sketch::SymmetryDetectOptions::tolGaps)
        .def_rw("axis_steps", &Sketch::SymmetryDetectOptions::axisSteps)
        .def_rw("axis_polish_steps", &Sketch::SymmetryDetectOptions::axisPolishSteps)
        .def_rw("axis_dedup_rad", &Sketch::SymmetryDetectOptions::axisDedupRad);

    m.def("detect_symmetry", [](
        nb::ndarray<float, nb::ndim<1>> data_x,
        nb::ndarray<float, nb::ndim<1>> data_y,
        Sketch::SymmetryDetectOptions opt
    ) -> Sketch::SymmetryReport {
        auto n = data_x.shape(0);
        if (data_y.shape(0) != n) {
            throw std::invalid_argument("detect_symmetry: data_x and data_y must have the same length");
        }
        std::vector<Sketch::SymPoint> pts;
        pts.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            pts.push_back({static_cast<double>(data_x.data()[i]), static_cast<double>(data_y.data()[i])});
        }
        return Sketch::detectSymmetry(pts, opt);
    },
    nb::arg("data_x"), nb::arg("data_y"), nb::arg("options") = Sketch::SymmetryDetectOptions{},
    "Report-only rotational/reflection symmetry detection on raw stroke points. "
    "Does not affect fit_stroke_multitree.");

    m.def("make_front_entry", [](std::string expr_x, std::string expr_y, float shape_err, int nodes,
                                 std::vector<float> curve_x, std::vector<float> curve_y,
                                 std::vector<float> region_err) -> Sketch::FrontEntry {
        Sketch::FrontEntry e;
        e.exprX = std::move(expr_x);
        e.exprY = std::move(expr_y);
        e.shapeErr = static_cast<Operon::Scalar>(shape_err);
        e.nodes = nodes;
        e.curveX.assign(curve_x.begin(), curve_x.end());
        e.curveY.assign(curve_y.begin(), curve_y.end());
        e.regionErr.assign(region_err.begin(), region_err.end());
        return e;
    },
    nb::arg("expr_x"), nb::arg("expr_y"), nb::arg("shape_err"), nb::arg("nodes"),
    nb::arg("curve_x") = std::vector<float>{}, nb::arg("curve_y") = std::vector<float>{},
    nb::arg("region_err") = std::vector<float>{},
    "Build a FrontEntry (tests / tooling that feed select_pool a synthetic pool).");

    m.def("resample_normalize", [](nb::ndarray<float, nb::ndim<1>> x, nb::ndarray<float, nb::ndim<1>> y,
                                   int n) -> nb::tuple {
        if (y.shape(0) != x.shape(0)) {
            throw std::invalid_argument("resample_normalize: x and y must have equal length");
        }
        std::vector<Operon::Scalar> ox, oy;
        if (!Sketch::resampleAndNormalizeStroke(x.data(), y.data(), static_cast<int>(x.shape(0)), n, ox, oy)) {
            throw std::invalid_argument("resample_normalize: invalid stroke (the bridge returns -1)");
        }
        return nb::make_tuple(ox, oy);
    },
    nb::arg("x"), nb::arg("y"), nb::arg("n") = 120,
    "The product bridge's arc-length resample + bbox-normalize (the frame select_pool expects).");

    m.def("select_pool", [](
        std::vector<float> data_x, std::vector<float> data_y,
        std::vector<Sketch::FrontEntry> pooled, std::vector<int> sources,
        std::vector<Sketch::MultitreeConfig> members,
        std::optional<std::vector<float>> point_s2, std::optional<std::vector<float>> point_nu,
        int ruler_mode, bool noise_floor_select, double adequacy_band, bool member_log
    ) -> nb::dict {
        if (data_x.size() != data_y.size()) {
            throw std::invalid_argument("select_pool: data_x and data_y must have equal length");
        }
        if (sources.size() != pooled.size()) {
            throw std::invalid_argument("select_pool: sources must have one entry per pooled member");
        }
        Sketch::NoiseModel noise;
        if (point_s2.has_value()) {
            noise.s2.assign(point_s2->begin(), point_s2->end());
            if (point_nu.has_value()) {
                if (point_nu->size() != point_s2->size()) {
                    throw std::invalid_argument("select_pool: point_nu must match point_s2 in length");
                }
                noise.nu.assign(point_nu->begin(), point_nu->end());
            } else {
                noise.nu.assign(noise.s2.size(), Operon::Scalar{0});
            }
        }
        std::vector<Operon::Scalar> dx(data_x.begin(), data_x.end());
        std::vector<Operon::Scalar> dy(data_y.begin(), data_y.end());
        Sketch::SelectionOptions opt;
        opt.rulerMode = ruler_mode;
        opt.noiseFloorSelect = noise_floor_select;
        opt.adequacyBand = adequacy_band;
        opt.memberLog = member_log;
        DTable dtable;
        Sketch::SelectionResult r;
        {
            nb::gil_scoped_release release;
            r = Sketch::selectFromPool<DTable>(dtable, std::move(pooled), sources, members, dx, dy, noise, opt);
        }
        nb::dict out;
        out["status"] = r.status;
        out["champion_index"] = static_cast<int>(r.championIndex);
        out["champion_rule"] = r.championRule;
        out["n_pooled"] = r.nPooled;
        out["n_deduped"] = r.nDeduped;
        out["scorer_built"] = r.scorerBuilt;
        out["champion_inserted_off_skyline"] = r.championInsertedOffSkyline;
        out["n_nml_abs_ok"] = r.nNmlAbsOk;
        out["n_nml_abs_missing"] = r.nNmlAbsMissing;
        std::vector<float> se, r2, dl;
        for (auto const& st : r.stats) { se.push_back(st.shapeErr); r2.push_back(st.r2); dl.push_back(st.dlBits); }
        out["shape_err"] = se;
        out["r2"] = r2;
        out["dl_bits"] = dl;
        out["skyline"] = nb::cast(std::move(r.skyline));
        return out;
    },
    nb::arg("data_x"), nb::arg("data_y"), nb::arg("pooled"), nb::arg("sources"), nb::arg("members"),
    nb::arg("point_s2") = nb::none(), nb::arg("point_nu") = nb::none(),
    nb::arg("ruler_mode") = 2, nb::arg("noise_floor_select") = true,
    nb::arg("adequacy_band") = 0.0, nb::arg("member_log") = false,
    "The app's post-search selection (Sketch::selectFromPool). data_x/y = "
    "the stroke in the fit frame (resample_normalize); pooled = FrontEntry list; sources = member "
    "index per entry (-1 = C2 route / non-two-tree); members = the portfolio MultitreeConfigs "
    "(member 0 builds the two-tree scorer). Returns status (0 ok, -2 empty), champion_index into "
    "skyline, champion_rule (0 accuracy end, 1 adequacy band, 2 net band), skyline, and the emitted "
    "shape_err/r2/dl_bits (dl_bits is display-only).");
}
