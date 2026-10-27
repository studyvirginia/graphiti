#pragma once

#include "sketch/multitree_gp.hpp"
#include "sketch/net_band_select.hpp"
#include "sketch/one_ruler.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef SKETCH_U1B
#define SKETCH_U1B 1
#endif

namespace Sketch {

struct SelectionOptions {
    int rulerMode{2};
    bool noiseFloorSelect{true};
    double adequacyBand{0.0};
    bool memberLog{false};
};

struct SelectionCandidateStats {
    float shapeErr{0};
    float r2{0};
    float dlBits{0};
};

struct SelectionResult {
    int status{0};
    std::vector<FrontEntry> skyline;
    std::vector<SelectionCandidateStats> stats;
    std::size_t championIndex{0};
    int championRule{0};
    std::size_t nPooled{0};
    std::size_t nDeduped{0};
    bool scorerBuilt{false};
    bool championInsertedOffSkyline{false};
    std::size_t nNmlAbsOk{0};
    std::size_t nNmlAbsMissing{0};
};

inline bool resampleAndNormalizeStroke(
    const float* in_x, const float* in_y, int n_in, int target_n,
    std::vector<Operon::Scalar>& out_x, std::vector<Operon::Scalar>& out_y)
{
    if (n_in < 4 || target_n < 4) return false;

    for (int i = 0; i < n_in; ++i) {
        if (!std::isfinite(in_x[i]) || !std::isfinite(in_y[i])) return false;
    }

    std::vector<double> cum(n_in, 0.0);
    for (int i = 1; i < n_in; ++i) {
        double dx = in_x[i] - in_x[i - 1];
        double dy = in_y[i] - in_y[i - 1];
        cum[i] = cum[i - 1] + std::hypot(dx, dy);
    }
    double total_len = cum.back();
    if (total_len <= 1e-4) return false;

    std::vector<double> rx(target_n), ry(target_n);
    int seg = 0;
    for (int i = 0; i < target_n; ++i) {
        double target_d = (static_cast<double>(i) / (target_n - 1)) * total_len;
        while (seg < n_in - 2 && cum[seg + 1] < target_d) {
            seg++;
        }
        double d0 = cum[seg];
        double d1 = cum[seg + 1];
        double t = (d1 > d0) ? (target_d - d0) / (d1 - d0) : 0.0;
        rx[i] = in_x[seg] + t * (in_x[seg + 1] - in_x[seg]);
        ry[i] = in_y[seg] + t * (in_y[seg + 1] - in_y[seg]);
    }

    double min_x = *std::min_element(rx.begin(), rx.end());
    double max_x = *std::max_element(rx.begin(), rx.end());
    double min_y = *std::min_element(ry.begin(), ry.end());
    double max_y = *std::max_element(ry.begin(), ry.end());

    double cx = (min_x + max_x) * 0.5;
    double cy = (min_y + max_y) * 0.5;
    double span = std::max({max_x - min_x, max_y - min_y, 1e-4});

    out_x.resize(target_n);
    out_y.resize(target_n);
    for (int i = 0; i < target_n; ++i) {
        out_x[i] = static_cast<Operon::Scalar>((rx[i] - cx) / span);
        out_y[i] = static_cast<Operon::Scalar>((ry[i] - cy) / span);
    }
    return true;
}

inline float candidateR2(std::vector<Operon::Scalar> const& dx, std::vector<Operon::Scalar> const& dy,
                         std::vector<Operon::Scalar> const& cvx, std::vector<Operon::Scalar> const& cvy)
{
    int const N = static_cast<int>(dx.size());
    int const M = static_cast<int>(std::min(cvx.size(), cvy.size()));
    if (N == 0 || M == 0) { return 0.0f; }
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < N; ++i) { mx += dx[static_cast<std::size_t>(i)]; my += dy[static_cast<std::size_t>(i)]; }
    mx /= N; my /= N;
    double ssTot = 0.0, ssRes = 0.0;
    for (int i = 0; i < N; ++i) {
        double const px = dx[static_cast<std::size_t>(i)], py = dy[static_cast<std::size_t>(i)];
        ssTot += (px - mx) * (px - mx) + (py - my) * (py - my);
        double best = std::numeric_limits<double>::infinity();
        for (int j = 0; j < M; ++j) {
            double const ex = px - cvx[static_cast<std::size_t>(j)], ey = py - cvy[static_cast<std::size_t>(j)];
            double const d = ex * ex + ey * ey;
            if (d < best) { best = d; }
        }
        ssRes += best;
    }
    if (!(ssTot > 1e-12)) { return 0.0f; }
    double const r2 = 1.0 - ssRes / ssTot;
    if (!std::isfinite(r2)) { return 0.0f; }
    return static_cast<float>(std::max(0.0, std::min(1.0, r2)));
}

inline auto candidateStats(FrontEntry const& e, std::vector<Operon::Scalar> const& dx,
                           std::vector<Operon::Scalar> const& dy) -> SelectionCandidateStats
{
    SelectionCandidateStats s;
    s.shapeErr = static_cast<float>(e.shapeErr);
    s.r2 = candidateR2(dx, dy, e.curveX, e.curveY);
    if (std::isfinite(e.nmlAbsTotal)) {
        constexpr double kLn2 = 0.6931471805599453;
        s.dlBits = static_cast<float>(e.nmlAbsTotal / kLn2);
    } else {
        float complexity_bits = static_cast<float>(e.nodes) * 3.58f;
        float error_bits = std::max(0.0f, -std::log2(std::max(1e-5f, s.shapeErr)) * 10.0f);
        s.dlBits = complexity_bits + error_bits;
    }
    return s;
}

inline auto selectionErrBucket(FrontEntry const& e) -> long long
{
    double const q = static_cast<double>(e.shapeErr) / 1e-7;
    if (!(q < 9.2e18)) { return std::numeric_limits<long long>::max(); }
    if (q < -9.2e18) { return std::numeric_limits<long long>::min(); }
    return std::llround(q);
}

template <typename DTable>
auto selectFromPool(DTable const& dtable,
                    std::vector<FrontEntry> pooled,
                    std::vector<int> const& pooledSrc,
                    std::vector<MultitreeConfig> const& members,
                    std::vector<Operon::Scalar> const& dx,
                    std::vector<Operon::Scalar> const& dy,
                    NoiseModel const& nativeNoise,
                    SelectionOptions const& opt) -> SelectionResult
{
    SelectionResult out;
    if (pooled.empty()) { out.status = -2; return out; }

    int const rulerMode = opt.rulerMode;
    bool const oneRuler = rulerMode == 1;
    bool const twoTreeRuler = rulerMode == 2;
    bool const k1Log = opt.memberLog;
    if (oneRuler || twoTreeRuler || k1Log) {
        std::vector<Operon::Scalar> const& s2 = nativeNoise.s2;
        bool anyNonTwoTree = false;
        for (std::size_t i = 0; i < pooled.size(); ++i) {
            int const srcI = i < pooledSrc.size() ? pooledSrc[i] : -2;
            if (srcI < 0) { anyNonTwoTree = true; break; }
        }
        std::unique_ptr<MultitreeGP<DTable>> ttScorer;
        if (((twoTreeRuler && anyNonTwoTree) || k1Log) && !members.empty()) {
            auto scorerCfg = members[0];
            scorerCfg.annealResolution = 0.0;
            ttScorer = std::make_unique<MultitreeGP<DTable>>(dtable, dx, dy, scorerCfg, nativeNoise);
            out.scorerBuilt = true;
        }
        for (std::size_t i = 0; i < pooled.size(); ++i) {
            auto& fe = pooled[i];
            int const srcI = i < pooledSrc.size() ? pooledSrc[i] : -2;
            double const ruler = (oneRuler || k1Log)
                ? oneRulerSiCisShapeNll(fe.curveX, fe.curveY, dx, dy, s2)
                : std::numeric_limits<double>::quiet_NaN();
            double const tt = (ttScorer != nullptr && (k1Log || srcI < 0) && fe.curveX.size() >= 2)
                ? static_cast<double>(ttScorer->frontObjectiveFromGrid(fe.curveX, fe.curveY))
                : std::numeric_limits<double>::quiet_NaN();
            if (k1Log) {
                int const src = srcI;
                char const* tag = "?";
                if (src == -1) {
                    tag = nativeNoise.empty() ? "route:pointloss+order+extra" : "route:whitened_pointloss+order+extra";
                } else if (src >= 0 && static_cast<std::size_t>(src) < members.size()) {
                    auto const& mc = members[static_cast<std::size_t>(src)];
                    tag = mc.useCisFitness ? (mc.useSiCisFitness ? "mt:si_cis_nll+revcov+order" : "mt:cis_nll+revcov+order")
                                           : "mt:pointloss+order";
                }
                std::fprintf(stderr, "K1M src=%d tag=%s nodes=%d err=%.9g ruler=%.9g tt=%.9g ruler_on=%d x=%s y=%s\n",
                             src, tag, fe.nodes, static_cast<double>(fe.shapeErr), ruler, tt, rulerMode,
                             fe.exprX.c_str(), fe.exprY.c_str());
            }
            if (twoTreeRuler && srcI < 0) {
                fe.shapeErr = std::isnan(tt) ? std::numeric_limits<Operon::Scalar>::max()
                                             : static_cast<Operon::Scalar>(tt);
            }
            if (oneRuler) {
                fe.shapeErr = std::isfinite(ruler) ? static_cast<Operon::Scalar>(ruler)
                                                   : std::numeric_limits<Operon::Scalar>::quiet_NaN();
            }
        }
        if (oneRuler) {
            pooled.erase(std::remove_if(pooled.begin(), pooled.end(),
                                        [](FrontEntry const& e) { return !std::isfinite(e.shapeErr); }),
                         pooled.end());
            if (pooled.empty()) { out.status = -2; return out; }
        }
    }
    out.nPooled = pooled.size();

    auto const activeErr = [](const auto& e) { return e.shapeErr; };

    std::vector<FrontEntry> deduped;
    {
        std::unordered_map<std::string, std::size_t> pos;
        for (auto& fe : pooled) {
            std::string key = fe.exprX + "||" + fe.exprY;
            auto const it = pos.find(key);
            if (it == pos.end()) {
                pos.emplace(std::move(key), deduped.size());
                deduped.push_back(std::move(fe));
            } else {
                auto& kept = deduped[it->second];
                bool const better = activeErr(fe) < activeErr(kept)
                    || (activeErr(fe) == activeErr(kept) && fe.nodes < kept.nodes);
                if (better) { kept = std::move(fe); }
            }
        }
    }
    out.nDeduped = deduped.size();

    std::sort(deduped.begin(), deduped.end(), [](const auto& a, const auto& b) {
        auto const aq = selectionErrBucket(a), bq = selectionErrBucket(b);
        if (aq != bq) { return aq < bq; }
        return a.nodes < b.nodes;
    });

#if SKETCH_U1B
    deduped.erase(std::remove_if(deduped.begin(), deduped.end(),
                                 [](FrontEntry const& e) { return !std::isfinite(e.shapeErr); }),
                  deduped.end());
    if (deduped.empty()) { out.status = -2; return out; }
#endif

    bool const t13Champ = Sketch::sketchT13On(Sketch::sketchT13Mask(), 0);
    std::size_t t13ArgminIdx = std::numeric_limits<std::size_t>::max();
    if (t13Champ) {
        double bestAbs = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < deduped.size(); ++i) {
            double const v = deduped[i].nmlAbsTotal;
            if (std::isfinite(v)) {
                ++out.nNmlAbsOk;
                if (v < bestAbs) { bestAbs = v; t13ArgminIdx = i; }
            } else {
                ++out.nNmlAbsMissing;
            }
        }
    }

    std::vector<FrontEntry>& skyline = out.skyline;
    int best_nodes = std::numeric_limits<int>::max();
    std::size_t t13ChampSkylineIdx = 0;
    for (std::size_t i = 0; i < deduped.size(); ++i) {
        auto& fe = deduped[i];
        int const nodes = fe.nodes;
        bool const isNewMin = nodes < best_nodes;
        bool const isT13Champ = t13Champ && (i == t13ArgminIdx);
        if (isNewMin || isT13Champ) {
            skyline.push_back(std::move(fe));
            if (isNewMin) { best_nodes = nodes; }
            if (isT13Champ) {
                t13ChampSkylineIdx = skyline.size() - 1;
                if (!isNewMin) { out.championInsertedOffSkyline = true; }
            }
        }
    }

    out.stats.reserve(skyline.size());
    for (auto const& e : skyline) { out.stats.push_back(candidateStats(e, dx, dy)); }

    std::size_t championIndex = 0;
    double const adequacyBand = opt.adequacyBand;
    if (adequacyBand > 0.0 && !skyline.empty()) {
        championIndex = 0;
        auto const errOf = [](const FrontEntry& e) -> double { return e.shapeErr; };
        double const bestErr = errOf(skyline.front());
        double const thresh = bestErr * (1.0 + adequacyBand);
        for (std::size_t i = 1; i < skyline.size(); ++i) {
            double const ei = errOf(skyline[i]);
            bool const adequate = (std::abs(ei - bestErr) <= 1e-7) || (ei <= thresh);
            if (!adequate) { break; }
            if (skyline[i].nodes < skyline[championIndex].nodes) { championIndex = i; }
        }
        if (championIndex != 0) { out.championRule = 1; }
    }
    if (opt.noiseFloorSelect && !nativeNoise.s2.empty() && !skyline.empty()
        && !skyline.front().regionErr.empty()) {
        int const R = static_cast<int>(skyline.front().regionErr.size());
        int const P = static_cast<int>(nativeNoise.s2.size());
        std::vector<double> sigma2(static_cast<std::size_t>(R), 0.0);
        std::vector<int> scnt(static_cast<std::size_t>(R), 0);
        for (int i = 0; i < P; ++i) {
            int k = static_cast<int>(static_cast<long long>(i) * R / P);
            if (k > R - 1) { k = R - 1; }
            sigma2[static_cast<std::size_t>(k)] += static_cast<double>(nativeNoise.s2[static_cast<std::size_t>(i)]);
            ++scnt[static_cast<std::size_t>(k)];
        }
        for (int k = 0; k < R; ++k) {
            if (scnt[static_cast<std::size_t>(k)] > 0) {
                sigma2[static_cast<std::size_t>(k)] /= scnt[static_cast<std::size_t>(k)];
            }
        }
        std::vector<double> const logVarR;
        std::vector<std::vector<double>> regionErr(skyline.size());
        std::vector<int> nds(skyline.size());
        std::vector<double> serr(skyline.size());
        bool anyRegion = false;
        for (std::size_t m = 0; m < skyline.size(); ++m) {
            auto const& re = skyline[m].regionErr;
            if (static_cast<int>(re.size()) == R) {
                regionErr[m].assign(re.begin(), re.end());
                anyRegion = true;
            }
            nds[m] = skyline[m].nodes;
            serr[m] = static_cast<double>(skyline[m].shapeErr);
        }
        std::vector<double> s2sorted = sigma2;
        std::sort(s2sorted.begin(), s2sorted.end());
        double const s2Median = s2sorted.empty() ? 1e9 : s2sorted[s2sorted.size() / 2];
        double s2Ceiling = 0.1;
        if (sf2On(sf2Mask(), 13) && dx.size() >= 2 && dx.size() == dy.size()) {
            double mx = 0.0, my = 0.0;
            for (std::size_t i = 0; i < dx.size(); ++i) { mx += dx[i]; my += dy[i]; }
            mx /= static_cast<double>(dx.size()); my /= static_cast<double>(dx.size());
            double cxx = 0.0, cyy = 0.0, cxy = 0.0;
            for (std::size_t i = 0; i < dx.size(); ++i) {
                double const a = dx[i] - mx, b = dy[i] - my;
                cxx += a * a; cyy += b * b; cxy += a * b;
            }
            cxx /= static_cast<double>(dx.size()); cyy /= static_cast<double>(dx.size()); cxy /= static_cast<double>(dx.size());
            s2Ceiling = 0.5 * (cxx + cyy) + std::sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);
        }
        bool const s2Plausible = std::isfinite(s2Median) && s2Median < s2Ceiling;
        if (anyRegion && s2Plausible) {
            auto const res = netBandSelect(regionErr, nds, serr, sigma2, P, logVarR);
            if (res.pickIndex >= 0) {
                championIndex = static_cast<std::size_t>(res.pickIndex);
                out.championRule = 2;
            }
        }
    }
    if (t13Champ && t13ArgminIdx != std::numeric_limits<std::size_t>::max()) {
        championIndex = t13ChampSkylineIdx;
        out.championRule = 3;
        double const champBits = out.stats[championIndex].dlBits;
        for (std::size_t i = 0; i < skyline.size(); ++i) {
            if (std::isfinite(skyline[i].nmlAbsTotal)) {
                out.stats[i].dlBits = static_cast<float>(out.stats[i].dlBits - champBits);
            }
        }
    }
    out.championIndex = championIndex;
    return out;
}

}
