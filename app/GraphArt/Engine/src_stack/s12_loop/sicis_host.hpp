#pragma once
#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif
#include "sketch/cis_scoring.hpp"
#include "sketch/sigma_v2_native.hpp"

namespace sicis {

inline bool enabled()
{
    static int const v = [] { char const* e = std::getenv("SICIS"); return (e && e[0] == '1') ? 1 : 0; }();
    return v != 0;
}

inline bool strokeEnergies(float const* x, float const* y, int n, Sketch::CisPerCellEnergyInput& out)
{
    out = Sketch::CisPerCellEnergyInput{};
    sketch::SigmaV2Model const* model = sketch::loadSigmaV2ModelCached();
    if (!model) { return false; }
    try {
        sketch::SigmaV2CellOutput const co = sketch::sigmaV2Infer(*model, x, y, n);
        out.perCellEnergy.resize(co.lnE.size());
        for (std::size_t c = 0; c < co.lnE.size(); ++c) { out.perCellEnergy[c] = std::exp(co.lnE[c]); }
        out.memberEnergy.assign(co.memberLnE.size(), std::vector<double>(co.lnE.size()));
        for (std::size_t m = 0; m < co.memberLnE.size(); ++m) {
            for (std::size_t c = 0; c < co.lnE.size(); ++c) { out.memberEnergy[m][c] = std::exp(co.memberLnE[m][c]); }
        }
        { static const bool fr = !std::getenv("NOISE_FRAME") || std::atoi(std::getenv("NOISE_FRAME")) != 0;
          if (fr) { std::array<double, 2> m{}; double sc = 1.0; (void)sketch::sigmaV2ShapeSpaceRaw(x, y, n, static_cast<int>(model->L), m, sc);
            double const f = sc * sc; for (double& v : out.perCellEnergy) { v *= f; } for (auto& mm : out.memberEnergy) { for (double& v : mm) { v *= f; } } } }
        { static const bool lf = std::getenv("NOISE_LEVELFIX") && std::atoi(std::getenv("NOISE_LEVELFIX")) != 0;
          static constexpr double kLevelFix[7] = {0.32, 0.27, 0.21, 0.21, 0.21, 0.21, 0.21};
          static constexpr int kLevelCells[7] = {24, 12, 6, 3, 1, 1, 1};
          if (lf && out.perCellEnergy.size() == 48) { int c = 0;
            for (int l = 0; l < 7; ++l) { for (int j = 0; j < kLevelCells[l]; ++j, ++c) { out.perCellEnergy[c] *= kLevelFix[l]; for (auto& mm : out.memberEnergy) { mm[c] *= kLevelFix[l]; } } } } }
        if (char const* ns = std::getenv("NOISE_SCALE")) { double const f = std::atof(ns); if (f > 0.0 && f != 1.0) {
            for (double& v : out.perCellEnergy) { v *= f; } for (auto& m : out.memberEnergy) { for (double& v : m) { v *= f; } } } }
        out.present = true;
    } catch (std::exception const& e) {
        std::fprintf(stderr, "[sicis] ZI5 inference failed: %s\n", e.what());
        return false;
    }
    return true;
}

inline double score(std::vector<double> const& res48, Sketch::CisPerCellEnergyInput const& energies)
{
    Sketch::CisScoringConfig cfg;
    cfg.useCisScoring = true;
    cfg.exactLikelihood = 1;
    cfg.scoreDcCell = true;
    Sketch::CisScoringResult const r = Sketch::evaluateCisScoring(res48, energies.perCellEnergy, {}, cfg, energies);
    return r.ok ? r.fitness : 1.0e30;
}

struct Fast48Ctx { bool ok = false; std::size_t M = 0; double pc[8][48], lp[8][48], src[8][48], sper[48]; };
inline bool fastEnabled() { static const bool v = [] { const char* e = std::getenv("SICIS_FAST"); return !(e && e[0] == '0'); }(); return v; }
inline double scoreGeneric48(float const* res48, Sketch::CisPerCellEnergyInput const& energies)
{
    std::vector<double> r(res48, res48 + 48);
    for (double v : r) { if (!std::isfinite(v)) { return 1.0e30; } }
    return score(r, energies);
}
inline double nan48() { return std::numeric_limits<double>::quiet_NaN(); }
inline double scoreFast48(float const* res48, Sketch::CisPerCellEnergyInput const& en)
{
    thread_local Fast48Ctx C; thread_local bool have = false;
    auto const& me = en.memberEnergy;
    if (!en.present || en.perCellEnergy.size() != 48 || me.empty() || me.size() > 8 || Sketch::cisEnsembleModeFromEnv() == 0 || Sketch::cisJointCoarseFromEnv() != 0) { return nan48(); }
    bool same = have && C.M == me.size();
    if (same) { for (std::size_t m = 0; m < C.M && same; ++m) { if (me[m].size() != 48 || std::memcmp(me[m].data(), C.src[m], 48 * sizeof(double)) != 0) { same = false; } } }
    if (!same) {
        C.ok = true; C.M = me.size();
        for (std::size_t m = 0; m < C.M; ++m) {
            if (me[m].size() != 48) { C.ok = false; break; }
            for (int c = 0; c < 48; ++c) {
                C.src[m][c] = me[m][c]; double const p = std::max(me[m][c], 0.0);
                if (!(p > 0.0) || !std::isfinite(p)) { C.ok = false; }
                C.pc[m][c] = p; C.lp[m][c] = C.ok ? std::log(p) : 0.0;
            }
            if (!C.ok) { break; }
        }
        have = true;
    }
    if (!C.ok) { return nan48(); }
    double x[48];
    for (int i = 0; i < 48; ++i) { x[i] = (double)res48[i]; if (!std::isfinite(x[i])) { return 1.0e30; } }
    constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
    double e[48]; int ne = 0; double cur[48], nxt[48]; int n = 48; for (int i = 0; i < 48; ++i) { cur[i] = x[i]; }
    while (n >= 2) {
        int const pairs = n / 2;
        for (int i = 0; i < pairs; ++i) { double const a = cur[2 * i], b = cur[2 * i + 1]; nxt[i] = (a + b) * kInvSqrt2; double const d = (a - b) * kInvSqrt2; e[ne++] = d * d; }
        int na = pairs; if (n % 2 == 1) { nxt[pairs] = cur[n - 1]; na = pairs + 1; }
        for (int i = 0; i < na; ++i) { cur[i] = nxt[i]; }
        n = na;
    }
    for (int i = 0; i < n; ++i) { e[ne++] = cur[i] * cur[i]; }
    if (ne != 48) { return nan48(); }
    for (int c = 0; c < 48; ++c) { if (!std::isfinite(e[c]) || e[c] < 0.0) { return nan48(); } }
    static const double kLog2Pi = std::log(2.0 * 3.14159265358979323846);
    double shapeM[8], fitM[8];
    for (std::size_t m = 0; m < C.M; ++m) {
        double S = 0.0, D = 0.0;
        for (int c = 0; c < 48; ++c) { S += e[c] / C.pc[m][c]; D += 1.0; }
        if (!std::isfinite(S) || S == 0.0) { return nan48(); }
        double const sHat = S / D, logS = std::log(S) - std::log(D);
        double shape = 0.0;
        for (int c = 0; c < 48; ++c) { double const y = (e[c] / C.pc[m][c]) / sHat; shape += 0.5 * y + 0.5 * 1.0 * (kLog2Pi + (logS + C.lp[m][c])); }
        double const scale = 0.5 * D * (sHat - 1.0 - logS);
        if (std::isnan(shape) || std::isnan(scale)) { return nan48(); }
        shapeM[m] = shape; fitM[m] = shape + scale;
    }
    auto lse = [&](double const* v) {
        std::size_t const M = C.M; double mn = v[0]; for (std::size_t i = 0; i < M; ++i) { mn = std::min(mn, v[i]); }
        if (mn == -std::numeric_limits<double>::infinity()) { return -std::numeric_limits<double>::infinity(); }
        if (!std::isfinite(mn)) { return mn; }
        double sum = 0.0; for (std::size_t i = 0; i < M; ++i) { double const d = v[i] - mn; sum += std::isfinite(d) ? std::exp(-d) : 0.0; }
        if (!(sum > 0.0)) { return std::numeric_limits<double>::infinity(); }
        return mn - std::log(sum) + std::log(static_cast<double>(M)); };
    double const sh = lse(shapeM), tot = lse(fitM), sc = tot - sh;
    if (std::isnan(sh) || std::isnan(sc)) { return 1.0e30; }
    double const fit = sh + sc;
    return std::isfinite(fit) ? fit : 1.0e30;
}

#if defined(__aarch64__)
inline bool scoreBatch48(float const* res, int n, Sketch::CisPerCellEnergyInput const& en, double* out, bool* fb)
{
    thread_local Fast48Ctx C; thread_local bool have = false;
    auto const& me = en.memberEnergy;
    if (!en.present || en.perCellEnergy.size() != 48 || me.empty() || me.size() > 8 || Sketch::cisEnsembleModeFromEnv() == 0 || Sketch::cisJointCoarseFromEnv() != 0) { return false; }
    bool same = have && C.M == me.size();
    if (same) { for (std::size_t m = 0; m < C.M && same; ++m) { if (me[m].size() != 48 || std::memcmp(me[m].data(), C.src[m], 48 * sizeof(double)) != 0) { same = false; } } }
    if (!same) {
        C.ok = true; C.M = me.size();
        for (std::size_t m = 0; m < C.M; ++m) {
            if (me[m].size() != 48) { C.ok = false; break; }
            for (int c = 0; c < 48; ++c) {
                C.src[m][c] = me[m][c]; double const p = std::max(me[m][c], 0.0);
                if (!(p > 0.0) || !std::isfinite(p)) { C.ok = false; }
                C.pc[m][c] = p; C.lp[m][c] = C.ok ? std::log(p) : 0.0;
            }
            if (!C.ok) { break; }
        }
        have = true;
    }
    if (!C.ok) { return false; }
    static const double kLog2Pi = std::log(2.0 * 3.14159265358979323846);
    static const double kLogD = std::log(48.0);
    constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
    std::size_t const M = C.M;
    static const double lmTab[9] = {0.0, std::log(1.0), std::log(2.0), std::log(3.0), std::log(4.0), std::log(5.0), std::log(6.0), std::log(7.0), std::log(8.0)};
    double const logM = lmTab[M];
    alignas(16) double eT[48][4]; alignas(16) double q[8][48][4];
    for (int g = 0; g < n; g += 4) {
        int const nb = std::min(4, n - g);
        bool laneOk[4];
        for (int b = 0; b < 4; ++b) {
            laneOk[b] = false;
            double e[48];
            if (b < nb) {
                float const* r = res + (std::size_t)(g + b) * 48; double x[48]; bool fin = true;
                for (int i = 0; i < 48; ++i) { x[i] = (double)r[i]; if (!std::isfinite(x[i])) { fin = false; break; } }
                if (!fin) { out[g + b] = 1.0e30; fb[g + b] = false; for (int c = 0; c < 48; ++c) { eT[c][b] = 1.0; } continue; }
                int ne = 0; double cur[48], nxt[48]; int nn = 48; for (int i = 0; i < 48; ++i) { cur[i] = x[i]; }
                while (nn >= 2) {
                    int const pairs = nn / 2;
                    for (int i = 0; i < pairs; ++i) { double const a = cur[2 * i], bb = cur[2 * i + 1]; nxt[i] = (a + bb) * kInvSqrt2; double const d = (a - bb) * kInvSqrt2; e[ne++] = d * d; }
                    int na = pairs; if (nn % 2 == 1) { nxt[pairs] = cur[nn - 1]; na = pairs + 1; }
                    for (int i = 0; i < na; ++i) { cur[i] = nxt[i]; }
                    nn = na;
                }
                for (int i = 0; i < nn; ++i) { e[ne++] = cur[i] * cur[i]; }
                bool good = ne == 48; for (int c = 0; c < 48 && good; ++c) { if (!std::isfinite(e[c]) || e[c] < 0.0) { good = false; } }
                if (!good) { fb[g + b] = true; for (int c = 0; c < 48; ++c) { eT[c][b] = 1.0; } continue; }
                fb[g + b] = false; laneOk[b] = true;
                for (int c = 0; c < 48; ++c) { eT[c][b] = e[c]; }
            } else { for (int c = 0; c < 48; ++c) { eT[c][b] = 1.0; } }
        }
        float64x2_t S[8][2];
        for (std::size_t m = 0; m < M; ++m) { S[m][0] = vdupq_n_f64(0.0); S[m][1] = vdupq_n_f64(0.0); }
        for (int c = 0; c < 48; ++c) {
            float64x2_t const e0 = vld1q_f64(&eT[c][0]), e1 = vld1q_f64(&eT[c][2]);
            for (std::size_t m = 0; m < M; ++m) {
                float64x2_t const pc = vdupq_n_f64(C.pc[m][c]);
                float64x2_t const q0 = vdivq_f64(e0, pc), q1 = vdivq_f64(e1, pc);
                vst1q_f64(&q[m][c][0], q0); vst1q_f64(&q[m][c][2], q1);
                S[m][0] = vaddq_f64(S[m][0], q0); S[m][1] = vaddq_f64(S[m][1], q1);
            }
        }
        alignas(16) double Sv[8][4], sHat[8][4], logS[8][4];
        for (std::size_t m = 0; m < M; ++m) { vst1q_f64(&Sv[m][0], S[m][0]); vst1q_f64(&Sv[m][2], S[m][1]); }
        for (int b = 0; b < nb; ++b) {
            if (!laneOk[b]) { for (std::size_t m = 0; m < M; ++m) { sHat[m][b] = 1.0; logS[m][b] = 0.0; } continue; }
            for (std::size_t m = 0; m < M; ++m) {
                if (!std::isfinite(Sv[m][b]) || Sv[m][b] == 0.0) { fb[g + b] = true; laneOk[b] = false; sHat[m][b] = 1.0; logS[m][b] = 0.0; continue; }
                sHat[m][b] = Sv[m][b] / 48.0; logS[m][b] = std::log(Sv[m][b]) - kLogD;
            }
            if (!laneOk[b]) { for (std::size_t m = 0; m < M; ++m) { sHat[m][b] = 1.0; logS[m][b] = 0.0; } }
        }
        for (int b = nb; b < 4; ++b) { for (std::size_t m = 0; m < M; ++m) { sHat[m][b] = 1.0; logS[m][b] = 0.0; } }
        float64x2_t shp[8][2];
        for (std::size_t m = 0; m < M; ++m) { shp[m][0] = vdupq_n_f64(0.0); shp[m][1] = vdupq_n_f64(0.0); }
        float64x2_t const half = vdupq_n_f64(0.5), vk = vdupq_n_f64(kLog2Pi);
        for (int c = 0; c < 48; ++c) {
            for (std::size_t m = 0; m < M; ++m) {
                float64x2_t const lp = vdupq_n_f64(C.lp[m][c]);
                float64x2_t const sh0 = vld1q_f64(&sHat[m][0]), sh1 = vld1q_f64(&sHat[m][2]);
                float64x2_t const ls0 = vld1q_f64(&logS[m][0]), ls1 = vld1q_f64(&logS[m][2]);
                float64x2_t const y0 = vdivq_f64(vld1q_f64(&q[m][c][0]), sh0), y1 = vdivq_f64(vld1q_f64(&q[m][c][2]), sh1);
                float64x2_t const K0 = vaddq_f64(vk, vaddq_f64(ls0, lp)), K1 = vaddq_f64(vk, vaddq_f64(ls1, lp));
                shp[m][0] = vaddq_f64(shp[m][0], vaddq_f64(vmulq_f64(half, y0), vmulq_f64(half, K0)));
                shp[m][1] = vaddq_f64(shp[m][1], vaddq_f64(vmulq_f64(half, y1), vmulq_f64(half, K1)));
            }
        }
        alignas(16) double shpv[8][4]; for (std::size_t m = 0; m < M; ++m) { vst1q_f64(&shpv[m][0], shp[m][0]); vst1q_f64(&shpv[m][2], shp[m][1]); }
        for (int b = 0; b < nb; ++b) {
            if (!laneOk[b]) { continue; }
            double shapeM[8], fitM[8]; bool bad = false;
            for (std::size_t m = 0; m < M; ++m) {
                double const scale = 0.5 * 48.0 * (sHat[m][b] - 1.0 - logS[m][b]);
                if (std::isnan(shpv[m][b]) || std::isnan(scale)) { bad = true; break; }
                shapeM[m] = shpv[m][b]; fitM[m] = shpv[m][b] + scale;
            }
            if (bad) { fb[g + b] = true; continue; }
            auto lse = [&](double const* v) {
                double mn = v[0]; for (std::size_t i = 0; i < M; ++i) { mn = std::min(mn, v[i]); }
                if (mn == -std::numeric_limits<double>::infinity()) { return -std::numeric_limits<double>::infinity(); }
                if (!std::isfinite(mn)) { return mn; }
                double sum = 0.0; for (std::size_t i = 0; i < M; ++i) { double const d = v[i] - mn; sum += std::isfinite(d) ? (d == 0.0 ? 1.0 : std::exp(-d)) : 0.0; }
                if (!(sum > 0.0)) { return std::numeric_limits<double>::infinity(); }
                return mn - std::log(sum) + logM; };
            double const sh = lse(shapeM), tot = lse(fitM), sc = tot - sh;
            if (std::isnan(sh) || std::isnan(sc)) { out[g + b] = 1.0e30; fb[g + b] = false; continue; }
            double const fit = sh + sc;
            out[g + b] = std::isfinite(fit) ? fit : 1.0e30; fb[g + b] = false;
        }
    }
    return true;
}
#endif

inline double score48(float const* res48, Sketch::CisPerCellEnergyInput const& energies)
{
    if (fastEnabled()) { double const f = scoreFast48(res48, energies); if (!std::isnan(f)) { return f; } }
    return scoreGeneric48(res48, energies);
}

}
