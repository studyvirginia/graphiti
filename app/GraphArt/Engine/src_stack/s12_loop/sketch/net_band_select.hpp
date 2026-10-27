#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace Sketch {

struct NetBandConfig {
    double z = 1.96;
    std::vector<int> scales{1, 2, 4, 8, 16};
};

struct NetBandRow {
    int index;
    int nodes;
    double shapeErr;
    bool adequate;
    double worstRatio;
};

struct NetBandResult {
    int pickIndex;
    int accuracyIndex;
    std::vector<NetBandRow> rows;
};

namespace detail {

inline double chi2UpperFactor(double m, double z) {
    m = std::max(m, 1.0);
    double const a = 2.0 / (9.0 * m);
    double const v = std::pow(1.0 - a + z * std::sqrt(a), 3.0);
    return std::max(1.0, v);
}

inline std::vector<double> aggToK(const std::vector<double>& v, int K) {
    int const R = static_cast<int>(v.size());
    if (K >= R) return v;
    std::vector<double> out(static_cast<std::size_t>(K), 0.0);
    std::vector<int> cnt(static_cast<std::size_t>(K), 0);
    double sumAll = 0.0;
    for (int r = 0; r < R; ++r) {
        int const k = r * K / R;
        out[static_cast<std::size_t>(k)] += v[static_cast<std::size_t>(r)];
        ++cnt[static_cast<std::size_t>(k)];
        sumAll += v[static_cast<std::size_t>(r)];
    }
    for (int k = 0; k < K; ++k) {
        if (cnt[static_cast<std::size_t>(k)] > 0) out[static_cast<std::size_t>(k)] /= cnt[static_cast<std::size_t>(k)];
        else out[static_cast<std::size_t>(k)] = sumAll / R;
    }
    return out;
}

}

inline std::vector<std::vector<double>> netBands(const std::vector<double>& sigma2,
                                                 const std::vector<double>& logVar,
                                                 int nData, const NetBandConfig& cfg) {
    std::vector<std::vector<double>> out;
    out.reserve(cfg.scales.size());
    bool const haveLv = (logVar.size() == sigma2.size());
    for (int K : cfg.scales) {
        auto s2K = detail::aggToK(sigma2, K);
        std::vector<double> lvK;
        if (haveLv) lvK = detail::aggToK(logVar, K);
        int const blocks = static_cast<int>(s2K.size());
        std::vector<double> band(static_cast<std::size_t>(blocks), 0.0);
        for (int i = 0; i < blocks; ++i) {
            double epi = 1.0;
            if (haveLv) {
                double const sdLog = std::sqrt(std::exp(std::clamp(lvK[static_cast<std::size_t>(i)], -8.0, 4.0)));
                epi = std::pow(10.0, 2.0 * cfg.z * sdLog);
            }
            double const edgeLo = static_cast<double>(i) / K;
            double const edgeHi = static_cast<double>(i + 1) / K;
            double const cnt = std::max(1.0, std::nearbyint((edgeHi - edgeLo) * nData));
            double const samp = detail::chi2UpperFactor(cnt, cfg.z);
            band[static_cast<std::size_t>(i)] = s2K[static_cast<std::size_t>(i)] * epi * samp;
        }
        out.push_back(std::move(band));
    }
    return out;
}

inline NetBandResult netBandSelect(
        const std::vector<std::vector<double>>& regionErr,
        const std::vector<int>& nodes, const std::vector<double>& shapeErr,
        const std::vector<double>& sigma2, int nData,
        const std::vector<double>& logVar = {}, const NetBandConfig& cfg = {}) {
    auto const bands = netBands(sigma2, logVar, nData, cfg);

    std::vector<NetBandRow> rows;
    std::set<std::pair<int, double>> seen;
    for (int i = 0; i < static_cast<int>(regionErr.size()); ++i) {
        if (regionErr[i].empty()) continue;
        double const keyErr = std::nearbyint(shapeErr[i] * 1e6) / 1e6;
        if (!seen.insert({nodes[i], keyErr}).second) continue;
        bool adequate = true;
        double worst = 0.0;
        for (std::size_t si = 0; si < cfg.scales.size(); ++si) {
            int const K = cfg.scales[si];
            auto const reK = detail::aggToK(std::vector<double>(regionErr[i].begin(), regionErr[i].end()), K);
            auto const& up = bands[si];
            int const n = std::min(static_cast<int>(reK.size()), static_cast<int>(up.size()));
            for (int k = 0; k < n; ++k) {
                double const b = std::max(up[static_cast<std::size_t>(k)], 1e-30);
                worst = std::max(worst, reK[static_cast<std::size_t>(k)] / b);
                if (reK[static_cast<std::size_t>(k)] > up[static_cast<std::size_t>(k)]) adequate = false;
            }
        }
        rows.push_back({i, nodes[i], shapeErr[i], adequate, worst});
    }
    std::stable_sort(rows.begin(), rows.end(), [](const NetBandRow& a, const NetBandRow& b) {
        if (a.nodes != b.nodes) return a.nodes < b.nodes;
        return a.shapeErr < b.shapeErr;
    });

    int pick = -1;
    for (const auto& r : rows) if (r.adequate) { pick = r.index; break; }
    int acc = -1;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& r : rows) if (r.shapeErr < best) { best = r.shapeErr; acc = r.index; }
    return {pick, acc, rows};
}

}
