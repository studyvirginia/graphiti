#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Sketch {

struct WeightedVarProConfig {
    bool useInverseVariance{false};
    double minWeight{0.05};
    double maxWeight{20.0};
    double ridgeFraction{1e-3};
};

struct WeightedVarProResult {
    std::vector<double> coefficients;
    std::vector<double> weights;
    double ridge{0.0};
    bool ok{false};
};

inline auto solveWeightedVarPro(std::vector<double> const& phi, std::size_t rows,
                                std::size_t cols, std::vector<double> const& target,
                                std::vector<double> const& variance = {},
                                WeightedVarProConfig cfg = {}) -> WeightedVarProResult
{
    WeightedVarProResult out;
    if (rows == 0 || cols == 0 || phi.size() != rows * cols || target.size() != rows ||
        !(cfg.minWeight > 0.0) || !(cfg.maxWeight >= cfg.minWeight) ||
        !(cfg.ridgeFraction >= 0.0) || !std::isfinite(cfg.ridgeFraction)) {
        return out;
    }

    out.weights.assign(rows, 1.0);
    if (cfg.useInverseVariance) {
        if (variance.size() != rows) { return out; }
        double sum = 0.0;
        for (std::size_t i = 0; i < rows; ++i) {
            if (!(variance[i] > 0.0) || !std::isfinite(variance[i])) { return out; }
            out.weights[i] = 1.0 / variance[i];
            sum += out.weights[i];
        }
        double const mean = sum / static_cast<double>(rows);
        if (!(mean > 0.0) || !std::isfinite(mean)) { return out; }
        for (double& w : out.weights) {
            w = std::clamp(w / mean, cfg.minWeight, cfg.maxWeight);
        }
    }

    std::vector<double> gram(cols * cols, 0.0), rhs(cols, 0.0);
    for (std::size_t i = 0; i < rows; ++i) {
        double const w = out.weights[i];
        if (!std::isfinite(target[i])) { return out; }
        for (std::size_t r = 0; r < cols; ++r) {
            double const pr = phi[i * cols + r];
            if (!std::isfinite(pr)) { return out; }
            rhs[r] += w * pr * target[i];
            for (std::size_t c = 0; c < cols; ++c) {
                double const pc = phi[i * cols + c];
                if (!std::isfinite(pc)) { return out; }
                gram[r * cols + c] += w * pr * pc;
            }
        }
    }
    double trace = 0.0;
    for (std::size_t k = 0; k < cols; ++k) { trace += gram[k * cols + k]; }
    out.ridge = cfg.ridgeFraction * trace;
    if (!std::isfinite(out.ridge)) { return out; }
    double const diagonal = (out.ridge > 0.0) ? out.ridge : std::numeric_limits<double>::epsilon();
    for (std::size_t k = 0; k < cols; ++k) { gram[k * cols + k] += diagonal; }

    std::vector<double> aug(cols * (cols + 1), 0.0);
    for (std::size_t r = 0; r < cols; ++r) {
        for (std::size_t c = 0; c < cols; ++c) { aug[r * (cols + 1) + c] = gram[r * cols + c]; }
        aug[r * (cols + 1) + cols] = rhs[r];
    }
    for (std::size_t col = 0; col < cols; ++col) {
        std::size_t pivot = col;
        double best = std::abs(aug[col * (cols + 1) + col]);
        for (std::size_t r = col + 1; r < cols; ++r) {
            double const candidate = std::abs(aug[r * (cols + 1) + col]);
            if (candidate > best) { best = candidate; pivot = r; }
        }
        if (!(best > 0.0) || !std::isfinite(best)) { return out; }
        if (pivot != col) {
            for (std::size_t c = col; c <= cols; ++c) {
                std::swap(aug[col * (cols + 1) + c], aug[pivot * (cols + 1) + c]);
            }
        }
        double const d = aug[col * (cols + 1) + col];
        for (std::size_t r = 0; r < cols; ++r) {
            if (r == col) { continue; }
            double const f = aug[r * (cols + 1) + col] / d;
            for (std::size_t c = col; c <= cols; ++c) {
                aug[r * (cols + 1) + c] -= f * aug[col * (cols + 1) + c];
            }
        }
    }
    out.coefficients.resize(cols);
    for (std::size_t r = 0; r < cols; ++r) {
        out.coefficients[r] = aug[r * (cols + 1) + cols] / aug[r * (cols + 1) + r];
        if (!std::isfinite(out.coefficients[r])) { out.coefficients.clear(); return out; }
    }
    out.ok = true;
    return out;
}

}
