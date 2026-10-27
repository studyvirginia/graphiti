#pragma once

#include "fp16_hot_path.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Sketch {

struct Fp16GramConfig {
    std::size_t rows{0};
    std::size_t columns{0};
    bool enabled{false};
    bool useFp16Storage{true};
    std::size_t threadgroupWidth{0};
};

struct Fp16GramBatch {
    std::vector<double> basis;
    std::vector<double> residual;
};

struct Fp16GramResult {
    bool ok{false};
    bool usedFallback{false};
    std::vector<double> gram;
    std::vector<double> rhs;
};

inline auto validFp16Gram(Fp16GramConfig const& cfg, Fp16GramBatch const& batch) -> bool {
    if (!cfg.enabled || cfg.rows == 0U || cfg.columns == 0U || cfg.rows > 4096U ||
        cfg.columns > 64U || batch.basis.size() != cfg.rows * cfg.columns ||
        batch.residual.size() != cfg.rows) { return false; }
    if (cfg.useFp16Storage) {
        for (double value : batch.basis) {
            std::uint16_t encoded = 0;
            if (!fp16::encodeFinite(static_cast<float>(value), encoded)) { return false; }
        }
        for (double value : batch.residual) {
            std::uint16_t encoded = 0;
            if (!fp16::encodeFinite(static_cast<float>(value), encoded)) { return false; }
        }
    } else {
        for (double value : batch.basis) { if (!std::isfinite(value)) { return false; } }
        for (double value : batch.residual) { if (!std::isfinite(value)) { return false; } }
    }
    return true;
}

inline auto accumulateFp16GramCpu(Fp16GramConfig const& cfg, Fp16GramBatch const& batch)
    -> Fp16GramResult {
    Fp16GramResult out;
    if (!validFp16Gram(cfg, batch)) { out.usedFallback = true; return out; }
    out.gram.assign(cfg.columns * cfg.columns, 0.0);
    out.rhs.assign(cfg.columns, 0.0);
    for (std::size_t r = 0; r < cfg.columns; ++r) {
        for (std::size_t c = 0; c < cfg.columns; ++c) {
            for (std::size_t i = 0; i < cfg.rows; ++i) {
                double br = 0.0, bc = 0.0;
                if (cfg.useFp16Storage) {
                    fp16::roundTrip(batch.basis[i * cfg.columns + r], br);
                    fp16::roundTrip(batch.basis[i * cfg.columns + c], bc);
                } else {
                    br = batch.basis[i * cfg.columns + r];
                    bc = batch.basis[i * cfg.columns + c];
                }
                out.gram[r * cfg.columns + c] += br * bc;
            }
        }
        for (std::size_t i = 0; i < cfg.rows; ++i) {
            double br = 0.0, rr = 0.0;
            if (cfg.useFp16Storage) {
                fp16::roundTrip(batch.basis[i * cfg.columns + r], br);
                fp16::roundTrip(batch.residual[i], rr);
            } else {
                br = batch.basis[i * cfg.columns + r];
                rr = batch.residual[i];
            }
            out.rhs[r] += br * rr;
        }
    }
    out.ok = true;
    return out;
}

class MetalFp16GramRuntime {
public:
    MetalFp16GramRuntime();
    ~MetalFp16GramRuntime();
    MetalFp16GramRuntime(MetalFp16GramRuntime const&) = delete;
    auto operator=(MetalFp16GramRuntime const&) -> MetalFp16GramRuntime& = delete;
    auto available() const -> bool;
    auto submit(Fp16GramConfig const& cfg, Fp16GramBatch const& batch) -> Fp16GramResult;
private:
    void* impl_{nullptr};
};

}
