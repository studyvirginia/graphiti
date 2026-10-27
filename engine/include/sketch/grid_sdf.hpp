#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#include <operon/operon/core/types.hpp>

namespace Sketch {

inline auto directFinalistEvalMode() -> bool& {
    static thread_local bool mode = false;
    return mode;
}

struct GridSdf {
    int resolution{128};
    Operon::Scalar x0{-0.55f};
    Operon::Scalar y0{-0.55f};
    Operon::Scalar x1{0.55f};
    Operon::Scalar y1{0.55f};
    Operon::Scalar h{0.0f};
    Operon::Scalar invH{0.0f};
    std::vector<Operon::Scalar> dist;
    std::vector<uint8_t> region;
    bool ok{false};
};

inline void buildGridSdf(std::vector<Operon::Scalar> const& dataX,
                         std::vector<Operon::Scalar> const& dataY,
                         std::size_t noiseK,
                         int resolution,
                         GridSdf& grid)
{
    grid.ok = false;
    auto const N = dataX.size();
    if (N < 2 || dataY.size() != N || resolution < 2) {
        return;
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (!std::isfinite(dataX[i]) || !std::isfinite(dataY[i])) {
            return;
        }
    }

    grid.resolution = resolution;
    grid.x0 = Operon::Scalar{-0.55f};
    grid.y0 = Operon::Scalar{-0.55f};
    grid.x1 = Operon::Scalar{0.55f};
    grid.y1 = Operon::Scalar{0.55f};
    grid.h = (grid.x1 - grid.x0) / static_cast<Operon::Scalar>(resolution - 1);
    grid.invH = static_cast<Operon::Scalar>(1.0) / grid.h;

    std::size_t const totalNodes = static_cast<std::size_t>(resolution) * static_cast<std::size_t>(resolution);
    grid.dist.resize(totalNodes);
    grid.region.resize(totalNodes);

    auto regionOf = [noiseK, N](std::size_t i) -> uint8_t {
        if (noiseK == 0) return 0;
        std::size_t const denom = (N > 1) ? (N - 1) : 1;
        return static_cast<uint8_t>(std::min(i * noiseK / denom, noiseK - 1));
    };

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    for (int r = 0; r < resolution; ++r) {
        float const yr = static_cast<float>(grid.y0) + static_cast<float>(r) * static_cast<float>(grid.h);
        float32x4_t const v_yr = vdupq_n_f32(yr);
        int c = 0;
        for (; c + 3 < resolution; c += 4) {
            float const xc0 = static_cast<float>(grid.x0) + static_cast<float>(c) * static_cast<float>(grid.h);
            float const xc1 = xc0 + static_cast<float>(grid.h);
            float const xc2 = xc1 + static_cast<float>(grid.h);
            float const xc3 = xc2 + static_cast<float>(grid.h);
            float32x4_t const v_xc = { xc0, xc1, xc2, xc3 };
            float32x4_t v_min_d2 = vdupq_n_f32(1e30f);
            uint32_t bestI[4] = {0, 0, 0, 0};

            for (std::size_t i = 0; i < N; ++i) {
                float32x4_t const v_xi = vdupq_n_f32(static_cast<float>(dataX[i]));
                float32x4_t const v_yi = vdupq_n_f32(static_cast<float>(dataY[i]));
                float32x4_t const v_dx = vsubq_f32(v_xc, v_xi);
                float32x4_t const v_dy = vsubq_f32(v_yr, v_yi);
                float32x4_t const v_d2 = vmlaq_f32(vmulq_f32(v_dx, v_dx), v_dy, v_dy);

                uint32x4_t const mask = vcltq_f32(v_d2, v_min_d2);
                uint32_t m[4];
                vst1q_u32(m, mask);
                if (m[0]) { bestI[0] = static_cast<uint32_t>(i); }
                if (m[1]) { bestI[1] = static_cast<uint32_t>(i); }
                if (m[2]) { bestI[2] = static_cast<uint32_t>(i); }
                if (m[3]) { bestI[3] = static_cast<uint32_t>(i); }
                v_min_d2 = vminq_f32(v_min_d2, v_d2);
            }

            float32x4_t const v_dist = vsqrtq_f32(v_min_d2);
            std::size_t const baseIdx = static_cast<std::size_t>(r) * resolution + c;
            vst1q_f32(reinterpret_cast<float*>(&grid.dist[baseIdx]), v_dist);
            grid.region[baseIdx]     = regionOf(bestI[0]);
            grid.region[baseIdx + 1] = regionOf(bestI[1]);
            grid.region[baseIdx + 2] = regionOf(bestI[2]);
            grid.region[baseIdx + 3] = regionOf(bestI[3]);
        }
        for (; c < resolution; ++c) {
            float const xc = static_cast<float>(grid.x0) + static_cast<float>(c) * static_cast<float>(grid.h);
            float min_d2 = 1e30f;
            std::size_t best_i = 0;
            for (std::size_t i = 0; i < N; ++i) {
                float const dx = xc - static_cast<float>(dataX[i]);
                float const dy = yr - static_cast<float>(dataY[i]);
                float const d2 = dx * dx + dy * dy;
                if (d2 < min_d2) {
                    min_d2 = d2;
                    best_i = i;
                }
            }
            std::size_t const idx = static_cast<std::size_t>(r) * resolution + c;
            grid.dist[idx] = static_cast<Operon::Scalar>(std::sqrt(min_d2));
            grid.region[idx] = regionOf(best_i);
        }
    }
#else
    for (int r = 0; r < resolution; ++r) {
        Operon::Scalar const yr = grid.y0 + static_cast<Operon::Scalar>(r) * grid.h;
        for (int c = 0; c < resolution; ++c) {
            Operon::Scalar const xc = grid.x0 + static_cast<Operon::Scalar>(c) * grid.h;
            Operon::Scalar min_d2 = std::numeric_limits<Operon::Scalar>::max();
            std::size_t best_i = 0;
            for (std::size_t i = 0; i < N; ++i) {
                Operon::Scalar const dx = xc - dataX[i];
                Operon::Scalar const dy = yr - dataY[i];
                Operon::Scalar const d2 = dx * dx + dy * dy;
                if (d2 < min_d2) {
                    min_d2 = d2;
                    best_i = i;
                }
            }
            std::size_t const idx = static_cast<std::size_t>(r) * resolution + c;
            grid.dist[idx] = std::sqrt(min_d2);
            grid.region[idx] = regionOf(best_i);
        }
    }
#endif

    grid.ok = true;
}

template <typename NoiseT>
inline void buildGridSdf(std::vector<Operon::Scalar> const& dataX,
                         std::vector<Operon::Scalar> const& dataY,
                         NoiseT const& noise,
                         int resolution,
                         GridSdf& grid)
{
    buildGridSdf(dataX, dataY, noise.size(), resolution, grid);
}

inline auto queryGridSdf(GridSdf const& grid,
                         Operon::Scalar px,
                         Operon::Scalar py,
                         uint8_t& outRegion) -> Operon::Scalar
{
    Operon::Scalar const u = (px - grid.x0) * grid.invH;
    Operon::Scalar const v = (py - grid.y0) * grid.invH;
    int const c = std::clamp(static_cast<int>(u), 0, grid.resolution - 2);
    int const r = std::clamp(static_cast<int>(v), 0, grid.resolution - 2);
    Operon::Scalar const su = std::clamp(u - static_cast<Operon::Scalar>(c), Operon::Scalar{0}, Operon::Scalar{1});
    Operon::Scalar const sv = std::clamp(v - static_cast<Operon::Scalar>(r), Operon::Scalar{0}, Operon::Scalar{1});

    int const idx = r * grid.resolution + c;
    Operon::Scalar const d00 = grid.dist[idx];
    Operon::Scalar const d10 = grid.dist[idx + 1];
    Operon::Scalar const d01 = grid.dist[idx + grid.resolution];
    Operon::Scalar const d11 = grid.dist[idx + grid.resolution + 1];

    Operon::Scalar const d0 = d00 + su * (d10 - d00);
    Operon::Scalar const d1 = d01 + su * (d11 - d01);
    Operon::Scalar const d = d0 + sv * (d1 - d0);

    int const cr = std::clamp(static_cast<int>(std::round(u)), 0, grid.resolution - 1);
    int const rr = std::clamp(static_cast<int>(std::round(v)), 0, grid.resolution - 1);
    outRegion = grid.region[rr * grid.resolution + cr];

    return d * d;
}

inline auto queryGridSdf(GridSdf const& grid,
                         Operon::Scalar px,
                         Operon::Scalar py) -> Operon::Scalar
{
    uint8_t dummy = 0;
    return queryGridSdf(grid, px, py, dummy);
}

inline void buildGridSdf(std::vector<Operon::Scalar> const& dataX,
                         std::vector<Operon::Scalar> const& dataY,
                         int resolution,
                         GridSdf& grid)
{
    buildGridSdf(dataX, dataY, std::size_t{0}, resolution, grid);
}

}
