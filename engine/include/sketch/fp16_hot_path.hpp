#pragma once


#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace Sketch::fp16 {

constexpr float kMaxFinite = 65504.0F;

inline auto bits(float value) -> std::uint32_t {
    std::uint32_t out = 0;
    std::memcpy(&out, &value, sizeof(out));
    return out;
}

inline auto fromBits(std::uint32_t value) -> float {
    float out = 0.0F;
    std::memcpy(&out, &value, sizeof(out));
    return out;
}

inline auto encodeFinite(float value, std::uint16_t& out) -> bool {
    if (!std::isfinite(value) || std::fabs(value) > kMaxFinite) { return false; }
    std::uint32_t const u = bits(value);
    std::uint32_t const sign = (u >> 16U) & 0x8000U;
    std::uint32_t const exponent = (u >> 23U) & 0xffU;
    std::uint32_t mantissa = u & 0x7fffffU;
    if (exponent == 0U) {
        out = static_cast<std::uint16_t>(sign);
        return true;
    }
    int const halfExp = static_cast<int>(exponent) - 127 + 15;
    if (halfExp <= 0) {
        int shift = 14 - halfExp;
        if (shift >= 24) { out = static_cast<std::uint16_t>(sign); return true; }
        mantissa |= 0x800000U;
        std::uint32_t halfMant = mantissa >> static_cast<unsigned>(shift);
        std::uint32_t const rem = mantissa & ((1U << static_cast<unsigned>(shift)) - 1U);
        if (rem > (1U << static_cast<unsigned>(shift - 1)) ||
            (rem == (1U << static_cast<unsigned>(shift - 1)) && (halfMant & 1U) != 0U)) {
            ++halfMant;
        }
        out = static_cast<std::uint16_t>(sign | halfMant);
        return true;
    }
    if (halfExp >= 31) { return false; }
    std::uint32_t halfMant = mantissa >> 13U;
    std::uint32_t const rem = mantissa & 0x1fffU;
    if (rem > 0x1000U || (rem == 0x1000U && (halfMant & 1U) != 0U)) { ++halfMant; }
    if (halfMant == 0x400U) {
        halfMant = 0U;
        if (halfExp + 1 >= 31) { return false; }
        out = static_cast<std::uint16_t>(sign |
            (static_cast<std::uint32_t>(halfExp + 1) << 10U));
        return true;
    }
    out = static_cast<std::uint16_t>(sign |
        (static_cast<std::uint32_t>(halfExp) << 10U) | halfMant);
    return true;
}

inline auto decode(std::uint16_t value) -> float {
    std::uint32_t const sign = (static_cast<std::uint32_t>(value & 0x8000U)) << 16U;
    std::uint32_t const exponent = (value >> 10U) & 0x1fU;
    std::uint32_t const mantissa = value & 0x3ffU;
    if (exponent == 0U) {
        if (mantissa == 0U) { return fromBits(sign); }
        float result = std::ldexp(static_cast<float>(mantissa), -24);
        return sign != 0U ? -result : result;
    }
    if (exponent == 31U) { return std::numeric_limits<float>::quiet_NaN(); }
    std::uint32_t const fb = sign |
        (static_cast<std::uint32_t>(exponent - 15 + 127) << 23U) | (mantissa << 13U);
    return fromBits(fb);
}

inline auto roundTrip(double value, double& restored) -> bool {
    std::uint16_t encoded = 0;
    if (!encodeFinite(static_cast<float>(value), encoded)) { return false; }
    restored = static_cast<double>(decode(encoded));
    return std::isfinite(restored);
}

inline auto winnerIndex(std::vector<double> const& scores) -> std::size_t {
    if (scores.empty()) { return scores.size(); }
    return static_cast<std::size_t>(std::distance(scores.begin(),
        std::min_element(scores.begin(), scores.end())));
}

inline auto winnerStable(std::vector<double> const& scores) -> bool {
    std::vector<double> rounded;
    rounded.reserve(scores.size());
    for (double score : scores) {
        double value = 0.0;
        if (!roundTrip(score, value)) { return false; }
        rounded.push_back(value);
    }
    std::size_t const originalWinner = winnerIndex(scores);
    std::size_t const roundedWinner = winnerIndex(rounded);
    if (originalWinner != roundedWinner) { return false; }
    for (std::size_t i = 0; i < rounded.size(); ++i) {
        if (i != originalWinner && rounded[i] == rounded[originalWinner] &&
            scores[i] != scores[originalWinner]) {
            return false;
        }
    }
    return true;
}

}
