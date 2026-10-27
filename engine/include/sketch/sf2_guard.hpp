#pragma once
#include <atomic>
#include <cstdlib>

namespace Sketch {

inline constexpr int kSf2Default = 0x69FF9F;

inline auto sf2MaskFromEnv() -> int
{
    char const* v = std::getenv("SKETCH_SF2");
    return (v != nullptr) ? static_cast<int>(std::strtol(v, nullptr, 0)) : kSf2Default;
}
inline auto sf2MaskRef() -> std::atomic<int>&
{
    static std::atomic<int> m{sf2MaskFromEnv()};
    return m;
}
inline auto sf2Mask() -> int { return sf2MaskRef().load(std::memory_order_relaxed); }
inline void setSf2Mask(int m) { sf2MaskRef().store(m < 0 ? sf2MaskFromEnv() : m, std::memory_order_relaxed); }
inline auto sf2On(int mask, int bit) -> bool { return ((mask >> bit) & 1) != 0; }

}
