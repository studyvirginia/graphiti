#pragma once
#include <atomic>
#include <cstdlib>

namespace Sketch {

inline constexpr int kSketchT13Default = 3;

inline auto sketchT13MaskFromEnv() -> int
{
    char const* v = std::getenv("SKETCH_T13");
    return (v != nullptr) ? static_cast<int>(std::strtol(v, nullptr, 0)) : kSketchT13Default;
}
inline auto sketchT13MaskRef() -> std::atomic<int>&
{
    static std::atomic<int> m{sketchT13MaskFromEnv()};
    return m;
}
inline auto sketchT13Mask() -> int { return sketchT13MaskRef().load(std::memory_order_relaxed); }
inline void setSketchT13Mask(int m) { sketchT13MaskRef().store(m < 0 ? sketchT13MaskFromEnv() : m, std::memory_order_relaxed); }
inline auto sketchT13On(int mask, int bit) -> bool { return ((mask >> bit) & 1) != 0; }

inline auto nmlDedupeEnabledFromEnv() -> bool
{
    static bool const v = [] {
        char const* e = std::getenv("SKETCH_NML_DEDUPE");
        return (e == nullptr) || (std::atoi(e) != 0);
    }();
    return v;
}

}
