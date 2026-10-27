#pragma once
#include <atomic>
#include <cstdlib>

namespace Sketch {

inline constexpr int kOp2Default = 7;

inline auto op2MaskFromEnv() -> int
{
    char const* v = std::getenv("SKETCH_OP2");
    return (v != nullptr) ? static_cast<int>(std::strtol(v, nullptr, 0)) : kOp2Default;
}
inline auto op2MaskRef() -> std::atomic<int>&
{
    static std::atomic<int> m{op2MaskFromEnv()};
    return m;
}
inline auto op2Mask() -> int { return op2MaskRef().load(std::memory_order_relaxed); }
inline void setOp2Mask(int m) { op2MaskRef().store(m < 0 ? op2MaskFromEnv() : m, std::memory_order_relaxed); }
inline auto op2On(int mask, int bit) -> bool { return ((mask >> bit) & 1) != 0; }

inline auto op2DiagEnabled() -> bool
{
    static bool const v = std::getenv("SKETCH_OP2_DIAG") != nullptr;
    return v;
}

}
