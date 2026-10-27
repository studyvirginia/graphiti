#pragma once
#include <cstdint>

namespace h156 {
inline bool g_on = false;
inline uint32_t g_touch = 0;

inline void touch_u32(uint32_t* p, uint32_t cur) {
  if (__atomic_load_n(p, __ATOMIC_RELAXED) != cur) __atomic_store_n(p, cur, __ATOMIC_RELAXED);
}
}
