#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>
#include "final_pick_h.hpp"

namespace sp {
constexpr int NMAX = 24;
constexpr int KMAX = 96;
constexpr int NPT = 48;
struct Tab {
  std::vector<uint16_t> lev;
  std::vector<double> val;
  Tab() {
    struct E { double p; int idx; };
    std::vector<E> e; e.reserve((size_t)NMAX * NMAX * (KMAX + 1));
    lev.assign((size_t)(NMAX + 1) * (NMAX + 1) * (KMAX + 1), 0);
    for (int nx = 1; nx <= NMAX; nx++) for (int ny = 1; ny <= NMAX; ny++) for (int k = 0; k <= KMAX; k++)
      e.push_back({fp::search_pen(k, NPT, nx, ny), (int)(((size_t)nx * (NMAX + 1) + ny) * (KMAX + 1) + k)});
    std::sort(e.begin(), e.end(), [](const E& a, const E& b) { return a.p < b.p; });
    val.assign(1, 0.0);
    for (size_t i = 0; i < e.size(); i++) { if (i == 0 || e[i].p != e[i - 1].p) val.push_back(e[i].p); lev[(size_t)e[i].idx] = (uint16_t)(val.size() - 1); }
  }
};
inline const Tab& tab() { static const Tab t; return t; }
inline uint16_t level(int nx, int ny, int k) {
  nx = nx < 1 ? 1 : nx > NMAX ? NMAX : nx; ny = ny < 1 ? 1 : ny > NMAX ? NMAX : ny; k = k < 0 ? 0 : k > KMAX ? KMAX : k;
  return tab().lev[((size_t)nx * (NMAX + 1) + ny) * (KMAX + 1) + k];
}
inline double pen_of_level(int lv) { return tab().val[(size_t)lv]; }
}
