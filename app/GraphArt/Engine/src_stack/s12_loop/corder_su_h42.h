#pragma once
#include <cstdint>
#include <cstring>
#include "corder_h14.h"

namespace h42 {
constexpr int MAXN = 512, MAXQ = 64, WN = MAXN / 64;

struct Undo { int x, pend; unsigned char held, prod; };

struct OrdScratch {
  int pend[MAXN + 1];
  unsigned char held[MAXN + 1], prod[MAXN + 1];
  int EA[MAXN + 1], EB[MAXN + 1];
  int L[MAXN + 2];
  Undo lg[3 * MAXN + 8];
  uint64_t fb[MAXQ][WN];
  uint64_t T[WN + 1];
  int rp[MAXQ], re[MAXQ];
  unsigned char valid[MAXQ];
};
static thread_local OrdScratch os;

static inline bool bin_op(int op) { return op >= 3 && op <= 6; }

template <bool COMMIT>
static inline void sim(OrdScratch& s, const h14::Tab& C, int q, int n, int& peakR, int& endR) {
  const int NN = n;
  const int* opnd = C.opnd;
  const int r = C.rr[q];
  const int re = opnd[r] < 0 ? r : NN;
  int len = 0;
  for (int t = C.fo[q]; t < C.fo[q + 1]; t++) { const int x = C.fl[t]; s.L[len] = x; len += (int)(s.prod[x] ^ 1); }
  s.L[len] = NN;
  int nl = 0, cur = 0, peak = 0;

  auto dec = [&](int a) {
    if (!COMMIT) s.lg[nl++] = Undo{a, s.pend[a], s.held[a], s.prod[a]};
    else s.T[a >> 6] |= 1ull << (a & 63);
    const int z = (--s.pend[a]) == 0;
    cur -= (int)(s.held[a] & (unsigned char)z);
    s.held[a] &= (unsigned char)(z ^ 1);
  };

  dec(re);
  for (int t = 0; t < len; t++) {
    const int x = s.L[t];
    dec(s.EA[x]);
    dec(s.EB[x]);
    if (!COMMIT) s.lg[nl++] = Undo{x, s.pend[x], s.held[x], s.prod[x]};
    else s.T[x >> 6] |= 1ull << (x & 63);
    s.prod[x] = 1;
    const int nx = s.L[t + 1];
    const int mu = (s.EA[nx] == x) + (s.EB[nx] == x);
    const int h = s.pend[x] > mu;
    s.held[x] = (unsigned char)h;
    cur += h;
    peak = cur > peak ? cur : peak;
  }
  peakR = peak; endR = cur;
  if (!COMMIT) {
    for (int i = nl - 1; i >= 0; --i) { const Undo& u = s.lg[i]; s.pend[u.x] = u.pend; s.held[u.x] = u.held; s.prod[u.x] = u.prod; }
  }
}

static bool order(const h14::Tab& C, int* out) {
  const int n = C.n, Q = C.Q;
  if (Q < 1 || Q > MAXQ || n < 1 || n > MAXN) return false;
  if (Q == 1) { out[0] = 0; return true; }
  OrdScratch& s = os;
  const int W = (n + 63) >> 6;
  for (int x = 0; x <= n; x++) { s.held[x] = 0; s.prod[x] = 0; }
  for (int w = 0; w <= WN; w++) s.T[w] = 0;

  int rem[MAXQ]; int nr = Q;
  for (int q = 0; q < Q; q++) { rem[q] = q; s.valid[q] = 0; }
  for (int step = 0; step < Q; step++) {
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      if (!s.valid[q]) { int p, e; sim<false>(s, C, q, n, p, e); s.rp[q] = p; s.re[q] = e; s.valid[q] = 1; }
    }
    int ba = 0;
    int bp = s.rp[rem[0]], be = s.re[rem[0]];
    for (int a = 1; a < nr; a++) {
      const int q = rem[a];
      if (s.rp[q] < bp || (s.rp[q] == bp && s.re[q] < be)) { bp = s.rp[q]; be = s.re[q]; ba = a; }
    }
    const int c = rem[ba];
    for (int w = 0; w <= WN; w++) s.T[w] = 0;
    { int p, e; sim<true>(s, C, c, n, p, e); }
    out[step] = c;
    for (int a = ba; a + 1 < nr; a++) rem[a] = rem[a + 1];
    --nr;
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      uint64_t hit = 0;
      for (int w = 0; w < W; w++) hit |= s.fb[q][w] & s.T[w];
      if (hit) s.valid[q] = 0;
    }
  }
  return true;
}
}
