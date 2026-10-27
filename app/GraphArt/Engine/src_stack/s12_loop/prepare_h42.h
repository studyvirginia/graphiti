#pragma once
#include <cstdint>
#include <cstring>
#include "corder_su_h42.h"

namespace h42 {
struct PrepBuf { int opnd[MAXN + 1]; int rep[MAXN]; };

struct PrepScratch {
  uint64_t key[2 * MAXN];
  int val[2 * MAXN];
  uint32_t kst[2 * MAXN];
  uint32_t kgen;
  uint32_t vst[MAXN];
  uint32_t vgen;
  uint32_t rst[MAXN];
  uint32_t rgen;
};
static thread_local PrepScratch ps;

static inline void next_gen(uint32_t& g, uint32_t* arr, int len) {
  if (++g == 0) { for (int i = 0; i < len; i++) arr[i] = 0; g = 1; }
}

static inline int build_F_st(const h14::Tab& C, int r, int* out, uint64_t* fbq, int* pend, const int* EA, const int* EB, int NN) {
  if (C.opnd[r] >= 0) return 0;
  const uint32_t g = ps.vgen, rg = ps.rgen;
  int stk[MAXN]; unsigned char ph[MAXN];
  int sp = 1, cnt = 0;
  stk[0] = r; ph[0] = 0;
  while (sp > 0) {
    const int t = sp - 1, x = stk[t];
    if (ph[t] == 0) {
      ph[t] = 1;
      const int c = EA[x];
      if (c != NN && ps.vst[c] != g) { stk[sp] = c; ph[sp] = 0; ++sp; }
      continue;
    }
    if (ph[t] == 1) {
      ph[t] = 2;
      const int c = EB[x];
      if (c != NN && ps.vst[c] != g) { stk[sp] = c; ph[sp] = 0; ++sp; }
      continue;
    }
    out[cnt++] = x; ps.vst[x] = g; --sp;
    fbq[x >> 6] |= 1ull << (x & 63);
    const int inc = ps.rst[x] != rg;
    ps.rst[x] = rg;
    pend[EA[x]] += inc;
    pend[EB[x]] += inc;
  }
  return cnt;
}

static void prepare(const View& V, PrepBuf& P) {
  OrdScratch& s = os;
  const int n = V.n();
  const int NN = n;
  for (int i = 0; i < n; i++) {
    const int o = V.op(i);
    P.opnd[i] = (i < V.nS) ? 64 + i : (o == 2 ? 900 : (o == 16 ? 901 : -1));
    P.rep[i] = i;
  }
  P.opnd[n] = 0;

  int cap = 16; while (cap < 2 * n) cap <<= 1;
  const uint64_t mask = (uint64_t)cap - 1;
  next_gen(ps.kgen, ps.kst, 2 * MAXN);
  for (int i = V.nS; i < n; i++) {
    const int o = V.op(i);
    if (o == 2 || o == 16) continue;
    int a = P.rep[V.a(i)], b = (o >= 3 && o <= 6) ? P.rep[V.b(i)] : -1;
    if ((o == 3 || o == 5) && a > b) { const int t = a; a = b; b = t; }
    const uint64_t k = ((uint64_t)(o + 1) << 40) | ((uint64_t)(a + 1) << 20) | (uint64_t)(b + 1);
    uint64_t p = ((k * 0x9E3779B97F4A7C15ull) >> 32) & mask;
    for (;;) {
      if (ps.kst[p] != ps.kgen) { ps.kst[p] = ps.kgen; ps.key[p] = k; ps.val[p] = i; break; }
      if (ps.key[p] == k) { P.rep[i] = ps.val[p]; break; }
      p = (p + 1) & mask;
    }
  }

  h14::Tab& C = h14::tab;
  const int* rep = P.rep;
  const int* opnd = P.opnd;
  C.n = n; C.Q = V.Q; C.opnd = opnd;
  for (int x = 0; x < n; x++) {
    if (rep[x] != x || opnd[x] >= 0) { C.OP[x] = -1; C.CA[x] = -1; C.CB[x] = -1; }
    else {
      const int o = V.op(x);
      C.OP[x] = o; C.CA[x] = rep[V.a(x)]; C.CB[x] = (o >= 3 && o <= 6) ? rep[V.b(x)] : -1;
    }
    const int A = C.CA[x], B = C.CB[x];
    s.EA[x] = (A >= 0 && opnd[A] < 0) ? A : NN;
    s.EB[x] = (bin_op(C.OP[x]) && B >= 0 && opnd[B] < 0) ? B : NN;
  }
  s.EA[NN] = NN; s.EB[NN] = NN;
  for (int x = 0; x <= n; x++) s.pend[x] = 0;
  next_gen(ps.rgen, ps.rst, MAXN);
  const int W = (n + 63) >> 6;
  int fill = 0;
  for (int q = 0; q < V.Q; q++) {
    C.rr[q] = rep[V.root(q)];
    C.fo[q] = fill;
    for (int w = 0; w < W; w++) s.fb[q][w] = 0;
    next_gen(ps.vgen, ps.vst, MAXN);
    fill += build_F_st(C, C.rr[q], C.fl + fill, s.fb[q], s.pend, s.EA, s.EB, NN);
  }
  C.fo[V.Q] = fill;
  for (int q = 0; q < V.Q; q++) { const int r = C.rr[q]; s.pend[opnd[r] < 0 ? r : NN]++; }
}
}
