#pragma once
#include "compile_tree_h12.h"
#include "corder_h14.h"

static bool compile_tree_corder_h14(const View& V, Tree& T) {
  if (V.Q > 17 || V.n() > h14::MAXN) return compile_tree_corder_h12(V, T);
  T.Q = V.Q; T.nth = V.nth; T.nS = V.nS; T.prog.clear(); T.sp.clear(); T.th.clear();
  if (V.nS > NSV || V.nth > 32) return false;
  T.sp.reserve((size_t)V.nS);
  for (int i = 0; i < V.nS; i++) { int o = V.op(i), a = V.a(i), b = V.b(i); if (a > 31 || b > 31) return false; T.sp.push_back((uint16_t)(o | (a << 6) | (b << 11))); }
  const int n = V.n();
  int opnd[h14::MAXN + 1], rep[h14::MAXN];
  for (int i = 0; i < n; i++) {
    const int o = V.op(i);
    opnd[i] = (i < V.nS) ? 64 + i : (o == 2 ? 900 : (o == 16 ? 901 : -1));
    rep[i] = i;
  }
  opnd[n] = 0;
  H14_P(h14::PR_HDR);
  {
    int cap = 16; while (cap < 2 * n) cap <<= 1;
    uint64_t key[2 * h14::MAXN]; int val[2 * h14::MAXN];
    memset(key, 0, sizeof(uint64_t) * (size_t)cap);
    const uint64_t mask = (uint64_t)cap - 1;
    for (int i = V.nS; i < n; i++) {
      const int o = V.op(i);
      if (o == 2 || o == 16) continue;
      int a = rep[V.a(i)], b = (o >= 3 && o <= 6) ? rep[V.b(i)] : -1;
      if ((o == 3 || o == 5) && a > b) { const int t = a; a = b; b = t; }
      const uint64_t k = ((uint64_t)(o + 1) << 40) | ((uint64_t)(a + 1) << 20) | (uint64_t)(b + 1);
      uint64_t p = ((k * 0x9E3779B97F4A7C15ull) >> 32) & mask;
      for (;;) {
        if (key[p] == 0) { key[p] = k; val[p] = i; break; }
        if (key[p] == k) { rep[i] = val[p]; break; }
        p = (p + 1) & mask;
      }
    }
  }
  H14_P(h14::PR_CSE);
  h14::Tab& C = h14::tab;
  h14::build_tables(C, V, rep, opnd, n);
  int ord[h14::MAXQ];
  for (int q = 0; q < V.Q; q++) ord[q] = q;
  if (V.Q > 1 && !h14::corder_h14(C, ord)) return compile_tree_corder_h12(V, T);

  int ns[h14::MAXM], as[h14::MAXM], bs[h14::MAXM]; const h14::Seq S{ns, as, bs};
  int reads[h14::MAXN], lastread[h14::MAXN];
  for (int x = 0; x < n; x++) { reads[x] = 0; lastread[x] = -1; }
  uint64_t D[h14::NWD]; memset(D, 0, sizeof D);
  int segend[h14::MAXQ]; int m = 0;
  for (int qi = 0; qi < V.Q; qi++) { m = h14::emit_q(C, ord[qi], D, S, reads, lastread, m); segend[qi] = m; }
  for (int k = 0; k < m; k++) {
    if (ns[k] == n) ns[k] = -1;
    if (bs[k] == n) bs[k] = -1;
  }
  H14_P(h14::PR_EMIT);

  int pos[h14::MAXN]; unsigned char chain[h14::MAXN];
  for (int x = 0; x < n; x++) { pos[x] = -1; chain[x] = 0; }
  for (int k = 0; k < m; k++) if (ns[k] >= 0) pos[ns[k]] = k;
  for (int k = 0; k + 1 < m; k++) {
    const int x = ns[k];
    if (x < 0 || reads[x] == 0) continue;
    const int r1 = (as[k + 1] == x) + (bs[k + 1] == x);
    if (r1 == reads[x]) chain[x] = 1;
  }
  H14_P(h14::PR_PREP);

  int slot[h14::MAXN], freeS[h14::MAXM]; int nf = 0, nsl = 0;
  for (int x = 0; x < n; x++) slot[x] = -1;
  uint16_t capt[18]; for (int q = 0; q < 18; q++) capt[q] = 1022;
  int qof[h14::MAXM], newend[h14::MAXM + 1];
  newend[0] = 0;
  { int q = 0; for (int k = 0; k < m; k++) { qof[k] = q; if (k + 1 == segend[q]) q++; } }
  auto code = [&](int x, int k) -> int {
    if (opnd[x] >= 0) return opnd[x];
    if (chain[x] && pos[x] == k - 1) return 1022;
    return slot[x];
  };
  T.prog.clear(); T.prog.reserve((size_t)m);
  for (int k = 0; k < m; k++) {
    const int node = ns[k], A0 = as[k], B0 = bs[k];
    const int A = code(A0, k), B = B0 >= 0 ? code(B0, k) : 0;
    if (A < 0 || B < 0) return false;
    if (opnd[A0] < 0 && slot[A0] >= 0 && lastread[A0] == k) { freeS[nf++] = slot[A0]; slot[A0] = -2; }
    if (B0 >= 0 && B0 != A0 && opnd[B0] < 0 && slot[B0] >= 0 && lastread[B0] == k) { freeS[nf++] = slot[B0]; slot[B0] = -2; }
    int d = 63;
    if (node >= 0 && reads[node] > 0 && !chain[node]) {
      if (nf > 0) d = freeS[--nf]; else d = nsl++;
      slot[node] = d;
    }
    if (node < 0) { capt[qof[k]] = (uint16_t)A; newend[k + 1] = (int)T.prog.size(); continue; }
    T.prog.push_back((uint32_t)C.OP[node] | (uint32_t)d << 6 | (uint32_t)A << 12 | (uint32_t)B << 22);
    newend[k + 1] = (int)T.prog.size();
  }
  H14_P(h14::PR_SLOT);

  for (int q = 0; q < 36; q++) T.segs[q] = 0;
  T.segs[0] = 0; for (int q = 0; q < V.Q; q++) T.segs[q + 1] = (uint16_t)newend[segend[q]];
  for (int q = 0; q < 18; q++) T.segs[18 + q] = 1022;
  for (int q = 0; q < V.Q; q++) T.segs[18 + q] = capt[q];
  T.nslot = nsl; if (nsl > MAXSLOT) return false;
  T.th.reserve((size_t)V.nth);
  for (int s = 0; s < V.nth; s++) T.th.push_back(V.th(s));
  H14_P(h14::PR_SEGS);
  return true;
}
