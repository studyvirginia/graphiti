#pragma once
#include <cstdint>
#include <cstring>

namespace su30 {
constexpr int MAXN = 512, MAXQ = 64;

struct G { int n, Q; const int* opnd; const int* OP; const int* CA; const int* CB; const int* rr; };

static inline bool bin_op(int op) { return op >= 3 && op <= 6; }

static int build_F(const G& g, int r, int* out, int* vis, int stamp) {
  int stx[MAXN], sts[MAXN];
  int sp = 1, cnt = 0;
  stx[0] = r; sts[0] = 0;
  while (sp > 0) {
    const int t = sp - 1, x = stx[t];
    if (sts[t] == 0) {
      if (g.opnd[x] >= 0 || vis[x] == stamp) { --sp; continue; }
      sts[t] = 1; stx[sp] = g.CA[x]; sts[sp] = 0; ++sp; continue;
    }
    if (sts[t] == 1 && bin_op(g.OP[x])) { sts[t] = 2; stx[sp] = g.CB[x]; sts[sp] = 0; ++sp; continue; }
    out[cnt++] = x; vis[x] = stamp; --sp;
  }
  return cnt;
}

struct St { int pend[MAXN]; unsigned char held[MAXN], prod[MAXN]; int live; };
struct Ulog { int x, pend; unsigned char held, prod; };

static int sim(St& s, const G& g, const int* L, int len, int r, bool commit, int& endLive, Ulog* lg) {
  int cur = s.live, peak = cur, nl = 0;
  auto touch = [&](int x) { if (!commit) lg[nl++] = Ulog{x, s.pend[x], s.held[x], s.prod[x]}; };
  if (g.opnd[r] < 0) {
    touch(r);
    if (--s.pend[r] == 0 && s.held[r]) { s.held[r] = 0; --cur; }
  }
  for (int t = 0; t < len; t++) {
    const int x = L[t];
    const int A = g.CA[x];
    if (g.opnd[A] < 0) {
      touch(A);
      if (--s.pend[A] == 0 && s.held[A]) { s.held[A] = 0; --cur; }
    }
    if (bin_op(g.OP[x])) {
      const int B = g.CB[x];
      if (B >= 0 && g.opnd[B] < 0) {
        touch(B);
        if (--s.pend[B] == 0 && s.held[B]) { s.held[B] = 0; --cur; }
      }
    }
    touch(x);
    s.prod[x] = 1;
    const int nx = (t + 1 < len) ? L[t + 1] : -1;
    int mu = 0;
    if (nx >= 0) { if (g.CA[nx] == x) ++mu; if (bin_op(g.OP[nx]) && g.CB[nx] == x) ++mu; }
    if (s.pend[x] > 0 && !(nx >= 0 && mu == s.pend[x])) { s.held[x] = 1; ++cur; }
    if (cur > peak) peak = cur;
  }
  endLive = cur;
  if (commit) s.live = cur;
  else for (int i = nl - 1; i >= 0; --i) { const Ulog& u = lg[i]; s.pend[u.x] = u.pend; s.held[u.x] = u.held; s.prod[u.x] = u.prod; }
  return peak;
}

static bool order(const G& g, int* out) {
  const int n = g.n, Q = g.Q;
  if (Q < 1 || Q > MAXQ || n < 1 || n > MAXN) return false;
  if (Q == 1) { out[0] = 0; return true; }
  static thread_local int fl[MAXQ * MAXN], buf[MAXQ * MAXN];
  static thread_local int fo[MAXQ + 1];
  static thread_local St s;
  static thread_local Ulog lg[4 * MAXN + 8];
  int vis[MAXN];
  for (int x = 0; x < n; x++) vis[x] = 0;
  int fill = 0;
  for (int q = 0; q < Q; q++) { fo[q] = fill; fill += build_F(g, g.rr[q], fl + fill, vis, q + 1); }
  fo[Q] = fill;

  unsigned char inR[MAXN];
  for (int x = 0; x < n; x++) { inR[x] = 0; s.pend[x] = 0; s.held[x] = 0; s.prod[x] = 0; }
  s.live = 0;
  for (int t = 0; t < fill; t++) inR[fl[t]] = 1;
  for (int x = 0; x < n; x++) {
    if (!inR[x]) continue;
    const int A = g.CA[x];
    if (g.opnd[A] < 0) s.pend[A]++;
    if (bin_op(g.OP[x])) { const int B = g.CB[x]; if (B >= 0 && g.opnd[B] < 0) s.pend[B]++; }
  }
  for (int q = 0; q < Q; q++) { const int r = g.rr[q]; if (g.opnd[r] < 0) s.pend[r]++; }

  int rem[MAXQ]; int nr = Q;
  for (int q = 0; q < Q; q++) rem[q] = q;
  for (int step = 0; step < Q; step++) {
    int lb[MAXQ], le[MAXQ]; int bp = 0;
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      lb[a] = bp;
      for (int t = fo[q]; t < fo[q + 1]; t++) { const int x = fl[t]; if (!s.prod[x]) buf[bp++] = x; }
      le[a] = bp;
    }
    int bestP = 0x7fffffff, bestE = 0x7fffffff, ba = 0;
    for (int a = 0; a < nr; a++) {
      int endL = 0;
      const int p = sim(s, g, buf + lb[a], le[a] - lb[a], g.rr[rem[a]], false, endL, lg);
      if (p < bestP || (p == bestP && endL < bestE)) { bestP = p; bestE = endL; ba = a; }
    }
    int endL = 0;
    sim(s, g, buf + lb[ba], le[ba] - lb[ba], g.rr[rem[ba]], true, endL, lg);
    out[step] = rem[ba];
    for (int a = ba; a + 1 < nr; a++) rem[a] = rem[a + 1];
    --nr;
  }
  return true;
}

static bool order_static(const G& g, int* out) {
  const int n = g.n, Q = g.Q;
  if (Q < 1 || Q > MAXQ || n < 1 || n > MAXN) return false;
  static thread_local int lst[MAXN];
  int need[MAXN], vis[MAXN];
  for (int x = 0; x < n; x++) { need[x] = 0; vis[x] = 0; }
  int key[MAXQ];
  for (int q = 0; q < Q; q++) {
    const int cnt = build_F(g, g.rr[q], lst, vis, q + 1);
    for (int t = 0; t < cnt; t++) {
      const int x = lst[t];
      const int A = g.CA[x];
      const int na = (g.opnd[A] >= 0) ? 0 : need[A];
      if (bin_op(g.OP[x])) {
        const int B = g.CB[x];
        const int nb = (B >= 0 && g.opnd[B] < 0) ? need[B] : 0;
        need[x] = (na == nb) ? na + 1 : (na > nb ? na : nb);
      } else {
        need[x] = na > 1 ? na : 1;
      }
    }
    key[q] = (g.opnd[g.rr[q]] >= 0) ? 0 : need[g.rr[q]];
  }
  for (int q = 0; q < Q; q++) out[q] = q;
  for (int i = 1; i < Q; i++) {
    const int v = out[i];
    int j = i - 1;
    while (j >= 0 && key[out[j]] < key[v]) { out[j + 1] = out[j]; --j; }
    out[j + 1] = v;
  }
  return true;
}

static int peak_of(const G& g, const int* ord, int* mOut) {
  const int n = g.n, Q = g.Q;
  static thread_local int sN[MAXN + MAXQ + 8], sA[MAXN + MAXQ + 8], sB[MAXN + MAXQ + 8];
  unsigned char done[MAXN];
  for (int x = 0; x < n; x++) done[x] = 0;
  int m = 0;
  for (int qi = 0; qi < Q; qi++) {
    const int r = g.rr[ord[qi]];
    const int m0 = m;
    int stx[MAXN], sts[MAXN]; int sp = 1; stx[0] = r; sts[0] = 0;
    while (sp > 0) {
      const int t = sp - 1, x = stx[t];
      if (sts[t] == 0) { if (g.opnd[x] >= 0 || done[x]) { --sp; continue; } sts[t] = 1; stx[sp] = g.CA[x]; sts[sp] = 0; ++sp; continue; }
      if (sts[t] == 1 && bin_op(g.OP[x])) { sts[t] = 2; stx[sp] = g.CB[x]; sts[sp] = 0; ++sp; continue; }
      sN[m] = x; sA[m] = g.CA[x]; sB[m] = bin_op(g.OP[x]) ? g.CB[x] : -1; ++m; done[x] = 1; --sp;
    }
    if (m == m0) { sN[m] = -1; sA[m] = r; sB[m] = -1; ++m; }
  }
  int rd[MAXN], lr[MAXN], pos[MAXN];
  for (int x = 0; x < n; x++) { rd[x] = 0; lr[x] = -1; pos[x] = -1; }
  for (int k = 0; k < m; k++) {
    const int a = sA[k], b = sB[k];
    if (a >= 0 && g.opnd[a] < 0) { rd[a]++; lr[a] = k; }
    if (b >= 0 && g.opnd[b] < 0) { rd[b]++; lr[b] = k; }
    if (sN[k] >= 0) pos[sN[k]] = k;
  }
  int chv[MAXN + MAXQ + 8];
  for (int k = 0; k < m; k++) {
    const int z = sN[k];
    if (z < 0 || rd[z] == 0) { chv[k] = 0; continue; }
    const int cnt = (k + 1 < m) ? ((sA[k + 1] == z) + (sB[k + 1] == z)) : 0;
    chv[k] = (cnt != rd[z]);
  }
  int live = 0, mx = 0;
  for (int k = 0; k < m; k++) {
    const int za = sA[k], zb = sB[k];
    const int fa = (za >= 0 && g.opnd[za] < 0 && lr[za] == k) ? chv[pos[za]] : 0;
    const int fb = (zb >= 0 && zb != za && g.opnd[zb] < 0 && lr[zb] == k) ? chv[pos[zb]] : 0;
    live += chv[k] - fa - fb;
    if (live > mx) mx = live;
  }
  if (mOut) *mOut = m;
  return mx;
}
}
