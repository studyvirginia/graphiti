#pragma once
#include <cstdint>
#include <cstring>
#include "../s02_nf/sk_nf_h2b.hpp"

#ifdef H14_PROF
void h14_probe(int id);
#define H14_P(id) h14_probe(id)
#else
#define H14_P(id) ((void)0)
#endif

namespace h14 {
constexpr int MAXN = 512, MAXQ = 64, MAXM = MAXN + MAXQ + 4, NWD = MAXN / 64;
constexpr int INF = 1 << 30, NOCUT = 0x7fffffff;
enum { PR_TABLES = 0, PR_STEP, PR_TRIAL_EMIT, PR_TRIAL_SCAN, PR_COMMIT,
       PR_HDR, PR_CSE, PR_EMIT, PR_PREP, PR_SLOT, PR_SEGS, PR_NSEC };

struct Tab {
  int n, Q;
  const int* opnd;
  int OP[MAXN], CA[MAXN], CB[MAXN];
  int rr[MAXQ];
  int fo[MAXQ + 1];
  int fl[MAXQ * MAXN];
  int gl[MAXQ * MAXN];
  int gb[MAXQ], ge[MAXQ];
};
static thread_local Tab tab;

static inline bool tst(const uint64_t* b, int x) { return (b[x >> 6] >> (x & 63)) & 1ull; }
static inline void put(uint64_t* b, int x) { b[x >> 6] |= 1ull << (x & 63); }

static inline int build_F(const Tab& C, int r, int* out, unsigned char* vis) {
  int stx[MAXN], sts[MAXN];
  int sp = 1, cnt = 0;
  stx[0] = r; sts[0] = 0;
  while (sp > 0) {
    const int t = sp - 1, x = stx[t];
    if (sts[t] == 0) {
      if (C.opnd[x] >= 0 || vis[x]) { --sp; continue; }
      sts[t] = 1; stx[sp] = C.CA[x]; sts[sp] = 0; ++sp; continue;
    }
    if (sts[t] == 1 && C.OP[x] >= 3 && C.OP[x] <= 6) { sts[t] = 2; stx[sp] = C.CB[x]; sts[sp] = 0; ++sp; continue; }
    out[cnt++] = x; vis[x] = 1; --sp;
  }
  return cnt;
}

static inline void build_tables(Tab& C, const skn::View& V, const int* rep, const int* opnd, int n) {
  C.n = n; C.Q = V.Q; C.opnd = opnd;
  for (int x = 0; x < n; x++) {
    if (rep[x] != x || opnd[x] >= 0) { C.OP[x] = -1; C.CA[x] = -1; C.CB[x] = -1; continue; }
    const int o = V.op(x);
    C.OP[x] = o; C.CA[x] = rep[V.a(x)]; C.CB[x] = (o >= 3 && o <= 6) ? rep[V.b(x)] : -1;
  }
  unsigned char vis[MAXN];
  int fill = 0;
  for (int q = 0; q < V.Q; q++) {
    C.rr[q] = rep[V.root(q)];
    C.fo[q] = fill;
    memset(vis, 0, (size_t)n);
    fill += build_F(C, C.rr[q], C.fl + fill, vis);
  }
  C.fo[V.Q] = fill;
}

struct Seq { int* sn; int* sa; int* sb; };

static inline int emit_q(const Tab& C, int q, uint64_t* D, const Seq& S, int* rd, int* lr, int m) {
  const int n = C.n, m0 = m;
  const int e = C.fo[q + 1];
  for (int t = C.fo[q]; t < e; t++) {
    const int x = C.fl[t];
    if (tst(D, x)) continue;
    const int A = C.CA[x], B = (C.OP[x] >= 3 && C.OP[x] <= 6) ? C.CB[x] : n;
    S.sn[m] = x; S.sa[m] = A; S.sb[m] = B;
    rd[A]++; lr[A] = m;
    if (B != n) { rd[B]++; lr[B] = m; }
    ++m;
    put(D, x);
  }
  if (m == m0) { const int r = C.rr[q]; S.sn[m] = n; S.sa[m] = r; S.sb[m] = n; rd[r]++; lr[r] = m; ++m; }
  return m;
}

static inline int emit_g(const Tab& C, int b, int e, int r, uint64_t* Dt, const Seq& S, int* rd, int* lr, int m) {
  const int n = C.n, m0 = m;
  const int* gl = C.gl;
  for (int t = b; t < e; t++) {
    const int x = gl[t];
    const uint64_t done = (Dt[x >> 6] >> (x & 63)) & 1ull;
    const int inc = (int)(done ^ 1ull);
    const int A = C.CA[x], Bn = C.CB[x];
    const int hasB = (Bn >= 0);
    const int Bi = hasB ? Bn : n;
    const int incB = inc & hasB;
    S.sn[m] = x; S.sa[m] = A; S.sb[m] = Bi;
    rd[A] += inc; lr[A] = inc ? m : lr[A];
    rd[Bi] += incB; lr[Bi] = incB ? m : lr[Bi];
    m += inc;
    Dt[x >> 6] |= (uint64_t)inc << (x & 63);
  }
  if (m == m0) { S.sn[m] = n; S.sa[m] = r; S.sb[m] = n; rd[r]++; lr[r] = m; ++m; }
  return m;
}

static inline int scan_two(const Tab& C, const Seq& S, const int* rd, const int* lr, int m, int lim, int thr,
                           int* chv, int* posOf) {
  const int* opnd = C.opnd;
  for (int k = 0; k < m; k++) {
    const int zn = S.sn[k];
    const int rz = rd[zn];
    const int cnt = (S.sa[k + 1] == zn) + (S.sb[k + 1] == zn);
    chv[k] = (rz > 0) & (cnt != rz);
    posOf[zn] = k;
  }
  int live = 0, mx = 0;
  for (int k = 0; k < lim; k++) {
    const int za = S.sa[k], zb = S.sb[k];
    const int fa = (opnd[za] < 0) & (lr[za] == k) & chv[posOf[za]];
    const int fb = (zb != za) & (opnd[zb] < 0) & (lr[zb] == k) & chv[posOf[zb]];
    live += chv[k] - fa - fb;
    mx = live > mx ? live : mx;
    if (mx * 1000 + m >= thr) return NOCUT;
  }
  return mx;
}

static bool corder_h14(Tab& C, int* out) {
  const int Q = C.Q, n = C.n;
  if (Q < 2 || Q > MAXQ || n < 1 || n > MAXN) return false;
  const Tab& CT = C;
  int sn[MAXM + 1], sa[MAXM + 1], sb[MAXM + 1]; const Seq S{sn, sa, sb};
  uint64_t D[NWD]; memset(D, 0, sizeof D);
  int rdP[MAXN + 1], lrP[MAXN + 1], rdT[MAXN + 1], lrT[MAXN + 1], chv[MAXM + 1], posOf[MAXN + 1];
  for (int x = 0; x <= n; x++) { rdP[x] = 0; lrP[x] = -1; posOf[x] = 0; }
  int rem[MAXQ]; int nr = Q; for (int q = 0; q < Q; q++) rem[q] = q;
  int mP = 0, no = 0;
  H14_P(PR_TABLES);
  while (nr > 0) {
    int gf = 0;
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      C.gb[q] = gf;
      const int e = CT.fo[q + 1];
      for (int t = CT.fo[q]; t < e; t++) { const int x = CT.fl[t]; if (!tst(D, x)) C.gl[gf++] = x; }
      C.ge[q] = gf;
    }
    H14_P(PR_STEP);
    int bc = INF, bi = -1;
    for (int i = 0; i < nr; i++) {
      memcpy(rdT, rdP, sizeof(int) * (size_t)(n + 1));
      memcpy(lrT, lrP, sizeof(int) * (size_t)(n + 1));
      uint64_t Dt[NWD]; memcpy(Dt, D, sizeof Dt);
      int m = mP;
      { const int q = rem[i]; m = emit_g(CT, C.gb[q], C.ge[q], CT.rr[q], Dt, S, rdT, lrT, m); }
      for (int j = 0; j < nr; j++) if (j != i) { const int q = rem[j]; m = emit_g(CT, C.gb[q], C.ge[q], CT.rr[q], Dt, S, rdT, lrT, m); }
      S.sn[m] = n; S.sa[m] = n; S.sb[m] = n;
      H14_P(PR_TRIAL_EMIT);
      const int mx = scan_two(CT, S, rdT, lrT, m, m, bc, chv, posOf);
      H14_P(PR_TRIAL_SCAN);
      if (mx == NOCUT) continue;
      const int c = mx * 1000 + m;
      if (c < bc) { bc = c; bi = i; }
    }
    out[no++] = rem[bi];
    mP = emit_q(CT, rem[bi], D, S, rdP, lrP, mP);
    for (int a = bi; a + 1 < nr; a++) rem[a] = rem[a + 1];
    --nr;
    H14_P(PR_COMMIT);
  }
  return true;
}
}
