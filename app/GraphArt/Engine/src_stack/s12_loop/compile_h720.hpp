#pragma once
#include "compile_memo_i17_h156.h"
#include "h313a_view.hpp"
#include "h313b.hpp"

#ifndef C308_TS
#define C308_TS(k) ((void)0)
#endif
#ifndef H720_V
#define H720_V 42
#endif
#ifndef H720_ATOM
#define H720_ATOM 0
#endif
namespace c720 {
#ifdef H720_COUNT
static long cn_trial = 0, cn_tlen = 0, cn_tnl = 0, cn_fill = 0, cn_n = 0, cn_pc = 0, cn_rec = 0, cn_steps = 0, cn_commit_len = 0;
#define H720_CNT(x, v) ((x) += (v))
#else
#define H720_CNT(x, v) ((void)0)
#endif
inline bool g_a2 = false;
constexpr int kN = 512, kQ = 64;
using h14::NWD;
using f38::rec_fields;
using f38::Scratch;

struct ColT { int ia[h14::MAXQ * (h14::MAXN + 1)]; int ib[h14::MAXQ * (h14::MAXN + 1)]; int io[h14::MAXQ + 1]; int pos[h14::MAXN]; };
static thread_local ColT ct;

static inline int build_F_c(const h14::Tab& C, int r, int* out, uint64_t* fbq, int* pend, const int* EA, const int* EB, int NN, int* ia, int* ib) {
  using namespace h42;
  if (C.opnd[r] >= 0) { ia[0] = -1; ib[0] = -1; return 0; }
  const uint32_t g = ps.vgen, rg = ps.rgen;
#if H720_V & 8
  int stk[h14::MAXN]; unsigned char ph[h14::MAXN];
  int sp = 0, cnt = 0, x = r; unsigned p = 0;
  for (;;) {
    if (p == 0) {
      const int c = EA[x];
      if (c != NN && ps.vst[c] != g) { stk[sp] = x; ph[sp] = 1; ++sp; x = c; continue; }
    }
    if (p <= 1) {
      const int c = EB[x];
      if (c != NN && ps.vst[c] != g) { stk[sp] = x; ph[sp] = 2; ++sp; x = c; p = 0; continue; }
    }
    out[cnt] = x; ps.vst[x] = g;
    ct.pos[x] = cnt;
    ia[cnt] = EA[x] == NN ? -1 : ct.pos[EA[x]];
    ib[cnt] = EB[x] == NN ? -1 : ct.pos[EB[x]];
    ++cnt;
    fbq[x >> 6] |= 1ull << (x & 63);
    const int inc = ps.rst[x] != rg;
    ps.rst[x] = rg;
    pend[EA[x]] += inc;
    pend[EB[x]] += inc;
    if (sp == 0) break;
    --sp; x = stk[sp]; p = ph[sp];
  }
  ia[cnt] = -1; ib[cnt] = -1;
  return cnt;
#else
  int stk[h14::MAXN]; unsigned char ph[h14::MAXN];
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
    out[cnt] = x; ps.vst[x] = g; --sp;
    ct.pos[x] = cnt;
    ia[cnt] = EA[x] == NN ? -1 : ct.pos[EA[x]];
    ib[cnt] = EB[x] == NN ? -1 : ct.pos[EB[x]];
    ++cnt;
    fbq[x >> 6] |= 1ull << (x & 63);
    const int inc = ps.rst[x] != rg;
    ps.rst[x] = rg;
    pend[EA[x]] += inc;
    pend[EB[x]] += inc;
  }
  ia[cnt] = -1; ib[cnt] = -1;
  return cnt;
#endif
}

static void prepare(const View& V, h42::PrepBuf& P) {
  using namespace h42;
  OrdScratch& s = os;
  const int n = V.n();
  const int NN = n;
  int cap = 16; while (cap < 2 * n) cap <<= 1;
  const uint64_t mask = (uint64_t)cap - 1;
  next_gen(ps.kgen, ps.kst, 2 * kN);
  h14::Tab& C = h14::tab;
  int* rep = P.rep; int* opnd = P.opnd;
  C.n = n; C.Q = V.Q; C.opnd = opnd;
  for (int i = 0; i < n; i++) {
    const int o = V.op(i);
    rep[i] = i;
    s.pend[i] = 0;
    s.EA[i] = NN; s.EB[i] = NN;
    C.OP[i] = -1; C.CA[i] = -1; C.CB[i] = -1;
    if (i < V.nS) { opnd[i] = 64 + i; continue; }
    if (o == 2 || o == 16) { opnd[i] = o == 2 ? 900 : 901; continue; }
    opnd[i] = -1;
    const int ra = rep[V.a(i)], rb = (o >= 3 && o <= 6) ? rep[V.b(i)] : -1;
    int a = ra, b = rb;
    if ((o == 3 || o == 5) && a > b) { const int t = a; a = b; b = t; }
    const uint64_t k = ((uint64_t)(o + 1) << 40) | ((uint64_t)(a + 1) << 20) | (uint64_t)(b + 1);
    uint64_t p = ((k * 0x9E3779B97F4A7C15ull) >> 32) & mask;
    for (;;) {
      if (ps.kst[p] != ps.kgen) { ps.kst[p] = ps.kgen; ps.key[p] = k; ps.val[p] = i; break; }
      if (ps.key[p] == k) { rep[i] = ps.val[p]; break; }
      p = (p + 1) & mask;
    }
    if (rep[i] != i) continue;
    C.OP[i] = o; C.CA[i] = ra; C.CB[i] = rb;
    s.EA[i] = opnd[ra] < 0 ? ra : NN;
    s.EB[i] = (o >= 3 && o <= 6 && opnd[rb] < 0) ? rb : NN;
  }
  opnd[n] = 0;
  s.EA[NN] = NN; s.EB[NN] = NN; s.pend[NN] = 0;
  next_gen(ps.rgen, ps.rst, kN);
  C308_TS(23);
  const int W = (n + 63) >> 6;
  int fill = 0, io = 0;
  for (int q = 0; q < V.Q; q++) {
    C.rr[q] = rep[V.root(q)];
    C.fo[q] = fill;
    ct.io[q] = io;
    for (int w = 0; w < W; w++) s.fb[q][w] = 0;
    next_gen(ps.vgen, ps.vst, kN);
    const int c = build_F_c(C, C.rr[q], C.fl + fill, s.fb[q], s.pend, s.EA, s.EB, NN, ct.ia + io, ct.ib + io);
    fill += c; io += c + 1;
  }
  C.fo[V.Q] = fill; H720_CNT(cn_fill, fill); H720_CNT(cn_n, n); H720_CNT(cn_rec, 1);
  for (int q = 0; q < V.Q; q++) { const int r = C.rr[q]; s.pend[opnd[r] < 0 ? r : NN]++; }
}

#if H720_V & 4
static inline void trial(const h42::OrdScratch& s, const h14::Tab& C, int q, int& peakR, int& endR) {
  constexpr int DUM = h14::MAXN + 1;
  int w[DUM + 1];
  int L[h14::MAXN + 1];
  const int b = C.fo[q], len = C.fo[q + 1] - b;
  const int* fl = C.fl + b;
  const int* ia = ct.ia + ct.io[q];
  const int* ib = ct.ib + ct.io[q];
  int nl = 0;
  H720_CNT(cn_trial, 1); H720_CNT(cn_tlen, len);
  for (int t = 0; t < len; t++) { const int x = fl[t]; w[t] = (s.pend[x] << 1) | (int)s.held[x]; L[nl] = t; nl += (int)(s.prod[x] ^ 1); }
  w[DUM] = (1 << 30) << 1;
  int cur = 0, peak = 0;
  auto dec = [&](int a) {
    const int i = a < 0 ? DUM : a;
    int v = w[i] - 2;
    const int z = v < 2;
    cur -= v & z;
    w[i] = v & ~z;
  };
  dec(len > 0 ? len - 1 : -1);
  for (int i = 0; i < nl; i++) {
    const int t = L[i];
    dec(ia[t]);
    dec(ib[t]);
    int mu = 0;
    if (i + 1 < nl) { const int nx = L[i + 1]; mu = (ia[nx] == t) + (ib[nx] == t); }
    const int h = (w[t] >> 1) > mu;
    w[t] = (w[t] & ~1) | h;
    cur += h;
    peak = cur > peak ? cur : peak;
  }
  H720_CNT(cn_tnl, nl);
  peakR = peak; endR = cur;
}

#else
static inline void trial(const h42::OrdScratch& s, const h14::Tab& C, int q, int& peakR, int& endR) {
  constexpr int DUM = h14::MAXN + 1;
  int pe[DUM + 1]; unsigned char he[DUM + 1];
  int L[h14::MAXN + 1];
  const int b = C.fo[q], len = C.fo[q + 1] - b;
  const int* fl = C.fl + b;
  const int* ia = ct.ia + ct.io[q];
  const int* ib = ct.ib + ct.io[q];
  int nl = 0;
  for (int t = 0; t < len; t++) { const int x = fl[t]; pe[t] = s.pend[x]; he[t] = s.held[x]; L[nl] = t; nl += (int)(s.prod[x] ^ 1); }
  pe[DUM] = 1 << 30; he[DUM] = 0;
  int cur = 0, peak = 0;
  auto dec = [&](int a) {
    const int i = a < 0 ? DUM : a;
    const int z = (--pe[i]) == 0;
    cur -= (int)(he[i] & (unsigned char)z);
    he[i] &= (unsigned char)(z ^ 1);
  };
  dec(len > 0 ? len - 1 : -1);
  for (int i = 0; i < nl; i++) {
    const int t = L[i];
    dec(ia[t]);
    dec(ib[t]);
    int mu = 0;
    if (i + 1 < nl) { const int nx = L[i + 1]; mu = (ia[nx] == t) + (ib[nx] == t); }
    const int h = pe[t] > mu;
    he[t] = (unsigned char)h;
    cur += h;
    peak = cur > peak ? cur : peak;
  }
  peakR = peak; endR = cur;
}

#endif
static inline void cap_node(Scratch& scr, const h14::Tab& C, int x, int m) {
  const int A = C.CA[x], op = C.OP[x], B = (op >= 3 && op <= 6) ? C.CB[x] : -1;
  scr.sn[m] = x; scr.sa[m] = A; scr.sb[m] = B;
  scr.reads[A]++; scr.lastread[A] = m;
  if (B >= 0) { scr.reads[B]++; scr.lastread[B] = m; }
  scr.pos[x] = m; scr.chain[x] = 0; scr.slot[x] = -1;
}
static inline void cap_mov(Scratch& scr, const h14::Tab& C, int q, int qi, int m) {
  const int r = C.rr[q];
  scr.sn[m] = -1; scr.sa[m] = r; scr.sb[m] = -1;
  scr.reads[r]++; scr.lastread[r] = m; scr.qof[m] = qi;
}

static inline int commit_cap(h42::OrdScratch& s, const h14::Tab& C, Scratch& scr, int q, int qi, int n, int m) {
  const int NN = n;
  const int* opnd = C.opnd;
  const int r = C.rr[q];
  const int re = opnd[r] < 0 ? r : NN;
  int len = 0;
  for (int t = C.fo[q]; t < C.fo[q + 1]; t++) { const int x = C.fl[t]; s.L[len] = x; len += (int)(s.prod[x] ^ 1); }
  s.L[len] = NN;
  int cur = 0;
  auto dec = [&](int a) {
    s.T[a >> 6] |= 1ull << (a & 63);
    const int z = (--s.pend[a]) == 0;
    cur -= (int)(s.held[a] & (unsigned char)z);
    s.held[a] &= (unsigned char)(z ^ 1);
  };
  dec(re);
  for (int t = 0; t < len; t++) {
    const int x = s.L[t];
    dec(s.EA[x]);
    dec(s.EB[x]);
    s.T[x >> 6] |= 1ull << (x & 63);
    s.prod[x] = 1;
    const int nx = s.L[t + 1];
    const int mu = (s.EA[nx] == x) + (s.EB[nx] == x);
    const int h = s.pend[x] > mu;
    s.held[x] = (unsigned char)h;
    cur += h;
    cap_node(scr, C, x, m + t);
  }
  H720_CNT(cn_commit_len, len); H720_CNT(cn_steps, 1);
  if (len == 0) { cap_mov(scr, C, q, qi, m); return m + 1; }
  return m + len;
}

static int order(const h14::Tab& C, int* out, Scratch& scr) {
  using namespace h42;
  const int n = C.n, Q = C.Q;
  if (Q < 2 || Q > kQ || n < 1 || n > kN) return -1;
  OrdScratch& s = os;
  int m = 0;
  const int W = (n + 63) >> 6;
  for (int x = 0; x <= n; x++) { s.held[x] = 0; s.prod[x] = 0; }
  for (int w = 0; w <= WN; w++) s.T[w] = 0;
#if H720_V & 2
  int rem[kQ]; int nr = Q;
  int inv[kQ]; int ni = Q;
  for (int q = 0; q < Q; q++) { rem[q] = q; inv[q] = q; }
  for (int step = 0; step < Q; step++) {
    if (nr > 1) for (int i = 0; i < ni; i++) { const int q = inv[i]; int p, e; C308_TS(15); trial(s, C, q, p, e); C308_TS(16); s.rp[q] = p; s.re[q] = e; }
    int ba = 0;
    int bk = nr > 1 ? s.rp[rem[0]] * 65536 + s.re[rem[0]] : 0;
    for (int a = 1; a < nr; a++) {
      const int q = rem[a];
      const int k = s.rp[q] * 65536 + s.re[q];
      const bool lt = k < bk;
      bk = lt ? k : bk; ba = lt ? a : ba;
    }
    const int c = rem[ba];
    for (int w = 0; w <= WN; w++) s.T[w] = 0;
    C308_TS(17);
    m = commit_cap(s, C, scr, c, step, n, m);
    C308_TS(18);
    scr.segend[step] = m;
    out[step] = c;
    for (int a = ba; a + 1 < nr; a++) rem[a] = rem[a + 1];
    --nr;
    ni = 0;
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      uint64_t hit = 0;
      for (int w = 0; w < W; w++) hit |= s.fb[q][w] & s.T[w];
      inv[ni] = q; ni += (int)(hit != 0);
    }
  }
#else
  int rem[kQ]; int nr = Q;
  for (int q = 0; q < Q; q++) { rem[q] = q; s.valid[q] = 0; }
  for (int step = 0; step < Q; step++) {
    for (int a = 0; a < nr; a++) {
      const int q = rem[a];
      if (!s.valid[q]) { int p, e; C308_TS(15); trial(s, C, q, p, e); C308_TS(16); s.rp[q] = p; s.re[q] = e; s.valid[q] = 1; }
    }
    int ba = 0;
    int bp = s.rp[rem[0]], be = s.re[rem[0]];
    for (int a = 1; a < nr; a++) {
      const int q = rem[a];
      if (s.rp[q] < bp || (s.rp[q] == bp && s.re[q] < be)) { bp = s.rp[q]; be = s.re[q]; ba = a; }
    }
    const int c = rem[ba];
    for (int w = 0; w <= WN; w++) s.T[w] = 0;
    C308_TS(17);
    m = commit_cap(s, C, scr, c, step, n, m);
    C308_TS(18);
    scr.segend[step] = m;
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
#endif
  return m;
}

template <class E>
static inline void segs_emit2(uint32_t* sg, const Tree& T, const uint32_t* segmap, const E& e) {
  const uint32_t e1022 = e.enc(1022);
  for (int q = 0; q < 18; q++) sg[q] = (q <= T.Q) ? segmap[q] : (uint32_t)T.segs[q];
  for (int q = 0; q < 18; q++) sg[18 + q] = q < T.Q ? e.enc(T.segs[18 + q]) : e1022;
}

struct ScanB {
  uint32_t SL, TT, TWO;
  inline uint32_t enc(uint32_t c) const {
    uint32_t r = (1u << 29) | (c - 64);
    r = c == 2000 ? TWO : r;
    r = c == 901 ? ((3u << 29) | (TT + 48)) : r;
    r = c == 900 ? (2u << 29) : r;
    r = c == 1022 ? ((1u << 31) | (3u << 29) | (TT + 48)) : r;
    r = c < 64 ? c * SL : r;
    return r;
  }
  inline uint32_t dscale(uint32_t d) const { return d * SL; }
};
struct GnB {
  uint32_t CD48, SB, NS, TWO;
  inline uint32_t enc(uint32_t c) const {
    uint32_t r = (1u << 30) | (SB + (c - 64));
    r = c == 2000 ? TWO : r;
    r = c == 901 ? ((1u << 30) | (SB + NS)) : r;
    r = c == 900 ? (1u << 29) : r;
    r = c == 1022 ? (1u << 31) : r;
    r = c < 64 ? c * CD48 : r;
    return r;
  }
  inline uint32_t dscale(uint32_t d) const { return d * CD48; }
};
static inline bool stage_impl(const View& V, Tree& T, cmi6::TreeEnc& te) {
  Scratch& scr = f38::scr;
  const int Q = V.Q, n = V.n();
  if (Q > 17 || n > h14::MAXN || (Q > 1 && n < 1)) { f38::Grid g; return f42::compile_fused_stage(V, T, te, g, 3, nullptr); }
  if (Q < 0) return false;
  T.Q = V.Q; T.nth = V.nth; T.nS = V.nS; T.prog.clear(); T.sp.clear(); T.th.clear();
  if (V.nS > NSV || V.nth > 32) return false;
#if H720_V & 32
  { T.sp.resize((size_t)V.nS); uint16_t* spo = T.sp.data(); int bad = 0;
    for (int i = 0; i < V.nS; i++) { const int o = V.op(i), a = V.a(i), b = V.b(i); bad |= (a > 31) | (b > 31); spo[i] = (uint16_t)(o | (a << 6) | (b << 11)); }
    if (bad) return false; }
#else
  T.sp.reserve((size_t)V.nS);
  for (int i = 0; i < V.nS; i++) { int o = V.op(i), a = V.a(i), b = V.b(i); if (a > 31 || b > 31) return false; T.sp.push_back((uint16_t)(o | (a << 6) | (b << 11))); }
#endif

  C308_TS(4);
  h42::PrepBuf P;
  c720::prepare(V, P);
  C308_TS(5);
  const int* opnd = P.opnd;
  int ord[h14::MAXQ];
  for (int q = 0; q < Q; q++) ord[q] = q;
  h14::Tab& C = h14::tab;
  for (int x = 0; x <= n; x++) { scr.reads[x] = 0; scr.lastread[x] = -1; }
  int m = 0;
  if (Q > 1) {
    m = c720::order(C, ord, scr);
    if (m < 0) { f38::Grid g; return f42::compile_fused_stage(V, T, te, g, 3, nullptr); }
  } else if (Q == 1) {
    for (int t = C.fo[0]; t < C.fo[1]; t++) cap_node(scr, C, C.fl[t], m++);
    if (m == 0) cap_mov(scr, C, 0, 0, m++);
    scr.segend[0] = m;
  }
  C308_TS(7);
#if H720_V & 32
  for (int k = 0; k + 1 < m; k++) {
    const int x = scr.sn[k], xi = x < 0 ? 0 : x;
    const int rd = scr.reads[xi];
    const int r1 = (scr.sa[k + 1] == x) + (scr.sb[k + 1] == x);
    const int c = (int)(x >= 0) & (int)(rd != 0) & (int)(r1 == rd);
    scr.chain[xi] = (unsigned char)(scr.chain[xi] | c);
  }
#else
  for (int k = 0; k + 1 < m; k++) {
    const int x = scr.sn[k];
    if (x < 0 || scr.reads[x] == 0) continue;
    const int r1 = (scr.sa[k + 1] == x) + (scr.sb[k + 1] == x);
    if (r1 == scr.reads[x]) scr.chain[x] = 1;
  }
#endif
  uint16_t capt[18]; for (int q = 0; q < 18; q++) capt[q] = 1022;
  C308_TS(8);
  auto code = [&](int x, int k) -> int {
    if (opnd[x] >= 0) return opnd[x];
    if (scr.chain[x] && scr.pos[x] == k - 1) return 1022;
    return scr.slot[x];
  };
  int nsl = 0, nf = 0, pc = 0;
  bool decFail = false;
  scr.newend[0] = 0;
  for (int k = 0; k < m; k++) {
    const int node = scr.sn[k], A0 = scr.sa[k], B0 = scr.sb[k];
#if H720_V & 16
    const int B0c = B0 < 0 ? 0 : B0;
    auto code2 = [&](int x) -> int { const int o = opnd[x]; const int cs = (int)(scr.chain[x] != 0) & (int)(scr.pos[x] == k - 1); const int sl = scr.slot[x]; return o >= 0 ? o : (cs ? 1022 : sl); };
    const int A = code2(A0);
    int B = code2(B0c); B = B0 >= 0 ? B : 0;
    if ((A | B) < 0) return false;
    { const int sA = scr.slot[A0]; const int cA = (int)(opnd[A0] < 0) & (int)(sA >= 0) & (int)(scr.lastread[A0] == k); scr.freeS[nf] = sA; nf += cA; scr.slot[A0] = cA ? -2 : sA; }
    { const int sB = scr.slot[B0c]; const int cB = (int)(B0 >= 0) & (int)(B0 != A0) & (int)(opnd[B0c] < 0) & (int)(sB >= 0) & (int)(scr.lastread[B0c] == k); scr.freeS[nf] = sB; nf += cB; scr.slot[B0c] = cB ? -2 : sB; }
    int d;
    { const int nd = node < 0 ? 0 : node;
      const int cN = (int)(node >= 0) & (int)(scr.reads[nd] > 0) & (int)(scr.chain[nd] == 0);
      const int has = nf > 0;
      const int pop = scr.freeS[has ? nf - 1 : 0];
      d = cN ? (has ? pop : nsl) : 63;
      nf -= cN & has; nsl += cN & (has ^ 1);
      scr.slot[nd] = cN ? d : scr.slot[nd]; }
#else
    const int A = code(A0, k), B = B0 >= 0 ? code(B0, k) : 0;
    if (A < 0 || B < 0) return false;
    if (opnd[A0] < 0 && scr.slot[A0] >= 0 && scr.lastread[A0] == k) { scr.freeS[nf++] = scr.slot[A0]; scr.slot[A0] = -2; }
    if (B0 >= 0 && B0 != A0 && opnd[B0] < 0 && scr.slot[B0] >= 0 && scr.lastread[B0] == k) { scr.freeS[nf++] = scr.slot[B0]; scr.slot[B0] = -2; }
    int d = 63;
    if (node >= 0 && scr.reads[node] > 0 && !scr.chain[node]) {
      d = nf > 0 ? scr.freeS[--nf] : nsl++;
      scr.slot[node] = d;
    }
#endif
    if (node < 0) { capt[scr.qof[k]] = (uint16_t)A; scr.newend[k + 1] = pc; continue; }
    const uint32_t op = (uint32_t)C.OP[node];
    scr.pk[pc] = op | (uint32_t)d << 6 | (uint32_t)A << 12 | (uint32_t)B << 22;
    if (!rec_fields(op, (uint32_t)d, (uint32_t)A, (uint32_t)B, scr.rec[pc])) decFail = true;
    ++pc;
    scr.newend[k + 1] = pc;
  }
  if (nsl > MAXSLOT) return false;
  C308_TS(9); H720_CNT(cn_pc, pc);
  T.prog.assign(scr.pk, scr.pk + pc);
#if H720_V & 32
  { T.segs[0] = 0;
    for (int q = 1; q < 18; q++) T.segs[q] = q <= Q ? (uint16_t)scr.newend[scr.segend[q - 1 < Q ? q - 1 : 0]] : (uint16_t)0;
    for (int q = 0; q < 18; q++) T.segs[18 + q] = capt[q]; }
#else
  for (int q = 0; q < 36; q++) T.segs[q] = 0;
  for (int q = 0; q < Q; q++) T.segs[q + 1] = (uint16_t)scr.newend[scr.segend[q]];
  for (int q = 0; q < 18; q++) T.segs[18 + q] = 1022;
  for (int q = 0; q < Q; q++) T.segs[18 + q] = capt[q];
#endif
  T.nslot = nsl;
#if H720_V & 32
  T.th.resize((size_t)V.nth); if (V.nth) std::memcpy(T.th.data(), V.theta, sizeof(float) * (size_t)V.nth);
#else
  T.th.reserve((size_t)V.nth);
  for (int s = 0; s < V.nth; s++) T.th.push_back(V.th(s));
#endif
  C308_TS(10);
  if (decFail) { f38::Grid g; return f42::compile_fused_stage(V, T, te, g, 3, nullptr); }

  const int SS = T.nth == 0 ? 2 : 16;
  #if H720_V & 1
  const enc34::ScanE es0 = enc34::make_scan(T, SS); const ScanB es{es0.SL, es0.TT, es0.TWO};
#else
  const enc34::ScanE es = enc34::make_scan(T, SS);
#endif
  const bool gin = T.nth >= 1 && T.nth <= GN_TH_MAX;
  const int CD = T.nth <= 4 ? 1 + (T.nth + 1) / 2 : 3;
  const uint32_t NSC = (uint32_t)T.nS + 2;
  const uint32_t SB = (uint32_t)T.nslot * (uint32_t)CD * 48;
  const bool gok = gin && SB < 65536;
  #if H720_V & 1
  const enc34::GnE eg0 = enc34::make_gn(T, CD); const GnB eg{eg0.CD48, eg0.SB, eg0.NS, eg0.TWO};
#else
  const enc34::GnE eg = enc34::make_gn(T, CD);
#endif
  const enc34::Dec* R = scr.rec;
  uint32_t* os = scr.ob; uint32_t* og = scr.ob2;
  size_t k = 0, out = 0; int qp = 0;
  const int Qt = T.Q;
  while (k < (size_t)pc) {
    while (qp <= Qt && (size_t)T.segs[qp] == k) scr.segmap[qp++] = (uint32_t)out;
    const enc34::Dec& a = R[k];
    bool fuse = false;
    if (k + 1 < (size_t)pc && a.zo == 0 && a.d63) {
      const enc34::Dec& b = R[k + 1];
      fuse = b.zo >= 7 && b.zo <= 15 && b.A == 1022 && !(qp <= Qt && (size_t)T.segs[qp] == k + 1);
    }
    const enc34::Dec& z = fuse ? R[k + 1] : a;
    if (cmi6::g_b1) {
      uint32_t* s2 = os + 2 * out;
      s2[0] = (uint32_t)z.zo | ((uint32_t)z.d << 8); s2[1] = h313b::map_code(a.A, (uint32_t)T.nS) | (h313b::map_code(a.B, (uint32_t)T.nS) << 16);
#ifdef H313B_NEG
      if (a.B == 64) s2[1] = h313b::map_code(a.A, (uint32_t)T.nS) | (65u << 16);
#endif
      k += fuse ? 2 : 1;
      ++out;
      continue;
    }
    uint32_t* s = os + 4 * out;
    s[0] = z.zo; s[1] = z.d63 ? 0xFFFFFFFFu : es.dscale(z.d); s[2] = es.enc(a.A); s[3] = es.enc(a.B);
    if (gok) {
      uint32_t* gg = og + 4 * out;
      gg[0] = z.zo; gg[1] = z.d63 ? 0xFFFFFFFFu : eg.dscale(z.d); gg[2] = eg.enc(a.A); gg[3] = eg.enc(a.B);
    }
    k += fuse ? 2 : 1;
    ++out;
  }
  while (qp <= Qt && (size_t)T.segs[qp] == (size_t)pc) scr.segmap[qp++] = (uint32_t)out;
  if (qp != Qt + 1) { f38::Grid g; return f42::compile_fused_stage(V, T, te, g, 3, nullptr); }
  C308_TS(11);
  if (!gok) std::memset(te.gg, 0, sizeof te.gg);
  if (cmi6::g_b1) { if (cmi6::g_a1) { te.tkp = os; te.ntkp = (uint32_t)(2 * out); } else te.tk.assign(os, os + 2 * out); }
  else if (cmi6::g_a1) { te.sp4p = os; te.nsp4p = (uint32_t)(4 * out); } else te.sp4.assign(os, os + 4 * out);
  segs_emit2(te.sg, T, scr.segmap, es);
  te.s_ok = 1;
  if (gok) {
    if (cmi6::g_b1) {} else if (cmi6::g_a1) { te.gp4p = og; te.ngp4p = (uint32_t)(4 * out); } else te.gp4.assign(og, og + 4 * out);
    segs_emit2(te.gg, T, scr.segmap, eg);
    te.ghw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
    te.gper = (uint32_t)((size_t)SB + 2 * (size_t)CD * NSC);
    te.g_ok = 1;
  }
  C308_TS(12);
  return true;
}

static inline void fix_b1(const Tree& T, cmi6::TreeEnc& te) {
  if (!cmi6::g_b1 || (te.sp4.empty() && te.gp4.empty())) return;
  bool ok = true;
  h313b::g_fix.fetch_add(1, std::memory_order_relaxed);
  te.tk.clear(); te.tkp = nullptr; te.ntkp = 0;
  if (te.s_ok) {
    ok = h313b::contract_scan(te.sp4.data(), te.sp4.size() / 4, h313b::scan_p(T.nth, T.nslot, T.nS), te.tk);
    if (ok && te.g_ok) {
      std::vector<uint32_t> g(te.tk.size() * 2);
      h313b::expand_gn(te.tk.data(), te.tk.size() / 2, h313b::gn_p(te.ghw), g.data());
      ok = g.size() == te.gp4.size() && std::memcmp(g.data(), te.gp4.data(), g.size() * 4) == 0;
    }
  } else if (te.g_ok) ok = h313b::contract_gn(te.gp4.data(), te.gp4.size() / 4, h313b::gn_p(te.ghw), te.tk);
  if (!ok) { fprintf(stderr, "H313b: fallback tree not representable in the shared stream (Q %d nth %d nS %d nslot %d)\n", T.Q, T.nth, T.nS, T.nslot); abort(); }
  te.sp4.clear(); te.gp4.clear();
}
static inline bool stage(const View& V, Tree& T, cmi6::TreeEnc& te) {
  const bool r = stage_impl(V, T, te);
  if (r) fix_b1(T, te);
  return r;
}

static inline uint32_t publish_idx(cc14::Handle& H, uint64_t h, cc14::Ent X, const uint32_t* kw, size_t nk, const uint32_t* pw, size_t np, const uint16_t* spw, size_t ns) {
  using cc14::Shared; cc14::Shared& S = H.S;
  uint32_t idx, koff, poff, soff;
  C308_TS(13);
  if (!cc14::take_block(S.nent, H.r_ent_next, H.r_ent_end, 1, cc14::RSV_ENT, Shared::MAXENT, idx)) { ++H.drops; return 0; }
  if (!cc14::take_block(S.nkey, H.r_key_next, H.r_key_end, (uint32_t)nk, cc14::RSV_WORDS, Shared::MAXKEY, koff)) { ++H.drops; return 0; }
  if (!cc14::take_block(S.nprog, H.r_prog_next, H.r_prog_end, (uint32_t)np, cc14::RSV_WORDS, Shared::MAXPROG, poff)) { ++H.drops; return 0; }
  if (!cc14::take_block(S.nsp, H.r_sp_next, H.r_sp_end, (uint32_t)ns, cc14::RSV_WORDS, Shared::MAXSP, soff)) { ++H.drops; return 0; }
  C308_TS(19);
  X.touch = h156::g_touch; X.key_off = koff; X.nkey = (uint32_t)nk; X.prog_off = poff; X.nprog = (uint32_t)np; X.sp_off = soff; X.nsp = (uint32_t)ns;
  cc14::Ent& D = S.ent[idx];
  D = X;
  C308_TS(20);
#ifdef H720_DIAG3
  if (0)
#endif
  { std::memcpy(S.key + koff, kw, sizeof(uint32_t) * nk);
  if (np) std::memcpy(S.prog + poff, pw, sizeof(uint32_t) * np);
  if (ns) std::memcpy(S.sp + soff, spw, sizeof(uint16_t) * ns); }
  C308_TS(21);
  uint32_t p = (uint32_t)(h & (Shared::TB - 1));
  for (;;) {
    uint32_t expected = 0;
#ifdef H720_DIAG
    if (H720_DIAG == 1) { if (S.tab[p].compare_exchange_strong(expected, idx + 1, std::memory_order_relaxed, std::memory_order_relaxed)) { C308_TS(22); return idx + 1; } }
    else if (S.tab[p].load(std::memory_order_relaxed) == 0) { S.tab[p].store(idx + 1, std::memory_order_relaxed); expected = 0; C308_TS(22); return idx + 1; }
    else expected = S.tab[p].load(std::memory_order_relaxed);
    if (0)
#endif
#if H720_ATOM & 1
    if (S.tab[p].compare_exchange_strong(expected, idx + 1, std::memory_order_release, std::memory_order_relaxed)) return idx + 1;
    std::atomic_thread_fence(std::memory_order_acquire);
#else
    if (S.tab[p].compare_exchange_strong(expected, idx + 1, std::memory_order_release, std::memory_order_acquire)) { C308_TS(22); return idx + 1; }
#endif
    const cc14::Ent& E = S.ent[expected - 1];
    if (E.h == h && E.nkey == nk && std::memcmp(S.key + E.key_off, kw, sizeof(uint32_t) * nk) == 0) return expected;
    p = (p + 1) & (Shared::TB - 1);
  }
}

static inline void tok_to_vec(cmi6::TreeEnc& te) {
  if (te.tkp) { te.tk.assign(te.tkp, te.tkp + te.ntkp); te.tkp = nullptr; te.ntkp = 0; }
  if (te.sp4p) { te.sp4.assign(te.sp4p, te.sp4p + te.nsp4p); te.sp4p = nullptr; te.nsp4p = 0; }
  if (te.gp4p) { te.gp4.assign(te.gp4p, te.gp4p + te.ngp4p); te.gp4p = nullptr; te.ngp4p = 0; }
}

static inline bool compile_enc(cmi6::Handle& H, const View& V, Tree& T, cmi6::TreeEnc& te, double& enc_ms, const uint64_t* hpre = nullptr) {
  const size_t nk = (size_t)(V.theta - V.w);
  const uint64_t h = hpre ? *hpre : cc14::hkey(V.w, nk);
  bool hit = false; uint32_t hidx = 0;
  uint32_t p = (uint32_t)(h & (cc14::Shared::TB - 1));
  for (;;) {
    const uint32_t v = H.C.S.tab[p].load(std::memory_order_acquire);
    if (v == 0) break;
    const cc14::Ent& E = H.C.S.ent[v - 1];
    if (E.h == h && E.nkey == nk && std::memcmp(H.C.S.key + E.key_off, V.w, sizeof(uint32_t) * nk) == 0) {
      ++H.C.hits;
      if (h156::g_on) h156::touch_u32(&H.C.S.ent[v - 1].touch, h156::g_touch);
      if (!E.ok) { ++H.C.neg_hits; return false; }
      C308_TS(0);
      T.Q = E.Q; T.nth = E.nth; T.nS = E.nS; T.nslot = E.nslot;
      if (g_a2) { T.progp = H.C.S.prog + E.prog_off; T.nprogp = E.nprog; T.segsp = E.segs;
#ifdef H313A_NEG
        if (E.nprog) T.nprogp = E.nprog - 1;
#endif
      }
      else {
        T.prog.assign(H.C.S.prog + E.prog_off, H.C.S.prog + E.prog_off + E.nprog);
        std::memcpy(T.segs, E.segs, sizeof T.segs);
      }
      T.sp.assign(H.C.S.sp + E.sp_off, H.C.S.sp + E.sp_off + E.nsp);
      T.th.resize((size_t)V.nth);
      for (int s = 0; s < V.nth; s++) T.th[s] = V.th(s);
      hit = true; hidx = v;
      C308_TS(1);
      break;
    }
    p = (p + 1) & (cc14::Shared::TB - 1);
  }
  if (!hit) {
    C308_TS(0);
    ++H.C.misses;
#ifndef H720_PF
#define H720_PF 0
#endif
    if (H720_PF & 1) __builtin_prefetch(&H.C.S.tab[p], 1, 3);
    if (H720_PF & 2) {
      const uint32_t* e = (const uint32_t*)&H.C.S.ent[H.C.r_ent_next]; __builtin_prefetch(e, 1, 3); __builtin_prefetch(e + 16, 1, 3); __builtin_prefetch(e + 32, 1, 3);
      __builtin_prefetch(H.C.S.key + H.C.r_key_next, 1, 3); __builtin_prefetch(H.C.S.key + H.C.r_key_next + 16, 1, 3);
      __builtin_prefetch(H.C.S.prog + H.C.r_prog_next, 1, 3); __builtin_prefetch(H.C.S.sp + H.C.r_sp_next, 1, 3);
    }
    if (H720_PF & 4) {
      const char* e = (const char*)&H.C.S.ent[H.C.r_ent_next]; __builtin_prefetch(e, 1, 3); __builtin_prefetch(e + 64, 1, 3); __builtin_prefetch(e + 127, 1, 3);
      const char* k = (const char*)(H.C.S.key + H.C.r_key_next); for (size_t o = 0; o < nk * 4; o += 64) __builtin_prefetch(k + o, 1, 3); __builtin_prefetch(k + nk * 4 - 1, 1, 3);
      const char* pg = (const char*)(H.C.S.prog + H.C.r_prog_next); const size_t pb = (size_t)V.n() * 4; for (size_t o = 0; o < pb; o += 64) __builtin_prefetch(pg + o, 1, 3);
      const char* sp = (const char*)(H.C.S.sp + H.C.r_sp_next); __builtin_prefetch(sp, 1, 3); __builtin_prefetch(sp + (size_t)V.nS * 2, 1, 3);
    }
    const bool ok = stage(V, T, te);
    cc14::Ent X; std::memset(&X, 0, sizeof X);
    X.h = h; X.ok = ok ? 1 : 0;
    if (ok) {
      X.Q = T.Q; X.nth = T.nth; X.nS = T.nS; X.nslot = T.nslot;
      std::memcpy(X.segs, T.segs, sizeof T.segs);
      hidx = publish_idx(H.C, h, X, V.w, nk, T.prog.data(), T.prog.size(), T.sp.data(), T.sp.size());
    } else {
      H.C.publish(h, X, V.w, nk, nullptr, 0, nullptr, 0);
      return false;
    }
    C308_TS(13);
  }
  const bool inscr = !hit && (cmi6::g_b1 ? te.tkp != nullptr : te.sp4p != nullptr);
  uint32_t idx = hidx - 1;
  if (hidx == 0 && !cmi6::locate(H.C.S, V.w, nk, h, idx)) {
    if (hit) { h313a::own_prog(T); enc34::encode_tree_h34(T, te); fix_b1(T, te); }
    if (inscr) tok_to_vec(te);
    H.enc_local++;
    C308_TS(hit ? 2 : 14);
    return true;
  }
  cmi6::Slot& s = H.M.at(idx);
  const uint32_t st = __atomic_load_n(&s.st, __ATOMIC_ACQUIRE);
  if (st == 2) { cmi6::load_enc(s, H.M.arena, te); H.enc_hit++; C308_TS(hit ? 2 : 14); return true; }
  if (hit) { h313a::own_prog(T); enc34::encode_tree_h34(T, te); fix_b1(T, te); }
  H.enc_miss++;
  uint32_t zero = 0;
  if (st == 0 && __atomic_compare_exchange_n(&s.st, &zero, 1u, false, (H720_ATOM & 2) ? __ATOMIC_RELAXED : __ATOMIC_ACQ_REL, (H720_ATOM & 2) ? __ATOMIC_RELAXED : __ATOMIC_ACQUIRE)) {
    const uint64_t n4 = cmi6::g_b1 ? te.ntkw() : te.nsp(), ng = cmi6::g_b1 ? 0 : te.ngp();
    const bool gt856 = cmi6::g_gt && cmi6::g_gtfn && cmi6::g_a1 && !cmi6::g_b1;
    static thread_local cmi6::GtScr gts856;
    const uint64_t ngt = gt856 ? cmi6::gt_make(T, gts856) : 0;
    uint8_t cm913 = 0; uint32_t tw913 = 0;
    if (cmi6::g_h913) { cm913 = cmi6::h913_cm(te, T.Q); tw913 = cmi6::tabw(cm913, (uint8_t)T.Q); }
    const uint64_t need = n4 + ng + ngt + tw913;
    const uint64_t off = H.M.bump.fetch_add(need, (H720_ATOM & 4) ? std::memory_order_relaxed : std::memory_order_seq_cst);
    if (off + need <= cmi6::ARENA_WORDS) {
      if (n4) std::memcpy(H.M.arena + off, cmi6::g_b1 ? te.ptk() : te.psp(), n4 * sizeof(uint32_t));
      if (ng) std::memcpy(H.M.arena + off + n4, te.pgp(), ng * sizeof(uint32_t));
      if (gt856) cmi6::gt_put(H.M.arena + off + n4 + ng, gts856);
      s.s_ok = te.s_ok; s.g_ok = te.g_ok; s.ghw = te.ghw; s.gper = te.gper;
      uint8_t cm = 0;
      if (cmi6::g_a3 && (unsigned)T.Q <= 17u) { if (cmi6::tab_canon(te.sg, T.Q)) cm |= 1; if (cmi6::tab_canon(te.gg, T.Q)) cm |= 2; }
      if (cmi6::g_h913) { cm = cm913; cmi6::h913_put(H.M.arena + off + n4 + ng + ngt, te, cm, T.Q); }
      else {
      if (cm & 1) cmi6::tab_pack(s.sg, te.sg, T.Q); else std::memcpy(s.sg, te.sg, sizeof s.sg);
      if (cm & 2) cmi6::tab_pack(s.gg, te.gg, T.Q); else std::memcpy(s.gg, te.gg, sizeof s.gg);
      }
      if (gt856) cm |= 4;
      s.cmode = cm; s.Qs = (uint8_t)T.Q;
      s.sp_off = (uint32_t)off; s.nsp4 = (uint32_t)n4; s.gp_off = (uint32_t)(off + n4); s.ngp4 = (uint32_t)ng;
      H.M.nst.fetch_add(1, std::memory_order_relaxed);
      __atomic_store_n(&s.st, 2u, __ATOMIC_RELEASE);
      if (gt856) te.gtp = H.M.arena + off + n4 + ng;
      if (inscr) { if (cmi6::g_b1) te.tkp = H.M.arena + off; else { te.sp4p = H.M.arena + off; if (te.gp4p) te.gp4p = H.M.arena + off + n4; } }
    } else {
      __atomic_store_n(&s.st, 3u, __ATOMIC_RELEASE);
      if (inscr) tok_to_vec(te);
    }
  } else if (inscr) tok_to_vec(te);
  C308_TS(hit ? 2 : 14);
  return true;
}
}
