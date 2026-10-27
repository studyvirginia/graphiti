#pragma once
#include "../s02_nf/sk_nf_h2b.hpp"
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cstdio>

namespace nfc {
using namespace skn;

#ifdef NFC_COUNT
struct Counters { uint64_t v[96]; };
inline int g_mm_n = 0; inline uint64_t g_mm[4096]; inline int g_tb_n = 0; inline uint64_t g_tb[4096];
inline Counters g_cnt;
#define CNT(i, k) (nfc::g_cnt.v[i] += (k))
enum { C_TE_CALL, C_TE_PROBE, C_TE_HIT, C_MON_CALL, C_MON_PROBE, C_MON_HIT, C_MON_KEYLEN, C_FAC_CALL, C_FAC_PROBE, C_FAC_HIT,
       C_SUM_CALL, C_SUM_PROBE, C_SUM_HIT, C_SORTFAC_CALL, C_SORTFAC_N, C_SORTFAC_MOVES, C_MONMUL, C_TEBIN, C_TEBIN_EARLY,
       C_NORM_CALL, C_NORM_N, C_NORM_SORTED, C_NORM_N_LE5, C_NORM_N_LE10, C_NORM_N_LE20, C_NORM_N_LE30, C_NORM_N_GT30,
       C_NORM_EQ2, C_NORM_EQ3P, C_NORM_MERGES, C_TREES, C_TABLE_CLEARS, C_TEVAL_NODES, C_TEVAL_LIBM, C_NTH, C_NI, C_EMIT_MONO, C_EMIT_FAC, C_EMIT_SUM, C_RUN3_FALLBACK, C_N24, C_HASHWORDS, C_MM_REPEAT, C_MM_IDENT, C_TB_REPEAT, C_CUBE_NV1, C_CUBE_NV2, C_CUBE_NV3, C_CUBE_BIG, C_SQ_NV2, C_SQ_NV3, C_SQ_NV4, C_SQ_NV5, C_SQ_BIG, C_NTE, C_NMON, C_NFAC, C_NSUM, C_NVP, C_INSTR };
#else
#define CNT(i, k) ((void)0)
#endif

struct KE { uint64_t h; int32_t coef, mono; };

struct Ctx {
  TE te[MAXTE]; int nte;
  KE kbuf[MAXP2];
  FAC fac[MAXFAC]; int nfac;
  MON mon[MAXMONO]; int nmon;
  int32_t mp[MAXMP]; uint64_t mph[MAXMP]; int nmp;
  SUMO sum[MAXSUM]; int nsum;
  PR sp[MAXSP]; int nsp;
  PR vp[MAXPAIR]; int nvp;
  uint64_t tte[TABN], tfac[TABN], tmon[TABN], tsum[TABN];
  uint32_t gen;
  bool ovf;
  int32_t scr[MAXSCR]; uint64_t shbuf[MAXSCR];
  PR scrA[MAXP2], scrB[MAXP2], scrC[MAXP2];
  struct Ins { uint32_t mask; uint16_t a, b; uint8_t op; };
  Ins is_[MAXINS], iv_[MAXINS];
  int ni, nsS, nvS, nth;
  bool have_one; uint16_t one_ref;
  int32_t slot_te[MAXTH]; uint8_t slot_freq[MAXTH];
  uint32_t freqTE[MAXTE];
  double tevals[MAXTE];
  uint32_t tgen[MAXTE]; uint32_t vgen = 0;
  uint32_t mm[MAXMONO], mf[MAXFAC], ms[MAXSUM], mt[MAXTE]; uint32_t mbase = 0;
  uint64_t binh[16], unh[16], facbh[8][16];
  Ctx() { memset(tgen, 0, sizeof tgen); memset(mm, 0, sizeof mm); memset(mf, 0, sizeof mf); memset(ms, 0, sizeof ms); memset(mt, 0, sizeof mt); memset(freqTE, 0, sizeof freqTE); for (int o = 0; o < 16; ++o) { binh[o] = mx(H_BIN, o); unh[o] = mx(H_UN, o); }
          const uint64_t bs[5] = {H_T, H_FN, H_RS, H_RM, H_COLL}; for (int k = 0; k < 5; ++k) for (int o = 0; o < 16; ++o) facbh[k][o] = mx(bs[k], o);
          gen = 0; memset(tte, 0, sizeof tte); memset(tfac, 0, sizeof tfac);
          memset(tmon, 0, sizeof tmon); memset(tsum, 0, sizeof tsum); }
  void bump_gen() {
    if (++gen == 0) { memset(tte, 0, sizeof tte); memset(tfac, 0, sizeof tfac);
                      memset(tmon, 0, sizeof tmon); memset(tsum, 0, sizeof tsum); gen = 1; }
  }

  int te_leaf(int pos) {
    if (nte >= MAXTE) { ovf = true; return COEF_LIN; }
    te[nte] = TE{TE_LEAF, 0, -1, -1, pos, H_LEAF};
    return nte++;
  }
  int te_intern(uint8_t kind, uint8_t op, int32_t a, int32_t b, uint64_t h) {
    uint32_t p = (uint32_t)(h & (TABN - 1));
    CNT(C_TE_CALL, 1);
    for (;;) {
      uint64_t e = tte[p];
      if (!cur(e, gen)) break;
      const TE& x = te[(uint32_t)e - 1];
      if (x.kind == kind && x.op == op && x.a == a && x.b == b) { CNT(C_TE_HIT, 1); return (int)(uint32_t)e - 1; }
      CNT(C_TE_PROBE, 1);
      p = (p + 1) & (TABN - 1);
    }
    if (nte >= MAXTE) { ovf = true; return COEF_LIN; }
    te[nte] = TE{kind, op, a, b, -1, h};
    tte[p] = ((uint64_t)gen << 32) | (uint32_t)(nte + 1);
    return nte++;
  }
  int te_bin(uint8_t op, int a, int b) {
    CNT(C_TEBIN, 1);
#ifdef NFC_COUNT
    { uint64_t key = ((uint64_t)op << 56) | ((uint64_t)(op == OP_ADD || op == OP_MUL ? std::min(a, b) : a) << 28) | (uint64_t)(op == OP_ADD || op == OP_MUL ? std::max(a, b) : b); bool f = false;
      for (int i = 0; i < g_tb_n; i++) if (g_tb[i] == key) { f = true; break; }
      if (f) CNT(C_TB_REPEAT, 1); else if (g_tb_n < 4096) g_tb[g_tb_n++] = key; }
#endif
    if (op == OP_MUL) { if (a == COEF_FIX) { CNT(C_TEBIN_EARLY, 1); return b; } if (b == COEF_FIX) { CNT(C_TEBIN_EARLY, 1); return a; } }
    if (op == OP_DIV && b == COEF_FIX) { CNT(C_TEBIN_EARLY, 1); return a; }
    if ((op == OP_ADD || op == OP_MUL) && (te[b].h < te[a].h || (te[b].h == te[a].h && b < a))) { int t = a; a = b; b = t; }
    uint64_t h = mx(mx(binh[op], te[a].h), te[b].h);
    return te_intern(TE_BIN, op, a, b, h);
  }
  int te_un(uint8_t op, int a) { return te_intern(TE_UN, op, a, -1, mx(unh[op], te[a].h)); }
  int te_neg(int a) { return te_intern(TE_NEG, 0, a, -1, mx(H_NEG, te[a].h)); }

  void sort_fac(int32_t* f, int n) const {
    CNT(C_SORTFAC_CALL, 1); CNT(C_SORTFAC_N, n);
    for (int i = 1; i < n; ++i) {
      int32_t x = f[i]; uint64_t hx = fac[x].h; int j = i - 1;
      while (j >= 0 && (fac[f[j]].h > hx || (fac[f[j]].h == hx && f[j] > x))) { CNT(C_SORTFAC_MOVES, 1); f[j + 1] = f[j]; --j; }
      f[j + 1] = x;
    }
  }
  int mon_intern_h(const int32_t* f, const uint64_t* fh, int n, uint64_t h) {
    uint32_t p = (uint32_t)(h & (TABN - 1));
    CNT(C_MON_CALL, 1); CNT(C_MON_KEYLEN, n);
    for (;;) {
      uint64_t e = tmon[p];
      if (!cur(e, gen)) break;
      const MON& M = mon[(uint32_t)e - 1];
      if (M.len == n && M.h == h && eq32(mp + M.off, f, n)) { CNT(C_MON_HIT, 1); return (int)(uint32_t)e - 1; }
      CNT(C_MON_PROBE, 1);
      p = (p + 1) & (TABN - 1);
    }
    if (nmon >= MAXMONO || nmp + n > MAXMP) { ovf = true; return 0; }
    for (int k = 0; k < n; ++k) { mp[nmp + k] = f[k]; mph[nmp + k] = fh[k]; }
    mon[nmon] = MON{nmp, n, h};
    nmp += n;
    tmon[p] = ((uint64_t)gen << 32) | (uint32_t)(nmon + 1);
    return nmon++;
  }
  int mon_single(int f) {
    int32_t x = f; uint64_t fh = fac[f].h;
    return mon_intern_h(&x, &fh, 1, mx(H_EMPTY, fh));
  }
  int mon_mul(int m1, int m2) {
    CNT(C_MONMUL, 1);
#ifdef NFC_COUNT
    { if (mon[m1].len == 0 || mon[m2].len == 0) CNT(C_MM_IDENT, 1);
      uint64_t key = m1 < m2 ? ((uint64_t)m1 << 32 | m2) : ((uint64_t)m2 << 32 | m1); bool f = false;
      for (int i = 0; i < g_mm_n; i++) if (g_mm[i] == key) { f = true; break; }
      if (f) CNT(C_MM_REPEAT, 1); else if (g_mm_n < 4096) g_mm[g_mm_n++] = key; }
#endif
    const int n1 = mon[m1].len, n2 = mon[m2].len;
    if (n1 + n2 > MAXSCR) { ovf = true; return 0; }
    if (n1 == 0) return m2;
    if (n2 == 0) return m1;
    const int32_t* a = mp + mon[m1].off; const uint64_t* ah = mph + mon[m1].off;
    const int32_t* b = mp + mon[m2].off; const uint64_t* bh = mph + mon[m2].off;
    uint64_t* sh = shbuf; uint64_t h = H_EMPTY;
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
      bool ta = ah[i] < bh[j] || (ah[i] == bh[j] && a[i] <= b[j]);
      if (ta) { scr[k] = a[i]; sh[k] = ah[i]; ++i; } else { scr[k] = b[j]; sh[k] = bh[j]; ++j; }
      h = mx(h, sh[k]); ++k;
    }
    while (i < n1) { scr[k] = a[i]; sh[k] = ah[i]; ++i; h = mx(h, sh[k]); ++k; }
    while (j < n2) { scr[k] = b[j]; sh[k] = bh[j]; ++j; h = mx(h, sh[k]); ++k; }
    return mon_intern_h(scr, sh, k, h);
  }
  int fac_intern(uint8_t kind, uint8_t op, int32_t arg, uint64_t argh) {
    uint64_t h = (kind == FAC_T) ? H_T : mx(facbh[kind][op], argh);
    uint32_t p = (uint32_t)(h & (TABN - 1));
    CNT(C_FAC_CALL, 1);
    for (;;) {
      uint64_t e = tfac[p];
      if (!cur(e, gen)) break;
      const FAC& x = fac[(uint32_t)e - 1];
      if (x.kind == kind && x.op == op && x.arg == arg) { CNT(C_FAC_HIT, 1); return (int)(uint32_t)e - 1; }
      CNT(C_FAC_PROBE, 1);
      p = (p + 1) & (TABN - 1);
    }
    if (nfac >= MAXFAC) { ovf = true; return 0; }
    fac[nfac] = FAC{kind, op, arg, h};
    tfac[p] = ((uint64_t)gen << 32) | (uint32_t)(nfac + 1);
    return nfac++;
  }
  int sum_intern(const PR* pr, int n) {
    uint64_t h = H_SUM;
    for (int k = 0; k < n; ++k) h = mx(mx(h, te[pr[k].coef].h), mon[pr[k].mono].h);
    uint32_t p = (uint32_t)(h & (TABN - 1));
    CNT(C_SUM_CALL, 1);
    for (;;) {
      uint64_t e = tsum[p];
      if (!cur(e, gen)) break;
      const SUMO& S = sum[(uint32_t)e - 1];
      if (S.len == n && S.h == h && eqpr(sp + S.off, pr, n)) { CNT(C_SUM_HIT, 1); return (int)(uint32_t)e - 1; }
      CNT(C_SUM_PROBE, 1);
      p = (p + 1) & (TABN - 1);
    }
    if (nsum >= MAXSUM || nsp + n > MAXSP) { ovf = true; return 0; }
    for (int k = 0; k < n; ++k) sp[nsp + k] = pr[k];
    sum[nsum] = SUMO{nsp, n, h};
    nsp += n;
    tsum[p] = ((uint64_t)gen << 32) | (uint32_t)(nsum + 1);
    return nsum++;
  }
  int norm(PR* p, int n, bool top, bool run3 = false) {
    if (n <= 1) return n;
    CNT(C_NORM_CALL, 1); CNT(C_NORM_N, n);
    CNT(n <= 5 ? C_NORM_N_LE5 : n <= 10 ? C_NORM_N_LE10 : n <= 20 ? C_NORM_N_LE20 : n <= 30 ? C_NORM_N_LE30 : C_NORM_N_GT30, 1);
    if (n < 24 || (n <= 64 && !run3)) {
      kbuf[0].h = mon[p[0].mono].h; kbuf[0].coef = p[0].coef; kbuf[0].mono = p[0].mono;
      for (int i = 1; i < n; ++i) {
        KE x; x.h = mon[p[i].mono].h; x.coef = p[i].coef; x.mono = p[i].mono;
        int j = i;
        while (j > 0 && (x.h < kbuf[j - 1].h || (x.h == kbuf[j - 1].h && x.mono < kbuf[j - 1].mono))) { kbuf[j] = kbuf[j - 1]; --j; }
        kbuf[j] = x;
      }
      if (n >= 24) {
        bool r3 = false;
        for (int k = 2; k < n; ++k) if (kbuf[k].mono == kbuf[k - 2].mono) { r3 = true; break; }
        if (r3) {
          for (int k = 0; k < n; ++k) { kbuf[k].h = mon[p[k].mono].h; kbuf[k].coef = p[k].coef; kbuf[k].mono = p[k].mono; }
          std::sort(kbuf, kbuf + n, [](const KE& a, const KE& b) {
            return a.h < b.h || (a.h == b.h && a.mono < b.mono);
          });
        }
      }
    } else if (n <= MAXP2) {
      for (int k = 0; k < n; ++k) { kbuf[k].h = mon[p[k].mono].h; kbuf[k].coef = p[k].coef; kbuf[k].mono = p[k].mono; }
      std::sort(kbuf, kbuf + n, [](const KE& a, const KE& b) {
        return a.h < b.h || (a.h == b.h && a.mono < b.mono);
      });
    } else {
      std::sort(p, p + n, [this](const PR& a, const PR& b) {
        uint64_t ha = mon[a.mono].h, hb = mon[b.mono].h;
        return ha < hb || (ha == hb && a.mono < b.mono);
      });
      for (int k = 0; k < n; ++k) { kbuf[k].coef = p[k].coef; kbuf[k].mono = p[k].mono; }
    }
#ifdef NFC_COUNT
    { int run = 1; for (int k = 1; k <= n; ++k) { if (k < n && kbuf[k].mono == kbuf[k - 1].mono) ++run; else { if (run == 2) CNT(C_NORM_EQ2, 1); else if (run >= 3) CNT(C_NORM_EQ3P, 1); run = 1; } } }
#endif
    int w = 0;
    for (int k = 0; k < n; ++k) {
      if (w > 0 && p[w - 1].mono == kbuf[k].mono)
        p[w - 1].coef = top ? COEF_LIN : te_bin(OP_ADD, p[w - 1].coef, kbuf[k].coef);
      else { p[w].coef = kbuf[k].coef; p[w].mono = kbuf[k].mono; ++w; }
    }
    return w;
  }
  int fn_mono_sum(uint8_t op, int sid) {
    if (ovf) return 0;
    int f = fac_intern(FAC_FN, op, sid, sum[sid].h);
    return mon_single(f);
  }
  int fn_mono(uint8_t op, const PR* args, int n) {
    int sid = sum_intern(args, n);
    if (ovf) return 0;
    int f = fac_intern(FAC_FN, op, sid, sum[sid].h);
    return mon_single(f);
  }
};

struct NFOut { const PR* cols; int q; };

inline int pivot_sum(Ctx& C, const PR* V, int nv, int p, PR* A) {
  int k = 0;
  A[k++] = PR{COEF_FIX, V[p].mono};
  for (int j = 0; j < nv; ++j) if (j != p) A[k++] = PR{C.te_bin(OP_DIV, V[j].coef, V[p].coef), V[j].mono};
  return C.norm(A, k, false);
}
static inline int nf_coll_bits() { static const int v = getenv("NF_COLL") ? atoi(getenv("NF_COLL")) : 0; return v; }
inline PR collapse(Ctx& C, const PR* V, int nv) {
  int k = pivot_sum(C, V, nv, 0, C.scrC);
  int sid = C.sum_intern(C.scrC, k);
  int f = C.fac_intern(FAC_COLL, 0, sid, C.sum[sid].h);
  return PR{V[0].coef, C.mon_single(f)};
}

#ifdef NFC_NI
#define NFC_NOINL __attribute__((noinline))
#else
#define NFC_NOINL
#endif
NFC_NOINL static int op_pow(Ctx& C, const PR* V, int nv, int power, PR* d) {
  int cnt = 1;
  if (power == 3) CNT(nv == 1 ? C_CUBE_NV1 : nv == 2 ? C_CUBE_NV2 : nv == 3 ? C_CUBE_NV3 : C_CUBE_BIG, 1); else CNT(nv == 1 ? C_CUBE_NV1 + 0 : nv == 2 ? C_SQ_NV2 : nv == 3 ? C_SQ_NV3 : nv == 4 ? C_SQ_NV4 : nv == 5 ? C_SQ_NV5 : C_SQ_BIG, nv == 1 ? 0 : 1);
  bool done = false;
  if (nv == 1) { d[0] = PR{C.te_bin(OP_MUL, V[0].coef, V[0].coef), C.mon_mul(V[0].mono, V[0].mono)};
                 if (power == 3) { d[0] = PR{C.te_bin(OP_MUL, d[0].coef, V[0].coef), C.mon_mul(d[0].mono, V[0].mono)}; }
                 done = true; }
  long nq = (power == 2) ? (long)nv * (nv + 1) / 2 : (long)nv * (nv + 1) * (nv + 2) / 6;
  if (!done && nq <= 16 && !(nf_coll_bits() & 1)) {
    PR* A = C.scrA; int na = 0;
    for (int a = 0; a < nv; ++a)
      for (int b = 0; b < nv; ++b) {
        if (b < a) { A[na] = A[b * nv + a]; ++na; }
        else A[na++] = PR{C.te_bin(OP_MUL, V[a].coef, V[b].coef), C.mon_mul(V[a].mono, V[b].mono)};
      }
    if (power == 2) { na = C.norm(A, na, false); for (int k = 0; k < na; ++k) d[k] = A[k]; cnt = na; }
    else {
      PR* B = C.scrB; int nb = 0;
      for (int a = 0; a < na; ++a) {
        const int ra = a / nv, rb = a - ra * nv;
        for (int b = 0; b < nv; ++b) {
          if (rb < ra) { B[nb] = B[(rb * nv + ra) * nv + b]; ++nb; }
          else B[nb++] = PR{C.te_bin(OP_MUL, A[a].coef, V[b].coef), C.mon_mul(A[a].mono, V[b].mono)};
        }
      }
      nb = C.norm(B, nb, false, true);
      for (int k = 0; k < nb; ++k) d[k] = B[k];
      cnt = nb;
    }
    done = true;
  }
  if (!done) {
    PR cb = collapse(C, V, nv);
    d[0] = PR{C.te_bin(OP_MUL, cb.coef, cb.coef), C.mon_mul(cb.mono, cb.mono)};
    if (power == 3) d[0] = PR{C.te_bin(OP_MUL, d[0].coef, cb.coef), C.mon_mul(d[0].mono, cb.mono)};
    cnt = 1;
  }
  return cnt;
}
static inline int find_const(const PR* V, int nv) { int ki = -1; for (int k = 0; k < nv; ++k) if (V[k].mono == 0) ki = k; return ki; }
NFC_NOINL static int op_sincos(Ctx& C, const PR* V, int nv, int o, PR* d) {
  const int ki = find_const(V, nv);
  int cnt = 1;
  if (ki >= 0) {
    int kc = V[ki].coef;
    int k2 = 0; PR* X = C.scrC;
    for (int k = 0; k < nv; ++k) if (k != ki) X[k2++] = V[k];
    const int sidX = C.sum_intern(X, k2);
    int mS = C.fn_mono_sum(OP_SIN, sidX);
    int mC = C.fn_mono_sum(OP_COS, sidX);
    int cc = C.te_un(OP_COS, kc), sc = C.te_un(OP_SIN, kc);
    if (o == OP_SIN) { d[0] = PR{cc, mS}; d[1] = PR{sc, mC}; }
    else { d[0] = PR{cc, mC}; d[1] = PR{C.te_neg(sc), mS}; }
    cnt = C.norm(d, 2, false);
  } else d[0] = PR{2, C.fn_mono((uint8_t)o, V, nv)};
  return cnt;
}
NFC_NOINL static int op_exp(Ctx& C, const PR* V, int nv, PR* d) {
  const int ki = find_const(V, nv);
  if (ki >= 0) {
    int k2 = 0; PR* X = C.scrC;
    for (int k = 0; k < nv; ++k) if (k != ki) X[k2++] = V[k];
    d[0] = PR{C.te_un(OP_EXP, V[ki].coef), C.fn_mono(OP_EXP, X, k2)};
  } else d[0] = PR{2, C.fn_mono(OP_EXP, V, nv)};
  return 1;
}
NFC_NOINL static int op_log(Ctx& C, const PR* V, int nv, PR* d) {
  const int ki = find_const(V, nv);
  int cnt = 1;
  int p = (ki >= 0) ? ki : 0;
  if (ki < 0 && V[0].coef == 2 && nv >= 1) {
    int k2 = pivot_sum(C, V, nv, 0, C.scrC);
    d[0] = PR{2, C.fn_mono(OP_LOG, C.scrC, k2)};
  } else {
    int k2 = pivot_sum(C, V, nv, p, C.scrC);
    int m = C.fn_mono(OP_LOG, C.scrC, k2);
    d[0] = PR{C.te_un(OP_LOG, V[p].coef), 0};
    d[1] = PR{2, m};
    cnt = C.norm(d, 2, false);
  }
  return cnt;
}
NFC_NOINL static int op_root(Ctx& C, const PR* V, int nv, int o, PR* d) {
  if (nv == 1) d[0] = PR{2, C.fn_mono((uint8_t)o, V, 1)};
  else {
    int k2 = pivot_sum(C, V, nv, 0, C.scrC);
    d[0] = PR{C.te_un((uint8_t)o, V[0].coef), C.fn_mono((uint8_t)o, C.scrC, k2)};
  }
  return 1;
}
NFC_NOINL static int op_addsub(Ctx& C, const PR* L, int nl, const PR* R, int nr, bool sub, PR* d) {
  int k = 0;
  for (int a = 0; a < nl; ++a) d[k++] = L[a];
  for (int b = 0; b < nr; ++b) {
    PR x = R[b];
    if (sub) x.coef = C.te_neg(x.coef);
    d[k++] = x;
  }
  return C.norm(d, k, false);
}
NFC_NOINL static int op_mul(Ctx& C, const PR* L, int nl, const PR* R, int nr, PR* d) {
  PR L1, R1;
  if ((nf_coll_bits() & 2) && (nl > 1 || nr > 1)) {
    if (nl > 1) { L1 = collapse(C, L, nl); L = &L1; nl = 1; }
    if (nr > 1) { R1 = collapse(C, R, nr); R = &R1; nr = 1; }
  } else if (nl > 1 && nr > 1) {
    if (nl <= nr) { L1 = collapse(C, L, nl); L = &L1; nl = 1; }
    else { R1 = collapse(C, R, nr); R = &R1; nr = 1; }
  }
  if (nl * nr > MAXP2) return -1;
  int k = 0;
  for (int a = 0; a < nl; ++a)
    for (int b = 0; b < nr; ++b) d[k++] = PR{C.te_bin(OP_MUL, L[a].coef, R[b].coef), C.mon_mul(L[a].mono, R[b].mono)};
  return C.norm(d, k, false);
}
NFC_NOINL static int op_div(Ctx& C, const PR* L, int nl, const PR* R, int nr, bool rtf, PR* d) {
  int cnt;
  if (rtf) {
    for (int a = 0; a < nl; ++a) d[a] = PR{C.te_bin(OP_DIV, L[a].coef, R[0].coef), L[a].mono};
    cnt = nl;
  } else {
    PR den = (nr == 1) ? R[0] : collapse(C, R, nr);
    int rm = C.mon_single(C.fac_intern(FAC_RMONO, 0, den.mono, C.mon[den.mono].h));
    int k = 0;
    for (int a = 0; a < nl; ++a) d[k++] = PR{C.te_bin(OP_DIV, L[a].coef, den.coef), C.mon_mul(L[a].mono, rm)};
    cnt = C.norm(d, k, false);
  }
  return cnt;
}

inline bool nf_tree(Ctx& C, const uint8_t* op, const float* c, int n, NFOut* out) {
  (void)c;
  if (n < 1 || n > MAXN) return false;
  C.bump_gen();
  CNT(C_TREES, 1);
#ifdef NFC_COUNT
  g_mm_n = 0; g_tb_n = 0;
#endif
  C.nte = 3; C.nfac = 1; C.nmon = 0; C.nmp = 0; C.nsum = 0; C.nsp = 0; C.nvp = 0; C.ovf = false;
  C.te[0] = TE{TE_LIN, 0, -1, -1, -1, H_LIN};
  C.te[1] = TE{TE_FIX, 0, -1, -1, -1, H_FIX};
  C.te[2] = TE{TE_ONEV, 0, -1, -1, -1, H_ONEV};
  C.fac[0] = FAC{FAC_T, 0, -1, H_T};
  C.nfac = 1;
  { uint64_t z = 0; C.mon_intern_h(nullptr, &z, 0, H_EMPTY); }
  int monoT = C.mon_single(0);

  int lc[MAXN], rc[MAXN], sp = 0;
  int stk[MAXN + 2]; int* st = stk + 2;
  bool tf[MAXN];
  static const int8_t ARITY[16] = {-1, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -1};
  for (int i = 0; i < n; ++i) {
    const uint8_t o = op[i];
    const int ar = o < 16 ? ARITY[o] : -1;
    if (ar < 0 || sp < ar) return false;
    const int Rr = st[sp - 1], Ll = st[sp - 2];
    lc[i] = ar == 2 ? Ll : ar == 1 ? Rr : -1;
    rc[i] = ar == 2 ? Rr : -1;
    tf[i] = ar == 0 ? (o == OP_CONST) : ar == 1 ? tf[Rr] : (tf[Ll] && tf[Rr]);
    sp += 1 - ar;
    st[sp - 1] = i;
  }
  if (sp != 1) return false;

  PR* val[MAXN]; int vl[MAXN];
  for (int i = 0; i < n; ++i) {
    uint8_t o = op[i];
    if (C.nvp + MAXP2 >= MAXPAIR) return false;
    PR* d = C.vp + C.nvp;
    int cnt = 1;
    if (o == OP_CONST) d[0] = PR{C.te_leaf(i), 0};
    else if (o == OP_T) d[0] = PR{C.te_leaf(i), monoT};
    else if (tf[i]) {
      if (o >= OP_SIN) d[0] = PR{C.te_un(o, val[lc[i]][0].coef), 0};
      else d[0] = PR{C.te_bin(o, val[lc[i]][0].coef, val[rc[i]][0].coef), 0};
    } else if (o >= OP_SIN) {
      const PR* V = val[lc[i]]; int nv = vl[lc[i]];
      if (o == OP_SQ || o == OP_CUBE) cnt = op_pow(C, V, nv, o == OP_SQ ? 2 : 3, d);
      else if (o == OP_SIN || o == OP_COS) cnt = op_sincos(C, V, nv, o, d);
      else if (o == OP_EXP) cnt = op_exp(C, V, nv, d);
      else if (o == OP_LOG) cnt = op_log(C, V, nv, d);
      else cnt = op_root(C, V, nv, o, d);
    } else if (o == OP_ADD || o == OP_SUB) {
      cnt = op_addsub(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], o == OP_SUB, d);
    } else if (o == OP_MUL) {
      cnt = op_mul(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], d);
      if (cnt < 0) return false;
    } else {
      cnt = op_div(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], tf[rc[i]], d);
    }
    if (C.ovf) return false;
    val[i] = d; vl[i] = cnt; C.nvp += cnt;
  }
  CNT(C_NTE, C.nte); CNT(C_NMON, C.nmon); CNT(C_NFAC, C.nfac); CNT(C_NSUM, C.nsum); CNT(C_NVP, C.nvp);
  const PR* root = val[n - 1];
  int q = vl[n - 1];
  if (q > MAXQ) return false;
  nf_add_const(root, q);
  out->cols = root; out->q = q;
  return true;
}

struct Prog {
  uint32_t w[WMAX];
  uint32_t* out = nullptr;
  int nw, nabs, nS, nV, Q, nth;
  bool wide;
};

namespace detail {
inline uint16_t push_ins(Ctx& C, uint8_t op, uint16_t a, uint16_t b, uint32_t mask, bool vec) {
  if (C.ni >= MAXINS) { C.ovf = true; return 0; }
  ++C.ni;
  if (vec) { int i = C.nvS++; C.iv_[i] = Ctx::Ins{mask, a, b, op}; return (uint16_t)(0x8000 | i); }
  int i = C.nsS++; C.is_[i] = Ctx::Ins{mask, a, b, op}; return (uint16_t)i;
}
inline const Ctx::Ins& ins_at(const Ctx& C, uint16_t r) { return (r & 0x8000) ? C.iv_[r & 0x7FFF] : C.is_[r]; }
inline uint16_t push_bin(Ctx& C, uint8_t op, uint16_t a, uint16_t b) {
  const Ctx::Ins& A = ins_at(C, a); const Ctx::Ins& B = ins_at(C, b);
  return push_ins(C, op, a, b, A.mask | B.mask, ((a | b) & 0x8000) != 0);
}
inline uint16_t push_un(Ctx& C, uint8_t op, uint16_t a) {
  const Ctx::Ins& A = ins_at(C, a);
  return push_ins(C, op, a, 0, A.mask, (a & 0x8000) != 0);
}
inline uint16_t ensure_one(Ctx& C) {
  if (!C.have_one) { C.one_ref = push_ins(C, I_ONE, 0, 0, 0, false); C.have_one = true; }
  return C.one_ref;
}
uint16_t emit_mono(Ctx& C, int m);
inline uint16_t emit_sum(Ctx& C, int s, bool freq);
inline uint16_t emit_fac(Ctx& C, int f) {
  CNT(C_EMIT_FAC, 1);
  if (C.mf[f] >= C.mbase) return (uint16_t)(C.mf[f] - C.mbase);
  const FAC& F = C.fac[f];
  uint16_t r;
  if (F.kind == FAC_T) r = push_ins(C, I_T, 0, 0, 0, true);
  else if (F.kind == FAC_FN) {
    uint16_t s = emit_sum(C, F.arg, F.op == OP_SIN || F.op == OP_COS);
    r = push_un(C, F.op, s);
  } else if (F.kind == FAC_COLL) { r = emit_sum(C, F.arg, false);
  } else if (F.kind == FAC_RSUM) { uint16_t one = ensure_one(C); uint16_t s = emit_sum(C, F.arg, false); r = push_bin(C, OP_DIV, one, s); }
  else { uint16_t one = ensure_one(C); uint16_t m = emit_mono(C, F.arg); r = push_bin(C, OP_DIV, one, m); }
  C.mf[f] = C.mbase + r;
  return r;
}
inline uint16_t emit_mono(Ctx& C, int m) {
  CNT(C_EMIT_MONO, 1);
  if (C.mm[m] >= C.mbase) return (uint16_t)(C.mm[m] - C.mbase);
  const MON& M = C.mon[m];
  uint16_t r;
  if (M.len == 0) r = ensure_one(C);
  else {
    r = emit_fac(C, C.mp[M.off]);
    for (int k = 1; k < M.len; ++k) { uint16_t f = emit_fac(C, C.mp[M.off + k]); r = push_bin(C, OP_MUL, r, f); }
  }
  C.mm[m] = C.mbase + r;
  return r;
}
inline uint16_t emit_term(Ctx& C, PR p) {
  uint16_t mi = emit_mono(C, p.mono);
  if (p.coef == COEF_LIN || p.coef == COEF_FIX) return mi;
  if (C.mt[p.coef] < C.mbase) {
    if (C.nth >= MAXTH) { C.ovf = true; return mi; }
    int slot = C.nth++;
    C.slot_te[slot] = p.coef;
    C.slot_freq[slot] = C.freqTE[p.coef] == C.mbase;
    C.mt[p.coef] = C.mbase + push_ins(C, I_THETA, (uint16_t)slot, 0, 1u << slot, false);
  }
  return push_bin(C, OP_MUL, (uint16_t)(C.mt[p.coef] - C.mbase), mi);
}
inline uint16_t emit_sum(Ctx& C, int s, bool freq) {
  if (C.ms[s] >= C.mbase && !freq) return (uint16_t)(C.ms[s] - C.mbase);
  const SUMO& S = C.sum[s];
  uint16_t r = emit_term(C, C.sp[S.off]);
  for (int k = 1; k < S.len; ++k) { uint16_t t = emit_term(C, C.sp[S.off + k]); r = push_bin(C, OP_ADD, r, t); }
  C.ms[s] = C.mbase + r;
  return r;
}
}

inline double te_val(Ctx& C, const float* c, int id) {
  if (C.tgen[id] == C.vgen) return C.tevals[id];
  const TE& E = C.te[id];
  CNT(C_TEVAL_NODES, 1);
  double v = 0;
  switch (E.kind) {
    case TE_LIN: v = 0; break;
    case TE_FIX: v = 1; break;
    case TE_ONEV: v = 1; break;
    case TE_LEAF: v = (double)c[E.leaf]; break;
    case TE_NEG: v = -te_val(C, c, E.a); break;
    case TE_UN: {
      double x = te_val(C, c, E.a);
      if (E.op != OP_SQ && E.op != OP_CUBE) CNT(C_TEVAL_LIBM, 1);
      switch (E.op) {
        case OP_SIN: v = std::sin(x); break;
        case OP_COS: v = std::cos(x); break;
        case OP_EXP: v = std::exp(x); break;
        case OP_LOG: v = std::log(std::fabs(x)); break;
        case OP_SQ: v = x * x; break;
        case OP_CUBE: v = x * x * x; break;
        case OP_SQRT: v = std::sqrt(std::fabs(x)); break;
        case OP_CBRT: v = std::cbrt(x); break;
        default: v = x;
      }
      break;
    }
    default: {
      double x = te_val(C, c, E.a), y = te_val(C, c, E.b);
      switch (E.op) {
        case OP_ADD: v = x + y; break;
        case OP_SUB: v = x - y; break;
        case OP_MUL: v = x * y; break;
        default: v = x / y;
      }
    }
  }
  C.tgen[id] = C.vgen; C.tevals[id] = v;
  return v;
}

#ifdef NFC_ABL
inline int g_abl = 99;
#define ABL(n) if (g_abl == (n)) return true;
#else
#define ABL(n)
#endif
inline bool compile(Ctx& C, const NFOut& nf, const float* c, Prog& P) {
  C.mbase += 65536;
  if (C.mbase == 0) { memset(C.mm, 0, sizeof C.mm); memset(C.mf, 0, sizeof C.mf); memset(C.ms, 0, sizeof C.ms); memset(C.mt, 0, sizeof C.mt);
                      memset(C.freqTE, 0, sizeof C.freqTE); C.mbase = 65536; }
  C.ni = 0; C.nsS = 0; C.nvS = 0; C.nth = 0; C.have_one = false; C.ovf = false;
  for (int f = 0; f < C.nfac; ++f)
    if (C.fac[f].kind == FAC_FN && (C.fac[f].op == OP_SIN || C.fac[f].op == OP_COS)) {
      const SUMO& S = C.sum[C.fac[f].arg];
      for (int k = 0; k < S.len; ++k) { int cf = C.sp[S.off + k].coef; if (cf >= 2) C.freqTE[cf] = C.mbase; }
    }
  uint16_t roots[MAXQ];
  int Q = nf.q;
  for (int k = 0; k < Q; ++k) roots[k] = detail::emit_mono(C, nf.cols[k].mono);
  if (C.ovf) return false;
  if (++C.vgen == 0) { memset(C.tgen, 0, sizeof C.tgen); C.vgen = 1; }
  const int n = C.ni;
  if (n > 8191 || C.nth > 32) return false;
  const int nS = C.nsS, nV = C.nvS;
  const bool wide = C.nth > 8;
  uint32_t* w = P.out ? P.out : P.w; int k = 0;
  w[k++] = (1u << 24) | ((uint32_t)Q << 16) | ((uint32_t)C.nth << 8) | (wide ? 1u : 0u);
  w[k++] = (uint32_t)nS | ((uint32_t)nV << 16);
  uint32_t fmask = 0;
  for (int s = 0; s < C.nth; ++s) if (C.slot_freq[s]) fmask |= 1u << s;
  w[k++] = fmask;
  const int rbase = k, nr = (Q + 1) / 2;
  uint32_t* ins = w + rbase + nr;
  uint32_t* mk = ins + n;
  const int mw = wide ? n : (n + 3) / 4;
  if (!wide) mk[mw - 1] = 0;
  uint8_t* mk8 = (uint8_t*)mk;
  auto res = [nS](uint32_t r) -> uint32_t { return (r & 0x8000u) ? nS + (r & 0x7FFFu) : r; };
  for (int i = 0; i < nS; ++i) {
    const Ctx::Ins& I = C.is_[i];
    ins[i] = (uint32_t)I.op | (res(I.a) << 6) | (res(I.b) << 19);
    if (wide) mk[i] = I.mask; else mk8[i] = (uint8_t)(I.mask & 0xFFu);
  }
  for (int i = 0; i < nV; ++i) {
    const Ctx::Ins& I = C.iv_[i];
    ins[nS + i] = (uint32_t)I.op | (res(I.a) << 6) | (res(I.b) << 19);
    if (wide) mk[nS + i] = I.mask; else mk8[nS + i] = (uint8_t)(I.mask & 0xFFu);
  }
  for (int q = 0; q < Q; q += 2) {
    uint32_t lo = res(roots[q]);
    uint32_t hi = (q + 1 < Q) ? res(roots[q + 1]) : 0u;
    w[rbase + q / 2] = lo | (hi << 16);
  }
  k = rbase + nr + n + mw;
  P.nabs = k;
  for (int s = 0; s < C.nth; ++s) {
    double vd = te_val(C, c, C.slot_te[s]);
    float v = (float)vd;
    uint32_t bits; memcpy(&bits, &v, 4);
    w[k++] = bits;
  }
  P.nw = k; P.nS = nS; P.nV = nV; P.Q = Q; P.nth = C.nth; P.wide = wide;
  return true;
}

struct Child {
  uint32_t koff = 0, klen = 0;
  uint32_t roff[2] = {0, 0}, rlen[2] = {0, 0};
  uint8_t ok[2] = {0, 0};
  uint64_t r0 = 0, r1 = 0, kh = 0;
};
constexpr size_t CHUNK = 64;
constexpr size_t CHUNK_WORDS = CHUNK * 2 * (1 + (size_t)WMAX + 1);
struct Arena { uint32_t* w; size_t cap; size_t n; };
static inline uint64_t sm(uint64_t x) {
  x += 0x9e3779b97f4a7c15ull; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}
static inline uint64_t hash_words(const uint32_t* w, size_t n) {
  uint64_t h = 0x243f6a8885a308d3ull ^ n;
  for (size_t i = 0; i < n; i++) h = sm(h ^ w[i]);
  return h;
}
static inline void need(const Arena& A, size_t k) {
  if (A.n + k > A.cap) { fprintf(stderr, "nfc: arena overflow\n"); abort(); }
}
static inline bool tree_flat(Ctx& C, Prog& P, const skn::Cand& cd, int t, Arena& A, Child& ch, uint32_t& hs, uint32_t& hn) {
  const size_t t0 = A.n;
  const int n = cd.len[t];
  bool good = false;
  NFOut out;
  need(A, 1 + (size_t)WMAX);
  P.out = A.w + A.n + 1;
  if (nf_tree(C, cd.op[t], cd.c[t], n, &out)) good = compile(C, out, cd.c[t], P);
  if (good && P.Q <= 16) {
    const uint32_t nw = (uint32_t)P.nw;
    A.w[A.n++] = nw;
    const size_t r0 = A.n;
    A.n += nw;
    ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw;
    hs = (uint32_t)r0; hn = nw;
    return true;
  }
  A.n = t0;
  const size_t s0 = A.n;
  need(A, 2 + 2 * (size_t)n);
  A.w[A.n++] = 0xFFFFFFFFu; A.w[A.n++] = (uint32_t)n;
  for (int j = 0; j < n; j++) { uint32_t b; memcpy(&b, &cd.c[t][j], 4); A.w[A.n++] = cd.op[t][j]; A.w[A.n++] = b; }
  hs = (uint32_t)s0; hn = (uint32_t)(A.n - s0);
  ch.rlen[t] = 0;
  return false;
}
static inline void hash4(const uint32_t* const w[4], const uint32_t n[4], uint64_t r[4]) {
  int id[4] = {0, 1, 2, 3};
  auto cs = [&](int x, int y) { if (n[id[y]] < n[id[x]]) { int t = id[x]; id[x] = id[y]; id[y] = t; } };
  cs(0, 1); cs(2, 3); cs(0, 2); cs(1, 3); cs(1, 2);
  const uint32_t* w0 = w[id[0]]; const uint32_t* w1 = w[id[1]]; const uint32_t* w2 = w[id[2]]; const uint32_t* w3 = w[id[3]];
  const uint32_t n0 = n[id[0]], n1 = n[id[1]], n2 = n[id[2]], n3 = n[id[3]];
  uint64_t h0 = 0x243f6a8885a308d3ull ^ n0, h1 = 0x243f6a8885a308d3ull ^ n1, h2 = 0x243f6a8885a308d3ull ^ n2, h3 = 0x243f6a8885a308d3ull ^ n3;
  uint32_t i = 0;
  for (; i < n0; i++) { h0 = sm(h0 ^ w0[i]); h1 = sm(h1 ^ w1[i]); h2 = sm(h2 ^ w2[i]); h3 = sm(h3 ^ w3[i]); }
  for (; i < n1; i++) { h1 = sm(h1 ^ w1[i]); h2 = sm(h2 ^ w2[i]); h3 = sm(h3 ^ w3[i]); }
  for (; i < n2; i++) { h2 = sm(h2 ^ w2[i]); h3 = sm(h3 ^ w3[i]); }
  for (; i < n3; i++) { h3 = sm(h3 ^ w3[i]); }
  r[id[0]] = h0; r[id[1]] = h1; r[id[2]] = h2; r[id[3]] = h3;
}
static void chunk_flat(Ctx& C, Prog& P, const skn::Cand* const* cands, const uint8_t* valid, int b0, int b1, Arena& A, Child* out) {
  A.n = 0;
  uint32_t hs[CHUNK][2], hn[CHUNK][2];
  for (int i = b0; i < b1; i++) {
    Child& ch = out[i];
    ch = Child();
    const int ii = i - b0;
    hs[ii][0] = hs[ii][1] = 0; hn[ii][0] = hn[ii][1] = 0;
    if (!valid[i]) { ch.koff = (uint32_t)A.n; continue; }
    const skn::Cand& cd = *cands[i];
    ch.koff = (uint32_t)A.n;
    ch.ok[0] = tree_flat(C, P, cd, 0, A, ch, hs[ii][0], hn[ii][0]) ? 1 : 0;
    ch.ok[1] = tree_flat(C, P, cd, 1, A, ch, hs[ii][1], hn[ii][1]) ? 1 : 0;
    ch.klen = (uint32_t)(A.n - ch.koff);
  }
  const int cnt = (b1 - b0) * 2;
  for (int g = 0; g < cnt; g += 4) {
    const uint32_t* w[4]; uint32_t n[4]; uint64_t r[4];
    for (int k = 0; k < 4; k++) { int c = g + k; if (c < cnt) { w[k] = A.w + hs[c >> 1][c & 1]; n[k] = hn[c >> 1][c & 1]; } else { w[k] = A.w; n[k] = 0; } }
    for (int k = 0; k < 4; k++) CNT(C_HASHWORDS, n[k]);
    hash4(w, n, r);
    for (int k = 0; k < 4; k++) { int c = g + k; if (c < cnt) { int i = b0 + (c >> 1); if (valid[i]) (c & 1 ? out[i].r1 : out[i].r0) = r[k]; } }
  }
  for (int i = b0; i < b1; i++) if (valid[i]) { Child& ch = out[i]; ch.kh = sm(ch.r0 ^ sm(ch.r1 + ch.klen)); }
}

}
