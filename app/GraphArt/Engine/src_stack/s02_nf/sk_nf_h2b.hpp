#pragma once
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <atomic>
#include <cstdlib>

namespace skn {

enum : uint8_t {
  OP_NOP = 0, OP_CONST = 1, OP_T = 2, OP_ADD = 3, OP_SUB = 4, OP_MUL = 5, OP_DIV = 6,
  OP_SIN = 7, OP_COS = 8, OP_EXP = 9, OP_LOG = 10, OP_SQ = 11, OP_CUBE = 12,
  OP_SQRT = 13, OP_CBRT = 14
};
enum : uint8_t { I_THETA = 1, I_T = 2, I_NEG = 15, I_ONE = 16 };

struct Cand { uint8_t len[2]; uint8_t op[2][24]; float c[2][24]; };

constexpr int MAXN = 24;
constexpr int MAXTE = 1024, MAXFAC = 512, MAXMONO = 1024, MAXSUM = 512;
constexpr int MAXMP = 1 << 15, MAXSP = 1 << 15, MAXPAIR = 1 << 17, MAXP2 = 1 << 14;
constexpr int TABN = 2048;
constexpr int MAXINS = 512, MAXTH = 32, MAXQ = 255;
constexpr int MAXSCR = 512;
constexpr int WMAX = 2 + 128 + MAXINS + MAXINS + MAXTH + 8;
static_assert(MAXTE <= TABN / 2 && MAXMONO <= TABN / 2 && MAXFAC <= TABN / 2 && MAXSUM <= TABN / 2, "intern table load <= 1/2");

constexpr int COEF_LIN = 0, COEF_FIX = 1;
enum : uint8_t { TE_LIN = 0, TE_FIX = 1, TE_ONEV = 2, TE_LEAF = 3, TE_BIN = 4, TE_UN = 5, TE_NEG = 6 };
enum : uint8_t { FAC_T = 0, FAC_FN = 1, FAC_RSUM = 2, FAC_RMONO = 3, FAC_COLL = 4 };

struct TE { uint8_t kind, op; int32_t a, b, leaf; uint64_t h; };
struct FAC { uint8_t kind, op; int32_t arg; uint64_t h; };
struct MON { int32_t off, len; uint64_t h; };
struct SUMO { int32_t off, len; uint64_t h; };
struct PR { int32_t coef, mono; };
static inline bool eqpr(const PR* a, const PR* b, int n) { for (int k = 0; k < n; ++k) if (a[k].coef != b[k].coef || a[k].mono != b[k].mono) return false; return true; }

static constexpr uint64_t H_LIN = 0x1001, H_FIX = 0x1002, H_ONEV = 0x1007, H_LEAF = 0x1003, H_BIN = 0x1004,
                          H_UN = 0x1005, H_NEG = 0x1006, H_T = 0x2001, H_FN = 0x2002,
                          H_RS = 0x2003, H_RM = 0x2004, H_COLL = 0x2005, H_EMPTY = 0x3001, H_SUM = 0x3002;

static inline uint64_t mx(uint64_t h, uint64_t x) {
  h = (h ^ (x + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2))) * 0xff51afd7ed558ccdull;
  return h ^ (h >> 33);
}
static inline bool cur(uint64_t e, uint32_t g) { return (uint32_t)(e >> 32) == g; }
static inline bool eq32(const int32_t* a, const int32_t* b, int n) { for (int k = 0; k < n; ++k) if (a[k] != b[k]) return false; return true; }

struct KE { uint64_t h; int32_t coef, mono; };

struct Ctx {
  TE te[MAXTE]; int nte;
  KE kbuf[MAXP2];
  FAC fac[MAXFAC]; int nfac;
  MON mon[MAXMONO]; int nmon;
  int32_t mp[MAXMP]; int nmp;
  SUMO sum[MAXSUM]; int nsum;
  PR sp[MAXSP]; int nsp;
  PR vp[MAXPAIR]; int nvp;
  uint64_t tte[TABN], tfac[TABN], tmon[TABN], tsum[TABN];
  uint32_t gen;
  bool ovf;
  int32_t scr[MAXSCR];
  PR scrA[MAXP2], scrB[MAXP2], scrC[MAXP2];
  uint8_t iop[MAXINS]; int32_t ia[MAXINS], ib[MAXINS];
  uint8_t icls[MAXINS]; uint32_t imask[MAXINS];
  int ni, nth, one_instr;
  int32_t slot_te[MAXTH]; uint8_t slot_freq[MAXTH];
  uint8_t freqTE[MAXTE];
  double tevals[MAXTE];
  uint8_t tdone[MAXTE];
  int32_t mm[MAXMONO], mf[MAXFAC], ms[MAXSUM], mt[MAXTE];
  Ctx() { gen = 0; memset(tte, 0, sizeof tte); memset(tfac, 0, sizeof tfac);
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
    for (;;) {
      uint64_t e = tte[p];
      if (!cur(e, gen)) break;
      const TE& x = te[(uint32_t)e - 1];
      if (x.kind == kind && x.op == op && x.a == a && x.b == b) return (int)(uint32_t)e - 1;
      p = (p + 1) & (TABN - 1);
    }
    if (nte >= MAXTE) { ovf = true; return COEF_LIN; }
    te[nte] = TE{kind, op, a, b, -1, h};
    tte[p] = ((uint64_t)gen << 32) | (uint32_t)(nte + 1);
    return nte++;
  }
  int te_bin(uint8_t op, int a, int b) {
    if (op == OP_MUL) { if (a == COEF_FIX) return b; if (b == COEF_FIX) return a; }
    if (op == OP_DIV && b == COEF_FIX) return a;
    if ((op == OP_ADD || op == OP_MUL) && (te[b].h < te[a].h || (te[b].h == te[a].h && b < a))) { int t = a; a = b; b = t; }
    uint64_t h = mx(mx(mx(H_BIN, op), te[a].h), te[b].h);
    return te_intern(TE_BIN, op, a, b, h);
  }
  int te_un(uint8_t op, int a) { return te_intern(TE_UN, op, a, -1, mx(mx(H_UN, op), te[a].h)); }
  int te_neg(int a) { return te_intern(TE_NEG, 0, a, -1, mx(H_NEG, te[a].h)); }

  void sort_fac(int32_t* f, int n) const {
    for (int i = 1; i < n; ++i) {
      int32_t x = f[i]; uint64_t hx = fac[x].h; int j = i - 1;
      while (j >= 0 && (fac[f[j]].h > hx || (fac[f[j]].h == hx && f[j] > x))) { f[j + 1] = f[j]; --j; }
      f[j + 1] = x;
    }
  }
  int mon_intern(const int32_t* f, int n) {
    uint64_t h = H_EMPTY;
    for (int k = 0; k < n; ++k) h = mx(h, fac[f[k]].h);
    uint32_t p = (uint32_t)(h & (TABN - 1));
    for (;;) {
      uint64_t e = tmon[p];
      if (!cur(e, gen)) break;
      const MON& M = mon[(uint32_t)e - 1];
      if (M.len == n && M.h == h && eq32(mp + M.off, f, n)) return (int)(uint32_t)e - 1;
      p = (p + 1) & (TABN - 1);
    }
    if (nmon >= MAXMONO || nmp + n > MAXMP) { ovf = true; return 0; }
    for (int k = 0; k < n; ++k) mp[nmp + k] = f[k];
    mon[nmon] = MON{nmp, n, h};
    nmp += n;
    tmon[p] = ((uint64_t)gen << 32) | (uint32_t)(nmon + 1);
    return nmon++;
  }
  int mon_single(int f) { int32_t x = f; return mon_intern(&x, 1); }
  int mon_mul(int m1, int m2) {
    int n1 = mon[m1].len, n2 = mon[m2].len;
    if (n1 + n2 > MAXSCR) { ovf = true; return 0; }
    for (int k = 0; k < n1; ++k) scr[k] = mp[mon[m1].off + k];
    for (int k = 0; k < n2; ++k) scr[n1 + k] = mp[mon[m2].off + k];
    sort_fac(scr, n1 + n2);
    return mon_intern(scr, n1 + n2);
  }
  int fac_intern(uint8_t kind, uint8_t op, int32_t arg, uint64_t argh) {
    uint64_t base = kind == FAC_FN ? H_FN : kind == FAC_RSUM ? H_RS : kind == FAC_RMONO ? H_RM : kind == FAC_COLL ? H_COLL : H_T;
    uint64_t h = (kind == FAC_T) ? H_T : mx(mx(base, op), argh);
    uint32_t p = (uint32_t)(h & (TABN - 1));
    for (;;) {
      uint64_t e = tfac[p];
      if (!cur(e, gen)) break;
      const FAC& x = fac[(uint32_t)e - 1];
      if (x.kind == kind && x.op == op && x.arg == arg) return (int)(uint32_t)e - 1;
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
    for (;;) {
      uint64_t e = tsum[p];
      if (!cur(e, gen)) break;
      const SUMO& S = sum[(uint32_t)e - 1];
      if (S.len == n && S.h == h && eqpr(sp + S.off, pr, n)) return (int)(uint32_t)e - 1;
      p = (p + 1) & (TABN - 1);
    }
    if (nsum >= MAXSUM || nsp + n > MAXSP) { ovf = true; return 0; }
    for (int k = 0; k < n; ++k) sp[nsp + k] = pr[k];
    sum[nsum] = SUMO{nsp, n, h};
    nsp += n;
    tsum[p] = ((uint64_t)gen << 32) | (uint32_t)(nsum + 1);
    return nsum++;
  }
  int norm(PR* p, int n, bool top) {
    if (n <= 1) return n;
    if (n <= MAXP2) {
      for (int k = 0; k < n; ++k) { kbuf[k].h = mon[p[k].mono].h; kbuf[k].coef = p[k].coef; kbuf[k].mono = p[k].mono; }
      std::sort(kbuf, kbuf + n, [](const KE& a, const KE& b) {
        return a.h < b.h || (a.h == b.h && a.mono < b.mono);
      });
      for (int k = 0; k < n; ++k) { p[k].coef = kbuf[k].coef; p[k].mono = kbuf[k].mono; }
    } else {
      std::sort(p, p + n, [this](const PR& a, const PR& b) {
        uint64_t ha = mon[a.mono].h, hb = mon[b.mono].h;
        return ha < hb || (ha == hb && a.mono < b.mono);
      });
    }
    int w = 0;
    for (int k = 0; k < n; ++k) {
      if (w > 0 && p[w - 1].mono == p[k].mono)
        p[w - 1].coef = top ? COEF_LIN : te_bin(OP_ADD, p[w - 1].coef, p[k].coef);
      else p[w++] = p[k];
    }
    return w;
  }
  int fn_mono(uint8_t op, const PR* args, int n) {
    int sid = sum_intern(args, n);
    if (ovf) return 0;
    int f = fac_intern(FAC_FN, op, sid, sum[sid].h);
    return mon_single(f);
  }
};

struct NFOut { const PR* cols; int q; };
static inline int free_offset() { static const int v = getenv("FREE_OFFSET") ? atoi(getenv("FREE_OFFSET")) : 1; return v; }
static inline void nf_add_const(const PR* root, int& q) {
  if (!free_offset() || q < 1 || q >= 16) return;
  for (int k = 0; k < q; ++k) if (root[k].mono == 0) return;
  const_cast<PR*>(root)[q] = PR{COEF_LIN, 0}; ++q;
}


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

inline bool nf_tree(Ctx& C, const uint8_t* op, const float* c, int n, NFOut* out) {
  (void)c;
  if (n < 1 || n > MAXN) return false;
  C.bump_gen();
  C.nte = 3; C.nfac = 1; C.nmon = 0; C.nmp = 0; C.nsum = 0; C.nsp = 0; C.nvp = 0; C.ovf = false;
  C.te[0] = TE{TE_LIN, 0, -1, -1, -1, H_LIN};
  C.te[1] = TE{TE_FIX, 0, -1, -1, -1, H_FIX};
  C.te[2] = TE{TE_ONEV, 0, -1, -1, -1, H_ONEV};
  C.fac[0] = FAC{FAC_T, 0, -1, H_T};
  C.nfac = 1;
  C.mon_intern(nullptr, 0);
  int monoT = C.mon_single(0);

  int lc[MAXN], rc[MAXN], st[MAXN], sp = 0;
  bool tf[MAXN];
  for (int i = 0; i < n; ++i) {
    uint8_t o = op[i];
    lc[i] = rc[i] = -1;
    if (o == OP_CONST) { tf[i] = true; st[sp++] = i; }
    else if (o == OP_T) { tf[i] = false; st[sp++] = i; }
    else if (o >= OP_SIN && o <= OP_CBRT) {
      if (sp < 1) return false;
      int a = st[--sp]; lc[i] = a; tf[i] = tf[a]; st[sp++] = i;
    } else if (o >= OP_ADD && o <= OP_DIV) {
      if (sp < 2) return false;
      int R = st[--sp], L = st[--sp];
      lc[i] = L; rc[i] = R; tf[i] = tf[L] && tf[R]; st[sp++] = i;
    } else return false;
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
      int ki = -1;
      for (int k = 0; k < nv; ++k) if (V[k].mono == 0) ki = k;
      if (o == OP_SQ || o == OP_CUBE) {
        int power = (o == OP_SQ) ? 2 : 3;
        bool done = false;
        if (nv == 1) { d[0] = PR{C.te_bin(OP_MUL, V[0].coef, V[0].coef), C.mon_mul(V[0].mono, V[0].mono)};
                       if (power == 3) { d[0] = PR{C.te_bin(OP_MUL, d[0].coef, V[0].coef), C.mon_mul(d[0].mono, V[0].mono)}; }
                       done = true; }
        long nq = (power == 2) ? (long)nv * (nv + 1) / 2 : (long)nv * (nv + 1) * (nv + 2) / 6;
        if (!done && nq <= 16 && !(nf_coll_bits() & 1)) {
          PR* A = C.scrA; int na = 0;
          for (int a = 0; a < nv; ++a)
            for (int b = 0; b < nv; ++b)
              A[na++] = PR{C.te_bin(OP_MUL, V[a].coef, V[b].coef), C.mon_mul(V[a].mono, V[b].mono)};
          if (power == 2) { na = C.norm(A, na, false); for (int k = 0; k < na; ++k) d[k] = A[k]; cnt = na; }
          else {
            PR* B = C.scrB; int nb = 0;
            for (int a = 0; a < na; ++a)
              for (int b = 0; b < nv; ++b)
                B[nb++] = PR{C.te_bin(OP_MUL, A[a].coef, V[b].coef), C.mon_mul(A[a].mono, V[b].mono)};
            nb = C.norm(B, nb, false);
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
      } else if (o == OP_SIN || o == OP_COS) {
        if (ki >= 0) {
          int kc = V[ki].coef;
          int k2 = 0; PR* X = C.scrC;
          for (int k = 0; k < nv; ++k) if (k != ki) X[k2++] = V[k];
          int mS = C.fn_mono(OP_SIN, X, k2);
          int mC = C.fn_mono(OP_COS, X, k2);
          int cc = C.te_un(OP_COS, kc), sc = C.te_un(OP_SIN, kc);
          if (o == OP_SIN) { d[0] = PR{cc, mS}; d[1] = PR{sc, mC}; }
          else { d[0] = PR{cc, mC}; d[1] = PR{C.te_neg(sc), mS}; }
          cnt = C.norm(d, 2, false);
        } else d[0] = PR{2, C.fn_mono((uint8_t)o, V, nv)};
      } else if (o == OP_EXP) {
        if (ki >= 0) {
          int k2 = 0; PR* X = C.scrC;
          for (int k = 0; k < nv; ++k) if (k != ki) X[k2++] = V[k];
          d[0] = PR{C.te_un(OP_EXP, V[ki].coef), C.fn_mono(OP_EXP, X, k2)};
        } else d[0] = PR{2, C.fn_mono(OP_EXP, V, nv)};
      } else if (o == OP_LOG) {
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
      } else {
        if (nv == 1) d[0] = PR{2, C.fn_mono((uint8_t)o, V, 1)};
        else {
          int k2 = pivot_sum(C, V, nv, 0, C.scrC);
          d[0] = PR{C.te_un((uint8_t)o, V[0].coef), C.fn_mono((uint8_t)o, C.scrC, k2)};
        }
      }
    } else if (o == OP_ADD || o == OP_SUB) {
      const PR* L = val[lc[i]]; int nl = vl[lc[i]];
      const PR* R = val[rc[i]]; int nr = vl[rc[i]];
      int k = 0;
      for (int a = 0; a < nl; ++a) d[k++] = L[a];
      for (int b = 0; b < nr; ++b) {
        PR x = R[b];
        if (o == OP_SUB) x.coef = C.te_neg(x.coef);
        d[k++] = x;
      }
      cnt = C.norm(d, k, false);
    } else if (o == OP_MUL) {
      PR L1, R1; const PR* L = val[lc[i]]; int nl = vl[lc[i]];
      const PR* R = val[rc[i]]; int nr = vl[rc[i]];
      if ((nf_coll_bits() & 2) && (nl > 1 || nr > 1)) {
        if (nl > 1) { L1 = collapse(C, L, nl); L = &L1; nl = 1; }
        if (nr > 1) { R1 = collapse(C, R, nr); R = &R1; nr = 1; }
      } else if (nl > 1 && nr > 1) {
        if (nl <= nr) { L1 = collapse(C, L, nl); L = &L1; nl = 1; }
        else { R1 = collapse(C, R, nr); R = &R1; nr = 1; }
      }
      if (nl * nr > MAXP2) return false;
      int k = 0;
      for (int a = 0; a < nl; ++a)
        for (int b = 0; b < nr; ++b) d[k++] = PR{C.te_bin(OP_MUL, L[a].coef, R[b].coef), C.mon_mul(L[a].mono, R[b].mono)};
      cnt = C.norm(d, k, false);
    } else {
      const PR* L = val[lc[i]]; int nl = vl[lc[i]];
      const PR* R = val[rc[i]]; int nr = vl[rc[i]];
      if (tf[rc[i]]) {
        for (int a = 0; a < nl; ++a) d[a] = PR{C.te_bin(OP_DIV, L[a].coef, R[0].coef), L[a].mono};
        cnt = nl;
      } else {
        PR den = (nr == 1) ? R[0] : collapse(C, R, nr);
        int rm = C.mon_single(C.fac_intern(FAC_RMONO, 0, den.mono, C.mon[den.mono].h));
        int k = 0;
        for (int a = 0; a < nl; ++a) d[k++] = PR{C.te_bin(OP_DIV, L[a].coef, den.coef), C.mon_mul(L[a].mono, rm)};
        cnt = C.norm(d, k, false);
      }
    }
    if (C.ovf) return false;
    val[i] = d; vl[i] = cnt; C.nvp += cnt;
  }
  const PR* root = val[n - 1];
  int q = vl[n - 1];
  if (q > MAXQ) return false;
  nf_add_const(root, q);
  out->cols = root; out->q = q;
  return true;
}

struct Prog {
  uint32_t w[WMAX];
  int nw, nabs, nS, nV, Q, nth;
  bool wide;
  double thd[MAXTH];
  float theta[MAXTH];
  uint64_t h1, h2;
};

namespace detail {
inline int push(Ctx& C, uint8_t op, int a, int b) {
  if (C.ni >= MAXINS) { C.ovf = true; return 0; }
  const int i = C.ni;
  C.iop[i] = op; C.ia[i] = a; C.ib[i] = b;
  if (op == I_T) { C.icls[i] = 1; C.imask[i] = 0; }
  else if (op == I_THETA) { C.icls[i] = 0; C.imask[i] = 1u << a; }
  else if (op == I_ONE) { C.icls[i] = 0; C.imask[i] = 0; }
  else if (op == OP_ADD || op == OP_SUB || op == OP_MUL || op == OP_DIV) {
    C.icls[i] = C.icls[a] | C.icls[b]; C.imask[i] = C.imask[a] | C.imask[b];
  } else { C.icls[i] = C.icls[a]; C.imask[i] = C.imask[a]; }
  return C.ni++;
}
inline int ensure_one(Ctx& C) { if (C.one_instr < 0) C.one_instr = push(C, I_ONE, 0, 0); return C.one_instr; }
int emit_mono(Ctx& C, int m);
inline int emit_sum(Ctx& C, int s, bool freq);
inline int emit_fac(Ctx& C, int f) {
  if (C.mf[f] >= 0) return C.mf[f];
  const FAC& F = C.fac[f];
  int r;
  if (F.kind == FAC_T) r = push(C, I_T, 0, 0);
  else if (F.kind == FAC_FN) {
    int s = emit_sum(C, F.arg, F.op == OP_SIN || F.op == OP_COS);
    r = push(C, F.op, s, 0);
  } else if (F.kind == FAC_COLL) { r = emit_sum(C, F.arg, false);
  } else if (F.kind == FAC_RSUM) { int one = ensure_one(C); int s = emit_sum(C, F.arg, false); r = push(C, OP_DIV, one, s); }
  else { int one = ensure_one(C); int m = emit_mono(C, F.arg); r = push(C, OP_DIV, one, m); }
  C.mf[f] = r;
  return r;
}
inline int emit_mono(Ctx& C, int m) {
  if (C.mm[m] >= 0) return C.mm[m];
  const MON& M = C.mon[m];
  int r;
  if (M.len == 0) r = ensure_one(C);
  else {
    r = emit_fac(C, C.mp[M.off]);
    for (int k = 1; k < M.len; ++k) { int f = emit_fac(C, C.mp[M.off + k]); r = push(C, OP_MUL, r, f); }
  }
  C.mm[m] = r;
  return r;
}
inline int emit_term(Ctx& C, PR p, bool freq) {
  (void)freq;
  int mi = emit_mono(C, p.mono);
  if (p.coef == COEF_LIN || p.coef == COEF_FIX) return mi;
  if (C.mt[p.coef] < 0) {
    if (C.nth >= MAXTH) { C.ovf = true; return mi; }
    int slot = C.nth++;
    C.slot_te[slot] = p.coef;
    C.slot_freq[slot] = C.freqTE[p.coef];
    C.mt[p.coef] = push(C, I_THETA, slot, 0);
  }
  return push(C, OP_MUL, C.mt[p.coef], mi);
}
inline int emit_sum(Ctx& C, int s, bool freq) {
  if (C.ms[s] >= 0 && !freq) return C.ms[s];
  const SUMO& S = C.sum[s];
  int r = emit_term(C, C.sp[S.off], freq);
  for (int k = 1; k < S.len; ++k) { int t = emit_term(C, C.sp[S.off + k], freq); r = push(C, OP_ADD, r, t); }
  C.ms[s] = r;
  return r;
}
}

inline double te_val(Ctx& C, const float* c, int id) {
  if (C.tdone[id]) return C.tevals[id];
  const TE& E = C.te[id];
  double v = 0;
  switch (E.kind) {
    case TE_LIN: v = 0; break;
    case TE_FIX: v = 1; break;
    case TE_ONEV: v = 1; break;
    case TE_LEAF: v = (double)c[E.leaf]; break;
    case TE_NEG: v = -te_val(C, c, E.a); break;
    case TE_UN: {
      double x = te_val(C, c, E.a);
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
  C.tdone[id] = 1; C.tevals[id] = v;
  return v;
}

inline bool compile(Ctx& C, const NFOut& nf, const float* c, Prog& P) {
  int nobj_m = C.nmon, nobj_f = C.nfac, nobj_s = C.nsum, nobj_t = C.nte;
  for (int i = 0; i < nobj_m; ++i) C.mm[i] = -1;
  for (int i = 0; i < nobj_f; ++i) C.mf[i] = -1;
  for (int i = 0; i < nobj_s; ++i) C.ms[i] = -1;
  for (int i = 0; i < nobj_t; ++i) { C.mt[i] = -1; C.freqTE[i] = 0; }
  C.ni = 0; C.nth = 0; C.one_instr = -1; C.ovf = false;
  for (int f = 0; f < C.nfac; ++f)
    if (C.fac[f].kind == FAC_FN && (C.fac[f].op == OP_SIN || C.fac[f].op == OP_COS)) {
      const SUMO& S = C.sum[C.fac[f].arg];
      for (int k = 0; k < S.len; ++k) { int cf = C.sp[S.off + k].coef; if (cf >= 2) C.freqTE[cf] = 1; }
    }
  int roots[MAXQ];
  int Q = nf.q;
  for (int k = 0; k < Q; ++k) roots[k] = detail::emit_mono(C, nf.cols[k].mono);
  if (C.ovf) return false;
  for (int i = 0; i < nobj_t; ++i) C.tdone[i] = 0;
  int newid[MAXINS]; int nS = 0, nV = 0;
  for (int i = 0; i < C.ni; ++i) if (!C.icls[i]) newid[i] = nS++;
  for (int i = 0; i < C.ni; ++i) if (C.icls[i]) newid[i] = nS + nV++;
  int n = C.ni;
  uint8_t op2[MAXINS]; int a2[MAXINS], b2[MAXINS]; uint32_t m2[MAXINS];
  for (int i = 0; i < n; ++i) {
    int j = newid[i];
    uint8_t o = C.iop[i];
    op2[j] = o; m2[j] = C.imask[i];
    if (o == I_THETA) { a2[j] = C.ia[i]; b2[j] = 0; }
    else if (o == I_T || o == I_ONE) { a2[j] = 0; b2[j] = 0; }
    else if (o == OP_ADD || o == OP_SUB || o == OP_MUL || o == OP_DIV) { a2[j] = newid[C.ia[i]]; b2[j] = newid[C.ib[i]]; }
    else { a2[j] = newid[C.ia[i]]; b2[j] = 0; }
  }
  if (n > 8191 || C.nth > 32) return false;
  bool wide = C.nth > 8;
  uint32_t* w = P.w; int k = 0;
  w[k++] = (1u << 24) | ((uint32_t)Q << 16) | ((uint32_t)C.nth << 8) | (wide ? 1u : 0u);
  w[k++] = (uint32_t)nS | ((uint32_t)nV << 16);
  uint32_t fmask = 0;
  for (int s = 0; s < C.nth; ++s) if (C.slot_freq[s]) fmask |= 1u << s;
  w[k++] = fmask;
  for (int q = 0; q < Q; q += 2) {
    uint32_t lo = (uint32_t)newid[roots[q]];
    uint32_t hi = (q + 1 < Q) ? (uint32_t)newid[roots[q + 1]] : 0u;
    w[k++] = lo | (hi << 16);
  }
  for (int i = 0; i < n; ++i) w[k++] = (uint32_t)op2[i] | ((uint32_t)a2[i] << 6) | ((uint32_t)b2[i] << 19);
  if (wide) { for (int i = 0; i < n; ++i) w[k++] = m2[i]; }
  else {
    for (int i = 0; i < n; i += 4) {
      uint32_t x = 0;
      for (int j = 0; j < 4 && i + j < n; ++j) x |= (m2[i + j] & 0xFFu) << (8 * j);
      w[k++] = x;
    }
  }
  P.nabs = k;
  for (int s = 0; s < C.nth; ++s) {
    double vd = te_val(C, c, C.slot_te[s]);
    P.thd[s] = vd;
    float v = (float)vd;
    P.theta[s] = v;
    uint32_t bits; memcpy(&bits, &v, 4);
    w[k++] = bits;
  }
  P.nw = k; P.nS = nS; P.nV = nV; P.Q = Q; P.nth = C.nth; P.wide = wide;
  uint64_t h1 = 0x243F6A8885A308D3ull, h2 = 0x13198A2E03707344ull;
  for (int i = 0; i < P.nabs; ++i) { h1 = mx(h1, w[i]); h2 = mx(h2 ^ 0xA5A5A5A5ull, (uint64_t)w[i] * 0x9E3779B1ull + i); }
  P.h1 = h1; P.h2 = h2;
  return true;
}

struct View {
  const uint32_t* w;
  int Q, nth, nS, nV; bool wide;
  uint32_t freq; const uint32_t* roots_w; const uint32_t* ins; const uint32_t* masks; const uint32_t* theta;
  explicit View(const uint32_t* p) : w(p) {
    Q = (p[0] >> 16) & 0xFF; nth = (p[0] >> 8) & 0xFF; wide = p[0] & 1;
    nS = p[1] & 0xFFFF; nV = p[1] >> 16;
    freq = p[2];
    roots_w = p + 3;
    int nr = (Q + 1) / 2;
    ins = roots_w + nr;
    int n = nS + nV;
    masks = ins + n;
    int mw = wide ? n : (n + 3) / 4;
    theta = masks + mw;
  }
  int n() const { return nS + nV; }
  int root(int q) const { uint32_t x = roots_w[q / 2]; return (q & 1) ? (int)(x >> 16) : (int)(x & 0xFFFF); }
  uint8_t op(int i) const { return ins[i] & 0x3F; }
  int a(int i) const { return (ins[i] >> 6) & 0x1FFF; }
  int b(int i) const { return (ins[i] >> 19) & 0x1FFF; }
  uint32_t mask(int i) const { return wide ? masks[i] : ((masks[i >> 2] >> (8 * (i & 3))) & 0xFF); }
  float th(int s) const { float v; memcpy(&v, theta + s, 4); return v; }
};


constexpr int MAXNODES = 256;
struct PNode { uint8_t kind, op; uint16_t a, b; int32_t leaf; };

inline void replay_nodes(const PNode* nd, int pn_len, const uint16_t* slot_root, int nth, const float* c, Prog& P) {
  double v[MAXNODES];
  for (int i = 0; i < pn_len; ++i) {
    const PNode& N = nd[i];
    double r;
    switch (N.kind) {
      case TE_LIN: r = 0; break;
      case TE_FIX: r = 1; break;
      case TE_ONEV: r = 1; break;
      case TE_LEAF: r = (double)c[N.leaf]; break;
      case TE_NEG: r = -v[N.a]; break;
      case TE_UN: {
        double x = v[N.a];
        switch (N.op) {
          case OP_SIN: r = std::sin(x); break;
          case OP_COS: r = std::cos(x); break;
          case OP_EXP: r = std::exp(x); break;
          case OP_LOG: r = std::log(std::fabs(x)); break;
          case OP_SQ: r = x * x; break;
          case OP_CUBE: r = x * x * x; break;
          case OP_SQRT: r = std::sqrt(std::fabs(x)); break;
          case OP_CBRT: r = std::cbrt(x); break;
          default: r = x;
        }
        break;
      }
      default: {
        double x = v[N.a], y = v[N.b];
        switch (N.op) {
          case OP_ADD: r = x + y; break;
          case OP_SUB: r = x - y; break;
          case OP_MUL: r = x * y; break;
          default: r = x / y;
        }
      }
    }
    v[i] = r;
  }
  for (int s = 0; s < nth; ++s) {
    double vd = v[slot_root[s]];
    P.thd[s] = vd;
    float f = (float)vd;
    P.theta[s] = f;
    memcpy(P.w + P.nabs + s, &f, 4);
  }
}

struct DagBuilder {
  int32_t lmap[MAXTE]; uint32_t lgen[MAXTE]; uint32_t lg;
  std::vector<PNode> out;
  DagBuilder() : lg(0) { memset(lgen, 0, sizeof lgen); out.reserve(MAXNODES); }
  bool emit_te(const Ctx& C, int id) {
    if (lgen[id] == lg) return true;
    const TE& T = C.te[id];
    PNode nd{T.kind, T.op, 0, 0, T.leaf};
    if (T.kind == TE_NEG || T.kind == TE_UN) {
      if (!emit_te(C, T.a)) return false;
      nd.a = (uint16_t)lmap[T.a];
    } else if (T.kind == TE_BIN || (T.kind != TE_LIN && T.kind != TE_FIX && T.kind != TE_ONEV && T.kind != TE_LEAF)) {
      if (!emit_te(C, T.a)) return false;
      if (!emit_te(C, T.b)) return false;
      nd.a = (uint16_t)lmap[T.a]; nd.b = (uint16_t)lmap[T.b];
    }
    if (out.size() >= (size_t)MAXNODES) return false;
    lmap[id] = (int32_t)out.size(); lgen[id] = lg;
    out.push_back(nd);
    return true;
  }
  bool build(const Ctx& C, uint16_t* slot_root) {
    lg++;
    if (lg == 0) { memset(lgen, 0, sizeof lgen); lg = 1; }
    out.clear();
    for (int s = 0; s < C.nth; ++s) {
      if (!emit_te(C, C.slot_te[s])) return false;
      slot_root[s] = (uint16_t)lmap[C.slot_te[s]];
    }
    return true;
  }
};

struct NFMemo {
  struct Key { uint8_t len; uint8_t op[MAXN]; };
  struct Ent {
    Key k;
    uint32_t aw_off, pn_off;
    uint16_t nabs, nS, nV, nth, pn_len;
    uint16_t slot_root[MAXTH];
    uint8_t status, Q, wide;
    uint64_t h1, h2;
  };
  static constexpr uint32_t TB = 1u << 18;
  static constexpr uint32_t MAXENT = 1u << 17;
  static constexpr size_t MAXAW = 1u << 24;
  static constexpr size_t MAXPN = 1u << 22;

  std::vector<uint32_t> tab;
  std::vector<Ent> ent;
  std::vector<uint32_t> aw;
  std::vector<PNode> pn;
  DagBuilder DB;
  long hits, misses, flushes;
  NFMemo() : tab(TB, 0), hits(0), misses(0), flushes(0) {
    ent.reserve(MAXENT); aw.reserve(MAXAW / 8); pn.reserve(MAXPN / 8);
  }
  void prefault() {
    ent.resize(MAXENT); ent.clear();
    aw.resize(MAXAW / 8); aw.clear();
    pn.resize(MAXPN / 8); pn.clear();
  }
  void flush() {
    std::fill(tab.begin(), tab.end(), 0u);
    ent.clear(); aw.clear(); pn.clear();
    ++flushes;
  }
  static inline uint64_t hkey(const Key& k) {
    uint64_t h = 0x243F6A8885A308D3ull ^ k.len;
    for (int i = 0; i < MAXN; i += 8) { uint64_t x; memcpy(&x, k.op + i, 8); h = mx(h, x); }
    return h;
  }
  bool nf_compile(Ctx& C, const uint8_t* op, const float* c, int n, Prog& P, int maxq, int* st = nullptr) {
    Key k; memset(&k, 0, sizeof k);
    k.len = (uint8_t)n;
    if (n >= 0 && n <= MAXN) memcpy(k.op, op, (size_t)n);
    const uint64_t h = hkey(k);
    uint32_t p = (uint32_t)(h & (TB - 1));
    while (tab[p] != 0) {
      const Ent& E = ent[tab[p] - 1];
      if (memcmp(&E.k, &k, sizeof k) == 0) {
        ++hits;
        if (st) *st = E.status == 2 ? 0 : (E.status == 0 ? 1 : 2);
        if (E.status != 2 || E.Q > maxq) return false;
        P.nabs = E.nabs; P.nS = E.nS; P.nV = E.nV; P.Q = E.Q; P.nth = E.nth; P.wide = E.wide != 0;
        P.nw = E.nabs + E.nth;
        memcpy(P.w, aw.data() + E.aw_off, sizeof(uint32_t) * E.nabs);
        P.h1 = E.h1; P.h2 = E.h2;
        replay_nodes(pn.data() + E.pn_off, E.pn_len, E.slot_root, E.nth, c, P);
        return true;
      }
      p = (p + 1) & (TB - 1);
    }
    ++misses;
    if (ent.size() >= MAXENT || aw.size() > MAXAW || pn.size() > MAXPN) { flush(); p = (uint32_t)(h & (TB - 1)); }
    NFOut out;
    bool good = false;
    uint8_t status = 0;
    if (nf_tree(C, op, c, n, &out)) {
      good = compile(C, out, c, P);
      status = good ? 2 : 1;
    }
    if (st) *st = good ? 0 : (status == 0 ? 1 : 2);
    ent.emplace_back(); Ent* E = &ent.back();
    E->k = k; E->status = status; E->aw_off = 0; E->pn_off = 0; E->nabs = E->nS = E->nV = E->nth = E->pn_len = 0;
    E->Q = 0; E->wide = 0; E->h1 = E->h2 = 0;
    if (good) {
      E->nabs = (uint16_t)P.nabs; E->nS = (uint16_t)P.nS; E->nV = (uint16_t)P.nV; E->Q = (uint8_t)P.Q;
      E->nth = (uint16_t)P.nth; E->wide = P.wide ? 1 : 0; E->h1 = P.h1; E->h2 = P.h2;
      if (DB.build(C, E->slot_root)) {
        E->pn_len = (uint16_t)DB.out.size();
        E->aw_off = (uint32_t)aw.size();
        E->pn_off = (uint32_t)pn.size();
        aw.insert(aw.end(), P.w, P.w + P.nabs);
        pn.insert(pn.end(), DB.out.begin(), DB.out.end());
      } else { ent.pop_back(); E = nullptr; }
    }
    if (E) tab[p] = (uint32_t)ent.size();
    return good && P.Q <= maxq;
  }
};

struct NFShared {
  static constexpr uint32_t TB = 1u << 19;
  static constexpr uint32_t MAXENT = 1u << 18;
  static constexpr size_t MAXAW = (size_t)1 << 25;
  static constexpr size_t MAXPN = (size_t)1 << 23;
  using Ent = NFMemo::Ent;
  std::atomic<uint32_t>* tab;
  Ent* ent;
  uint32_t* aw;
  PNode* pn;
  std::atomic<uint32_t> nent, naw, npn;
  NFShared() : tab(new std::atomic<uint32_t>[TB]()) {
    ent = (Ent*)calloc(MAXENT, sizeof(Ent));
    aw = (uint32_t*)calloc(MAXAW, sizeof(uint32_t));
    pn = (PNode*)calloc(MAXPN, sizeof(PNode));
    nent.store(0); naw.store(0); npn.store(0);
  }
  ~NFShared() { delete[] tab; free(ent); free(aw); free(pn); }
  void prefault(size_t bytes) {
    memset(ent, 0, std::min(bytes, sizeof(Ent) * (size_t)MAXENT));
    memset(aw, 0, std::min(bytes, sizeof(uint32_t) * MAXAW));
    memset(pn, 0, std::min(bytes, sizeof(PNode) * MAXPN));
  }
  NFShared(const NFShared&) = delete;
  NFShared& operator=(const NFShared&) = delete;
};

struct NFHandle {
  NFShared& S;
  DagBuilder DB;
  long hits, misses, drops, bigskip;
  explicit NFHandle(NFShared& s) : S(s), hits(0), misses(0), drops(0), bigskip(0) {}

  void publish(uint64_t h, const NFMemo::Ent& X, const uint32_t* words, int nabs, const PNode* nodes, int npn_) {
    uint32_t idx = S.nent.fetch_add(1, std::memory_order_relaxed);
    if (idx >= NFShared::MAXENT) { ++drops; return; }
    uint32_t aoff = S.naw.fetch_add((uint32_t)nabs, std::memory_order_relaxed);
    if ((size_t)aoff + (size_t)nabs > NFShared::MAXAW) { ++drops; return; }
    uint32_t poff = S.npn.fetch_add((uint32_t)npn_, std::memory_order_relaxed);
    if ((size_t)poff + (size_t)npn_ > NFShared::MAXPN) { ++drops; return; }
    NFMemo::Ent& D = S.ent[idx];
    D = X;
    D.aw_off = aoff; D.pn_off = poff;
    if (nabs) memcpy(S.aw + aoff, words, sizeof(uint32_t) * (size_t)nabs);
    if (npn_) memcpy(S.pn + poff, nodes, sizeof(PNode) * (size_t)npn_);
    uint32_t p = (uint32_t)(h & (NFShared::TB - 1));
    for (;;) {
      uint32_t expected = 0;
      if (S.tab[p].compare_exchange_strong(expected, idx + 1, std::memory_order_release, std::memory_order_acquire)) return;
      if (memcmp(&S.ent[expected - 1].k, &X.k, sizeof(X.k)) == 0) return;
      p = (p + 1) & (NFShared::TB - 1);
    }
  }

  bool nf_compile(Ctx& C, const uint8_t* op, const float* c, int n, Prog& P, int maxq, int* st = nullptr) {
    NFMemo::Key k; memset(&k, 0, sizeof k);
    k.len = (uint8_t)n;
    if (n >= 0 && n <= MAXN) memcpy(k.op, op, (size_t)n);
    const uint64_t h = NFMemo::hkey(k);
    uint32_t p = (uint32_t)(h & (NFShared::TB - 1));
    for (;;) {
      uint32_t v = S.tab[p].load(std::memory_order_acquire);
      if (v == 0) break;
      const NFMemo::Ent& E = S.ent[v - 1];
      if (memcmp(&E.k, &k, sizeof k) == 0) {
        ++hits;
        if (st) *st = E.status == 2 ? 0 : (E.status == 0 ? 1 : 2);
        if (E.status != 2 || E.Q > maxq) return false;
        P.nabs = E.nabs; P.nS = E.nS; P.nV = E.nV; P.Q = E.Q; P.nth = E.nth; P.wide = E.wide != 0;
        P.nw = E.nabs + E.nth;
        memcpy(P.w, S.aw + E.aw_off, sizeof(uint32_t) * E.nabs);
        P.h1 = E.h1; P.h2 = E.h2;
        replay_nodes(S.pn + E.pn_off, E.pn_len, E.slot_root, E.nth, c, P);
        return true;
      }
      p = (p + 1) & (NFShared::TB - 1);
    }
    ++misses;
    NFOut out;
    bool good = false;
    uint8_t status = 0;
    if (nf_tree(C, op, c, n, &out)) {
      good = compile(C, out, c, P);
      status = good ? 2 : 1;
    }
    if (st) *st = good ? 0 : (status == 0 ? 1 : 2);
    NFMemo::Ent X; memset(&X, 0, sizeof X);
    X.k = k; X.status = status;
    bool cache = true;
    if (good) {
      X.nabs = (uint16_t)P.nabs; X.nS = (uint16_t)P.nS; X.nV = (uint16_t)P.nV; X.Q = (uint8_t)P.Q;
      X.nth = (uint16_t)P.nth; X.wide = P.wide ? 1 : 0; X.h1 = P.h1; X.h2 = P.h2;
      if (DB.build(C, X.slot_root)) {
        X.pn_len = (uint16_t)DB.out.size();
        publish(h, X, P.w, P.nabs, DB.out.data(), (int)DB.out.size());
      } else { cache = false; ++bigskip; }
    } else publish(h, X, nullptr, 0, nullptr, 0);
    (void)cache;
    return good && P.Q <= maxq;
  }
};

}
