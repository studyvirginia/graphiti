#pragma once
#include <cstdint>
#include <cstring>
#include <vector>
#include "compile_memo_i13_h156.h"

namespace enc34 {

struct Dec { uint16_t zo, d, A, B; uint8_t d63, fuse; };

struct Scratch {
  std::vector<Dec> dec;
  std::vector<uint32_t> ob, ob2;
  uint32_t segmap[18];
  size_t nout = 0;
};

enum PlanResult { PLAN_OK = 0, PLAN_DECODE_FAIL = 1, PLAN_FALLBACK = 2 };

struct DecTab { uint8_t cls[64], zN[64], zB[64], zBoth[64]; uint16_t bothB[64]; };
static constexpr DecTab make_dectab() {
  DecTab t{};
  for (uint32_t op = 0; op < 64; op++) {
    t.cls[op] = op >= 3 && op <= 6 ? 0 : (op >= 7 && op <= 15) ? 1 : 2;
    t.zN[op] = (uint8_t)(op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4);
    t.zB[op] = (uint8_t)(op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5);
    t.zBoth[op] = (uint8_t)(op == 5 ? 11 : op == 3 ? 0 : 255);
    t.bothB[op] = (uint16_t)(op == 5 ? 901 : 2000);
  }
  return t;
}
static constexpr DecTab kDec = make_dectab();

static inline bool decode_one(uint32_t w, Dec& x) {
  const uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22;
  const uint32_t cls = kDec.cls[op];
  const uint32_t cs = (a == 1022 && b == 1022) ? 2u : (b == 1022 ? 1u : 0u);
  const uint32_t zbin = cs == 0 ? kDec.zN[op] : cs == 1 ? kDec.zB[op] : kDec.zBoth[op];
  if (cls == 0 && zbin == 255) return false;
  x.zo = (uint16_t)(cls == 0 ? zbin : cls == 1 ? op : 6u);
  x.A = (uint16_t)(cls == 0 && cs == 1 ? b : a);
  x.B = (uint16_t)(cls != 0 ? 901u : cs == 0 ? b : cs == 1 ? a : kDec.bothB[op]);
  x.d = (uint16_t)d; x.d63 = d == 63; x.fuse = 0;
  return true;
}

static inline int plan(const Tree& T, Scratch& S) {
  if (T.Q < 0 || T.Q > 17) return PLAN_FALLBACK;
  const size_t n = T.prog.size();
  if (S.dec.size() < n) S.dec.resize(n);
  Dec* D = S.dec.data();
  size_t dec_end = 0;
  auto ensure = [&](size_t i) -> bool {
    while (dec_end <= i) { if (!decode_one(T.prog[dec_end], D[dec_end])) return false; dec_end++; }
    return true;
  };
  const int Q = T.Q;
  for (int q = 1; q <= Q; q++) if (T.segs[q] < T.segs[q - 1]) return PLAN_FALLBACK;
  size_t k = 0, out = 0; int qp = 0;
  while (k < n) {
    if (!ensure(k)) return PLAN_DECODE_FAIL;
    while (qp <= Q && (size_t)T.segs[qp] == k) S.segmap[qp++] = (uint32_t)out;
    bool fuse = false;
    if (k + 1 < n && D[k].zo == 0 && D[k].d63) {
      if (!ensure(k + 1)) return PLAN_DECODE_FAIL;
      fuse = D[k + 1].zo >= 7 && D[k + 1].zo <= 15 && D[k + 1].A == 1022 && !(qp <= Q && (size_t)T.segs[qp] == k + 1);
    }
    if (fuse) { D[k].fuse = 1; k += 2; }
    else { k += 1; }
    out++;
  }
  while (qp <= Q && (size_t)T.segs[qp] == n) S.segmap[qp++] = (uint32_t)out;
  if (qp != Q + 1) return PLAN_FALLBACK;
  S.nout = out;
  return PLAN_OK;
}

struct ScanE {
  uint32_t SL, TT, TWO;
  inline uint32_t enc(uint32_t c) const {
    if (c < 64) return c * SL;
    if (c == 1022) return (1u << 31) | (3u << 29) | (TT + 48);
    if (c == 900) return 2u << 29;
    if (c == 901) return (3u << 29) | (TT + 48);
    if (c == 2000) return TWO;
    return (1u << 29) | (c - 64);
  }
  inline uint32_t dscale(uint32_t d) const { return d * SL; }
};

struct GnE {
  uint32_t CD48, SB, NS, TWO;
  inline uint32_t enc(uint32_t c) const {
    if (c < 64) return c * CD48;
    if (c == 1022) return 1u << 31;
    if (c == 900) return 1u << 29;
    if (c == 901) return (1u << 30) | (SB + NS);
    if (c == 2000) return TWO;
    return (1u << 30) | (SB + (c - 64));
  }
  inline uint32_t dscale(uint32_t d) const { return d * CD48; }
};

template <class E>
static inline void emit(std::vector<uint32_t>& outv, Scratch& S, size_t n, const E& e) {
  const Dec* D = S.dec.data();
  const size_t nout = S.nout;
  if (S.ob.size() < 4 * n) S.ob.resize(4 * n);
  uint32_t* o = S.ob.data();
  size_t w = 0;
  for (size_t k = 0; k < n;) {
    const Dec& a = D[k];
    if (a.fuse) {
      const Dec& b = D[k + 1];
      o[w] = b.zo; o[w + 1] = b.d63 ? 0xFFFFFFFFu : e.dscale(b.d); o[w + 2] = e.enc(a.A); o[w + 3] = e.enc(a.B);
      k += 2;
    } else {
      o[w] = a.zo; o[w + 1] = a.d63 ? 0xFFFFFFFFu : e.dscale(a.d); o[w + 2] = e.enc(a.A); o[w + 3] = e.enc(a.B);
      k += 1;
    }
    w += 4;
  }
  outv.assign(o, o + 4 * nout);
}

static inline void emit_both(Scratch& S, size_t n, const ScanE& es, const GnE* eg, std::vector<uint32_t>& outS, std::vector<uint32_t>* outG) {
  const Dec* D = S.dec.data();
  const size_t nout = S.nout;
  if (S.ob.size() < 4 * n) S.ob.resize(4 * n);
  if (eg && S.ob2.size() < 4 * n) S.ob2.resize(4 * n);
  uint32_t* os = S.ob.data();
  uint32_t* og = S.ob2.data();
  size_t w = 0;
  for (size_t k = 0; k < n;) {
    const Dec& a = D[k];
    const Dec& z = a.fuse ? D[k + 1] : a;
    uint32_t* s = os + w;
    s[0] = z.zo; s[1] = z.d63 ? 0xFFFFFFFFu : es.dscale(z.d); s[2] = es.enc(a.A); s[3] = es.enc(a.B);
    if (eg) {
      uint32_t* g = og + w;
      g[0] = z.zo; g[1] = z.d63 ? 0xFFFFFFFFu : eg->dscale(z.d); g[2] = eg->enc(a.A); g[3] = eg->enc(a.B);
    }
    k += a.fuse ? 2 : 1;
    w += 4;
  }
  outS.assign(os, os + 4 * nout);
  if (eg) outG->assign(og, og + 4 * nout);
}

template <class E>
static inline void segs_emit(uint32_t* sg, const Tree& T, const uint32_t* segmap, const E& e) {
  for (int q = 0; q < 18; q++) sg[q] = (q <= T.Q) ? segmap[q] : (uint32_t)T.segs[q];
  for (int q = 0; q < 18; q++) sg[18 + q] = e.enc(T.segs[18 + q]);
}

static inline ScanE make_scan(const Tree& T, int SS) {
  const uint32_t SL = (uint32_t)SS * 48;
  const int nsc = T.nS <= 4 ? 4 : T.nS <= 8 ? 8 : T.nS <= 16 ? 16 : 32;
  const uint32_t TT = (uint32_t)T.nslot * SL + (uint32_t)(SS * nsc);
  return ScanE{SL, TT, (3u << 29) | (TT + 49)};
}
static inline GnE make_gn(const Tree& T, int CD) {
  const uint32_t SB = (uint32_t)T.nslot * (uint32_t)CD * 48;
  return GnE{(uint32_t)CD * 48, SB, (uint32_t)T.nS, (1u << 30) | (SB + (uint32_t)T.nS + 1)};
}

static inline void scan_h34(const Tree& T, Scratch& S, int SS, std::vector<uint32_t>& p4, uint32_t* segs) {
  const ScanE e = make_scan(T, SS);
  emit(p4, S, T.prog.size(), e);
  segs_emit(segs, T, S.segmap, e);
}

static inline bool gn_h34(const Tree& T, Scratch& S, int CD, std::vector<uint32_t>& p4, uint32_t* segs, uint32_t& hw, size_t& per) {
  const uint32_t NSC = (uint32_t)T.nS + 2, SB = (uint32_t)T.nslot * (uint32_t)CD * 48;
  if (SB >= 65536) return false;
  const GnE e = make_gn(T, CD);
  emit(p4, S, T.prog.size(), e);
  segs_emit(segs, T, S.segmap, e);
  hw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
  per = (size_t)SB + 2 * (size_t)CD * NSC;
  return true;
}

static inline void encode_tree_h34(const Tree& tr, cmi6::TreeEnc& te, Scratch& S) {
  std::memset(te.sg, 0, sizeof te.sg); std::memset(te.gg, 0, sizeof te.gg);
  if (plan(tr, S) != PLAN_OK) { cmi6::encode_tree_i6(tr, te); return; }
  const int SS = tr.nth == 0 ? 2 : 16;
  const ScanE es = make_scan(tr, SS);
  const bool gin = tr.nth >= 1 && tr.nth <= GN_TH_MAX;
  const int CD = tr.nth <= 4 ? 1 + (tr.nth + 1) / 2 : 3;
  const uint32_t SB = (uint32_t)tr.nslot * (uint32_t)CD * 48;
  const bool gok = gin && SB < 65536;
  const GnE eg = make_gn(tr, CD);
  emit_both(S, tr.prog.size(), es, gok ? &eg : nullptr, te.sp4, gok ? &te.gp4 : nullptr);
  segs_emit(te.sg, tr, S.segmap, es);
  te.s_ok = 1;
  if (gok) {
    segs_emit(te.gg, tr, S.segmap, eg);
    const uint32_t NSC = (uint32_t)tr.nS + 2;
    te.ghw = NSC | (SB << 8) | ((uint32_t)tr.nth << 24);
    te.gper = (uint32_t)((size_t)SB + 2 * (size_t)CD * NSC);
    te.g_ok = 1;
  }
}

static inline void encode_tree_h34(const Tree& tr, cmi6::TreeEnc& te) {
  static thread_local Scratch S;
  encode_tree_h34(tr, te, S);
}

}
