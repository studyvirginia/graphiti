#ifndef COFIT_NO_X
#define COFIT_H5_1 1
#define COFIT_H5_2 1
#define COFIT_H5_3 1
#define COFIT_H5_4 1
#endif
#ifndef COFIT_NO_X
#define COFIT_X2 1
#define COFIT_X3 1
#define COFIT_X4 1
#define COFIT_RI_TG 1
#endif
#pragma once
#include <functional>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <arm_neon.h>
#include "fjn_h950.hpp"
#ifdef COFIT_FMA
#define CF_AM(acc_, a_, b_) vfmaq_f32(acc_, a_, b_)
#define CF_SM(acc_, a_, b_) vfmsq_f32(acc_, a_, b_)
#define CF_DOT3(a0_, b0_, a1_, b1_, a2_, b2_) vfmaq_f32(vfmaq_f32(vmulq_f32(a0_, b0_), a1_, b1_), a2_, b2_)
#define CF_TERM(xx_, ss_, gg_, sel_) vfmaq_f32(sel_, vmulq_f32(xx_, ss_), gg_)
#define CF_SACC(acc_, a_, b_) acc_ = std::fmaf(a_, b_, acc_)
#define CF_SSUB(acc_, a_, b_) acc_ = std::fmaf(-(a_), b_, acc_)
#else
#define CF_AM(acc_, a_, b_) vaddq_f32(acc_, vmulq_f32(a_, b_))
#define CF_SM(acc_, a_, b_) vsubq_f32(acc_, vmulq_f32(a_, b_))
#define CF_DOT3(a0_, b0_, a1_, b1_, a2_, b2_) vaddq_f32(vaddq_f32(vmulq_f32(a0_, b0_), vmulq_f32(a1_, b1_)), vmulq_f32(a2_, b2_))
#define CF_TERM(xx_, ss_, gg_, sel_) vaddq_f32(vmulq_f32(vmulq_f32(xx_, ss_), gg_), sel_)
#define CF_SACC(acc_, a_, b_) acc_ += (a_) * (b_)
#define CF_SSUB(acc_, a_, b_) acc_ -= (a_) * (b_)
#endif
#pragma clang fp contract(off)
#include <mach/mach_time.h>
#ifdef COFIT_PROF
static thread_local double g_ptn[8]; static thread_local double g_op[2][16], g_opn[2][16];
struct PTN { int id; uint64_t t; PTN(int i):id(i),t(mach_absolute_time()){} ~PTN(){ g_ptn[id]+= (double)(mach_absolute_time()-t); } };
#define PROFN(i) PTN _ptn##i(i)
#else
#define PROFN(i)
#endif
#ifdef COFIT_PROF
static thread_local long g_cnt[4];
#define COFIT_CNT(i) (g_cnt[i]++)
#else
#define COFIT_CNT(i)
#endif
#ifndef CF_LEANF
#define CF_LEANF false
#endif
#ifdef COFIT_LEANP
#define CF_P2DUAL false
#else
#define CF_P2DUAL true
#endif
#ifdef COFIT_CBRT2
#ifndef COFIT_LANEFIX
#define COFIT_LANEFIX 1
#endif
#endif
namespace cofit {
#define PHASE(i)
#ifdef COFIT_H5_1
#define COFIT_CPYROWS(A_, B_, n_, nb4_) { const int nv_ = (nb4_) / 4 + 1; for (int i_ = 0; i_ < (n_); i_++) { float* d_ = (A_) + i_ * NS; const float* s_ = (B_) + i_ * NS; \
    for (int v_ = 0; v_ < nv_; v_++) { vst1q_f32(d_ + 4 * v_, vld1q_f32(s_ + 4 * v_)); __asm__ volatile("" ::: ); } } }
#else
#define COFIT_CPYROWS(A_, B_, n_, nb4_) { for (int i_ = 0; i_ < (n_); i_++) std::memcpy((A_) + i_ * NS, (B_) + i_ * NS, (size_t)(nb4_) * 4 + 16); }
#endif
static const int MAXC = 16;
static const int NS = 40;
static const int MAXN = 38;
#ifdef COFIT_H5_2
#define COFIT_ABSLT(x_, c_) vcaltq_f32((x_), (c_))
#else
#define COFIT_ABSLT(x_, c_) vcltq_f32(vabsq_f32(x_), (c_))
#endif
inline void set_ftz() {
#if defined(__aarch64__)
  uint64_t f; __asm__ volatile("mrs %0, fpcr" : "=r"(f)); f |= (1ull << 24); __asm__ volatile("msr fpcr, %0" : : "r"(f));
#endif
}
struct Tree {
  int Q = 0, nth = 0, nS = 0, nslot = 0; std::vector<uint16_t> mskD, mskA, mskB; uint16_t mskOut[18]; double coneFrac = 1.0;
  std::vector<uint32_t> tokP, tokD; uint32_t sgP[36], sgD[36]; uint32_t SBp = 0, SBd = 0; std::vector<uint16_t> sp; bool ok = false;
#ifdef COFIT_VREUSE
  std::vector<int32_t> vro; int nvr = 0;
#endif
};
struct TreeIn { int Q, nth, nS, nslot; const uint32_t* prog; size_t nprog; const uint16_t* sp; size_t nsp; const uint16_t* segs; };
inline void mark_acc(std::vector<uint32_t>& tok, const uint32_t* sg, int Q) {
  bool live = false;
  for (int q = Q - 1; q >= 0; q--) {
    if (sg[18 + q] >> 31) live = true;
    for (int64_t i = (int64_t)sg[q + 1] - 1; i >= (int64_t)sg[q]; i--) {
      uint32_t& zo = tok[4 * (size_t)i]; zo = (zo & 0xFFu) | (live ? 0x100u : 0u); live = false;
      if (tok[4 * (size_t)i + 2] >> 31) live = true;
      if (tok[4 * (size_t)i + 3] >> 31) live = true;
    }
  }
}
#ifdef COFIT_X6
inline void mark_accdst(std::vector<uint32_t>& tok, const uint32_t* sg, int Q, uint32_t span, int* nprod = nullptr, int* nskip = nullptr) {
  const size_t nt = tok.size() / 4; std::vector<char> bad(nt, 0), isprod(nt, 0); bool okAll = true; int64_t p = -1; const uint32_t NONE = 0xFFFFFFFFu;
  auto ovl = [&](uint32_t a, uint32_t b) { return a < b + span && b < a + span; };
  auto verify = [&](int64_t pi, int64_t endEx, int64_t self) {
    const uint32_t dp = tok[4 * (size_t)pi + 1]; if (dp == NONE) return;
    for (int64_t t = pi + 1; t < endEx; t++) { const uint32_t dt = tok[4 * (size_t)t + 1]; if (dt != NONE && ovl(dp, dt)) { bad[(size_t)pi] = 1; return; } }
    if (self >= 0) { const uint32_t ds = tok[4 * (size_t)self + 1]; if (ds != NONE && ds != dp && ovl(dp, ds)) bad[(size_t)pi] = 1; } };
  for (int q = 0; q < Q; q++) {
    for (int64_t i = sg[q]; i < (int64_t)sg[q + 1]; i++) {
      const uint32_t wa = tok[4 * (size_t)i + 2], wb = tok[4 * (size_t)i + 3];
      if ((wa >> 31) || (wb >> 31)) { if (p < 0) okAll = false; else verify(p, i, i); }
      if ((tok[4 * (size_t)i] >> 8) & 1u) { p = i; isprod[(size_t)i] = 1; }
    }
    if (sg[18 + q] >> 31) { if (p < 0) okAll = false; else verify(p, (int64_t)sg[q + 1], -1); }
  }
  for (size_t i = 0; i < nt; i++) if (isprod[i] && tok[4 * i + 1] != NONE) { if (nprod) (*nprod)++; if (okAll && !bad[i]) { tok[4 * i] |= 0x200u; if (nskip) (*nskip)++; } }
}
#endif
inline bool g_cone_full = false;
inline void cone_masks(Tree& t, bool) {
  const int nth = t.nth, ND = nth + 1, C = 1 + ND, nS = t.nS, nsc2 = nS + 2; const uint32_t SB = t.SBd; const uint32_t SS = 48u * (uint32_t)C;
  std::vector<uint16_t> scal((size_t)nS + 2, 0); for (int j = 0; j < nS; j++) { const uint16_t w = t.sp[(size_t)j]; if ((w & 63u) == 1u) scal[(size_t)j] = (uint16_t)(1u << (((w >> 6) & 31u) + 1)); }
  std::vector<uint16_t> slotm((size_t)t.nslot + 1, 0); uint16_t acc = 0; const size_t nt = t.tokD.size() / 4; t.mskD.assign(nt, 0); t.mskA.assign(nt, 0); t.mskB.assign(nt, 0);
  auto opm = [&](uint32_t w) -> uint16_t { if (w >> 31) return acc; const uint32_t of = w & 0x1FFFFFFFu; const bool vec = ((w >> 30) & 1u) == 0u;
    if ((w >> 29) & 1u) return (uint16_t)(1u << (C - 1));
    if (!vec) return of >= SB && of - SB < (uint32_t)nS + 2 ? scal[of - SB] : (uint16_t)0;
    const uint32_t sl = of / SS; return sl < slotm.size() ? slotm[sl] : (uint16_t)0; };
  (void)nsc2;
  for (int q = 0; q < t.Q; q++) {
    for (uint32_t i = t.sgD[q]; i < t.sgD[q + 1]; i++) { const uint32_t dst = t.tokD[4 * i + 1]; const uint16_t ma = opm(t.tokD[4 * i + 2]), mb = opm(t.tokD[4 * i + 3]); const uint16_t m = ma | mb;
      t.mskD[i] = m; t.mskA[i] = ma; t.mskB[i] = mb; acc = m; if (dst != 0xFFFFFFFFu) { const uint32_t sl = dst / SS; if (sl < slotm.size()) slotm[sl] = m; } }
    t.mskOut[q] = opm(t.sgD[18 + q]);
  }
  if (g_cone_full)
  { const uint16_t full = (uint16_t)(((1u << C) - 1u) & ~1u); for (size_t i = 0; i < nt; i++) { t.mskD[i] = full; t.mskA[i] = full; t.mskB[i] = full; } for (int q = 0; q < t.Q; q++) t.mskOut[q] = full; }
  double num = 0, den = 0; for (size_t i = 0; i < nt; i++) { num += 1 + __builtin_popcount(t.mskD[i]); den += C; } t.coneFrac = den > 0 ? num / den : 1.0;
}
inline const bool g_lf8h = getenv("COFIT_F8H") && atoi(getenv("COFIT_F8H")) == 1;
static inline uint32_t lf8h_w(uint32_t v, uint32_t cd, uint32_t k) { return (v >> 31) ? v : (((v >> 30) & 1u) ? v + k : (((v >> 29) & 1u) ? v : v * cd)); }
inline bool build_tree(Tree& t, const TreeIn& in) {
  t.Q = in.Q; t.nth = in.nth; t.nS = in.nS; t.nslot = in.nslot; t.sp.assign(in.sp, in.sp + in.nsp);
  struct TT { int nth, nS, nslot; } tt{in.nth, in.nS, in.nslot};
  const int ND = in.nth + 1, C = 1 + ND;
  uint32_t hwP, hwD;
  if (g_lf8h) { hwD = 0;
    t.ok = fjn::enc_tree(tt, 48u, std::make_pair(in.prog, in.nprog), in.segs, t.tokP, t.sgP, hwP, t.SBp);
    if (t.ok) { const uint32_t cd = (uint32_t)C; t.SBd = t.SBp * cd; t.ok = t.SBd < 65536u;
      if (t.ok) { const uint32_t k8 = t.SBd - t.SBp; const size_t nw = t.tokP.size(); t.tokD.resize(nw); const uint32_t* s_ = t.tokP.data(); uint32_t* d_ = t.tokD.data();
        for (size_t i = 0; i < nw; i += 4) { d_[i] = s_[i]; d_[i + 1] = s_[i + 1] == 0xFFFFFFFFu ? s_[i + 1] : s_[i + 1] * cd; d_[i + 2] = lf8h_w(s_[i + 2], cd, k8); d_[i + 3] = lf8h_w(s_[i + 3], cd, k8); }
        for (int q = 0; q < 18; q++) t.sgD[q] = t.sgP[q]; for (int q = 18; q < 36; q++) t.sgD[q] = lf8h_w(t.sgP[q], cd, k8); } } }
  else
  t.ok = fjn::enc_tree(tt, 48u, std::make_pair(in.prog, in.nprog), in.segs, t.tokP, t.sgP, hwP, t.SBp) &&
         fjn::enc_tree(tt, 48u * (uint32_t)C, std::make_pair(in.prog, in.nprog), in.segs, t.tokD, t.sgD, hwD, t.SBd);
  if (t.ok) { mark_acc(t.tokP, t.sgP, in.Q); mark_acc(t.tokD, t.sgD, in.Q);
#ifdef COFIT_X6
    mark_accdst(t.tokP, t.sgP, in.Q, 48u); mark_accdst(t.tokD, t.sgD, in.Q, 48u * (uint32_t)C);
#endif
    cone_masks(t, hwD > 0);
#ifdef COFIT_VREUSE
    { const size_t nt = t.tokP.size() / 4; t.vro.assign(nt, -1); int off = 0;
      for (size_t i = 0; i < nt; i++) { const uint32_t zw = t.tokP[4 * i], zo = zw & 0xFFu; const bool dead = !((zw >> 8) & 1u) && t.tokP[4 * i + 1] == 0xFFFFFFFFu; if (dead) continue;
        if (zo == 7 || zo == 8) { t.vro[i] = off; off += 96; } else if (zo == 9 || zo == 10 || zo == 14) { t.vro[i] = off; off += 48; } }
      t.nvr = off; }
#endif
  }
  return t.ok && C <= MAXC;
}
struct Out { float Fx = INFINITY, Fy = INFINITY, stalled = 0, it = 0, trials = 0, dual = 0, plain = 0, solves = 0; float th[2][32]; float a[2][16]; float u[10]; float traj[14]; float res48[96]; int ok = 0; };
struct Ctx { const float* yv96; float flr; bool closed; int maxit; const fjn::WarpMat* WM; };
static inline int ppos(int p) { const int i = p / 3, c = p % 3; return 12 * (i >> 2) + 4 * c + (i & 3); }
#ifdef COFIT_RI_TR
#define COFIT_RI_SC
#define COFIT_RI_EX
#define COFIT_RI_LG
#endif
template <bool CHK>
static inline bool vsincos_t(float32x4_t x, float32x4_t& s, float32x4_t& c) {
#ifdef COFIT_SCALAR_LIBM
  return false;
#endif
#ifdef COFIT_RI_SC
  if constexpr (CHK) { if (!(vmaxvq_f32(vabsq_f32(x)) < 8192.0f)) return false; }
  const float32x4_t magic = vdupq_n_f32(12582912.0f), tq = vfmaq_f32(magic, x, vdupq_n_f32(0.63661977236758134f)), q = vsubq_f32(tq, magic);
  const int32x4_t qi = vreinterpretq_s32_f32(tq);
#else
  if constexpr (CHK) { const uint32x4_t ok = COFIT_ABSLT(x, vdupq_n_f32(8192.0f)); if (vminvq_u32(ok) == 0) return false; }
  const float32x4_t q = vrndnq_f32(vmulq_f32(x, vdupq_n_f32(0.63661977236758134f)));
#endif
  float32x4_t r = vfmsq_f32(x, q, vdupq_n_f32(1.5707963705062866f));
  r = vfmsq_f32(r, q, vdupq_n_f32(-4.371138828673793e-08f));
  r = vfmsq_f32(r, q, vdupq_n_f32(-1.7151245100058819e-15f));
  const float32x4_t z = vmulq_f32(r, r);
  float32x4_t ps = vfmaq_f32(vdupq_n_f32(8.3321608736e-3f), vdupq_n_f32(-1.9515295891e-4f), z);
  ps = vfmaq_f32(vdupq_n_f32(-1.6666654611e-1f), ps, z);
  const float32x4_t sr = vfmaq_f32(r, vmulq_f32(r, z), ps);
  float32x4_t pc = vfmaq_f32(vdupq_n_f32(-1.388731625493765e-3f), vdupq_n_f32(2.443315711809948e-5f), z);
  pc = vfmaq_f32(vdupq_n_f32(4.166664568298827e-2f), pc, z);
  const float32x4_t cr = vfmaq_f32(vfmsq_f32(vdupq_n_f32(1.0f), z, vdupq_n_f32(0.5f)), vmulq_f32(z, z), pc);
#ifdef COFIT_RI_SC
  const int32x4_t t31 = vshlq_n_s32(qi, 31), t30 = vshlq_n_s32(qi, 30);
  const uint32x4_t swap = vreinterpretq_u32_s32(vshrq_n_s32(t31, 31));
  const float32x4_t sa = vbslq_f32(swap, cr, sr), ca = vbslq_f32(swap, sr, cr);
  const uint32x4_t sb = vdupq_n_u32(0x80000000u);
  const uint32x4_t sneg = vandq_u32(vreinterpretq_u32_s32(t30), sb), cneg = vandq_u32(vreinterpretq_u32_s32(veorq_s32(t30, t31)), sb);
  s = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(sa), sneg));
  c = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(ca), cneg));
#else
  const int32x4_t qi = vcvtq_s32_f32(q);
  const uint32x4_t swap = vceqq_s32(vandq_s32(qi, vdupq_n_s32(1)), vdupq_n_s32(1));
  const float32x4_t sa = vbslq_f32(swap, cr, sr), ca = vbslq_f32(swap, sr, cr);
  const uint32x4_t sneg = vceqq_s32(vandq_s32(qi, vdupq_n_s32(2)), vdupq_n_s32(2));
  const uint32x4_t cneg = vceqq_s32(vandq_s32(vaddq_s32(qi, vdupq_n_s32(1)), vdupq_n_s32(2)), vdupq_n_s32(2));
  const uint32x4_t sb = vdupq_n_u32(0x80000000u);
  s = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(sa), vandq_u32(sneg, sb)));
  c = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(ca), vandq_u32(cneg, sb)));
#endif
  return true;
}
static inline bool vsincos(float32x4_t x, float32x4_t& s, float32x4_t& c) { return vsincos_t<true>(x, s, c); }
static inline void vsincos_k(float32x4_t x, float32x4_t& s, float32x4_t& c) { vsincos_t<false>(x, s, c); }
template <bool CHK>
static inline bool vexp_t(float32x4_t x, float32x4_t& e) {
#ifdef COFIT_SCALAR_LIBM
  return false;
#endif
#ifdef COFIT_RI_EX
  if constexpr (CHK) { if (!(vmaxvq_f32(vabsq_f32(x)) < 80.0f)) return false; }
  const float32x4_t magic = vdupq_n_f32(12582912.0f), tn = vfmaq_f32(magic, x, vdupq_n_f32(1.44269504088896341f)), n = vsubq_f32(tn, magic);
#else
  if constexpr (CHK) { const uint32x4_t ok = COFIT_ABSLT(x, vdupq_n_f32(80.0f)); if (vminvq_u32(ok) == 0) return false; }
  const float32x4_t n = vrndnq_f32(vmulq_f32(x, vdupq_n_f32(1.44269504088896341f)));
#endif
  float32x4_t r = vfmsq_f32(x, n, vdupq_n_f32(0.693359375f));
  r = vfmsq_f32(r, n, vdupq_n_f32(-2.12194440e-4f));
  float32x4_t p = vdupq_n_f32(1.9875691500e-4f);
  p = vfmaq_f32(vdupq_n_f32(1.3981999507e-3f), p, r);
  p = vfmaq_f32(vdupq_n_f32(8.3334519073e-3f), p, r);
  p = vfmaq_f32(vdupq_n_f32(4.1665795894e-2f), p, r);
  p = vfmaq_f32(vdupq_n_f32(1.6666665459e-1f), p, r);
  p = vfmaq_f32(vdupq_n_f32(5.0000001201e-1f), p, r);
  const float32x4_t z = vmulq_f32(r, r);
  const float32x4_t y = vaddq_f32(vfmaq_f32(r, p, z), vdupq_n_f32(1.0f));
#ifdef COFIT_RI_EX
  e = vreinterpretq_f32_s32(vaddq_s32(vreinterpretq_s32_f32(y), vshlq_n_s32(vreinterpretq_s32_f32(tn), 23)));
#else
  const int32x4_t ni = vcvtq_s32_f32(n);
  e = vreinterpretq_f32_s32(vaddq_s32(vreinterpretq_s32_f32(y), vshlq_n_s32(ni, 23)));
#endif
  return true;
}
static inline bool vexp(float32x4_t x, float32x4_t& e) { return vexp_t<true>(x, e); }
static inline void vexp_k(float32x4_t x, float32x4_t& e) { vexp_t<false>(x, e); }
template <bool CHK>
static inline bool vlog_t(float32x4_t x, float32x4_t& l) {
#ifdef COFIT_SCALAR_LIBM
  return false;
#endif
  const uint32x4_t xb = vreinterpretq_u32_f32(x);
#ifdef COFIT_RI_LG
  if constexpr (CHK) { const uint32x4_t i0 = vsubq_u32(xb, vdupq_n_u32(0x00800000u)); if (vmaxvq_u32(i0) >= 0x7f000000u) return false; }
  const int32x4_t ii = vreinterpretq_s32_u32(vsubq_u32(xb, vdupq_n_u32(0x3f3504f3u)));
  int32x4_t e = vshrq_n_s32(ii, 23);
  float32x4_t m = vreinterpretq_f32_u32(vaddq_u32(vandq_u32(vreinterpretq_u32_s32(ii), vdupq_n_u32(0x007fffffu)), vdupq_n_u32(0x3f3504f3u)));
#else
  if constexpr (CHK) { const uint32x4_t ok = vandq_u32(vcgeq_u32(xb, vdupq_n_u32(0x00800000u)), vcltq_u32(xb, vdupq_n_u32(0x7f800000u))); if (vminvq_u32(ok) == 0) return false; }
  int32x4_t e = vsubq_s32(vreinterpretq_s32_u32(vshrq_n_u32(xb, 23)), vdupq_n_s32(127));
  float32x4_t m = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(xb, vdupq_n_u32(0x007fffffu)), vdupq_n_u32(0x3f800000u)));
  const uint32x4_t big = vcgeq_f32(m, vdupq_n_f32(1.41421356237f));
  m = vbslq_f32(big, vmulq_f32(m, vdupq_n_f32(0.5f)), m); e = vaddq_s32(e, vandq_s32(vreinterpretq_s32_u32(big), vdupq_n_s32(1)));
#endif
  const float32x4_t f = vsubq_f32(m, vdupq_n_f32(1.0f)), ef = vcvtq_f32_s32(e), z = vmulq_f32(f, f);
  float32x4_t p = vdupq_n_f32(7.0376836292e-2f);
  p = vfmaq_f32(vdupq_n_f32(-1.1514610310e-1f), p, f); p = vfmaq_f32(vdupq_n_f32(1.1676998740e-1f), p, f); p = vfmaq_f32(vdupq_n_f32(-1.2420140846e-1f), p, f);
  p = vfmaq_f32(vdupq_n_f32(1.4249322787e-1f), p, f); p = vfmaq_f32(vdupq_n_f32(-1.6668057665e-1f), p, f); p = vfmaq_f32(vdupq_n_f32(2.0000714765e-1f), p, f);
  p = vfmaq_f32(vdupq_n_f32(-2.4999993993e-1f), p, f); p = vfmaq_f32(vdupq_n_f32(3.3333331174e-1f), p, f);
  float32x4_t y = vmulq_f32(vmulq_f32(f, z), p);
  y = vfmaq_f32(y, ef, vdupq_n_f32(-2.12194440e-4f));
  y = vfmsq_f32(y, z, vdupq_n_f32(0.5f));
  l = vfmaq_f32(vaddq_f32(f, y), ef, vdupq_n_f32(0.693359375f));
  return true;
}
static inline bool vlog(float32x4_t x, float32x4_t& l) { return vlog_t<true>(x, l); }
static inline void vlog_k(float32x4_t x, float32x4_t& l) { vlog_t<false>(x, l); }
static inline bool vcbrt(float32x4_t a, float32x4_t& y) {
#ifdef COFIT_SCALAR_LIBM
  return false;
#endif
  float32x4_t l, e; if (!vlog(a, l)) return false; if (!vexp(vmulq_f32(l, vdupq_n_f32(0.33333333333f)), e)) return false;
  y = vdivq_f32(vfmaq_f32(vmulq_f32(vdupq_n_f32(2.0f), e), a, vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(e, e))), vdupq_n_f32(3.0f));
  return true;
}
static inline void vcbrt_k(float32x4_t a, float32x4_t& y) {
  float32x4_t l, e; vlog_k(a, l); vexp_k(vmulq_f32(l, vdupq_n_f32(0.33333333333f)), e);
  y = vdivq_f32(vfmaq_f32(vmulq_f32(vdupq_n_f32(2.0f), e), a, vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(e, e))), vdupq_n_f32(3.0f));
}
inline void vm_cbrt4(const float* x, float* o) { float32x4_t Y; if (!vcbrt(vld1q_f32(x), Y)) { for (int i = 0; i < 4; i++) o[i] = std::pow(x[i], 1.0f / 3.0f); return; } vst1q_f32(o, Y); }
inline void vm_sincos4(const float* x, float* s, float* c) { float32x4_t S, Cc; float32x4_t xv = vld1q_f32(x); if (!vsincos(xv, S, Cc)) { for (int i = 0; i < 4; i++) { s[i] = std::sin(x[i]); c[i] = std::cos(x[i]); } return; } vst1q_f32(s, S); vst1q_f32(c, Cc); }
inline void vm_exp4(const float* x, float* o) { float32x4_t E; if (!vexp(vld1q_f32(x), E)) { for (int i = 0; i < 4; i++) o[i] = std::exp(x[i]); return; } vst1q_f32(o, E); }
inline void vm_log4(const float* x, float* o) { float32x4_t L; if (!vlog(vld1q_f32(x), L)) { for (int i = 0; i < 4; i++) o[i] = std::log(x[i]); return; } vst1q_f32(o, L); }
static inline float dot48(const float* a, const float* b) {
  float32x4_t L[4];
  for (int g = 0; g < 4; g++) {
    const float* pa = a + 12 * g; const float* pb = b + 12 * g;
    L[g] = CF_DOT3(vld1q_f32(pa), vld1q_f32(pb), vld1q_f32(pa + 4), vld1q_f32(pb + 4), vld1q_f32(pa + 8), vld1q_f32(pb + 8));
  }
  const float32x4_t X = vpaddq_f32(L[0], L[1]), Y = vpaddq_f32(L[2], L[3]), Z = vpaddq_f32(X, Y);
  const float32x2_t W = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  return vget_lane_f32(W, 0) + vget_lane_f32(W, 1);
}
static inline float red16v(const float32x4_t* l) {
  const float32x4_t X = vpaddq_f32(l[0], l[1]), Y = vpaddq_f32(l[2], l[3]), Z = vpaddq_f32(X, Y);
  const float32x2_t W = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  return vget_lane_f32(W, 0) + vget_lane_f32(W, 1);
}
static inline float32x4_t dot_acc(const float* a, const float* b) {
  float32x4_t s0 = vmulq_f32(vld1q_f32(a), vld1q_f32(b)), s1 = vmulq_f32(vld1q_f32(a + 4), vld1q_f32(b + 4)), s2 = vmulq_f32(vld1q_f32(a + 8), vld1q_f32(b + 8));
  for (int v = 3; v < 12; v += 3) { s0 = vfmaq_f32(s0, vld1q_f32(a + 4 * v), vld1q_f32(b + 4 * v)); s1 = vfmaq_f32(s1, vld1q_f32(a + 4 * v + 4), vld1q_f32(b + 4 * v + 4)); s2 = vfmaq_f32(s2, vld1q_f32(a + 4 * v + 8), vld1q_f32(b + 4 * v + 8)); }
  return vaddq_f32(vaddq_f32(s0, s1), s2);
}
#ifdef COFIT_EXACT_GRAM
static inline float dot_fast(const float* a, const float* b) { float hs = 0; for (int p = 0; p < 48; p++) { const int q = ppos(p); hs = std::fmaf(a[q], b[q], hs); } return hs; }
#else
static inline float dot_fast(const float* a, const float* b) { return vaddvq_f32(dot_acc(a, b)); }
#endif
static inline float32x4_t dot4k(const float* a, const float* b0, const float* b1, const float* b2, const float* b3) {
  float32x4_t av[12]; for (int i = 0; i < 12; i++) av[i] = vld1q_f32(a + 4 * i);
  const float* bb[4] = {b0, b1, b2, b3}; float32x2_t W[4];
  for (int d = 0; d < 4; d++) {
    float32x4_t L[4];
    for (int g = 0; g < 4; g++) L[g] = CF_DOT3(av[3 * g], vld1q_f32(bb[d] + 12 * g), av[3 * g + 1], vld1q_f32(bb[d] + 12 * g + 4), av[3 * g + 2], vld1q_f32(bb[d] + 12 * g + 8));
    const float32x4_t X = vpaddq_f32(L[0], L[1]), Yv = vpaddq_f32(L[2], L[3]), Z = vpaddq_f32(X, Yv);
    W[d] = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  }
  return vpaddq_f32(vcombine_f32(W[0], W[1]), vcombine_f32(W[2], W[3]));
}
#ifdef COFIT_G22
template <int RA, int CB>
static inline void gram_tile(const float* a0, const float* a1, const float* b0, const float* b1, float* out) {
  float32x4_t s[2][2][3];
#pragma clang loop unroll(full)
  for (int v = 0; v < 12; v++) {
    const float32x4_t A0 = vld1q_f32(a0 + 4 * v), A1 = RA > 1 ? vld1q_f32(a1 + 4 * v) : A0, B0 = vld1q_f32(b0 + 4 * v), B1 = CB > 1 ? vld1q_f32(b1 + 4 * v) : B0;
    const int k = v % 3;
    if (v < 3) { s[0][0][k] = vmulq_f32(A0, B0); if (CB > 1) s[0][1][k] = vmulq_f32(A0, B1); if (RA > 1) { s[1][0][k] = vmulq_f32(A1, B0); if (CB > 1) s[1][1][k] = vmulq_f32(A1, B1); } }
    else { s[0][0][k] = vfmaq_f32(s[0][0][k], A0, B0); if (CB > 1) s[0][1][k] = vfmaq_f32(s[0][1][k], A0, B1); if (RA > 1) { s[1][0][k] = vfmaq_f32(s[1][0][k], A1, B0); if (CB > 1) s[1][1][k] = vfmaq_f32(s[1][1][k], A1, B1); } }
  }
  float32x4_t T[2][2];
  for (int r = 0; r < RA; r++) for (int c = 0; c < CB; c++) T[r][c] = vaddq_f32(vaddq_f32(s[r][c][0], s[r][c][1]), s[r][c][2]);
  for (int r = 0; r < RA; r++) for (int c = 0; c < CB; c++) out[r * 2 + c] = vaddvq_f32(T[r][c]);
}
#endif
struct PlainRes { float rss; float a[16]; float Y[48]; float E[16 * 48]; };
static inline void copy_plain(PlainRes& d, const PlainRes& s, int Q) { d.rss = s.rss; std::memcpy(d.a, s.a, sizeof d.a); std::memcpy(d.Y, s.Y, sizeof d.Y); std::memcpy(d.E, s.E, (size_t)Q * 48 * 4); }
struct Scratch {
  std::vector<float> Rv, Rv2;
#ifdef COFIT_VREUSE
  std::vector<float> VR[4];
#endif
  alignas(16) float outq[16 * MAXC * 48]; alignas(16) float COL[26 * 48];
  alignas(16) float ST[9 * 48]; float Skey[48 * 9]; bool stv = false;
  alignas(16) float Hd[NS * NS], Ab[NS * NS], A[NS * NS], bsel[NS * NS];
  alignas(16) float A2[NS * NS]; alignas(16) float xv2[NS + 8], dg2[NS + 8];
  alignas(16) float GV[NS + 8], DEL[NS + 8], scl[NS + 8], xv[NS + 8], dg[NS + 8], Gd[NS + 8], xs[NS + 8]; alignas(16) float FV[4 * (NS + 8)]; alignas(16) float GR[28 * 28]; alignas(16) float zero48[48];
  Scratch() { std::memset(Hd, 0, sizeof Hd); std::memset(Ab, 0, sizeof Ab); std::memset(A, 0, sizeof A); std::memset(GV, 0, sizeof GV); std::memset(scl, 0, sizeof scl); std::memset(xv, 0, sizeof xv); std::memset(FV, 0, sizeof FV); std::memset(GR, 0, sizeof GR); std::memset(zero48, 0, sizeof zero48); std::memset(DEL, 0, sizeof DEL); std::memset(dg, 0, sizeof dg); std::memset(Gd, 0, sizeof Gd); std::memset(xs, 0, sizeof xs); std::memset(bsel, 0, sizeof bsel); std::memset(A2, 0, sizeof A2); std::memset(xv2, 0, sizeof xv2); std::memset(dg2, 0, sizeof dg2); std::memset(outq, 0, sizeof outq); std::memset(COL, 0, sizeof COL); }
};
struct Opd { const float* base; int cs; bool vec; };
#ifdef COFIT_LANEFIX
template <class F>
static inline float32x4_t cf_fixlanes(float32x4_t r, float32x4_t m, uint32x4_t okm, F f) {
  float rr[4], mm[4]; uint32_t om[4]; vst1q_f32(rr, r); vst1q_f32(mm, m); vst1q_u32(om, okm);
  for (int i = 0; i < 4; i++) if (!om[i]) rr[i] = f(mm[i]);
  return vld1q_f32(rr);
}
static inline uint32x4_t cf_normal_mask(float32x4_t am) { return vcltq_u32(vsubq_u32(vreinterpretq_u32_f32(am), vdupq_n_u32(0x00800000u)), vdupq_n_u32(0x7f000000u)); }
#ifdef COFIT_CBRT2
static inline uint32x4_t cf_cbrt_mask(float32x4_t am) { return vcltq_u32(vsubq_u32(vreinterpretq_u32_f32(am), vdupq_n_u32(0x0704EC3Du)), vdupq_n_u32(0x7B4097CEu - 0x0704EC3Du)); }
static inline float32x4_t cf_cbrt_fast(float32x4_t a) {
  const float32x4_t fi = vcvtq_f32_u32(vreinterpretq_u32_f32(a));
  float32x4_t y = vreinterpretq_f32_u32(vaddq_u32(vcvtq_u32_f32(vmulq_f32(fi, vdupq_n_f32(0.33333333f))), vdupq_n_u32(0x2A5137A0u)));
  const float32x4_t two = vdupq_n_f32(2.0f);
  const float32x4_t y3 = vmulq_f32(vmulq_f32(y, y), y);
  y = vmulq_f32(y, vdivq_f32(vfmaq_f32(y3, two, a), vfmaq_f32(a, two, y3)));
  const float32x4_t y2 = vmulq_f32(y, y);
  const float32x4_t r = vfmsq_f32(a, y2, y), den = vfmaq_f32(a, two, vmulq_f32(y2, y));
  return vfmaq_f32(y, y, vdivq_f32(r, den));
}
#endif
#endif
#define J12 _Pragma("clang loop unroll(full)") for (int j = 0; j < 12; j++)
#ifdef COFIT_VREUSE
#ifdef COFIT_H5_5
#error COFIT_VREUSE is built for the non-H5_5 exec_tok signature
#endif
template <bool AV, bool BV, bool ACC, bool DUAL = true, int VRM = 0>
#else
template <bool AV, bool BV, bool ACC, bool DUAL = true>
#endif
#ifdef COFIT_H5_5
static inline void exec_tok(uint32_t zo, const float* pa, const float* pb, uint64_t cs2, float* acc, float* dst, uint64_t m2, float32x4_t* nfp_in) {
  const int csa = (int)(uint32_t)cs2, csb = (int)(uint32_t)(cs2 >> 32); const uint32_t mA = DUAL ? (uint32_t)m2 : 0u, mB = DUAL ? (uint32_t)(m2 >> 32) : 0u;
  float32x4_t* const nfp = DUAL ? nfp_in : nullptr;
  float32x4_t a0[12], b0[12];
  { const float32x4_t sa = vdupq_n_f32(pa[0]); J12 { a0[j] = AV ? vld1q_f32(pa + 4 * j) : sa; }
    if (zo != 6) { const float32x4_t sb = vdupq_n_f32(pb[0]); J12 { b0[j] = BV ? vld1q_f32(pb + 4 * j) : sb; } } }
#else
#ifdef COFIT_VREUSE
static inline void exec_tok(uint32_t zo, int C, const float* pa, int csa, const float* pb, int csb, float* acc, float* dst, uint32_t mA, uint32_t mB, float32x4_t* nfp, float* svp = nullptr) {
#else
static inline void exec_tok(uint32_t zo, int C, const float* pa, int csa, const float* pb, int csb, float* acc, float* dst, uint32_t mA, uint32_t mB, float32x4_t* nfp) {
#endif
  (void)C;
  float32x4_t a0[12], b0[12];
  { const float32x4_t sa = vdupq_n_f32(pa[0]), sb = vdupq_n_f32(pb[0]); J12 { a0[j] = AV ? vld1q_f32(pa + 4 * j) : sa; b0[j] = BV ? vld1q_f32(pb + 4 * j) : sb; } }
#endif
  const float32x4_t zz = vdupq_n_f32(0.0f); float32x4_t nfl = zz;
#define LDA(c, j) (HA ? (AV ? vld1q_f32(pa + (c) * csa + 4 * (j)) : sac) : zz)
#define LDB(c, j) (HB ? (BV ? vld1q_f32(pb + (c) * csb + 4 * (j)) : sbc) : zz)
#define ST(c, j, val) { const float32x4_t _v = (val); if (ACC) vst1q_f32(acc + (c) * 48 + 4 * (j), _v); if (dst) vst1q_f32(dst + (c) * 48 + 4 * (j), _v); }
#define ST0(j, val) { const float32x4_t _v = (val); if (nfp) nfl = vfmaq_f32(nfl, _v, zz); if (ACC) vst1q_f32(acc + 4 * (j), _v); if (dst) vst1q_f32(dst + 4 * (j), _v); }
#define FOR_TANG(...) for (uint32_t um_ = (mA | mB) & ~1u; um_; um_ &= um_ - 1) { const int c = __builtin_ctz(um_); const bool ha_ = (mA >> c) & 1u, hb_ = (mB >> c) & 1u; \
    const float32x4_t sac = (AV || !ha_) ? zz : vdupq_n_f32(pa[c * csa]); const float32x4_t sbc = (BV || !hb_) ? zz : vdupq_n_f32(pb[c * csb]); \
    if (ha_ && hb_) { constexpr bool HA = true, HB = true; __VA_ARGS__ } else if (ha_) { constexpr bool HA = true, HB = false; __VA_ARGS__ } else { constexpr bool HA = false, HB = true; __VA_ARGS__ } }
  switch (zo) {
    case 0: {
      FOR_TANG(J12 { const float32x4_t A = LDA(c, j), B = LDB(c, j); ST(c, j, HA ? (HB ? vfmaq_f32(vmulq_f32(a0[j], B), A, b0[j]) : vmulq_f32(A, b0[j])) : vmulq_f32(a0[j], B)); })
      J12 { ST0(j, vmulq_f32(a0[j], b0[j])); } } break;
    case 1: J12 { ST0(j, vaddq_f32(a0[j], b0[j])); }
      FOR_TANG(J12 { ST(c, j, HA ? (HB ? vaddq_f32(LDA(c, j), LDB(c, j)) : LDA(c, j)) : LDB(c, j)); }) break;
    case 2: J12 { ST0(j, vsubq_f32(a0[j], b0[j])); }
      FOR_TANG(J12 { ST(c, j, HA ? (HB ? vsubq_f32(LDA(c, j), LDB(c, j)) : LDA(c, j)) : vnegq_f32(LDB(c, j))); }) break;
    case 3: J12 { ST0(j, vnegq_f32(vsubq_f32(a0[j], b0[j]))); }
      FOR_TANG(J12 { ST(c, j, HA ? (HB ? vnegq_f32(vsubq_f32(LDA(c, j), LDB(c, j))) : vnegq_f32(LDA(c, j))) : LDB(c, j)); }) break;
    case 4: {
      float32x4_t v[12]; J12 { v[j] = vdivq_f32(a0[j], b0[j]); }
      FOR_TANG(J12 { ST(c, j, vdivq_f32(HA ? (HB ? CF_SM(LDA(c, j), v[j], LDB(c, j)) : LDA(c, j)) : vnegq_f32(vmulq_f32(v[j], LDB(c, j))), b0[j])); })
      J12 { ST0(j, v[j]); } } break;
    case 5: {
      float32x4_t v[12]; J12 { v[j] = vdivq_f32(b0[j], a0[j]); }
      FOR_TANG(J12 { ST(c, j, vdivq_f32(HB ? (HA ? CF_SM(LDB(c, j), v[j], LDA(c, j)) : LDB(c, j)) : vnegq_f32(vmulq_f32(v[j], LDA(c, j))), a0[j])); })
      J12 { ST0(j, v[j]); } } break;
    case 6: J12 { ST0(j, a0[j]); }
      FOR_TANG(J12 { ST(c, j, LDA(c, j)); }) break;
    default: {
      float32x4_t v[12], d[12]; bool tg_done = false;
#ifdef COFIT_VREUSE
      if constexpr (VRM == 2) {
        if (svp) {
          if (zo == 7 || zo == 8) { J12 { v[j] = vld1q_f32(svp + 4 * j); d[j] = vld1q_f32(svp + 48 + 4 * j); } tg_done = true; }
          else if (zo == 9) { J12 { v[j] = vld1q_f32(svp + 4 * j); d[j] = v[j]; } tg_done = true; }
          else if (zo == 10) { J12 { v[j] = vld1q_f32(svp + 4 * j); d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(a0[j], b0[j])); } tg_done = true; }
          else if (zo == 14) { J12 { v[j] = vld1q_f32(svp + 4 * j); d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), v[j]), v[j])); } tg_done = true; } } }
#endif
#ifdef COFIT_RI_TG
#ifndef COFIT_RI_TGM
#define COFIT_RI_TGM 0xFF80
#endif
      if (!tg_done && zo >= 11 && zo <= 15 && ((COFIT_RI_TGM >> zo) & 1)) {
        bool okc = true;
#ifdef COFIT_CBRT2
        if (zo == 14) { uint32x4_t mxu = vdupq_n_u32(0); J12 { mxu = vmaxq_u32(mxu, vsubq_u32(vreinterpretq_u32_f32(vabsq_f32(vmulq_f32(a0[j], b0[j]))), vdupq_n_u32(0x0704EC3Du))); } okc = vmaxvq_u32(mxu) < 0x7B4097CEu - 0x0704EC3Du; }
#else
        if (zo == 14) { uint32x4_t mxu = vdupq_n_u32(0); J12 { mxu = vmaxq_u32(mxu, vsubq_u32(vreinterpretq_u32_f32(vabsq_f32(vmulq_f32(a0[j], b0[j]))), vdupq_n_u32(0x00800000u))); } okc = vmaxvq_u32(mxu) < 0x7f000000u; }
#endif
        if (okc) {
          if (zo == 11) { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); v[j] = vmulq_f32(m, m); d[j] = vmulq_f32(vdupq_n_f32(2.0f), m); } }
          else if (zo == 12) { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); v[j] = vmulq_f32(vmulq_f32(m, m), m); d[j] = vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), m), m); } }
          else if (zo == 13) { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); const float32x4_t s = vsqrtq_f32(vabsq_f32(m)); v[j] = s; const float32x4_t half = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(vdupq_n_f32(0.5f)))); if (mA | mB) d[j] = vdivq_f32(half, s); } }
#ifdef COFIT_CBRT2
          else if (zo == 14) { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); const float32x4_t y = cf_cbrt_fast(vabsq_f32(m));
            const float32x4_t vv = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(y))); v[j] = vv; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vv), vv)); } }
#endif
          else if (zo == 14) { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); float32x4_t y; vcbrt_k(vabsq_f32(m), y); const float32x4_t vv = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(y))); v[j] = vv; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vv), vv)); } }
          else { J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]); v[j] = vnegq_f32(m); d[j] = vdupq_n_f32(-1.0f); } }
          tg_done = true; } }
      if (!tg_done && (zo == 7 || zo == 8 || zo == 9 || zo == 10) && ((COFIT_RI_TGM >> zo) & 1)) {
        float32x4_t mv[12]; float32x4_t mxf = vdupq_n_f32(0.0f); uint32x4_t mxu = vdupq_n_u32(0);
        if (zo == 10) { J12 { mv[j] = vabsq_f32(vmulq_f32(a0[j], b0[j])); mxu = vmaxq_u32(mxu, vsubq_u32(vreinterpretq_u32_f32(mv[j]), vdupq_n_u32(0x00800000u))); } }
        else { J12 { mv[j] = vmulq_f32(a0[j], b0[j]); mxf = vmaxq_f32(mxf, vabsq_f32(mv[j])); } }
        const bool allok = zo == 10 ? vmaxvq_u32(mxu) < 0x7f000000u : vmaxvq_f32(mxf) < (zo == 9 ? 80.0f : 8192.0f);
        if (allok) {
          if (zo == 7) { J12 { float32x4_t s_, c_; vsincos_k(mv[j], s_, c_); v[j] = s_; d[j] = c_; } }
          else if (zo == 8) { J12 { float32x4_t s_, c_; vsincos_k(mv[j], s_, c_); v[j] = c_; d[j] = vnegq_f32(s_); } }
          else if (zo == 9) { J12 { float32x4_t e_; vexp_k(mv[j], e_); v[j] = e_; d[j] = e_; } }
          else { J12 { float32x4_t l_; vlog_k(mv[j], l_); v[j] = l_; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(a0[j], b0[j])); } }
          tg_done = true; } }
#endif
#ifdef COFIT_PROF2
      if (!tg_done) g_fbt[(mA | mB) ? 1 : 0][zo]++;
#endif
      if (!tg_done) J12 { const float32x4_t m = vmulq_f32(a0[j], b0[j]);
#ifdef COFIT_LANEFIX
        if (zo == 7 || zo == 8) { const uint32x4_t okm = COFIT_ABSLT(m, vdupq_n_f32(8192.0f)); float32x4_t s, c_; vsincos_k(vbslq_f32(okm, m, zz), s, c_);
          if (vminvq_u32(okm) == 0) { float mm[4], sv[4], cv[4]; uint32_t om[4]; vst1q_f32(mm, m); vst1q_f32(sv, s); vst1q_f32(cv, c_); vst1q_u32(om, okm);
            for (int i = 0; i < 4; i++) if (!om[i]) { sv[i] = std::sin(mm[i]); cv[i] = std::cos(mm[i]); } s = vld1q_f32(sv); c_ = vld1q_f32(cv); }
          if (zo == 7) { v[j] = s; d[j] = c_; } else { v[j] = c_; d[j] = vnegq_f32(s); } continue; }
        if (zo == 9) { const uint32x4_t okm = COFIT_ABSLT(m, vdupq_n_f32(80.0f)); float32x4_t e; vexp_k(vbslq_f32(okm, m, zz), e);
          if (vminvq_u32(okm) == 0) e = cf_fixlanes(e, m, okm, [](float x) { return std::exp(x); });
          v[j] = e; d[j] = e; continue; }
        if (zo == 10) { const float32x4_t am = vabsq_f32(m); const uint32x4_t okm = cf_normal_mask(am); float32x4_t l; vlog_k(vbslq_f32(okm, am, vdupq_n_f32(1.0f)), l);
          if (vminvq_u32(okm) == 0) l = cf_fixlanes(l, am, okm, [](float x) { return x == 0.0f ? -INFINITY : (x > 3.4028235e38f || x != x) ? x : std::log(x); });
          v[j] = l; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), m); continue; }
#ifdef COFIT_CBRT2
        if (zo == 14) { const float32x4_t am = vabsq_f32(m); const uint32x4_t okf = cf_cbrt_mask(am); float32x4_t y = cf_cbrt_fast(vbslq_f32(okf, am, vdupq_n_f32(1.0f))); uint32x4_t okm = okf;
          if (vminvq_u32(okf) == 0) { okm = cf_normal_mask(am); if (vmaxvq_u32(vbicq_u32(okm, okf))) { float32x4_t yo; vcbrt_k(vbslq_f32(okm, am, vdupq_n_f32(1.0f)), yo); y = vbslq_f32(okf, y, yo); } }
          float32x4_t vv = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(y)));
          if (vminvq_u32(okm) == 0) vv = cf_fixlanes(vv, m, okm, [](float x) { const float ax = std::fabs(x); return (ax == 0.0f || !(ax < INFINITY)) ? x : std::copysign(std::pow(ax, 1.0f / 3.0f), x); });
          v[j] = vv; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vv), vv)); continue; }
#endif
        if (zo == 14) { const float32x4_t am = vabsq_f32(m); const uint32x4_t okm = cf_normal_mask(am); float32x4_t y; vcbrt_k(vbslq_f32(okm, am, vdupq_n_f32(1.0f)), y);
          float32x4_t vv = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(y)));
          if (vminvq_u32(okm) == 0) vv = cf_fixlanes(vv, m, okm, [](float x) { const float ax = std::fabs(x); return (ax == 0.0f || !(ax < INFINITY)) ? x : std::copysign(std::pow(ax, 1.0f / 3.0f), x); });
          v[j] = vv; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vv), vv)); continue; }
#endif
        switch (zo) {
          case 7: case 8: { float32x4_t s, c_; if (vsincos(m, s, c_)) { if (zo == 7) { v[j] = s; d[j] = c_; } else { v[j] = c_; d[j] = vnegq_f32(s); } } else { float mm[4], vv[4], dd[4]; vst1q_f32(mm, m); for (int i = 0; i < 4; i++) { vv[i] = zo == 7 ? std::sin(mm[i]) : std::cos(mm[i]); dd[i] = zo == 7 ? std::cos(mm[i]) : -std::sin(mm[i]); } v[j] = vld1q_f32(vv); d[j] = vld1q_f32(dd); } } break;
          case 9: { float32x4_t e; if (!vexp(m, e)) { float mm[4], vv[4]; vst1q_f32(mm, m); for (int i = 0; i < 4; i++) vv[i] = std::exp(mm[i]); e = vld1q_f32(vv); } v[j] = e; d[j] = e; } break;
          case 10: { float32x4_t l; const float32x4_t am = vabsq_f32(m); if (!vlog(am, l)) { float mm[4], vv[4]; vst1q_f32(mm, am); for (int i = 0; i < 4; i++) vv[i] = std::log(mm[i]); l = vld1q_f32(vv); } v[j] = l; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), m); } break;
          case 11: v[j] = vmulq_f32(m, m); d[j] = vmulq_f32(vdupq_n_f32(2.0f), m); break;
          case 12: v[j] = vmulq_f32(vmulq_f32(m, m), m); d[j] = vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), m), m); break;
          case 13: { const float32x4_t s = vsqrtq_f32(vabsq_f32(m)); v[j] = s; const float32x4_t half = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(vdupq_n_f32(0.5f)))); if (mA | mB) d[j] = vdivq_f32(half, s); } break;
          case 14: { float32x4_t y; const float32x4_t am = vabsq_f32(m);
            if (vcbrt(am, y)) { const float32x4_t vv = vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(m), vdupq_n_u32(0x80000000u)), vreinterpretq_u32_f32(y))); v[j] = vv; if (mA | mB) d[j] = vdivq_f32(vdupq_n_f32(1.0f), vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vv), vv)); break; }
            float mm[4], vv[4], dd[4]; vst1q_f32(mm, m); for (int i = 0; i < 4; i++) { vv[i] = std::copysign(std::pow(std::fabs(mm[i]), 1.0f / 3.0f), mm[i]); dd[i] = 1.0f / (3.0f * vv[i] * vv[i]); } v[j] = vld1q_f32(vv); d[j] = vld1q_f32(dd); } break;
          default: v[j] = vnegq_f32(m); d[j] = vdupq_n_f32(-1.0f); break;
        } }
#ifdef COFIT_VREUSE
      if constexpr (VRM == 1) { if (svp) { J12 { vst1q_f32(svp + 4 * j, v[j]); } if (zo == 7 || zo == 8) { J12 { vst1q_f32(svp + 48 + 4 * j, d[j]); } } } }
#endif
      if (nfp) J12 { nfl = vfmaq_f32(nfl, d[j], zz); }
      FOR_TANG(J12 { ST(c, j, vmulq_f32(d[j], HA ? (HB ? vfmaq_f32(vmulq_f32(a0[j], LDB(c, j)), LDA(c, j), b0[j]) : vmulq_f32(LDA(c, j), b0[j])) : vmulq_f32(a0[j], LDB(c, j)))); })
      J12 { ST0(j, v[j]); } } break;
  }
#undef LDA
#undef LDB
#undef ST
#undef ST0
#undef FOR_TANG
  if (nfp) *nfp = vaddq_f32(*nfp, nfl);
}
#ifdef COFIT_VR_CNT
struct VrCnt { std::atomic<uint64_t> dualtok[16], reused[16], plaintok2[16], pairs, duals, dualsreuse; VrCnt() { for (int i = 0; i < 16; i++) dualtok[i] = reused[i] = plaintok2[i] = 0; pairs = duals = dualsreuse = 0; }
  ~VrCnt() { fprintf(stderr, "VRCNT pairs %llu dualwalks %llu dualwalks_with_reuse %llu\n", (unsigned long long)pairs, (unsigned long long)duals, (unsigned long long)dualsreuse);
    for (int z = 0; z < 16; z++) if (dualtok[z] || plaintok2[z]) fprintf(stderr, "VRCNT op %d dual_tokens %llu dual_tokens_reused %llu plain_cand_pair_tokens %llu\n", z, (unsigned long long)dualtok[z], (unsigned long long)reused[z], (unsigned long long)plaintok2[z]); } };
inline VrCnt g_vrc;
#endif
#ifdef COFIT_VREUSE
static bool run_tree(const Tree& T, bool dual, const float* th, const float* tsh, int C, Scratch& S, float* outq, bool cone = true, float* vrsv = nullptr) {
#else
static bool run_tree(const Tree& T, bool dual, const float* th, const float* tsh, int C, Scratch& S, float* outq, bool cone = true) {
#endif
  PHASE(16);
  const uint32_t* tok = dual ? T.tokD.data() : T.tokP.data(); const uint32_t* sg = dual ? T.sgD : T.sgP; const uint32_t SB = dual ? T.SBd : T.SBp;
  const int nS = T.nS, nsc2 = nS + 2; { const size_t need = (size_t)SB + (size_t)C * nsc2 + 8; if (S.Rv.size() < need) S.Rv.resize(need); } float* R = S.Rv.data(); alignas(16) float acc[MAXC * 48];
  std::memset(R + SB, 0, (size_t)C * nsc2 * 4);
  for (int j = 0; j < nS; j++) { const uint16_t w = T.sp[(size_t)j]; const uint32_t op = w & 63u, a = (w >> 6) & 31u; R[SB + j] = (op == 1u) ? th[a] : 1.0f;
    if (op == 1u && (int)a + 1 < C) R[SB + ((int)a + 1) * nsc2 + j] = 1.0f; }
  R[SB + nS] = 1.0f; R[SB + nS + 1] = 2.0f;
  std::memset(acc, 0, (size_t)C * 48 * 4);
#ifdef COFIT_CONE
  const bool lim = dual && cone;
#else
  const bool lim = false; (void)cone;
#endif
#ifdef COFIT_CONE_NOGUARD
  constexpr bool NFCHK_ = false;
#else
  const bool NFCHK_ = lim;
#endif
  const uint32_t fullm = (uint32_t)(((1u << C) - 1u) & ~1u); float32x4_t nfv = vdupq_n_f32(0.0f);
#ifdef COFIT_X6
  const float* accp = acc;
  auto fetch = [&](uint32_t w) -> Opd {
    if (w >> 31) return Opd{accp, 48, true};
#else
  auto fetch = [&](uint32_t w) -> Opd {
    if (w >> 31) return Opd{acc, 48, true};
#endif
    const uint32_t of = w & 0x1FFFFFFFu; const bool vec = ((w >> 30) & 1u) == 0u; const float* bp = ((w >> 29) & 1u) ? tsh : R;
    return Opd{bp + of, vec ? 48 : nsc2, vec}; };
  for (int q = 0; q < T.Q; q++) {
    PHASE(dual ? 4 : 1);
    for (uint32_t i = sg[q]; i < sg[q + 1]; i++) {
      const uint32_t zw = tok[4 * i], zo = zw & 0xFFu, accn = (zw >> 8) & 1u, dst = tok[4 * i + 1], wa = tok[4 * i + 2], wb = tok[4 * i + 3];
      if (!accn && dst == 0xFFFFFFFFu) continue;
#ifdef COFIT_PROF2
      g_tokn[dual ? 1 : 0][zo > 15 ? 15 : zo]++;
#ifdef COFIT_PROF_OPT
      const uint64_t _ot0 = mach_absolute_time();
#endif
#endif
#ifdef COFIT_PROF_OPS
      const uint64_t _t0 = mach_absolute_time();
#endif
#ifdef COFIT_VREUSE
      float* svt = (dual && vrsv && T.vro[i] >= 0) ? vrsv + T.vro[i] : nullptr;
#ifdef COFIT_VR_CNT
      if (dual) { g_vrc.dualtok[zo > 15 ? 15 : zo]++; if (svt) g_vrc.reused[zo > 15 ? 15 : zo]++; }
#endif
#endif
      const Opd oa = fetch(wa), ob = fetch(wb); float* dp = (dst != 0xFFFFFFFFu) ? R + dst : nullptr; const uint32_t mA = !dual ? 0u : (lim ? (uint32_t)T.mskA[i] : fullm), mB = !dual ? 0u : (lim ? (uint32_t)T.mskB[i] : fullm);
#ifdef COFIT_X6
      const bool skAcc = (zw >> 9) & 1u;
      switch ((oa.vec ? 4 : 0) | (ob.vec ? 2 : 0) | ((accn && !skAcc) ? 1 : 0)) {
#else
      switch ((oa.vec ? 4 : 0) | (ob.vec ? 2 : 0) | (accn ? 1 : 0)) {
#endif
#ifdef COFIT_H5_5
#ifdef COFIT_LEANP
#define EX(i, A_, B_, K_) case i: if (dual) exec_tok<A_, B_, K_, true>(zo, oa.base, ob.base, (uint64_t)(uint32_t)oa.cs | ((uint64_t)(uint32_t)ob.cs << 32), acc, dp, (uint64_t)mA | ((uint64_t)mB << 32), NFCHK_ ? &nfv : nullptr); \
    else exec_tok<A_, B_, K_, CF_LEANF>(zo, oa.base, ob.base, (uint64_t)(uint32_t)oa.cs | ((uint64_t)(uint32_t)ob.cs << 32), acc, dp, 0ull, nullptr); break;
#else
#define EX(i, A_, B_, K_) case i: exec_tok<A_, B_, K_>(zo, oa.base, ob.base, (uint64_t)(uint32_t)oa.cs | ((uint64_t)(uint32_t)ob.cs << 32), acc, dp, (uint64_t)mA | ((uint64_t)mB << 32), NFCHK_ ? &nfv : nullptr); break;
#endif
#else
#ifdef COFIT_VREUSE
#define EX(i, A_, B_, K_) case i: if (svt) exec_tok<A_, B_, K_, true, 2>(zo, C, oa.base, oa.cs, ob.base, ob.cs, acc, dp, mA, mB, NFCHK_ ? &nfv : nullptr, svt); else exec_tok<A_, B_, K_>(zo, C, oa.base, oa.cs, ob.base, ob.cs, acc, dp, mA, mB, NFCHK_ ? &nfv : nullptr, nullptr); break;
#else
#define EX(i, A_, B_, K_) case i: exec_tok<A_, B_, K_>(zo, C, oa.base, oa.cs, ob.base, ob.cs, acc, dp, mA, mB, NFCHK_ ? &nfv : nullptr); break;
#endif
#endif
        EX(0, false, false, false) EX(1, false, false, true) EX(2, false, true, false) EX(3, false, true, true) EX(4, true, false, false) EX(5, true, false, true) EX(6, true, true, false) EX(7, true, true, true)
#undef EX
      }
#ifdef COFIT_X6
      if (accn) accp = skAcc ? (const float*)dp : (const float*)acc;
#endif
#if defined(COFIT_PROF2) && defined(COFIT_PROF_OPT)
      g_tokt[dual ? 1 : 0][zo > 15 ? 15 : zo] += mach_absolute_time() - _ot0;
#endif
#ifdef COFIT_PROF_OPS
      g_op[dual?1:0][zo>15?15:zo] += (double)(mach_absolute_time()-_t0); g_opn[dual?1:0][zo>15?15:zo] += 1;
#endif
    }
    PHASE(17);
    const Opd ox = fetch(sg[18 + q]);
    const uint32_t mo = !dual ? 0u : (lim ? (uint32_t)T.mskOut[q] : fullm);
    for (int c = 0; c < C; c++) {
#ifdef COFIT_X2
      if (dual && c == 0) continue;
#endif
      if (c > 0 && dual && !((mo >> c) & 1u)) continue; float* o = outq + (size_t)q * C * 48 + c * 48; if (ox.vec) std::memcpy(o, ox.base + c * ox.cs, 48 * 4); else { const float32x4_t b = vdupq_n_f32(ox.base[c * ox.cs]); for (int j = 0; j < 12; j++) vst1q_f32(o + 4 * j, b); } }
  }
  if (NFCHK_ && vmaxvq_u32(vmvnq_u32(vceqq_f32(nfv, nfv))) != 0) return false;
  return true;
}
#ifdef COFIT_VREUSE
static void run_tree2(const Tree& T, const float* th0, const float* th1, const float* ts0, const float* ts1, Scratch& S, float* out0, float* out1, float* vr0 = nullptr, float* vr1 = nullptr) {
#else
static void run_tree2(const Tree& T, const float* th0, const float* th1, const float* ts0, const float* ts1, Scratch& S, float* out0, float* out1) {
#endif
  const uint32_t* tok = T.tokP.data(); const uint32_t* sg = T.sgP; const uint32_t SB = T.SBp;
  const int nS = T.nS, nsc2 = nS + 2; const size_t need = (size_t)SB + (size_t)nsc2 + 8; if (S.Rv.size() < need) S.Rv.resize(need); if (S.Rv2.size() < need) S.Rv2.resize(need);
  float* R0 = S.Rv.data(); float* R1 = S.Rv2.data(); alignas(16) float acc0[48], acc1[48];
  std::memset(R0 + SB, 0, (size_t)nsc2 * 4); std::memset(R1 + SB, 0, (size_t)nsc2 * 4);
  for (int j = 0; j < nS; j++) { const uint16_t w = T.sp[(size_t)j]; const uint32_t op = w & 63u, a = (w >> 6) & 31u; R0[SB + j] = (op == 1u) ? th0[a] : 1.0f; R1[SB + j] = (op == 1u) ? th1[a] : 1.0f; }
  R0[SB + nS] = 1.0f; R0[SB + nS + 1] = 2.0f; R1[SB + nS] = 1.0f; R1[SB + nS + 1] = 2.0f;
  std::memset(acc0, 0, sizeof acc0); std::memset(acc1, 0, sizeof acc1);
  for (int q = 0; q < T.Q; q++) {
    for (uint32_t i = sg[q]; i < sg[q + 1]; i++) {
      const uint32_t zw = tok[4 * i], zo = zw & 0xFFu, accn = (zw >> 8) & 1u, dst = tok[4 * i + 1], wa = tok[4 * i + 2], wb = tok[4 * i + 3];
      if (!accn && dst == 0xFFFFFFFFu) continue;
      const bool ca = wa >> 31, cb = wb >> 31; const uint32_t ofa = wa & 0x1FFFFFFFu, ofb = wb & 0x1FFFFFFFu;
      const bool va = ca || ((wa >> 30) & 1u) == 0u, vb = cb || ((wb >> 30) & 1u) == 0u; const bool ta = !ca && ((wa >> 29) & 1u), tb = !cb && ((wb >> 29) & 1u);
      const float* a0 = ca ? acc0 : (ta ? ts0 : R0) + ofa; const float* a1 = ca ? acc1 : (ta ? ts1 : R1) + ofa;
      const float* b0 = cb ? acc0 : (tb ? ts0 : R0) + ofb; const float* b1 = cb ? acc1 : (tb ? ts1 : R1) + ofb;
#ifdef COFIT_VREUSE
      float* sv0 = (vr0 && T.vro[i] >= 0) ? vr0 + T.vro[i] : nullptr; float* sv1 = sv0 ? vr1 + T.vro[i] : nullptr;
#ifdef COFIT_VR_CNT
      g_vrc.plaintok2[zo > 15 ? 15 : zo]++;
#endif
#endif
      const int csa = va ? 48 : nsc2, csb = vb ? 48 : nsc2; float* d0 = (dst != 0xFFFFFFFFu) ? R0 + dst : nullptr; float* d1 = (dst != 0xFFFFFFFFu) ? R1 + dst : nullptr;
      switch ((va ? 4 : 0) | (vb ? 2 : 0) | (accn ? 1 : 0)) {
#ifdef COFIT_H5_5
#define EX2(i, A_, B_, K_) case i: { const uint64_t cs2_ = (uint64_t)(uint32_t)csa | ((uint64_t)(uint32_t)csb << 32); exec_tok<A_, B_, K_, CF_P2DUAL>(zo, a0, b0, cs2_, acc0, d0, 0ull, nullptr); exec_tok<A_, B_, K_, CF_P2DUAL>(zo, a1, b1, cs2_, acc1, d1, 0ull, nullptr); } break;
#else
#ifdef COFIT_VREUSE
#define EX2(i, A_, B_, K_) case i: if (sv0) { exec_tok<A_, B_, K_, true, 1>(zo, 1, a0, csa, b0, csb, acc0, d0, 0u, 0u, nullptr, sv0); exec_tok<A_, B_, K_, true, 1>(zo, 1, a1, csa, b1, csb, acc1, d1, 0u, 0u, nullptr, sv1); } else { exec_tok<A_, B_, K_>(zo, 1, a0, csa, b0, csb, acc0, d0, 0u, 0u, nullptr, nullptr); exec_tok<A_, B_, K_>(zo, 1, a1, csa, b1, csb, acc1, d1, 0u, 0u, nullptr, nullptr); } break;
#else
#define EX2(i, A_, B_, K_) case i: exec_tok<A_, B_, K_>(zo, 1, a0, csa, b0, csb, acc0, d0, 0u, 0u, nullptr); exec_tok<A_, B_, K_>(zo, 1, a1, csa, b1, csb, acc1, d1, 0u, 0u, nullptr); break;
#endif
#endif
        EX2(0, false, false, false) EX2(1, false, false, true) EX2(2, false, true, false) EX2(3, false, true, true) EX2(4, true, false, false) EX2(5, true, false, true) EX2(6, true, true, false) EX2(7, true, true, true)
#undef EX2
      }
    }
    const uint32_t wo = sg[18 + q]; const bool co = wo >> 31; const uint32_t ofo = wo & 0x1FFFFFFFu; const bool vo = co || ((wo >> 30) & 1u) == 0u; const bool to_ = !co && ((wo >> 29) & 1u);
    const float* o0 = co ? acc0 : (to_ ? ts0 : R0) + ofo; const float* o1 = co ? acc1 : (to_ ? ts1 : R1) + ofo;
    float* p0 = out0 + (size_t)q * 48; float* p1 = out1 + (size_t)q * 48;
    if (vo) { std::memcpy(p0, o0, 48 * 4); std::memcpy(p1, o1, 48 * 4); }
    else { const float32x4_t x0 = vdupq_n_f32(o0[0]), x1 = vdupq_n_f32(o1[0]); for (int j = 0; j < 12; j++) { vst1q_f32(p0 + 4 * j, x0); vst1q_f32(p1 + 4 * j, x1); } }
  }
}
static inline void prep_warp(const fjn::WarpMat& W, Scratch& S) {
  if (S.stv && !std::memcmp(S.Skey, W.S, sizeof S.Skey)) return;
  std::memcpy(S.Skey, W.S, sizeof S.Skey); S.stv = true;
  for (int k = 0; k < 9; k++) for (int p = 0; p < 48; p++) S.ST[k * 48 + ppos(p)] = W.S[p * 9 + k];
}
static inline void tgrid(const Scratch& S, const float* u, float t0, float* t) {
  float u2[9]; float U = 0; for (int k = 0; k < 9; k++) { u2[k] = u[k] * u[k]; U += u2[k]; }
  const float32x4_t Uv = vdupq_n_f32(U), t0v = vdupq_n_f32(t0);
  for (int j = 0; j < 12; j++) { float32x4_t s = vdupq_n_f32(0.0f); for (int k = 0; k < 9; k++) s = vfmaq_f32(s, vld1q_f32(S.ST + k * 48 + 4 * j), vdupq_n_f32(u2[k])); vst1q_f32(t + 4 * j, vaddq_f32(t0v, vdivq_f32(s, Uv))); }
}
struct Mgs { PlainRes* r; int Q; bool bad; float on[16]; float nxt = 0.0f, rssv = 0.0f; bool hasRss = false; alignas(16) float t[16][16]; alignas(16) float tY[16]; };
static void plain_post(const Tree& T, const float* yv, PlainRes& r, Mgs& m);
static void plain_walk(const Tree& T, const float* th, const float* t48, const float* yv, Scratch& S, PlainRes& r, Mgs& m) {
  float* E = r.E; { PROFN(0); run_tree(T, false, th, t48, 1, S, E); }
#ifdef COFIT_PROF_REP
  { static thread_local int ctr = 0; if ((ctr++ & 7) == 0) { uint64_t t0 = mach_absolute_time(); for (int r_ = 0; r_ < 16; r_++) run_tree(T, false, th, t48, 1, S, E); g_prep_t += mach_absolute_time() - t0; g_prep_n += 16; g_prep_tok += 16 * (uint64_t)(T.tokP.size() / 4); } }
#endif
  plain_post(T, yv, r, m);
}
static void plain_post(const Tree& T, const float* yv, PlainRes& r, Mgs& m) {
  PHASE(2);
  float* E = r.E;
  const int Q = T.Q; m.r = &r; m.Q = Q;
  uint32x4_t fin = vdupq_n_u32(~0u); { const float32x4_t inf = vdupq_n_f32(INFINITY); for (int i = 0; i < Q * 12; i++) fin = vandq_u32(fin, COFIT_ABSLT(vld1q_f32(E + 4 * i), inf)); }
  m.bad = vminvq_u32(fin) == 0; m.hasRss = false;
  std::memcpy(r.Y, yv, 48 * 4);
  for (int k = 0; k < Q; k++) m.on[k] = std::sqrt(dot48(E + k * 48, E + k * 48));
  for (int k = 0; k < Q; k++) { std::memset(m.t[k], 0, sizeof m.t[k]); m.t[k][k] = 1.0f; } std::memset(m.tY, 0, sizeof m.tY);
}
#ifdef COFIT_X3
static inline float32x4_t dot4k_r(const float32x4_t* av, const float* b0, const float* b1, const float* b2, const float* b3) {
  const float* bb[4] = {b0, b1, b2, b3}; float32x2_t W[4];
  for (int d = 0; d < 4; d++) {
    float32x4_t L[4];
    for (int g = 0; g < 4; g++) L[g] = CF_DOT3(av[3 * g], vld1q_f32(bb[d] + 12 * g), av[3 * g + 1], vld1q_f32(bb[d] + 12 * g + 4), av[3 * g + 2], vld1q_f32(bb[d] + 12 * g + 8));
    const float32x4_t X = vpaddq_f32(L[0], L[1]), Yv = vpaddq_f32(L[2], L[3]), Z = vpaddq_f32(X, Yv);
    W[d] = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  }
  return vpaddq_f32(vcombine_f32(W[0], W[1]), vcombine_f32(W[2], W[3]));
}
#ifdef COFIT_DOTN
template <int NB>
static inline float32x4_t dot4k_rn(const float32x4_t* av, const float* b0, const float* b1, const float* b2, const float* b3) {
  const float* bb[4] = {b0, b1, b2, b3}; float32x2_t W[4];
#pragma clang loop unroll(full)
  for (int d = 0; d < NB; d++) {
    float32x4_t L[4];
    for (int g = 0; g < 4; g++) L[g] = CF_DOT3(av[3 * g], vld1q_f32(bb[d] + 12 * g), av[3 * g + 1], vld1q_f32(bb[d] + 12 * g + 4), av[3 * g + 2], vld1q_f32(bb[d] + 12 * g + 8));
    const float32x4_t X = vpaddq_f32(L[0], L[1]), Yv = vpaddq_f32(L[2], L[3]), Z = vpaddq_f32(X, Yv);
    W[d] = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  }
#pragma clang loop unroll(full)
  for (int d = NB; d < 4; d++) W[d] = W[0];
  return vpaddq_f32(vcombine_f32(W[0], W[1]), vcombine_f32(W[2], W[3]));
}
static inline float32x4_t dot4k_rs(int nb, const float32x4_t* av, const float* b0, const float* b1, const float* b2, const float* b3) {
  switch (nb) { case 1: return dot4k_rn<1>(av, b0, b1, b2, b3); case 2: return dot4k_rn<2>(av, b0, b1, b2, b3); case 3: return dot4k_rn<3>(av, b0, b1, b2, b3); default: return dot4k_rn<4>(av, b0, b1, b2, b3); }
}
#endif
static inline float dot48_r(const float32x4_t* a) {
  float32x4_t L[4];
  for (int g = 0; g < 4; g++) L[g] = CF_DOT3(a[3 * g], a[3 * g], a[3 * g + 1], a[3 * g + 1], a[3 * g + 2], a[3 * g + 2]);
  const float32x4_t X = vpaddq_f32(L[0], L[1]), Y = vpaddq_f32(L[2], L[3]), Z = vpaddq_f32(X, Y);
  const float32x2_t W = vadd_f32(vget_low_f32(Z), vget_high_f32(Z));
  return vget_lane_f32(W, 0) + vget_lane_f32(W, 1);
}
static inline void mgs_step(Mgs& m, int k, float tau) {
  float* E = m.r->E; float* Y = m.r->Y; const int Q = m.Q; float (*t)[16] = m.t; float* tY = m.tY;
  float* ck = E + k * 48;
  const float nrm = k == 0 ? m.on[0] : m.nxt; const bool keep = nrm > tau * m.on[k]; const float inv = keep ? 1.0f / nrm : 0.0f; const float32x4_t iv = vdupq_n_f32(inv);
#ifdef COFIT_MGST
  const int nvt = (Q + 3) >> 2;
#else
  const int nvt = 4;
#endif
  for (int v = 0; v < nvt; v++) vst1q_f32(t[k] + 4 * v, vmulq_f32(vld1q_f32(t[k] + 4 * v), iv));
  float32x4_t ckv[12]; for (int v = 0; v < 12; v++) { ckv[v] = vmulq_f32(vld1q_f32(ck + 4 * v), iv); vst1q_f32(ck + 4 * v, ckv[v]); }
  float32x4_t tk[4]; for (int v = 0; v < nvt; v++) tk[v] = vld1q_f32(t[k] + 4 * v);
  const int cnt = Q - k;
  for (int m0 = 0; m0 < cnt; m0 += 4) {
    const int nb = std::min(4, cnt - m0);
#ifdef COFIT_H5_4
#define H5CJ(q_) ({ const int j_ = k + 1 + m0 + ((q_) < nb ? (q_) : 0); (j_ < Q) ? E + j_ * 48 : Y; })
#define H5TJ(q_) ({ const int j_ = k + 1 + m0 + ((q_) < nb ? (q_) : 0); (j_ < Q) ? t[j_] : tY; })
#ifdef COFIT_DOTN
    const float32x4_t rr4 = dot4k_rs(nb, ckv, H5CJ(0), H5CJ(1), H5CJ(2), H5CJ(3));
#else
    const float32x4_t rr4 = dot4k_r(ckv, H5CJ(0), H5CJ(1), H5CJ(2), H5CJ(3));
#endif
    auto h5body = [&](auto QI) { constexpr int q = decltype(QI)::value; const float rr = vgetq_lane_f32(rr4, q); float* const cjq = H5CJ(q); float* const tjq = H5TJ(q);
#define H5Q_CJ cjq
#define H5Q_TJ tjq
#else
    float* cj[4]; float* tj[4];
    for (int q = 0; q < 4; q++) { const int j = k + 1 + m0 + (q < nb ? q : 0); cj[q] = (j < Q) ? E + j * 48 : Y; tj[q] = (j < Q) ? t[j] : tY; }
    const float32x4_t rr4 = dot4k_r(ckv, cj[0], cj[1], cj[2], cj[3]); float rrv[4]; vst1q_f32(rrv, rr4);
    for (int q = 0; q < nb; q++) { const float rr = rrv[q];
#define H5Q_CJ cj[q]
#define H5Q_TJ tj[q]
#endif
      const float32x4_t rv = vdupq_n_f32(rr), nrv = vdupq_n_f32(-rr);
      for (int v = 0; v < nvt; v++) vst1q_f32(H5Q_TJ + 4 * v, vfmaq_f32(vld1q_f32(H5Q_TJ + 4 * v), nrv, tk[v]));
      if (m0 == 0 && q == 0) {
        float32x4_t nv[12]; for (int v = 0; v < 12; v++) { nv[v] = CF_SM(vld1q_f32(H5Q_CJ + 4 * v), rv, ckv[v]); vst1q_f32(H5Q_CJ + 4 * v, nv[v]); }
        const float d2 = dot48_r(nv); if (k + 1 < Q) m.nxt = std::sqrt(d2); else { m.rssv = d2; m.hasRss = true; }
      } else
      for (int v = 0; v < 12; v++) vst1q_f32(H5Q_CJ + 4 * v, CF_SM(vld1q_f32(H5Q_CJ + 4 * v), rv, ckv[v])); };
#ifdef COFIT_H5_4
    h5body(std::integral_constant<int, 0>()); if (nb > 1) h5body(std::integral_constant<int, 1>()); if (nb > 2) h5body(std::integral_constant<int, 2>()); if (nb > 3) h5body(std::integral_constant<int, 3>());
#endif
  }
}
#undef H5Q_CJ
#undef H5Q_TJ
#else
static inline void mgs_step(Mgs& m, int k, float tau) {
  float* E = m.r->E; float* Y = m.r->Y; const int Q = m.Q; float (*t)[16] = m.t; float* tY = m.tY;
  float* ck = E + k * 48; const float nrm = std::sqrt(dot48(ck, ck)); const bool keep = nrm > tau * m.on[k]; const float inv = keep ? 1.0f / nrm : 0.0f; const float32x4_t iv = vdupq_n_f32(inv);
  for (int v = 0; v < 4; v++) vst1q_f32(t[k] + 4 * v, vmulq_f32(vld1q_f32(t[k] + 4 * v), iv));
  for (int v = 0; v < 12; v++) vst1q_f32(ck + 4 * v, vmulq_f32(vld1q_f32(ck + 4 * v), iv));
  float32x4_t tk[4]; for (int v = 0; v < 4; v++) tk[v] = vld1q_f32(t[k] + 4 * v);
  const int cnt = Q - k;
  for (int m0 = 0; m0 < cnt; m0 += 4) {
    const int nb = std::min(4, cnt - m0); float* cj[4]; float* tj[4];
    for (int q = 0; q < 4; q++) { const int j = k + 1 + m0 + (q < nb ? q : 0); cj[q] = (j < Q) ? E + j * 48 : Y; tj[q] = (j < Q) ? t[j] : tY; }
    const float32x4_t rr4 = dot4k(ck, cj[0], cj[1], cj[2], cj[3]); float rrv[4]; vst1q_f32(rrv, rr4);
    for (int q = 0; q < nb; q++) { const float rr = rrv[q]; const float32x4_t rv = vdupq_n_f32(rr), nrv = vdupq_n_f32(-rr);
      for (int v = 0; v < 4; v++) vst1q_f32(tj[q] + 4 * v, vfmaq_f32(vld1q_f32(tj[q] + 4 * v), nrv, tk[v]));
      for (int v = 0; v < 12; v++) vst1q_f32(cj[q] + 4 * v, CF_SM(vld1q_f32(cj[q] + 4 * v), rv, vld1q_f32(ck + 4 * v))); }
  }
}
#endif
static inline void mgs_finish(Mgs& m) {
#ifdef COFIT_X3
  PlainRes& r = *m.r; const float rss = m.hasRss ? m.rssv : dot48(r.Y, r.Y);
#else
  PlainRes& r = *m.r; const float rss = dot48(r.Y, r.Y);
#endif
  for (int q = 0; q < 16; q++) r.a[q] = (q < m.Q) ? -m.tY[q] : 0.0f;
  r.rss = (m.bad || !std::isfinite(rss)) ? INFINITY : rss;
}
struct Cand { alignas(16) float thx[32], thy[32], u[10]; alignas(16) float t[48]; PlainRes rx, ry; };
template <int KB>
static inline void ldl_trailing(int n, int j0, float* __restrict A, const float* __restrict FV, const float* const* ckp) {
  for (int j = j0; j < n; j++) { float* cj = A + j * NS; float32x4_t cv[KB]; for (int kk = 0; kk < KB; kk++) cv[kk] = vdupq_n_f32(ckp[kk][j]);
    for (int i = j; i < n; i += 4) { float32x4_t v = vld1q_f32(cj + i);
#pragma clang loop unroll(full)
      for (int kk = 0; kk < KB; kk++) v = CF_SM(v, vld1q_f32(FV + kk * (NS + 8) + i), cv[kk]);
      vst1q_f32(cj + i, v); } }
}
static inline bool solve_ldl(int n, float* __restrict A, float* __restrict x, float* __restrict FV, float* __restrict xs) {
  bool okc = true;
  for (int k0 = 0; k0 < n; k0 += 4) {
    const int kb = std::min(4, n - k0); const float* ckp[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int kk = 0; kk < kb; kk++) {
      const int k = k0 + kk; const float* ck = A + k * NS; ckp[kk] = ck; float* fv = FV + kk * (NS + 8); const float piv = ck[k]; okc = okc && (piv > 1e-12f) && std::isfinite(piv);
      const float id = 1.0f / std::max(piv, 1e-12f); const float cb = x[k] * id; const float32x4_t idv = vdupq_n_f32(id), cbv = vdupq_n_f32(cb);
      for (int i = k + 1; i < n; i += 4) { const float32x4_t a = vld1q_f32(ck + i); vst1q_f32(fv + i, vmulq_f32(a, idv)); vst1q_f32(x + i, CF_SM(vld1q_f32(x + i), a, cbv)); }
      for (int j = k + 1; j < k0 + kb; j++) { float* cj = A + j * NS; const float32x4_t ckj = vdupq_n_f32(ck[j]);
        for (int i = j; i < n; i += 4) vst1q_f32(cj + i, CF_SM(vld1q_f32(cj + i), vld1q_f32(fv + i), ckj)); }
    }
    switch (kb) { case 4: ldl_trailing<4>(n, k0 + kb, A, FV, ckp); break; case 3: ldl_trailing<3>(n, k0 + kb, A, FV, ckp); break; case 2: ldl_trailing<2>(n, k0 + kb, A, FV, ckp); break; default: ldl_trailing<1>(n, k0 + kb, A, FV, ckp); break; }
  }
  for (int i = n - 1; i >= 0; i--) { const float* ci = A + i * NS; float acc = x[i]; for (int k = n - 1; k > i; k--) CF_SSUB(acc, ci[k], xs[k]); xs[i] = acc / ci[i]; }
  for (int i = 0; i < n; i++) x[i] = xs[i];
  return okc;
}
#ifdef COFIT_RI_LDL
template <int NB>
static inline void ldl_step_ri(int n, int k, int b0, float* __restrict A, float* __restrict x, float id, float cb) {
  const float* ck = A + k * NS; float32x4_t fv[NB]; const float32x4_t idv = vdupq_n_f32(id), cbv = vdupq_n_f32(cb);
#pragma clang loop unroll(full)
  for (int t = 0; t < NB; t++) { const float32x4_t a = vld1q_f32(ck + 4 * (b0 + t)); fv[t] = vmulq_f32(a, idv); vst1q_f32(x + 4 * (b0 + t), vfmsq_f32(vld1q_f32(x + 4 * (b0 + t)), a, cbv)); }
  for (int j = k + 1; j < n; j++) { float* cj = A + j * NS + 4 * b0; const float32x4_t ckj = vld1q_dup_f32(ck + j);
#pragma clang loop unroll(full)
    for (int t = 0; t < NB; t++) vst1q_f32(cj + 4 * t, vfmsq_f32(vld1q_f32(cj + 4 * t), fv[t], ckj)); }
}
static inline bool solve_ldl_ri(int n, float* __restrict A, float* __restrict x, float* __restrict xs, float* __restrict ids) {
  bool okc = true; const int nbv = (n + 3) >> 2;
  for (int k = 0; k < n; k++) {
    const float piv = A[k * NS + k]; okc = okc && (piv > 1e-12f) && std::isfinite(piv);
    const float id = 1.0f / std::max(piv, 1e-12f); ids[k] = id; const float xk = x[k]; xs[k] = xk; const int b0 = (k + 1) >> 2;
    if (b0 >= nbv) break;
    switch (nbv - b0) { case 1: ldl_step_ri<1>(n, k, b0, A, x, id, xk * id); break; case 2: ldl_step_ri<2>(n, k, b0, A, x, id, xk * id); break; case 3: ldl_step_ri<3>(n, k, b0, A, x, id, xk * id); break; default: ldl_step_ri<4>(n, k, b0, A, x, id, xk * id); break; }
  }
#ifndef RI_NOBS
  for (int i = n - 1; i >= 0; i--) { const float* ci = A + i * NS; float acc = xs[i]; for (int k = n - 1; k > i; k--) acc = std::fmaf(-ci[k], xs[k], acc); xs[i] = acc * ids[i]; }
  for (int i = 0; i < n; i++) x[i] = xs[i];
#endif
  return okc;
}
#endif
template <int NB>
static inline void gj_step(int n, int k, float* __restrict A, float* __restrict x, float ip, float xk) {
  const float* rk = A + k * NS; float32x4_t rv[NB > 0 ? NB : 1];
#pragma clang loop unroll(full)
  for (int b = 0; b < NB; b++) rv[b] = vld1q_f32(rk + k + 1 + 4 * b);
#ifdef COFIT_GJ1
  auto row_ = [&](int i) { float* ri = A + i * NS; const float f = ri[k] * ip; const float32x4_t fv = vdupq_n_f32(f);
#pragma clang loop unroll(full)
    for (int b = 0; b < NB; b++) { const int j = k + 1 + 4 * b; vst1q_f32(ri + j, vfmsq_f32(vld1q_f32(ri + j), fv, rv[b])); }
    x[i] = std::fmaf(-f, xk, x[i]); };
  if (k + 1 < n) row_(k + 1);
  for (int i = 0; i < k; i++) row_(i);
  for (int i = k + 2; i < n; i++) row_(i);
#else
  for (int i = 0; i < n; i++) { if (i == k) continue; float* ri = A + i * NS; const float f = ri[k] * ip; const float32x4_t fv = vdupq_n_f32(f);
#pragma clang loop unroll(full)
    for (int b = 0; b < NB; b++) { const int j = k + 1 + 4 * b; vst1q_f32(ri + j, vfmsq_f32(vld1q_f32(ri + j), fv, rv[b])); }
    x[i] = std::fmaf(-f, xk, x[i]); }
#endif
}
static inline bool solve_gj(int n, float* __restrict A, float* __restrict x, float* __restrict dg) {
  bool okc = true;
  for (int k = 0; k < n; k++) {
    const float pk = A[k * NS + k]; okc = okc && (pk > 1e-12f) && std::isfinite(pk); const float ip = 1.0f / std::max(pk, 1e-12f); dg[k] = pk;
    const int nb = (n - k - 1 + 3) / 4; const float xk = x[k];
    switch (nb) { case 0: gj_step<0>(n, k, A, x, ip, xk); break; case 1: gj_step<1>(n, k, A, x, ip, xk); break; case 2: gj_step<2>(n, k, A, x, ip, xk); break; case 3: gj_step<3>(n, k, A, x, ip, xk); break; default: gj_step<4>(n, k, A, x, ip, xk); break; }
  }
  for (int i = 0; i < n; i++) x[i] = x[i] / std::max(dg[i], 1e-12f);
  return okc;
}
#ifdef COFIT_SV2
static inline bool solve_gj2(int n, float* __restrict A0, float* __restrict x0, float* __restrict d0, float* __restrict A1, float* __restrict x1, float* __restrict d1, bool& ok1) {
  bool ok0 = true; ok1 = true;
  for (int k = 0; k < n; k++) {
    const float p0 = A0[k * NS + k], p1 = A1[k * NS + k]; ok0 = ok0 && (p0 > 1e-12f) && std::isfinite(p0); ok1 = ok1 && (p1 > 1e-12f) && std::isfinite(p1);
    const float i0 = 1.0f / std::max(p0, 1e-12f), i1 = 1.0f / std::max(p1, 1e-12f); d0[k] = p0; d1[k] = p1;
    const int nb = (n - k - 1 + 3) / 4; const float xk0 = x0[k], xk1 = x1[k];
    switch (nb) { case 0: gj_step<0>(n, k, A0, x0, i0, xk0); gj_step<0>(n, k, A1, x1, i1, xk1); break; case 1: gj_step<1>(n, k, A0, x0, i0, xk0); gj_step<1>(n, k, A1, x1, i1, xk1); break;
      case 2: gj_step<2>(n, k, A0, x0, i0, xk0); gj_step<2>(n, k, A1, x1, i1, xk1); break; case 3: gj_step<3>(n, k, A0, x0, i0, xk0); gj_step<3>(n, k, A1, x1, i1, xk1); break;
      default: gj_step<4>(n, k, A0, x0, i0, xk0); gj_step<4>(n, k, A1, x1, i1, xk1); break; }
  }
  for (int i = 0; i < n; i++) { x0[i] = x0[i] / std::max(d0[i], 1e-12f); x1[i] = x1[i] / std::max(d1[i], 1e-12f); }
  return ok0;
}
#endif
static void fit_pair(const Ctx& cx, const Tree& X, const Tree& Y, const float* thx0, const float* thy0, int Qc, Scratch& S, Out& out) {
  PHASE(0);
  const fjn::WarpMat& W = *cx.WM; const int nx = X.nth, ny = Y.nth, nk = cx.closed ? fjn::KW : fjn::NW, n = nx + ny + nk;
  const float tau = (float)Qc * 9.0f * std::ldexp(1.0f, -24);
  prep_warp(W, S);
  alignas(16) float yvp[96]; for (int p = 0; p < 48; p++) { yvp[ppos(p)] = cx.yv96[p]; yvp[48 + ppos(p)] = cx.yv96[48 + p]; }
  const float* yv = yvp; const float flr = cx.flr; const int maxit = cx.maxit;
  Cand C[2];
  float thA[2][32] = {}; alignas(16) float uA[10]; alignas(16) float tA[48]; float Fcur = INFINITY, lamS = 1e-3f, Lc = 0; float rssA[2] = {INFINITY, INFINITY};
  PlainRes acc[2];
  for (int s_ = 0; s_ < 2; s_++) { acc[s_].rss = 0.0f; std::memset(acc[s_].a, 0, sizeof acc[s_].a); std::memset(acc[s_].Y, 0, sizeof acc[s_].Y); }
#ifdef COFIT_SWAPCOPY
  PlainRes* accp[2] = {&acc[0], &acc[1]}; PlainRes* rxp[2] = {&C[0].rx, &C[1].rx}; PlainRes* ryp[2] = {&C[0].ry, &C[1].ry};
#define COFIT_ACC(s_) (*accp[s_])
#define COFIT_CRX(h_) (*rxp[h_])
#define COFIT_CRY(h_) (*ryp[h_])
#else
#define COFIT_ACC(s_) acc[s_]
#define COFIT_CRX(h_) C[h_].rx
#define COFIT_CRY(h_) C[h_].ry
#endif
  for (int j = 0; j < nx; j++) thA[0][j] = thx0[j]; for (int j = 0; j < ny; j++) thA[1][j] = thy0[j];
  for (int k = 0; k < 9; k++) uA[k] = W.u0[k]; uA[9] = 0.0f;
  for (int h = 0; h < 2; h++) { std::memcpy(C[h].thx, thA[0], 128); std::memcpy(C[h].thy, thA[1], 128); std::memcpy(C[h].u, uA, 40); if (cx.closed) C[h].u[9] = (float)h * 0.125f; }
  for (int h = 0; h < 2; h++) tgrid(S, C[h].u, C[h].u[9], C[h].t);
  unsigned mode = 0, rungp = 0, scanr = 0, bestk = 0; bool done = false, init = true; float bestF = INFINITY; int it = 0;
  float plain = 0, dual = 0, trials = 0, solves = 0, stalled = 0;
  float pred[2]; bool fail[2]; float Delh[2][32 + 32 + 10];
  Out& o = out; o = Out(); std::memset(o.traj, 0, sizeof o.traj);
  const int nb4 = (n + 3) & ~3;
  static const bool lazy_ = getenv("COFIT_LAZY") && atoi(getenv("COFIT_LAZY")) != 0;
  bool pend1 = false;
#ifdef COFIT_VREUSE
#ifndef COFIT_FUSE2
#error COFIT_VREUSE needs COFIT_FUSE2 (the saves are made by run_tree2)
#endif
  bool vr_ok = false; int vr_ph = 0;
#ifdef COFIT_VR_CNT
  g_vrc.pairs++;
#endif
#endif
  auto solve_h_lam = [&](int h) {
        const float lamh = h ? Lc * std::exp2((float)(2u * rungp + 1u)) : Lc; bool okc = true;
        float* x = S.xv; for (int jb = 0; jb < nb4 + 4; jb += 4) vst1q_f32(x + jb, vmulq_f32(vld1q_f32(S.scl + jb), vld1q_f32(S.GV + jb)));
        float* A = S.A;
        COFIT_CPYROWS(A, S.Ab, n, nb4);
        { PROFN(6); PHASE(11);
#ifdef COFIT_RI_LDL
          if (n <= 16) { for (int i = 0; i < n; i++) A[i * NS + i] = S.Gd[i] + ((S.scl[i] > 0.0f) ? lamh : 1.0f); okc = solve_ldl_ri(n, A, x, S.xs, S.dg); }
#else
          if (n <= 16) { for (int i = 0; i < n; i++) A[i * NS + i] = S.Gd[i] + ((S.scl[i] > 0.0f) ? lamh : 1.0f); okc = solve_gj(n, A, x, S.dg); }
#endif
        else { for (int i = 0; i < n; i++) A[i * NS + i] += ((S.scl[i] > 0.0f) ? lamh : 1.0f); okc = solve_ldl(n, A, x, S.FV, S.xs); } }
        PHASE(12);
        for (int i = n; i < nb4 + 4; i++) x[i] = 0.0f;
        { float32x4_t pl[4] = {vdupq_n_f32(0.0f), vdupq_n_f32(0.0f), vdupq_n_f32(0.0f), vdupq_n_f32(0.0f)}; const float32x4_t lv = vdupq_n_f32(lamh), z0 = vdupq_n_f32(0.0f);
          for (int ib = 0; ib < nb4; ib += 4) { const float32x4_t xx = vld1q_f32(x + ib), ss = vld1q_f32(S.scl + ib), gg = vld1q_f32(S.GV + ib);
            const float32x4_t tl = vmulq_f32(vmulq_f32(lv, xx), xx); const float32x4_t term = CF_TERM(xx, ss, gg, vbslq_f32(vcgtq_f32(ss, z0), tl, z0));
            pl[(ib >> 2) & 3] = vaddq_f32(pl[(ib >> 2) & 3], term); }
          pred[h] = red16v(pl); }
        fail[h] = !okc;
        const bool f_ = !okc;
#ifdef COFIT_VSP
        for (int ib = 0; ib < nb4; ib += 4) vst1q_f32(Delh[h] + ib, vmulq_f32(vld1q_f32(S.scl + ib), vld1q_f32(x + ib)));
        { static const int32_t lane_[4] = {0, 1, 2, 3}; const int32x4_t lv = vld1q_s32(lane_); const float32x4_t z0 = vdupq_n_f32(0.0f);
          const int32x4_t nxv = vdupq_n_s32(f_ ? 0 : nx), nyv = vdupq_n_s32(f_ ? 0 : ny);
          for (int jb = 0; jb < 32; jb += 4) { const int32x4_t jv = vaddq_s32(lv, vdupq_n_s32(jb));
            const uint32x4_t mx = vcltq_s32(jv, nxv), my = vcltq_s32(jv, nyv);
            vst1q_f32(C[h].thx + jb, vaddq_f32(vld1q_f32(thA[0] + jb), vbslq_f32(mx, vld1q_f32(Delh[h] + jb), z0)));
            vst1q_f32(C[h].thy + jb, vaddq_f32(vld1q_f32(thA[1] + jb), vbslq_f32(my, vld1q_f32(Delh[h] + nx + jb), z0))); } }
#else
        for (int i = 0; i < n; i++) Delh[h][i] = S.scl[i] * x[i];
        for (int j = 0; j < 32; j++) { C[h].thx[j] = thA[0][j] + ((j < nx && !f_) ? Delh[h][j] : 0.0f); C[h].thy[j] = thA[1][j] + ((j < ny && !f_) ? Delh[h][nx + j] : 0.0f); }
#endif
        for (int k = 0; k < 10; k++) { float d_ = 0; if (k < nk && !f_) d_ = Delh[h][nx + ny + k]; float u_ = uA[k] + d_; if (k == 9) u_ = std::min(std::max(u_, 0.0f), 1.0f); C[h].u[k] = u_; }
        tgrid(S, C[h].u, C[h].u[9], C[h].t);
      };
#ifdef COFIT_X4
  auto& solve_h = solve_h_lam;
#else
  std::function<void(int)> solve_h;
#endif
  while (!done) {
    if (mode == 1) {  PROFN(3);
      for (int i = 0; i < n; i++) std::memset(S.Hd + i * NS, 0, (size_t)(nb4 + 4) * 4);
      std::memset(S.GV, 0, (size_t)(nb4 + 4) * 4); std::memset(S.DEL, 0, (size_t)(nb4 + 4) * 4);
      float UU = 0; for (int k = 0; k < 9; k++) CF_SACC(UU, uA[k], uA[k]);
      const float t0v = uA[9];
      alignas(16) float tdiff[48]; for (int j = 0; j < 12; j++) vst1q_f32(tdiff + 4 * j, vsubq_f32(vld1q_f32(tA + 4 * j), vdupq_n_f32(t0v)));
      for (int s = 0; s < 2; s++) {
        const Tree& T = s ? Y : X; const int nth = T.nth, ND = nth + 1, Cd = 1 + ND; const float* th = thA[s];
        alignas(16) float tsh[MAXC * 48]; std::memcpy(tsh, tA, 48 * 4); for (int c = 1; c < Cd; c++) { const float32x4_t z = vdupq_n_f32((c == Cd - 1) ? 1.0f : 0.0f); for (int j = 0; j < 12; j++) vst1q_f32(tsh + c * 48 + 4 * j, z); }
        float* outq = S.outq; bool coneOk = true; { PROFN(1);
#ifdef COFIT_VREUSE
          float* vrs_ = (vr_ok && T.nvr) ? S.VR[2 * vr_ph + s].data() : nullptr;
#ifdef COFIT_VR_CNT
          if (s == 0) { g_vrc.duals++; if (vr_ok) g_vrc.dualsreuse++; }
#endif
          if (!run_tree(T, true, th, tsh, Cd, S, outq, true, vrs_)) { coneOk = false; COFIT_CNT(1); run_tree(T, true, th, tsh, Cd, S, outq, false, vrs_); }
#else
          if (!run_tree(T, true, th, tsh, Cd, S, outq, true)) { coneOk = false; COFIT_CNT(1); run_tree(T, true, th, tsh, Cd, S, outq, false); }
#endif
 COFIT_CNT(0); }
        alignas(16) float V[MAXC * 48];
        PHASE(5);
        { PROFN(7); const float32x4_t inf = vdupq_n_f32(INFINITY);
          for (int d = 0; d < ND; d++) { float32x4_t v[12]; for (int j = 0; j < 12; j++) v[j] = vdupq_n_f32(0.0f);
#ifdef COFIT_X2
            for (int q = 0; q < T.Q; q++) { if (coneOk && !((T.mskOut[q] >> (d + 1)) & 1u)) continue; const float32x4_t aq = vdupq_n_f32(COFIT_ACC(s).a[q]); const float* src = outq + (size_t)q * Cd * 48 + (d + 1) * 48;
              for (int j = 0; j < 12; j++) v[j] = CF_AM(v[j], aq, vld1q_f32(src + 4 * j)); }
            { uint32x4_t fin = COFIT_ABSLT(v[0], inf); for (int j = 1; j < 12; j++) fin = vandq_u32(fin, COFIT_ABSLT(v[j], inf));
              if (vminvq_u32(fin) == 0) {
                for (int j = 0; j < 12; j++) v[j] = vdupq_n_f32(0.0f);
#endif
            for (int q = 0; q < T.Q; q++) { if (coneOk && !((T.mskOut[q] >> (d + 1)) & 1u)) continue; const float32x4_t aq = vdupq_n_f32(COFIT_ACC(s).a[q]); const float* src = outq + (size_t)q * Cd * 48 + (d + 1) * 48;
              for (int j = 0; j < 12; j++) { const float32x4_t t_ = vmulq_f32(aq, vld1q_f32(src + 4 * j)); const uint32x4_t fm = COFIT_ABSLT(t_, inf); v[j] = vaddq_f32(v[j], vreinterpretq_f32_u32(vandq_u32(fm, vreinterpretq_u32_f32(t_)))); } }
#ifdef COFIT_X2
              } }
#endif
            for (int j = 0; j < 12; j++) vst1q_f32(V + d * 48 + 4 * j, v[j]); } }
        PHASE(6);
        const int nsd = nth; const int mc = nsd + nk; float* COL = S.COL;
        for (int ic = 0; ic < mc; ic++) {
          float* cc = COL + ic * 48;
          if (ic < nsd) std::memcpy(cc, V + ic * 48, 48 * 4);
          else { const int k_ = ic - nsd; if (k_ < 9) { const float cf = 2.0f * uA[k_] / UU; const float32x4_t cfv = vdupq_n_f32(cf);
              for (int j = 0; j < 12; j++) vst1q_f32(cc + 4 * j, vmulq_f32(vld1q_f32(V + nsd * 48 + 4 * j), vmulq_f32(cfv, vsubq_f32(vld1q_f32(S.ST + k_ * 48 + 4 * j), vld1q_f32(tdiff + 4 * j))))); }
            else std::memcpy(cc, V + nsd * 48, 48 * 4); }
        }
        PHASE(7);
        for (int c0 = 0; c0 < mc; c0 += 4) {
          const int nb = std::min(4, mc - c0); float* cp[4]; for (int m = 0; m < 4; m++) cp[m] = COL + (c0 + (m < nb ? m : 0)) * 48;
          PHASE(7); float on[4]; for (int m = 0; m < nb; m++) on[m] = std::sqrt(dot48(cp[m], cp[m]));
          for (int m = 0; m < nb; m++) { const int ic = c0 + m; const int var = (ic < nsd) ? (s ? nx + ic : ic) : (nx + ny + ic - nsd); CF_SACC(S.DEL[var], on[m], on[m]); }
#ifdef COFIT_RI_CGS
          { alignas(16) float Cm[16][4];
            for (int k = 0; k < T.Q; k++) { vst1q_f32(Cm[k], dot4k(COFIT_ACC(s).E + k * 48, cp[0], cp[1], cp[2], cp[3])); }
#ifdef COFIT_RI_CGS2
            { int m = 0; for (; m + 1 < nb; m += 2) { float32x4_t v[12], w[12]; for (int j = 0; j < 12; j++) { v[j] = vld1q_f32(cp[m] + 4 * j); w[j] = vld1q_f32(cp[m + 1] + 4 * j); }
                for (int k = 0; k < T.Q; k++) { const float* qk = COFIT_ACC(s).E + k * 48; const float32x4_t rv = vdupq_n_f32(Cm[k][m]), rw = vdupq_n_f32(Cm[k][m + 1]);
                  for (int j = 0; j < 12; j++) { const float32x4_t q_ = vld1q_f32(qk + 4 * j); v[j] = CF_SM(v[j], rv, q_); w[j] = CF_SM(w[j], rw, q_); } }
                for (int j = 0; j < 12; j++) { vst1q_f32(cp[m] + 4 * j, v[j]); vst1q_f32(cp[m + 1] + 4 * j, w[j]); } }
              if (m < nb) { float32x4_t v[12]; for (int j = 0; j < 12; j++) v[j] = vld1q_f32(cp[m] + 4 * j);
                for (int k = 0; k < T.Q; k++) { const float* qk = COFIT_ACC(s).E + k * 48; const float32x4_t rv = vdupq_n_f32(Cm[k][m]);
                  for (int j = 0; j < 12; j++) v[j] = CF_SM(v[j], rv, vld1q_f32(qk + 4 * j)); }
                for (int j = 0; j < 12; j++) vst1q_f32(cp[m] + 4 * j, v[j]); } } }
#else
            for (int m = 0; m < nb; m++) { float32x4_t v[12]; for (int j = 0; j < 12; j++) v[j] = vld1q_f32(cp[m] + 4 * j);
              for (int k = 0; k < T.Q; k++) { const float* qk = COFIT_ACC(s).E + k * 48; const float32x4_t rv = vdupq_n_f32(Cm[k][m]);
                for (int j = 0; j < 12; j++) v[j] = CF_SM(v[j], rv, vld1q_f32(qk + 4 * j)); }
              for (int j = 0; j < 12; j++) vst1q_f32(cp[m] + 4 * j, v[j]); } }
#endif
#else
          PHASE(18);
#ifdef COFIT_TQK
          for (int k = 0; k < T.Q; k++) { const float* qk = COFIT_ACC(s).E + k * 48; float32x4_t qv[12];
#pragma clang loop unroll(full)
            for (int v = 0; v < 12; v++) qv[v] = vld1q_f32(qk + 4 * v);
#ifdef COFIT_DOTN
            const float32x4_t rr4 = dot4k_rs(nb, qv, cp[0], cp[1], cp[2], cp[3]); float rrv[4]; vst1q_f32(rrv, rr4);
#else
            const float32x4_t rr4 = dot4k_r(qv, cp[0], cp[1], cp[2], cp[3]); float rrv[4]; vst1q_f32(rrv, rr4);
#endif
            for (int m = 0; m < nb; m++) { const float32x4_t rv = vdupq_n_f32(rrv[m]);
#pragma clang loop unroll(full)
              for (int j = 0; j < 12; j++) vst1q_f32(cp[m] + 4 * j, CF_SM(vld1q_f32(cp[m] + 4 * j), rv, qv[j])); } }
#else
          for (int k = 0; k < T.Q; k++) { const float* qk = COFIT_ACC(s).E + k * 48; const float32x4_t rr4 = dot4k(qk, cp[0], cp[1], cp[2], cp[3]); float rrv[4]; vst1q_f32(rrv, rr4);
            for (int m = 0; m < nb; m++) { const float32x4_t rv = vdupq_n_f32(rrv[m]);
              for (int j = 0; j < 12; j++) vst1q_f32(cp[m] + 4 * j, CF_SM(vld1q_f32(cp[m] + 4 * j), rv, vld1q_f32(qk + 4 * j))); } }
#endif
#endif
          PHASE(19);
          for (int m = 0; m < nb; m++) { const float nrm = std::sqrt(dot48(cp[m], cp[m])); if (!(nrm > tau * on[m])) std::memset(cp[m], 0, 48 * 4); }
        }
#if defined(COFIT_X5) && !defined(COFIT_TILE_GRAM) && !defined(COFIT_EXACT_GRAM)
        { PROFN(5); PHASE(8);
          const float* Yg = COFIT_ACC(s).Y;
          for (int ic = 0; ic < mc; ic++) {
            const int i_ = (ic < nsd) ? (s ? nx + ic : ic) : (nx + ny + ic - nsd); const float* ca = COL + ic * 48;
            float32x4_t av[12];
#pragma clang loop unroll(full)
            for (int v = 0; v < 12; v++) av[v] = vld1q_f32(ca + 4 * v);
            const int ncol = ic + 2;
            for (int j0 = 0; j0 < ncol; j0 += 4) {
              const int nbc = std::min(4, ncol - j0); const float* bp[4]; for (int m = 0; m < 4; m++) { const int idx = j0 + (m < nbc ? m : 0); bp[m] = (idx <= ic) ? COL + idx * 48 : Yg; }
              float32x4_t Tm[4];
#pragma clang loop unroll(full)
              for (int m = 0; m < 4; m++) {
                float32x4_t s0 = vmulq_f32(av[0], vld1q_f32(bp[m])), s1 = vmulq_f32(av[1], vld1q_f32(bp[m] + 4)), s2 = vmulq_f32(av[2], vld1q_f32(bp[m] + 8));
#pragma clang loop unroll(full)
                for (int v = 3; v < 12; v += 3) { s0 = vfmaq_f32(s0, av[v], vld1q_f32(bp[m] + 4 * v)); s1 = vfmaq_f32(s1, av[v + 1], vld1q_f32(bp[m] + 4 * v + 4)); s2 = vfmaq_f32(s2, av[v + 2], vld1q_f32(bp[m] + 4 * v + 8)); }
                Tm[m] = vaddq_f32(vaddq_f32(s0, s1), s2); }
              float hv[4]; vst1q_f32(hv, vpaddq_f32(vpaddq_f32(Tm[0], Tm[1]), vpaddq_f32(Tm[2], Tm[3])));
              for (int m = 0; m < nbc; m++) { const int idx = j0 + m; const float hs = hv[m];
                if (idx <= ic) { const int j_ = (idx < nsd) ? (s ? nx + idx : idx) : (nx + ny + idx - nsd); S.Hd[i_ * NS + j_] += hs; if (i_ != j_) S.Hd[j_ * NS + i_] += hs; }
                else S.GV[i_] += hs; }
            }
          }
        }
#elif defined(COFIT_G22) && !defined(COFIT_TILE_GRAM) && !defined(COFIT_EXACT_GRAM)
        { PROFN(5); PHASE(8);
          int vi[28]; for (int ic = 0; ic < mc; ic++) vi[ic] = (ic < nsd) ? (s ? nx + ic : ic) : (nx + ny + ic - nsd);
          const float* Yg = COFIT_ACC(s).Y;
          for (int ic0 = 0; ic0 < mc; ic0 += 2) {
            const int ra = std::min(2, mc - ic0); const float* a0 = COL + ic0 * 48; const float* a1 = ra > 1 ? a0 + 48 : a0; float o[4];
            for (int jc0 = 0; jc0 < ic0 + ra; jc0 += 2) {
              const int cb = std::min(2, ic0 + ra - jc0); const float* b0 = COL + jc0 * 48; const float* b1 = cb > 1 ? b0 + 48 : b0;
              if (ra == 2) { if (cb == 2) gram_tile<2, 2>(a0, a1, b0, b1, o); else gram_tile<2, 1>(a0, a1, b0, b1, o); } else { if (cb == 2) gram_tile<1, 2>(a0, a1, b0, b1, o); else gram_tile<1, 1>(a0, a1, b0, b1, o); }
              for (int r = 0; r < ra; r++) for (int c = 0; c < cb; c++) { const int ic = ic0 + r, jc = jc0 + c; if (jc > ic) continue;
                const float hs = o[r * 2 + c]; const int i_ = vi[ic], j_ = vi[jc]; S.Hd[i_ * NS + j_] += hs; if (i_ != j_) S.Hd[j_ * NS + i_] += hs; }
            }
            if (ra == 2) gram_tile<2, 1>(a0, a1, Yg, Yg, o); else gram_tile<1, 1>(a0, a1, Yg, Yg, o);
            for (int r = 0; r < ra; r++) S.GV[vi[ic0 + r]] += o[r * 2];
          }
        }
#elif !defined(COFIT_TILE_GRAM) || defined(COFIT_EXACT_GRAM)
        { PROFN(5); PHASE(8);
        for (int ic = 0; ic < mc; ic++) {
          const int i_ = (ic < nsd) ? (s ? nx + ic : ic) : (nx + ny + ic - nsd); const float* ca = COL + ic * 48;
          for (int jc = 0; jc <= ic; jc++) {
            const int j_ = (jc < nsd) ? (s ? nx + jc : jc) : (nx + ny + jc - nsd);
            const float hs = dot_fast(ca, COL + jc * 48);
            S.Hd[i_ * NS + j_] += hs; if (i_ != j_) S.Hd[j_ * NS + i_] += hs; }
          S.GV[i_] += dot_fast(ca, COFIT_ACC(s).Y);
        }
        }
#else
        { PROFN(5);
          const int R = mc + 1; const float* rp[28]; for (int r = 0; r < R; r++) rp[r] = (r < mc) ? COL + r * 48 : COFIT_ACC(s).Y; for (int r = R; r < 28; r++) rp[r] = S.zero48;
          const int nbk = (R + 3) / 4; float* G = S.GR;
          for (int ib = 0; ib < nbk; ib++) for (int jb = 0; jb <= ib; jb++) {
            const float* a0p = rp[4 * ib]; const float* a1p = rp[4 * ib + 1]; const float* a2p = rp[4 * ib + 2]; const float* a3p = rp[4 * ib + 3];
            const float* b0p = rp[4 * jb]; const float* b1p = rp[4 * jb + 1]; const float* b2p = rp[4 * jb + 2]; const float* b3p = rp[4 * jb + 3];
            const float32x4_t z = vdupq_n_f32(0.0f);
            float32x4_t g00 = z, g01 = z, g02 = z, g03 = z, g10 = z, g11 = z, g12 = z, g13 = z, g20 = z, g21 = z, g22 = z, g23 = z, g30 = z, g31 = z, g32 = z, g33 = z;
#ifdef COFIT_H5_3
#pragma clang loop unroll(full)
#else
#pragma clang loop unroll(disable)
#endif
            for (int v = 0; v < 12; v++) {
              const float32x4_t A0 = vld1q_f32(a0p + 4 * v), A1 = vld1q_f32(a1p + 4 * v), A2 = vld1q_f32(a2p + 4 * v), A3 = vld1q_f32(a3p + 4 * v);
              const float32x4_t B0 = vld1q_f32(b0p + 4 * v), B1 = vld1q_f32(b1p + 4 * v), B2 = vld1q_f32(b2p + 4 * v), B3 = vld1q_f32(b3p + 4 * v);
#define FM(R_, C_) g##R_##C_ = vfmaq_f32(g##R_##C_, A##R_, B##C_);
              FM(0, 0) FM(0, 1) FM(0, 2) FM(0, 3) FM(1, 0) FM(1, 1) FM(1, 2) FM(1, 3) FM(2, 0) FM(2, 1) FM(2, 2) FM(2, 3) FM(3, 0) FM(3, 1) FM(3, 2) FM(3, 3)
#undef FM
            }
            vst1q_f32(G + (4 * ib + 0) * 28 + 4 * jb, vpaddq_f32(vpaddq_f32(g00, g01), vpaddq_f32(g02, g03)));
            vst1q_f32(G + (4 * ib + 1) * 28 + 4 * jb, vpaddq_f32(vpaddq_f32(g10, g11), vpaddq_f32(g12, g13)));
            vst1q_f32(G + (4 * ib + 2) * 28 + 4 * jb, vpaddq_f32(vpaddq_f32(g20, g21), vpaddq_f32(g22, g23)));
            vst1q_f32(G + (4 * ib + 3) * 28 + 4 * jb, vpaddq_f32(vpaddq_f32(g30, g31), vpaddq_f32(g32, g33)));
          }
          for (int ic = 0; ic < mc; ic++) {
            const int i_ = (ic < nsd) ? (s ? nx + ic : ic) : (nx + ny + ic - nsd);
            for (int jc = 0; jc <= ic; jc++) { const int j_ = (jc < nsd) ? (s ? nx + jc : jc) : (nx + ny + jc - nsd); const float hs = G[ic * 28 + jc];
              S.Hd[i_ * NS + j_] += hs; if (i_ != j_) S.Hd[j_ * NS + i_] += hs; }
            S.GV[i_] += G[mc * 28 + ic]; }
        }
#endif
      }
      PHASE(9);
      for (int i = 0; i < n; i++) { const float d_ = S.DEL[i]; S.scl[i] = (d_ > 0.0f && std::isfinite(d_)) ? 1.0f / std::sqrt(d_) : 0.0f; }
      for (int i = n; i < nb4 + 4; i++) { S.scl[i] = 0.0f; S.GV[i] = 0.0f; }
      for (int i = 0; i < n; i++) { const float32x4_t si = vdupq_n_f32(S.scl[i]); for (int jb = 0; jb < nb4; jb += 4) vst1q_f32(S.Ab + i * NS + jb, vmulq_f32(vmulq_f32(si, vld1q_f32(S.scl + jb)), vld1q_f32(S.Hd + i * NS + jb)));
        S.Gd[i] = (S.Hd[i * NS + i] * S.scl[i]) * S.scl[i]; }
      Lc = lamS; dual += 1; rungp = 0; mode = 2;
    }
    if (mode == 2) {  PROFN(4); PHASE(10);
      solves += 1;
#ifdef COFIT_SV2
      if (n <= 16) {
        float lam2[2]; bool ok2[2]; float* xs2[2] = {S.xv, S.xv2}; float* As2[2] = {S.A, S.A2}; float* dg2[2] = {S.dg, S.dg2};
        for (int h = 0; h < 2; h++) { lam2[h] = h ? Lc * std::exp2((float)(2u * rungp + 1u)) : Lc; float* x = xs2[h]; float* A = As2[h];
          for (int jb = 0; jb < nb4 + 4; jb += 4) vst1q_f32(x + jb, vmulq_f32(vld1q_f32(S.scl + jb), vld1q_f32(S.GV + jb)));
          COFIT_CPYROWS(A, S.Ab, n, nb4);
          for (int i = 0; i < n; i++) A[i * NS + i] = S.Gd[i] + ((S.scl[i] > 0.0f) ? lam2[h] : 1.0f); }
        { PROFN(6); ok2[0] = solve_gj2(n, As2[0], xs2[0], dg2[0], As2[1], xs2[1], dg2[1], ok2[1]); }
        for (int h = 0; h < 2; h++) { const float lamh = lam2[h]; float* x = xs2[h]; const bool okc = ok2[h];
          for (int i = n; i < nb4 + 4; i++) x[i] = 0.0f;
          { float32x4_t pl[4] = {vdupq_n_f32(0.0f), vdupq_n_f32(0.0f), vdupq_n_f32(0.0f), vdupq_n_f32(0.0f)}; const float32x4_t lv = vdupq_n_f32(lamh), z0 = vdupq_n_f32(0.0f);
            for (int ib = 0; ib < nb4; ib += 4) { const float32x4_t xx = vld1q_f32(x + ib), ss = vld1q_f32(S.scl + ib), gg = vld1q_f32(S.GV + ib);
              const float32x4_t tl = vmulq_f32(vmulq_f32(lv, xx), xx); const float32x4_t term = CF_TERM(xx, ss, gg, vbslq_f32(vcgtq_f32(ss, z0), tl, z0));
              pl[(ib >> 2) & 3] = vaddq_f32(pl[(ib >> 2) & 3], term); }
            pred[h] = red16v(pl); }
          fail[h] = !okc;
          for (int i = 0; i < n; i++) Delh[h][i] = S.scl[i] * x[i];
          const bool f_ = !okc;
          for (int j = 0; j < 32; j++) { C[h].thx[j] = thA[0][j] + ((j < nx && !f_) ? Delh[h][j] : 0.0f); C[h].thy[j] = thA[1][j] + ((j < ny && !f_) ? Delh[h][nx + j] : 0.0f); }
          for (int k = 0; k < 10; k++) { float d_ = 0; if (k < nk && !f_) d_ = Delh[h][nx + ny + k]; float u_ = uA[k] + d_; if (k == 9) u_ = std::min(std::max(u_, 0.0f), 1.0f); C[h].u[k] = u_; }
          tgrid(S, C[h].u, C[h].u[9], C[h].t); }
      } else
#endif
#ifndef COFIT_X4
      solve_h = solve_h_lam;
#endif
      solve_h(0); pend1 = lazy_ && !init; if (!pend1) solve_h(1);
      mode = 0;
    }
    if (mode == 0) {
      plain += 1;
      { PROFN(2); Mgs M[4];
        #ifdef COFIT_FUSE2
        const bool lz = false, skip0 = false;
#ifdef COFIT_VREUSE
        for (int s_ = 0; s_ < 2; s_++) { const Tree& T_ = s_ ? Y : X; for (int h_ = 0; h_ < 2; h_++) if (S.VR[2 * h_ + s_].size() < (size_t)T_.nvr + 4) S.VR[2 * h_ + s_].resize((size_t)T_.nvr + 4); }
        { PROFN(0); run_tree2(X, C[0].thx, C[1].thx, C[0].t, C[1].t, S, COFIT_CRX(0).E, COFIT_CRX(1).E, X.nvr ? S.VR[0].data() : nullptr, X.nvr ? S.VR[2].data() : nullptr); run_tree2(Y, C[0].thy, C[1].thy, C[0].t, C[1].t, S, COFIT_CRY(0).E, COFIT_CRY(1).E, Y.nvr ? S.VR[1].data() : nullptr, Y.nvr ? S.VR[3].data() : nullptr); }
        vr_ok = true;
#else
        { PROFN(0); run_tree2(X, C[0].thx, C[1].thx, C[0].t, C[1].t, S, COFIT_CRX(0).E, COFIT_CRX(1).E); run_tree2(Y, C[0].thy, C[1].thy, C[0].t, C[1].t, S, COFIT_CRY(0).E, COFIT_CRY(1).E); }
#endif
        plain_post(X, yv, COFIT_CRX(0), M[0]); plain_post(Y, yv + 48, COFIT_CRY(0), M[1]); plain_post(X, yv, COFIT_CRX(1), M[2]); plain_post(Y, yv + 48, COFIT_CRY(1), M[3]);
#else
        const bool lz = lazy_ && (!init || !cx.closed);
        const bool skip0 = lazy_ && !init && fail[0];
        if (!skip0) { plain_walk(X, C[0].thx, C[0].t, yv, S, COFIT_CRX(0), M[0]); plain_walk(Y, C[0].thy, C[0].t, yv + 48, S, COFIT_CRY(0), M[1]); }
        else { M[0].Q = 0; M[1].Q = 0; }
        if (!lz) { plain_walk(X, C[1].thx, C[1].t, yv, S, COFIT_CRX(1), M[2]); plain_walk(Y, C[1].thy, C[1].t, yv + 48, S, COFIT_CRY(1), M[3]); }
#endif
        PHASE(3);
        const int kmax = std::max(X.Q, Y.Q); const int nq = lz ? 2 : 4;
        for (int k = 0; k < kmax; k++) for (int q = 0; q < nq; q++) if (k < M[q].Q) mgs_step(M[q], k, tau);
        for (int q = 0; q < nq; q++) if (!(skip0 && q < 2)) mgs_finish(M[q]);
        if (skip0) { COFIT_CRX(0).rss = COFIT_CRY(0).rss = INFINITY; }
        if (lz && !init) { const float F0_ = COFIT_CRX(0).rss + COFIT_CRY(0).rss;
          if (fail[0] || !(F0_ < Fcur)) { if (pend1) { solve_h(1); pend1 = false; }
            if (fail[1]) { COFIT_CRX(1).rss = COFIT_CRY(1).rss = INFINITY; } else { plain_walk(X, C[1].thx, C[1].t, yv, S, COFIT_CRX(1), M[2]); plain_walk(Y, C[1].thy, C[1].t, yv + 48, S, COFIT_CRY(1), M[3]);
            for (int k = 0; k < kmax; k++) for (int q = 2; q < 4; q++) if (k < M[q].Q) mgs_step(M[q], k, tau);
            for (int q = 2; q < 4; q++) mgs_finish(M[q]); } }
          else COFIT_CRX(1).rss = COFIT_CRY(1).rss = INFINITY; }
        else if (lz) COFIT_CRX(1).rss = COFIT_CRY(1).rss = INFINITY; }
      PHASE(13);
      const float F0 = COFIT_CRX(0).rss + COFIT_CRY(0).rss, F1 = COFIT_CRX(1).rss + COFIT_CRY(1).rss; int pick = -1;
      if (init) {
        if (cx.closed && scanr < 4) {
          if (std::isfinite(F0) && F0 < bestF) { bestF = F0; bestk = 2 * scanr; }
          if (std::isfinite(F1) && F1 < bestF) { bestF = F1; bestk = 2 * scanr + 1; }
          scanr++;
          for (int h = 0; h < 2; h++) { C[h].u[9] = (scanr < 4) ? (float)(2 * scanr + h) * 0.125f : (float)bestk * 0.125f; tgrid(S, C[h].u, C[h].u[9], C[h].t); }
          if (scanr == 4 && !std::isfinite(bestF)) done = true;
          pick = -2;
        } else if (std::isfinite(F0)) pick = 0; else done = true;
      } else { if (!fail[0] && F0 < Fcur) pick = 0; else if (!fail[1] && F1 < Fcur) pick = 1; }
      if (done) break;
      if (pick >= 0) {
        const int ph = pick; const float Fn = ph == 0 ? F0 : F1;
#ifdef COFIT_VREUSE
        vr_ph = ph;
#endif
        std::memcpy(thA[0], C[ph].thx, 128); std::memcpy(thA[1], C[ph].thy, 128); std::memcpy(uA, C[ph].u, 40); std::memcpy(tA, C[ph].t, 48 * 4);
#ifdef COFIT_SWAPCOPY
        rssA[0] = rxp[ph]->rss; rssA[1] = ryp[ph]->rss; std::swap(accp[0], rxp[ph]); std::swap(accp[1], ryp[ph]);
#else
        copy_plain(COFIT_ACC(0), COFIT_CRX(ph), X.Q); copy_plain(COFIT_ACC(1), COFIT_CRY(ph), Y.Q); rssA[0] = COFIT_CRX(ph).rss; rssA[1] = COFIT_CRY(ph).rss;
#endif
        if (!init) {
          const float lamacc = (ph == 0) ? Lc : Lc * std::exp2((float)(2u * rungp + 1u)); const float pr = pred[ph];
          const float rho = std::min(std::max((Fcur - Fn) / std::max(pr, 1e-30f), -1e3f), 1e3f); const float c_ = 2.0f * rho - 1.0f;
          lamS = lamacc * std::max(1.0f / 3.0f, 1.0f - c_ * c_ * c_);
          if (it < 13) o.traj[it] = Fn;
        } else o.traj[13] = Fn;
        const float dF = Fcur - Fn; Fcur = Fn;
        if (init) { init = false; mode = 1; }
        else { it++; if (dF < flr || it >= maxit) done = true; else mode = 1; }
      } else if (pick == -2) {
      } else {
        Lc = Lc * std::exp2((float)(4u * rungp + 3u)); rungp++; trials += 2;
        if (rungp >= 6) { done = true; stalled = 1; } else mode = 2;
      }
    }
  }
  PHASE(14);
  o.Fx = rssA[0]; o.Fy = rssA[1]; o.stalled = stalled; o.it = (float)it; o.trials = trials; o.dual = dual; o.plain = plain; o.solves = solves;
  for (int j = 0; j < 32; j++) { o.th[0][j] = j < nx ? thA[0][j] : 0.0f; o.th[1][j] = j < ny ? thA[1][j] : 0.0f; }
  for (int j = 0; j < 16; j++) { o.a[0][j] = COFIT_ACC(0).a[j]; o.a[1][j] = COFIT_ACC(1).a[j]; } for (int j = 0; j < 10; j++) o.u[j] = uA[j];
  for (int p = 0; p < 48; p++) { o.res48[p] = COFIT_ACC(0).Y[ppos(p)]; o.res48[48 + p] = COFIT_ACC(1).Y[ppos(p)]; }
  o.ok = 1;
#ifdef COFIT_PROF2
  PHASE(15);
  g_pairs++; g_pn[std::min(n, 39)]++; g_pq[0][X.Q]++; g_pq[1][Y.Q]++; g_npl += (uint64_t)plain; g_ndu += (uint64_t)dual; g_nso += (uint64_t)solves; g_tr += (uint64_t)trials; if (n > 16) g_nsolve_ldl += (uint64_t)solves;
#endif
}
#ifdef COFIT_PROF2
static inline void prof_dump() {
  static const char* nm[24] = {"setup","plain_walk_tokens","plain_post","mgs","dual_walk_tokens","Vbuild","colbuild","tailproj","gram_grad","mode1_misc(zero,scl,Ab)","solve_prep(x,cpyrows,diag)","solve_core","solve_post(pred,Delh,tgrid)","accept_bookkeeping","end","after(outside fit_pair)","run_tree_setup","run_tree_outcopy","tail_projection","tail_recheck_norm","","","",""};
  fprintf(stderr, "COFITPROF pairs %llu plain %llu dual %llu solves %llu trials %llu ldl_solves %llu\n", (unsigned long long)g_pairs, (unsigned long long)g_npl, (unsigned long long)g_ndu, (unsigned long long)g_nso, (unsigned long long)g_tr, (unsigned long long)g_nsolve_ldl);
  for (int i = 0; i < 20; i++) fprintf(stderr, "COFITPROF phase %d %s ticks %llu cnt %llu\n", i, nm[i], (unsigned long long)g_pt[i], (unsigned long long)g_pcnt[i]);
  for (int d = 0; d < 2; d++) for (int z = 0; z < 16; z++) fprintf(stderr, "COFITPROF tok %s op %d n %llu ticks %llu\n", d ? "dual" : "plain", z, (unsigned long long)g_tokn[d][z], (unsigned long long)g_tokt[d][z]);
  for (int d = 0; d < 2; d++) for (int z = 7; z < 16; z++) fprintf(stderr, "COFITPROF fb %s op %d n %llu\n", d ? "dual" : "plain", z, (unsigned long long)g_fbt[d][z]);
  for (int i = 0; i < 40; i++) if (g_pn[i]) fprintf(stderr, "COFITPROF n %d pairs %llu\n", i, (unsigned long long)g_pn[i]);
  for (int sd = 0; sd < 2; sd++) for (int i = 0; i < 17; i++) if (g_pq[sd][i]) fprintf(stderr, "COFITPROF Q side %d Q %d pairs %llu\n", sd, i, (unsigned long long)g_pq[sd][i]);
  fprintf(stderr, "COFITPROF rep_ticks %llu walks %llu toks %llu\n", (unsigned long long)g_prep_t, (unsigned long long)g_prep_n, (unsigned long long)g_prep_tok);
  fprintf(stderr, "COFITPROF runone_ticks %llu scorer_ticks %llu\n", (unsigned long long)g_runone_t, (unsigned long long)g_scorer_t);
}
#endif
}
#pragma clang fp contract(on)
