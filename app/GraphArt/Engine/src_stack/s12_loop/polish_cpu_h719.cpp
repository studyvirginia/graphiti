#include "polish_cpu_h719.hpp"
#include <Accelerate/Accelerate.h>
#include <arm_neon.h>
#include <cmath>
#include <cstring>
#include <algorithm>

#include <cstdlib>
#ifdef H840_PROF
#include <mach/mach_time.h>
static inline uint64_t h840_t() { uint64_t v; asm volatile("mrs %0, cntvct_el0" : "=r"(v)); return v; }
static uint64_t g840_acc[16], g840_cnt[16]; static uint64_t g840_cpc[4][2]; static uint64_t g840_op_t[2][16], g840_op_n[2][16], g840_res_t[2];
#define P840(i, t0) do { const uint64_t _t = h840_t(); g840_acc[i] += _t - (t0); g840_cnt[i]++; (t0) = _t; } while (0)
#define P840_START(t0) uint64_t t0 = h840_t()
#else
#define P840(i, t0) do {} while (0)
#define P840_START(t0) do {} while (0)
#endif
#include "h641_exec.hpp"
inline int g842p_bits() { static const int v = [] { const char* e = getenv("H842"); return e ? atoi(e) : 0; }(); return v; }
inline int g840_bits() { static const int v = [] { const char* e = getenv("H840"); if (!e) return 0; const int x = atoi(e); return x == 1 ? 255 : x; }(); return v; }
inline void g840_copy(float* __restrict d, const float* __restrict s, size_t nf) {
  for (size_t k = 0; k < nf; k += 8) { const float32x4_t a = vld1q_f32(s + k), b = vld1q_f32(s + k + 4); vst1q_f32(d + k, a); vst1q_f32(d + k + 4, b); }
}
#ifndef G644_ON_DEFINED
#define G644_ON_DEFINED
inline bool g644_on() { static const bool v = getenv("H644") && atoi(getenv("H644")) != 0; return v; }
#endif
inline void g840_pf1(const PolInputs& in, uint32_t tree) {
  __builtin_prefetch(in.hdr + 4 * (size_t)tree);
  const char* sg = (const char*)(in.segs + 36 * (size_t)tree); __builtin_prefetch(sg); __builtin_prefetch(sg + 64); __builtin_prefetch(sg + 143);
  const char* th = (const char*)(in.th0 + 32 * (size_t)tree); __builtin_prefetch(th); __builtin_prefetch(th + 64);
  const char* tv = (const char*)(in.tvec + 48 * (size_t)tree); __builtin_prefetch(tv); __builtin_prefetch(tv + 64); __builtin_prefetch(tv + 128);
}
inline void g840_pf2(const PolInputs& in, uint32_t tree, int Q) {
  const uint32_t* h = in.hdr + 4 * (size_t)tree;
  const uint32_t* sg = in.segs + 36 * (size_t)tree;
  const uint32_t x = h[0], nS = h[1] >> 16, z = h[2], ni = sg[Q];
  const char* t = (const char*)(in.tok + 4 * (size_t)x);
  for (uint32_t b = 0; b < 16 * ni; b += 64) __builtin_prefetch(t + b);
  if (ni) __builtin_prefetch(t + 16 * ni - 1);
  const char* sp = (const char*)(in.sprog + z);
  __builtin_prefetch(sp); if (nS > 32) __builtin_prefetch(sp + 2 * nS - 1);
}
namespace {
inline float p331_taud() { static const float v = getenv("H331_TAUD") ? (float)atof(getenv("H331_TAUD")) : 1.0f; return v; }
constexpr int P240_NP = 48;
constexpr int P240_MAXQ = 16;
constexpr int P240_MAXC = 25;
constexpr float P240_INF_ = INFINITY;
constexpr int P240_MAXN = 24, P240_MAXN1 = P240_MAXN + 1;

inline float P240_dot48(const float* x, const float* y) {
  float32x4_t a0 = vdupq_n_f32(0.f), a1 = vdupq_n_f32(0.f), a2 = vdupq_n_f32(0.f);
  for (int k = 0; k < P240_NP; k += 12) {
    a0 = vfmaq_f32(a0, vld1q_f32(x + k), vld1q_f32(y + k));
    a1 = vfmaq_f32(a1, vld1q_f32(x + k + 4), vld1q_f32(y + k + 4));
    a2 = vfmaq_f32(a2, vld1q_f32(x + k + 8), vld1q_f32(y + k + 8));
  }
  return vaddvq_f32(vaddq_f32(vaddq_f32(a0, a1), a2));
}
inline void P240_axpy48(float* y, float r, const float* x) {
  const float32x4_t vr = vdupq_n_f32(r);
  for (int k = 0; k < P240_NP; k += 4) vst1q_f32(y + k, vfmsq_f32(vld1q_f32(y + k), vr, vld1q_f32(x + k)));
}
inline void P240_scal48(float* y, float s) {
  const float32x4_t vs = vdupq_n_f32(s);
  for (int k = 0; k < P240_NP; k += 4) vst1q_f32(y + k, vmulq_f32(vld1q_f32(y + k), vs));
}
inline bool P240_fin(float x) { return std::isfinite(x); }
inline float P240_unv(int op, float m) {
  switch (op) {
    case 7: return sinf(m);
    case 8: return cosf(m);
    case 9: return expf(m);
    case 10: return logf(fabsf(m));
    case 11: return m * m;
    case 12: return m * m * m;
    case 13: return sqrtf(fabsf(m));
    case 14: return copysignf(powf(fabsf(m), 1.0f / 3.0f), m);
    default: return -m;
  }
}
inline float P240_und(int op, float m, float r) {
  switch (op) {
    case 7: return cosf(m);
    case 8: return -sinf(m);
    case 9: return r;
    case 10: return 1.0f / m;
    case 11: return 2.0f * m;
    case 12: return 3.0f * m * m;
    case 13: return copysignf(0.5f, m) / r;
    case 14: return 1.0f / (3.0f * r * r);
    default: return -1.0f;
  }
}
}

PolishWork::PolishWork() {}

const float* PolishWork::load(uint32_t w, float* buf, int C) {
  if (w >> 31) return accp_;
  const uint32_t off = w & 0x1FFFFFFFu;
  if ((w >> 30) & 1u) {
    for (int c = 0; c < C; c++) { const float32x4_t v = vdupq_n_f32(S_[(size_t)c * nsc2_ + (off - SB_)]); for (int p = 0; p < P240_NP; p += 4) vst1q_f32(buf + c * P240_NP + p, v); }
    return buf;
  }
  if ((w >> 29) & 1u) return TS_.data() + off;
  return R_.data() + vslot(off);
}

void PolishWork::exec(const uint32_t* w4, int C) {
  if (g421_on()) { exec421(w4, C); return; }
  if (g417_on()) { exec417(w4, C); return; }
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = load(w4[2], a_.data(), C);
  const float* b = load(w4[3], b_.data(), C);
  float* acc = (accp_ == acc_.data()) ? acc2_.data() : acc_.data();
  if (op == 0u) {
    for (int p = 0; p < P240_NP; p += 4) vst1q_f32(acc + p, vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < P240_NP; p += 4)
        vst1q_f32(acc + c * P240_NP + p, vfmaq_f32(vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + c * P240_NP + p)), vld1q_f32(a + c * P240_NP + p), vld1q_f32(b + p)));
  } else if (op < 4u) {
    for (int c = 0; c < C; c++)
      for (int p = 0; p < P240_NP; p += 4) {
        const float32x4_t va = vld1q_f32(a + c * P240_NP + p), vb = vld1q_f32(b + c * P240_NP + p);
        vst1q_f32(acc + c * P240_NP + p, op == 1u ? vaddq_f32(va, vb) : op == 2u ? vsubq_f32(va, vb) : vsubq_f32(vb, va));
      }
  } else if (op == 4u) {
    for (int p = 0; p < P240_NP; p += 4) vst1q_f32(acc + p, vdivq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < P240_NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), bc = vld1q_f32(b + c * P240_NP + p), ac = vld1q_f32(a + c * P240_NP + p);
        vst1q_f32(acc + c * P240_NP + p, vdivq_f32(vfmsq_f32(ac, q0, bc), vld1q_f32(b + p)));
      }
  } else if (op == 5u) {
    for (int p = 0; p < P240_NP; p += 4) vst1q_f32(acc + p, vdivq_f32(vld1q_f32(b + p), vld1q_f32(a + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < P240_NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), ac = vld1q_f32(a + c * P240_NP + p), bc = vld1q_f32(b + c * P240_NP + p);
        vst1q_f32(acc + c * P240_NP + p, vdivq_f32(vfmsq_f32(bc, q0, ac), vld1q_f32(a + p)));
      }
  } else if (op >= 7u && op <= 15u) {
    alignas(16) float m[P240_NP], d[P240_NP]; int n = P240_NP;
    for (int p = 0; p < P240_NP; p += 4) vst1q_f32(m + p, vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    switch (op) {
      case 7: vvsincosf(acc, d, m, &n); break;
      case 8: { vvsincosf(d, acc, m, &n); for (int p = 0; p < P240_NP; p += 4) vst1q_f32(d + p, vnegq_f32(vld1q_f32(d + p))); } break;
      case 9: vvexpf(acc, m, &n); std::memcpy(d, acc, sizeof d); break;
      case 10: { alignas(16) float am[P240_NP]; for (int p = 0; p < P240_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(am + p, vabsq_f32(vm)); vst1q_f32(d + p, vdivq_f32(vdupq_n_f32(1.0f), vm)); } vvlogf(acc, am, &n); } break;
      case 11: for (int p = 0; p < P240_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vm, vm)); vst1q_f32(d + p, vmulq_f32(vdupq_n_f32(2.0f), vm)); } break;
      case 12: for (int p = 0; p < P240_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vmulq_f32(vm, vm), vm)); vst1q_f32(d + p, vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vm), vm)); } break;
      case 13: { alignas(16) float am[P240_NP]; for (int p = 0; p < P240_NP; p += 4) vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vvsqrtf(acc, am, &n);
                 for (int p = 0; p < P240_NP; p++) d[p] = copysignf(0.5f, m[p]) / acc[p]; } break;
      case 14: { alignas(16) float am[P240_NP], e3[P240_NP]; for (int p = 0; p < P240_NP; p += 4) { vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vst1q_f32(e3 + p, vdupq_n_f32(1.0f / 3.0f)); }
                 vvpowf(acc, e3, am, &n); for (int p = 0; p < P240_NP; p++) { acc[p] = copysignf(acc[p], m[p]); d[p] = 1.0f / (3.0f * acc[p] * acc[p]); } } break;
      default: for (int p = 0; p < P240_NP; p += 4) { vst1q_f32(acc + p, vnegq_f32(vld1q_f32(m + p))); vst1q_f32(d + p, vdupq_n_f32(-1.0f)); } break;
    }
    for (int c = 1; c < C; c++)
      for (int p = 0; p < P240_NP; p += 4) {
        const float32x4_t a0 = vld1q_f32(a + p), b0 = vld1q_f32(b + p), ac = vld1q_f32(a + c * P240_NP + p), bc = vld1q_f32(b + c * P240_NP + p);
        { const float32x4_t sm = vfmaq_f32(vmulq_f32(a0, bc), ac, b0), r = vmulq_f32(vld1q_f32(d + p), sm);
          vst1q_f32(acc + c * P240_NP + p, g_h790 ? vbslq_f32(vceqq_f32(sm, vdupq_n_f32(0.f)), vdupq_n_f32(0.f), r) : r); }
      }
  } else {
    std::memcpy(acc, a, sizeof(float) * (size_t)C * P240_NP);
  }
  accp_ = acc;
  if (dst != 0xFFFFFFFFu) std::memcpy(R_.data() + vslot(dst), acc, sizeof(float) * (size_t)C * P240_NP);
}

void PolishWork::exec417(const uint32_t* w4, int C) {
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = nullptr; const float* b = nullptr; float as[P240_MAXC], bs[P240_MAXC]; bool sa = false, sb = false;
  auto res = [&](uint32_t w, const float*& p, float* s, bool& sc) {
    if (w >> 31) { p = accp_; return; }
    const uint32_t off = w & 0x1FFFFFFFu;
    if ((w >> 30) & 1u) { sc = true; for (int c = 0; c < C; c++) s[c] = S_[(size_t)c * nsc2_ + (off - SB_)]; return; }
    if ((w >> 29) & 1u) { p = TS_.data() + off; return; }
    p = R_.data() + vslot(off);
  };
  res(w4[2], a, as, sa); res(w4[3], b, bs, sb);
  float* acc = (accp_ == acc_.data()) ? acc2_.data() : acc_.data();
  float* R = dst != 0xFFFFFFFFu ? R_.data() + vslot(dst) : nullptr;
  if (R && R != a && R != b) acc = R;
  if (sa) { if (sb) H417_body<P240_NP, true, true>(op, a, as, b, bs, acc, C); else H417_body<P240_NP, true, false>(op, a, as, b, bs, acc, C); }
  else { if (sb) H417_body<P240_NP, false, true>(op, a, as, b, bs, acc, C); else H417_body<P240_NP, false, false>(op, a, as, b, bs, acc, C); }
  accp_ = acc;
  if (R && acc != R) std::memcpy(R, acc, sizeof(float) * (size_t)C * P240_NP);
}

template <bool SA, bool SBB>
__attribute__((noinline)) static void H840_unary(uint32_t op, const float* a, const float* asp, const float* b, const float* bsp, float* acc, int C, size_t ns) {
  float as[P240_MAXC], bs[P240_MAXC];
  if (SA) for (int c = 0; c < C; c++) as[c] = asp[(size_t)c * ns];
  if (SBB) for (int c = 0; c < C; c++) bs[c] = bsp[(size_t)c * ns];
  H641X_body<P240_NP, SA, SBB>(op, a, as, b, bs, acc, C);
}
template <bool SA, bool SBB>
static inline void H840_simple(uint32_t op, const float* a, const float* asp, const float* b, const float* bsp, float* acc, int C, size_t ns) {
  constexpr int NP = P240_NP;
#define LA(c, p) (SA ? vdupq_n_f32(asp[(size_t)(c) * ns]) : vld1q_f32(a + (c) * NP + (p)))
#define LB(c, p) (SBB ? vdupq_n_f32(bsp[(size_t)(c) * ns]) : vld1q_f32(b + (c) * NP + (p)))
  if (op == 0u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t a0[4], b0[4];
      for (int u = 0; u < 4; u++) { a0[u] = LA(0, p + 4 * u); b0[u] = LB(0, p + 4 * u); vst1q_f32(acc + p + 4 * u, vmulq_f32(a0[u], b0[u])); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vfmaq_f32(vmulq_f32(a0[u], LB(c, p + 4 * u)), LA(c, p + 4 * u), b0[u]));
    }
  } else if (op < 4u) {
    if (op == 1u) { for (int c = 0; c < C; c++) for (int p = 0; p < NP; p += 4) vst1q_f32(acc + c * NP + p, vaddq_f32(LA(c, p), LB(c, p))); }
    else if (op == 2u) { for (int c = 0; c < C; c++) for (int p = 0; p < NP; p += 4) vst1q_f32(acc + c * NP + p, vsubq_f32(LA(c, p), LB(c, p))); }
    else { for (int c = 0; c < C; c++) for (int p = 0; p < NP; p += 4) { const float32x4_t va = LA(c, p), vb = LB(c, p); vst1q_f32(acc + c * NP + p, vsubq_f32(vb, va)); } }
  } else if (op == 4u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t q0[4], b0[4];
      for (int u = 0; u < 4; u++) { b0[u] = LB(0, p + 4 * u); q0[u] = vdivq_f32(LA(0, p + 4 * u), b0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LA(c, p + 4 * u), q0[u], LB(c, p + 4 * u)), b0[u]));
    }
  } else if (op == 5u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t q0[4], a0[4];
      for (int u = 0; u < 4; u++) { a0[u] = LA(0, p + 4 * u); q0[u] = vdivq_f32(LB(0, p + 4 * u), a0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LB(c, p + 4 * u), q0[u], LA(c, p + 4 * u)), a0[u]));
    }
  } else {
    if (SA) { for (int c = 0; c < C; c++) for (int p = 0; p < NP; p += 4) vst1q_f32(acc + c * NP + p, LA(c, p)); }
    else g840_copy(acc, a, (size_t)C * NP);
  }
#undef LA
#undef LB
}

void PolishWork::exec421(const uint32_t* w4, int C) {
  P840_START(xt_);
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = nullptr; const float* b = nullptr; float as[P240_MAXC], bs[P240_MAXC]; bool sa = false, sb = false;
  auto res = [&](uint32_t w, const float*& p, float* s, bool& sc) {
    if (w >> 31) { p = accp_; return; }
    const uint32_t off = w & 0x1FFFFFFFu;
    if ((w >> 30) & 1u) { sc = true; for (int c = 0; c < C; c++) s[c] = S_[(size_t)c * nsc2_ + (off - SB_)]; return; }
    if ((w >> 29) & 1u) { p = TS_.data() + off; return; }
    p = R_.data() + vslot(off);
  };
  res(w4[2], a, as, sa); res(w4[3], b, bs, sb);
  float* acc = (accp_ == acc_.data()) ? acc2_.data() : acc_.data();
  float* R = dst != 0xFFFFFFFFu ? R_.data() + vslot(dst) : nullptr;
  if (R && R != a && R != b) acc = R;
#ifdef H840_PROF
  { const uint64_t _t = h840_t(); g840_res_t[C == 1 ? 0 : 1] += _t - xt_; xt_ = _t; }
#endif
  if (sa) { if (sb) H641X_body<P240_NP, true, true>(op, a, as, b, bs, acc, C); else H641X_body<P240_NP, true, false>(op, a, as, b, bs, acc, C); }
  else { if (sb) H641X_body<P240_NP, false, true>(op, a, as, b, bs, acc, C); else H641X_body<P240_NP, false, false>(op, a, as, b, bs, acc, C); }
#ifdef H840_PROF
  { const uint64_t _t = h840_t(); g840_op_t[C == 1 ? 0 : 1][op & 15] += _t - xt_; g840_op_n[C == 1 ? 0 : 1][op & 15]++; xt_ = _t; }
#endif
  accp_ = acc;
  if (R && acc != R) { if (g840_bits() & 1) g840_copy(R, acc, (size_t)C * P240_NP); else std::memcpy(R, acc, sizeof(float) * (size_t)C * P240_NP);
#ifdef H840_PROF
    g840_cpc[C == 1 ? 0 : 1][0]++; g840_cpc[C == 1 ? 0 : 1][1] += C;
#endif
  }
}

void PolishWork::exec840(const uint32_t* w4, int C) {
  const uint32_t op = w4[0], dst = w4[1];
  const float *a = nullptr, *b = nullptr, *asp = nullptr, *bsp = nullptr; bool sa = false, sb = false;
  auto res = [&](uint32_t w, const float*& p, const float*& sp, bool& sc) {
    if (w >> 31) { p = accp_; return; }
    const uint32_t off = w & 0x1FFFFFFFu;
    if ((w >> 30) & 1u) { sc = true; sp = S_.data() + (off - SB_); return; }
    if ((w >> 29) & 1u) { p = TS_.data() + off; return; }
    p = R_.data() + vslot(off);
  };
  res(w4[2], a, asp, sa); res(w4[3], b, bsp, sb);
  float* acc = (accp_ == acc_.data()) ? acc2_.data() : acc_.data();
  float* R = dst != 0xFFFFFFFFu ? R_.data() + vslot(dst) : nullptr;
  if (R && R != a && R != b) acc = R;
  const size_t ns = nsc2_;
  if (op >= 7u && op <= 15u) {
    if (sa) { if (sb) H840_unary<true, true>(op, a, asp, b, bsp, acc, C, ns); else H840_unary<true, false>(op, a, asp, b, bsp, acc, C, ns); }
    else { if (sb) H840_unary<false, true>(op, a, asp, b, bsp, acc, C, ns); else H840_unary<false, false>(op, a, asp, b, bsp, acc, C, ns); }
  } else {
    if (sa) { if (sb) H840_simple<true, true>(op, a, asp, b, bsp, acc, C, ns); else H840_simple<true, false>(op, a, asp, b, bsp, acc, C, ns); }
    else { if (sb) H840_simple<false, true>(op, a, asp, b, bsp, acc, C, ns); else H840_simple<false, false>(op, a, asp, b, bsp, acc, C, ns); }
  }
  accp_ = acc;
  if (R && acc != R) g840_copy(R, acc, (size_t)C * P240_NP);
}

void PolishWork::eval(const PolInputs& in, const Tree& T, const uint32_t* sg, const float* TH, const float* tvec, int C, int Q) {
  P840_START(pt_);
  const int pb_ = C == 1 ? 0 : 8; (void)pb_;
  cstride_ = (uint32_t)C;
  { const size_t need = ((size_t)per_ / (cd_ * 48u) + 2) * (size_t)C * 48u + 64; if (R_.size() < need) R_.resize(need); }
  const int b842 = g842p_bits();
  if (b842 & 1) {
    for (uint32_t j = 0; j < T.nS; ++j) {
      const uint32_t w = in.sprog[T.z + j];
      const uint32_t op = w & 63u, a = (w >> 6) & 31u;
      S_[j] = (op == 1u) ? TH[a] : 1.0f;
    }
    if (C > sb842_) {
      for (int c = sb842_; c < C; c++) {
        for (uint32_t j = 0; j < T.nS; ++j) {
          const uint32_t w = in.sprog[T.z + j];
          const uint32_t op = w & 63u, a = (w >> 6) & 31u;
          S_[(size_t)c * nsc2_ + j] = (op == 1u && a == (uint32_t)(c - 1)) ? 1.0f : 0.0f;
        }
        S_[(size_t)c * nsc2_ + T.nS] = 0.0f; S_[(size_t)c * nsc2_ + T.nS + 1] = 0.0f;
      }
      sb842_ = C;
    }
    S_[T.nS] = 1.0f; S_[T.nS + 1] = 2.0f;
    if (!ts0_842_) { std::memcpy(TS_.data(), tvec, sizeof(float) * P240_NP); ts0_842_ = true; }
    if (C > tsz842_) { std::memset(TS_.data() + (size_t)tsz842_ * P240_NP, 0, sizeof(float) * P240_NP * (size_t)(C - tsz842_)); tsz842_ = C; }
  } else {
  for (uint32_t j = 0; j < T.nS; ++j) {
    const uint32_t w = in.sprog[T.z + j];
    const uint32_t op = w & 63u, a = (w >> 6) & 31u;
    S_[j] = (op == 1u) ? TH[a] : 1.0f;
    for (int c = 1; c < C; c++) S_[(size_t)c * nsc2_ + j] = (op == 1u && a == (uint32_t)(c - 1)) ? 1.0f : 0.0f;
  }
  for (int c = 0; c < C; c++) {
    S_[(size_t)c * nsc2_ + T.nS] = c == 0 ? 1.0f : 0.0f;
    S_[(size_t)c * nsc2_ + T.nS + 1] = c == 0 ? 2.0f : 0.0f;
  }
  std::memcpy(TS_.data(), tvec, sizeof(float) * P240_NP);
  std::memset(TS_.data() + P240_NP, 0, sizeof(float) * P240_NP * (size_t)(C - 1));
  }
  bool zacc842 = true;
  if ((b842 & 2) && !g644_on() && sg[1] > sg[0]) { const uint32_t* w0 = in.tok + 4 * (size_t)(T.x + sg[0]); if (!((w0[2] | w0[3]) >> 31)) zacc842 = false; }
  if (zacc842) std::memset(acc_.data(), 0, sizeof(float) * P240_NP * (size_t)C);
#ifdef H842_CNT
  { static struct Cn { long n = 0, sk = 0; ~Cn() { fprintf(stderr, "H842_CNT %s evals %ld acc-zero skipped %ld\n", "p", n, sk); } } cn; cn.n++; if (!zacc842) cn.sk++; }
#endif
  accp_ = acc_.data();
  if (g644_on()) { eval644(in, T, sg, C, Q); return; }
  P840(pb_ + 0, pt_);
  const bool h840_x = (g840_bits() & 16) && g421_on();
  for (int q = 0; q < Q; q++) {
    for (uint32_t i = sg[q]; i < sg[q + 1]; i++) {
      if (h840_x) exec840(in.tok + 4 * (T.x + i), C);
      else exec(in.tok + 4 * (T.x + i), C);
    }
    P840(pb_ + 1, pt_);
    float* dq = rdst_ ? rdst_ + (size_t)q * rstr_ : roots_.data() + (size_t)q * P240_MAXC * P240_NP;
    const float* src = load(sg[18 + q], dq, C);
    if (src != dq) { if (g840_bits() & 1) g840_copy(dq, src, (size_t)C * P240_NP); else std::memcpy(dq, src, sizeof(float) * (size_t)C * P240_NP); }
#ifdef H840_PROF
    g840_cpc[2 + (C == 1 ? 0 : 1)][0]++; g840_cpc[2 + (C == 1 ? 0 : 1)][1] += C;
#endif
    P840(pb_ + 2, pt_);
  }
}

void PolishWork::eval644(const PolInputs& in, const Tree& T, const uint32_t* sg, int C, int Q) {
  K644& K = k644_[C == 1 ? 0 : 1];
  if (K.tok != in.tok || K.sg != sg || K.C != C || K.Q != Q || K.x != T.x || K.cd != cd_ || K.nsc2 != nsc2_ || K.SB != SB_ || K.R != R_.data() ||
      K.S != S_.data() || K.TS != TS_.data() || K.acc != acc_.data() || K.acc2 != acc2_.data() || K.h313 != false) {
    K.tok = in.tok; K.sg = sg; K.C = C; K.Q = Q; K.x = T.x; K.cd = cd_; K.nsc2 = nsc2_; K.SB = SB_; K.R = R_.data(); K.S = S_.data();
    K.TS = TS_.data(); K.acc = acc_.data(); K.acc2 = acc2_.data(); K.h313 = false;
    K.rec.clear(); K.beg.assign((size_t)Q + 1, 0); K.root.assign((size_t)Q, nullptr); K.rsc.assign((size_t)Q, 0);
    float* accp = acc_.data();
    auto res = [&](uint32_t w, const float*& p, const float*& s, uint8_t& sc) {
      p = nullptr; s = nullptr; sc = 0;
      if (w >> 31) { p = accp; return; }
      const uint32_t off = w & 0x1FFFFFFFu;
      if ((w >> 30) & 1u) { sc = 1; s = S_.data() + (off - SB_); return; }
      if ((w >> 29) & 1u) { p = TS_.data() + off; return; }
      p = R_.data() + vslot(off);
    };
    for (int q = 0; q < Q; q++) {
      K.beg[q] = (uint32_t)K.rec.size();
      for (uint32_t i = sg[q]; i < sg[q + 1]; i++) {
        const uint32_t* w4; uint32_t w4b[4]; (void)w4b;
        w4 = in.tok + 4 * (T.x + i);
        const uint32_t op = w4[0], dst = w4[1];
        const float *a, *as, *b, *bs; uint8_t sa, sb;
        res(w4[2], a, as, sa); res(w4[3], b, bs, sb);
        float* acc = (accp == acc_.data()) ? acc2_.data() : acc_.data();
        float* R = dst != 0xFFFFFFFFu ? R_.data() + vslot(dst) : nullptr;
        if (R && R != a && R != b) acc = R;
        K.rec.push_back(R644{op, sa, sb, sa ? as : a, sb ? bs : b, acc, (R && acc != R) ? R : nullptr});
        accp = acc;
      }
      const uint32_t w = sg[18 + q];
      const float *p, *s; uint8_t sc; res(w, p, s, sc);
      K.root[q] = sc ? s : p; K.rsc[q] = sc;
    }
    K.beg[Q] = (uint32_t)K.rec.size();
    K.accp_end = accp;
  }
  const size_t ns = nsc2_;
  for (int q = 0; q < Q; q++) {
    for (uint32_t i = K.beg[q]; i < K.beg[q + 1]; i++) {
      const R644& r = K.rec[i];
      float as[P240_MAXC], bs[P240_MAXC];
      if (r.sa) for (int c = 0; c < C; c++) as[c] = r.a[(size_t)c * ns];
      if (r.sb) for (int c = 0; c < C; c++) bs[c] = r.b[(size_t)c * ns];
      if (r.sa) { if (r.sb) H641X_body<P240_NP, true, true>(r.op, r.a, as, r.b, bs, r.acc, C); else H641X_body<P240_NP, true, false>(r.op, r.a, as, r.b, bs, r.acc, C); }
      else { if (r.sb) H641X_body<P240_NP, false, true>(r.op, r.a, as, r.b, bs, r.acc, C); else H641X_body<P240_NP, false, false>(r.op, r.a, as, r.b, bs, r.acc, C); }
      if (r.cp) std::memcpy(r.cp, r.acc, sizeof(float) * (size_t)C * P240_NP);
    }
    float* dq = roots_.data() + (size_t)q * P240_MAXC * P240_NP;
    if (K.rsc[q]) {
      for (int c = 0; c < C; c++) { const float32x4_t v = vdupq_n_f32(K.root[q][(size_t)c * ns]); for (int p = 0; p < P240_NP; p += 4) vst1q_f32(dq + c * P240_NP + p, v); }
    } else if (K.root[q] != dq) std::memcpy(dq, K.root[q], sizeof(float) * (size_t)C * P240_NP);
  }
  accp_ = K.accp_end;
}

namespace {
void P240_mgs(int n, float* cols, float* tr, int L, float tau, bool dual, const float* on, float& pred) {
  pred = 0.0f;
  for (int k = 0; k < n - 1; k++) {
    float* ck = cols + k * P240_NP;
    const float nrm = sqrtf(P240_dot48(ck, ck));
    const bool keep = nrm > tau * on[k];
    const float inv = keep ? 1.0f / nrm : 0.0f;
    for (int l = 0; l < L; l++) tr[k * L + l] *= inv;
    P240_scal48(ck, inv);
    float rj[P240_MAXN1 + P240_MAXQ + 1];
    for (int j = k + 1; j < n; j++) rj[j] = P240_dot48(ck, cols + j * P240_NP);
    for (int j = k + 1; j < n; j++) {
      const float r = rj[j];
      for (int l = 0; l < L; l++) tr[j * L + l] = fmaf(-r, tr[k * L + l], tr[j * L + l]);
      if (dual && j == n - 1) pred = fmaf(r, r, pred);
      P240_axpy48(cols + j * P240_NP, r, ck);
    }
  }
}
}

namespace {
struct P240_PassS { float E[P240_MAXQ * P240_NP]; float Y[P240_NP]; float a[P240_MAXQ]; float rss; };
}
void PolishWork::polish(const PolInputs& in, uint32_t tree, int Q, float tau, float flr, int itmax, PolResult& out) {
  const uint32_t* h = in.hdr + 4 * tree;
  Tree T{h[0], h[1] >> 16, h[2], h[3] & 255u, (h[3] >> 8) & 0xFFFFu, h[3] >> 24};
  const uint32_t* sg = in.segs + 36 * tree;
  nsc2_ = T.nsc2; SB_ = T.SB; per_ = in.per;
  cd_ = (T.nth <= 4) ? 1 + (T.nth + 1) / 2 : 3;
  if (S_.size() < (size_t)P240_MAXC * T.nsc2 + 8) S_.resize((size_t)P240_MAXC * T.nsc2 + 8);
  const size_t ts_sz842 = TS_.size();
  if (TS_.size() < (size_t)P240_MAXC * P240_NP) TS_.resize((size_t)P240_MAXC * P240_NP);
  if (acc_.size() < (size_t)P240_MAXC * P240_NP) acc_.resize((size_t)P240_MAXC * P240_NP);
  if (acc2_.size() < (size_t)P240_MAXC * P240_NP) acc2_.resize((size_t)P240_MAXC * P240_NP);
  if (roots_.size() < (size_t)P240_MAXQ * P240_MAXC * P240_NP) roots_.resize((size_t)P240_MAXQ * P240_MAXC * P240_NP);
  if (a_.size() < (size_t)P240_MAXC * P240_NP) a_.resize((size_t)P240_MAXC * P240_NP);
  if (b_.size() < (size_t)P240_MAXC * P240_NP) b_.resize((size_t)P240_MAXC * P240_NP);
  sb842_ = 1; ts0_842_ = false;
  if (TS_.size() != ts_sz842) tsz842_ = 1;
#ifdef H842_POISON
  { std::memset(S_.data(), 0xFF, sizeof(float) * S_.size()); std::memset(acc_.data(), 0xFF, sizeof(float) * acc_.size()); std::memset(acc2_.data(), 0xFF, sizeof(float) * acc2_.size()); std::memset(TS_.data(), 0xFF, sizeof(float) * P240_NP); }
#endif
  const float* tvec = in.tvec + 48 * tree;
  const int nth = (int)T.nth;
  const int Cd = nth + 1;

  auto plain = [&](const float* TH, P240_PassS& P) {
    const bool direct = (g840_bits() & 2) && !g644_on() && (g641_on() || (g643_bits() & 1));
    if (direct) { rdst_ = P.E; rstr_ = P240_NP; }
    eval(in, T, sg, TH, tvec, 1, Q);
    rdst_ = nullptr;
    P840_START(pp_);
    if (g641_on() || (g643_bits() & 1)) {
      float tr[(P240_MAXQ + 1) * P240_MAXQ], on[P240_MAXQ + 1]; float* col[P240_MAXQ + 1];
      bool bad = false;
      for (int q = 0; q < Q; q++) {
        const float* rq = direct ? (const float*)(P.E + q * P240_NP) : roots_.data() + (size_t)q * P240_MAXC * P240_NP;
        float* e = P.E + q * P240_NP; col[q] = e;
        uint32x4_t nf = vdupq_n_u32(0);
        for (int p = 0; p < P240_NP; p += 4) { const float32x4_t v = vld1q_f32(rq + p); if (!direct) vst1q_f32(e + p, v); nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(v, vdupq_n_f32(INFINITY)))); }
        if (vmaxvq_u32(nf)) bad = true;
        on[q] = sqrtf(H64_dot48(e, e));
      }
      std::memcpy(P.Y, in.yv, sizeof(float) * P240_NP); col[Q] = P.Y;
      for (int k = 0; k < (Q + 1) * Q; k++) tr[k] = 0.0f;
      for (int q = 0; q < Q; q++) tr[q * Q + q] = 1.0f;
      float pred;
      P840(3, pp_);
      H64_mgs(Q + 1, col, tr, Q, tau, false, on, pred);
      P840(4, pp_);
      for (int q = 0; q < Q; q++) P.a[q] = -tr[Q * Q + q];
      const float rss = H64_dot48(P.Y, P.Y);
      P.rss = (bad || !P240_fin(rss)) ? P240_INF_ : rss;
      return;
    }
    float cols[(P240_MAXQ + 1) * P240_NP], tr[(P240_MAXQ + 1) * P240_MAXQ], on[P240_MAXQ + 1];
    bool bad = false;
    for (int q = 0; q < Q; q++) {
      const float* rq = roots_.data() + (size_t)q * P240_MAXC * P240_NP;
      std::memcpy(cols + q * P240_NP, rq, sizeof(float) * P240_NP);
      uint32x4_t nf = vdupq_n_u32(0);
      for (int p = 0; p < P240_NP; p += 4) nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(vld1q_f32(rq + p), vdupq_n_f32(INFINITY))));
      if (vmaxvq_u32(nf)) bad = true;
    }
    std::memcpy(cols + Q * P240_NP, in.yv, sizeof(float) * P240_NP);
    for (int k = 0; k < Q; k++) on[k] = sqrtf(P240_dot48(cols + k * P240_NP, cols + k * P240_NP));
    for (int k = 0; k < (Q + 1) * Q; k++) tr[k] = 0.0f;
    for (int q = 0; q < Q; q++) tr[q * Q + q] = 1.0f;
    float pred;
    P240_mgs(Q + 1, cols, tr, Q, tau, false, on, pred);
    for (int q = 0; q < Q; q++) { std::memcpy(P.E + q * P240_NP, cols + q * P240_NP, sizeof(float) * P240_NP); P.a[q] = -tr[Q * Q + q]; }
    std::memcpy(P.Y, cols + Q * P240_NP, sizeof(float) * P240_NP);
    const float rss = P240_dot48(P.Y, P.Y);
    P.rss = (bad || !P240_fin(rss)) ? P240_INF_ : rss;
  };
  auto dual = [&](const float* TH, const P240_PassS& P, float* d, float& pred) {
    eval(in, T, sg, TH, tvec, Cd, Q);
    P840_START(pd_);
    if (g643_bits()) {
      float V[(P240_MAXN1) * P240_NP], tr[P240_MAXN1 * P240_MAXN], on[P240_MAXN1]; float* col[P240_MAXN1];
      if (g643_bits() & 2) H643_dual_setup(Q, nth, roots_.data(), P240_MAXC * P240_NP, P.a, P.E, V, on);
      else {
        for (int k = 0; k < (nth + 1) * P240_NP; k++) V[k] = 0.0f;
        for (int q = 0; q < Q; q++) {
          const float aq = P.a[q];
          const float* rq = roots_.data() + (size_t)q * P240_MAXC * P240_NP;
          const float32x4_t va = vdupq_n_f32(aq), vinf = vdupq_n_f32(INFINITY);
          for (int j = 0; j < nth; j++)
            for (int p = 0; p < P240_NP; p += 4) {
              const float32x4_t t_ = vmulq_f32(va, vld1q_f32(rq + (j + 1) * P240_NP + p));
              const float32x4_t tz = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(t_), vcaltq_f32(t_, vinf)));
              vst1q_f32(V + j * P240_NP + p, vaddq_f32(vld1q_f32(V + j * P240_NP + p), tz));
            }
        }
        for (int j = 0; j < nth; j++) on[j] = sqrtf(P240_dot48(V + j * P240_NP, V + j * P240_NP));
        for (int q = 0; q < Q; q++) {
          const float* qk = P.E + q * P240_NP;
          for (int j = 0; j < nth; j++) { const float r = P240_dot48(qk, V + j * P240_NP); P240_axpy48(V + j * P240_NP, r, qk); }
        }
      }
      std::memcpy(V + nth * P240_NP, P.Y, sizeof(float) * P240_NP);
      for (int k = 0; k < (nth + 1) * nth; k++) tr[k] = 0.0f;
      for (int j = 0; j < nth; j++) tr[j * nth + j] = 1.0f;
      for (int j = 0; j <= nth; j++) col[j] = V + j * P240_NP;
      P840(5, pd_);
      H64_mgs(nth + 1, col, tr, nth, tau * p331_taud(), true, on, pred);
      P840(6, pd_);
      const float rss = P240_dot48(V + nth * P240_NP, V + nth * P240_NP);
      for (int j = 0; j < nth; j++) d[j] = -tr[nth * nth + j];
      if (!P240_fin(rss)) pred = -1.0f;
      P840(7, pd_);
      return;
    }
    float V[(P240_MAXN1) * P240_NP], tr[P240_MAXN1 * P240_MAXN], on[P240_MAXN1];
    for (int k = 0; k < (nth + 1) * P240_NP; k++) V[k] = 0.0f;
    for (int q = 0; q < Q; q++) {
      const float aq = P.a[q];
      const float* rq = roots_.data() + (size_t)q * P240_MAXC * P240_NP;
      const float32x4_t va = vdupq_n_f32(aq), vinf = vdupq_n_f32(INFINITY);
      for (int j = 0; j < nth; j++)
        for (int p = 0; p < P240_NP; p += 4) {
          const float32x4_t t_ = vmulq_f32(va, vld1q_f32(rq + (j + 1) * P240_NP + p));
          const float32x4_t tz = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(t_), vcaltq_f32(t_, vinf)));
          vst1q_f32(V + j * P240_NP + p, vaddq_f32(vld1q_f32(V + j * P240_NP + p), tz));
        }
    }
    for (int j = 0; j < nth; j++) on[j] = sqrtf(P240_dot48(V + j * P240_NP, V + j * P240_NP));
    for (int q = 0; q < Q; q++) {
      const float* qk = P.E + q * P240_NP;
      for (int j = 0; j < nth; j++) { const float r = P240_dot48(qk, V + j * P240_NP); P240_axpy48(V + j * P240_NP, r, qk); }
    }
    std::memcpy(V + nth * P240_NP, P.Y, sizeof(float) * P240_NP);
    for (int k = 0; k < (nth + 1) * nth; k++) tr[k] = 0.0f;
    for (int j = 0; j < nth; j++) tr[j * nth + j] = 1.0f;
    P240_mgs(nth + 1, V, tr, nth, tau * p331_taud(), true, on, pred);
    const float rss = P240_dot48(V + nth * P240_NP, V + nth * P240_NP);
    for (int j = 0; j < nth; j++) d[j] = -tr[nth * nth + j];
    if (!P240_fin(rss)) pred = -1.0f;
  };

  float th[32] = {0}, cand0[32] = {0}, cand1[32] = {0}, d[32] = {0};
  for (int j = 0; j < nth; j++) th[j] = in.th0[32 * tree + j];
  P240_PassS PB[3];
  int ic = 0;
  P240_PassS* cur = &PB[0];
  uint32_t steps = 0, npp = 1, ndp = 0;
  plain(th, PB[0]);
  float rcur = PB[0].rss;
  const float rss0 = PB[0].rss;
  bool done = !P240_fin(rcur) || itmax <= 0;
  if (!done) {
    bool need_dual = true;
    float f = 1.0f;
    int hp = 0, it = 0;
    while (true) {
      if (need_dual) {
        ndp++;
        float pred;
        dual(th, *cur, d, pred);
        if (!(pred >= flr)) break;
        hp = 0; f = 1.0f;
      }
      for (int j = 0; j < nth; j++) { cand0[j] = th[j] + f * d[j]; cand1[j] = th[j] + 0.5f * f * d[j]; }
      npp++;
      const int i0 = (ic + 1) % 3, i1 = (ic + 2) % 3;
      plain(cand0, PB[i0]);
      const float r0 = PB[i0].rss;
      float r1 = r0;
      int pk = -1;
      if (r0 <= rcur) pk = 0;
      else if (hp < 2) { plain(cand1, PB[i1]); r1 = PB[i1].rss; if (r1 <= rcur) pk = 1; }
      if (pk >= 0) {
        rcur = pk ? r1 : r0;
        for (int j = 0; j < nth; j++) th[j] = pk ? cand1[j] : cand0[j];
        ic = pk ? i1 : i0;
        cur = &PB[ic];
        steps++; it++;
        if (it >= itmax || !P240_fin(rcur)) break;
        need_dual = true;
        continue;
      }
      hp++; f *= 0.25f;
      if (hp >= 3) break;
      need_dual = false;
    }
  }
  std::memset(out.th, 0, sizeof out.th);
  for (int j = 0; j < nth; j++) out.th[j] = th[j];
  std::memset(out.a, 0, sizeof out.a);
  for (int q = 0; q < Q; q++) out.a[q] = cur->a[q];
  out.rss0 = rss0; out.rcur = rcur;
  out.steps = steps; out.npp = npp; out.ndp = ndp;
}

void PolishWork::roots_debug(const PolInputs& in, uint32_t tree, const float* TH, int C, int Q, float* out) {
  const uint32_t* h = in.hdr + 4 * tree;
  Tree T{h[0], h[1] >> 16, h[2], h[3] & 255u, (h[3] >> 8) & 0xFFFFu, h[3] >> 24};
  nsc2_ = T.nsc2; SB_ = T.SB;
  S_.assign((size_t)P240_MAXC * T.nsc2 + 8, 0.0f);
  per_ = in.per; cd_ = (T.nth <= 4) ? 1 + (T.nth + 1) / 2 : 3; R_.assign(((size_t)in.per / (cd_ * 48u) + 2) * (size_t)(P240_MAXC) * 48u + 64, 0.0f);
  TS_.assign((size_t)P240_MAXC * P240_NP, 0.0f); acc_.assign((size_t)P240_MAXC * P240_NP, 0.0f); acc2_.assign((size_t)P240_MAXC * P240_NP, 0.0f);
  roots_.assign((size_t)P240_MAXQ * P240_MAXC * P240_NP, 0.0f); a_.assign((size_t)P240_MAXC * P240_NP, 0.0f); b_.assign((size_t)P240_MAXC * P240_NP, 0.0f);
  sb842_ = 1; tsz842_ = 1; ts0_842_ = false;
  eval(in, T, in.segs + 36 * tree, TH, in.tvec + 48 * tree, C, Q);
  for (int q = 0; q < Q; q++) std::memcpy(out + (size_t)q * C * P240_NP, roots_.data() + (size_t)q * P240_MAXC * P240_NP, sizeof(float) * C * P240_NP);
}
