#pragma once
#include "h417_exec.hpp"
inline bool g421_on() { static const bool v = getenv("H421") && atoi(getenv("H421")) != 0; return v; }
template <int NP, bool SA, bool SBB>
inline void H421_body(uint32_t op, const float* a, const float* as, const float* b, const float* bs, float* acc, int C) {
#define LA(c, p) (SA ? vdupq_n_f32(as[(c)]) : vld1q_f32(a + (c) * NP + (p)))
#define LB(c, p) (SBB ? vdupq_n_f32(bs[(c)]) : vld1q_f32(b + (c) * NP + (p)))
  if (op == 0u) {
    for (int p = 0; p < NP; p += 4) vst1q_f32(acc + p, vmulq_f32(LA(0, p), LB(0, p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < NP; p += 4)
        vst1q_f32(acc + c * NP + p, vfmaq_f32(vmulq_f32(LA(0, p), LB(c, p)), LA(c, p), LB(0, p)));
  } else if (op < 4u) {
    for (int c = 0; c < C; c++)
      for (int p = 0; p < NP; p += 4) {
        const float32x4_t va = LA(c, p), vb = LB(c, p);
        vst1q_f32(acc + c * NP + p, op == 1u ? vaddq_f32(va, vb) : op == 2u ? vsubq_f32(va, vb) : vsubq_f32(vb, va));
      }
  } else if (op == 4u) {
    for (int p = 0; p < NP; p += 4) vst1q_f32(acc + p, vdivq_f32(LA(0, p), LB(0, p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), bc = LB(c, p), ac = LA(c, p);
        vst1q_f32(acc + c * NP + p, vdivq_f32(vfmsq_f32(ac, q0, bc), LB(0, p)));
      }
  } else if (op == 5u) {
    for (int p = 0; p < NP; p += 4) vst1q_f32(acc + p, vdivq_f32(LB(0, p), LA(0, p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), ac = LA(c, p), bc = LB(c, p);
        vst1q_f32(acc + c * NP + p, vdivq_f32(vfmsq_f32(bc, q0, ac), LA(0, p)));
      }
  } else if (op >= 7u && op <= 15u) {
    alignas(16) float m[NP], d[NP]; int n = NP;
    for (int p = 0; p < NP; p += 4) vst1q_f32(m + p, vmulq_f32(LA(0, p), LB(0, p)));
    if (C == 1) {
      switch (op) {
        case 7: vvsinf(acc, m, &n); break;
        case 8: vvsincosf(d, acc, m, &n); break;
        case 9: vvexpf(acc, m, &n); break;
        case 10: { alignas(16) float am[NP]; for (int p = 0; p < NP; p += 4) vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vvlogf(acc, am, &n); } break;
        case 11: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vm, vm)); } break;
        case 12: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vmulq_f32(vm, vm), vm)); } break;
        case 13: { alignas(16) float am[NP]; for (int p = 0; p < NP; p += 4) vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vvsqrtf(acc, am, &n); } break;
        case 14: { alignas(16) float am[NP], e3[NP]; for (int p = 0; p < NP; p += 4) { vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vst1q_f32(e3 + p, vdupq_n_f32(1.0f / 3.0f)); }
                   vvpowf(acc, e3, am, &n); for (int p = 0; p < NP; p++) acc[p] = copysignf(acc[p], m[p]); } break;
        default: for (int p = 0; p < NP; p += 4) vst1q_f32(acc + p, vnegq_f32(vld1q_f32(m + p))); break;
      }
      return;
    }
    switch (op) {
      case 7: vvsincosf(acc, d, m, &n); break;
      case 8: { vvsincosf(d, acc, m, &n); for (int p = 0; p < NP; p += 4) vst1q_f32(d + p, vnegq_f32(vld1q_f32(d + p))); } break;
      case 9: vvexpf(acc, m, &n); std::memcpy(d, acc, sizeof d); break;
      case 10: { alignas(16) float am[NP]; for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(am + p, vabsq_f32(vm)); vst1q_f32(d + p, vdivq_f32(vdupq_n_f32(1.0f), vm)); } vvlogf(acc, am, &n); } break;
      case 11: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vm, vm)); vst1q_f32(d + p, vmulq_f32(vdupq_n_f32(2.0f), vm)); } break;
      case 12: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vmulq_f32(vm, vm), vm)); vst1q_f32(d + p, vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vm), vm)); } break;
      case 13: { alignas(16) float am[NP]; for (int p = 0; p < NP; p += 4) vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vvsqrtf(acc, am, &n);
                 for (int p = 0; p < NP; p++) d[p] = copysignf(0.5f, m[p]) / acc[p]; } break;
      case 14: { alignas(16) float am[NP], e3[NP]; for (int p = 0; p < NP; p += 4) { vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vst1q_f32(e3 + p, vdupq_n_f32(1.0f / 3.0f)); }
                 vvpowf(acc, e3, am, &n); for (int p = 0; p < NP; p++) { acc[p] = copysignf(acc[p], m[p]); d[p] = 1.0f / (3.0f * acc[p] * acc[p]); } } break;
      default: for (int p = 0; p < NP; p += 4) { vst1q_f32(acc + p, vnegq_f32(vld1q_f32(m + p))); vst1q_f32(d + p, vdupq_n_f32(-1.0f)); } break;
    }
    for (int c = 1; c < C; c++)
      for (int p = 0; p < NP; p += 4) {
        const float32x4_t a0 = LA(0, p), b0 = LB(0, p), ac = LA(c, p), bc = LB(c, p);
        vst1q_f32(acc + c * NP + p, vmulq_f32(vld1q_f32(d + p), vfmaq_f32(vmulq_f32(a0, bc), ac, b0)));
      }
  } else {
    if (SA) { for (int c = 0; c < C; c++) for (int p = 0; p < NP; p += 4) vst1q_f32(acc + c * NP + p, LA(c, p)); }
    else std::memcpy(acc, a, sizeof(float) * (size_t)C * NP);
  }
#undef LA
#undef LB
}
