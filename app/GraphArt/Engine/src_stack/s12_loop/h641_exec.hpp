#pragma once
#include "h421_exec.hpp"
inline bool g641_on() { static const bool v = getenv("H641") && atoi(getenv("H641")) != 0; return v; }
inline bool g642_on() { static const bool v = getenv("H642") && atoi(getenv("H642")) != 0; return v; }
inline int g643_bits() { static const int v = getenv("H643") ? atoi(getenv("H643")) : 0; return v; }

template <int NP, bool SA, bool SBB>
inline void H642_body(uint32_t op, const float* a, const float* as, const float* b, const float* bs, float* acc, int C) {
  if (C == 1 || !(op == 0u || op == 4u || op == 5u || (op >= 7u && op <= 15u))) { H421_body<NP, SA, SBB>(op, a, as, b, bs, acc, C); return; }
#define LA(c, p) (SA ? vdupq_n_f32(as[(c)]) : vld1q_f32(a + (c) * NP + (p)))
#define LB(c, p) (SBB ? vdupq_n_f32(bs[(c)]) : vld1q_f32(b + (c) * NP + (p)))
  constexpr int U = 4;
  static_assert(NP % (4 * U) == 0, "NP");
  if (op == 0u) {
    for (int p = 0; p < NP; p += 4 * U) {
      float32x4_t a0[U], b0[U];
      for (int u = 0; u < U; u++) { a0[u] = LA(0, p + 4 * u); b0[u] = LB(0, p + 4 * u); vst1q_f32(acc + p + 4 * u, vmulq_f32(a0[u], b0[u])); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < U; u++)
          vst1q_f32(acc + c * NP + p + 4 * u, vfmaq_f32(vmulq_f32(a0[u], LB(c, p + 4 * u)), LA(c, p + 4 * u), b0[u]));
    }
  } else if (op == 4u) {
    for (int p = 0; p < NP; p += 4 * U) {
      float32x4_t q0[U], b0[U];
      for (int u = 0; u < U; u++) { b0[u] = LB(0, p + 4 * u); q0[u] = vdivq_f32(LA(0, p + 4 * u), b0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < U; u++)
          vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LA(c, p + 4 * u), q0[u], LB(c, p + 4 * u)), b0[u]));
    }
  } else if (op == 5u) {
    for (int p = 0; p < NP; p += 4 * U) {
      float32x4_t q0[U], a0[U];
      for (int u = 0; u < U; u++) { a0[u] = LA(0, p + 4 * u); q0[u] = vdivq_f32(LB(0, p + 4 * u), a0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < U; u++)
          vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LB(c, p + 4 * u), q0[u], LA(c, p + 4 * u)), a0[u]));
    }
  } else {
    alignas(16) float m[NP], d[NP]; int n = NP;
    for (int p = 0; p < NP; p += 4) vst1q_f32(m + p, vmulq_f32(LA(0, p), LB(0, p)));
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
    for (int p = 0; p < NP; p += 4 * U) {
      float32x4_t a0[U], b0[U], dd[U];
      for (int u = 0; u < U; u++) { a0[u] = LA(0, p + 4 * u); b0[u] = LB(0, p + 4 * u); dd[u] = vld1q_f32(d + p + 4 * u); }
      for (int c = 1; c < C; c++)
        for (int u = 0; u < U; u++)
          vst1q_f32(acc + c * NP + p + 4 * u, vmulq_f32(dd[u], vfmaq_f32(vmulq_f32(a0[u], LB(c, p + 4 * u)), LA(c, p + 4 * u), b0[u])));
    }
  }
#undef LA
#undef LB
}
template <int NP, bool SA, bool SBB>
inline void H641X_body(uint32_t op, const float* a, const float* as, const float* b, const float* bs, float* acc, int C) {
  if (g642_on()) H642_body<NP, SA, SBB>(op, a, as, b, bs, acc, C); else H421_body<NP, SA, SBB>(op, a, as, b, bs, acc, C);
}

inline float H64_dot48r(const float32x4_t* x, const float32x4_t* y) {
  float32x4_t a0 = vdupq_n_f32(0.f), a1 = vdupq_n_f32(0.f), a2 = vdupq_n_f32(0.f);
  for (int u = 0; u < 12; u += 3) { a0 = vfmaq_f32(a0, x[u], y[u]); a1 = vfmaq_f32(a1, x[u + 1], y[u + 1]); a2 = vfmaq_f32(a2, x[u + 2], y[u + 2]); }
  return vaddvq_f32(vaddq_f32(vaddq_f32(a0, a1), a2));
}
inline float H64_dot48(const float* x, const float* y) {
  float32x4_t a0 = vdupq_n_f32(0.f), a1 = vdupq_n_f32(0.f), a2 = vdupq_n_f32(0.f);
  for (int k = 0; k < 48; k += 12) {
    a0 = vfmaq_f32(a0, vld1q_f32(x + k), vld1q_f32(y + k));
    a1 = vfmaq_f32(a1, vld1q_f32(x + k + 4), vld1q_f32(y + k + 4));
    a2 = vfmaq_f32(a2, vld1q_f32(x + k + 8), vld1q_f32(y + k + 8));
  }
  return vaddvq_f32(vaddq_f32(vaddq_f32(a0, a1), a2));
}
inline void H641_mgs(int n, float* const* col, float* tr, int L, float tau, bool dual, const float* on, float& pred) {
  pred = 0.0f;
  for (int k = 0; k < n - 1; k++) {
    float* ck = col[k];
    const float nrm = sqrtf(H64_dot48(ck, ck));
    const bool keep = nrm > tau * on[k];
    const float inv = keep ? 1.0f / nrm : 0.0f;
    for (int l = 0; l < L; l++) tr[k * L + l] *= inv;
    { const float32x4_t vs = vdupq_n_f32(inv); for (int p = 0; p < 48; p += 4) vst1q_f32(ck + p, vmulq_f32(vld1q_f32(ck + p), vs)); }
    float rj[48];
    for (int j = k + 1; j < n; j++) rj[j] = H64_dot48(ck, col[j]);
    for (int j = k + 1; j < n; j++) {
      const float r = rj[j];
      for (int l = 0; l < L; l++) tr[j * L + l] = fmaf(-r, tr[k * L + l], tr[j * L + l]);
      if (dual && j == n - 1) pred = fmaf(r, r, pred);
      const float32x4_t vr = vdupq_n_f32(r); float* cj = col[j];
      for (int p = 0; p < 48; p += 4) vst1q_f32(cj + p, vfmsq_f32(vld1q_f32(cj + p), vr, vld1q_f32(ck + p)));
    }
  }
}
inline void H643_mgs(int n, float* const* col, float* tr, int L, float tau, bool dual, const float* on, float& pred) {
  pred = 0.0f;
  if (n < 2) return;
  float nrm;
  { const float* c0 = col[0]; float32x4_t x[12]; for (int u = 0; u < 12; u++) x[u] = vld1q_f32(c0 + 4 * u); nrm = sqrtf(H64_dot48r(x, x)); }
  for (int k = 0; k < n - 1; k++) {
    float* ck = col[k];
    const bool keep = nrm > tau * on[k];
    const float inv = keep ? 1.0f / nrm : 0.0f;
    for (int l = 0; l < L; l++) tr[k * L + l] *= inv;
    float32x4_t q[12];
    { const float32x4_t vs = vdupq_n_f32(inv); for (int u = 0; u < 12; u++) { q[u] = vmulq_f32(vld1q_f32(ck + 4 * u), vs); vst1q_f32(ck + 4 * u, q[u]); } }
    for (int j = k + 1; j < n; j++) {
      float* cj = col[j];
      float32x4_t x[12];
      for (int u = 0; u < 12; u++) x[u] = vld1q_f32(cj + 4 * u);
      const float r = H64_dot48r(q, x);
      for (int l = 0; l < L; l++) tr[j * L + l] = fmaf(-r, tr[k * L + l], tr[j * L + l]);
      if (dual && j == n - 1) pred = fmaf(r, r, pred);
      const float32x4_t vr = vdupq_n_f32(r);
      for (int u = 0; u < 12; u++) { x[u] = vfmsq_f32(x[u], vr, q[u]); vst1q_f32(cj + 4 * u, x[u]); }
      if (j == k + 1 && k + 1 < n - 1) nrm = sqrtf(H64_dot48r(x, x));
    }
  }
}
inline void H64_mgs(int n, float* const* col, float* tr, int L, float tau, bool dual, const float* on, float& pred) {
  if (g643_bits() & 1) H643_mgs(n, col, tr, L, tau, dual, on, pred); else H641_mgs(n, col, tr, L, tau, dual, on, pred);
}
inline void H643_dual_setup(int Q, int nth, const float* roots, int rstride, const float* a, const float* E, float* V, float* on) {
  const float32x4_t vinf = vdupq_n_f32(INFINITY);
  for (int j = 0; j < nth; j++) {
    float32x4_t v[12];
    for (int u = 0; u < 12; u++) v[u] = vdupq_n_f32(0.0f);
    for (int q = 0; q < Q; q++) {
      const float32x4_t va = vdupq_n_f32(a[q]);
      const float* rq = roots + (size_t)q * rstride + (j + 1) * 48;
      for (int u = 0; u < 12; u++) {
        const float32x4_t t_ = vmulq_f32(va, vld1q_f32(rq + 4 * u));
        const float32x4_t tz = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(t_), vcaltq_f32(t_, vinf)));
        v[u] = vaddq_f32(v[u], tz);
      }
    }
    on[j] = sqrtf(H64_dot48r(v, v));
    for (int q = 0; q < Q; q++) {
      const float* qk = E + q * 48;
      float32x4_t e[12];
      for (int u = 0; u < 12; u++) e[u] = vld1q_f32(qk + 4 * u);
      const float r = H64_dot48r(e, v);
      const float32x4_t vr = vdupq_n_f32(r);
      for (int u = 0; u < 12; u++) v[u] = vfmsq_f32(v[u], vr, e[u]);
    }
    for (int u = 0; u < 12; u++) vst1q_f32(V + j * 48 + 4 * u, v[u]);
  }
}
