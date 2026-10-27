#include "gn_cpu_h940.hpp"
#include <Accelerate/Accelerate.h>
#include <arm_neon.h>
#include <cmath>
#include <cstring>
#include <algorithm>

#include <cstdlib>
#include <map>
#include <cstdio>
#ifdef H836_PROF
#include <mach/mach_time.h>
static inline uint64_t h836_t() { uint64_t v; asm volatile("mrs %0, cntvct_el0" : "=r"(v)); return v; }
static uint64_t g836_acc[16], g836_cnt[16]; static uint64_t g836_cpc[4][2]; static uint64_t g836_op_t[2][16], g836_op_n[2][16], g836_res_t[2];
#define P836(i, t0) do { const uint64_t _t = h836_t(); g836_acc[i] += _t - (t0); g836_cnt[i]++; (t0) = _t; } while (0)
#define P836_START(t0) uint64_t t0 = h836_t()
#else
#define P836(i, t0) do {} while (0)
#define P836_START(t0) do {} while (0)
#endif
#include "h641_exec.hpp"
inline int g940_bits() { static const int v = [] { const char* e = getenv("H940"); return e ? atoi(e) : 0; }(); return v; }
inline int g941_bits() { static const int v = [] { const char* e = getenv("H941"); return e ? atoi(e) : 0; }(); return v; }
static long g941_st[8];
struct G941St { ~G941St() { if (getenv("H941_STATS")) fprintf(stderr, "H941_STATS masked evals %ld fallback(g2/g3) %ld fallback(g1) %ld | H940 fused initial passes %ld wasted(ndp==0) %ld | sum npp %ld sum ndp %ld\n", g941_st[0], g941_st[1], g941_st[2], g941_st[4], g941_st[5], g941_st[6], g941_st[7]); } };
static G941St g941_stp;
inline int g942_bits() { static const int v = [] { const char* e = getenv("H942"); return e ? atoi(e) : 0; }(); return v; }
inline int g842g_bits() { static const int v = [] { const char* e = getenv("H842"); return e ? atoi(e) : 0; }(); return v; }
inline int g836_bits() { static const int v = [] { const char* e = getenv("H836"); if (!e) return 0; const int x = atoi(e); return x == 1 ? 255 : x; }(); return v; }
inline void g836_copy(float* __restrict d, const float* __restrict s, size_t nf) {
  for (size_t k = 0; k < nf; k += 8) { const float32x4_t a = vld1q_f32(s + k), b = vld1q_f32(s + k + 4); vst1q_f32(d + k, a); vst1q_f32(d + k + 4, b); }
}
#ifndef G644_ON_DEFINED
#define G644_ON_DEFINED
inline bool g644_on() { static const bool v = getenv("H644") && atoi(getenv("H644")) != 0; return v; }
#endif
inline void g836_pf1(const GNInputs& in, uint32_t tree) {
  __builtin_prefetch(in.hdr + 4 * (size_t)tree);
  const char* sg = (const char*)(in.segs + 36 * (size_t)tree); __builtin_prefetch(sg); __builtin_prefetch(sg + 64); __builtin_prefetch(sg + 143);
  const char* th = (const char*)(in.th0 + 32 * (size_t)tree); __builtin_prefetch(th); __builtin_prefetch(th + 64);
}
inline void g836_pf2(const GNInputs& in, uint32_t tree, int Q) {
  const uint32_t* h = in.hdr + 4 * (size_t)tree;
  const uint32_t* sg = in.segs + 36 * (size_t)tree;
  const uint32_t x = h[0], nS = h[1] >> 16, z = h[2], ni = sg[Q];
  const char* t = (const char*)(in.tok + 4 * (size_t)x);
  for (uint32_t b = 0; b < 16 * ni; b += 64) __builtin_prefetch(t + b);
  if (ni) __builtin_prefetch(t + 16 * ni - 1);
  const char* sp = (const char*)(in.sprog + z);
  __builtin_prefetch(sp); if (nS > 32) __builtin_prefetch(sp + 2 * nS - 1);
}
#ifdef H941_CENSUS
struct G941Cell { long items = 0, polishes = 0, duals = 0, tok = 0, tanglanes = 0, livelanes = 0, thetafree = 0, unary = 0, unary_live = 0, unary_lanes = 0, unary_livelanes = 0;
                  long roots = 0, rootlanes = 0, rootlive = 0, rbw = 0, nfd = 0, nfd_risk = 0, guardfail_g2 = 0, tokmismatch = 0; };
static std::map<std::pair<int,int>, G941Cell> g941_cells;
static int g941_nf = 0, g941_risk = 0;
static inline void h941_chk(uint32_t op, const float* acc0) {
  uint32x4_t nf = vdupq_n_u32(0); const float32x4_t vinf = vdupq_n_f32(INFINITY);
  for (int p = 0; p < 48; p += 4) nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(vld1q_f32(acc0 + p), vinf)));
  if (vmaxvq_u32(nf)) g941_nf = 1;
  if (op == 10) { for (int p = 0; p < 48; p++) if (acc0[p] < -80.0f) g941_risk = 1; }
  else if (op == 13 || op == 14) { for (int p = 0; p < 48; p++) if (acc0[p] == 0.0f) g941_risk = 1; }
}
struct G941Rep { ~G941Rep() {
  long T_tok = 0, T_lane = 0, T_live = 0, T_unl = 0, T_unlive = 0, T_rl = 0, T_rlive = 0, T_duals = 0, T_nf = 0, T_risk = 0, T_free = 0, T_rbw = 0, T_mis = 0, T_items = 0;
  fprintf(stderr, "H941CENSUS cell(Q,nth) items duals tokens kappa_tok(live tangent lanes/dense) thetafree_tok_share kappa_unary kappa_root guardfail_g2 risk_d rbw tokmismatch\n");
  for (auto& kv : g941_cells) { const G941Cell& c = kv.second;
    fprintf(stderr, "H941CELL %2d %2d %7ld %7ld %9ld %.4f %.4f %.4f %.4f %ld %ld %ld %ld\n", kv.first.first, kv.first.second, c.items, c.duals, c.tok, c.tanglanes ? (double)c.livelanes / c.tanglanes : 0.0,
      c.tok ? (double)c.thetafree / c.tok : 0.0, c.unary_lanes ? (double)c.unary_livelanes / c.unary_lanes : 0.0, c.rootlanes ? (double)c.rootlive / c.rootlanes : 0.0, c.guardfail_g2, c.nfd_risk, c.rbw, c.tokmismatch);
    T_tok += c.tok; T_lane += c.tanglanes; T_live += c.livelanes; T_unl += c.unary_lanes; T_unlive += c.unary_livelanes; T_rl += c.rootlanes; T_rlive += c.rootlive; T_duals += c.duals; T_nf += c.guardfail_g2; T_risk += c.nfd_risk; T_free += c.thetafree; T_rbw += c.rbw; T_mis += c.tokmismatch; T_items += c.items; }
  fprintf(stderr, "H941TOTAL items %ld duals %ld tokens %ld kappa_tok %.4f (live %ld / dense %ld) unary kappa %.4f thetafree_tok_share %.4f kappa_root %.4f guardfail_g2 %ld (%.3f%% of duals) risk_d %ld rbw %ld tokmismatch %ld\n",
    T_items, T_duals, T_tok, (double)T_live / T_lane, T_live, T_lane, T_unl ? (double)T_unlive / T_unl : 0.0, (double)T_free / T_tok, (double)T_rlive / T_rl, T_nf, 100.0 * T_nf / T_duals, T_risk, T_rbw, T_mis); } };
static G941Rep g941_rep;
#endif
#ifdef H940_ZCHECK
struct G940Z { long evals = 0, roots = 0, badevals = 0, badroots = 0, opev[16] = {0}, badop[16] = {0};
  ~G940Z() { fprintf(stderr, "ZCHECK dual evals %ld roots %ld | mismatching evals %ld roots %ld\n", evals, roots, badevals, badroots);
    for (int o = 0; o < 16; o++) if (opev[o]) fprintf(stderr, "  op %2d present in %ld evals, mismatching evals containing it %ld\n", o, opev[o], badop[o]); } };
static G940Z g940_z;
#endif
namespace {
inline float g331_taud() { static const float v = getenv("H331_TAUD") ? (float)atof(getenv("H331_TAUD")) : 1.0f; return v; }
constexpr int G250_NP = 48;
constexpr int G250_MAXQ = 16;
constexpr int G250_MAXC = 25;
constexpr float G250_INF_ = INFINITY;
constexpr int G250_MAXN = 24, G250_MAXN1 = G250_MAXN + 1;

inline float G250_dot48(const float* x, const float* y) {
  float32x4_t a0 = vdupq_n_f32(0.f), a1 = vdupq_n_f32(0.f), a2 = vdupq_n_f32(0.f);
  for (int k = 0; k < G250_NP; k += 12) {
    a0 = vfmaq_f32(a0, vld1q_f32(x + k), vld1q_f32(y + k));
    a1 = vfmaq_f32(a1, vld1q_f32(x + k + 4), vld1q_f32(y + k + 4));
    a2 = vfmaq_f32(a2, vld1q_f32(x + k + 8), vld1q_f32(y + k + 8));
  }
  return vaddvq_f32(vaddq_f32(vaddq_f32(a0, a1), a2));
}
inline void G250_axpy48(float* y, float r, const float* x) {
  const float32x4_t vr = vdupq_n_f32(r);
  for (int k = 0; k < G250_NP; k += 4) vst1q_f32(y + k, vfmsq_f32(vld1q_f32(y + k), vr, vld1q_f32(x + k)));
}
inline void G250_scal48(float* y, float s) {
  const float32x4_t vs = vdupq_n_f32(s);
  for (int k = 0; k < G250_NP; k += 4) vst1q_f32(y + k, vmulq_f32(vld1q_f32(y + k), vs));
}
inline bool G250_fin(float x) { return std::isfinite(x); }
inline float G250_unv(int op, float m) {
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
inline float G250_und(int op, float m, float r) {
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

GNWork::GNWork() {}

const float* GNWork::load(uint32_t w, float* buf, int C) {
  if (w >> 31) return accp_;
  const uint32_t off = w & 0x1FFFFFFFu;
  if ((w >> 30) & 1u) {
    for (int c = 0; c < C; c++) { const float32x4_t v = vdupq_n_f32(S_[(size_t)c * nsc2_ + (off - SB_)]); for (int p = 0; p < G250_NP; p += 4) vst1q_f32(buf + c * G250_NP + p, v); }
    return buf;
  }
  if ((w >> 29) & 1u) return TS_.data() + off;
  return R_.data() + vslot(off);
}

void GNWork::exec(const uint32_t* w4, int C) {
  if (g421_on()) { exec421(w4, C); return; }
  if (g417_on()) { exec417(w4, C); return; }
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = load(w4[2], a_.data(), C);
  const float* b = load(w4[3], b_.data(), C);
  float* acc = (accp_ == acc_.data()) ? acc2_.data() : acc_.data();
  if (op == 0u) {
    for (int p = 0; p < G250_NP; p += 4) vst1q_f32(acc + p, vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < G250_NP; p += 4)
        vst1q_f32(acc + c * G250_NP + p, vfmaq_f32(vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + c * G250_NP + p)), vld1q_f32(a + c * G250_NP + p), vld1q_f32(b + p)));
  } else if (op < 4u) {
    for (int c = 0; c < C; c++)
      for (int p = 0; p < G250_NP; p += 4) {
        const float32x4_t va = vld1q_f32(a + c * G250_NP + p), vb = vld1q_f32(b + c * G250_NP + p);
        vst1q_f32(acc + c * G250_NP + p, op == 1u ? vaddq_f32(va, vb) : op == 2u ? vsubq_f32(va, vb) : vsubq_f32(vb, va));
      }
  } else if (op == 4u) {
    for (int p = 0; p < G250_NP; p += 4) vst1q_f32(acc + p, vdivq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < G250_NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), bc = vld1q_f32(b + c * G250_NP + p), ac = vld1q_f32(a + c * G250_NP + p);
        vst1q_f32(acc + c * G250_NP + p, vdivq_f32(vfmsq_f32(ac, q0, bc), vld1q_f32(b + p)));
      }
  } else if (op == 5u) {
    for (int p = 0; p < G250_NP; p += 4) vst1q_f32(acc + p, vdivq_f32(vld1q_f32(b + p), vld1q_f32(a + p)));
    for (int c = 1; c < C; c++)
      for (int p = 0; p < G250_NP; p += 4) {
        const float32x4_t q0 = vld1q_f32(acc + p), ac = vld1q_f32(a + c * G250_NP + p), bc = vld1q_f32(b + c * G250_NP + p);
        vst1q_f32(acc + c * G250_NP + p, vdivq_f32(vfmsq_f32(bc, q0, ac), vld1q_f32(a + p)));
      }
  } else if (op >= 7u && op <= 15u) {
    alignas(16) float m[G250_NP], d[G250_NP]; int n = G250_NP;
    for (int p = 0; p < G250_NP; p += 4) vst1q_f32(m + p, vmulq_f32(vld1q_f32(a + p), vld1q_f32(b + p)));
    switch (op) {
      case 7: vvsincosf(acc, d, m, &n); break;
      case 8: { vvsincosf(d, acc, m, &n); for (int p = 0; p < G250_NP; p += 4) vst1q_f32(d + p, vnegq_f32(vld1q_f32(d + p))); } break;
      case 9: vvexpf(acc, m, &n); std::memcpy(d, acc, sizeof d); break;
      case 10: { alignas(16) float am[G250_NP]; for (int p = 0; p < G250_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(am + p, vabsq_f32(vm)); vst1q_f32(d + p, vdivq_f32(vdupq_n_f32(1.0f), vm)); } vvlogf(acc, am, &n); } break;
      case 11: for (int p = 0; p < G250_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vm, vm)); vst1q_f32(d + p, vmulq_f32(vdupq_n_f32(2.0f), vm)); } break;
      case 12: for (int p = 0; p < G250_NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vmulq_f32(vm, vm), vm)); vst1q_f32(d + p, vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vm), vm)); } break;
      case 13: { alignas(16) float am[G250_NP]; for (int p = 0; p < G250_NP; p += 4) vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vvsqrtf(acc, am, &n);
                 for (int p = 0; p < G250_NP; p++) d[p] = copysignf(0.5f, m[p]) / acc[p]; } break;
      case 14: { alignas(16) float am[G250_NP], e3[G250_NP]; for (int p = 0; p < G250_NP; p += 4) { vst1q_f32(am + p, vabsq_f32(vld1q_f32(m + p))); vst1q_f32(e3 + p, vdupq_n_f32(1.0f / 3.0f)); }
                 vvpowf(acc, e3, am, &n); for (int p = 0; p < G250_NP; p++) { acc[p] = copysignf(acc[p], m[p]); d[p] = 1.0f / (3.0f * acc[p] * acc[p]); } } break;
      default: for (int p = 0; p < G250_NP; p += 4) { vst1q_f32(acc + p, vnegq_f32(vld1q_f32(m + p))); vst1q_f32(d + p, vdupq_n_f32(-1.0f)); } break;
    }
    for (int c = 1; c < C; c++)
      for (int p = 0; p < G250_NP; p += 4) {
        const float32x4_t a0 = vld1q_f32(a + p), b0 = vld1q_f32(b + p), ac = vld1q_f32(a + c * G250_NP + p), bc = vld1q_f32(b + c * G250_NP + p);
        { const float32x4_t sm = vfmaq_f32(vmulq_f32(a0, bc), ac, b0), r = vmulq_f32(vld1q_f32(d + p), sm);
          vst1q_f32(acc + c * G250_NP + p, g_h790 ? vbslq_f32(vceqq_f32(sm, vdupq_n_f32(0.f)), vdupq_n_f32(0.f), r) : r); }
      }
  } else {
    std::memcpy(acc, a, sizeof(float) * (size_t)C * G250_NP);
  }
  accp_ = acc;
  if (dst != 0xFFFFFFFFu) std::memcpy(R_.data() + vslot(dst), acc, sizeof(float) * (size_t)C * G250_NP);
}

void GNWork::exec417(const uint32_t* w4, int C) {
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = nullptr; const float* b = nullptr; float as[G250_MAXC], bs[G250_MAXC]; bool sa = false, sb = false;
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
  if (sa) { if (sb) H417_body<G250_NP, true, true>(op, a, as, b, bs, acc, C); else H417_body<G250_NP, true, false>(op, a, as, b, bs, acc, C); }
  else { if (sb) H417_body<G250_NP, false, true>(op, a, as, b, bs, acc, C); else H417_body<G250_NP, false, false>(op, a, as, b, bs, acc, C); }
  accp_ = acc;
  if (R && acc != R) std::memcpy(R, acc, sizeof(float) * (size_t)C * G250_NP);
}

template <bool SA, bool SBB>
__attribute__((noinline)) static void H836_unary(uint32_t op, const float* a, const float* asp, const float* b, const float* bsp, float* acc, int C, size_t ns) {
  float as[G250_MAXC], bs[G250_MAXC];
  if (SA) for (int c = 0; c < C; c++) as[c] = asp[(size_t)c * ns];
  if (SBB) for (int c = 0; c < C; c++) bs[c] = bsp[(size_t)c * ns];
  H641X_body<G250_NP, SA, SBB>(op, a, as, b, bs, acc, C);
}
template <bool SA, bool SBB>
static inline void H836_simple(uint32_t op, const float* a, const float* asp, const float* b, const float* bsp, float* acc, int C, size_t ns) {
  constexpr int NP = G250_NP;
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
    else g836_copy(acc, a, (size_t)C * NP);
  }
#undef LA
#undef LB
}

void GNWork::exec421(const uint32_t* w4, int C) {
  P836_START(xt_);
  const uint32_t op = w4[0], dst = w4[1];
  const float* a = nullptr; const float* b = nullptr; float as[G250_MAXC], bs[G250_MAXC]; bool sa = false, sb = false;
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
#ifdef H836_PROF
  { const uint64_t _t = h836_t(); g836_res_t[C == 1 ? 0 : 1] += _t - xt_; xt_ = _t; }
#endif
  if (sa) { if (sb) H641X_body<G250_NP, true, true>(op, a, as, b, bs, acc, C); else H641X_body<G250_NP, true, false>(op, a, as, b, bs, acc, C); }
  else { if (sb) H641X_body<G250_NP, false, true>(op, a, as, b, bs, acc, C); else H641X_body<G250_NP, false, false>(op, a, as, b, bs, acc, C); }
#ifdef H836_PROF
  { const uint64_t _t = h836_t(); g836_op_t[C == 1 ? 0 : 1][op & 15] += _t - xt_; g836_op_n[C == 1 ? 0 : 1][op & 15]++; xt_ = _t; }
#endif
  accp_ = acc;
  if (R && acc != R) { if (g836_bits() & 1) g836_copy(R, acc, (size_t)C * G250_NP); else std::memcpy(R, acc, sizeof(float) * (size_t)C * G250_NP);
#ifdef H836_PROF
    g836_cpc[C == 1 ? 0 : 1][0]++; g836_cpc[C == 1 ? 0 : 1][1] += C;
#endif
  }
}

void GNWork::exec836(const uint32_t* w4, int C) {
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
    if (sa) { if (sb) H836_unary<true, true>(op, a, asp, b, bsp, acc, C, ns); else H836_unary<true, false>(op, a, asp, b, bsp, acc, C, ns); }
    else { if (sb) H836_unary<false, true>(op, a, asp, b, bsp, acc, C, ns); else H836_unary<false, false>(op, a, asp, b, bsp, acc, C, ns); }
  } else {
    if (sa) { if (sb) H836_simple<true, true>(op, a, asp, b, bsp, acc, C, ns); else H836_simple<true, false>(op, a, asp, b, bsp, acc, C, ns); }
    else { if (sb) H836_simple<false, true>(op, a, asp, b, bsp, acc, C, ns); else H836_simple<false, false>(op, a, asp, b, bsp, acc, C, ns); }
  }
  accp_ = acc;
  if (R && acc != R) g836_copy(R, acc, (size_t)C * G250_NP);
}

alignas(16) static const float g941_zrow[48] = {0};
#define H941_PA(c) (((am >> ((c) - 1)) & 1u) ? a + (c) * NP : g941_zrow)
#define H941_PB(c) (((bm >> ((c) - 1)) & 1u) ? b + (c) * NP : g941_zrow)
template <bool SA, bool SBB>
__attribute__((noinline)) static bool H941_simple(uint32_t op, const float* a, const float* asp, uint32_t am, const float* b, const float* bsp, uint32_t bm, uint32_t rm, float* acc, int C, size_t ns) {
  constexpr int NP = G250_NP; (void)C;
  const float32x4_t vinf = vdupq_n_f32(INFINITY); uint32x4_t ok = vdupq_n_u32(~0u);
#define LA0(p) (SA ? vdupq_n_f32(asp[0]) : vld1q_f32(a + (p)))
#define LB0(p) (SBB ? vdupq_n_f32(bsp[0]) : vld1q_f32(b + (p)))
#define LAC(p) (SA ? sa_ : vld1q_f32(pa + (p)))
#define LBC(p) (SBB ? sb_ : vld1q_f32(pb + (p)))
#define H941_COMP_SETUP const int c = __builtin_ctz(mm) + 1; const float* pa = SA ? nullptr : H941_PA(c); const float* pb = SBB ? nullptr : H941_PB(c); \
    const float32x4_t sa_ = SA ? vdupq_n_f32(asp[(size_t)c * ns]) : vdupq_n_f32(0.f), sb_ = SBB ? vdupq_n_f32(bsp[(size_t)c * ns]) : vdupq_n_f32(0.f)
  if (op == 0u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t a0[4], b0[4];
      for (int u = 0; u < 4; u++) { a0[u] = LA0(p + 4 * u); b0[u] = LB0(p + 4 * u); const float32x4_t r = vmulq_f32(a0[u], b0[u]); vst1q_f32(acc + p + 4 * u, r); ok = vandq_u32(ok, vcaltq_f32(r, vinf)); }
      for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP;
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vfmaq_f32(vmulq_f32(a0[u], LBC(p + 4 * u)), LAC(p + 4 * u), b0[u])); }
    }
  } else if (op < 4u) {
    for (int p = 0; p < NP; p += 4) {
      const float32x4_t va = LA0(p), vb = LB0(p); const float32x4_t r = op == 1u ? vaddq_f32(va, vb) : op == 2u ? vsubq_f32(va, vb) : vsubq_f32(vb, va);
      vst1q_f32(acc + p, r); ok = vandq_u32(ok, vcaltq_f32(r, vinf));
    }
    for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP;
      for (int p = 0; p < NP; p += 4) { const float32x4_t va = LAC(p), vb = LBC(p); vst1q_f32(acc + c * NP + p, op == 1u ? vaddq_f32(va, vb) : op == 2u ? vsubq_f32(va, vb) : vsubq_f32(vb, va)); } }
  } else if (op == 4u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t q0[4], b0[4];
      for (int u = 0; u < 4; u++) { b0[u] = LB0(p + 4 * u); q0[u] = vdivq_f32(LA0(p + 4 * u), b0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); ok = vandq_u32(ok, vcaltq_f32(q0[u], vinf)); }
      for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP;
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LAC(p + 4 * u), q0[u], LBC(p + 4 * u)), b0[u])); }
    }
  } else if (op == 5u) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t q0[4], a0[4];
      for (int u = 0; u < 4; u++) { a0[u] = LA0(p + 4 * u); q0[u] = vdivq_f32(LB0(p + 4 * u), a0[u]); vst1q_f32(acc + p + 4 * u, q0[u]); ok = vandq_u32(ok, vcaltq_f32(q0[u], vinf)); }
      for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP;
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vdivq_f32(vfmsq_f32(LBC(p + 4 * u), q0[u], LAC(p + 4 * u)), a0[u])); }
    }
  } else {
    for (int p = 0; p < NP; p += 4) vst1q_f32(acc + p, LA0(p));
    for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP; (void)pb; (void)sb_;
      for (int p = 0; p < NP; p += 4) vst1q_f32(acc + c * NP + p, LAC(p)); }
  }
  return vminvq_u32(ok) != 0;
}
template <bool SA, bool SBB>
__attribute__((noinline)) static bool H941_unary(uint32_t op, const float* a, const float* asp, uint32_t am, const float* b, const float* bsp, uint32_t bm, uint32_t& rm, uint32_t full, float* acc, int C, size_t ns) {
  constexpr int NP = G250_NP; (void)C;
  const float32x4_t vinf = vdupq_n_f32(INFINITY); uint32x4_t ok = vdupq_n_u32(~0u);
  alignas(16) float m[NP], d[NP]; int n = NP;
  for (int p = 0; p < NP; p += 4) vst1q_f32(m + p, vmulq_f32(LA0(p), LB0(p)));
  switch (op) {
    case 7: vvsincosf(acc, d, m, &n); break;
    case 8: { vvsincosf(d, acc, m, &n); for (int p = 0; p < NP; p += 4) vst1q_f32(d + p, vnegq_f32(vld1q_f32(d + p))); } break;
    case 9: vvexpf(acc, m, &n); std::memcpy(d, acc, sizeof d); break;
    case 10: { alignas(16) float am_[NP]; for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(am_ + p, vabsq_f32(vm)); vst1q_f32(d + p, vdivq_f32(vdupq_n_f32(1.0f), vm)); } vvlogf(acc, am_, &n); } break;
    case 11: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vm, vm)); vst1q_f32(d + p, vmulq_f32(vdupq_n_f32(2.0f), vm)); } break;
    case 12: for (int p = 0; p < NP; p += 4) { const float32x4_t vm = vld1q_f32(m + p); vst1q_f32(acc + p, vmulq_f32(vmulq_f32(vm, vm), vm)); vst1q_f32(d + p, vmulq_f32(vmulq_f32(vdupq_n_f32(3.0f), vm), vm)); } break;
    case 13: { alignas(16) float am_[NP]; for (int p = 0; p < NP; p += 4) vst1q_f32(am_ + p, vabsq_f32(vld1q_f32(m + p))); vvsqrtf(acc, am_, &n);
               for (int p = 0; p < NP; p++) d[p] = copysignf(0.5f, m[p]) / acc[p]; } break;
    case 14: { alignas(16) float am_[NP], e3[NP]; for (int p = 0; p < NP; p += 4) { vst1q_f32(am_ + p, vabsq_f32(vld1q_f32(m + p))); vst1q_f32(e3 + p, vdupq_n_f32(1.0f / 3.0f)); }
               vvpowf(acc, e3, am_, &n); for (int p = 0; p < NP; p++) { acc[p] = copysignf(acc[p], m[p]); d[p] = 1.0f / (3.0f * acc[p] * acc[p]); } } break;
    default: for (int p = 0; p < NP; p += 4) { vst1q_f32(acc + p, vnegq_f32(vld1q_f32(m + p))); vst1q_f32(d + p, vdupq_n_f32(-1.0f)); } break;
  }
  for (int p = 0; p < NP; p += 4) ok = vandq_u32(ok, vcaltq_f32(vld1q_f32(acc + p), vinf));
  if (!vminvq_u32(ok)) return false;
  if (op == 10u || op == 13u || op == 14u) {
    uint32x4_t okd = vdupq_n_u32(~0u);
    for (int p = 0; p < NP; p += 4) okd = vandq_u32(okd, vcaltq_f32(vld1q_f32(d + p), vinf));
    if (!vminvq_u32(okd)) rm = full;
  }
  if (rm) {
    for (int p = 0; p < NP; p += 16) {
      float32x4_t a0[4], b0[4], dd[4];
      for (int u = 0; u < 4; u++) { a0[u] = LA0(p + 4 * u); b0[u] = LB0(p + 4 * u); dd[u] = vld1q_f32(d + p + 4 * u); }
      for (uint32_t mm = rm; mm; mm &= mm - 1) { H941_COMP_SETUP;
        for (int u = 0; u < 4; u++) vst1q_f32(acc + c * NP + p + 4 * u, vmulq_f32(dd[u], vfmaq_f32(vmulq_f32(a0[u], LBC(p + 4 * u)), LAC(p + 4 * u), b0[u]))); }
    }
  }
  return true;
}
#undef LA0
#undef LB0
#undef LAC
#undef LBC
#undef H941_COMP_SETUP

#define H941_SID(off) (cd_ == 1u ? (off) / 48u : cd_ == 2u ? (off) / 96u : cd_ == 3u ? (off) / 144u : (off) / (cd_ * 48u))
inline bool GNWork::opmask941(uint32_t w, uint32_t& m) {
  if (w >> 31) { m = accm_; return true; }
  const uint32_t off = w & 0x1FFFFFFFu;
  if ((w >> 30) & 1u) { const uint32_t sx = off - SB_; if (sx >= smask_.size()) return false; m = smask_[sx]; return true; }
  if ((w >> 29) & 1u) { m = 0u; return true; }
  const uint32_t sid = H941_SID(off); if (sid >= slotst_.size()) return false;
  const uint64_t e = slotst_[sid]; if ((uint32_t)(e >> 32) != gen941_) return false;
  m = (uint32_t)e; return true;
}
bool GNWork::exec940(const uint32_t* w4, int C) {
  const uint32_t op = w4[0], dst = w4[1];
  uint32_t am, bm;
  if (!opmask941(w4[2], am) || !opmask941(w4[3], bm)) return false;
  uint32_t rm = (op == 6u || op > 15u) ? am : (am | bm);
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
  bool ok;
  if (op >= 7u && op <= 15u) {
    if (sa) { if (sb) ok = H941_unary<true, true>(op, a, asp, am, b, bsp, bm, rm, full941_, acc, C, ns); else ok = H941_unary<true, false>(op, a, asp, am, b, bsp, bm, rm, full941_, acc, C, ns); }
    else { if (sb) ok = H941_unary<false, true>(op, a, asp, am, b, bsp, bm, rm, full941_, acc, C, ns); else ok = H941_unary<false, false>(op, a, asp, am, b, bsp, bm, rm, full941_, acc, C, ns); }
  } else {
    if (sa) { if (sb) ok = H941_simple<true, true>(op, a, asp, am, b, bsp, bm, rm, acc, C, ns); else ok = H941_simple<true, false>(op, a, asp, am, b, bsp, bm, rm, acc, C, ns); }
    else { if (sb) ok = H941_simple<false, true>(op, a, asp, am, b, bsp, bm, rm, acc, C, ns); else ok = H941_simple<false, false>(op, a, asp, am, b, bsp, bm, rm, acc, C, ns); }
  }
  if (!ok) return false;
#ifdef H941_POISON
  for (int c = 1; c < C; c++) if (!((rm >> (c - 1)) & 1u)) std::memset(acc + c * G250_NP, 0xFF, sizeof(float) * G250_NP);
#endif
  accp_ = acc; accm_ = rm;
  if (dst != 0xFFFFFFFFu) {
    const uint32_t sid = H941_SID(dst); if (sid >= slotst_.size()) return false;
    slotst_[sid] = ((uint64_t)gen941_ << 32) | rm;
    if (acc != R) {
      g836_copy(R, acc, G250_NP);
      for (uint32_t mm = rm; mm; mm &= mm - 1) { const int c = __builtin_ctz(mm) + 1; g836_copy(R + c * G250_NP, acc + c * G250_NP, G250_NP); }
#ifdef H941_POISON
      for (int c = 1; c < C; c++) if (!((rm >> (c - 1)) & 1u)) std::memset(R + c * G250_NP, 0xFF, sizeof(float) * G250_NP);
#endif
    }
  }
  return true;
}

void GNWork::eval(const GNInputs& in, const Tree& T, const uint32_t* sg, const float* TH, const float* tvec, int C, int Q) {
  rootm_valid_ = false;
  if (C > 1 && C == (int)T.nth + 1 && (g941_bits() & 4) && (g941_bits() & 1) && !g644_on() && !h313b::g_on && (g836_bits() & 16) && g421_on() && !g_h790 && T.nth <= 24 && Q <= 16) {
    bool fin = true; for (uint32_t j = 0; j < T.nth; j++) if (!(fabsf(TH[j]) < INFINITY)) fin = false;
    if (!fin) g941_st[2]++;
    else {
      if (!smask_ok_) {
        smask_.resize(T.nS + 2);
        for (uint32_t j = 0; j < T.nS; j++) { const uint32_t w = in.sprog[T.z + j]; smask_[j] = ((w & 63u) == 1u) ? (1u << ((w >> 6) & 31u)) : 0u; }
        smask_[T.nS] = smask_[T.nS + 1] = 0u; smask_ok_ = true; full941_ = (1u << T.nth) - 1u;
      }
      { const size_t nslot = (size_t)per_ / (cd_ * 48u) + 2; if (slotst_.size() < nslot) slotst_.resize(nslot, 0); }
      gen941_++; accm_ = 0;
      g941_st[0]++;
      if (eval_impl(in, T, sg, TH, tvec, C, Q, true)) { rootm_valid_ = true; return; }
      g941_st[1]++;
    }
  }
  eval_impl(in, T, sg, TH, tvec, C, Q, false);
}

bool GNWork::eval_impl(const GNInputs& in, const Tree& T, const uint32_t* sg, const float* TH, const float* tvec, int C, int Q, bool masked) {
  P836_START(pt_);
  const int pb_ = C == 1 ? 0 : 8; (void)pb_;
  cstride_ = (uint32_t)C;
  { const size_t need = ((size_t)per_ / (cd_ * 48u) + 2) * (size_t)C * 48u + 64; if (R_.size() < need) R_.resize(need); }
  const int b842 = g842g_bits();
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
    if (!ts0_842_) { std::memcpy(TS_.data(), tvec, sizeof(float) * G250_NP); ts0_842_ = true; }
    if (C > tsz842_) { std::memset(TS_.data() + (size_t)tsz842_ * G250_NP, 0, sizeof(float) * G250_NP * (size_t)(C - tsz842_)); tsz842_ = C; }
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
  std::memcpy(TS_.data(), tvec, sizeof(float) * G250_NP);
  std::memset(TS_.data() + G250_NP, 0, sizeof(float) * G250_NP * (size_t)(C - 1));
  }
  bool zacc842 = true;
  if ((b842 & 2) && !g644_on() && !h313b::g_on && sg[1] > sg[0]) { const uint32_t* w0 = in.tok + 4 * (size_t)(T.x + sg[0]); if (!((w0[2] | w0[3]) >> 31)) zacc842 = false; }
  if (zacc842) std::memset(acc_.data(), 0, sizeof(float) * G250_NP * (size_t)C);
#ifdef H842_CNT
  { static struct Cn { long n = 0, sk = 0; ~Cn() { fprintf(stderr, "H842_CNT %s evals %ld acc-zero skipped %ld\n", "g", n, sk); } } cn; cn.n++; if (!zacc842) cn.sk++; }
#endif
  accp_ = acc_.data();
  if (g644_on()) { eval644(in, T, sg, C, Q); return true; }
  P836(pb_ + 0, pt_);
  const bool h836_x = (g836_bits() & 16) && g421_on();
  for (int q = 0; q < Q; q++) {
    for (uint32_t i = sg[q]; i < sg[q + 1]; i++) {
      if (masked) { if (!exec940(in.tok + 4 * (T.x + i), C)) return false; }
      else if (h313b::g_on) { uint32_t w4[4]; h313b::expand_gn(in.tok + 2 * (size_t)(T.x + i), 1, h313b::GnP{cd_ * 48u, T.SB, T.nS}, w4); exec(w4, C); }
      else if (h836_x) { exec836(in.tok + 4 * (T.x + i), C);
#ifdef H941_CENSUS
        if (C > 1) h941_chk(in.tok[4 * (T.x + i)], accp_);
#endif
      }
      else exec(in.tok + 4 * (T.x + i), C);
    }
    P836(pb_ + 1, pt_);
    float* dq = rdst_ ? rdst_ + (size_t)q * rstr_ : roots_.data() + (size_t)q * G250_MAXC * G250_NP;
    if (masked && !opmask941(sg[18 + q], rootm_[q])) return false;
    const float* src = load(sg[18 + q], dq, C);
    if (src != dq) {
      if (masked && (g941_bits() & 2)) {
        g836_copy(dq, src, G250_NP);
        for (uint32_t mm = rootm_[q]; mm; mm &= mm - 1) { const int c = __builtin_ctz(mm) + 1; g836_copy(dq + c * G250_NP, src + c * G250_NP, G250_NP); }
#ifdef H941_POISON
        for (int c = 1; c < C; c++) if (!((rootm_[q] >> (c - 1)) & 1u)) std::memset(dq + c * G250_NP, 0xFF, sizeof(float) * G250_NP);
#endif
      }
      else if (g836_bits() & 1) g836_copy(dq, src, (size_t)C * G250_NP); else std::memcpy(dq, src, sizeof(float) * (size_t)C * G250_NP);
    }
#ifdef H836_PROF
    g836_cpc[2 + (C == 1 ? 0 : 1)][0]++; g836_cpc[2 + (C == 1 ? 0 : 1)][1] += C;
#endif
    P836(pb_ + 2, pt_);
  }
  return true;
}

void GNWork::eval644(const GNInputs& in, const Tree& T, const uint32_t* sg, int C, int Q) {
  K644& K = k644_[C == 1 ? 0 : 1];
  if (K.tok != in.tok || K.sg != sg || K.C != C || K.Q != Q || K.x != T.x || K.cd != cd_ || K.nsc2 != nsc2_ || K.SB != SB_ || K.R != R_.data() ||
      K.S != S_.data() || K.TS != TS_.data() || K.acc != acc_.data() || K.acc2 != acc2_.data() || K.h313 != h313b::g_on) {
    K.tok = in.tok; K.sg = sg; K.C = C; K.Q = Q; K.x = T.x; K.cd = cd_; K.nsc2 = nsc2_; K.SB = SB_; K.R = R_.data(); K.S = S_.data();
    K.TS = TS_.data(); K.acc = acc_.data(); K.acc2 = acc2_.data(); K.h313 = h313b::g_on;
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
        if (h313b::g_on) { h313b::expand_gn(in.tok + 2 * (size_t)(T.x + i), 1, h313b::GnP{cd_ * 48u, T.SB, T.nS}, w4b); w4 = w4b; } else w4 = in.tok + 4 * (T.x + i);
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
      float as[G250_MAXC], bs[G250_MAXC];
      if (r.sa) for (int c = 0; c < C; c++) as[c] = r.a[(size_t)c * ns];
      if (r.sb) for (int c = 0; c < C; c++) bs[c] = r.b[(size_t)c * ns];
      if (r.sa) { if (r.sb) H641X_body<G250_NP, true, true>(r.op, r.a, as, r.b, bs, r.acc, C); else H641X_body<G250_NP, true, false>(r.op, r.a, as, r.b, bs, r.acc, C); }
      else { if (r.sb) H641X_body<G250_NP, false, true>(r.op, r.a, as, r.b, bs, r.acc, C); else H641X_body<G250_NP, false, false>(r.op, r.a, as, r.b, bs, r.acc, C); }
      if (r.cp) std::memcpy(r.cp, r.acc, sizeof(float) * (size_t)C * G250_NP);
    }
    float* dq = roots_.data() + (size_t)q * G250_MAXC * G250_NP;
    if (K.rsc[q]) {
      for (int c = 0; c < C; c++) { const float32x4_t v = vdupq_n_f32(K.root[q][(size_t)c * ns]); for (int p = 0; p < G250_NP; p += 4) vst1q_f32(dq + c * G250_NP + p, v); }
    } else if (K.root[q] != dq) std::memcpy(dq, K.root[q], sizeof(float) * (size_t)C * G250_NP);
  }
  accp_ = K.accp_end;
}

inline void H941_dual_setup(int Q, int nth, const float* roots, int rstride, const float* a, const float* E, float* V, float* on, const uint32_t* rootm) {
  const float32x4_t vinf = vdupq_n_f32(INFINITY);
  for (int j = 0; j < nth; j++) {
    float32x4_t v[12];
    for (int u = 0; u < 12; u++) v[u] = vdupq_n_f32(0.0f);
    for (int q = 0; q < Q; q++) {
      if (!((rootm[q] >> j) & 1u)) continue;
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
template <bool TRIM>
inline void H942_mgs(int n, float* const* col, float* tr, int L, int Lp, float tau, bool dual, const float* on, float& pred) {
  pred = 0.0f;
  if (n < 2) return;
  const int G = Lp >> 2;
  float nrm;
  { const float* c0 = col[0]; float32x4_t x[12]; for (int u = 0; u < 12; u++) x[u] = vld1q_f32(c0 + 4 * u); nrm = sqrtf(H64_dot48r(x, x)); }
  for (int k = 0; k < n - 1; k++) {
    float* ck = col[k];
    const bool keep = nrm > tau * on[k];
    const float inv = keep ? 1.0f / nrm : 0.0f;
    const bool fin_inv = TRIM && (fabsf(inv) < INFINITY);
    float* trk = tr + k * Lp;
    const int Gk = TRIM ? (((k >> 2) + 1) < G ? ((k >> 2) + 1) : G) : G;
    { const float32x4_t vi = vdupq_n_f32(inv); const int g1 = fin_inv ? Gk : G; for (int g = 0; g < g1; g++) vst1q_f32(trk + 4 * g, vmulq_f32(vld1q_f32(trk + 4 * g), vi)); }
    float32x4_t q[12];
    { const float32x4_t vs = vdupq_n_f32(inv); for (int u = 0; u < 12; u++) { q[u] = vmulq_f32(vld1q_f32(ck + 4 * u), vs); vst1q_f32(ck + 4 * u, q[u]); } }
    for (int j = k + 1; j < n; j++) {
      float* cj = col[j];
      float32x4_t x[12];
      for (int u = 0; u < 12; u++) x[u] = vld1q_f32(cj + 4 * u);
      const float r = H64_dot48r(q, x);
      { float* trj = tr + j * Lp; const float32x4_t vnr = vdupq_n_f32(-r);
        const int g1 = (TRIM && fin_inv && fabsf(r) < INFINITY) ? Gk : G;
        for (int g = 0; g < g1; g++) vst1q_f32(trj + 4 * g, vfmaq_f32(vld1q_f32(trj + 4 * g), vnr, vld1q_f32(trk + 4 * g))); }
      if (dual && j == n - 1) pred = fmaf(r, r, pred);
      const float32x4_t vr = vdupq_n_f32(r);
      for (int u = 0; u < 12; u++) { x[u] = vfmsq_f32(x[u], vr, q[u]); vst1q_f32(cj + 4 * u, x[u]); }
      if (j == k + 1 && k + 1 < n - 1) nrm = sqrtf(H64_dot48r(x, x));
    }
  }
}
inline void H942_tr_init(float* tr, int L, int Lp) {
  const int G = Lp >> 2; const float32x4_t z = vdupq_n_f32(0.f);
  const float oh[4][4] = {{1.f, 0.f, 0.f, 0.f}, {0.f, 1.f, 0.f, 0.f}, {0.f, 0.f, 1.f, 0.f}, {0.f, 0.f, 0.f, 1.f}};
  for (int r = 0; r <= L; r++) for (int g = 0; g < G; g++) vst1q_f32(tr + r * Lp + 4 * g, (r < L && (r >> 2) == g) ? vld1q_f32(oh[r & 3]) : z);
}
namespace {
void G250_mgs(int n, float* cols, float* tr, int L, float tau, bool dual, const float* on, float& pred) {
  pred = 0.0f;
  for (int k = 0; k < n - 1; k++) {
    float* ck = cols + k * G250_NP;
    const float nrm = sqrtf(G250_dot48(ck, ck));
    const bool keep = nrm > tau * on[k];
    const float inv = keep ? 1.0f / nrm : 0.0f;
    for (int l = 0; l < L; l++) tr[k * L + l] *= inv;
    G250_scal48(ck, inv);
    float rj[G250_MAXN1 + G250_MAXQ + 1];
    for (int j = k + 1; j < n; j++) rj[j] = G250_dot48(ck, cols + j * G250_NP);
    for (int j = k + 1; j < n; j++) {
      const float r = rj[j];
      for (int l = 0; l < L; l++) tr[j * L + l] = fmaf(-r, tr[k * L + l], tr[j * L + l]);
      if (dual && j == n - 1) pred = fmaf(r, r, pred);
      G250_axpy48(cols + j * G250_NP, r, ck);
    }
  }
}
}

namespace {
struct G250_PassS { float E[G250_MAXQ * G250_NP]; float Y[G250_NP]; float a[G250_MAXQ]; float rss; };
}
void GNWork::polish(const GNInputs& in, uint32_t tree, int Q, float tau, float flr, int itmax, GNResult& out) {
  const uint32_t* h = in.hdr + 4 * tree;
  Tree T{h[0], h[1] >> 16, h[2], h[3] & 255u, (h[3] >> 8) & 0xFFFFu, h[3] >> 24};
  const uint32_t* sg = in.segs + 36 * tree;
  nsc2_ = T.nsc2; SB_ = T.SB; per_ = in.per;
  cd_ = (T.nth <= 4) ? 1 + (T.nth + 1) / 2 : 3;
  if (S_.size() < (size_t)G250_MAXC * T.nsc2 + 8) S_.resize((size_t)G250_MAXC * T.nsc2 + 8);
  const size_t ts_sz842 = TS_.size();
  if (TS_.size() < (size_t)G250_MAXC * G250_NP) TS_.resize((size_t)G250_MAXC * G250_NP);
  if (acc_.size() < (size_t)G250_MAXC * G250_NP) acc_.resize((size_t)G250_MAXC * G250_NP);
  if (acc2_.size() < (size_t)G250_MAXC * G250_NP) acc2_.resize((size_t)G250_MAXC * G250_NP);
  if (roots_.size() < (size_t)G250_MAXQ * G250_MAXC * G250_NP) roots_.resize((size_t)G250_MAXQ * G250_MAXC * G250_NP);
  if (a_.size() < (size_t)G250_MAXC * G250_NP) a_.resize((size_t)G250_MAXC * G250_NP);
  if (b_.size() < (size_t)G250_MAXC * G250_NP) b_.resize((size_t)G250_MAXC * G250_NP);
  smask_ok_ = false; rootm_valid_ = false;
  sb842_ = 1; ts0_842_ = false;
  if (TS_.size() != ts_sz842) tsz842_ = 1;
#ifdef H842_POISON
  { std::memset(S_.data(), 0xFF, sizeof(float) * S_.size()); std::memset(acc_.data(), 0xFF, sizeof(float) * acc_.size()); std::memset(acc2_.data(), 0xFF, sizeof(float) * acc2_.size()); std::memset(TS_.data(), 0xFF, sizeof(float) * G250_NP); }
#endif
  const float* tvec = in.tlin;
  const int nth = (int)T.nth;
  const int Cd = nth + 1;
#ifdef H941_CENSUS
  g941_cells[{Q, nth}].items++;
#endif

  const void* tang_of = nullptr;
  auto plain_impl = [&](const float* TH, G250_PassS& P, int Cev) {
    tang_of = nullptr;
    const bool direct = (g836_bits() & 2) && !g644_on() && (g641_on() || (g643_bits() & 1)) && Cev == 1;
    if (direct) { rdst_ = P.E; rstr_ = G250_NP; }
    eval(in, T, sg, TH, tvec, Cev, Q);
    rdst_ = nullptr;
    if (Cev > 1) tang_of = &P;
    P836_START(pp_);
    if (g641_on() || (g643_bits() & 1)) {
      float tr[(G250_MAXQ + 1) * G250_MAXQ], on[G250_MAXQ + 1]; float* col[G250_MAXQ + 1];
      bool bad = false;
      for (int q = 0; q < Q; q++) {
        const float* rq = direct ? (const float*)(P.E + q * G250_NP) : roots_.data() + (size_t)q * G250_MAXC * G250_NP;
        float* e = P.E + q * G250_NP; col[q] = e;
        uint32x4_t nf = vdupq_n_u32(0);
        for (int p = 0; p < G250_NP; p += 4) { const float32x4_t v = vld1q_f32(rq + p); if (!direct) vst1q_f32(e + p, v); nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(v, vdupq_n_f32(INFINITY)))); }
        if (vmaxvq_u32(nf)) bad = true;
        on[q] = sqrtf(H64_dot48(e, e));
      }
      std::memcpy(P.Y, in.yv, sizeof(float) * G250_NP); col[Q] = P.Y;
      float pred;
      if ((g942_bits() & 1) && (g643_bits() & 1) && !(g942_bits() & 8)) {
        const int Lp = (Q + 3) & ~3; H942_tr_init(tr, Q, Lp);
#ifdef H942_DBG
        { float cc[(G250_MAXQ + 1) * G250_NP], trr[(G250_MAXQ + 1) * G250_MAXQ]; float* cp[G250_MAXQ + 1]; for (int c = 0; c <= Q; c++) { std::memcpy(cc + c * G250_NP, col[c], 4 * G250_NP); cp[c] = cc + c * G250_NP; }
          for (int k = 0; k < (Q + 1) * Q; k++) trr[k] = 0.0f; for (int q = 0; q < Q; q++) trr[q * Q + q] = 1.0f; float pr;
          H64_mgs(Q + 1, cp, trr, Q, tau, false, on, pr);
          float c2[(G250_MAXQ + 1) * G250_NP]; float* cq[G250_MAXQ + 1]; for (int c = 0; c <= Q; c++) { std::memcpy(c2 + c * G250_NP, col[c], 4 * G250_NP); cq[c] = c2 + c * G250_NP; }
          float t2[(G250_MAXQ + 1) * G250_MAXQ]; std::memcpy(t2, tr, sizeof t2); float p2;
          H942_mgs<false>(Q + 1, cq, t2, Q, Lp, tau, false, on, p2);
          static int nd = 0; bool bad2 = false; for (int q = 0; q < Q; q++) if (std::memcmp(&trr[Q * Q + q], &t2[Q * Lp + q], 4)) bad2 = true;
          for (int c = 0; c <= Q; c++) if (std::memcmp(cc + c * G250_NP, c2 + c * G250_NP, 4 * G250_NP)) bad2 = true;
          if (bad2 && nd++ < 5) { fprintf(stderr, "H942_DBG mismatch Q %d Lp %d\n", Q, Lp); for (int c = 0; c <= Q; c++) for (int p = 0; p < G250_NP; p++) if (std::memcmp(cc + c * G250_NP + p, c2 + c * G250_NP + p, 4)) { unsigned u1, u2; std::memcpy(&u1, cc + c * G250_NP + p, 4); std::memcpy(&u2, c2 + c * G250_NP + p, 4); fprintf(stderr, "  col %d p %d ref %08x new %08x\n", c, p, u1, u2); break; } for (int q = 0; q < Q; q++) { unsigned u1, u2; std::memcpy(&u1, &trr[Q * Q + q], 4); std::memcpy(&u2, &t2[Q * Lp + q], 4); fprintf(stderr, "  a[%d] ref %08x new %08x\n", q, u1, u2); } } }
#endif
        P836(3, pp_);
        if (g942_bits() & 2) H942_mgs<true>(Q + 1, col, tr, Q, Lp, tau, false, on, pred); else H942_mgs<false>(Q + 1, col, tr, Q, Lp, tau, false, on, pred);
        P836(4, pp_);
        for (int q = 0; q < Q; q++) P.a[q] = -tr[Q * Lp + q];
      } else {
      for (int k = 0; k < (Q + 1) * Q; k++) tr[k] = 0.0f;
      for (int q = 0; q < Q; q++) tr[q * Q + q] = 1.0f;
      P836(3, pp_);
      H64_mgs(Q + 1, col, tr, Q, tau, false, on, pred);
      P836(4, pp_);
      for (int q = 0; q < Q; q++) P.a[q] = -tr[Q * Q + q];
      }
      const float rss = H64_dot48(P.Y, P.Y);
      P.rss = (bad || !G250_fin(rss)) ? G250_INF_ : rss;
      return;
    }
    float cols[(G250_MAXQ + 1) * G250_NP], tr[(G250_MAXQ + 1) * G250_MAXQ], on[G250_MAXQ + 1];
    bool bad = false;
    for (int q = 0; q < Q; q++) {
      const float* rq = roots_.data() + (size_t)q * G250_MAXC * G250_NP;
      std::memcpy(cols + q * G250_NP, rq, sizeof(float) * G250_NP);
      uint32x4_t nf = vdupq_n_u32(0);
      for (int p = 0; p < G250_NP; p += 4) nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(vld1q_f32(rq + p), vdupq_n_f32(INFINITY))));
      if (vmaxvq_u32(nf)) bad = true;
    }
    std::memcpy(cols + Q * G250_NP, in.yv, sizeof(float) * G250_NP);
    for (int k = 0; k < Q; k++) on[k] = sqrtf(G250_dot48(cols + k * G250_NP, cols + k * G250_NP));
    for (int k = 0; k < (Q + 1) * Q; k++) tr[k] = 0.0f;
    for (int q = 0; q < Q; q++) tr[q * Q + q] = 1.0f;
    float pred;
    G250_mgs(Q + 1, cols, tr, Q, tau, false, on, pred);
    for (int q = 0; q < Q; q++) { std::memcpy(P.E + q * G250_NP, cols + q * G250_NP, sizeof(float) * G250_NP); P.a[q] = -tr[Q * Q + q]; }
    std::memcpy(P.Y, cols + Q * G250_NP, sizeof(float) * G250_NP);
    const float rss = G250_dot48(P.Y, P.Y);
    P.rss = (bad || !G250_fin(rss)) ? G250_INF_ : rss;
  };
  auto plain = [&](const float* TH, G250_PassS& P) { plain_impl(TH, P, 1); };
#ifdef H941_CENSUS
  auto h941_census = [&](int Qc) {
    G941Cell& cell = g941_cells[{Qc, nth}];
    const uint32_t nS = T.nS; uint32_t smask[512]; if (nS + 2 > 512) return;
    for (uint32_t j = 0; j < nS; j++) { const uint32_t w = in.sprog[T.z + j]; smask[j] = ((w & 63u) == 1u) ? (1u << ((w >> 6) & 31u)) : 0u; }
    smask[nS] = smask[nS + 1] = 0;
    const uint32_t NS = 0xFFFFFFFFu;
    static thread_local std::vector<uint32_t> slotm; const size_t nslot = (size_t)per_ / (cd_ * 48u) + 2; slotm.assign(nslot, NS);
    uint32_t accm = 0; long ntok = 0;
    auto opm = [&](uint32_t w) -> uint32_t {
      if (w >> 31) return accm;
      const uint32_t off = w & 0x1FFFFFFFu;
      if ((w >> 30) & 1u) { const uint32_t s = off - SB_; return s < nS + 2 ? smask[s] : 0xFFFFFFFFu; }
      if ((w >> 29) & 1u) return 0u;
      const size_t sl = off / (cd_ * 48u); if (sl >= nslot) return 0xFFFFFFFFu;
      if (slotm[sl] == NS) cell.rbw++;
      return slotm[sl] == NS ? 0u : slotm[sl];
    };
    const uint32_t full = nth >= 32 ? 0xFFFFFFFFu : ((1u << nth) - 1u);
    for (uint32_t i = sg[0]; i < sg[Qc]; i++) {
      const uint32_t* w4 = in.tok + 4 * (size_t)(T.x + i);
      const uint32_t op = w4[0], dst = w4[1];
      const uint32_t ma = opm(w4[2]), mb = opm(w4[3]);
      const uint32_t rm = (ma | mb) & full;
      if (dst != 0xFFFFFFFFu) { const size_t sl = dst / (cd_ * 48u); if (sl < nslot) slotm[sl] = rm; }
      accm = rm; ntok++;
      cell.tanglanes += nth; cell.livelanes += __builtin_popcount(rm);
      if (rm == 0) cell.thetafree++;
      if (op >= 7u && op <= 15u) { cell.unary++; cell.unary_lanes += nth; cell.unary_livelanes += __builtin_popcount(rm); }
    }
    if (ntok != (long)(sg[Qc] - sg[0])) cell.tokmismatch++;
    cell.tok += ntok; cell.duals++;
    for (int q = 0; q < Qc; q++) { const uint32_t rmq = opm(sg[18 + q]) & full; cell.roots++; cell.rootlanes += nth; cell.rootlive += __builtin_popcount(rmq); }
    if (g941_nf) cell.guardfail_g2++;
    if (g941_risk) cell.nfd_risk++;
  };
#endif
  auto dual = [&](const float* TH, const G250_PassS& P, float* d, float& pred) {
#ifdef H941_CENSUS
    g941_nf = 0; g941_risk = 0;
#endif
    if (tang_of != &P) eval(in, T, sg, TH, tvec, Cd, Q);
    tang_of = nullptr;
    const bool sk1 = (g941_bits() & 1) && rootm_valid_;
#ifdef H941_CENSUS
    h941_census(Q);
#endif
#ifdef H940_ZCHECK
    {
      alignas(16) float z1[G250_MAXQ * G250_NP];
      rdst_ = z1; rstr_ = G250_NP;
      eval(in, T, sg, TH, tvec, 1, Q);
      rdst_ = nullptr;
      uint32_t ops = 0; for (uint32_t i = sg[0]; i < sg[Q]; i++) ops |= 1u << (in.tok[4 * (T.x + i)] & 15);
      g940_z.evals++;
      bool anybad = false;
      for (int q = 0; q < Q; q++) {
        g940_z.roots++;
        if (std::memcmp(z1 + q * G250_NP, roots_.data() + (size_t)q * G250_MAXC * G250_NP, sizeof(float) * G250_NP) != 0) {
          g940_z.badroots++; anybad = true;
          if (g940_z.badroots <= 5) fprintf(stderr, "ZCHECK mismatch tree %u q %d Q %d nth %d opmask %x\n", 0u, q, Q, nth, ops);
        }
      }
      if (anybad) { g940_z.badevals++; for (int o = 0; o < 16; o++) if ((ops >> o) & 1) g940_z.badop[o]++; }
      for (int o = 0; o < 16; o++) if ((ops >> o) & 1) g940_z.opev[o]++;
    }
#endif
    P836_START(pd_);
    if (g643_bits()) {
      float V[(G250_MAXN1) * G250_NP], tr[G250_MAXN1 * G250_MAXN], on[G250_MAXN1]; float* col[G250_MAXN1];
      if (g643_bits() & 2) { if (sk1) H941_dual_setup(Q, nth, roots_.data(), G250_MAXC * G250_NP, P.a, P.E, V, on, rootm_); else H643_dual_setup(Q, nth, roots_.data(), G250_MAXC * G250_NP, P.a, P.E, V, on); }
      else {
        for (int k = 0; k < (nth + 1) * G250_NP; k++) V[k] = 0.0f;
        for (int q = 0; q < Q; q++) {
          const float aq = P.a[q];
          const float* rq = roots_.data() + (size_t)q * G250_MAXC * G250_NP;
          const float32x4_t va = vdupq_n_f32(aq), vinf = vdupq_n_f32(INFINITY);
          for (int j = 0; j < nth; j++) {
            if (sk1 && !((rootm_[q] >> j) & 1u)) continue;
            for (int p = 0; p < G250_NP; p += 4) {
              const float32x4_t t_ = vmulq_f32(va, vld1q_f32(rq + (j + 1) * G250_NP + p));
              const float32x4_t tz = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(t_), vcaltq_f32(t_, vinf)));
              vst1q_f32(V + j * G250_NP + p, vaddq_f32(vld1q_f32(V + j * G250_NP + p), tz));
            }
          }
        }
        for (int j = 0; j < nth; j++) on[j] = sqrtf(G250_dot48(V + j * G250_NP, V + j * G250_NP));
        for (int q = 0; q < Q; q++) {
          const float* qk = P.E + q * G250_NP;
          for (int j = 0; j < nth; j++) { const float r = G250_dot48(qk, V + j * G250_NP); G250_axpy48(V + j * G250_NP, r, qk); }
        }
      }
      std::memcpy(V + nth * G250_NP, P.Y, sizeof(float) * G250_NP);
      const bool v942 = (g942_bits() & 1) && (g643_bits() & 1) && !(g942_bits() & 4); const int Lp = (nth + 3) & ~3;
      if (v942) H942_tr_init(tr, nth, Lp);
      else { for (int k = 0; k < (nth + 1) * nth; k++) tr[k] = 0.0f;
        for (int j = 0; j < nth; j++) tr[j * nth + j] = 1.0f; }
      for (int j = 0; j <= nth; j++) col[j] = V + j * G250_NP;
      P836(5, pd_);
      if (v942) { if (g942_bits() & 2) H942_mgs<true>(nth + 1, col, tr, nth, Lp, tau * g331_taud(), true, on, pred); else H942_mgs<false>(nth + 1, col, tr, nth, Lp, tau * g331_taud(), true, on, pred); }
      else H64_mgs(nth + 1, col, tr, nth, tau * g331_taud(), true, on, pred);
      P836(6, pd_);
      const float rss = G250_dot48(V + nth * G250_NP, V + nth * G250_NP);
      for (int j = 0; j < nth; j++) d[j] = v942 ? -tr[nth * Lp + j] : -tr[nth * nth + j];
      if (!G250_fin(rss)) pred = -1.0f;
      P836(7, pd_);
      return;
    }
    float V[(G250_MAXN1) * G250_NP], tr[G250_MAXN1 * G250_MAXN], on[G250_MAXN1];
    for (int k = 0; k < (nth + 1) * G250_NP; k++) V[k] = 0.0f;
    for (int q = 0; q < Q; q++) {
      const float aq = P.a[q];
      const float* rq = roots_.data() + (size_t)q * G250_MAXC * G250_NP;
      const float32x4_t va = vdupq_n_f32(aq), vinf = vdupq_n_f32(INFINITY);
      for (int j = 0; j < nth; j++) {
        if (sk1 && !((rootm_[q] >> j) & 1u)) {
#ifdef H941_DBG
          { static long nbad = 0; for (int p = 0; p < G250_NP; p++) { const float v = rq[(j + 1) * G250_NP + p]; if (v != 0.0f && std::isfinite(v)) { if (nbad++ < 10) fprintf(stderr, "H941_DBG dead root comp nonzero: Q %d nth %d q %d j %d p %d v %g rootm %x\n", Q, nth, q, j, p, v, rootm_[q]); break; } } }
#endif
          continue;
        }
        for (int p = 0; p < G250_NP; p += 4) {
          const float32x4_t t_ = vmulq_f32(va, vld1q_f32(rq + (j + 1) * G250_NP + p));
          const float32x4_t tz = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(t_), vcaltq_f32(t_, vinf)));
          vst1q_f32(V + j * G250_NP + p, vaddq_f32(vld1q_f32(V + j * G250_NP + p), tz));
        }
      }
    }
    for (int j = 0; j < nth; j++) on[j] = sqrtf(G250_dot48(V + j * G250_NP, V + j * G250_NP));
    for (int q = 0; q < Q; q++) {
      const float* qk = P.E + q * G250_NP;
      for (int j = 0; j < nth; j++) { const float r = G250_dot48(qk, V + j * G250_NP); G250_axpy48(V + j * G250_NP, r, qk); }
    }
    std::memcpy(V + nth * G250_NP, P.Y, sizeof(float) * G250_NP);
    for (int k = 0; k < (nth + 1) * nth; k++) tr[k] = 0.0f;
    for (int j = 0; j < nth; j++) tr[j * nth + j] = 1.0f;
    G250_mgs(nth + 1, V, tr, nth, tau * g331_taud(), true, on, pred);
    const float rss = G250_dot48(V + nth * G250_NP, V + nth * G250_NP);
    for (int j = 0; j < nth; j++) d[j] = -tr[nth * nth + j];
    if (!G250_fin(rss)) pred = -1.0f;
  };

  float th[32] = {0}, cand0[32] = {0}, cand1[32] = {0}, d[32] = {0};
  for (int j = 0; j < nth; j++) th[j] = in.th0[32 * tree + j];
  G250_PassS PB[3];
  int ic = 0;
  G250_PassS* cur = &PB[0];
  uint32_t steps = 0, npp = 1, ndp = 0;
  const bool fused_init = (g940_bits() & 1) && nth >= 1 && itmax > 0;
  plain_impl(th, PB[0], fused_init ? Cd : 1);
  float rcur = PB[0].rss;
  const float rss0 = PB[0].rss;
  bool done = !G250_fin(rcur) || itmax <= 0;
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
      const int cev = ((g940_bits() & 2) && nth >= 1 && it + 1 < itmax) ? Cd : 1;
      plain_impl(cand0, PB[i0], cev);
      const float r0 = PB[i0].rss;
      float r1 = r0;
      int pk = -1;
      if (r0 <= rcur) pk = 0;
      else if (hp < 2) { plain_impl(cand1, PB[i1], cev); r1 = PB[i1].rss; if (r1 <= rcur) pk = 1; }
      if (pk >= 0) {
        rcur = pk ? r1 : r0;
        for (int j = 0; j < nth; j++) th[j] = pk ? cand1[j] : cand0[j];
        ic = pk ? i1 : i0;
        cur = &PB[ic];
        steps++; it++;
        if (it >= itmax || !G250_fin(rcur)) break;
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
  if (fused_init) { g941_st[4]++; if (ndp == 0) g941_st[5]++; }
  g941_st[6] += npp; g941_st[7] += ndp;
}

void GNWork::roots_debug(const GNInputs& in, uint32_t tree, const float* TH, int C, int Q, float* out) {
  const uint32_t* h = in.hdr + 4 * tree;
  Tree T{h[0], h[1] >> 16, h[2], h[3] & 255u, (h[3] >> 8) & 0xFFFFu, h[3] >> 24};
  nsc2_ = T.nsc2; SB_ = T.SB;
  S_.assign((size_t)G250_MAXC * T.nsc2 + 8, 0.0f);
  per_ = in.per; cd_ = (T.nth <= 4) ? 1 + (T.nth + 1) / 2 : 3; R_.assign(((size_t)in.per / (cd_ * 48u) + 2) * (size_t)(G250_MAXC) * 48u + 64, 0.0f);
  TS_.assign((size_t)G250_MAXC * G250_NP, 0.0f); acc_.assign((size_t)G250_MAXC * G250_NP, 0.0f); acc2_.assign((size_t)G250_MAXC * G250_NP, 0.0f);
  roots_.assign((size_t)G250_MAXQ * G250_MAXC * G250_NP, 0.0f); a_.assign((size_t)G250_MAXC * G250_NP, 0.0f); b_.assign((size_t)G250_MAXC * G250_NP, 0.0f);
  sb842_ = 1; tsz842_ = 1; ts0_842_ = false;
  eval(in, T, in.segs + 36 * tree, TH, in.tlin, C, Q);
  for (int q = 0; q < Q; q++) std::memcpy(out + (size_t)q * C * G250_NP, roots_.data() + (size_t)q * G250_MAXC * G250_NP, sizeof(float) * C * G250_NP);
}
