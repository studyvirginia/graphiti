#pragma once
#include <cstdint>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <arm_neon.h>
#include <Accelerate/Accelerate.h>

#ifndef S290_VF
#define S290_VF 1
#endif
#ifndef S290_VFMASK
#define S290_VFMASK 0xFFu
#endif
#ifndef S290_GS
#define S290_GS 1
#endif
namespace s250b {

struct Geo { int S, LPS, P, T, E; };
inline Geo geo_for(int Q, bool s2) {
  if (s2) return {2, 16, 3, 32, 96};
  if (Q <= 5) return {16, 8, 6, 128, 768};
  return {16, 16, 3, 256, 768};
}

struct Inputs {
  const uint32_t* hdr;
  const uint32_t* prog;
  const uint16_t* sprog;
  const float* bred;
  const float* sob;
  const float* yv;
  const uint32_t* segs;
  const float* ctab;
};

template<typename R> struct Work {
  Geo g{};
  int Q = 0, nslot = 0, nsc = 0;
  uint32_t TT = 0;
  std::vector<R> NS;
  std::vector<uint32_t> ix1, ix2;
  std::vector<R> A, B, ACC, ACC2, C;
  R* accp = nullptr;
  std::vector<R> M, D, ybuf;
  std::vector<R> cap;
  std::vector<R> cols;

  void prepare(int Q_, bool s2, int nslot_, int nsc_) {
    Q = Q_; nslot = nslot_; nsc = nsc_; g = geo_for(Q, s2);
    const int E = g.E, T = g.T, P = g.P, S = g.S, LPS = g.LPS;
    TT = (uint32_t)(nslot * E + S * nsc);
    NS.assign((size_t)nslot * E + (size_t)S * nsc + 50, (R)0);
    ix1.resize(E); ix2.resize(E);
    for (int p = 0; p < P; p++)
      for (int l = 0; l < T; l++) {
        const int e = p * T + l, sl = l / LPS, lq = l % LPS;
        ix1[e] = (uint32_t)(nslot * E + sl * nsc);
        ix2[e] = (uint32_t)(TT + lq * P + p);
      }
    A.resize(E); B.resize(E); ACC.resize(E); ACC2.resize(E); C.resize(E); M.resize(E); D.resize(E); ybuf.resize((size_t)S * 48);
    accp = ACC.data();
    cap.resize((size_t)Q * E); cols.resize((size_t)Q * S * 48);
  }

  inline const R* fetch(uint32_t raw, R* buf) const {
    const uint32_t cl = (raw >> 29) & 3u, of = raw & 0x1FFFFFFFu;
    const R* N = NS.data(); const int E = g.E, T = g.T, P = g.P, S = g.S, LPS = g.LPS;
    if (cl == 0u) return N + of;
    if (cl == 1u) {
      const R* sv = N + of + (size_t)nslot * E;
      for (int sl = 0; sl < S; sl++) { const R v = sv[(size_t)sl * nsc]; for (int lq = 0; lq < LPS; lq++) buf[sl * LPS + lq] = v; }
      for (int p = 1; p < P; p++) std::memcpy(buf + p * T, buf, sizeof(R) * T);
      return buf;
    }
    if (cl == 2u) {
      const R* tb = N + TT + of;
      for (int p = 0; p < P; p++) {
        R* row = buf + p * T;
        for (int lq = 0; lq < LPS; lq++) row[lq] = tb[lq * P + p];
        for (int sl = 1; sl < S; sl++) std::memcpy(row + sl * LPS, row, sizeof(R) * LPS);
      }
      return buf;
    }
    const R v = N[of]; for (int e = 0; e < E; e++) buf[e] = v;
    return buf;
  }
};

template<typename R> inline void exec_instr(Work<R>& W, uint32_t x, uint32_t y, uint32_t zraw, uint32_t wraw, const float* ctab) {
  const int E = W.g.E;
  const R* A = (zraw >> 31) ? W.accp : W.fetch(zraw, W.A.data());
  const R* B = W.fetch(wraw, W.B.data());
  R* ACC = (W.accp == W.ACC.data()) ? W.ACC2.data() : W.ACC.data();
  if (x == 0u) {
    for (int e = 0; e < E; e++) ACC[e] = A[e] * B[e];
  } else if (x < 4u) {
    const R cy = (R)ctab[4 * x + 1], cz = (R)ctab[4 * x + 2];
    for (int e = 0; e < E; e++) ACC[e] = cy * std::fma(cz, B[e], A[e]);
  } else {
    const uint32_t op = x & 255u;
    if (S290_VF && std::is_same<R, float>::value && ((op >= 7u && op <= 10u) || op == 13u || op == 14u) && ((S290_VFMASK >> (op - 7u)) & 1u)) {
      float* m = (float*)W.M.data(); float* d = (float*)W.D.data(); float* acc = (float*)ACC; int n = E;
      const float* a = (const float*)A; const float* b = (const float*)B;
      for (int e = 0; e < E; e += 4) {
        const float32x4_t vm = vmulq_f32(vld1q_f32(a + e), vld1q_f32(b + e));
        vst1q_f32(m + e, (op == 7u || op == 8u || op == 9u) ? vm : vabsq_f32(vm));
      }
      switch (op) {
        case 7: vvsincosf(acc, d, m, &n); break;
        case 8: vvsincosf(d, acc, m, &n); break;
        case 9: vvexpf(acc, m, &n); break;
        case 10: vvlogf(acc, m, &n); break;
        case 13: vvsqrtf(acc, m, &n); break;
        default: { for (int e = 0; e < E; e += 4) vst1q_f32(d + e, vdupq_n_f32(1.0f / 3.0f)); vvpowf(acc, d, m, &n);
                   for (int e = 0; e < E; e++) acc[e] = std::copysign(acc[e], a[e] * b[e]); } break;
      }
    } else switch (op) {
      case 4: for (int e = 0; e < E; e++) ACC[e] = A[e] / B[e]; break;
      case 5: for (int e = 0; e < E; e++) ACC[e] = B[e] / A[e]; break;
      case 7: for (int e = 0; e < E; e++) ACC[e] = std::sin(A[e] * B[e]); break;
      case 8: for (int e = 0; e < E; e++) ACC[e] = std::cos(A[e] * B[e]); break;
      case 9: for (int e = 0; e < E; e++) ACC[e] = std::exp(A[e] * B[e]); break;
      case 10: for (int e = 0; e < E; e++) ACC[e] = std::log(std::fabs(A[e] * B[e])); break;
      case 11: for (int e = 0; e < E; e++) { const R m = A[e] * B[e]; ACC[e] = m * m; } break;
      case 12: for (int e = 0; e < E; e++) { const R m = A[e] * B[e]; ACC[e] = m * m * m; } break;
      case 13: for (int e = 0; e < E; e++) ACC[e] = std::sqrt(std::fabs(A[e] * B[e])); break;
      case 14: for (int e = 0; e < E; e++) { const R m = A[e] * B[e]; ACC[e] = std::copysign(std::pow(std::fabs(m), (R)1 / (R)3), m); } break;
      case 15: for (int e = 0; e < E; e++) ACC[e] = -(A[e] * B[e]); break;
      case 16: case 17: case 18: case 19: {
        const uint32_t cw = ((x >> 8) & 0x7FFFu) | (((x >> 23) & 3u) << 29);
        const R* Cc = W.fetch(cw, W.C.data());
        if (op == 16u) { for (int e = 0; e < E; e++) ACC[e] = (A[e] * B[e]) * Cc[e]; }
        else { const int k = (int)op - 16; const R ky = (R)ctab[4 * k + 1], kz = (R)ctab[4 * k + 2];
               for (int e = 0; e < E; e++) ACC[e] = ky * std::fma(kz, Cc[e], A[e] * B[e]); }
        break;
      }
      default: std::memcpy(ACC, A, sizeof(R) * E); break;
    }
  }
  W.accp = ACC;
  if (y != 0xFFFFFFFFu) { R* N = W.NS.data() + y; std::memcpy(N, ACC, sizeof(R) * E); }
}

template<typename R> inline void gs_start_ref(int Q, R* cs, R* y, R tau, R& rss_out, bool& bad) {
  bad = false;
  R on[16];
  for (int k = 0; k < Q; k++) for (int i = 0; i < 48; i++) if (!std::isfinite(cs[k * 48 + i])) bad = true;
  for (int k = 1; k < Q; k++) { R n = 0; for (int i = 0; i < 48; i++) n += cs[k*48+i]*cs[k*48+i]; on[k] = std::sqrt(n); }
  for (int k = 0; k < Q; k++) {
    R* ck = cs + k * 48;
    R nn = 0; for (int i = 0; i < 48; i++) nn += ck[i] * ck[i];
    R inv;
    if (k == 0) inv = nn > 0 ? (R)1 / std::sqrt(nn) : (R)0;
    else        inv = (std::sqrt(nn) > tau * on[k]) ? (R)1 / std::sqrt(nn) : (R)0;
    for (int i = 0; i < 48; i++) ck[i] *= inv;
    for (int j = k + 1; j < Q; j++) { R* cj = cs + j * 48; R r = 0; for (int i = 0; i < 48; i++) r += ck[i] * cj[i];
      for (int i = 0; i < 48; i++) cj[i] -= r * ck[i]; }
    R r = 0; for (int i = 0; i < 48; i++) r += ck[i] * y[i];
    for (int i = 0; i < 48; i++) y[i] -= r * ck[i];
  }
  R rss = 0; for (int i = 0; i < 48; i++) rss += y[i] * y[i];
  rss_out = rss;
}

inline float s290_dot48(const float* x, const float* y) {
  float32x4_t a0 = vdupq_n_f32(0.f), a1 = vdupq_n_f32(0.f), a2 = vdupq_n_f32(0.f);
  for (int k = 0; k < 48; k += 12) {
    a0 = vfmaq_f32(a0, vld1q_f32(x + k), vld1q_f32(y + k));
    a1 = vfmaq_f32(a1, vld1q_f32(x + k + 4), vld1q_f32(y + k + 4));
    a2 = vfmaq_f32(a2, vld1q_f32(x + k + 8), vld1q_f32(y + k + 8));
  }
  return vaddvq_f32(vaddq_f32(vaddq_f32(a0, a1), a2));
}
inline void gs_start_f(int Q, float* cs, float* y, float tau, float& rss_out, bool& bad) {
  uint32x4_t nf = vdupq_n_u32(0);
  for (int i = 0; i < Q * 48; i += 4) nf = vorrq_u32(nf, vmvnq_u32(vcaltq_f32(vld1q_f32(cs + i), vdupq_n_f32(INFINITY))));
  bad = vmaxvq_u32(nf) != 0;
  float on[16];
  for (int k = 1; k < Q; k++) on[k] = std::sqrt(s290_dot48(cs + k * 48, cs + k * 48));
  for (int k = 0; k < Q; k++) {
    float* ck = cs + k * 48;
    const float nn = s290_dot48(ck, ck);
    float inv;
    if (k == 0) inv = nn > 0 ? 1.0f / std::sqrt(nn) : 0.0f;
    else        inv = (std::sqrt(nn) > tau * on[k]) ? 1.0f / std::sqrt(nn) : 0.0f;
    const float32x4_t vi = vdupq_n_f32(inv);
    for (int i = 0; i < 48; i += 4) vst1q_f32(ck + i, vmulq_f32(vld1q_f32(ck + i), vi));
    float r[17];
    for (int j = k + 1; j < Q; j++) r[j] = s290_dot48(ck, cs + j * 48);
    r[Q] = s290_dot48(ck, y);
    for (int j = k + 1; j <= Q; j++) {
      float* cj = (j < Q) ? cs + j * 48 : y; const float32x4_t vr = vdupq_n_f32(r[j]);
      for (int i = 0; i < 48; i += 4) vst1q_f32(cj + i, vfmsq_f32(vld1q_f32(cj + i), vr, vld1q_f32(ck + i)));
    }
  }
  rss_out = s290_dot48(y, y);
}
template<typename R> inline void gs_start(int Q, R* cs, R* y, R tau, R& rss_out, bool& bad) {
  if constexpr (S290_GS && std::is_same<R, float>::value) gs_start_f(Q, cs, y, tau, rss_out, bad);
  else gs_start_ref(Q, cs, y, tau, rss_out, bad);
}

template<typename R> inline void run_tree(Work<R>& W, const Inputs& I, uint32_t tree, R tau, R* out16) {
  const uint32_t* h = I.hdr + 4 * (size_t)tree;
  const uint32_t hx = h[0], hy = h[1], hz = h[2], hw = h[3];
  const int nth = (int)(hw & 255u), nS = (int)(hy >> 16);
  const uint32_t thoff = hw >> 8;
  const Geo& g = W.g; const int E = g.E, S = g.S, LPS = g.LPS, P = g.P, T = g.T, Q = W.Q;
  R* N = W.NS.data(); R* SV = N + (size_t)W.nslot * E; const uint32_t TT = W.TT;
  std::fill(W.NS.begin(), W.NS.end(), (R)0);
  for (uint32_t i = 0; i < 50; i++) N[TT + i] = i < 48 ? (R)i * ((R)1 / (R)47) : (i == 48 ? (R)1 : (R)2);
  for (int ix = 0; ix < S * nS; ix++) {
    const int st = ix % S, i = ix / S;
    const uint16_t w = I.sprog[hz + i]; const uint32_t op = w & 63u, a = (w >> 6) & 31u; R x = (R)1;
    if (op == 1u) {
      const R u = (nth == 1) ? (R)st / (R)15 : (R)I.sob[st * 32 + a];
      const R bf = (R)I.bred[thoff + a];
      x = (bf - (R)12.566370614359172) + ((R)2 * (R)12.566370614359172) * u;
    }
    SV[st * W.nsc + i] = x;
  }
  const uint32_t* sg = I.segs + 36 * (size_t)tree;
  std::fill(W.ACC.begin(), W.ACC.end(), (R)0); W.accp = W.ACC.data();
  for (int k = 0; k < Q; k++) {
    for (uint32_t i = sg[k]; i < sg[k + 1]; i++) {
      const uint32_t* w = I.prog + 4 * ((size_t)hx + i);
      exec_instr(W, w[0], w[1], w[2], w[3], I.ctab);
    }
    const uint32_t cw = sg[18 + k];
    R* dst = W.cap.data() + (size_t)k * E;
    const R* src = (cw >> 31) ? W.accp : W.fetch(cw, dst);
    if (src != dst) std::memcpy(dst, src, sizeof(R) * E);
  }
  for (int k = 0; k < Q; k++)
    for (int p = 0; p < P; p++)
      for (int l = 0; l < T; l++) {
        const int s = l / LPS, lq = l % LPS;
        W.cols[((size_t)k * S + s) * 48 + lq * P + p] = W.cap[(size_t)k * E + p * T + l];
      }
  R* y = W.ybuf.data();
  for (int s = 0; s < S; s++) for (int pt = 0; pt < 48; pt++) y[s * 48 + pt] = I.yv[pt];
  R cs[16 * 48];
  for (int s = 0; s < S; s++) {
    for (int k = 0; k < Q; k++) for (int pt = 0; pt < 48; pt++) cs[k * 48 + pt] = W.cols[((size_t)k * S + s) * 48 + pt];
    R rss; bool bad;
    gs_start(Q, cs, y + s * 48, tau, rss, bad);
    out16[s] = (bad || !std::isfinite(rss)) ? std::numeric_limits<R>::infinity() : rss;
  }
}

}
