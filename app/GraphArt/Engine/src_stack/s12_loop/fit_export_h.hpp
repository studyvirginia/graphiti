#pragma once
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <algorithm>
#include "fjn_h950.hpp"
#include "ui_bridge_h.hpp"
#include "disp_simplify_h.hpp"

namespace fitx {

struct Node { uint8_t op; uint16_t a, b; float c; };
struct Side {
  int Q = 0, nth = 0; bool drv = false;
  std::vector<Node> n; std::vector<uint16_t> root; float a[16] = {0}; float th[32] = {0};
  std::vector<std::pair<uint16_t, uint8_t>> thref;
};
struct Rec {
  uint64_t key = 0; float err = INFINITY; int kind = 0;
  Side s[2]; float u[10] = {0};
};

inline bool on() { static const int v = (getenv("FIT_EXPORT") || getenv("FINAL_PICK")) ? 1 : 0; return uib::g_ring != nullptr || v; }

inline int g_pin_t = -1; inline int g_pin_op = -1;
template <class Cand>
inline uint64_t gkey(const Cand& g) {
  uint64_t h = 0x9E3779B97F4A7C15ull;
  auto mix = [&](uint64_t x) { h ^= x + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2); h *= 0xff51afd7ed558ccdull; h ^= h >> 33; };
  for (int t = 0; t < 2; t++) {
    if (t == g_pin_t) { mix(1u | 0x100u); mix((uint64_t)g_pin_op); continue; }
    mix(g.len[t] | 0x100u);
    for (int k = 0; k < g.len[t] && k < 24; k++) { mix(g.op[t][k]); if (g.op[t][k] == 1) { uint32_t b; std::memcpy(&b, &g.c[t][k], 4); mix(b); } }
  }
  return h;
}

inline bool tiny_on() { static const bool v = getenv("TINY_FIX") && atoi(getenv("TINY_FIX")) != 0; return v; }
template <class R>
inline bool build_side(const R& T, const float* th, const float* a, Side& S) {
  S = Side(); S.Q = T.Q; S.nth = T.nth;
  if (T.Q < 1 || T.Q > 16 || T.nth > 32 || (!T.prog && (T.nprog || !tiny_on())) || !T.segs) return false;
  for (int i = 0; i < T.nth; i++) S.th[i] = th[i];
  for (int q = 0; q < T.Q; q++) S.a[q] = a[q];
  auto add = [&](uint8_t op, int x, int y, float c) -> int { if (S.n.size() >= 4000) return -1; S.n.push_back({op, (uint16_t)std::max(x, 0), (uint16_t)std::max(y, 0), c}); return (int)S.n.size() - 1; };
  int tn = -1, one = -1; std::vector<int> sn((size_t)T.nS, -1); int slot[64]; for (int i = 0; i < 64; i++) slot[i] = -1;
  int chain = -1; bool bad = false;
  auto opnd = [&](uint32_t c) -> int {
    if (c == 1022) return chain;
    if (c < 64) return slot[c];
    if (c == 900) { if (tn < 0) tn = add(2, 0, 0, 0); return tn; }
    if (c == 901) { if (one < 0) one = add(1, 0, 0, 1.0f); return one; }
    const uint32_t i = c - 64; if (i >= (uint32_t)T.nS) return -1;
    if (sn[i] < 0) { const uint16_t w = T.sp[i]; const uint32_t op = w & 63u, ai = (w >> 6) & 31u; sn[i] = add(1, 0, 0, op == 1u ? S.th[ai] : 1.0f); if (op == 1u && sn[i] >= 0) S.thref.push_back({(uint16_t)sn[i], (uint8_t)ai}); }
    return sn[i]; };
  S.root.assign((size_t)T.Q, 0);
  for (int q = 0; q < T.Q && !bad; q++) {
    for (uint32_t i = T.segs[q]; i < T.segs[q + 1]; i++) {
      const uint32_t w = T.prog[i]; const uint32_t op = w & 63, d = (w >> 6) & 63, ca = (w >> 12) & 1023, cb = w >> 22;
      const int A = opnd(ca); if (A < 0) { bad = true; break; }
      int r;
      if (op >= 3 && op <= 6) { const int B = opnd(cb); if (B < 0) { bad = true; break; } r = add((uint8_t)op, A, B, 0); }
      else if (op >= 7 && op <= 15) r = add((uint8_t)op, A, 0, 0);
      else r = A;
      if (r < 0) { bad = true; break; }
      chain = r; if (d != 63) slot[d] = r;
    }
    if (bad) break;
    const uint32_t cw = T.segs[18 + q]; const int rr = cw == 1022 ? chain : opnd(cw);
    if (rr < 0) { bad = true; break; }
    S.root[(size_t)q] = (uint16_t)rr;
  }
  return !bad;
}

inline double eval_side(const Side& S, double t, double* v) {
  if (S.drv) return t;
  const size_t N = S.n.size();
  for (size_t i = 0; i < N; i++) {
    const Node& n = S.n[i]; double r;
    switch (n.op) {
      case 1: r = n.c; break;
      case 2: r = t; break;
      case 3: r = v[n.a] + v[n.b]; break;
      case 4: r = v[n.a] - v[n.b]; break;
      case 5: r = v[n.a] * v[n.b]; break;
      case 6: r = v[n.a] / v[n.b]; break;
      case 7: r = std::sin(v[n.a]); break;
      case 8: r = std::cos(v[n.a]); break;
      case 9: r = std::exp(v[n.a]); break;
      case 10: r = std::log(std::fabs(v[n.a])); break;
      case 11: r = v[n.a] * v[n.a]; break;
      case 12: r = v[n.a] * v[n.a] * v[n.a]; break;
      case 13: r = std::sqrt(std::fabs(v[n.a])); break;
      case 14: r = std::cbrt(v[n.a]); break;
      default: r = -v[n.a]; break;
    }
    v[i] = r;
  }
  double s = 0; for (int q = 0; q < S.Q; q++) s += (double)S.a[q] * v[S.root[(size_t)q]];
  return s;
}

inline void eval_side_vec(const Side& S, const double* t, int M, double* out, std::vector<double>& W) {
  if (S.drv) { for (int m = 0; m < M; m++) out[m] = t[m]; return; }
  const size_t N = S.n.size(); W.resize(N * (size_t)M + 1);
  for (size_t i = 0; i < N; i++) {
    const Node& n = S.n[i]; double* r = W.data() + i * (size_t)M; const double* A = W.data() + (size_t)n.a * M; const double* B = W.data() + (size_t)n.b * M;
    switch (n.op) {
      case 1: for (int m = 0; m < M; m++) r[m] = n.c; break;
      case 2: for (int m = 0; m < M; m++) r[m] = t[m]; break;
      case 3: for (int m = 0; m < M; m++) r[m] = A[m] + B[m]; break;
      case 4: for (int m = 0; m < M; m++) r[m] = A[m] - B[m]; break;
      case 5: for (int m = 0; m < M; m++) r[m] = A[m] * B[m]; break;
      case 6: for (int m = 0; m < M; m++) r[m] = A[m] / B[m]; break;
      case 7: for (int m = 0; m < M; m++) r[m] = std::sin(A[m]); break;
      case 8: for (int m = 0; m < M; m++) r[m] = std::cos(A[m]); break;
      case 9: for (int m = 0; m < M; m++) r[m] = std::exp(A[m]); break;
      case 10: for (int m = 0; m < M; m++) r[m] = std::log(std::fabs(A[m])); break;
      case 11: for (int m = 0; m < M; m++) r[m] = A[m] * A[m]; break;
      case 12: for (int m = 0; m < M; m++) r[m] = A[m] * A[m] * A[m]; break;
      case 13: for (int m = 0; m < M; m++) r[m] = std::sqrt(std::fabs(A[m])); break;
      case 14: for (int m = 0; m < M; m++) r[m] = std::cbrt(A[m]); break;
      default: for (int m = 0; m < M; m++) r[m] = -A[m]; break;
    }
  }
  for (int m = 0; m < M; m++) { double s = 0; for (int q = 0; q < S.Q; q++) s += (double)S.a[q] * W[(size_t)S.root[(size_t)q] * M + m]; out[m] = s; }
}

inline void ls_fit(Side& S, const double* t, const float* Y, double tau = 1e-5) {
  const int Q = S.Q, N = 48; std::vector<double> v(S.n.size() + 1); std::vector<std::vector<double>> E((size_t)Q, std::vector<double>(N)), T((size_t)Q, std::vector<double>((size_t)Q, 0.0));
  Side u = S; for (int q = 0; q < Q; q++) { for (int j = 0; j < 16; j++) u.a[j] = j == q ? 1.f : 0.f; for (int i = 0; i < N; i++) E[(size_t)q][(size_t)i] = eval_side(u, t[i], v.data()); T[(size_t)q][(size_t)q] = 1.0; }
  std::vector<double> y(Y, Y + N), ty((size_t)Q, 0.0); std::vector<char> keep((size_t)Q, 1);
  for (int k = 0; k < Q; k++) {
    double on = 0; for (int i = 0; i < N; i++) on += E[(size_t)k][(size_t)i] * E[(size_t)k][(size_t)i]; on = std::sqrt(on);
    for (int pass = 0; pass < 1; pass++) {}
    double nn = 0; for (int i = 0; i < N; i++) nn += E[(size_t)k][(size_t)i] * E[(size_t)k][(size_t)i]; const double nrm = std::sqrt(nn);
    if (!(nrm > tau * on) || !std::isfinite(nrm)) { keep[(size_t)k] = 0; continue; }
    for (int i = 0; i < N; i++) E[(size_t)k][(size_t)i] /= nrm; for (int c = 0; c < Q; c++) T[(size_t)k][(size_t)c] /= nrm;
    for (int j = k + 1; j < Q; j++) { double r = 0; for (int i = 0; i < N; i++) r += E[(size_t)k][(size_t)i] * E[(size_t)j][(size_t)i];
      for (int i = 0; i < N; i++) E[(size_t)j][(size_t)i] -= r * E[(size_t)k][(size_t)i]; for (int c = 0; c < Q; c++) T[(size_t)j][(size_t)c] -= r * T[(size_t)k][(size_t)c]; }
    double r = 0; for (int i = 0; i < N; i++) r += E[(size_t)k][(size_t)i] * y[(size_t)i];
    for (int i = 0; i < N; i++) y[(size_t)i] -= r * E[(size_t)k][(size_t)i]; for (int c = 0; c < Q; c++) ty[(size_t)c] -= r * T[(size_t)k][(size_t)c];
  }
  for (int c = 0; c < Q; c++) S.a[c] = (float)(-ty[(size_t)c]);
}

inline void warp_t(const float* u, int n, double* t) {
  static std::vector<double> kn; if (kn.empty()) { const int deg = 3, M = 10, ni = M - deg - 1; kn.assign(deg + 1, 0.0); for (int i = 1; i <= ni; i++) kn.push_back((double)i / (ni + 1)); for (int i = 0; i <= deg; i++) kn.push_back(1.0); }
  double U = 0; for (int k = 0; k < 9; k++) U += (double)u[k] * u[k];
  for (int j = 0; j < n; j++) {
    double B[10]; fjn::bspl_basis(kn, 3, 10, n > 1 ? (double)j / (n - 1) : 0.0, B);
    double s = 0; for (int k = 0; k < 9; k++) { double Sk = 0; for (int m = k + 1; m < 10; m++) Sk += B[m]; s += Sk * (double)u[k] * u[k]; }
    t[j] = (double)u[9] + (U > 0 ? s / U : 0.0);
  }
}

inline bool emit_tree(const Rec& r, int side, uib::Tree& T) {
  const Side& S = r.s[side]; int len = 0; bool ovf = false;
  auto put = [&](uint8_t op, float c) { if (len >= uib::MAXLEN) { ovf = true; return; } T.op[len] = op; T.c[len] = c; T.cid[len] = op == 1 ? (int16_t)(side * uib::MAXLEN + len) : (int16_t)-1; len++; };
  if (S.drv) { put(2, 0); T.len = (uint8_t)len; return !ovf; }
  struct Em { const Side& S; decltype(put)& put; bool& ovf;
    void go(int i) { if (ovf) return; const Node& n = S.n[(size_t)i];
      if (n.op == 1 || n.op == 2) { put(n.op, n.c); return; }
      if (n.op == 15) { go(n.a); put(1, -1.0f); put(5, 0); return; }
      go(n.a); if (n.op <= 6) go(n.b); put(n.op, 0); } } em{S, put, ovf};
  int terms = 0;
  for (int q = 0; q < S.Q && !ovf; q++) {
    const float a = S.a[q]; if (!(a != 0.0f) || !std::isfinite(a)) continue;
    const Node& rn = S.n[S.root[(size_t)q]];
    if (rn.op == 1) put(1, a * rn.c); else { put(1, a); em.go(S.root[(size_t)q]); put(5, 0); }
    if (terms++) put(3, 0);
  }
  if (!terms) put(1, 0.0f);
  T.len = (uint8_t)std::min(len, uib::MAXLEN);
  return !ovf;
}

inline float g_target[2][48];
inline const float* g_yv1 = nullptr;
inline bool display_pair(const Rec& r, uib::Tree T[2], uib::ds::Stats* st, double ink[2][120], bool simp = true) {
  for (int s = 0; s < 2; s++) { std::memset(&T[s], 0, sizeof T[s]); if (!emit_tree(r, s, T[s])) return false; if (uib::g_xf_on && !uib::to_xframe(T[s], s)) return false; }
  uib::ds::Stats lst; uib::ds::Stats& S = st ? *st : lst; S = uib::ds::Stats();
  if (!simp || !uib::ds::on()) return true;
  double x48[48], x120[120];
  if (r.kind == 1) { warp_t(r.u, 48, x48); warp_t(r.u, 120, x120); }
  else { for (int i = 0; i < 48; i++) { const double s = i / 47.0; x48[i] = uib::g_xf_on ? (s - uib::g_xf_b) / uib::g_xf_a : s; }
         for (int j = 0; j < 120; j++) { const double s = j / 119.0; x120[j] = uib::g_xf_on ? (s - uib::g_xf_b) / uib::g_xf_a : s; } }
  const float* Yp[2] = {nullptr, nullptr};
  if (r.kind == 1 && g_yv1) { Yp[0] = g_yv1; Yp[1] = g_yv1 + 48; } else if (r.kind == 2) { Yp[0] = g_target[0]; Yp[1] = g_target[1]; }
  double dummy[2][120]; uib::ds::simplify_pair(T, x48, x120, Yp, uib::g_xf_on ? uib::g_xf_drv : -1, S, ink ? ink : dummy);
  return true;
}

template <class Cand>
inline void fill_pair(const Rec& r, uib::Pair& P, const Cand& g) {
  uib::fill_from_cand(g, P);
  bool ok = true, disp = false, shown = false; double ink[2][120];
  if (uib::ds::on()) {
    uib::Tree D[2]; uib::ds::Stats st;
    if (display_pair(r, D, &st, ink)) { P.x = D[0]; P.y = D[1]; disp = true; shown = st.applied; }
  }
  if (!disp) for (int s = 0; s < 2; s++) { uib::Tree& T = s ? P.y : P.x; uib::Tree tmp; std::memset(&tmp, 0, sizeof tmp); if (emit_tree(r, s, tmp)) T = tmp; else ok = false; }
  uib::Fit& F = P.fit; std::memset(&F, 0, sizeof F);
  F.kind = (uint8_t)r.kind; F.ovf = ok ? 0 : 1; F.pad[0] = disp ? 1 : 0;
  for (int s = 0; s < 2; s++) { F.q[s] = (uint8_t)r.s[s].Q; F.nth[s] = (uint8_t)r.s[s].nth; F.drv[s] = r.s[s].drv ? 1 : 0;
    std::memcpy(F.th[s], r.s[s].th, sizeof F.th[s]); std::memcpy(F.a[s], r.s[s].a, sizeof F.a[s]); }
  std::memcpy(F.u, r.u, sizeof F.u);
  double tt[uib::PTS]; if (r.kind == 1) warp_t(r.u, uib::PTS, tt); else for (int j = 0; j < uib::PTS; j++) tt[j] = (double)j / (uib::PTS - 1);
  std::vector<double> v; v.resize(std::max(r.s[0].n.size(), r.s[1].n.size()) + 1);
  for (int j = 0; j < uib::PTS; j++) {
    if (shown) { P.xs[j] = (float)(uib::g_map[0] * ink[0][j] + uib::g_map[1]); P.ys[j] = (float)(uib::g_map[2] * ink[1][j] + uib::g_map[3]); continue; }
    P.xs[j] = (float)(uib::g_smap[0] * eval_side(r.s[0], tt[j], v.data()) + uib::g_smap[1]);
    P.ys[j] = (float)(uib::g_smap[2] * eval_side(r.s[1], tt[j], v.data()) + uib::g_smap[3]);
  }
  P.fitted = (uint8_t)r.kind;
}

inline float g_gate = INFINITY;
inline int topn(int ft = 0) { static const int v = getenv("FIT_TOPN") ? atoi(getenv("FIT_TOPN")) : 1024; static const int w = getenv("FIT_TOPN_FT") ? atoi(getenv("FIT_TOPN_FT")) : 0; return ft ? w : v; }
template <class F> inline float batch_gate(size_t n, F errat, int ft = 0) {
  std::vector<float> e; e.reserve(n); for (size_t k = 0; k < n; k++) { const float v = errat(k); if (std::isfinite(v)) e.push_back(v); }
  float g = g_gate; const int tn = topn(ft); if (tn <= 0) return g; const size_t N = (size_t)tn;
  if (e.size() > N) { std::nth_element(e.begin(), e.begin() + (long)(N - 1), e.end()); g = std::min(g, e[N - 1]); }
  return g;
}
inline double gate_margin() { static const double v = getenv("FIT_GATE") ? atof(getenv("FIT_GATE")) : 80.0; return v; }
inline double g_ms_proto = 0;
inline double g_ms = 0, g_ms_wb = 0, g_ms_commit = 0, g_ms_retain = 0;
struct Cap { uint64_t key = 0; float err = INFINITY; int kind = 0; std::shared_ptr<const Side> p[2]; size_t off = 0; float u[10] = {0}; };
struct Store {
  std::mutex mx; std::vector<float> ar; std::vector<Cap> cur; std::unordered_map<uint64_t, Cap> all; std::vector<uint8_t> got;
  std::unordered_map<uint64_t, std::weak_ptr<const Side>> intern; std::unordered_map<uint64_t, std::shared_ptr<const Side>> pc;
  long ncap = 0, nlook = 0, nhit = 0, nbuildfail = 0;
};
inline Store& st() { static Store s; return s; }
template <class R> inline uint64_t shash(const R& T) {
  uint64_t h = 0x84222325cbf29ce4ull;
  auto m = [&](uint64_t x) { h ^= x + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2); h *= 0xff51afd7ed558ccdull; h ^= h >> 31; };
  m((uint64_t)T.Q | ((uint64_t)T.nth << 8) | ((uint64_t)T.nS << 16));
  for (size_t i = 0; i < T.nprog; i++) m(T.prog[i]);
  for (size_t i = 0; i < T.nsp; i++) m(T.sp[i]);
  for (int i = 0; i < 36; i++) m(T.segs[i]);
  return h;
}
inline std::shared_ptr<const Side> drv_proto() {
  static const std::shared_ptr<const Side> d = [] { auto s = std::make_shared<Side>(); s->drv = true; s->Q = 1; return std::shared_ptr<const Side>(s); }();
  return d;
}
template <class R> inline std::shared_ptr<const Side> proto_of(const R& T, uint64_t h) {
  Store& S = st(); std::shared_ptr<const Side> r;
  auto ii = S.intern.find(h); if (ii != S.intern.end()) r = ii->second.lock();
  if (!r) { float z[32] = {0}; auto sd = std::make_shared<Side>(); if (build_side(T, z, z, *sd)) { r = sd; S.intern[h] = r; } else S.nbuildfail++; }
  return r;
}
inline void begin(int n) { Store& S = st(); std::lock_guard<std::mutex> lk(S.mx); S.got.assign((size_t)n, 0); S.pc.clear(); S.cur.clear(); }
inline bool was_captured(int c) { Store& S = st(); return c >= 0 && (size_t)c < S.got.size() && S.got[(size_t)c]; }
inline bool capture_p(int child, uint64_t key, float err, int kind, const std::shared_ptr<const Side>& p0, const std::shared_ptr<const Side>& p1, const float* thx, const float* ax, const float* thy, const float* ay, const float* u) {
  Store& S = st(); Cap c; c.key = key; c.err = err; c.kind = kind;
  c.p[0] = p0; c.p[1] = p1;
  if (!c.p[0] || !c.p[1]) return false;
  const float* th[2] = {thx, thy}; const float* a[2] = {ax, ay};
  size_t n = 0; for (int s = 0; s < 2; s++) if (!c.p[s]->drv) n += (size_t)c.p[s]->nth + (size_t)c.p[s]->Q;
  c.off = S.ar.size(); S.ar.resize(c.off + n); { float* d = S.ar.data() + c.off;
    for (int s = 0; s < 2; s++) if (!c.p[s]->drv) { std::memcpy(d, th[s], sizeof(float) * (size_t)c.p[s]->nth); d += c.p[s]->nth; std::memcpy(d, a[s], sizeof(float) * (size_t)c.p[s]->Q); d += c.p[s]->Q; } }
  if (u) std::memcpy(c.u, u, sizeof c.u);
  S.cur.push_back(std::move(c)); S.ncap++; if (child >= 0 && (size_t)child < S.got.size()) S.got[(size_t)child] = 1;
  return true;
}
template <class R>
inline bool capture(int child, uint64_t key, float err, int kind, const R* tx, uint64_t hx, const float* thx, const float* ax, const R* ty, uint64_t hy, const float* thy, const float* ay, const float* u) {
  return capture_p(child, key, err, kind, tx ? proto_of(*tx, hx) : drv_proto(), ty ? proto_of(*ty, hy) : drv_proto(), thx, ax, thy, ay, u);
}
inline void materialize(const Cap& c, Rec& r) {
  r = Rec(); r.key = c.key; r.err = c.err; r.kind = c.kind; std::memcpy(r.u, c.u, sizeof r.u);
  size_t o = c.off; const float* pl = st().ar.data();
  for (int s = 0; s < 2; s++) { Side& d = r.s[s]; d = *c.p[s]; if (d.drv) continue;
    for (int i = 0; i < d.nth; i++) d.th[i] = pl[o + (size_t)i]; o += (size_t)d.nth;
    for (int q = 0; q < d.Q; q++) d.a[q] = pl[o + (size_t)q]; o += (size_t)d.Q;
    for (auto& tr : d.thref) d.n[tr.first].c = d.th[tr.second]; }
}
inline void commit() {
  const auto t0 = std::chrono::steady_clock::now(); struct T_ { std::chrono::steady_clock::time_point t; ~T_() { g_ms_commit += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); } } tm_{t0};
  Store& S = st(); std::lock_guard<std::mutex> lk(S.mx); if (S.all.size() + S.cur.size() + 16 > (size_t)(S.all.bucket_count() * S.all.max_load_factor())) S.all.reserve(2 * (S.all.size() + S.cur.size()) + 16);
  for (Cap& c : S.cur) { auto it = S.all.find(c.key); if (it == S.all.end()) S.all.emplace(c.key, std::move(c)); else if (c.err <= it->second.err) it->second = std::move(c); }
  S.cur.clear(); S.pc.clear();
}
inline bool g_lazy_retain = false;
inline void retain(const std::vector<uint64_t>& popkeys, bool force = false) {
  if (g_lazy_retain && !force) { Store& S0 = st(); std::lock_guard<std::mutex> lk0(S0.mx); if (S0.all.size() < 2 * popkeys.size() + 4096) return; }
  const auto t0 = std::chrono::steady_clock::now(); struct T_ { std::chrono::steady_clock::time_point t; ~T_() { g_ms_retain += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); } } tm_{t0};
  Store& S = st(); std::lock_guard<std::mutex> lk(S.mx);
  std::unordered_map<uint64_t, Cap> keep; keep.reserve(popkeys.size() * 2);
  for (uint64_t k : popkeys) { auto it = S.all.find(k); if (it != S.all.end() && !keep.count(k)) keep.emplace(k, std::move(it->second)); }
  { std::vector<float> nar; nar.reserve(keep.size() * 24);
    for (auto& kv : keep) { Cap& c = kv.second; size_t n = 0; for (int s = 0; s < 2; s++) if (!c.p[s]->drv) n += (size_t)c.p[s]->nth + (size_t)c.p[s]->Q;
      const size_t no = nar.size(); nar.insert(nar.end(), S.ar.begin() + (long)c.off, S.ar.begin() + (long)(c.off + n)); c.off = no; }
    S.ar.swap(nar); }
  S.all.swap(keep);
  for (auto it = S.intern.begin(); it != S.intern.end();) { if (it->second.expired()) it = S.intern.erase(it); else ++it; }
}
inline void present_mask(const std::vector<uint64_t>& keys, size_t lo, size_t hi, uint8_t* out) { Store& S = st(); for (size_t i = lo; i < hi; i++) out[i - lo] = S.all.find(keys[i]) != S.all.end() ? 1 : 0; }
inline void reset() { Store& S = st(); std::lock_guard<std::mutex> lk(S.mx); g_gate = INFINITY; S.cur.clear(); S.all.clear(); S.ar.clear(); S.got.clear(); S.pc.clear(); S.intern.clear(); S.ncap = S.nlook = S.nhit = S.nbuildfail = 0; }
template <class Cand>
inline bool find(const Cand& g, Rec& out) {
  Store& S = st(); std::lock_guard<std::mutex> lk(S.mx); S.nlook++;
  auto it = S.all.find(gkey(g)); if (it == S.all.end()) return false; materialize(it->second, out); S.nhit++;
  if (out.kind == 2) {
    double tg[48]; for (int i = 0; i < 48; i++) tg[i] = (double)i / 47.0;
    for (int sd = 0; sd < 2; sd++) if (!out.s[sd].drv) ls_fit(out.s[sd], tg, g_target[sd]);
  }
  return true;
}
template <class Cand>
inline bool qnth(const Cand& g, int* q, int* nth, int need_kind = 0) {
  Store& S = st(); std::lock_guard<std::mutex> lk(S.mx);
  auto it = S.all.find(gkey(g)); if (it == S.all.end()) return false;
  if (need_kind && it->second.kind != need_kind) return false;
  for (int s = 0; s < 2; s++) { q[s] = it->second.p[s]->Q; nth[s] = it->second.p[s]->nth; }
  return true;
}

struct GeoPack { std::vector<uint32_t> proto, mrec, moff; std::vector<uint8_t> st; long npacked = 0, nunk = 0; size_t maxnodes = 0; };
inline void pack_geo(const std::vector<uint64_t>& keys, int nodecap, GeoPack& P) {
  Store& S = st(); std::lock_guard<std::mutex> lk(S.mx);
  const size_t nm = keys.size(); P.proto.clear(); P.mrec.clear(); P.moff.assign(nm, 0u); P.st.assign(nm, 0); P.npacked = P.nunk = 0; P.maxnodes = 0;
  P.proto.push_back(0u);
  std::unordered_map<const Side*, uint32_t> pm; pm.reserve(4096);
  auto proto_of_side = [&](const Side* sd) -> uint32_t {
    auto it = pm.find(sd); if (it != pm.end()) return it->second;
    const uint32_t off = (uint32_t)P.proto.size(); pm[sd] = off;
    const uint32_t n = sd->drv ? 0u : (uint32_t)sd->n.size();
    P.proto.push_back(n); P.proto.push_back((uint32_t)sd->Q); P.proto.push_back(sd->drv ? 1u : 0u); P.proto.push_back((uint32_t)sd->nth); P.proto.push_back(0u);
    std::vector<uint8_t> ts(n, 0); for (const auto& tr : sd->thref) if (tr.first < n) ts[tr.first] = (uint8_t)(tr.second + 1);
    for (uint32_t i = 0; i < n; i++) { const Node& nd = sd->n[i]; uint32_t cb; std::memcpy(&cb, &nd.c, 4);
      P.proto.push_back((uint32_t)nd.op | ((uint32_t)ts[i] << 8) | ((uint32_t)nd.a << 16)); P.proto.push_back((uint32_t)nd.b); P.proto.push_back(cb); }
    if (!sd->drv) for (int q = 0; q < sd->Q; q++) P.proto.push_back((uint32_t)sd->root[(size_t)q]);
    return off; };
  for (size_t r = 0; r < nm; r++) {
    auto it = S.all.find(keys[r]); if (it == S.all.end() || it->second.kind != 1) continue;
    const Cap& c = it->second; bool okc = true;
    for (int s = 0; s < 2; s++) { const Side& sd = *c.p[s]; if (sd.drv) continue; P.maxnodes = std::max(P.maxnodes, sd.n.size()); if ((int)sd.n.size() > nodecap || sd.Q < 1 || sd.Q > 16 || sd.nth > 32) okc = false; }
    if (!okc) { P.st[r] = 2; P.nunk++; continue; }
    P.st[r] = 1; P.npacked++;
    P.moff[r] = (uint32_t)P.mrec.size();
    P.mrec.push_back(proto_of_side(c.p[0].get())); P.mrec.push_back(proto_of_side(c.p[1].get()));
    uint32_t t0b; std::memcpy(&t0b, &c.u[9], 4); P.mrec.push_back(t0b);
    size_t n = 0; for (int s = 0; s < 2; s++) if (!c.p[s]->drv) n += (size_t)c.p[s]->nth + (size_t)c.p[s]->Q;
    const size_t o = P.mrec.size(); P.mrec.resize(o + n); if (n) std::memcpy(P.mrec.data() + o, S.ar.data() + c.off, n * 4);
  }
  P.mrec.push_back(0u);
}
}
