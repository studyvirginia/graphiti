#pragma once
#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
#include "fit_export_h.hpp"
#include "sicis_host.hpp"
namespace geop {
inline double g_sc[2] = {1.0, 1.0};

inline double geo_nll(const fitx::Rec& fr, const float* Y0, const float* Y1, const Sketch::CisPerCellEnergyInput& en, const float* nrm = nullptr, int M = 1024, float* cout = nullptr) {
  static const bool sgn_tan = !getenv("GEO_SIGN") || atoi(getenv("GEO_SIGN")) != 0;
  const double t0 = fr.kind == 1 ? (double)fr.u[9] : 0.0, t1 = t0 + 1.0;
  static const bool gfast = !getenv("GEO_FAST") || atoi(getenv("GEO_FAST")) != 0;
  thread_local std::vector<double> cx, cy, v, tt, W, Dp, Dc, best; thread_local std::vector<uint16_t> from16; thread_local std::vector<int> from;
  cx.resize((size_t)M); cy.resize((size_t)M);
  { size_t nn = std::max(fr.s[0].n.size(), fr.s[1].n.size()) + 2; v.assign(nn, 0.0); }
  if (gfast) {
    tt.resize((size_t)M); for (int j = 0; j < M; j++) tt[(size_t)j] = t0 + (t1 - t0) * j / (M - 1);
    fitx::eval_side_vec(fr.s[0], tt.data(), M, cx.data(), W); fitx::eval_side_vec(fr.s[1], tt.data(), M, cy.data(), W);
    for (int j = 0; j < M; j++) { cx[(size_t)j] = g_sc[0] * cx[(size_t)j]; cy[(size_t)j] = g_sc[1] * cy[(size_t)j];
      if (!std::isfinite(cx[(size_t)j]) || !std::isfinite(cy[(size_t)j])) return std::numeric_limits<double>::quiet_NaN(); }
  } else
  for (int j = 0; j < M; j++) { const double t = t0 + (t1 - t0) * j / (M - 1);
    cx[(size_t)j] = g_sc[0] * fitx::eval_side(fr.s[0], t, v.data()); cy[(size_t)j] = g_sc[1] * fitx::eval_side(fr.s[1], t, v.data());
    if (!std::isfinite(cx[(size_t)j]) || !std::isfinite(cy[(size_t)j])) return std::numeric_limits<double>::quiet_NaN(); }
  const int N = 48; const double INF = 1e300;
  int jp[48];
  if (gfast && M <= 65536) {
    Dp.resize((size_t)M); Dc.resize((size_t)M); best.resize((size_t)M); from16.resize((size_t)N * M);
    for (int p = 0; p < N; p++) {
      for (int j = 0; j < M; j++) { const double dx = g_sc[0] * (double)Y0[p] - cx[(size_t)j], dy = g_sc[1] * (double)Y1[p] - cy[(size_t)j]; best[(size_t)j] = dx * dx + dy * dy; }
      if (p == 0) { for (int j = 0; j < M; j++) Dp[(size_t)j] = best[(size_t)j]; continue; }
      double pm = INF; int pj = 0; uint16_t* fr16 = from16.data() + (size_t)p * M;
      for (int j = 0; j < M; j++) { const double d = Dp[(size_t)j]; if (d < pm) { pm = d; pj = j; } fr16[j] = (uint16_t)pj; Dc[(size_t)j] = pm + best[(size_t)j]; }
      Dp.swap(Dc);
    }
    int j = 0; for (int k = 1; k < M; k++) if (Dp[(size_t)k] < Dp[(size_t)j]) j = k;
    for (int p = N - 1; p >= 0; p--) { jp[p] = j; if (p > 0) j = from16[(size_t)p * M + j]; }
  } else {
  std::vector<double> D((size_t)N * M); best.resize((size_t)M); from.resize((size_t)N * M);
  for (int p = 0; p < N; p++) {
    for (int j = 0; j < M; j++) { const double dx = g_sc[0] * (double)Y0[p] - cx[(size_t)j], dy = g_sc[1] * (double)Y1[p] - cy[(size_t)j]; best[(size_t)j] = dx * dx + dy * dy; }
    if (p == 0) { for (int j = 0; j < M; j++) D[(size_t)j] = best[(size_t)j]; continue; }
    double pm = INF; int pj = 0;
    for (int j = 0; j < M; j++) { const double d = D[(size_t)(p - 1) * M + j]; if (d < pm) { pm = d; pj = j; } from[(size_t)p * M + j] = pj; D[(size_t)p * M + j] = pm + best[(size_t)j]; }
  }
  int j = 0; for (int k = 1; k < M; k++) if (D[(size_t)(N - 1) * M + k] < D[(size_t)(N - 1) * M + j]) j = k;
  for (int p = N - 1; p >= 0; p--) { jp[p] = j; if (p > 0) j = from[(size_t)p * M + j]; }
  }
  float c[48];
  for (int p = 0; p < N; p++) {
    double bd = INF, bzx = 0, bzy = 0, bcr = 0; const double px = g_sc[0] * Y0[p], py = g_sc[1] * Y1[p];
    for (int s = std::max(jp[p] - 1, 0); s <= std::min(jp[p], M - 2); s++) {
      const double ax = cx[(size_t)s], ay = cy[(size_t)s], ex = cx[(size_t)s + 1] - ax, ey = cy[(size_t)s + 1] - ay, L2 = ex * ex + ey * ey;
      double u = L2 > 0 ? ((px - ax) * ex + (py - ay) * ey) / L2 : 0.0; u = std::min(1.0, std::max(0.0, u));
      const double zx = px - (ax + u * ex), zy = py - (ay + u * ey), d = zx * zx + zy * zy; if (d < bd) { bd = d; bzx = zx; bzy = zy; bcr = ex * zy - ey * zx; } }
    const double h = std::hypot(bzx, bzy); c[p] = nrm ? (float)(bzx * nrm[p] + bzy * nrm[48 + p]) : (float)((sgn_tan ? bcr < 0 : bzx + bzy < 0) ? -h : h);
  }
  if (cout) for (int p = 0; p < N; p++) cout[p] = c[p];
  return sicis::score48(c, en);
}

inline int trig_extra_k(const fitx::Rec& fr, int M = 64) {
  const double t0 = fr.kind == 1 ? (double)fr.u[9] : 0.0, t1 = t0 + 1.0, PI = 3.14159265358979323846;
  struct TN { int side, idx; bool sin_; std::vector<double> arg; double span; };
  std::vector<TN> tn; int extra = 0;
  for (int s = 0; s < 2; s++) {
    const fitx::Side& S = fr.s[s]; if (S.drv) continue; const size_t N = S.n.size(); if (!N) continue;
    std::vector<double> v(N + 1), V((size_t)M * N);
    for (int m = 0; m < M; m++) { fitx::eval_side(S, t0 + (t1 - t0) * m / (M - 1), v.data()); for (size_t i = 0; i < N; i++) V[(size_t)m * N + i] = v[i]; }
    std::vector<char> has(N, 0);
    for (size_t i = 0; i < N; i++) { const fitx::Node& n = S.n[i]; const int op = n.op; if (op == 1 || op == 2) continue;
      bool h = has[n.a]; if (op >= 3 && op <= 6) h = h || has[n.b];
      if (op >= 7 && op <= 9) { double lo = 1e300, hi = -1e300, mx = 0; for (int m = 0; m < M; m++) { const double a = V[(size_t)m * N + n.a]; lo = std::min(lo, a); hi = std::max(hi, a); mx = std::max(mx, std::fabs(a)); }
        if (!std::isfinite(hi - lo) || !(hi - lo > 1e-9 * (1.0 + mx))) { has[i] = h; continue; }
        if (has[n.a]) extra++;
        h = true;
        if (op != 9) { TN x; x.side = s; x.idx = (int)i; x.sin_ = op == 7; x.span = hi - lo; x.arg.resize((size_t)M); for (int m = 0; m < M; m++) x.arg[(size_t)m] = V[(size_t)m * N + n.a]; tn.push_back(std::move(x)); } }
      has[i] = h; }
  }
  for (const TN& a : tn) {
    if (!(a.span < PI)) continue;
    bool paired = false;
    for (const TN& b : tn) { if (b.sin_ == a.sin_) continue;
      for (int sg = -1; sg <= 1 && !paired; sg += 2) { double lo = 1e300, hi = -1e300; for (int m = 0; m < M; m++) { const double d = a.arg[(size_t)m] - sg * b.arg[(size_t)m]; lo = std::min(lo, d); hi = std::max(hi, d); }
        if (hi - lo < PI / 4) paired = true; }
      if (paired) break; }
    if (!paired) extra++;
  }
  return extra;
}
inline bool freet_nonaffine(const fitx::Rec& fr, int M = 64) {
  if (fr.kind != 1) return false;
  const double t0 = (double)fr.u[9];
  for (int s = 0; s < 2; s++) {
    const fitx::Side& S = fr.s[s]; if (S.drv) continue; const size_t N = S.n.size(); if (!N) continue;
    std::vector<double> v(N + 1), V((size_t)M * N), tt((size_t)M);
    for (int m = 0; m < M; m++) { tt[(size_t)m] = t0 + (double)m / (M - 1); fitx::eval_side(S, tt[(size_t)m], v.data()); for (size_t i = 0; i < N; i++) V[(size_t)m * N + i] = v[i]; }
    for (size_t i = 0; i < N; i++) { const int op = S.n[i].op; if (op < 7 || op > 9) continue; const int a = S.n[i].a;
      double st = 0, sa = 0, stt = 0, sta = 0, mx = 0;
      for (int m = 0; m < M; m++) { const double x = tt[(size_t)m], y = V[(size_t)m * N + a]; if (!std::isfinite(y)) return true; st += x; sa += y; stt += x * x; sta += x * y; mx = std::max(mx, std::fabs(y)); }
      const double d = M * stt - st * st, b1 = d != 0 ? (M * sta - st * sa) / d : 0, b0 = (sa - b1 * st) / M;
      double dev = 0; for (int m = 0; m < M; m++) dev = std::max(dev, std::fabs(V[(size_t)m * N + a] - (b0 + b1 * tt[(size_t)m])));
      if (dev > 1e-6 * (1.0 + mx)) return true; }
  }
  return false;
}
inline bool aliased(const fitx::Rec& fr, int M = 48) {
  const double t0 = fr.kind == 1 ? (double)fr.u[9] : 0.0, lim = 3.14159265358979323846 * 47.0;
  for (int s = 0; s < 2; s++) {
    const fitx::Side& S = fr.s[s]; if (S.drv) continue; const size_t N = S.n.size(); if (!N) continue;
    static const bool aw = !getenv("ALIAS_WARP") || atoi(getenv("ALIAS_WARP")) != 0;
    std::vector<double> tt((size_t)M); if (aw && fr.kind == 1) fitx::warp_t(fr.u, M, tt.data()); else for (int m = 0; m < M; m++) tt[(size_t)m] = t0 + (double)m / (M - 1);
    std::vector<double> v(N + 1), V((size_t)M * N);
    for (int m = 0; m < M; m++) { fitx::eval_side(S, tt[(size_t)m], v.data()); for (size_t i = 0; i < N; i++) V[(size_t)m * N + i] = v[i]; }
    for (size_t i = 0; i < N; i++) { const int op = S.n[i].op; if (op != 7 && op != 8) continue; const int a = S.n[i].a;
      for (int m = 1; m < M; m++) { const double d = std::fabs(V[(size_t)m * N + a] - V[(size_t)(m - 1) * N + a]) * (M - 1); if (std::isfinite(d) && d >= lim) return true; } }
    static const bool an = !getenv("ALIAS_NODE") || atoi(getenv("ALIAS_NODE")) != 0;
    if (an) {
      const int R = 16; std::vector<double> W((size_t)(R + 1) * N);
      for (int m = 1; m < M; m++) {
        const double ta = tt[(size_t)m - 1], tb = tt[(size_t)m];
        for (int r = 0; r <= R; r++) { fitx::eval_side(S, ta + (tb - ta) * r / R, v.data()); for (size_t i = 0; i < N; i++) W[(size_t)r * N + i] = v[i]; }
        for (size_t i = 0; i < N; i++) {
          int ext = 0, last = 0;
          for (int r = 1; r <= R; r++) {
            const double y0 = W[(size_t)(r - 1) * N + i], y1 = W[(size_t)r * N + i];
            if (!std::isfinite(y0) || !std::isfinite(y1)) { ext = 0; break; }
            const double dy = y1 - y0; if (std::fabs(dy) <= 1e-9 * (1.0 + std::fabs(y1))) continue;
            const int sg = dy > 0 ? 1 : -1; if (last != 0 && sg != last) ext++; last = sg;
          }
          if (ext >= 2) return true;
        }
      }
    }
  }
  return false;
}
}
