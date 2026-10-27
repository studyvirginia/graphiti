#pragma once
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <functional>
#include "../s02_nf/sk_nf_h2b.hpp"
#include "sobol_tab.h"
using namespace skn;

enum { T_T = 1, T_ONE = 2 };
static const int NSV = 32, MAXSLOT = 40, NPT = 48;
static const int GN_TH_MAX = 24;

struct Tree { int Q, nth, nS, nslot; std::vector<uint32_t> prog; std::vector<uint16_t> sp; std::vector<float> th; std::vector<uint32_t> rec; uint16_t segs[36];
  const uint32_t* progp = nullptr; const uint16_t* segsp = nullptr; uint32_t nprogp = 0; };

static bool compile_tree(const View& V, Tree& T) {
  T.Q = V.Q; T.nth = V.nth; T.nS = V.nS; T.prog.clear(); T.sp.clear(); T.th.clear();
  if (V.nS > NSV || V.nth > 32) return false;
  for (int i = 0; i < V.nS; i++) { int o = V.op(i), a = V.a(i), b = V.b(i); if (a > 31 || b > 31) return false; T.sp.push_back((uint16_t)(o | (a << 6) | (b << 11))); }
  const int n = V.n(); std::vector<int> opnd(n, -1);
  for (int i = 0; i < n; i++) { int o = V.op(i); if (i < V.nS) opnd[i] = 64 + i; else if (o == 2) opnd[i] = 900; else if (o == 16) opnd[i] = 901; }
  std::vector<int> rep(n); for (int i = 0; i < n; i++) rep[i] = i;
  { std::map<std::tuple<int,int,int>, int> hc;
    for (int i = V.nS; i < n; i++) { int o = V.op(i); if (o == 2 || o == 16) continue; int a = rep[V.a(i)], b = (o >= 3 && o <= 6) ? rep[V.b(i)] : -1;
      if ((o == 3 || o == 5) && a > b) std::swap(a, b); auto key = std::make_tuple(o, a, b); auto it = hc.find(key); if (it == hc.end()) hc[key] = i; else rep[i] = it->second; } }
  auto CA = [&](int x) { return rep[V.a(x)]; }; auto CB = [&](int x) { return rep[V.b(x)]; };
  struct I { int node, op, A, B; };
  std::vector<int> corder(V.Q); for (int q = 0; q < V.Q; q++) corder[q] = q;
  auto peak = [&](const std::vector<int>& ord) {
    std::vector<I> in2; std::vector<char> dn(n, 0);
    std::function<void(int)> df = [&](int x) { x = rep[x]; if (opnd[x] >= 0 || dn[x]) return; int o = V.op(x);
      if (o >= 3 && o <= 6) { df(CA(x)); df(CB(x)); in2.push_back({x, o, CA(x), CB(x)}); } else { df(CA(x)); in2.push_back({x, o, CA(x), -1}); } dn[x] = 1; };
    for (int q : ord) { int r = rep[V.root(q)]; size_t b0 = in2.size(); df(r); if (in2.size() == b0) in2.push_back({-1, 0, r, -1}); }
    const int mm = (int)in2.size(); std::vector<int> rd(n, 0), lr(n, -1); for (int k = 0; k < mm; k++) { rd[in2[k].A]++; lr[in2[k].A] = k; if (in2[k].B >= 0) { rd[in2[k].B]++; lr[in2[k].B] = k; } }
    int live = 0, mx = 0; std::vector<char> held(n, 0);
    for (int k = 0; k < mm; k++) { auto& z = in2[k]; int ch = 0;
      if (z.node >= 0 && rd[z.node] > 0) { bool chained = (k + 1 < mm) && ((in2[k + 1].A == z.node) + (in2[k + 1].B == z.node) == rd[z.node]); if (!chained) ch = 1; }
      if (opnd[z.A] < 0 && held[z.A] && lr[z.A] == k) { held[z.A] = 0; live--; } if (z.B >= 0 && z.B != z.A && opnd[z.B] < 0 && held[z.B] && lr[z.B] == k) { held[z.B] = 0; live--; }
      if (ch) { held[z.node] = 1; live++; mx = std::max(mx, live); } }
    return mx * 1000 + mm; };
  if (V.Q > 1) { std::vector<int> ord, rem = corder;
    while (!rem.empty()) { int bi = 0, bc = 1 << 30; for (size_t i = 0; i < rem.size(); i++) { std::vector<int> tr = ord; tr.push_back(rem[i]); for (size_t j = 0; j < rem.size(); j++) if (j != i) tr.push_back(rem[j]); int c = peak(tr); if (c < bc) { bc = c; bi = (int)i; } }
      ord.push_back(rem[bi]); rem.erase(rem.begin() + bi); }
    corder = ord; }
  std::vector<I> ins; std::vector<char> done(n, 0); std::vector<int> segend;
  std::function<void(int)> dfs = [&](int x) { x = rep[x]; if (opnd[x] >= 0 || done[x]) return; int o = V.op(x);
    if (o >= 3 && o <= 6) { dfs(CA(x)); dfs(CB(x)); ins.push_back({x, o, CA(x), CB(x)}); } else { dfs(CA(x)); ins.push_back({x, o, CA(x), -1}); } done[x] = 1; };
  for (int qi = 0; qi < V.Q; qi++) { int r = rep[V.root(corder[qi])]; size_t before = ins.size(); dfs(r);
    if (ins.size() == before) ins.push_back({-1, 0, r, -1}); segend.push_back((int)ins.size()); }
  const int m = (int)ins.size(); std::vector<int> reads(n, 0), lastread(n, -1);
  for (int k = 0; k < m; k++) { auto& z = ins[k]; reads[z.A]++; lastread[z.A] = k; if (z.B >= 0) { reads[z.B]++; lastread[z.B] = k; } }
  std::vector<int> pos(n, -1); for (int k = 0; k < m; k++) if (ins[k].node >= 0) pos[ins[k].node] = k;
  std::vector<char> chain(n, 0);
  for (int k = 0; k + 1 < m; k++) { int x = ins[k].node; if (x < 0 || reads[x] == 0) continue; auto& z = ins[k + 1]; int r1 = (z.A == x) + (z.B == x); if (r1 == reads[x]) chain[x] = 1; }
  std::vector<int> slot(n, -1), freeS; int nsl = 0;
  std::vector<uint16_t> capt(17, 1022); std::vector<int> qof(m, -1), newend(m + 1, 0); { int q = 0; for (int k = 0; k < m; k++) { qof[k] = q; if (k + 1 == segend[q]) q++; } }
  auto code = [&](int x, int k) { if (opnd[x] >= 0) return opnd[x]; if (chain[x] && pos[x] == k - 1) return 1022; return slot[x]; };
  T.prog.clear();
  for (int k = 0; k < m; k++) { auto& z = ins[k]; int A = code(z.A, k), B = z.B >= 0 ? code(z.B, k) : 0; if (A < 0 || B < 0) return false;
    if (opnd[z.A] < 0 && slot[z.A] >= 0 && lastread[z.A] == k) { freeS.push_back(slot[z.A]); slot[z.A] = -2; }
    if (z.B >= 0 && z.B != z.A && opnd[z.B] < 0 && slot[z.B] >= 0 && lastread[z.B] == k) { freeS.push_back(slot[z.B]); slot[z.B] = -2; }
    int d = 63; if (z.node >= 0 && reads[z.node] > 0 && !chain[z.node]) { if (!freeS.empty()) { d = freeS.back(); freeS.pop_back(); } else d = nsl++; slot[z.node] = d; }
    if (z.node < 0) { capt[qof[k]] = (uint16_t)A; newend[k + 1] = (int)T.prog.size(); continue; }
    T.prog.push_back((uint32_t)z.op | (uint32_t)d << 6 | (uint32_t)A << 12 | (uint32_t)B << 22); newend[k + 1] = (int)T.prog.size(); }
  for (int q = 0; q < 36; q++) T.segs[q] = 0;
  T.segs[0] = 0; for (int q = 0; q < V.Q; q++) T.segs[q + 1] = (uint16_t)newend[segend[q]];
  for (int q = 0; q < 18; q++) T.segs[18 + q] = 1022;
  for (int q = 0; q < V.Q; q++) T.segs[18 + q] = capt[q];
  T.nslot = nsl; if (nsl > MAXSLOT) return false;
  for (int s = 0; s < V.nth; s++) T.th.push_back(V.th(s)); return true;
}

static bool zconvS(const Tree& T, std::vector<uint32_t>& p4, uint32_t* segs, int SS) {
  const uint32_t SL = SS * 48; const int nsc = T.nS <= 4 ? 4 : T.nS <= 8 ? 8 : T.nS <= 16 ? 16 : 32; const uint32_t TT = T.nslot * SL + SS * nsc;
  auto enc = [&](uint32_t c) -> uint32_t { if (c == 1022) return (1u << 31) | (3u << 29) | (TT + 48); if (c < 64) return c * SL; if (c == 900) return 2u << 29; if (c == 901) return (3u << 29) | (TT + 48); return (1u << 29) | (c - 64); };
  const uint32_t TWO = (3u << 29) | (TT + 49);
  p4.clear();
  for (auto w : T.prog) { uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22, zo, A, B;
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a); B = enc(901); } else if (op == 3) { zo = 0; A = enc(a); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b); B = enc(a); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a); B = enc(b); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a); B = enc(901); }
    else { zo = 6; A = enc(a); B = enc(901); }
    p4.push_back(zo); p4.push_back(d == 63 ? 0xFFFFFFFFu : d * SL); p4.push_back(A); p4.push_back(B); }
  for (int q = 0; q < 18; q++) segs[q] = T.segs[q];
  for (int q = 0; q < 18; q++) segs[18 + q] = enc(T.segs[18 + q]);
  { const size_t n = p4.size() / 4; std::vector<char> bnd(n + 1, 0); for (int q = 0; q <= T.Q; q++) if (T.segs[q] <= n) bnd[T.segs[q]] = 1;
    std::vector<uint32_t> o; std::vector<int> newidx(n + 1, 0); size_t k = 0;
    while (k < n) { newidx[k] = (int)(o.size() / 4); const uint32_t* a = &p4[k * 4];
      if (k + 1 < n && !bnd[k + 1] && a[0] == 0 && a[1] == 0xFFFFFFFFu) { const uint32_t* b = &p4[(k + 1) * 4];
        if (b[0] >= 7 && b[0] <= 15 && (b[2] >> 31)) { o.push_back(b[0]); o.push_back(b[1]); o.push_back(a[2]); o.push_back(a[3]); newidx[k + 1] = newidx[k]; k += 2; continue; } }
      o.insert(o.end(), a, a + 4); k++; }
    newidx[n] = (int)(o.size() / 4); for (int q = 0; q <= T.Q; q++) segs[q] = newidx[T.segs[q]]; p4 = o; }
  return true; }

static bool zgconv(const Tree& T, int CD, std::vector<uint32_t>& p4, uint32_t* segs, uint32_t& hw, size_t& per) {
  const uint32_t NSC = T.nS + 2, SB = T.nslot * CD * 48;
  auto enc = [&](uint32_t c) -> uint32_t { if (c == 1022) return 1u << 31; if (c < 64) return c * CD * 48; if (c == 900) return 1u << 29; if (c == 901) return (1u << 30) | (SB + T.nS); return (1u << 30) | (SB + (c - 64)); };
  const uint32_t TWO = (1u << 30) | (SB + T.nS + 1);
  p4.clear();
  for (auto w : T.prog) { uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22, zo, A, B;
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a); B = enc(901); } else if (op == 3) { zo = 0; A = enc(a); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b); B = enc(a); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a); B = enc(b); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a); B = enc(901); }
    else { zo = 6; A = enc(a); B = enc(901); }
    p4.push_back(zo); p4.push_back(d == 63 ? 0xFFFFFFFFu : d * CD * 48); p4.push_back(A); p4.push_back(B); }
  for (int q = 0; q < 18; q++) segs[q] = T.segs[q];
  for (int q = 0; q < 18; q++) segs[18 + q] = enc(T.segs[18 + q]);
  { const size_t n = p4.size() / 4; std::vector<char> bnd(n + 1, 0); for (int q = 0; q <= T.Q; q++) if (T.segs[q] <= n) bnd[T.segs[q]] = 1;
    std::vector<uint32_t> o; std::vector<int> newidx(n + 1, 0); size_t k = 0;
    while (k < n) { newidx[k] = (int)(o.size() / 4); const uint32_t* a = &p4[k * 4];
      if (k + 1 < n && !bnd[k + 1] && a[0] == 0 && a[1] == 0xFFFFFFFFu) { const uint32_t* b = &p4[(k + 1) * 4];
        if (b[0] >= 7 && b[0] <= 15 && (b[2] >> 31)) { o.push_back(b[0]); o.push_back(b[1]); o.push_back(a[2]); o.push_back(a[3]); newidx[k + 1] = newidx[k]; k += 2; continue; } }
      o.insert(o.end(), a, a + 4); k++; }
    newidx[n] = (int)(o.size() / 4); for (int q = 0; q <= T.Q; q++) segs[q] = newidx[T.segs[q]]; p4 = o; }
  if (SB >= 65536) return false; hw = NSC | (SB << 8) | ((uint32_t)T.nth << 24); per = (size_t)SB + 2 * CD * NSC; return true; }

static double eval_rss(const View& V, const double* th, const double* y) {
  int n = V.n(), Q = V.Q; std::vector<double> v(n), R((size_t)Q * (Q + 1), 0.0); double rs = 0;
  for (int p = 0; p < NPT; p++) { double t = p / 47.0;
    for (int i = 0; i < n; i++) { int o = V.op(i); double A = (o >= 3 && o <= 15 && o != 16) ? v[V.a(i)] : 0, B = (o >= 3 && o <= 6) ? v[V.b(i)] : 0; double x = 0;
      switch (o) { case 1: x = th[V.a(i)]; break; case 2: x = t; break; case 16: x = 1; break; case 3: x = A + B; break; case 4: x = A - B; break; case 5: x = A * B; break; case 6: x = A / B; break;
        case 7: x = sin(A); break; case 8: x = cos(A); break; case 9: x = exp(A); break; case 10: x = log(fabs(A)); break; case 11: x = A * A; break; case 12: x = A * A * A; break; case 13: x = sqrt(fabs(A)); break; case 14: x = cbrt(A); break; case 15: x = -A; break; }
      v[i] = x; }
    std::vector<double> row(Q + 1); for (int q = 0; q < Q; q++) { row[q] = v[V.root(q)]; if (!std::isfinite(row[q])) return INFINITY; } row[Q] = y[p];
    for (int k = 0; k < Q; k++) { double a = R[k * (Q + 1) + k], b = row[k], h = std::hypot(a, b); double c = h > 0 ? a / h : 1, s = h > 0 ? b / h : 0; R[k * (Q + 1) + k] = h;
      for (int j = k + 1; j <= Q; j++) { double rj = R[k * (Q + 1) + j], xj = row[j]; R[k * (Q + 1) + j] = c * rj + s * xj; row[j] = c * xj - s * rj; } }
    rs += row[Q] * row[Q]; }
  return rs; }
static void start_theta(const std::vector<float>& bred, int nth, int s, double* th) { const float P = 12.566370614359172f;
  for (int j = 0; j < nth; j++) { float u = nth == 1 ? (float)s / 15.0f : SOB_U[s][j]; th[j] = (double)((bred[j] - P) + (2.0f * P) * u); } }
static bool eval_cols(const View& V, const double* th, std::vector<double>& E) {
  const int n = V.n(), Q = V.Q; std::vector<double> v(n); E.assign((size_t)Q * NPT, 0);
  for (int p = 0; p < NPT; p++) { double t = p / 47.0;
    for (int i = 0; i < n; i++) { int o = V.op(i); double A = (o >= 3 && o <= 15) ? v[V.a(i)] : 0, B = (o >= 3 && o <= 6) ? v[V.b(i)] : 0, x = 0;
      switch (o) { case 1: x = th[V.a(i)]; break; case 2: x = t; break; case 16: x = 1; break; case 3: x = A + B; break; case 4: x = A - B; break; case 5: x = A * B; break; case 6: x = A / B; break;
        case 7: x = sin(A); break; case 8: x = cos(A); break; case 9: x = exp(A); break; case 10: x = log(fabs(A)); break; case 11: x = A * A; break; case 12: x = A * A * A; break; case 13: x = sqrt(fabs(A)); break; case 14: x = cbrt(A); break; case 15: x = -A; break; }
      v[i] = x; }
    for (int q = 0; q < Q; q++) { double x = v[V.root(q)]; if (!std::isfinite(x)) return false; E[(size_t)q * NPT + p] = x; } }
  return true; }
static double vp(const std::vector<double>& E, int Q, const double* y, double tau, std::vector<double>& a, std::vector<double>& r) {
  std::vector<double> Qm(E); std::vector<double> R((size_t)Q * Q, 0), z(y, y + NPT); std::vector<char> keep(Q, 0); a.assign(Q, 0);
  for (int k = 0; k < Q; k++) { double on = 0; for (int p = 0; p < NPT; p++) on += E[(size_t)k * NPT + p] * E[(size_t)k * NPT + p]; on = sqrt(on);
    for (int i = 0; i < k; i++) if (keep[i]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)i * NPT + p] * Qm[(size_t)k * NPT + p]; R[i * Q + k] = s; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] -= s * Qm[(size_t)i * NPT + p]; }
    double nn = 0; for (int p = 0; p < NPT; p++) nn += Qm[(size_t)k * NPT + p] * Qm[(size_t)k * NPT + p]; nn = sqrt(nn);
    keep[k] = nn > tau * on; R[k * Q + k] = keep[k] ? nn : 0; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] = keep[k] ? Qm[(size_t)k * NPT + p] / nn : 0; }
  std::vector<double> qy(Q, 0); for (int k = 0; k < Q; k++) if (keep[k]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)k * NPT + p] * z[p]; qy[k] = s; for (int p = 0; p < NPT; p++) z[p] -= s * Qm[(size_t)k * NPT + p]; }
  for (int k = Q - 1; k >= 0; k--) { if (!keep[k]) continue; double s = qy[k]; for (int j = k + 1; j < Q; j++) s -= R[k * Q + j] * a[j]; a[k] = s / R[k * Q + k]; }
  r = z; double rss = 0; for (int p = 0; p < NPT; p++) rss += z[p] * z[p]; return rss; }
struct RefOut { double rss0, rssf; int steps; };
static RefOut ref_gn(const View& V, std::vector<double> th, const double* y, double tau, double floor_) {
  RefOut o{INFINITY, INFINITY, 0}; const int m = V.nth, Q = V.Q; std::vector<double> E, a, r;
  if (!eval_cols(V, th.data(), E)) return o; double rss = vp(E, Q, y, tau, a, r); o.rss0 = rss;
  for (int it = 0; it < 3; it++) {
    std::vector<double> J((size_t)NPT * m);
    for (int j = 0; j < m; j++) { double hh = 1e-6 * std::max(1.0, fabs(th[j])); std::vector<double> t2 = th; t2[j] += hh; std::vector<double> E2; if (!eval_cols(V, t2.data(), E2)) { o.rssf = rss; return o; }
      std::vector<double> g(NPT, 0); for (int q = 0; q < Q; q++) if (a[q] != 0) for (int p = 0; p < NPT; p++) g[p] += a[q] * (E2[(size_t)q * NPT + p] - E[(size_t)q * NPT + p]) / hh;
      std::vector<double> at, pr; vp(E, Q, g.data(), tau, at, pr); for (int p = 0; p < NPT; p++) J[(size_t)p * m + j] = -pr[p]; }
    std::vector<double> Jc((size_t)m * NPT); for (int j = 0; j < m; j++) for (int p = 0; p < NPT; p++) Jc[(size_t)j * NPT + p] = -J[(size_t)p * m + j];
    std::vector<double> d, rr; vp(Jc, m, r.data(), tau, d, rr); double pred = rss; for (int p = 0; p < NPT; p++) pred -= rr[p] * rr[p];
    if (pred < floor_) break;
    bool acc = false; double f = 1;
    for (int hv = 0; hv <= 4 && !acc; hv++, f *= 0.5) { std::vector<double> tn = th; for (int j = 0; j < m; j++) tn[j] += f * d[j]; std::vector<double> E2, a2, r2;
      if (!eval_cols(V, tn.data(), E2)) continue; double rs2 = vp(E2, Q, y, tau, a2, r2); if (rs2 <= rss) { th = tn; E = E2; a = a2; r = r2; rss = rs2; acc = true; } }
    if (!acc) break; o.steps++; }
  o.rssf = rss; return o; }
