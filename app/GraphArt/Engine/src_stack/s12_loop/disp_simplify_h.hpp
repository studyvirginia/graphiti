#pragma once
#include "ui_bridge_h.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace uib { namespace ds {

inline bool on() { static const int v = getenv("DISPLAY_SIMPLIFY") ? atoi(getenv("DISPLAY_SIMPLIFY")) : 1; return v != 0; }
constexpr double DEAD = 1e-4;
inline bool sincos_on() { static const int v = getenv("DISPLAY_SINCOS") ? atoi(getenv("DISPLAY_SINCOS")) : 1; return v != 0; }
inline bool abs_on() { static const int v = getenv("DISPLAY_ABS") ? atoi(getenv("DISPLAY_ABS")) : 1; return v != 0; }
inline double ph_tol() { static const double v = getenv("DISPLAY_PHASE_TOL") ? atof(getenv("DISPLAY_PHASE_TOL")) : 1e-3; return v; }

struct Stats { bool applied = false; int mask = -1; int len0[2] = {0, 0}, len1[2] = {0, 0}; double rss0[2] = {0, 0}, rss1[2] = {0, 0}, dmax[2] = {0, 0}; double span = 0; bool haveY = false; };

struct Nd { uint8_t op; double c; int a, b; };
struct Tm { double k; int p; };
struct Lin { double c0 = 0; std::vector<Tm> t; };

struct Simp {
  std::vector<Nd> v; std::vector<std::string> ks; bool phase = true;
  int mk(uint8_t op, double c, int a = -1, int b = -1) { v.push_back({op, c, a, b}); ks.emplace_back(); return (int)v.size() - 1; }
  int cst(double c) { return mk(1, c); }
  const std::string& key(int n) {
    if (ks[(size_t)n].empty()) { const Nd d = v[(size_t)n]; std::string s; char b[48];
      if (d.op == 1) { snprintf(b, sizeof b, "%a", d.c); s = b; } else if (d.op == 2) s = "x";
      else { s = "(" + std::to_string((int)d.op) + " " + key(d.a); if (d.op <= 6) s += " " + key(d.b); s += ")"; }
      ks[(size_t)n] = s; }
    return ks[(size_t)n]; }
  double ev(int n, double x) const {
    const Nd& d = v[(size_t)n];
    if (d.op == 1) return d.c; if (d.op == 2) return x;
    if (d.op <= 6) return op_eval(d.op, ev(d.a, x), ev(d.b, x), x, 0);
    return op_eval(d.op, ev(d.a, x), 0, x, 0); }
  void addt(Lin& L, double k, int p) {
    if (k == 0) return; const std::string& s = key(p);
    for (size_t i = 0; i < L.t.size(); i++) if (key(L.t[i].p) == s) { L.t[i].k += k; if (L.t[i].k == 0) L.t.erase(L.t.begin() + (long)i); return; }
    L.t.push_back({k, p}); }
  static Lin scale(Lin L, double c) { if (c == 0) return Lin(); L.c0 *= c; for (auto& t : L.t) t.k *= c; return L; }
  Lin plus(const Lin& a, const Lin& b, double sg) { Lin r = a; r.c0 += sg * b.c0; for (auto& t : b.t) addt(r, sg * t.k, t.p); return r; }
  static bool single(const Lin& L) { return L.t.size() == 1 && L.c0 == 0; }
  Lin mono(double k, int p) { Lin r; addt(r, k, p); return r; }
  void sincos(Lin& L) {
    if (!sincos_on()) return;
    bool again = true;
    while (again) { again = false;
      for (size_t i = 0; i < L.t.size() && !again; i++) {
        const int pi0 = L.t[i].p, oi = v[(size_t)pi0].op; if (oi != 7 && oi != 8) continue;
        const int ai = v[(size_t)pi0].a; if (v[(size_t)ai].op == 1) continue;
        const std::string ka = key(ai);
        size_t j = L.t.size();
        for (size_t q = 0; q < L.t.size() && j == L.t.size(); q++) { const int pq = L.t[q].p; if (q != i && v[(size_t)pq].op + oi == 15 && key(v[(size_t)pq].a) == ka) j = q; }
        if (j == L.t.size()) continue;
        const double a = oi == 7 ? L.t[i].k : L.t[j].k, b = oi == 7 ? L.t[j].k : L.t[i].k;
        const double R = std::hypot(a, b), phi = std::atan2(a, b);
        Lin A = norm(ai); if (A.t.empty()) continue;
        A.c0 = std::remainder(A.c0 - phi, 2 * M_PI);
        if (A.t[0].k < 0) A = scale(A, -1);
        const int cn = mk(8, 0, fromLin(A));
        Lin N; N.c0 = L.c0;
        for (size_t q = 0; q < L.t.size(); q++) if (q != i && q != j) addt(N, L.t[q].k, L.t[q].p);
        addt(N, R, cn); L = N; again = true; } } }
  int fromLin(const Lin& L0) {
    Lin L = L0; sincos(L);
    int acc = -1;
    for (size_t i = 0; i < L.t.size(); i++) { double k = L.t[i].k; int p = L.t[i].p;
      if (v[(size_t)p].op == 9 && (std::fabs(k) < 1e-3 || std::fabs(k) > 1e3)) {
        const double lk = std::log(std::fabs(k)); const int z = v[(size_t)p].a; p = mk(9, 0, mk(lk < 0 ? 4 : 3, 0, z, cst(std::fabs(lk)))); k = k < 0 ? -1.0 : 1.0; }
      if (acc < 0) acc = k == 1.0 ? p : mk(5, 0, cst(k), p);
      else { const double kk = std::fabs(k); const int tm = kk == 1.0 ? p : mk(5, 0, cst(kk), p); acc = mk(k < 0 ? 4 : 3, 0, acc, tm); } }
    if (acc < 0) return cst(L.c0);
    if (L.c0 != 0) acc = mk(L.c0 < 0 ? 4 : 3, 0, acc, cst(std::fabs(L.c0)));
    return acc; }
  int mulP(int p, int q) { if (key(p) == key(q)) {
      if (abs_on() && v[(size_t)p].op == 13) return mk(15, 0, v[(size_t)p].a);
      if (abs_on() && v[(size_t)p].op == 15) return mk(11, 0, v[(size_t)p].a);
      return mk(11, 0, p); } if (key(q) < key(p)) std::swap(p, q); return mk(5, 0, p, q); }
  Lin trig(uint8_t op, Lin L) {
    double sg = 1;
    { const double h = M_PI / 2, m = std::round(L.c0 / h), e = L.c0 - m * h, tol = phase ? ph_tol() : 0;
      if ((L.c0 != 0 || m != 0) && std::fabs(e) <= tol && std::isfinite(m) && std::fabs(m) < 1e9) {
        const int mi = (int)((((long)m) % 4 + 4) % 4); L.c0 = 0;
        static const int so[2][4] = {{7, 8, 7, 8}, {8, 7, 8, 7}}; static const int ss[2][4] = {{1, 1, -1, -1}, {1, -1, -1, 1}};
        sg = ss[op == 8][mi]; op = (uint8_t)so[op == 8][mi]; } }
    if (L.t[0].k < 0) { L = scale(L, -1); if (op == 7) sg = -sg; }
    return mono(sg, mk(op, 0, fromLin(L))); }
  Lin norm(int n) {
    const Nd d = v[(size_t)n];
    switch (d.op) {
      case 1: { Lin r; r.c0 = d.c; return r; }
      case 2: return mono(1, mk(2, 0));
      case 3: case 4: { Lin a = norm(d.a), b = norm(d.b); return plus(a, b, d.op == 3 ? 1.0 : -1.0); }
      case 5: { Lin a = norm(d.a), b = norm(d.b);
        if (a.t.empty()) return scale(b, a.c0); if (b.t.empty()) return scale(a, b.c0);
        if (single(a) && single(b)) return mono(a.t[0].k * b.t[0].k, mulP(a.t[0].p, b.t[0].p));
        if (single(a)) return mono(a.t[0].k, mk(5, 0, a.t[0].p, fromLin(b)));
        if (single(b)) return mono(b.t[0].k, mk(5, 0, fromLin(a), b.t[0].p));
        return mono(1, mk(5, 0, fromLin(a), fromLin(b))); }
      case 6: { Lin a = norm(d.a), b = norm(d.b);
        if (b.t.empty()) { if (b.c0 != 0 && std::isfinite(1.0 / b.c0)) return scale(a, 1.0 / b.c0); return mono(1, mk(6, 0, fromLin(a), cst(b.c0))); }
        if (a.t.empty() && a.c0 == 0) return Lin();
        double k = 1; int num, den;
        if (single(b)) { k /= b.t[0].k; den = b.t[0].p; } else den = fromLin(b);
        if (a.t.empty()) { k *= a.c0; num = cst(1); } else if (single(a)) { k *= a.t[0].k; num = a.t[0].p; } else num = fromLin(a);
        return mono(k, mk(6, 0, num, den)); }
      default: break;
    }
    Lin a = norm(d.a); const uint8_t o = d.op;
    if (a.t.empty()) { const double r = op_eval(o, a.c0, 0, 0, 0); if (std::isfinite(r)) { Lin c; c.c0 = r; return c; } return mono(1, mk(o, 0, cst(a.c0))); }
    if (o == 7 || o == 8) return trig(o, a);
    if (abs_on() && o == 9 && a.t.size() == 1 && a.t[0].k == 1.0 && v[(size_t)a.t[0].p].op == 10) {
      const double ec = std::exp(a.c0); if (std::isfinite(ec) && ec > 1e-300) return mono(ec, mk(15, 0, v[(size_t)a.t[0].p].a)); }
    if (abs_on() && o == 10 && a.t.size() == 1 && a.c0 == 0 && v[(size_t)a.t[0].p].op == 9) {
      Lin r = norm(v[(size_t)a.t[0].p].a); r.c0 += std::log(std::fabs(a.t[0].k)); if (std::isfinite(r.c0)) return r; }
    if (abs_on() && o == 15) {
      if (single(a)) { const int p = a.t[0].p, po = v[(size_t)p].op; return mono(std::fabs(a.t[0].k), (po == 15 || po == 11 || po == 9) ? p : mk(15, 0, p)); }
      if (a.t[0].k < 0) a = scale(a, -1); return mono(1, mk(15, 0, fromLin(a))); }
    if (o == 9) {
      const double ec = std::exp(a.c0); Lin rest = a; rest.c0 = 0;
      if (std::isfinite(ec) && ec > 1e-300) return mono(ec, mk(9, 0, fromLin(rest)));
      return mono(1, mk(9, 0, fromLin(a))); }
    if (o == 10) { if (a.t[0].k < 0) a = scale(a, -1); return mono(1, mk(10, 0, fromLin(a))); }
    if (single(a)) {
      const double k = a.t[0].k; const int p = a.t[0].p;
      if (o == 11 && (v[(size_t)p].op == 7 || v[(size_t)p].op == 8)) {
        Lin u2 = scale(norm(v[(size_t)p].a), 2.0); Lin r = scale(trig(8, u2), (v[(size_t)p].op == 8 ? 0.5 : -0.5) * k * k); r.c0 += 0.5 * k * k; return r; }
      if (abs_on()) {
        const int po = v[(size_t)p].op;
        if (o == 11 && po == 13) return mono(k * k, mk(15, 0, v[(size_t)p].a));
        if (o == 13 && po == 11) return mono(std::sqrt(std::fabs(k)), mk(15, 0, v[(size_t)p].a));
        if (o == 11 && po == 15) return mono(k * k, mk(11, 0, v[(size_t)p].a));
        if (o == 12 && po == 14) { Lin r = scale(norm(v[(size_t)p].a), k * k * k); return r; }
        if (o == 14 && po == 12) { Lin r = scale(norm(v[(size_t)p].a), std::cbrt(k)); return r; } }
      if (o == 11) return mono(k * k, mk(11, 0, p));
      if (o == 12) return mono(k * k * k, mk(12, 0, p));
      if (o == 13) return mono(std::sqrt(std::fabs(k)), mk(13, 0, p));
      if (o == 14) return mono(std::cbrt(k), mk(14, 0, p)); }
    double sg = 1;
    if (a.t[0].k < 0) { a = scale(a, -1); if (o == 12 || o == 14) sg = -1; }
    return mono(sg, mk(o, 0, fromLin(a)));
  }
  int parse(const Tree& T) {
    std::vector<int> st;
    for (int k = 0; k < T.len; k++) { const uint8_t o = T.op[k];
      if (o == 1) st.push_back(cst(T.c[k])); else if (o == 2) st.push_back(mk(2, 0));
      else if (o >= 3 && o <= 6) { if (st.size() < 2) return -1; const int r = st.back(); st.pop_back(); const int l = st.back(); st.pop_back(); st.push_back(mk(o, 0, l, r)); }
      else if (o >= 7 && o <= 15) { if (st.empty()) return -1; const int a = st.back(); st.pop_back(); st.push_back(mk(o, 0, a)); }
      else return -1; }
    return st.size() == 1 ? st[0] : -1; }
  bool emit(int n, Tree& R, int side, int& len) const {
    const Nd& d = v[(size_t)n];
    if (d.op >= 3 && d.op <= 6) { if (!emit(d.a, R, side, len) || !emit(d.b, R, side, len)) return false; }
    else if (d.op >= 7) { if (!emit(d.a, R, side, len)) return false; }
    if (len >= 255) return false;
    R.op[len] = d.op; R.c[len] = d.op == 1 ? (float)d.c : 0.f; R.cid[len] = d.op == 1 ? (int16_t)(side * MAXLEN + len) : (int16_t)-1; len++; return true; }
  static void lsq(const std::vector<std::vector<double>>& C, const double* y, int N, std::vector<double>& coef, std::vector<char>& keep) {
    const size_t m = C.size(); std::vector<std::vector<double>> Q(m); std::vector<std::vector<double>> R(m, std::vector<double>(m, 0.0)); coef.assign(m, 0.0); keep.assign(m, 0);
    for (size_t k = 0; k < m; k++) { std::vector<double> w = C[k]; double on = 0; for (int i = 0; i < N; i++) on += w[(size_t)i] * w[(size_t)i]; on = std::sqrt(on);
      for (int pass = 0; pass < 2; pass++) for (size_t j = 0; j < k; j++) if (keep[j]) { double r = 0; for (int i = 0; i < N; i++) r += Q[j][(size_t)i] * w[(size_t)i]; R[j][k] += r; for (int i = 0; i < N; i++) w[(size_t)i] -= r * Q[j][(size_t)i]; }
      double nn = 0; for (int i = 0; i < N; i++) nn += w[(size_t)i] * w[(size_t)i]; nn = std::sqrt(nn);
      if (!(nn > 1e-5 * on) || !std::isfinite(nn)) continue;
      for (int i = 0; i < N; i++) w[(size_t)i] /= nn; Q[k] = w; R[k][k] = nn; keep[k] = 1; }
    std::vector<double> qy(m, 0.0); for (size_t k = 0; k < m; k++) if (keep[k]) for (int i = 0; i < N; i++) qy[k] += Q[k][(size_t)i] * y[i];
    for (size_t k = m; k-- > 0;) if (keep[k]) { double s = qy[k]; for (size_t j = k + 1; j < m; j++) if (keep[j]) s -= R[k][j] * coef[j]; coef[k] = s / R[k][k]; } }
  void dead(Lin& L, const double* x48, const double* Y, double span) {
    const int N = 48; std::vector<std::vector<double>> C(L.t.size(), std::vector<double>(N)); std::vector<char> alive(L.t.size(), 1); bool any = false;
    for (size_t j = 0; j < L.t.size(); j++) { double mx = 0; for (int i = 0; i < N; i++) { C[j][(size_t)i] = ev(L.t[j].p, x48[i]); mx = std::max(mx, std::fabs(L.t[j].k * C[j][(size_t)i])); if (!std::isfinite(C[j][(size_t)i])) return; }
      if (mx < DEAD * span) { alive[j] = 0; any = true; } }
    bool c0dead = L.c0 != 0 && std::fabs(L.c0) < DEAD * span; if (c0dead) any = true;
    if (!any) return;
    std::vector<std::vector<double>> cols; std::vector<size_t> idx; for (size_t j = 0; j < L.t.size(); j++) if (alive[j]) { cols.push_back(C[j]); idx.push_back(j); }
    const bool hasc = L.c0 != 0 && !c0dead; if (hasc) cols.emplace_back(N, 1.0);
    Lin R; std::vector<double> coef; std::vector<char> keep;
    if (!cols.empty()) lsq(cols, Y, N, coef, keep);
    if (cols.empty() || !std::all_of(coef.begin(), coef.end(), [](double c) { return std::isfinite(c); })) {
      for (size_t j = 0; j < L.t.size(); j++) if (alive[j]) R.t.push_back(L.t[j]); R.c0 = c0dead ? 0 : L.c0; L = R; return; }
    for (size_t q = 0; q < idx.size(); q++) if (keep[q]) R.t.push_back({coef[q], L.t[idx[q]].p});
    if (hasc && keep[idx.size()]) R.c0 = coef[idx.size()];
    L = R; }
};

inline bool side(const Tree& in, int sd, int mask, const double* x48, const double* Y, double span, bool isdrv, Tree& out) {
  Simp S; S.phase = (mask & 1) != 0;
  const int r = S.parse(in); if (r < 0) return false;
  Lin L = S.norm(r);
  if ((mask & 2) && Y) {
    S.dead(L, x48, Y, span);
    if (isdrv) { double mx = 0; for (int i = 0; i < 48; i++) { double a = L.c0; for (auto& t : L.t) a += t.k * S.ev(t.p, x48[i]); mx = std::max(mx, std::fabs(a - x48[i])); }
      if (mx < DEAD * span) { L = Lin(); L.t.push_back({1.0, S.mk(2, 0)}); } } }
  const int root = S.fromLin(L); int len = 0; Tree R; std::memset(&R, 0, sizeof R);
  if (!S.emit(root, R, sd, len)) return false;
  R.len = (uint8_t)len; out = R; return true;
}

inline bool simplify_pair(Tree T[2], const double* x48, const double* x120, const float* const Y[2], int drv, Stats& st, double ink[2][120]) {
  st = Stats(); double o120[2][120], o48[2][48];
  for (int s = 0; s < 2; s++) { st.len0[s] = T[s].len; for (int i = 0; i < 48; i++) o48[s][i] = eval_tree(T[s], x48[i]); for (int j = 0; j < 120; j++) o120[s][j] = eval_tree(T[s], x120[j]);
    for (int i = 0; i < 48; i++) if (!std::isfinite(o48[s][i])) return false; for (int j = 0; j < 120; j++) if (!std::isfinite(o120[s][j])) return false; }
  const bool haveY = Y && Y[0] && Y[1]; st.haveY = haveY;
  auto rng = [](const double* a, int n) { double lo = a[0], hi = a[0]; for (int i = 1; i < n; i++) { lo = std::min(lo, a[i]); hi = std::max(hi, a[i]); } return hi - lo; };
  double yv[2][48];
  for (int s = 0; s < 2; s++) for (int i = 0; i < 48; i++) yv[s][i] = haveY ? (s == drv ? x48[i] : (double)Y[s][i]) : o48[s][i];
  double span = std::max(rng(yv[0], 48), rng(yv[1], 48)); if (!(span > 0)) return false; st.span = span;
  double r0 = 0; for (int s = 0; s < 2; s++) { double q = 0; for (int i = 0; i < 48; i++) { const double d = yv[s][i] - o48[s][i]; q += d * d; } st.rss0[s] = q; r0 += q; }
  const double fl = 48 * (1.2e-7 * span) * (1.2e-7 * span);
  static const int masks[4] = {3, 1, 2, 0};
  for (int mi = 0; mi < 4; mi++) { const int m = masks[mi];
    Tree N[2]; double n120[2][120], r1 = 0, dm[2], rs[2]; bool ok = true;
    for (int s = 0; s < 2 && ok; s++) {
      if (!side(T[s], s, m, x48, haveY ? yv[s] : nullptr, span, s == drv, N[s])) ok = false;
      else if (N[s].len > T[s].len) N[s] = T[s];
      if (!ok) break;
      dm[s] = 0; for (int j = 0; j < 120; j++) { n120[s][j] = eval_tree(N[s], x120[j]); if (!std::isfinite(n120[s][j])) { ok = false; break; } dm[s] = std::max(dm[s], std::fabs(n120[s][j] - o120[s][j])); }
      if (!ok) break;
      rs[s] = 0; for (int i = 0; i < 48; i++) { const double y = eval_tree(N[s], x48[i]), d = yv[s][i] - y; rs[s] += d * d; } r1 += rs[s];
      if (!std::isfinite(r1) || dm[s] > ((m & 3) ? 1e-3 : 1e-6) * span) ok = false; }
    if (!ok || (haveY && r1 > 1.01 * r0 + fl)) continue;
    for (int s = 0; s < 2; s++) { st.len1[s] = N[s].len; st.rss1[s] = rs[s]; st.dmax[s] = dm[s] / span; T[s] = N[s]; std::memcpy(ink[s], n120[s], sizeof n120[s]); }
    st.applied = true; st.mask = m; return true; }
  return false;
}

}}
