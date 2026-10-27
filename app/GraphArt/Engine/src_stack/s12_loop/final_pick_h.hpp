#pragma once
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "sketch/cis_scoring.hpp"

namespace fp {

inline std::string expr_of(const uint8_t* op, const float* c, int len, int& pos) {
  if (pos < 0 || pos >= len) return "?";
  const int k = pos--;
  const uint8_t o = op[k];
  char b[64];
  if (o == 0) return "?";
  if (o == 1) { snprintf(b, sizeof b, "%.6g", (double)c[k]); return b; }
  if (o == 2) return "t";
  static const char* un[] = {"", "", "", "", "", "", "", "sin", "cos", "exp", "log", "sq", "cube", "sqrt", "cbrt", "abs"};
  static const char* bi[] = {"", "", "", "+", "-", "*", "/"};
  if (o >= 7 && o <= 15) { std::string a = expr_of(op, c, len, pos); return std::string(un[o]) + "(" + a + ")"; }
  if (o >= 3 && o <= 6) {
    std::string r = expr_of(op, c, len, pos);
    std::string l = expr_of(op, c, len, pos);
    return "(" + l + " " + bi[o] + " " + r + ")";
  }
  return "?";
}

template <class Cand>
inline std::string tree_str(const Cand& g, int t) {
  const int len = g.len[t];
  int pos = len - 1;
  return expr_of(g.op[t], g.c[t], len, pos);
}

struct Pick { int idx; int cx; double err; bool band; };

inline int pick_index(const std::vector<double>& err, const std::vector<int>& cx, int m_dof, double z, const char* reason_out[1], bool nll = false) {
  const int n = (int)err.size();
  int best = -1;
  for (int i = 0; i < n; i++) if (std::isfinite(err[i]) && (best < 0 || err[i] < err[best])) best = i;
  if (best < 0) { reason_out[0] = "none"; return -1; }
  std::vector<int> sky;
  for (int i = 0; i < n; i++) {
    if (!std::isfinite(err[i])) continue;
    bool dom = false;
    for (int j = 0; j < n && !dom; j++)
      if (j != i && std::isfinite(err[j]) && err[j] <= err[i] && cx[j] <= cx[i] && (err[j] < err[i] || cx[j] < cx[i])) dom = true;
    if (!dom) sky.push_back(i);
  }
  const double m = std::max(1.0, (double)m_dof);
  const double a = 2.0 / (9.0 * m);
  const double fac = std::max(1.0, std::pow(1.0 - a + z * std::sqrt(a), 3.0));
  const double band = nll ? err[best] + 0.5 * z * z : err[best] * fac;
  int pick = -1;
  for (int i : sky) if (err[i] <= band && (pick < 0 || cx[i] < cx[pick] || (cx[i] == cx[pick] && err[i] < err[pick]))) pick = i;
  reason_out[0] = pick >= 0 ? "simplest skyline member inside band" : "accuracy end (no member inside band)";
  return pick >= 0 ? pick : best;
}

namespace nml {

inline constexpr double kLn2 = 0.6931471805599453;
inline constexpr double kPi = 3.141592653589793;
inline constexpr int kFlagN = 1, kFlagS = 2;

struct LogRing {
  static double zero() { return -std::numeric_limits<double>::infinity(); }
  static double fromInt(long long k) { return k > 0 ? std::log((double)k) : zero(); }
  static double add(double a, double b) {
    if (a == zero()) return b;
    if (b == zero()) return a;
    const double m = std::max(a, b);
    return m + std::log1p(std::exp(std::min(a, b) - m));
  }
  static double mul(double a, double b) { return (a == zero() || b == zero()) ? zero() : a + b; }
  static double pairSame(double a) { return a == zero() ? zero() : a + add(a, 0.0) - kLn2; }
};

inline std::vector<double> canon_counts(int cmax, bool alpha8 = true) {
  using R = LogRing;
  std::vector<std::array<double, 4>> n((size_t)std::max(cmax, 1) + 1);
  for (auto& z : n) z.fill(R::zero());
  if (cmax < 1) return std::vector<double>(n.size(), R::zero());
  auto acc = [&](std::array<double, 4>& dst, int f, double v) { dst[f] = R::add(dst[f], v); };
  acc(n[1], kFlagN | kFlagS, R::fromInt(1));
  acc(n[1], kFlagN, R::fromInt(1));
  for (int c = 2; c <= cmax; ++c) {
    auto& out = n[c];
    for (int g = 0; g < 4; ++g) {
      const double v = n[c - 1][g];
      if (v == R::zero()) continue;
      const bool gN = g & kFlagN, gS = g & kFlagS;
      acc(out, (gN || gS) ? kFlagN : 0, v);
      acc(out, gS ? kFlagN : 0, v);
      for (int u = 0; u < 4; ++u) acc(out, 0, v);
      if (alpha8) { acc(out, gN ? kFlagN : 0, v); acc(out, gN ? kFlagN : 0, v); }
    }
    for (int a = 1; a <= c - 2; ++a) {
      const int b = c - 1 - a;
      const auto& A = n[a];
      const auto& B = n[b];
      for (int f = 0; f < 4; ++f) {
        if (A[f] == R::zero()) continue;
        for (int g = 0; g < 4; ++g) {
          if (B[g] == R::zero()) continue;
          const bool fN = f & kFlagN, fS = f & kFlagS, gN = g & kFlagN, gS = g & kFlagS;
          const double ordered = R::mul(A[f], B[g]);
          const int addF = ((fN && gN) ? kFlagN : 0) | ((fS || gS) ? kFlagS : 0);
          const int mulF = (fN || gN) ? kFlagN : 0;
          acc(out, mulF, ordered);
          if (!gN) acc(out, addF, ordered);
          if (a < b) { acc(out, addF, ordered); acc(out, mulF, ordered); }
          else if (a == b) {
            if (f < g) { acc(out, addF, ordered); acc(out, mulF, ordered); }
            else if (f == g) { const double pr = R::pairSame(A[f]); acc(out, addF, pr); acc(out, mulF, pr); }
          }
        }
      }
    }
  }
  std::vector<double> out(n.size(), R::zero());
  for (size_t c = 0; c < n.size(); ++c) for (double v : n[c]) out[c] = R::add(out[c], v);
  return out;
}

inline double struct_code_canon(int c) {
  static const bool alpha8 = [] { const char* v = std::getenv("STRUCT_ALPHABET"); return v == nullptr || std::strtol(v, nullptr, 0) != 0; }();
  static const std::vector<double> table = canon_counts(128, alpha8);
  if (c < 1) return std::numeric_limits<double>::infinity();
  if ((size_t)c >= table.size()) return std::numeric_limits<double>::infinity();
  int fl = 0;
  for (unsigned v = (unsigned)c; v > 1U; v >>= 1U) ++fl;
  return table[c] + (2.0 * fl + 1.0) * kLn2;
}

inline int const_leaves(const uint8_t* op, int len) {
  int k = 0;
  for (int i = 0; i < len; i++) if (op[i] == 1) k++;
  return k;
}

}

inline constexpr int kT13Default = 3;
inline int t13_mask() {
  static const int m = [] { const char* v = std::getenv("SKETCH_T13"); return v != nullptr ? (int)std::strtol(v, nullptr, 0) : kT13Default; }();
  return m;
}
inline bool t13_bit(int b) { return ((t13_mask() >> b) & 1) != 0; }

inline int pick_index_nml(const std::vector<double>& err, const std::vector<int>& cx,
                          const std::vector<int>& nx, const std::vector<int>& ny,
                          const std::vector<int>& k, int n_pts, std::vector<double>* total_out = nullptr,
                          bool all_members = false) {
  const int n = (int)err.size();
  std::vector<int> sky;
  for (int i = 0; i < n; i++) {
    if (!std::isfinite(err[i])) continue;
    bool dom = false;
    for (int j = 0; !all_members && j < n && !dom; j++)
      if (j != i && std::isfinite(err[j]) && err[j] <= err[i] && cx[j] <= cx[i] && (err[j] < err[i] || cx[j] < cx[i])) dom = true;
    if (!dom) sky.push_back(i);
  }
  if (total_out) total_out->assign(n, std::numeric_limits<double>::infinity());
  const double half_log = 0.5 * std::log((double)std::max(n_pts, 1) / (2.0 * nml::kPi));
  int pick = -1;
  double bestL = std::numeric_limits<double>::infinity();
  for (int i : sky) {
    const double L = err[i] + k[i] * half_log + nml::struct_code_canon(nx[i]) + nml::struct_code_canon(ny[i]);
    if (total_out) (*total_out)[i] = L;
    if (std::isfinite(L) && (pick < 0 || L < bestL || (L == bestL && err[i] < err[pick]))) { pick = i; bestL = L; }
  }
  return pick;
}

namespace nml {

struct Member {
  std::vector<double> theta;
  std::vector<int> roles;
  int nx = 0;
  int nodesX = 0, nodesY = 0;
  int nodesXReduced = 0, nodesYReduced = 0;
};

struct Ranges { double extent = 1.0, freq = 0.0, phase = 0.0, expRate = 0.0; bool expRateOk = false; };

using ResidualFn = std::function<bool(const std::vector<double>& theta, std::vector<double>& r48)>;

struct FullTerms {
  double nll = std::numeric_limits<double>::quiet_NaN();
  double paramTerm = std::numeric_limits<double>::quiet_NaN();
  double fisherTerm = std::numeric_limits<double>::quiet_NaN();
  double rangeTerm = std::numeric_limits<double>::quiet_NaN();
  double structTerm = std::numeric_limits<double>::quiet_NaN();
  double total = std::numeric_limits<double>::infinity();
  int k = 0, kId = 0, n = 0, nDead = 0, nCollinear = 0, nOutOfRangeId = 0;
  bool admissible = true, ok = false;
};

struct IdentLD { double logDet = 0; int rank = 0, nDead = 0, nCollinear = 0; std::vector<int> kept; bool ok = false; };
inline IdentLD identifiable_logdet(const std::vector<std::vector<double>>& cols, int n, double eps) {
  IdentLD out;
  const int k = (int)cols.size();
  if (k == 0) { out.ok = true; return out; }
  double tol = (double)std::max(n, k) * eps, maxNorm = 0;
  std::vector<double> norms(k, 0.0);
  for (int c = 0; c < k; c++) { double s = 0; for (double v : cols[c]) s += v * v; norms[c] = std::sqrt(s); maxNorm = std::max(maxNorm, norms[c]); }
  if (!(maxNorm > 0) || !std::isfinite(maxNorm)) { out.nDead = k; out.ok = true; return out; }
  std::vector<int> live;
  for (int c = 0; c < k; c++) { if (norms[c] > tol * maxNorm) live.push_back(c); else out.nDead++; }
  if (live.empty()) { out.ok = true; return out; }
  const int L = (int)live.size();
  std::vector<std::vector<double>> q(L, std::vector<double>(n));
  for (int c = 0; c < L; c++) for (int i = 0; i < n; i++) q[c][i] = cols[live[c]][i] / norms[live[c]];
  std::vector<char> used(L, 0);
  double r00 = 0;
  for (int step = 0; step < std::min(n, L); step++) {
    int piv = -1; double best = -1;
    for (int c = 0; c < L; c++) if (!used[c]) { double s = 0; for (double v : q[c]) s += v * v; if (s > best) { best = s; piv = c; } }
    if (piv < 0) break;
    const double rii = std::sqrt(best);
    if (step == 0) r00 = rii;
    if (!(rii > tol * r00) || !(rii > 0)) break;
    used[piv] = 1;
    for (int i = 0; i < n; i++) q[piv][i] /= rii;
    for (int c = 0; c < L; c++) if (!used[c]) { double d = 0; for (int i = 0; i < n; i++) d += q[piv][i] * q[c][i]; for (int i = 0; i < n; i++) q[c][i] -= d * q[piv][i]; }
    out.logDet += 2.0 * std::log(rii) + 2.0 * std::log(norms[live[piv]]);
    out.kept.push_back(live[piv]);
    out.rank++;
  }
  out.nCollinear = L - out.rank;
  std::sort(out.kept.begin(), out.kept.end());
  out.ok = true;
  return out;
}

inline bool full_terms(const Member& m, const ResidualFn& resid, const Sketch::CisPerCellEnergyInput& en,
                       const Ranges& rg, FullTerms& t, double fd_rel = 1e-3) {
  const int k = (int)m.theta.size();
  t.k = k;
  if ((int)m.roles.size() != k) return false;
  std::vector<double> r0;
  if (!resid(m.theta, r0) || r0.empty()) return false;
  Sketch::CisScoringConfig cfg;
  cfg.useCisScoring = true;
  cfg.exactLikelihood = 1;
  cfg.scoreDcCell = true;
  const Sketch::CisScoringResult cis = Sketch::evaluateCisScoring(r0, en.perCellEnergy, {}, cfg, en);
  if (!cis.ok || cis.cells.empty() || !(cis.sHat > 0) || !std::isfinite(cis.sHat)) return false;
  const int n = (int)cis.cells.size();
  t.n = n;
  std::vector<std::vector<double>> Jw(k, std::vector<double>(n));
  for (int c = 0; c < k; c++) {
    const double h = fd_rel * std::max(1.0, std::fabs(m.theta[c]));
    std::vector<double> tp = m.theta, tm = m.theta, rp, rm;
    tp[c] += h; tm[c] -= h;
    if (!resid(tp, rp) || !resid(tm, rm) || rp.size() != r0.size() || rm.size() != r0.size()) return false;
    std::vector<double> d(rp.size());
    for (size_t i = 0; i < d.size(); i++) d[i] = (rp[i] - rm[i]) / (2.0 * h);
    const std::vector<double> flat = Sketch::haarFlattenCoeffs(Sketch::haarForward(d));
    if ((int)flat.size() != n) return false;
    for (int r = 0; r < n; r++) {
      const double p = cis.cells[r].predictedEnergy;
      static const bool fw = !std::getenv("FISHER_WHITEN") || std::atoi(std::getenv("FISHER_WHITEN")) != 0;
      Jw[c][r] = flat[r] / std::sqrt(fw ? p : cis.sHat * p);
    }
  }
  const IdentLD id = identifiable_logdet(Jw, n, (double)FLT_EPSILON);
  if (!id.ok) return false;
  t.kId = id.rank; t.nDead = id.nDead; t.nCollinear = id.nCollinear;
  auto rangeOf = [&](int c) {
    switch (m.roles[c]) {
      case 1: return rg.freq;
      case 2: return rg.phase;
      case 3: return rg.expRateOk ? rg.expRate : rg.extent;
      default: return rg.extent;
    }
  };
  double sumLogRange = 0;
  std::vector<char> kept(k, 0);
  for (int c : id.kept) kept[c] = 1;
  for (int c : id.kept) sumLogRange += std::log(rangeOf(c));
  int nOut = 0;
  for (int c = 0; c < k; c++) {
    if (!kept[c]) continue;
    const double a = std::fabs(m.theta[c]);
    if (m.roles[c] == 1 && !(a <= rg.freq)) nOut++;
    if (m.roles[c] == 3 && rg.expRateOk && !(a <= 0.5 * rg.expRate)) nOut++;
  }
  t.nOutOfRangeId = nOut;
  t.admissible = (nOut == 0);
  const double n_d = (double)n, kId = (double)id.rank;
  t.nll = std::numeric_limits<double>::quiet_NaN();
  t.paramTerm = 0.5 * kId * std::log(n_d / (2.0 * kPi));
  t.fisherTerm = 0.5 * (id.logDet - kId * std::log(n_d));
  t.rangeTerm = sumLogRange;
  const int cxs = m.nodesXReduced > 0 ? m.nodesXReduced : m.nodesX;
  const int cys = m.nodesYReduced > 0 ? m.nodesYReduced : m.nodesY;
  t.structTerm = struct_code_canon(cxs) + struct_code_canon(cys);
  t.total = t.paramTerm + t.fisherTerm + t.rangeTerm + t.structTerm;
  t.ok = std::isfinite(t.total);
  if (!t.admissible) t.total = std::numeric_limits<double>::infinity();
  return true;
}

}

inline double pen_param(int k, int n_pts) { return 0.5 * k * std::log((double)std::max(n_pts, 1) / (2.0 * nml::kPi)); }
inline double pen_struct(int nx, int ny) { return nml::struct_code_canon(nx) + nml::struct_code_canon(ny); }
inline double search_pen(int k, int n_pts, int nx, int ny) { return pen_param(k, n_pts) + pen_struct(nx, ny); }
inline nml::FullTerms fallback_terms(int k, int n_pts, int nx, int ny) {
  nml::FullTerms t;
  t.k = k; t.n = n_pts;
  t.paramTerm = pen_param(k, n_pts);
  t.fisherTerm = 0.0; t.rangeTerm = 0.0;
  t.structTerm = pen_struct(nx, ny);
  t.total = t.paramTerm + t.structTerm;
  t.ok = false;
  return t;
}

inline int pick_index_nml_full(const std::vector<double>& err, const std::vector<int>& cx,
                               const std::vector<nml::FullTerms>& fulls, std::vector<double>* total_out = nullptr,
                               bool all_members = false) {
  const int n = (int)err.size();
  std::vector<int> sky;
  for (int i = 0; i < n; i++) {
    if (!std::isfinite(err[i])) continue;
    bool dom = false;
    for (int j = 0; !all_members && j < n && !dom; j++)
      if (j != i && std::isfinite(err[j]) && err[j] <= err[i] && cx[j] <= cx[i] && (err[j] < err[i] || cx[j] < cx[i])) dom = true;
    if (!dom) sky.push_back(i);
  }
  if (total_out) total_out->assign(n, std::numeric_limits<double>::infinity());
  int pick = -1;
  double bestL = std::numeric_limits<double>::infinity();
  for (int i : sky) {
    const double L = err[i] + fulls[i].total;
    if (total_out) (*total_out)[i] = L;
    if (std::isfinite(L) && (pick < 0 || L < bestL || (L == bestL && err[i] < err[pick]))) { pick = i; bestL = L; }
  }
  return pick;
}

}
