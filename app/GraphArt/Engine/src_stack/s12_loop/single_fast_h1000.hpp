#pragma once
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

struct SingleFastDecision {
  int axis = 0;
  double margin = 0;
  double eps = 0;
  double backtrack_x = 0, backtrack_y = 0;
  bool closed = false;
};

static double sf_worst_backtrack(const std::vector<double>& v, double eps) {
  if (v.size() < 3) return 0;
  const double dir = v.back() >= v.front() ? 1.0 : -1.0;
  double ext = v[0], worst = 0;
  for (size_t i = 1; i < v.size(); i++) {
    const double d = dir * (v[i] - ext);
    if (d >= 0) ext = v[i];
    else worst = std::max(worst, -d);
  }
  return worst / eps;
}

static double g_sf_a0 = 0, g_sf_a1 = 1;
static void single_fast_resample(const std::vector<double>& X, const std::vector<double>& Y, int axis, int npt,
                                 std::vector<double>& ox, std::vector<double>& oy) {
  const std::vector<double>& A = axis == 1 ? X : Y;
  const std::vector<double>& B = axis == 1 ? Y : X;
  std::vector<size_t> idx(A.size()); for (size_t i = 0; i < idx.size(); i++) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(), [&](size_t p, size_t q) { return A[p] < A[q]; });
  const double a0 = A[idx.front()], a1 = A[idx.back()];
  g_sf_a0 = a0; g_sf_a1 = a1;
  std::vector<double> ga(npt), gb(npt);
  size_t j = 0;
  for (int p = 0; p < npt; p++) {
    const double a = a0 + (a1 - a0) * p / (npt - 1);
    while (j + 2 < idx.size() && A[idx[j + 1]] < a) j++;
    const double da = A[idx[j + 1]] - A[idx[j]];
    const double w = da <= 0 ? 0.0 : std::min(1.0, std::max(0.0, (a - A[idx[j]]) / da));
    ga[p] = (double)p / (npt - 1); gb[p] = B[idx[j]] + w * (B[idx[j + 1]] - B[idx[j]]);
  }
  if (axis == 1) { ox = ga; oy = gb; } else { ox = gb; oy = ga; }
}

static SingleFastDecision single_fast_detect(const std::vector<double>& X, const std::vector<double>& Y) {
  SingleFastDecision D;
  const size_t n = X.size();
  if (n < 3 || Y.size() != n) return D;
  double x0 = X[0], x1 = X[0], y0 = Y[0], y1 = Y[0];
  for (size_t i = 0; i < n; i++) { x0 = std::min(x0, X[i]); x1 = std::max(x1, X[i]); y0 = std::min(y0, Y[i]); y1 = std::max(y1, Y[i]); }
  const double frac = getenv("SINGLE_FAST_EPS") ? atof(getenv("SINGLE_FAST_EPS")) : 0.02;
  const double bw = x1 - x0, bh = y1 - y0;
  D.eps = frac * std::max(bw, bh);
  if (D.eps <= 0) return D;
  double L = 0; for (size_t i = 1; i < n; i++) L += std::hypot(X[i] - X[i - 1], Y[i] - Y[i - 1]);
  D.closed = std::hypot(X[n - 1] - X[0], Y[n - 1] - Y[0]) <= 2.0 * D.eps || std::hypot(X[n - 1] - X[0], Y[n - 1] - Y[0]) <= L / std::sqrt((double)n);
  if (D.closed) return D;
  D.backtrack_x = sf_worst_backtrack(X, D.eps);
  D.backtrack_y = sf_worst_backtrack(Y, D.eps);
  const double mx = (1.0 - D.backtrack_x), my = (1.0 - D.backtrack_y);
  if (D.backtrack_x == 0 && D.backtrack_y == 0) { D.axis = bw >= bh ? 1 : 2; D.margin = 1.0; }
  else if (mx >= my && D.backtrack_x <= 1.0) { D.axis = 1; D.margin = mx; }
  else if (D.backtrack_y <= 1.0) { D.axis = 2; D.margin = my; }
  else { D.axis = 0; D.margin = std::max(mx, my); }
  return D;
}
