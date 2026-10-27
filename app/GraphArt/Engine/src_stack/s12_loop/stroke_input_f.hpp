#pragma once
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

namespace stroke_in {

inline bool load_raw(const char* path, std::vector<float>& xs, std::vector<float>& ys) {
  FILE* f = fopen(path, "r"); if (!f) return false;
  char line[512];
  while (fgets(line, sizeof line, f)) {
    double x, y;
    if (line[0] == '#' || line[0] == '\n') continue;
    if (sscanf(line, "%lf %lf", &x, &y) == 2) { xs.push_back((float)x); ys.push_back((float)y); }
  }
  fclose(f);
  return true;
}

inline double g_cx = 0, g_cy = 0, g_span = 1;
inline bool resample_normalize(const float* in_x, const float* in_y, int n_in, int target_n,
                               std::vector<double>& out_x, std::vector<double>& out_y) {
  if (n_in < 4 || target_n < 4) return false;
  for (int i = 0; i < n_in; ++i)
    if (!std::isfinite(in_x[i]) || !std::isfinite(in_y[i])) return false;
  std::vector<double> cum(n_in, 0.0);
  for (int i = 1; i < n_in; ++i)
    cum[i] = cum[i - 1] + std::hypot((double)in_x[i] - in_x[i - 1], (double)in_y[i] - in_y[i - 1]);
  const double total_len = cum.back();
  if (total_len <= 1e-4) return false;
  std::vector<double> rx(target_n), ry(target_n);
  int seg = 0;
  for (int i = 0; i < target_n; ++i) {
    const double target_d = ((double)i / (target_n - 1)) * total_len;
    while (seg < n_in - 2 && cum[seg + 1] < target_d) seg++;
    const double d0 = cum[seg], d1 = cum[seg + 1];
    const double t = (d1 > d0) ? (target_d - d0) / (d1 - d0) : 0.0;
    rx[i] = in_x[seg] + t * (in_x[seg + 1] - in_x[seg]);
    ry[i] = in_y[seg] + t * (in_y[seg + 1] - in_y[seg]);
  }
  const double min_x = *std::min_element(rx.begin(), rx.end()), max_x = *std::max_element(rx.begin(), rx.end());
  const double min_y = *std::min_element(ry.begin(), ry.end()), max_y = *std::max_element(ry.begin(), ry.end());
  const double cx = (min_x + max_x) * 0.5, cy = (min_y + max_y) * 0.5;
  const double span = std::max({max_x - min_x, max_y - min_y, 1e-4});
  g_cx = cx; g_cy = cy; g_span = span;
  out_x.resize(target_n); out_y.resize(target_n);
  for (int i = 0; i < target_n; ++i) {
    out_x[i] = (float)((rx[i] - cx) / span);
    out_y[i] = (float)((ry[i] - cy) / span);
  }
  return true;
}

}
