#pragma once
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

namespace fjn {
static const int NW = 9, KW = 10, NPTS = 48;

struct WarpMat { float S[NPTS * NW]; float u0[NW]; double Sd[NPTS * NW]; double u0d[NW]; };
inline void bspl_basis(const std::vector<double>& kn, int deg, int M, double x, double* B) {
  int span = deg; const int nk = (int)kn.size();
  for (int i = deg; i < nk - deg - 1; i++) if (x >= kn[i] && x < kn[i + 1]) { span = i; break; }
  if (x >= 1.0) span = nk - deg - 2;
  for (int m = 0; m < M; m++) B[m] = 0.0;
  double N[4] = {1, 0, 0, 0}, left[4], right[4];
  for (int j = 1; j <= deg; j++) {
    left[j] = x - kn[span + 1 - j]; right[j] = kn[span + j] - x; double saved = 0;
    for (int r = 0; r < j; r++) { const double den = right[r + 1] + left[j - r]; const double tmp = den != 0 ? N[r] / den : 0.0; N[r] = saved + right[r + 1] * tmp; saved = left[j - r] * tmp; }
    N[j] = saved;
  }
  for (int r = 0; r <= deg; r++) B[span - deg + r] = N[r];
}
inline void make_warp(WarpMat& W) {
  const int deg = 3, M = 10, ni = M - deg - 1;
  std::vector<double> kn(deg + 1, 0.0); for (int i = 1; i <= ni; i++) kn.push_back((double)i / (ni + 1)); for (int i = 0; i <= deg; i++) kn.push_back(1.0);
  double d0[NW], sd = 0; double gre[10];
  for (int m = 0; m < M; m++) { gre[m] = 0; for (int j = 1; j <= deg; j++) gre[m] += kn[m + j]; gre[m] /= deg; }
  for (int k = 0; k < NW; k++) { d0[k] = gre[k + 1] - gre[k]; sd += d0[k]; }
  for (int k = 0; k < NW; k++) { W.u0d[k] = std::sqrt(d0[k] / sd); W.u0[k] = (float)W.u0d[k]; }
  for (int i = 0; i < NPTS; i++) {
    double B[10]; bspl_basis(kn, deg, M, (double)i / (NPTS - 1), B);
    for (int k = 0; k < NW; k++) { double s = 0; for (int m = k + 1; m < M; m++) s += B[m]; W.Sd[i * NW + k] = s; W.S[i * NW + k] = (float)s; }
  }
}

inline int g_fuse = getenv("FJN_FUSE_HOST") ? atoi(getenv("FJN_FUSE_HOST")) : 0;
inline int geom_CD(int T) { if (g_fuse == 2) return 1 + (T + 2) / 2; const int ND = T + 1; const int SUB = ND <= 4 ? 1 : (ND + 3) / 4; const int CH = SUB == 1 ? (ND + 1) / 2 : 2; if (g_fuse && SUB > 1) return 1 + CH * SUB; return 1 + CH; }
inline int g_lay = getenv("FJN_LAY") ? atoi(getenv("FJN_LAY")) : 0;
inline int g_sloth = getenv("FJN_SLOTH_HOST") ? atoi(getenv("FJN_SLOTH_HOST")) : 0;
inline int g_slot2 = getenv("FJN_SLOT2_HOST") ? atoi(getenv("FJN_SLOT2_HOST")) : 0;
inline int g_v3ts = getenv("FJN_V3TS_HOST") ? atoi(getenv("FJN_V3TS_HOST")) : 0;
inline int g_v3col = getenv("FJN_V3COL_HOST") ? atoi(getenv("FJN_V3COL_HOST")) : 0;
inline int ovl_floats(int Q, int T) { if (g_v3col) return 0; const int NC = T + KW, NMAX = 2 * T + KW, NPKK = NMAX * (NMAX + 1) / 2; const int colf = g_sloth ? NC * 24 : NC * 48; if (g_slot2) return std::max(2 * NPKK, (colf + 3) / 4 * 4); return (2 * NPKK + colf + 3) / 4 * 4; }
inline int g_res48_ws = (getenv("FJN_RES48_HOST") || getenv("FJN_CISFUSE_HOST")) ? 288 : 0;
inline int g_v4 = getenv("FJN_V4_HOST") ? atoi(getenv("FJN_V4_HOST")) : 0;
inline int g_diet = getenv("FJN_DIET_HOST") ? atoi(getenv("FJN_DIET_HOST")) : 0;
inline int g_v5 = getenv("FJN_V5_HOST") ? atoi(getenv("FJN_V5_HOST")) : 0;
inline int g_v5qmax = getenv("FJN_V5QMAX_HOST") ? atoi(getenv("FJN_V5QMAX_HOST")) : 99;
inline int ws_floats(int Q, int T) {
  const int NPK = std::max(4, (T + 3) / 4 * 4), QP = std::max(4, (Q + 3) / 4 * 4), NC = T + KW, NMAX = 2 * T + KW, NPKK = NMAX * (NMAX + 1) / 2;
  int o = 2 * NPK + 4 * NPK + 16 + 32 + 2 * QP + 4 * QP + 48 + 2 + 4 + 32 + 2 * NMAX + 2 * NMAX + NMAX + NMAX + NPKK + (g_lay == 2 || g_lay == 3 ? 0 : 2 * NPKK) + (g_lay == 0 ? NC * 48 : 0) + 96 + 4 + (g_diet ? 0 : g_res48_ws - (((g_v4 & 1) && g_res48_ws) ? 96 : 0)) - (((g_v5 & 4) && Q <= g_v5qmax && g_res48_ws) ? 192 : 0) - ((g_v4 & 2) ? 48 : 0) - ((g_v4 & 4) ? 2 * NMAX : 0);
  return (o + 3) / 4 * 4;
}
inline int per_min(int T) { const int NC = T + KW, NMAX = 2 * T + KW; if (g_v3col) return std::max(NMAX * (NMAX + 1) / 2, 72);  return std::max(NMAX * (NMAX + 1) / 2, (NC * 48 + 1) / 2); }
inline int slotc_floats(int Q, int T) { return (g_v3ts ? 2 * 48 : 2 * geom_CD(T) * 48) + ws_floats(Q, T); }
inline int HSq(int Q) { return (Q + 1) * 48 + (Q + 1) * 16; }
inline int eb_stride(int Q, int T) { return 4 * HSq(Q) + (g_v3col ? (T + KW) * 48 : 0); }

struct Rec {
  std::vector<uint32_t> tokP, tokD;
  uint32_t hP[4], hD[4];
  uint32_t sgP[36], sgD[36];
  uint32_t perP = 0, perD = 0;
  int nth = 0, Q = 0, nS = 0;
};
template <class TreeT, class ProgFn, class SegFn>
inline bool enc_tree(const TreeT& T, uint32_t SS, ProgFn prog, SegFn segsf, std::vector<uint32_t>& p4, uint32_t* segs, uint32_t& hw, uint32_t& SBout) {
  const uint32_t NSC = (uint32_t)T.nS + 2, SB = (uint32_t)T.nslot * SS; SBout = SB;
  auto enc = [&](uint32_t c) -> uint32_t { if (c == 1022) return 1u << 31; if (c < 64) return c * SS; if (c == 900) return 1u << 29;
    if (c == 901) return (1u << 30) | (SB + (uint32_t)T.nS); return (1u << 30) | (SB + (c - 64)); };
  const uint32_t TWO = (1u << 30) | (SB + (uint32_t)T.nS + 1);
  p4.clear();
  const uint32_t* prg = prog.first; const size_t nprg = prog.second;
  for (size_t pi = 0; pi < nprg; pi++) { const uint32_t w = prg[pi]; uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22, zo, A, B;
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a); B = enc(901); } else if (op == 3) { zo = 0; A = enc(a); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b); B = enc(a); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a); B = enc(b); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a); B = enc(901); }
    else { zo = 6; A = enc(a); B = enc(901); }
    p4.push_back(zo); p4.push_back(d == 63 ? 0xFFFFFFFFu : d * SS); p4.push_back(A); p4.push_back(B); }
  const uint16_t* tsg = segsf;
  for (int q = 0; q < 18; q++) segs[q] = tsg[q];
  for (int q = 0; q < 18; q++) segs[18 + q] = enc(tsg[18 + q]);
  if (SB >= 65536) return false;
  hw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
  return true;
}
}
