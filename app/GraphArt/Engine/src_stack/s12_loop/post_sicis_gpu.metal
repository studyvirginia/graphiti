#include <metal_stdlib>
using namespace metal;

constant uint RS_NPT = 48u;
constant uint RS_MAXN = 96u;
constant uint RS_MAXQ = 16u;

static inline bool rs_side(device const uint* rec, uint o, device const float* thp, thread const float* tb,
                           device const float* ytg, thread float* z) {
  const uint p0 = rec[o], p1 = rec[o + 1u];
  const uint Q = (p0 >> 16) & 0xFFu;
  const uint nS = p1 & 0xFFFFu, nV = p1 >> 16;
  const uint n = nS + nV;
  if (n > RS_MAXN || Q > RS_MAXQ || Q == 0u) return false;
  const uint nr = (Q + 1u) / 2u;
  device const uint* roots_w = rec + o + 3u;
  device const uint* ins = roots_w + nr;
  float v[RS_MAXN];
  float E[RS_MAXQ * RS_NPT];
  for (uint p = 0u; p < RS_NPT; p++) {
    const float t = tb[p];
    for (uint i = 0u; i < n; i++) {
      const uint w = ins[i];
      const uint op = w & 0x3Fu;
      const uint a = (w >> 6) & 0x1FFFu, b = (w >> 19) & 0x1FFFu;
      const float A = (op >= 3u && op <= 15u) ? v[a] : 0.0f;
      const float B = (op >= 3u && op <= 6u) ? v[b] : 0.0f;
      float x = 0.0f;
      switch (op) {
        case 1u: x = thp[a]; break;
        case 2u: x = t; break;
        case 16u: x = 1.0f; break;
        case 3u: x = A + B; break;
        case 4u: x = A - B; break;
        case 5u: x = A * B; break;
        case 6u: x = A / B; break;
        case 7u: x = sin(A); break;
        case 8u: x = cos(A); break;
        case 9u: x = exp(A); break;
        case 10u: x = log(fabs(A)); break;
        case 11u: x = A * A; break;
        case 12u: x = A * A * A; break;
        case 13u: x = sqrt(fabs(A)); break;
        case 14u: x = copysign(pow(fabs(A), 1.0f / 3.0f), A); break;
        case 15u: x = -A; break;
        default: break;
      }
      v[i] = x;
    }
    for (uint q = 0u; q < Q; q++) {
      const uint r = (q & 1u) ? (roots_w[q / 2u] >> 16) : (roots_w[q / 2u] & 0xFFFFu);
      const float x = v[r];
      if (!isfinite(x)) return false;
      E[q * RS_NPT + p] = x;
    }
  }
  float Qm[RS_MAXQ * RS_NPT];
  bool keep[RS_MAXQ];
  for (uint i = 0u; i < Q * RS_NPT; i++) Qm[i] = E[i];
  const float tau = (float)Q * 9.0f * 5.9604644775390625e-8f;
  for (uint k = 0u; k < Q; k++) {
    float on = 0.0f;
    for (uint p = 0u; p < RS_NPT; p++) on += E[k * RS_NPT + p] * E[k * RS_NPT + p];
    on = sqrt(on);
    for (uint i = 0u; i < k; i++) {
      if (!keep[i]) continue;
      float s = 0.0f;
      for (uint p = 0u; p < RS_NPT; p++) s += Qm[i * RS_NPT + p] * Qm[k * RS_NPT + p];
      for (uint p = 0u; p < RS_NPT; p++) Qm[k * RS_NPT + p] -= s * Qm[i * RS_NPT + p];
    }
    float nn = 0.0f;
    for (uint p = 0u; p < RS_NPT; p++) nn += Qm[k * RS_NPT + p] * Qm[k * RS_NPT + p];
    nn = sqrt(nn);
    keep[k] = nn > tau * on;
    for (uint p = 0u; p < RS_NPT; p++) Qm[k * RS_NPT + p] = keep[k] ? Qm[k * RS_NPT + p] / nn : 0.0f;
  }
  for (uint p = 0u; p < RS_NPT; p++) z[p] = ytg[p];
  for (uint k = 0u; k < Q; k++) {
    if (!keep[k]) continue;
    float s = 0.0f;
    for (uint p = 0u; p < RS_NPT; p++) s += Qm[k * RS_NPT + p] * z[p];
    for (uint p = 0u; p < RS_NPT; p++) z[p] -= s * Qm[k * RS_NPT + p];
  }
  return true;
}

kernel void post_sicis_residual(device const float* th0 [[buffer(0)]],
                                device const float* th1 [[buffer(1)]],
                                device const float* tb [[buffer(2)]],
                                device const uint* rec [[buffer(3)]],
                                device const uint* off [[buffer(4)]],
                                device const float* yv [[buffer(5)]],
                                device float* res [[buffer(6)]],
                                device uint* bad [[buffer(7)]],
                                constant uint& n [[buffer(8)]],
                                uint gid [[thread_position_in_grid]]) {
  if (gid >= n) return;
  float tbl[RS_NPT];
  for (uint p = 0u; p < RS_NPT; p++) tbl[p] = tb[gid * RS_NPT + p];
  float z0[RS_NPT], z1[RS_NPT];
  bool ok = rs_side(rec, off[2u * gid], th0 + gid * 32u, tbl, yv, z0);
  if (ok) ok = rs_side(rec, off[2u * gid + 1u], th1 + gid * 32u, tbl, yv + RS_NPT, z1);
  if (!ok) { bad[gid] = 1u; for (uint p = 0u; p < RS_NPT; p++) res[gid * RS_NPT + p] = 0.0f; return; }
  bad[gid] = 0u;
  for (uint p = 0u; p < RS_NPT; p++) {
    const float zx = z0[p], zy = z1[p];
    const float h = sqrt(zx * zx + zy * zy);
    res[gid * RS_NPT + p] = (zx + zy < 0.0f) ? -h : h;
  }
}

kernel void post_sicis_residual_n(device const float* th0 [[buffer(0)]],
                                device const float* th1 [[buffer(1)]],
                                device const float* tb [[buffer(2)]],
                                device const uint* rec [[buffer(3)]],
                                device const uint* off [[buffer(4)]],
                                device const float* yv [[buffer(5)]],
                                device float* res [[buffer(6)]],
                                device uint* bad [[buffer(7)]],
                                constant uint& n [[buffer(8)]],
                                device const float* nrm [[buffer(9)]],
                                uint gid [[thread_position_in_grid]]) {
  if (gid >= n) return;
  float tbl[RS_NPT];
  for (uint p = 0u; p < RS_NPT; p++) tbl[p] = tb[gid * RS_NPT + p];
  float z0[RS_NPT], z1[RS_NPT];
  bool ok = rs_side(rec, off[2u * gid], th0 + gid * 32u, tbl, yv, z0);
  if (ok) ok = rs_side(rec, off[2u * gid + 1u], th1 + gid * 32u, tbl, yv + RS_NPT, z1);
  if (!ok) { bad[gid] = 1u; for (uint p = 0u; p < RS_NPT; p++) res[gid * RS_NPT + p] = 0.0f; return; }
  bad[gid] = 0u;
  for (uint p = 0u; p < RS_NPT; p++) {
    res[gid * RS_NPT + p] = z0[p] * nrm[p] + z1[p] * nrm[RS_NPT + p];
  }
}

kernel void post_sicis_score(device const float* res [[buffer(0)]],
                             device const float* cfc [[buffer(1)]],
                             device float* out [[buffer(2)]],
                             device uint* flag [[buffer(3)]],
                             constant uint& n [[buffer(4)]],
                             uint gid [[thread_position_in_grid]]) {
  if (gid >= n) return;
  const float kI = 0.70710678118654752f;
  device const float* rowp = res + gid * 48u;
  float cur[48];
  bool bad = false;
  for (uint p = 0u; p < 48u; p++) { const float v = rowp[p]; bad = bad || !isfinite(v); cur[p] = v; }
  if (bad) { out[gid] = 1.0e30f; flag[gid] = 1u; return; }
  float E[48];
  uint len = 48u; uint eo = 0u;
  for (uint lev = 0u; lev < 4u; lev++) {
    const uint h = len / 2u;
    for (uint i = 0u; i < h; i++) {
      const float x = cur[2u * i], y = cur[2u * i + 1u];
      const float d = (x - y) * kI;
      cur[i] = (x + y) * kI;
      E[eo + i] = d * d;
    }
    eo += h; len = h;
  }
  {
    const float x = cur[0], y = cur[1];
    const float d = (x - y) * kI;
    E[45] = d * d;
    const float x2 = (x + y) * kI, y2 = cur[2];
    const float d2 = (x2 - y2) * kI;
    E[46] = d2 * d2;
    const float dc = (x2 + y2) * kI;
    E[47] = dc * dc;
  }
  const uint M = (uint)cfc[0];
  if (M == 0u) { out[gid] = 1.0e30f; flag[gid] = 0u; return; }
  const bool joint = cfc[4] == 1.0f;
  float g[8];
  float gmin = 1.0e38f;
  for (uint m = 0u; m < M; m++) {
    float sm = 0.0f;
    if (!joint) {
      for (uint c = 0u; c < 48u; c++) sm = fma(E[c], cfc[8u + m * 48u + c], sm);
    } else {
      const uint bs[7] = {0u, 24u, 36u, 42u, 45u, 46u, 47u};
      const uint bl[7] = {24u, 12u, 6u, 3u, 1u, 1u, 1u};
      for (uint b = 0u; b < 7u; b++) {
        float S = 0.0f;
        for (uint c = bs[b]; c < bs[b] + bl[b]; c++) S = fma(E[c], cfc[8u + m * 48u + c], S);
        if (!(S > 0.0f)) { out[gid] = 1.0e30f; flag[gid] = 1u; return; }
        sm += 0.5f * (float)bl[b] * log(S);
      }
    }
    g[m] = sm + cfc[392u + m];
    gmin = min(gmin, g[m]);
  }
  bool gbad = false;
  for (uint m = 0u; m < M; m++) gbad = gbad || !isfinite(g[m]);
  if (gbad) { out[gid] = 1.0e30f; flag[gid] = 1u; return; }
  float sw = 0.0f;
  for (uint m = 0u; m < M; m++) sw += exp(-(g[m] - gmin));
  const float r = ((gmin - log(sw) + log((float)M)) + cfc[3]) + cfc[2];
  const bool rbad = !isfinite(r) || !(sw > 0.0f) || r <= -1.0e30f;
  out[gid] = rbad ? 1.0e30f : min(r, 1.0e30f);
  flag[gid] = rbad ? 1u : 0u;
}
