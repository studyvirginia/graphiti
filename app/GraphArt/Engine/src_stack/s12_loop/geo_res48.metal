#include <metal_stdlib>
using namespace metal;

constant uint GM = 1024u;
constant uint GNMAX = 128u;
constant uint GTG = 256u;

struct GeoArgs { uint base; uint nchunk; uint sgn_tan; uint pad; float gsx; float gsy; float pad2; float pad3; };

struct MinIdx { float v; uint i; };
static inline MinIdx mi_comb(MinIdx L, MinIdx R) { return (R.v < L.v) ? R : L; }

static inline float4 geo_side(device const uint* proto, uint po, device const float* pl, float4 tv, thread bool& bad, thread float4& asum) {
  asum = fabs(tv);
  const uint n = proto[po], Q = proto[po + 1u], drv = proto[po + 2u], nth = proto[po + 3u];
  if (drv != 0u) return tv;
  if (n > GNMAX || Q == 0u || Q > 16u) { bad = true; return float4(0.0f); }
  device const uint* nd = proto + po + 5u;
  device const uint* rt = nd + 3u * n;
  float4 v[GNMAX];
  for (uint i = 0u; i < n; i++) {
    const uint w0 = nd[3u * i], w1 = nd[3u * i + 1u], w2 = nd[3u * i + 2u];
    const uint op = w0 & 0xFFu, ts = (w0 >> 8) & 0xFFu, a = w0 >> 16, b = w1;
    const float c = ts != 0u ? pl[ts - 1u] : as_type<float>(w2);
    const float4 A = v[a], B = v[b];
    float4 x;
    switch (op) {
      case 1u: x = float4(c); break;
      case 2u: x = tv; break;
      case 3u: x = A + B; break;
      case 4u: x = A - B; break;
      case 5u: x = A * B; break;
      case 6u: x = A / B; break;
      case 7u: x = precise::sin(A); break;
      case 8u: x = precise::cos(A); break;
      case 9u: x = precise::exp(A); break;
      case 10u: x = precise::log(fabs(A)); break;
      case 11u: x = A * A; break;
      case 12u: x = A * A * A; break;
      case 13u: x = precise::sqrt(fabs(A)); break;
      case 14u: x = copysign(precise::pow(fabs(A), 1.0f / 3.0f), A); break;
      default: x = -A; break;
    }
    v[i] = x;
  }
  float4 s = float4(0.0f);
  asum = float4(0.0f);
  for (uint q = 0u; q < Q; q++) { const float4 term = pl[nth + q] * v[rt[q] & 0xFFFFu]; s += term; asum += fabs(term); }
  return s;
}

kernel void geo_res48(device const uint* proto [[buffer(0)]],
                      device const uint* mrec [[buffer(1)]],
                      device const uint* moff [[buffer(2)]],
                      device const float* Y [[buffer(3)]],
                      device float* res [[buffer(4)]],
                      device uint* badv [[buffer(5)]],
                      device ushort* scratch [[buffer(6)]],
                      constant GeoArgs& A [[buffer(7)]],
                      device float* condv [[buffer(8)]],
                      uint grp [[threadgroup_position_in_grid]],
                      uint tid [[thread_position_in_threadgroup]],
                      uint sg [[simdgroup_index_in_threadgroup]],
                      uint lane [[thread_index_in_simdgroup]]) {
  threadgroup float tx[GM];
  threadgroup float ty[GM];
  threadgroup MinIdx sgt[2][8];
  threadgroup float sgv[8];
  threadgroup uint sgi[8];
  threadgroup uint jps[48];
  threadgroup atomic_uint badflag;
  const uint m = A.base + grp;
  if (tid == 0u) atomic_store_explicit(&badflag, 0u, memory_order_relaxed);
  threadgroup_barrier(mem_flags::mem_threadgroup);
  const uint mo = moff[m];
  const uint p0 = mrec[mo], p1 = mrec[mo + 1u];
  const float t0 = as_type<float>(mrec[mo + 2u]);
  device const float* pl0 = (device const float*)(mrec + mo + 3u);
  const uint nth0 = proto[p0 + 3u], Q0 = proto[p0 + 1u], drv0 = proto[p0 + 2u];
  device const float* pl1 = pl0 + (drv0 != 0u ? 0u : nth0 + Q0);
  const uint j0 = 4u * tid;
  const float4 tv = t0 + float4(float(j0), float(j0 + 1u), float(j0 + 2u), float(j0 + 3u)) * (1.0f / float(GM - 1u));
  bool bad = false;
  float4 ax4, ay4;
  float4 cx = geo_side(proto, p0, pl0, tv, bad, ax4) * A.gsx;
  float4 cy = geo_side(proto, p1, pl1, tv, bad, ay4) * A.gsy;
  const float4 am = max(ax4 * A.gsx, ay4 * A.gsy);
  const float amax = max(max(am.x, am.y), max(am.z, am.w));
  bad = bad || !all(isfinite(cx)) || !all(isfinite(cy));
  if (bad) atomic_store_explicit(&badflag, 1u, memory_order_relaxed);
  tx[j0] = cx.x; tx[j0 + 1u] = cx.y; tx[j0 + 2u] = cx.z; tx[j0 + 3u] = cx.w;
  ty[j0] = cy.x; ty[j0 + 1u] = cy.y; ty[j0 + 2u] = cy.z; ty[j0 + 3u] = cy.w;
  { const float sm = simd_max(amax); if (lane == 0u) sgv[sg] = sm; }
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (tid == 0u) { float c = sgv[0]; for (uint s = 1u; s < 8u; s++) c = max(c, sgv[s]); condv[m] = c; }
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (atomic_load_explicit(&badflag, memory_order_relaxed) != 0u) {
    if (tid == 0u) badv[m] = 1u;
    if (tid < 48u) res[m * 48u + tid] = 0.0f;
    return;
  }
  device ushort* sc = scratch + (ulong)grp * (48ul * GM);
  float4 Dp = float4(0.0f);
  const float INF = INFINITY;
  for (uint p = 0u; p < 48u; p++) {
    const float px = A.gsx * Y[p], py = A.gsy * Y[48u + p];
    const float4 dx = px - cx, dy = py - cy;
    const float4 best = dx * dx + dy * dy;
    if (p == 0u) { Dp = best; continue; }
    MinIdx ag; ag.v = Dp.x; ag.i = j0;
    if (Dp.y < ag.v) { ag.v = Dp.y; ag.i = j0 + 1u; }
    if (Dp.z < ag.v) { ag.v = Dp.z; ag.i = j0 + 2u; }
    if (Dp.w < ag.v) { ag.v = Dp.w; ag.i = j0 + 3u; }
    MinIdx inc = ag;
    for (uint off = 1u; off < 32u; off <<= 1) {
      MinIdx o; o.v = simd_shuffle_up(inc.v, off); o.i = simd_shuffle_up(inc.i, off);
      if (lane >= off) inc = mi_comb(o, inc);
    }
    if (lane == 31u) sgt[p & 1u][sg] = inc;
    MinIdx ex; ex.v = simd_shuffle_up(inc.v, 1u); ex.i = simd_shuffle_up(inc.i, 1u);
    if (lane == 0u) { ex.v = INF; ex.i = 0u; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    MinIdx pre; pre.v = INF; pre.i = 0u;
    for (uint s = 0u; s < sg; s++) pre = mi_comb(pre, sgt[p & 1u][s]);
    MinIdx cur = mi_comb(pre, ex);
    ushort4 fr;
    float dv[4] = {Dp.x, Dp.y, Dp.z, Dp.w};
    float nd[4];
    for (uint k = 0u; k < 4u; k++) {
      if (dv[k] < cur.v) { cur.v = dv[k]; cur.i = j0 + k; }
      fr[k] = (ushort)cur.i;
      nd[k] = cur.v + best[k];
    }
    *(device ushort4*)(sc + p * GM + j0) = fr;
    Dp = float4(nd[0], nd[1], nd[2], nd[3]);
  }
  float mv = Dp.x; uint mj = j0;
  if (Dp.y < mv) { mv = Dp.y; mj = j0 + 1u; }
  if (Dp.z < mv) { mv = Dp.z; mj = j0 + 2u; }
  if (Dp.w < mv) { mv = Dp.w; mj = j0 + 3u; }
  const float smin = simd_min(mv);
  if (lane == 0u) sgv[sg] = smin;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  float gmin = sgv[0]; for (uint s = 1u; s < 8u; s++) gmin = min(gmin, sgv[s]);
  const uint cand = (mv == gmin) ? mj : 0xFFFFFFFFu;
  const uint sidx = simd_min(cand);
  if (lane == 0u) sgi[sg] = sidx;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (tid == 0u) {
    uint j = sgi[0]; for (uint s = 1u; s < 8u; s++) j = min(j, sgi[s]);
    if (!isfinite(gmin) || j >= GM) { badv[m] = 1u; for (uint p = 0u; p < 48u; p++) jps[p] = 0u; }
    else {
      badv[m] = 0u;
      for (uint p = 47u; ; p--) { jps[p] = j; if (p == 0u) break; j = sc[p * GM + j]; }
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (tid < 48u) {
    const uint p = tid;
    const uint jp = jps[p];
    const float px = A.gsx * Y[p], py = A.gsy * Y[48u + p];
    float bd = INF, bzx = 0.0f, bzy = 0.0f, bcr = 0.0f;
    const uint s0 = jp > 0u ? jp - 1u : 0u, s1 = min(jp, GM - 2u);
    for (uint s = s0; s <= s1; s++) {
      const float ax = tx[s], ay = ty[s], ex = tx[s + 1u] - ax, ey = ty[s + 1u] - ay, L2 = ex * ex + ey * ey;
      float u = L2 > 0.0f ? ((px - ax) * ex + (py - ay) * ey) / L2 : 0.0f;
      u = min(1.0f, max(0.0f, u));
      const float zx = px - (ax + u * ex), zy = py - (ay + u * ey), d = zx * zx + zy * zy;
      if (d < bd) { bd = d; bzx = zx; bzy = zy; bcr = ex * zy - ey * zx; }
    }
    const float h = precise::sqrt(bzx * bzx + bzy * bzy);
    const bool neg = A.sgn_tan != 0u ? (bcr < 0.0f) : (bzx + bzy < 0.0f);
    res[m * 48u + p] = neg ? -h : h;
  }
}
