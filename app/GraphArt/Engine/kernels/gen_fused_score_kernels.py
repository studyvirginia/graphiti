#!/usr/bin/env python3
import sys, os
G = 3
_IM2_B = int(os.environ.get('FJN_IM2', '0') or 0)
IM2_ON0_B, IM2_LAD_B = bool(_IM2_B & 4), bool(_IM2_B & 8)
AFF_B = bool(os.environ.get('FJN_AFF'))
GUARD_B = bool(os.environ.get('FJN_GUARD'))
STORE_B = bool(os.environ.get('FJN_STORE'))
CISFUSE_B = bool(os.environ.get('FJN_CISFUSE'))
RES48_B = bool(os.environ.get('FJN_RES48')) or CISFUSE_B
DBG = bool(os.environ.get('FJN_DBG'))
import os as _os
RS = int(_os.environ.get('FJN_RS', '0'))
UNPROJ = bool(_os.environ.get('FJN_UNPROJ'))
LDL = _os.environ.get('FJN_LDL') == '1'
LAY = 'ovl' if _os.environ.get('FJN_OVL') == '1' else ('lpr' if _os.environ.get('FJN_LPR') == '1' else ('pts' if _os.environ.get('FJN_PTSTREAM') == '1' else 'base'))
SIMDRED = _os.environ.get('FJN_SIMDRED') == '1'
ABL_SOLVE = 'true' if _os.environ.get('FJN_ABL') == '1' else 'false'
SLOTH = _os.environ.get('FJN_SLOTH') == '1'
SLOT2 = _os.environ.get('FJN_SLOT') == '2'
_RS_T = set(int(v) for v in _os.environ['FJN_RS_T'].split(',')) if _os.environ.get('FJN_RS_T') else None
HTAN = bool(_os.environ.get('FJN_H_TAN'))
HPLAIN = bool(_os.environ.get('FJN_H_PLAIN'))
FN = 'fast::' if _os.environ.get('FJN_FASTMATH') else ''
FUSE = bool(_os.environ.get('FJN_FUSE'))
KW = 10
NW = 9
UNV = {7: '{F}sin({x})', 8: '{F}cos({x})', 9: '{F}exp({x})', 10: '{F}log(fabs({x}))', 11: '({x})*({x})', 12: '({x})*({x})*({x})', 13: '{F}sqrt(fabs({x}))', 14: 'copysign({F}pow(fabs({x}), 1.0f/3.0f), {x})', 15: '-({x})'}
UND = {7: '{F}cos({a})', 8: '-{F}sin({a})', 9: '{r}', 10: '1.0f / ({a})', 11: '2.0f * ({a})', 12: '3.0f * ({a}) * ({a})', 13: 'copysign(0.5f, {a}) / ({r})', 14: '1.0f / (3.0f * ({r}) * ({r}))', 15: '-1.0f'}
CTAB = "constant float4 ctab[4] = {float4(1,0,0,0), float4(0,1,1,0), float4(0,1,-1,0), float4(0,-1,-1,0)};"
def HS(Q): return (Q + 1) * 48 + (Q + 1) * 16
def each(f): return " ".join(f(v) for v in range(G))

def geom(Q, T):
    ND = T + 1
    SUB = 1 if ND <= 4 else (ND + 3) // 4
    CH = ((ND + 1) // 2) if SUB == 1 else 2
    CD = 1 + CH
    if FUSE and SUB > 1: CD = 1 + CH * SUB
    return ND, SUB, CH, CD

class Lay:
    def __init__(s): s.o = 0; s.d = {}
    def add(s, name, n): s.d[name] = s.o; s.o += n
    def __getitem__(s, k): return s.d[k]

def layout(Q, T):
    NPK = max(4, (T + 3) // 4 * 4); QP = max(4, (Q + 3) // 4 * 4)
    NC = T + KW; NMAX = 2 * T + KW; NPK_ = NMAX * (NMAX + 1) // 2
    L = Lay()
    L.add('THA', 2 * NPK); L.add('THC', 4 * NPK); L.add('UA', 16); L.add('UC', 32)
    L.add('ACA', 2 * QP); L.add('ACC', 4 * QP); L.add('TA', 48)
    L.add('RSSA', 2); L.add('RSSC', 4); L.add('SC', 32)
    L.add('DEL', 2 * NMAX); L.add('XB', 2 * NMAX); L.add('SCL', NMAX); L.add('GV', NMAX)
    L.add('HP', NPK_)
    if LAY == 'lpr':
        L.add('RES', 96); L.add('SINK', 4)
        L.o = (L.o + 3) // 4 * 4
    elif LAY == 'ovl':
        L.add('RES', 96); L.o = (L.o + 3) // 4 * 4
        if RES48_B: L.add('RESC', 4 * 48); L.add('RESA', 96); L.o = (L.o + 3) // 4 * 4
        if not SLOT2: L.add('OVL', 0)
        if SLOT2:
            L.add('SINK', 4); L.o = (L.o + 3) // 4 * 4; L.add('OVL', 0); L.add('LP', 0); L.add('COL', 0)
            L.d['COL'] = L['LP']; L.o = L['LP'] + max(2 * NPK_, (NC * (24 if SLOTH else 48) + 3) // 4 * 4)
        else:
            L.add('LP', 2 * NPK_); L.add('COL', NC * (24 if SLOTH else 48))
        L.add('OVE', 0); L.ovn = (L.o - L['OVL'] + 3) // 4 * 4
        L.o = L['OVL']
    else:
        L.add('LP', 2 * NPK_); L.add('RES', 96)
        if LAY == 'base': L.add('COL', NC * 48)
        L.o = (L.o + 3) // 4 * 4
    if RES48_B and LAY != 'ovl': L.add('RESC', 4 * 48); L.add('RESA', 96); L.o = (L.o + 3) // 4 * 4
    return L, NPK, QP, NC, NMAX, NPK_

OUT_ORIG_B = """  if (lw == 0u) {{
    po[0] = WR[{L['RSSA']}u + 0u]; po[1] = WR[{L['RSSA']}u + 1u]; po[2] = WR[{L['SC']}u + 10u]; po[3] = (float)it; po[4] = WR[{L['SC']}u + 12u]; po[5] = WR[{L['SC']}u + 17u]; po[6] = WR[{L['SC']}u + 16u]; po[7] = WR[{L['SC']}u + 18u];
  }}"""
OUT_LANES_B = """  if (lw < 8u) po[lw] = (lw == 0u) ? WR[{L['RSSA']}u + 0u] : (lw == 1u) ? WR[{L['RSSA']}u + 1u] : (lw == 2u) ? WR[{L['SC']}u + 10u] : (lw == 3u) ? (float)it : (lw == 4u) ? WR[{L['SC']}u + 12u] : (lw == 5u) ? WR[{L['SC']}u + 17u] : (lw == 6u) ? WR[{L['SC']}u + 16u] : WR[{L['SC']}u + 18u];"""
SIMDRED_F0 = 'simd_shuffle(rsS0, (ushort)0) + simd_shuffle(rsS1, (ushort)0)'
SIMDRED_F1 = 'simd_shuffle(rsS0, (ushort)16) + simd_shuffle(rsS1, (ushort)16)'
CISFUSE_TXT = r"""  // ---- g_cisfuse: SI-CIS score of the accepted residual (cis_scoring.hpp evaluateCisScoring, exactLikelihood 1, DC cell, ensemble mixture)
  // cfc: [0] M, [1] check flag, [2] Khi, [3] Klo, [4] joint flag (CIS_JOINT_COARSE), [8 + m*48 + c] = 0.5/p_mc, [392 + m] = K_m - Kmin.  score = -log(mean_m exp(-fit_m)); per-cell fit_m = sum_c E_c h_mc + K_m; joint fit_m = sum_b (D_b/2) log(sum_{c in b} E_c h_mc) + K_m
  simdgroup_barrier(mem_flags::mem_threadgroup);
  if (cfc[1] != 0.0f) { for (uint i = lw; i < 96u; i += 32u) res48[pair * 97u + i] = WR[@RESA@u + i]; }
  {
    threadgroup float* RA_ = WR + @RESA@u; threadgroup float* P0_ = WR + @RESC@u; threadgroup float* P1_ = WR + @RESC@u + 48u; threadgroup float* E_ = WR + @RESC@u + 96u;
    const float kI_ = 0.70710678118654752f;
    bool nf_ = false;
    for (uint p = lw; p < 48u; p += 32u) {
      const float zx = RA_[p], zy = RA_[48u + p]; const float h = precise::sqrt(zx * zx + zy * zy);
      const float c = (zx + zy < 0.0f) ? -h : h; P0_[p] = c; nf_ = nf_ || ((as_type<uint>(c) & 0x7f800000u) == 0x7f800000u);
    }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    if (lw < 24u) { const float a = P0_[2u * lw], b = P0_[2u * lw + 1u]; const float d = (a - b) * kI_; P1_[lw] = (a + b) * kI_; E_[lw] = d * d; }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    if (lw < 12u) { const float a = P1_[2u * lw], b = P1_[2u * lw + 1u]; const float d = (a - b) * kI_; P0_[lw] = (a + b) * kI_; E_[24u + lw] = d * d; }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    if (lw < 6u) { const float a = P0_[2u * lw], b = P0_[2u * lw + 1u]; const float d = (a - b) * kI_; P1_[lw] = (a + b) * kI_; E_[36u + lw] = d * d; }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    if (lw < 3u) { const float a = P1_[2u * lw], b = P1_[2u * lw + 1u]; const float d = (a - b) * kI_; P0_[lw] = (a + b) * kI_; E_[42u + lw] = d * d; }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    if (lw == 0u) {   // n = 3: one butterfly + identity carry; n = 2: one butterfly; DC = last approximation
      const float a = P0_[0], b = P0_[1]; const float d = (a - b) * kI_; E_[45] = d * d; P1_[0] = (a + b) * kI_; P1_[1] = P0_[2];
      const float a2 = P1_[0], b2 = P1_[1]; const float d2 = (a2 - b2) * kI_; E_[46] = d2 * d2; const float dc = (a2 + b2) * kI_; E_[47] = dc * dc;
    }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    const uint M_ = (uint)cfc[0];
    float g_ = 1.0e38f; bool bad_ = nf_;
    if (lw < M_) {
      float sm = 0.0f;
      if (cfc[4] == 0.0f) { for (uint c = 0; c < 48u; c++) sm = fma(E_[c], cfc[8u + lw * 48u + c], sm); }
      else {   // CIS_JOINT_COARSE (S4): per Haar level block b, D_b = size, sum_b (D_b/2) log(S_mb/2), S_mb/2 = sum_c E_c h_mc
        uint c0_ = 0u;
        for (uint b_ = 0u; b_ < 7u; b_++) {
          const uint n_ = (b_ == 0u) ? 24u : (b_ == 1u) ? 12u : (b_ == 2u) ? 6u : (b_ == 3u) ? 3u : 1u;
          float sb_ = 0.0f; for (uint c = c0_; c < c0_ + n_; c++) sb_ = fma(E_[c], cfc[8u + lw * 48u + c], sb_);
          sm += (0.5f * (float)n_) * precise::log(sb_); c0_ += n_;
        }
      }
      g_ = sm + cfc[392u + lw]; bad_ = bad_ || ((as_type<uint>(g_) & 0x7f800000u) == 0x7f800000u);
    }
    const bool anyBad_ = simd_any(bad_) || (M_ == 0u);
    const float gmin_ = simd_min(anyBad_ ? 0.0f : g_);
    const float w_ = (lw < M_ && !anyBad_) ? precise::exp(-(g_ - gmin_)) : 0.0f;
    const float sw_ = simd_sum(w_);
    if (lw == 0u) {
      float r_ = ((gmin_ - precise::log(sw_) + precise::log((float)M_)) + cfc[3]) + cfc[2];
      const bool rbad_ = anyBad_ || ((as_type<uint>(r_) & 0x7f800000u) == 0x7f800000u) || !(sw_ > 0.0f);
      res48[pair * 97u + 96u] = rbad_ ? 1.0e30f : fmin(r_, 1.0e30f);
    }
  }
"""
def kern(Q, T):
    global RS
    RS_G = RS
    if _RS_T is not None: RS = 1 if T in _RS_T else 0
    try: return _kern(Q, T)
    finally: RS = RS_G

def _kern(Q, T):
    ND, SUB, CH, CD = geom(Q, T)
    L, NPK, QP, NC, NMAX, NPKK = layout(Q, T)
    WS = L.o; OVN = getattr(L, 'ovn', 0)
    if LAY in ('lpr', 'pts'): L.d['COL'] = WS
    COLW = 'CB[' if LAY in ('lpr', 'pts') else (f"((threadgroup half*)(WR + {L['COL']}u))[" if SLOTH else f"WR[{L['COL']}u + ")
    COLR = 'CB' if LAY in ('lpr', 'pts') else (f"((threadgroup const half*)(WR + {L['COL']}u))" if SLOTH else f"(WR + {L['COL']}u)")
    COLT = 'threadgroup const half*' if SLOTH else 'threadgroup const float*'
    if SLOTH or SLOT2: assert LAY == 'ovl'
    SINKW = f"WR[{L['SINK']}u + (e & 3u)]" if (LAY == 'lpr' or SLOT2) else f"WR[{L['LP']}u + (e & 255u)]"
    LPHX = 'R' if LAY == 'lpr' else f"WR + {L['LP']}u + hf * {NPKK}u"
    OVLX = f"max({OVN}u, 2u * per)" if LAY == 'ovl' else "2u * per"
    SLOTC = 2 * CD * 48 + WS
    M2 = (NC + 1) // 2
    def fetchB(w, X, comps):
        return (f"      {{ const uint of_ = {w} & 0x1FFFFFFFu; const bool vec_ = (({w} >> 30) & 1u) == 0u; threadgroup const float* bp_ = (({w} >> 29) & 1u) ? TSH : R; const uint b_ = of_ + (vec_ ? p0 : 0u); const uint cs_ = vec_ ? 48u : nsc2; const uint ps_ = vec_ ? 1u : 0u;\n"
                f"        " + " ".join(f"{X}{c}_{v} = bp_[b_ + {c}u * cs_ + {v}u * ps_];" for c in comps for v in range(G)) + " }\n")
    def fetchA(w, X, comps):
        return f"      if ({w} >> 31) {{ " + " ".join(f"{X}{c}_{v} = acc{c}_{v};" for c in comps for v in range(G)) + " }\n      else" + fetchB(w, X, comps)
    def step(comps, PRG):
        C = len(comps)
        s = f"    {{ const uint4 w = {PRG}[hx + i];\n"
        s += "      float " + ", ".join([f"a{c}_{v}" for c in comps for v in range(G)] + [f"b{c}_{v}" for c in comps for v in range(G)]) + ";\n"
        s += fetchA('w.z', 'a', comps) + fetchB('w.w', 'b', comps)
        s += "      if (w.x == 0u) { " + " ".join(f"acc{c}_{v} = fma(a{c}_{v}, b0_{v}, a0_{v} * b{c}_{v});" for c in comps if c > 0 for v in range(G)) + " " + " ".join(f"acc0_{v} = a0_{v} * b0_{v};" for v in range(G)) + " }\n"
        if AFF_B: s += "      else if (w.x < 4u) { " + " ".join(f"acc{c}_{v} = (w.x == 1u) ? (a{c}_{v} + b{c}_{v}) : (a{c}_{v} - b{c}_{v}) * ((w.x == 2u) ? 1.0f : -1.0f);" for c in comps for v in range(G)) + " }\n"
        else: s += "      else if (w.x < 4u) { const float4 k_ = ctab[w.x]; " + " ".join(f"acc{c}_{v} = k_.y * fma(k_.z, b{c}_{v}, a{c}_{v});" for c in comps for v in range(G)) + " }\n"
        s += "      else { switch (w.x) {\n"
        s += "      case 4: " + " ".join(f"acc0_{v} = a0_{v} / b0_{v};" for v in range(G)) + " " + " ".join(f"acc{c}_{v} = (a{c}_{v} - acc0_{v} * b{c}_{v}) / b0_{v};" for c in comps if c > 0 for v in range(G)) + " break;\n"
        s += "      case 5: " + " ".join(f"acc0_{v} = b0_{v} / a0_{v};" for v in range(G)) + " " + " ".join(f"acc{c}_{v} = (b{c}_{v} - acc0_{v} * a{c}_{v}) / a0_{v};" for c in comps if c > 0 for v in range(G)) + " break;\n"
        for op in range(7, 16):
            s += f"      case {op}: {{ " + " ".join(f"const float m{v} = a0_{v} * b0_{v};" for v in range(G)) + " " + " ".join(f"acc0_{v} = {UNV[op].format(x=f'm{v}', F=FN)};" for v in range(G))
            if C > 1: s += " { " + " ".join(f"const float d{v} = {UND[op].format(a=f'm{v}', r=f'acc0_{v}', F=FN)};" for v in range(G)) + " " + " ".join(f"acc{c}_{v} = d{v} * fma(a{c}_{v}, b0_{v}, a0_{v} * b{c}_{v});" for c in comps if c > 0 for v in range(G)) + " }"
            s += " } break;\n"
        s += "      default: " + " ".join(f"acc{c}_{v} = a{c}_{v};" for c in comps for v in range(G)) + " break; } }\n"
        def rw(c, v):
            hr = (HTAN and c > 0) or (HPLAIN and c == 0)
            return f"(isinf((float)(half)acc{c}_{v}) || isnan((float)(half)acc{c}_{v}) ? acc{c}_{v} : (float)(half)acc{c}_{v})" if hr else f"acc{c}_{v}"
        s += "      if (w.y != 0xFFFFFFFFu) { " + " ".join(f"R[w.y + {c * 48}u + p0 + {v}u] = {rw(c, v)};" for c in comps for v in range(G)) + " } }\n"
        return s

    def plain_pass():
        comps = [0]
        P = []
        P.append("  for (uint s = 0; s < 2u; s++) {")
        P.append("    const uint rec = s ? rc.y : rc.x; const uint4 h = hdr[2u * rec]; const uint nS = h.y & 0xFFFFu, Qs = h.y >> 16, nsc2 = h.w & 255u, SBp = (h.w >> 8) & 0xFFFFu; const uint hx = h.x;")
        P.append("    device const uint* sg = segs + rec * 72u; threadgroup const float* TH = WR + " + f"{L['THC']}u + (hf * 2u + s) * {NPK}u;")
        P.append("    for (uint j = lq; j < nS; j += 16u) { const ushort w = sprog[h.z + j]; const uint op = w & 63u, a = (w >> 6) & 31u; R[SBp + j] = (op == 1u) ? TH[a] : 1.0f; }")
        P.append("    if (lq == 0u) { R[SBp + nS] = 1.0f; R[SBp + nS + 1u] = 2.0f; }")
        P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append("    float acc0_0 = 0, acc0_1 = 0, acc0_2 = 0;")
        P.append("    float " + ", ".join(f"E{q}_{v} = 0" for q in range(Q) for v in range(G)) + ";")
        P.append(f"    for (uint q = 0; q < {Q}u; q++) {{")
        P.append(f"    for (uint i = sg[q]; i < sg[q + 1]; i++)\n" + step(comps, 'progP'))
        P.append(f"    {{ const uint cw = sg[18 + q]; float x0_0, x0_1, x0_2;\n" + fetchA('cw', 'x', comps))
        P.append("      if (q >= Qs) { x0_0 = 0.0f; x0_1 = 0.0f; x0_2 = 0.0f; }")
        P.append("      switch (q) { " + " ".join(f"case {q}: " + each(lambda v: f"E{q}_{v} = x0_{v};") + " break;" for q in range(Q)) + " }")
        P.append("    } }")
        P.append("    float Y_0 = yv[s * 48u + p0], Y_1 = yv[s * 48u + p0 + 1u], Y_2 = yv[s * 48u + p0 + 2u];")
        P.append("    bool bad = false;")
        for q in range(Q): P.append("    bad = bad || " + " || ".join(f"!isfinite(E{q}_{v})" for v in range(G)) + ";")
        cols = [f"E{q}" for q in range(Q)] + ["Y"]
        n = len(cols)
        P.append("    " + " ".join(f"float t{k} = (lq == {k}u) ? 1.0f : 0.0f;" for k in range(n - 1)) + " float tY = 0.0f;")
        P.append("    " + " ".join(f"const float on{k} = sqrt(hsum16(" + " + ".join(f"{cols[k]}_{v} * {cols[k]}_{v}" for v in range(G)) + "));" for k in range(n - 1)))
        tn = [f"t{k}" for k in range(n - 1)] + ["tY"]
        for k in range(n - 1):
            ck = cols[k]
            P.append("    {")
            if IM2_ON0_B and k == 0: P.append(f"     const float nrm = on0;")
            else: P.append(f"     const float nn = hsum16(" + " + ".join(f"{ck}_{v} * {ck}_{v}" for v in range(G)) + "); const float nrm = sqrt(nn);")
            P.append(f"     const bool keep = nrm > tau * on{k}; const float inv = keep ? 1.0f / nrm : 0.0f; {tn[k]} *= inv;")
            P.append("     " + each(lambda v: f"{ck}_{v} *= inv;"))
            for j in range(k + 1, n):
                cj = cols[j]
                P.append(f"     {{ const float r = hsum16(" + " + ".join(f"{ck}_{v} * {cj}_{v}" for v in range(G)) + f"); {tn[j]} = fma(-r, {tn[k]}, {tn[j]}); " + each(lambda v: f"{cj}_{v} -= r * {ck}_{v};") + " }")
            P.append("    }")
        P.append(f"    {{ device float* EBh = EB + ebase + (s * 2u + hf) * {HS(Q)}u + p0; " + " ".join(f"EBh[{k * 48 + v}u] = E{k}_{v};" for k in range(Q) for v in range(G)) + " " + " ".join(f"EBh[{Q * 48 + v}u] = Y_{v};" for v in range(G)) + " }")
        P.append("    const float rss = hsum16(" + " + ".join(f"Y_{v} * Y_{v}" for v in range(G)) + "); const bool badh = bad_half(bad);")
        P.append(f"    if (lq < {QP}u) WR[{L['ACC']}u + (hf * 2u + s) * {QP}u + lq] = (lq < {Q}u) ? -tY : 0.0f;")
        if RES48_B: P.append(f"    {{ const uint rb_ = {L['RESC']}u + (hf * 2u + s) * 48u + p0; WR[rb_] = Y_0; WR[rb_ + 1u] = Y_1; WR[rb_ + 2u] = Y_2; }}")
        P.append(f"    if (lq == 0u) WR[{L['RSSC']}u + hf * 2u + s] = (badh || !isfinite(rss)) ? INFINITY : rss;")
        if SIMDRED: P.append(f"    if (s == 0u) rsS0 = (badh || !isfinite(rss)) ? INFINITY : rss; else rsS1 = (badh || !isfinite(rss)) ? INFINITY : rss;")
        P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append("  }")
        return "\n".join(P)

    def dual_pass():
        comps = list(range(CD))
        P = []
        P.append("  for (uint s = 0; s < 2u; s++) {")
        P.append("    const uint rec = s ? rc.y : rc.x; const uint4 h = hdr[2u * rec + 1u]; const uint nS = h.y & 0xFFFFu, Qs = h.y >> 16, nsc2 = h.w & 255u, SBd = (h.w >> 8) & 0xFFFFu; const uint hx = h.x;")
        P.append("    device const uint* sg = segs + rec * 72u + 36u; threadgroup const float* TH = WR + " + f"{L['THA']}u + s * {NPK}u;")
        FU = FUSE and SUB > 1
        if SUB > 1:
            P.append("    float " + ", ".join(f"V{j}_{v} = 0" for j in range(ND) for v in range(G)) + ";")
        if SUB > 1 and not FU:
            P.append(f"    for (uint sp_ = 0; sp_ < {SUB}u; sp_++) {{")
        else:
            P.append("    const uint sp_ = 0u;")
        if FU: dirx = lambda c: f"({((c - 1) // CH) * 4}u + hf * {CH}u + {(c - 1) % CH}u)"
        else: dirx = lambda c: f"(sp_ * 4u + hf * {CH}u + {c - 1}u)" if SUB > 1 else f"(hf * {CH}u + {c - 1}u)"
        P.append("    for (uint j = lq; j < nS; j += 16u) { const ushort w = sprog[h.z + j]; const uint op = w & 63u, a = (w >> 6) & 31u; const float x = (op == 1u) ? TH[a] : 1.0f; R[SBd + j] = x; "
                 + " ".join(f"R[SBd + {c}u * nsc2 + j] = (op == 1u && a == {dirx(c)}) ? 1.0f : 0.0f;" for c in comps if c > 0) + " }")
        P.append("    if (lq < " + str(CD) + "u) { R[SBd + lq * nsc2 + nS] = lq == 0u ? 1.0f : 0.0f; R[SBd + lq * nsc2 + nS + 1u] = lq == 0u ? 2.0f : 0.0f; }")
        P.append("    for (uint i = lq; i < 48u; i += 16u) { TSH[i] = WR[" + str(L['TA']) + "u + i]; " + " ".join(f"TSH[{c * 48}u + i] = ({dirx(c)} == {T}u) ? 1.0f : 0.0f;" for c in comps if c > 0) + " }")
        P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append("    float " + ", ".join(f"acc{c}_{v} = 0" for c in comps for v in range(G)) + ";")
        NWV = CD - 1 if FU else CH
        P.append("    float " + ", ".join(f"Wv{j}_{v} = 0" for j in range(NWV) for v in range(G)) + ";")
        P.append(f"    for (uint q = 0; q < {Q}u; q++) {{")
        P.append(f"    for (uint i = sg[q]; i < sg[q + 1]; i++)\n" + step(comps, 'progD'))
        if GUARD_B: P.append("    if (q < Qs) {")
        P.append(f"    {{ const uint cw = sg[18 + q]; float " + ", ".join(f"x{c}_{v}" for c in comps for v in range(G)) + ";\n" + fetchA('cw', 'x', comps))
        P.append("      { const float aq = (q < Qs) ? WR[" + str(L['ACA']) + "u + s * " + str(QP) + "u + q] : 0.0f; " + " ".join(f"{{ const float t_ = aq * x{c}_{v}; Wv{c-1}_{v} += isfinite(t_) ? t_ : 0.0f; }}" for c in comps if c > 0 for v in range(G)) + " }")
        P.append("    } }" + (" }" if GUARD_B else ""))
        if FU:
            P.append("    " + " ".join(f"{{ const float o_ = simd_shuffle_xor(Wv{2 * k + j % 2}_{v}, 16); V{4 * k + j}_{v} = (hf == {j // 2}u) ? Wv{2 * k + j % 2}_{v} : o_; }}" for k in range(SUB) for j in range(4) if 4 * k + j < ND for v in range(G)))
        elif SUB > 1:
            P.append("    switch (sp_) { " + " ".join(f"case {k}: " + " ".join(f"{{ const float o_ = simd_shuffle_xor(Wv{j % 2}_{v}, 16); V{4 * k + j}_{v} = (hf == {j // 2}u) ? Wv{j % 2}_{v} : o_; }}" for j in range(4) if 4 * k + j < ND for v in range(G)) + " break;" for k in range(SUB)) + " }")
            P.append("    }")
        else:
            for j in range(ND):
                src = j % CH; fromhalf = j // CH
                P.append(f"    float " + ", ".join(f"V{j}_{v}" for v in range(G)) + "; { " + " ".join(f"{{ const float o_ = simd_shuffle_xor(Wv{src}_{v}, 16); V{j}_{v} = (hf == {fromhalf}u) ? Wv{src}_{v} : o_; }}" for v in range(G)) + " }")
        if LAY == 'ovl': P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append(f"    device const float* EBc = EB + ebase + (s * 2u + curh) * {HS(Q)}u + p0;")
        P.append("    const float UU = WR[" + str(L['SC']) + "u + 3u];")
        P.append("    const float t0v = WR[" + str(L['UA']) + "u + 9u];")
        for m in range(M2):
            P.append(f"    const uint idx{m} = {2 * m}u + hf; float c{m}_0, c{m}_1, c{m}_2;")
            if 2 * m < T:
                sel = (lambda v: f"(hf ? V{2*m+1}_{v} : V{2*m}_{v})") if 2 * m + 1 < T else (lambda v: f"V{2*m}_{v}")
                P.append(f"    if (idx{m} < {T}u) {{ " + " ".join(f"c{m}_{v} = {sel(v)};" for v in range(G)) + " }")
                P.append("    else {")
            else:
                P.append("    {")
            P.append(f"      const uint k_ = idx{m} - {T}u; " + " ".join(
                f"{{ const float dd_ = (k_ < 9u) ? 2.0f * WR[{L['UA']}u + k_] / UU * (Sm[(p0 + {v}u) * 9u + k_] - (TSH[p0 + {v}u] - t0v)) : 1.0f; c{m}_{v} = (idx{m} < {NC}u) ? V{T}_{v} * dd_ : 0.0f; }}" for v in range(G)) + " }")
        P.append("    " + " ".join(f"const float on{m} = sqrt(hsum16(" + " + ".join(f"c{m}_{v} * c{m}_{v}" for v in range(G)) + "));" for m in range(M2)))
        if UNPROJ:
            P.append("    " + " ".join(f"if (lq == 0u && idx{m} < {T}u && (s == 0u ? idx{m} < nx : idx{m} < ny)) WR[{L['DEL']}u + (s == 0u ? idx{m} : nx + idx{m})] += on{m} * on{m};" for m in range(M2)))
            P.append("    " + " ".join(f"if (lq == 0u && idx{m} >= {T}u && idx{m} < {T}u + nk) WR[{L['DEL']}u + nx + ny + idx{m} - {T}u] += on{m} * on{m};" for m in range(M2)))
        P.append(f"    for (uint k = 0; k < {Q}u; k++) {{ " + " ".join(f"const float q{v} = EBc[k * 48u + {v}u];" for v in range(G)))
        for m in range(M2):
            P.append(f"      {{ const float r = hsum16(" + " + ".join(f"q{v} * c{m}_{v}" for v in range(G)) + "); " + each(lambda v: f"c{m}_{v} -= r * q{v};") + " }")
        P.append("    }")
        for m in range(M2):
            P.append(f"    {{ const float nrm = sqrt(hsum16(" + " + ".join(f"c{m}_{v} * c{m}_{v}" for v in range(G)) + f")); const bool keep = nrm > tau * on{m}; if (!keep) {{ c{m}_0 = 0; c{m}_1 = 0; c{m}_2 = 0; }}")
            P.append(f"      if (idx{m} < {NC}u) {{ " + " ".join(f"{COLW}idx{m} * 48u + p0 + {v}u] = c{m}_{v};" for v in range(G)) + " } }")
        P.append("    if (hf == 0u) { " + " ".join(f"WR[{L['RES']}u + s * 48u + p0 + {v}u] = EBc[{Q * 48 + v}u];" for v in range(G)) + " }")
        P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append("    // ---- Gram / gradient contributions of this side (COL holds this side only)")
        P.append("    for (uint e = lw; e < nn; e += 32u) {")
        P.append("      uint i_ = (uint)((sqrt(8.0f * (float)e + 1.0f) - 1.0f) * 0.5f); while (i_ * (i_ + 1u) / 2u > e) i_--; while ((i_ + 1u) * (i_ + 2u) / 2u <= e) i_++;")
        P.append("      const uint j_ = e - i_ * (i_ + 1u) / 2u;")
        P.append("      const bool ai = (i_ < nx) ? (s == 0u) : (i_ < nx + ny) ? (s == 1u) : true; const bool aj = (j_ < nx) ? (s == 0u) : (j_ < nx + ny) ? (s == 1u) : true;")
        P.append("      if (ai && aj) {")
        P.append(f"        const uint ci = (i_ < nx) ? i_ : (i_ < nx + ny) ? (i_ - nx) : ({T}u + (i_ - nx - ny)); const uint cj = (j_ < nx) ? j_ : (j_ < nx + ny) ? (j_ - nx) : ({T}u + (j_ - nx - ny));")
        P.append(f"        {COLT} ca = {COLR} + ci * 48u; {COLT} cb = {COLR} + cj * 48u; float hs = 0;")
        P.append(f"        for (uint p = 0; p < 48u; p++) hs = fma(float(ca[p]), float(cb[p]), hs);")
        P.append(f"        WR[{L['HP']}u + e] += hs;")
        P.append(f"        for (uint gr = 0; gr < rp.z; gr++) {{ float hs2 = 0; for (uint p = 0; p < 48u; p++) hs2 = fma(float(ca[p]), float(cb[p]), hs2); {SINKW} = hs2; }} }}")
        P.append("    }")
        P.append("    for (uint i_ = lw; i_ < n; i_ += 32u) {")
        P.append("      const bool ai = (i_ < nx) ? (s == 0u) : (i_ < nx + ny) ? (s == 1u) : true;")
        P.append("      if (ai) {")
        P.append(f"        const uint ci = (i_ < nx) ? i_ : (i_ < nx + ny) ? (i_ - nx) : ({T}u + (i_ - nx - ny));")
        P.append(f"        {COLT} ca = {COLR} + ci * 48u; threadgroup const float* rr = WR + {L['RES']}u + s * 48u; float gs = 0;")
        P.append(f"        for (uint p = 0; p < 48u; p++) gs = fma(float(ca[p]), rr[p], gs);")
        P.append(f"        WR[{L['GV']}u + i_] += gs; }}")
        P.append("    }")
        P.append("    simdgroup_barrier(mem_flags::mem_threadgroup);")
        P.append("  }")
        return "\n".join(P)


    def gj_text(NR):
        P = []
        P.append("        { const uint base_ = hf * 16u; const bool live = lq < n; const float si = live ? WR[%du + lq] : 0.0f;" % L['SCL'])
        P.append("          float " + ", ".join(f"a{j}" for j in range(NR)) + ", bb, dg = 1.0f;")
        for j in range(NR):
            P.append(f"          {{ const bool lv = live && {j}u < n; const float sj = lv ? WR[{L['SCL']}u + {j}u] : 0.0f; const float h = lv ? WR[{L['HP']}u + ((lq >= {j}u) ? (lq * (lq + 1u) / 2u + {j}u) : ({j}u * ({j}u + 1u) / 2u + lq))] : 0.0f;")
            P.append(f"            a{j} = (lq == {j}u) ? (live ? (h * si * si + ((si > 0.0f) ? lamh : 1.0f)) : 1.0f) : si * sj * h; }}")
        P.append(f"          bb = live ? si * WR[{L['GV']}u + lq] : 0.0f;")
        GA = int(_os.environ.get('FJN_GJABL','0'))
        if GA: pass
        elif RS == 1:
          for k in range(NR):
            P.append(f"          {{ const float pk = simd_shuffle(a{k}, (ushort)(base_ + {k}u)); okc = okc && (pk > 1e-12f) && isfinite(pk); const float ip = 1.0f / max(pk, 1e-12f);")
            P.append(f"            const float f = (lq == {k}u) ? 0.0f : a{k} * ip; if (lq == {k}u) dg = pk;")
            for j in range(k + 1, NR):
                P.append(f"            a{j} = fma(-f, simd_shuffle(a{j}, (ushort)(base_ + {k}u)), a{j});")
            P.append(f"            bb = fma(-f, simd_shuffle(bb, (ushort)(base_ + {k}u)), bb); }}")
        else:
          P.append("          bool okl = true;")
          for k in range(NR):
            P.append(f"          {{ if (lq == {k}u) {{ okl = (a{k} > 1e-12f) && isfinite(a{k}); dg = a{k}; }} const float ip = 1.0f / max(a{k}, 1e-12f); const float f = (lq == {k}u) ? 0.0f : a{k};")
            for j in range(k + 1, NR):
                P.append(f"            a{j} = fma(-f, simd_shuffle(a{j} * ip, (ushort)(base_ + {k}u)), a{j});")
            P.append(f"            bb = fma(-f, simd_shuffle(bb * ip, (ushort)(base_ + {k}u)), bb); }}")
          P.append("          okc = okc && !bad_half(!okl);")
        P.append(f"          if (live) XBh[lq] = bb / max(dg, 1e-12f); }}")
        return "\n".join(P)

    plain = plain_pass(); dual = dual_pass()
    RSSDECL = "      float rsS0 = INFINITY, rsS1 = INFINITY;\n" if SIMDRED else ""
    ZDEL = f"for (uint i_ = lw; i_ < n; i_ += 32u) WR[{L['DEL']}u + i_] = 0.0f;" if UNPROJ else ''
    SCLF = (f"for (uint i_ = lw; i_ < n; i_ += 32u) {{ const float d_ = WR[{L['DEL']}u + i_]; WR[{L['SCL']}u + i_] = (d_ > 0.0f && isfinite(d_)) ? rsqrt(d_) : 0.0f; }}") if UNPROJ else \
            (f"for (uint i_ = lw; i_ < n; i_ += 32u) {{ const float d_ = WR[{L['HP']}u + i_ * (i_ + 1u) / 2u + i_]; WR[{L['SCL']}u + i_] = (d_ > 0.0f && isfinite(d_)) ? rsqrt(d_) : 0.0f; }}")
    CHOL = r"""      for (uint k = 0; k < n; k++) {
        const float piv = LPh[k * (k + 1u) / 2u + k]; okc = okc && (piv > 1e-12f) && isfinite(piv);
        const float lkk = sqrt(max(piv, 1e-12f)); const float il = 1.0f / lkk;
        simdgroup_barrier(mem_flags::mem_threadgroup);
        for (uint i_ = k + 1u + lq; i_ < n; i_ += 16u) LPh[i_ * (i_ + 1u) / 2u + k] *= il;
        if (lq == 0u) LPh[k * (k + 1u) / 2u + k] = lkk;
        simdgroup_barrier(mem_flags::mem_threadgroup);
        for (uint i_ = k + 1u + lq; i_ < n; i_ += 16u) { const float lik = LPh[i_ * (i_ + 1u) / 2u + k]; for (uint j_ = k + 1u; j_ <= i_; j_++) LPh[i_ * (i_ + 1u) / 2u + j_] -= lik * LPh[j_ * (j_ + 1u) / 2u + k]; }
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }
      // forward: L y = b (b in XBh), then back: L^T x = y
      for (uint k = 0; k < n; k++) {
        const float yk = XBh[k] / LPh[k * (k + 1u) / 2u + k];
        simdgroup_barrier(mem_flags::mem_threadgroup);
        if (lq == 0u) XBh[k] = yk;
        for (uint i_ = k + 1u + lq; i_ < n; i_ += 16u) XBh[i_] -= LPh[i_ * (i_ + 1u) / 2u + k] * yk;
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }
      for (uint kk = 0; kk < n; kk++) {
        const uint k = n - 1u - kk;
        const float xk = XBh[k] / LPh[k * (k + 1u) / 2u + k];
        simdgroup_barrier(mem_flags::mem_threadgroup);
        if (lq == 0u) XBh[k] = xk;
        for (uint i_ = lq; i_ < k; i_ += 16u) XBh[i_] -= LPh[k * (k + 1u) / 2u + i_] * xk;
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }
"""
    LDL_TXT = r"""      for (uint k = 0; k < n; k++) {
        const uint kk0 = k * (k + 1u) / 2u;
        const float piv = LPh[kk0 + k]; okc = okc && (piv > 1e-12f) && isfinite(piv);
        const float id = 1.0f / max(piv, 1e-12f); const float cb = XBh[k] * id;
        for (uint i_ = k + 1u + lq; i_ < n; i_ += 16u) {
          const uint ii = i_ * (i_ + 1u) / 2u; const float aik = LPh[ii + k]; const float f = aik * id;
          XBh[i_] -= aik * cb;
          for (uint j_ = k + 1u; j_ <= i_; j_++) LPh[ii + j_] -= f * LPh[j_ * (j_ + 1u) / 2u + k];
        }
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }
      for (uint kk = 0; kk < n; kk++) {
        const uint k = n - 1u - kk;
        const float xk = XBh[k] / LPh[k * (k + 1u) / 2u + k];
        if (lq == 0u) XBh[k] = xk;
        for (uint i_ = lq; i_ < k; i_ += 16u) XBh[i_] -= LPh[k * (k + 1u) / 2u + i_] * xk;
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }
"""
    SOLVE = LDL_TXT if LDL else CHOL
    GJ8 = gj_text(10); GJ12 = gj_text(12); GJ14 = gj_text(14); GJ16 = gj_text(16); RS_ON = 'true' if RS else 'false'
    NM = NMAX
    def Hidx(i, j): return f"(({i}) * (({i}) + 1u) / 2u + ({j}))"
    DBGF = 'true' if DBG else 'false'
    DBGG = (f"if (it == 0u) {{ for (uint i_ = lw; i_ < n; i_ += 32u) {{ po[40u + i_] = WR[{L['GV']}u + i_]; po[72u + i_] = WR[{L['HP']}u + i_ * (i_ + 1u) / 2u + i_]; }} }}") if DBG else ''
    DBGD = (f"if (it == 0u && rungp == 0u) {{ for (uint i_ = lq; i_ < n; i_ += 16u) {{ if (hf == 0u) po[8u + i_] = WR[{L['DEL']}u + i_]; else if (i_ < 24u) po[104u + i_] = WR[{L['DEL']}u + {NMAX}u + i_]; }} }}") if DBG else ''
    _ou = (OUT_LANES_B if STORE_B else OUT_ORIG_B)
    RES_FLUSH = (f"for (uint i = lw; i < 96u; i += 32u) res48[pair * 96u + i] = WR[{L['RESA']}u + i];\n  ") if (RES48_B and not CISFUSE_B) else ''
    if CISFUSE_B: RES_FLUSH = CISFUSE_TXT.replace('@RESA@', str(L['RESA'])).replace('@RESC@', str(L['RESC']))
    RES_ACC_LINE = (f"for (uint i = lw; i < 96u; i += 32u) WR[{L['RESA']}u + i] = WR[{L['RESC']}u + (ph * 2u + i / 48u) * 48u + i % 48u];\n        ") if RES48_B else ''
    RESP = ("device float* res48 [[buffer(15)]], " + ("constant float* cfc [[buffer(16)]], " if CISFUSE_B else "")) if RES48_B else ""
    return (f"""
kernel void fjn_q{Q}_t{T}(device const uint4* progP [[buffer(0)]], device const uint4* progD [[buffer(1)]], device const uint4* hdr [[buffer(2)]], device const ushort* sprog [[buffer(3)]],
    device const uint* segs [[buffer(4)]], constant float* yv [[buffer(5)]], device const uint2* pairs [[buffer(6)]], device const float* th0 [[buffer(7)]], device float* pout [[buffer(8)]],
    constant float4& tf [[buffer(9)]], constant uint4& rp [[buffer(14)]], device float* EB [[buffer(10)]], constant float* Sm [[buffer(11)]], device const uint* plist [[buffer(12)]], constant uint4& wq [[buffer(13)]],
    {RESP}threadgroup float* TGM [[threadgroup(0)]], uint tg [[threadgroup_position_in_grid]], uint l [[thread_index_in_threadgroup]], uint nsg [[simdgroups_per_threadgroup]]) {{
  const uint sid = l >> 5, lw = l & 31u, hf = lw >> 4, lq = lw & 15u, p0 = lq * 3u;
  const uint gi = tg * nsg + sid; if (gi >= wq.x) return;
  const uint per = wq.y; const bool closed = wq.z != 0u; const uint maxit = wq.w;
  const float tau = tf.x, flr = tf.y;
  const uint pair = plist[gi]; const uint2 rc = pairs[pair]; device float* po = pout + pair * 128u;
  const uint ovl = {OVLX};
  threadgroup float* SBm = TGM + sid * ({SLOTC}u + ovl);
  threadgroup float* TSH = SBm + hf * {CD * 48}u;
  threadgroup float* WR = SBm + {2 * CD * 48}u;
  threadgroup float* R = WR + {WS}u + hf * per;
  threadgroup float* CB = WR + {WS}u; (void)CB;
  const uint ebase = gi * {4 * HS(Q)}u;
  const uint nx = hdr[2u * rc.x].w >> 24, ny = hdr[2u * rc.y].w >> 24;
  const uint nk = closed ? {KW}u : {NW}u; const uint n = nx + ny + nk; const uint nn = n * (n + 1u) / 2u;
  // ---- init state
  for (uint j = lw; j < {NPK}u; j += 32u) {{ WR[{L['THA']}u + j] = (j < nx) ? th0[pair * 64u + j] : 0.0f; WR[{L['THA']}u + {NPK}u + j] = (j < ny) ? th0[pair * 64u + 32u + j] : 0.0f; }}
  if (lw < {NW}u) WR[{L['UA']}u + lw] = Sm[48u * 9u + lw];
  if (lw == 9u) WR[{L['UA']}u + 9u] = 0.0f;
  simdgroup_barrier(mem_flags::mem_threadgroup);
  // candidates = accepted (both halves), t grid from u
  for (uint j = lw; j < {2 * NPK}u; j += 32u) {{ const float v = WR[{L['THA']}u + j]; WR[{L['THC']}u + j] = v; WR[{L['THC']}u + {2 * NPK}u + j] = v; }}
  if (lw < 16u) {{ const float v = WR[{L['UA']}u + lw]; WR[{L['UC']}u + lw] = v; WR[{L['UC']}u + 16u + lw] = v; }}
  if (closed && lw < 2u) WR[{L['UC']}u + lw * 16u + 9u] = (float)lw * 0.125f;
  simdgroup_barrier(mem_flags::mem_threadgroup);
  if (lw == 0u) {{ float U_ = 0; for (uint k = 0; k < 9u; k++) U_ += WR[{L['UC']}u + k] * WR[{L['UC']}u + k]; WR[{L['SC']}u + 3u] = U_; }}
  simdgroup_barrier(mem_flags::mem_threadgroup);
  for (uint hh = 0; hh < 2u; hh++) {{
    if (hf == hh) {{ const float U_ = WR[{L['SC']}u + 3u]; const float t0c = WR[{L['UC']}u + hh * 16u + 9u];
      for (uint v = 0; v < 3u; v++) {{ float s_ = 0; for (uint k = 0; k < 9u; k++) {{ const float uk = WR[{L['UC']}u + hh * 16u + k]; s_ = fma(Sm[(p0 + v) * 9u + k], uk * uk, s_); }} TSH[p0 + v] = t0c + s_ / U_; }} }} }}
  simdgroup_barrier(mem_flags::mem_threadgroup);
  // SC: 0 F, 1 lam, 2 nu(unused), 3 UU (accepted u norm^2), 4 mode, 5 it, 6 pair index, 7 Lcur, 8 pred0, 9 pred1, 10 stalled, 11 curh, 12 trials, 13 fail0, 14 fail1, 15 dF last, 16 plain passes, 17 dual passes
  if (lw < 2u) WR[{L['RSSA']}u + lw] = INFINITY;
  if (lw == 0u) {{ WR[{L['SC']}u + 0u] = INFINITY; WR[{L['SC']}u + 1u] = 1e-3f; WR[{L['SC']}u + 4u] = 0.0f; WR[{L['SC']}u + 5u] = 0.0f; WR[{L['SC']}u + 10u] = 0.0f; WR[{L['SC']}u + 11u] = 0.0f; WR[{L['SC']}u + 12u] = 0.0f; WR[{L['SC']}u + 13u] = 0.0f; WR[{L['SC']}u + 14u] = 0.0f; WR[{L['SC']}u + 18u] = 0.0f; WR[{L['SC']}u + 16u] = 0.0f; WR[{L['SC']}u + 17u] = 0.0f; }}
  simdgroup_barrier(mem_flags::mem_threadgroup);
  uint mode = 0u, curh = 0u, it = 0u, rungp = 0u, scanr = 0u, bestk = 0u; bool done = false, init = true; float bestF = INFINITY;
  float Fcur = INFINITY;
  while (!done) {{
    if (mode == 1u) {{
      for (uint rdl = 0; rdl < rp.y; rdl++) {{
      for (uint e = lw; e < nn; e += 32u) WR[{L['HP']}u + e] = 0.0f;
      for (uint i_ = lw; i_ < n; i_ += 32u) WR[{L['GV']}u + i_] = 0.0f;{ZDEL}
      simdgroup_barrier(mem_flags::mem_threadgroup);
{dual}
      }}
      // Gram and gradient were accumulated per side inside the dual pass
      simdgroup_barrier(mem_flags::mem_threadgroup);
      {SCLF}
      {DBGG}
      if (lw == 0u) {{ WR[{L['SC']}u + 6u] = 0.0f; WR[{L['SC']}u + 7u] = WR[{L['SC']}u + 1u]; WR[{L['SC']}u + 17u] += 1.0f; }}
      simdgroup_barrier(mem_flags::mem_threadgroup);
      for (uint br = 0; br < (uint)tf.z; br++) simdgroup_barrier(mem_flags::mem_threadgroup);
      rungp = 0u; mode = 2u;
    }}
    if (mode == 2u) {{
      if (lw == 0u) WR[{L['SC']}u + 18u] += 1.0f;
      // ---- LM ladder pair: half 0 solves lambda_{{2p}}, half 1 lambda_{{2p+1}}; candidates set; t grids
      const float Lc = WR[{L['SC']}u + 7u]; const float lamh = hf ? Lc * exp2((float)(2u * rungp + 1u)) : Lc;
      threadgroup float* LPh = {LPHX}; threadgroup float* XBh = WR + {L['XB']}u + hf * {NMAX}u;
      bool okc = true;
      for (uint rsl = 0; rsl < rp.w; rsl++) {{
      okc = true;
      const bool rsd = {RS_ON} && n <= 16u;
      if (rsd) {{
        simdgroup_barrier(mem_flags::mem_threadgroup);
        if (n <= 10u) {{
{GJ8}
        }} else if (n <= 12u) {{
{GJ12}
        }} else if (n <= 14u) {{
{GJ14}
        }} else {{
{GJ16}
        }}
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }} else {{
      for (uint e = lq; e < nn; e += 16u) {{
        uint i_ = (uint)((sqrt(8.0f * (float)e + 1.0f) - 1.0f) * 0.5f); while (i_ * (i_ + 1u) / 2u > e) i_--; while ((i_ + 1u) * (i_ + 2u) / 2u <= e) i_++;
        const uint j_ = e - i_ * (i_ + 1u) / 2u; const float si = WR[{L['SCL']}u + i_], sj = WR[{L['SCL']}u + j_];
        float a_ = si * sj * WR[{L['HP']}u + e];
        if (i_ == j_) a_ += (si > 0.0f) ? lamh : 1.0f;
        LPh[e] = a_;
      }}
      for (uint i_ = lq; i_ < n; i_ += 16u) XBh[i_] = WR[{L['SCL']}u + i_] * WR[{L['GV']}u + i_];
      simdgroup_barrier(mem_flags::mem_threadgroup);
      okc = true;
      if ({ABL_SOLVE}) {{   // timing ablation only (FJN_ABL=1): Jacobi step in place of the Cholesky solves
        for (uint i_ = lq; i_ < n; i_ += 16u) XBh[i_] = XBh[i_] / max(LPh[i_ * (i_ + 1u) / 2u + i_], 1e-6f);
        simdgroup_barrier(mem_flags::mem_threadgroup);
      }} else {{
{SOLVE}
      }}
      }}
      }}
      // pred = x . b + lam |x_live|^2 ; delta = scl * x
      float pa = 0;
      for (uint i_ = lq; i_ < n; i_ += 16u) {{ const float x_ = XBh[i_]; const float s_ = WR[{L['SCL']}u + i_]; pa += x_ * s_ * WR[{L['GV']}u + i_] + ((s_ > 0.0f) ? lamh * x_ * x_ : 0.0f); }}
      pa = hsum16(pa);
      for (uint i_ = lq; i_ < n; i_ += 16u) WR[{L['DEL']}u + hf * {NMAX}u + i_] = WR[{L['SCL']}u + i_] * XBh[i_];
      if (lq == 0u) {{ WR[{L['SC']}u + 8u + hf] = pa; WR[{L['SC']}u + 13u + hf] = okc ? 0.0f : 1.0f; }}
      simdgroup_barrier(mem_flags::mem_threadgroup);
      {DBGD}
      // candidates of this half
      const bool fail = !okc;
      for (uint j = lq; j < {NPK}u; j += 16u) {{
        const float dx = (j < nx && !fail) ? WR[{L['DEL']}u + hf * {NMAX}u + j] : 0.0f; const float dy = (j < ny && !fail) ? WR[{L['DEL']}u + hf * {NMAX}u + nx + j] : 0.0f;
        WR[{L['THC']}u + (hf * 2u + 0u) * {NPK}u + j] = WR[{L['THA']}u + j] + dx; WR[{L['THC']}u + (hf * 2u + 1u) * {NPK}u + j] = WR[{L['THA']}u + {NPK}u + j] + dy; }}
      if (lq < 10u) {{ float d_ = 0; if (lq < nk && !fail) d_ = WR[{L['DEL']}u + hf * {NMAX}u + nx + ny + lq];
        float u_ = WR[{L['UA']}u + lq] + d_; if (lq == 9u) u_ = min(max(u_, 0.0f), 1.0f); WR[{L['UC']}u + hf * 16u + lq] = u_; }}
      simdgroup_barrier(mem_flags::mem_threadgroup);
      {{ float U_ = 0; for (uint k = 0; k < 9u; k++) {{ const float uk = WR[{L['UC']}u + hf * 16u + k]; U_ += uk * uk; }} const float t0c = WR[{L['UC']}u + hf * 16u + 9u];
        for (uint v = 0; v < 3u; v++) {{ float s_ = 0; for (uint k = 0; k < 9u; k++) {{ const float uk = WR[{L['UC']}u + hf * 16u + k]; s_ = fma(Sm[(p0 + v) * 9u + k], uk * uk, s_); }} TSH[p0 + v] = t0c + s_ / U_; }} }}
      simdgroup_barrier(mem_flags::mem_threadgroup);
      mode = 0u;
    }}
    if (mode == 0u) {{
{RSSDECL}      for (uint rpl = 0; rpl < rp.x; rpl++) {{
{plain}
      }}
      // ---- decision (uniform: every lane reads the same TG scalars)
      if (lw == 0u) WR[{L['SC']}u + 16u] += 1.0f;
      simdgroup_barrier(mem_flags::mem_threadgroup);
      {("const float F0 = " + SIMDRED_F0 + "; const float F1 = " + SIMDRED_F1 + ";") if SIMDRED else ("const float F0 = WR[" + str(L['RSSC']) + "u + 0u] + WR[" + str(L['RSSC']) + "u + 1u]; const float F1 = WR[" + str(L['RSSC']) + "u + 2u] + WR[" + str(L['RSSC']) + "u + 3u];")}
      int pick = -1;
      if (init) {{
        if (closed && scanr < 4u) {{
          // t0 start scan (closed strokes): 4 plain passes with t0 = (2r+h)/8 (identity warp, scan theta), then one pass with both halves at the best t0
          if (isfinite(F0) && F0 < bestF) {{ bestF = F0; bestk = 2u * scanr; }}
          if (isfinite(F1) && F1 < bestF) {{ bestF = F1; bestk = 2u * scanr + 1u; }}
          scanr++;
          if (lw < 2u) {{ const float t0n = (scanr < 4u) ? (float)(2u * scanr + lw) * 0.125f : (float)bestk * 0.125f; WR[{L['UC']}u + lw * 16u + 9u] = t0n; }}
          simdgroup_barrier(mem_flags::mem_threadgroup);
          {{ float U_ = 0; for (uint k = 0; k < 9u; k++) {{ const float uk = WR[{L['UC']}u + hf * 16u + k]; U_ += uk * uk; }} const float t0c = WR[{L['UC']}u + hf * 16u + 9u];
            for (uint v = 0; v < 3u; v++) {{ float s_ = 0; for (uint k = 0; k < 9u; k++) {{ const float uk = WR[{L['UC']}u + hf * 16u + k]; s_ = fma(Sm[(p0 + v) * 9u + k], uk * uk, s_); }} TSH[p0 + v] = t0c + s_ / U_; }} }}
          simdgroup_barrier(mem_flags::mem_threadgroup);
          if (scanr == 4u && !isfinite(bestF)) done = true;
          pick = -2;
        }} else if (isfinite(F0)) pick = 0; else {{ done = true; }}
      }} else {{
        const bool f0 = WR[{L['SC']}u + 13u] == 0.0f, f1 = WR[{L['SC']}u + 14u] == 0.0f;
        if (f0 && F0 < Fcur) pick = 0; else if (f1 && F1 < Fcur) pick = 1;
      }}
      if (done) break;
      if (pick == -2) {{ }} else if (pick >= 0) {{
        const uint ph = (uint)pick; const float Fn = (ph == 0u) ? F0 : F1;
        // accept: copy candidate ph -> accepted
        for (uint j = lw; j < {2 * NPK}u; j += 32u) WR[{L['THA']}u + j] = WR[{L['THC']}u + ph * {2 * NPK}u + j];
        for (uint j = lw; j < {2 * QP}u; j += 32u) WR[{L['ACA']}u + j] = WR[{L['ACC']}u + ph * {2 * QP}u + j];
        if (lw < 10u) WR[{L['UA']}u + lw] = WR[{L['UC']}u + ph * 16u + lw];
        for (uint i = lw; i < 48u; i += 32u) WR[{L['TA']}u + i] = TGM[sid * ({SLOTC}u + ovl) + ph * {CD * 48}u + i];
        if (lw < 2u) WR[{L['RSSA']}u + lw] = WR[{L['RSSC']}u + ph * 2u + lw];
        {RES_ACC_LINE}simdgroup_barrier(mem_flags::mem_threadgroup);
        if (lw == 0u) {{
          float U_ = 0; for (uint k = 0; k < 9u; k++) U_ += WR[{L['UA']}u + k] * WR[{L['UA']}u + k]; WR[{L['SC']}u + 3u] = U_;
          if (!init) {{
            const float lamacc = (ph == 0u) ? WR[{L['SC']}u + 7u] : WR[{L['SC']}u + 7u] * exp2((float)(2u * rungp + 1u));
            const float pred = WR[{L['SC']}u + 8u + ph];
            const float rho = clamp((Fcur - Fn) / max(pred, 1e-30f), -1e3f, 1e3f);
            const float c_ = 2.0f * rho - 1.0f;
            WR[{L['SC']}u + 1u] = lamacc * max(1.0f / 3.0f, 1.0f - c_ * c_ * c_);
            WR[{L['SC']}u + 15u] = Fcur - Fn;
          }}
          WR[{L['SC']}u + 0u] = Fn; WR[{L['SC']}u + 11u] = (float)ph;
          if (!{DBGF} && !init && it < 13u) {{ po[114u + it] = Fn; }}
          if (init) po[127u] = Fn;
        }}
        simdgroup_barrier(mem_flags::mem_threadgroup);
        curh = ph; Fcur = WR[{L['SC']}u + 0u];
        if (init) {{ init = false; mode = 1u; }}
        else {{
          it++;
          const float dF = WR[{L['SC']}u + 15u];
          if (dF < flr || it >= maxit) done = true; else mode = 1u;
        }}
      }} else {{
        // rejected pair: advance the ladder
        if (lw == 0u) {{ const float Lc = WR[{L['SC']}u + 7u]; WR[{L['SC']}u + 7u] = {'Lc * exp2((float)(4u * rungp + 3u))' if IM2_LAD_B else 'Lc * exp2((float)(2u * rungp + 1u)) * exp2((float)(2u * rungp + 2u))'}; }}
        simdgroup_barrier(mem_flags::mem_threadgroup);
        rungp++;
        if (lw == 0u) WR[{L['SC']}u + 12u] += 2.0f;
        if (rungp >= 6u) {{ done = true; if (lw == 0u) WR[{L['SC']}u + 10u] = 1.0f; }} else mode = 2u;
      }}
    }}
  }}
  // ---- outputs: pout[pair*128]: 0 Fx, 1 Fy, 2 stalled, 3 it, 4 trials, 5 dual passes, 6 plain passes, 8..39 thx, 40..71 thy, 72..87 ax, 88..103 ay, 104..113 u
@@OUT_BLOCK@@
  if (!{DBGF}) for (uint j = lw; j < 32u; j += 32u) {{ po[8u + j] = (j < {NPK}u) ? WR[{L['THA']}u + j] : 0.0f; po[40u + j] = (j < {NPK}u) ? WR[{L['THA']}u + {NPK}u + j] : 0.0f; }}
  if (!{DBGF}) for (uint j = lw; j < 16u; j += 32u) {{ po[72u + j] = (j < {QP}u) ? WR[{L['ACA']}u + j] : 0.0f; po[88u + j] = (j < {QP}u) ? WR[{L['ACA']}u + {QP}u + j] : 0.0f; }}
  if (lw < 10u) po[104u + lw] = WR[{L['UA']}u + lw];
{RES_FLUSH}}}
""").replace("@@OUT_BLOCK@@", eval(f'f"""{_ou}"""'))

def header():
    return "\n".join(["#include <metal_stdlib>", "using namespace metal;", CTAB,
        "inline float hsum16(float x) { x = quad_sum(x); x += simd_shuffle_xor(x, 8); x += simd_shuffle_xor(x, 4); return x; }",
        "inline bool bad_half(bool b) { uint m = (uint)b; m |= simd_shuffle_xor(m, 8); m |= simd_shuffle_xor(m, 4); m |= simd_shuffle_xor(m, 2); m |= simd_shuffle_xor(m, 1); return m != 0; }"])

if __name__ == "__main__":
    QS = [int(x) for x in sys.argv[1].split(',')]; TS = [int(x) for x in sys.argv[2].split(',')]
    if os.environ.get('FJN_TMAX_H'):
        with open(os.environ['FJN_TMAX_H'], 'w') as f: f.write(f"#define FJN_TMAX {max(TS)}\n")
    print(header())
    ONLY = None
    if os.environ.get('FJN_ONLY'): ONLY = set(tuple(int(v) for v in c.split(':')) for c in os.environ['FJN_ONLY'].split(','))
    for Q in QS:
        for T in TS:
            if ONLY is None or (Q, T) in ONLY: print(kern(Q, T))
