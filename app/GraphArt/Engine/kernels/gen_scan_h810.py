import sys, os
H805 = os.environ.get('H805') == '1'
H806 = os.environ.get('H806') == '1'
FMAX = os.environ.get('FMAX', '1') == '1'
FSQ, FDOT, FAX = ('a', 'a', '1') if FMAX else (os.environ.get('FSQ', 'p'), os.environ.get('FDOT', 'p'), os.environ.get('FAX', 'p'))
B1 = os.environ.get('B1') == '1'
BU = os.environ.get('BU', '0')
H622 = os.environ.get('H622', '1') != '0'
H333 = os.environ.get('H333') == '1'
if H806 and not H622: sys.exit('H806 requires H622=1 (the all-bad guard)')
if (H805 or H806) and H333: sys.exit('H805/H806 are not supported with H333=1')
UN = {7: 'sin({x})', 8: 'cos({x})', 9: 'exp({x})', 10: 'log(fabs({x}))', 11: '({x})*({x})', 12: '({x})*({x})*({x})', 13: 'sqrt(fabs({x}))', 14: 'copysign(pow(fabs({x}), 1.0f/3.0f), {x})', 15: '-({x})'}

def gen(kind, G, QS, sfx):
    NP = 24 if kind == '24' else 48
    SS = 16 if kind == '24' else 4
    LN = NP // G; NT = SS * LN; SL = SS * NP
    SH = {384: 1, 192: 2}[SL]
    V = range(G)
    def dot(a, b, kind_='S'):
        f_ = FSQ if kind_ == 'S' else FDOT
        if f_ == 'u': return " + ".join(f"simd_shuffle_xor({a(v)} * {b(v)}, 0)" for v in V)
        if f_ == 'p': return " + ".join(f"{a(v)} * {b(v)}" for v in V)
        order = list(range(G)) if f_ in 'ao' else list(range(G))[::-1]
        s_ = f"{a(order[0])} * {b(order[0])}"
        for v in order[1:]: s_ = f"fma({a(v)}, {b(v)}, {s_})"
        return s_
    def each(f): return " ".join(f(v) for v in V)
    def fetchB(w, X, nacc=False):
        if B1:
            cl = f"{w} >> 29" if nacc else f"({w} >> 29) & 3u"
            return (f"    {{ const uint cl_ = {cl}, of_ = {w} & 0x1FFFFFFFu; const uint b_ = (of_ >> (cl_ == 0u ? {SH}u : 0u)) + (cl_ == 0u ? l : cl_ == 1u ? sb : cl_ == 2u ? p0 : 0u - DIFF); const uint s_ = cl_ == 0u ? {NT}u : cl_ == 2u ? 1u : 0u;\n"
                    f"      " + each(lambda v: f"{X}{v} = NS[b_ + s_ * {v}u];") + " }\n")
        return (f"    {{ const uint cl_ = ({w} >> 29) & 3u, of_ = {w} & 0x1FFFFFFFu; const uint b_ = cl_ == 0u ? (of_ >> {SH}u) + l : cl_ == 1u ? of_ + sb : cl_ == 2u ? of_ + p0 : of_ - DIFF; const uint s_ = cl_ == 0u ? {NT}u : cl_ == 2u ? 1u : 0u;\n"
                f"      " + each(lambda v: f"{X}{v} = NS[b_ + s_ * {v}u];") + " }\n")
    def fetchA(w, X):
        return (f"    if ({w} >> 31) {{ " + each(lambda v: f"{X}{v} = acc{v};") + " }\n    else" + fetchB(w, X, True))
    def step(idx):
        s = f"    {{ const uint4 w = prog[h.x + {idx}];\n" if BU != '2' else f"    {{ const uint4 w0_ = prog[h.x + {idx}]; const uint4 w = uint4(simd_broadcast_first(w0_.x), simd_broadcast_first(w0_.y), simd_broadcast_first(w0_.z), simd_broadcast_first(w0_.w));\n"
        s += "    float " + ", ".join([f"a{v}" for v in V] + [f"b{v}" for v in V]) + ";\n"
        s += fetchA('w.z', 'a') + fetchB('w.w', 'b')
        s += "    if (w.x == 0u) { " + each(lambda v: f"acc{v} = a{v} * b{v};") + " }\n"
        s += "    else if (w.x < 4u) { const float4 c = ctab[w.x]; " + each(lambda v: f"acc{v} = c.y * fma(c.z, b{v}, a{v});") + " }\n"
        s += "    else { switch (w.x & 255u) { case 4: " + each(lambda v: f"acc{v} = a{v} / b{v};") + " break; case 5: " + each(lambda v: f"acc{v} = b{v} / a{v};") + " break;\n"
        for c, f in UN.items(): s += f"      case {c}: " + each(lambda v: f"{{ const float m_ = a{v} * b{v}; acc{v} = {f.format(x='m_')}; }}") + " break;\n"
        s += "      case 16: case 17: case 18: case 19: { const uint cw_ = ((w.x >> 8) & 0x7FFFu) | (((w.x >> 23) & 3u) << 29); float " + ", ".join(f"c{v}" for v in V) + ";\n" + fetchB('cw_', 'c') + "      if ((w.x & 255u) == 16u) { " + each(lambda v: f"acc{v} = (a{v} * b{v}) * c{v};") + " } else { const float4 k_ = ctab[(w.x & 255u) - 16u]; " + each(lambda v: f"acc{v} = k_.y * fma(k_.z, c{v}, a{v} * b{v});") + " } } break;\n"
        s += "      default: " + each(lambda v: f"acc{v} = a{v};") + " break; } }\n"
        s += f"    if (w.y != 0xFFFFFFFFu) {{ const uint d_ = w.y >> {SH}u; " + each(lambda v: f"NS[d_ + l + {NT*v}u] = acc{v};") + " } }\n"
        return s
    def runseg(i0, i1):
        return f"  for (uint i = {i0}; i < {i1}; i++) {{\n" + step('i') + "  }\n"
    hs = [LN >> i for i in range(1, 6) if LN >> i >= 1]
    out = []
    out.append(f"inline float hsum{sfx}(float x) {{ x = quad_sum(x); " + " ".join(f"x += simd_shuffle_xor(x, {k});" for k in hs if k >= 4) + " return x; }")
    out.append(f"inline bool badh{sfx}(bool b) {{ uint m = (uint)b; " + " ".join(f"m |= simd_shuffle_xor(m, {k});" for k in hs) + " return m != 0; }")
    for Q in QS:
        L = []
        name = ('rcscan24_q' if kind == '24' else 'rcscan4x_q') + str(Q)
        HSIG = ' device float* hand [[buffer(16)]], constant uint& hbase [[buffer(17)]],' if (H333 and kind != '24') else ''
        L.append(f'''kernel void {name}({'constant' if BU == '1' else 'device const'} uint4* prog [[buffer(0)]], device const uint4* hdr [[buffer(1)]], device const ushort* sprog [[buffer(2)]], device const float* bred [[buffer(3)]],
                   constant float* sob [[buffer(4)]], constant float* yv [[buffer(5)]], device float* out [[buffer(6)]], device const uint* treelist [[buffer(7)]],
                   device const uint* segs [[buffer(8)]], constant uint& nslot [[buffer(9)]], constant uint& nsc [[buffer(10)]], constant float4* ctab [[buffer(11)]], device float* outa [[buffer(12)]], constant float& tau [[buffer(13)]], threadgroup float* NS [[threadgroup(0)]],{HSIG}
                   uint tg [[threadgroup_position_in_grid]], uint l [[thread_index_in_threadgroup]], uint ln [[thread_index_in_simdgroup]]) {{
  const uint tree = treelist[tg]; const uint4 h = hdr[tree];
  const int nth = h.w & 255; const uint thoff = h.w >> 8; const int nS = h.y >> 16;
  threadgroup float* SV = NS + nslot * {SL}u; const uint TT = nslot * {SL}u + {SS}u * nsc;
  const uint DIFF = nslot * {768 - SL}u + {16 - SS}u * nsc;   // host class-3 offsets are TT(768 layout)+48/49
  device const uint* sg = segs + tree * 36;
  const uint sl = l / {LN}u, lq = l % {LN}u; const uint sb = nslot * {SL}u + sl * nsc; const uint p0 = TT + lq * {G}u;''')
        if kind == '24':
            L.append(f"  for (uint i = l; i < 50u; i += {NT}u) NS[TT + i] = i < 24u ? (float)(2u * i) * (1.0f / 47) : (i == 48u ? 1.0f : (i == 49u ? 2.0f : 0.0f));")
            L.append(f"  for (uint ix = l; ix < 16u * (uint)nS; ix += {NT}u) {{ const uint st = ix % 16u, i = ix / 16u; const ushort w = sprog[h.z + i]; const uint op = w & 63, a = (w >> 6) & 31; float x = 1.0f;")
        else:
            L.append("  threadgroup uint selS[4];")
            L.append("  if (l < 32u) { const float vl = (l < 16u) ? outa[tree * 16u + l] : INFINITY; uint r = 0;")
            L.append("    for (ushort j = 0; j < 16; j++) { const float vj = simd_shuffle(vl, j); r += (vj < vl || (vj == vl && (uint)j < l)) ? 1u : 0u; }")
            L.append("    if (l < 16u && r < 4u) selS[r] = l; }")
            L.append("  if (l < 16u) out[tree * 16u + l] = INFINITY;")
            L.append("  threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);")
            L.append(f"  for (uint i = l; i < 50u; i += {NT}u) NS[TT + i] = i < 48u ? i * (1.0f / 47) : (i == 48u ? 1.0f : 2.0f);")
            L.append(f"  for (uint ix = l; ix < 4u * (uint)nS; ix += {NT}u) {{ const uint st = selS[ix % 4u], i = ix / 4u; const ushort w = sprog[h.z + i]; const uint op = w & 63, a = (w >> 6) & 31; float x = 1.0f;")
        L.append("    if (op == 1u) { float u = (nth == 1) ? (float)st / 15.0f : sob[st*32 + a]; float bf = bred[thoff + a]; x = (bf - 12.566370614359172f) + (2.0f*12.566370614359172f)*u; }")
        L.append(f"    SV[(ix % {SS}u) * nsc + i] = x; }}")
        L.append("  threadgroup_barrier(mem_flags::mem_threadgroup);")
        L.append("  float " + ", ".join(f"acc{v} = 0" for v in V) + ";")
        if H806:
            L.append("  float " + ", ".join(f"c{q}_{v}" for q in range(Q) for v in V) + ";")
            L.append('  bool bad = false;')
            for q in range(Q):
                if q > 0: L.append("  if (!dead) {")
                L.append(runseg(f"sg[{q}]", f"sg[{q+1}]"))
                L.append(f"  {{ const uint cw = sg[{18+q}];\n" + fetchA('cw', f'c{q}_') + "  }")
                L.append("  bad = bad || " + " || ".join(f"!isfinite(c{q}_{v})" for v in V) + ";")
                if q <= Q - 2 and q == 0: L.append(f"  bool dead = simd_all(badh{sfx}(bad));")
                elif q <= Q - 2: L.append(f"  dead = simd_all(badh{sfx}(bad));")
                if q > 0: L.append("  }")
        else:
          for q in range(Q):
            L.append(runseg(f"sg[{q}]", f"sg[{q+1}]"))
            L.append("  float " + ", ".join(f"c{q}_{v}" for v in V) + ";")
            L.append(f"  {{ const uint cw = sg[{18+q}];\n" + fetchA('cw', f'c{q}_') + "  }")
        if kind == '24': L.append("  float " + ", ".join(f"c{Q}_{v} = yv[2u * (lq * {G}u + {v}u)]" for v in V) + ";")
        else: L.append("  float " + ", ".join(f"c{Q}_{v} = yv[lq * {G} + {v}]" for v in V) + ";")
        TRK = H333 and kind != '24'
        if not H806: L.append('  bool bad = false;')
        if TRK: L.append("  " + " ".join(f"float t{k} = (lq == {k}u) ? 1.0f : 0.0f;" for k in range(Q)) + " float tY = 0.0f;")
        tn = [f"t{k}" for k in range(Q)] + ["tY"]
        if not H806:
            for q in range(Q): L.append("  bad = bad || " + " || ".join(f"!isfinite(c{q}_{v})" for v in V) + ";")
        if not TRK and H622: L.append(f"  const bool badh = badh{sfx}(bad); float rss = INFINITY;\n  if (!simd_all(badh)) {{")
        L.append("  " + " ".join(f"const float on{k} = sqrt(hsum{sfx}(" + dot(lambda v: f"c{k}_{v}", lambda v: f"c{k}_{v}") + "));" for k in range(1, Q)))
        for k in range(Q):
            L.append("  {")
            if k == 0:
                L.append(f"   const float nn = hsum{sfx}(" + dot(lambda v: f"c{k}_{v}", lambda v: f"c{k}_{v}") + "); const float inv = nn > 0 ? rsqrt(nn) : 0.0f;")
            else:
                L.append(f"   const float nn = hsum{sfx}(" + dot(lambda v: f"c{k}_{v}", lambda v: f"c{k}_{v}") + f"); const float inv = (sqrt(nn) > tau * on{k}) ? rsqrt(nn) : 0.0f;")
            L.append("   " + (f"t{k} *= inv; " if TRK else "") + (each(lambda v: f"const float o{k}_{v} = c{k}_{v};") + " " if (FDOT == 'o' or FAX in '23') else "") + each(lambda v: f"c{k}_{v} *= inv;"))
            if H805: L.append("   if (!simd_all(inv == 0.0f)) {")
            for j in range(k+1, Q+1):
                dk = (lambda v: f"o{k}_{v}") if FDOT == 'o' else (lambda v: f"c{k}_{v}")
                axp = {'u': lambda v: f"c{j}_{v} -= simd_shuffle_xor(r * c{k}_{v}, 0);", 'p': lambda v: f"c{j}_{v} -= r * c{k}_{v};", '1': lambda v: f"c{j}_{v} = fma(-r, c{k}_{v}, c{j}_{v});", '2': lambda v: f"c{j}_{v} = fma(-inv, r * o{k}_{v}, c{j}_{v});", '3': lambda v: f"c{j}_{v} = fma(-(r * inv), o{k}_{v}, c{j}_{v});"}[FAX]
                L.append(f"   {{ const float r = hsum{sfx}(" + dot(dk, lambda v: f"c{j}_{v}", "R") + ")" + (" * inv" if FDOT == 'o' else "") + "; " + (f"{tn[j]} = fma(-r, t{k}, {tn[j]}); " if TRK else "") + each(axp) + " }")
            if H805: L.append("   }")
            L.append("  }")
        if TRK or not H622: L.append(f"  const float rss = hsum{sfx}(" + dot(lambda v: f"c{Q}_{v}", lambda v: f"c{Q}_{v}") + f"); const bool badh = badh{sfx}(bad);")
        else: L.append(f"  rss = hsum{sfx}(" + dot(lambda v: f"c{Q}_{v}", lambda v: f"c{Q}_{v}") + ");\n  }")
        if kind == '24': L.append('  if (lq == 0) outa[tree * 16 + sl] = (badh || !isfinite(rss)) ? INFINITY : rss;')
        else: L.append('  if (lq == 0) out[tree * 16 + selS[sl]] = (badh || !isfinite(rss)) ? INFINITY : rss;')
        if TRK:
            HSTR = (Q + 1) * 48 + 20
            L.append('  threadgroup float rssS[4]; if (lq == 0) rssS[sl] = (badh || !isfinite(rss)) ? INFINITY : rss;')
            L.append('  threadgroup_barrier(mem_flags::mem_threadgroup);')
            L.append('  int wsl = -1; { float bv = INFINITY; uint bi = 99u; for (uint k = 0; k < 4u; k++) { const float v_ = rssS[k]; const uint si = selS[k]; if (v_ < bv || (v_ == bv && isfinite(v_) && si < bi)) { bv = v_; bi = si; wsl = (int)k; } } }')
            L.append(f'  if ((int)sl == wsl) {{ device float* H_ = hand + hbase + tg * {HSTR}u; ' + ' '.join(f'H_[{q * 48}u + lq * {G}u + {v}u] = c{q}_{v};' for q in range(Q + 1) for v in V) + f' if (lq < {Q}u) H_[{(Q + 1) * 48}u + lq] = -tY; if (lq == 0u) H_[{(Q + 1) * 48 + 16}u] = rss; }}')

        L.append('}')
        out.append('\n'.join(L))
    return '\n'.join(out)

HI = list(range(6, 17)); LO = list(range(1, 6))
print('#include <metal_stdlib>\nusing namespace metal;')
print(gen('24', 3, HI, '24g3'))
print(gen('24', 6, LO, '24g6'))
print(gen('4x', 3, HI, '4xg3'))
print(gen('4x', 6, LO, '4xg6'))
