#!/usr/bin/env python3 -I
import re, subprocess, sys, os
mask = int(sys.argv[1]); here = os.path.dirname(os.path.abspath(__file__)) + '/'
t = subprocess.run([sys.executable, '-I', here + 'gen_scan_h810.py'], capture_output=True, text=True, check=True, env=os.environ).stdout
def sub(pat, rep, s):
    s2, n = re.subn(pat, rep, s)
    if n == 0: sys.exit('no match for ' + pat)
    return s2, n
if mask & 3:
    sg = bool(mask & 1)
    T = 'int' if sg else 'ushort'
    F = re.compile(r'const uint cl_ = ([^,]+), of_ = ([\w.]+) & 0x1FFFFFFFu; const uint b_ = \(of_ >> \(cl_ == 0u \? (\d)u : 0u\)\) \+ \(cl_ == 0u \? l : cl_ == 1u \? sb : cl_ == 2u \? p0 : 0u - DIFF\); const uint s_ = cl_ == 0u \? (\d+)u : cl_ == 2u \? 1u : 0u;')
    if sg:
        def rf(m): return (f'const int cl_ = (int)({m.group(1)}), of_ = (int)({m.group(2)} & 0x1FFFFFFFu); const int b_ = (of_ >> (cl_ == 0 ? {m.group(3)} : 0)) + (cl_ == 0 ? (int)l : cl_ == 1 ? (int)sb : cl_ == 2 ? (int)p0 : -(int)DIFF); const int s_ = cl_ == 0 ? {m.group(4)} : cl_ == 2 ? 1 : 0;')
        t, n = F.subn(rf, t); assert n > 0
        t, _ = sub(r'NS\[b_ \+ s_ \* (\d)u\]', r'NS[b_ + s_ * \1]', t)
        t, _ = sub(r'const uint d_ = w\.y >> (\d)u;', r'const int d_ = (int)(w.y >> \1u);', t)
        t, _ = sub(r'NS\[d_ \+ l \+ (\d+)u\]', r'NS[d_ + (int)l + \1]', t)
    else:
        def rf(m):
            if mask & 16:
                X = m.group(2); hi = f'(xh_.y >> 13)' if '& 3u' not in m.group(1) else f'((xh_.y >> 13) & 3)'
                return (f'const ushort2 xh_ = as_type<ushort2>({X}); const ushort cl_ = (ushort)({hi}), of_ = xh_.x; const ushort b_ = (ushort)((ushort)(of_ >> (cl_ == 0 ? {m.group(3)} : 0)) + (cl_ == 0 ? (ushort)l : cl_ == 1 ? (ushort)sb : cl_ == 2 ? (ushort)p0 : (ushort)(0u - DIFF))); const ushort s_ = cl_ == 0 ? (ushort){m.group(4)} : cl_ == 2 ? (ushort)1 : (ushort)0;')
            return (f'const ushort cl_ = (ushort)({m.group(1)}), of_ = (ushort)({m.group(2)} & 0x1FFFFFFFu); const ushort b_ = (ushort)((ushort)(of_ >> (cl_ == 0 ? {m.group(3)} : 0)) + (cl_ == 0 ? (ushort)l : cl_ == 1 ? (ushort)sb : cl_ == 2 ? (ushort)p0 : (ushort)(0u - DIFF))); const ushort s_ = cl_ == 0 ? (ushort){m.group(4)} : cl_ == 2 ? (ushort)1 : (ushort)0;')
        t, n = F.subn(rf, t); assert n > 0
        t, _ = sub(r'NS\[b_ \+ s_ \* (\d)u\]', r'NS[(ushort)(b_ + (ushort)(s_ * \1))]', t)
        t, _ = sub(r'const uint d_ = w\.y >> (\d)u;', r'const ushort d_ = (ushort)(w.y >> \1u);', t)
        t, _ = sub(r'NS\[d_ \+ l \+ (\d+)u\]', r'NS[(ushort)(d_ + (ushort)l + \1)]', t)
    assert 'const uint b_' not in t
if mask & 32:
    t, n = sub(r'const uint4 w = prog\[h\.x \+ ([a-z0-9]+)\];\n', r'const uint4 w = prog[h.x + \1]; const ushort2 xw_ = as_type<ushort2>(w.x);\n', t)
    t, _ = sub(r'if \(w\.x == 0u\) \{', 'if (xw_.x == 0) {', t)
    t, _ = sub(r'else if \(w\.x < 4u\) \{ const float4 c = ctab\[w\.x\];', 'else if (xw_.x < 4) { const float4 c = ctab[xw_.x];', t)
    t, _ = sub(r'switch \(w\.x & 255u\)', 'switch (xw_.x & 255)', t)
if mask & 64:
    t, n = sub(r'for \(uint i = (sg\[\d+\]); i < (sg\[\d+\]); i\+\+\)', r'for (ushort i = (ushort)\1; i < (ushort)\2; i++)', t)
if mask & 128:
    t, n = sub(r'(const uint sb = nslot \* \d+u \+ sl \* nsc; const uint p0 = TT \+ lq \* \dU?u;)'.replace('U?',''), r'\1 const ushort l16_ = (ushort)l, sb16_ = (ushort)sb, p016_ = (ushort)p0, df16_ = (ushort)(0u - DIFF);', t)
    for a, b in (('(ushort)l :', 'l16_ :'), ('(ushort)sb :', 'sb16_ :'), ('(ushort)p0 :', 'p016_ :'), ('(ushort)(0u - DIFF)))', 'df16_))'), ('(ushort)l + ', 'l16_ + ')):
        assert a in t, a
        t = t.replace(a, b)
if mask & 4:
    t, n = sub(r'for \(uint i = l; i < 50u; i \+= (\d+)u\) NS\[TT \+ i\] = i < 24u \? \(float\)\(2u \* i\) \* \(1\.0f / 47\) : \(i == 48u \? 1\.0f : \(i == 49u \? 2\.0f : 0\.0f\)\);',
               r'for (uint i = l; i < 26u; i += \1u) NS[TT + i] = i < 24u ? (float)(2u * i) * (1.0f / 47) : (i == 24u ? 1.0f : 2.0f);', t)
    t, n2 = sub(r'const uint DIFF = nslot \* 384u \+ 0u \* nsc;', 'const uint DIFF = nslot * 384u + 0u * nsc + 24u;', t)
    assert n == n2 == 17 - 0 or n == n2, (n, n2)
if mask & 8:
    t, n = re.subn(r'(  if \(!simd_all\(badh\)\) \{\n).*?(  rss = hsum)', r'\1\2', t, flags=re.S); assert n > 0, n
sys.stdout.write(t)
