
#pragma once

struct SelH33 {
  std::vector<uint64_t> hkey;
  std::vector<uint32_t> hpid, hpos;
  std::vector<uint64_t> tkey;
  std::vector<uint32_t> tent;
  std::vector<int> mstart, cur, fillA;
  struct PRec { uint32_t dst_slot, p, pos, j; };
  std::vector<PRec> plist;
};
static SelH33 g_sel_h33;

#ifdef H33_PROF
#define H33_T0 auto _h33t = std::chrono::steady_clock::now();
#define H33_LAP(name) { auto _h33n = std::chrono::steady_clock::now(); \
  fprintf(stderr, "h33 %s %.2f ms\n", name, std::chrono::duration<double, std::milli>(_h33n - _h33t).count()); _h33t = _h33n; }
#else
#define H33_T0
#define H33_LAP(name)
#endif

static void select_h33(evo::Pool& pool, evo::State& S, const std::vector<Meta>& pm, const std::vector<Meta>& cm,
                       const uint8_t* admit, std::vector<Meta>& nm, Times& T) {
  using namespace evo;
  H33_T0
  SelH33& H = g_sel_h33;
  const int NP = NPOP;
  const size_t HMAX = (size_t)2 * DECK + MIGR;
  if (H.hkey.size() != (size_t)NDEME * HMAX) {
    H.hkey.assign((size_t)NDEME * HMAX, 0); H.hpid.assign((size_t)NDEME * HMAX, 0); H.hpos.assign((size_t)NDEME * HMAX, 0);
  }
  H.mstart.assign(NDEME, 0); H.cur.assign(NDEME, 0); H.fillA.assign(NDEME, 0);

  pool.run(NDEME, [&](int d) { migrants(d, S.par.data(), S.migtab.data() + d * MIGR, S.sel[d].mtmp); });
  H33_LAP("migrants")
  pool.run(NDEME, [&](int d) {
    SelScratch& s = S.sel[d];
    s.pidx.clear();
    for (int i = d * DECK; i < (d + 1) * DECK; i++) s.pidx.push_back((uint32_t)i);
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      if (!admit[i]) continue;
      Ind& c = S.chi[i];
      float sum = 0.0f; for (int q = 0; q < NRES; q++) sum += c.res[q];
      c.err = sum; c.cx = (uint16_t)cplx(c.g);
      s.pidx.push_back((uint32_t)(NP + i));
    }
    H.mstart[d] = (int)s.pidx.size();
    int pd = (d + NDEME - 1) % NDEME;
    for (int k = 0; k < MIGR; k++) s.pidx.push_back(S.migtab[pd * MIGR + k]);
    const int m = (int)s.pidx.size();
    s.pe.resize(m); s.pc.resize(m); s.prank.resize(m); s.pcrowd.resize(m); s.surv.resize(m);
    for (int q = 0; q < m; q++) {
      uint32_t p = s.pidx[q];
      const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
      s.pe[q] = src.err; s.pc[q] = (uint16_t)src.cx;
    }
    rank_crowd(m, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    choose_survivors(m, s.pcrowd.data(), s.rs, s.surv.data(), s.tmp);
    const size_t base = (size_t)d * HMAX;
    for (int q = 0; q < m; q++) {
      const uint32_t pos = s.surv[q], p = s.pidx[pos];
      const Meta& mt = p < (uint32_t)NP ? pm[p] : cm[p - NP];
      H.hkey[base + q] = mt.kh; H.hpid[base + q] = p; H.hpos[base + q] = pos;
    }
  });

  H33_LAP("prep+rank")
  size_t cap = 1024; while (cap < 2 * (size_t)NP + 2) cap <<= 1;
  if (H.tkey.size() != cap) { H.tkey.assign(cap, 0); H.tent.assign(cap, 0); }
  else std::fill(H.tent.begin(), H.tent.end(), 0u);
  const uint64_t mask = (uint64_t)cap - 1;
  auto metaOf = [&](uint32_t p) -> const Meta& { return p < (uint32_t)NP ? pm[p] : cm[p - NP]; };
  auto place_key = [&](uint64_t h, uint32_t p, uint32_t f) -> bool {
    uint64_t i = h & mask;
    while (H.tent[i]) {
      if (H.tkey[i] == h && metaOf(H.hpid[H.tent[i] - 1]).key == metaOf(p).key) return false;
      i = (i + 1) & mask;
    }
    H.tkey[i] = h; H.tent[i] = f + 1; return true;
  };
  H.plist.clear(); H.plist.reserve((size_t)NP);
  auto record = [&](int dst, int slot, int j, uint32_t p, uint32_t pos) {
    H.plist.push_back(SelH33::PRec{(uint32_t)(dst * DECK + slot), p, pos, (uint32_t)j});
    if (p >= (uint32_t)NP) T.surv_child += 1;
  };
  const size_t PF = 8;
  for (int d = 0; d < NDEME; d++) {
    const size_t base = (size_t)d * HMAX;
    const int m = (int)S.sel[d].pidx.size();
    int k = 0, q = 0;
    for (; q < m && k < DECK; q++) {
      if (q + (int)PF < m) {
        const uint64_t ph = H.hkey[base + q + PF] & mask;
        __builtin_prefetch(&H.tkey[ph]); __builtin_prefetch(&H.tent[ph]);
      }
      const size_t f = base + q;
      const uint32_t p = H.hpid[f], pos = H.hpos[f];
      if (!place_key(H.hkey[f], p, (uint32_t)f)) { T.sel_skip += 1; if (pos >= (uint32_t)H.mstart[d]) T.sel_skip_mig += 1; continue; }
      record(d, k, d, p, pos); k++;
    }
    H.cur[d] = q; H.fillA[d] = k;
  }
  for (int d = 0; d < NDEME; d++) {
    int slot = H.fillA[d];
    for (int j = 0; j < NDEME && slot < DECK; j++) {
      const size_t base = (size_t)j * HMAX;
      const int mj = (int)S.sel[j].pidx.size();
      while (H.cur[j] < mj && slot < DECK) {
        const size_t f = base + (size_t)H.cur[j]++;
        const uint32_t p = H.hpid[f], pos = H.hpos[f];
        if (!place_key(H.hkey[f], p, (uint32_t)f)) { T.sel_skip += 1; if (pos >= (uint32_t)H.mstart[j]) T.sel_skip_mig += 1; continue; }
        record(d, slot, j, p, pos); slot++; T.xfill += 1;
      }
    }
    if (slot < DECK) { fprintf(stderr, "select: deme %d short by %d after both passes (impossible: parents are distinct)\n", d, DECK - slot); exit(1); }
  }

  H33_LAP("serial placement")
  const size_t NPL = H.plist.size(), CHK = 1024;
  pool.run((int)((NPL + CHK - 1) / CHK), [&](int c) {
    const size_t r1 = std::min(NPL, (size_t)(c + 1) * CHK);
    for (size_t r = (size_t)c * CHK; r < r1; r++) {
      if (r + 4 < NPL) {
        const uint32_t pn = H.plist[r + 4].p;
        if (pn < (uint32_t)NP) { __builtin_prefetch(&pm[pn]); __builtin_prefetch(&S.par[pn]); }
        else { __builtin_prefetch(&cm[pn - NP]); __builtin_prefetch(&S.chi[pn - NP]); }
      }
      const SelH33::PRec& q = H.plist[r];
      const uint32_t p = q.p;
      const Meta& mt = p < (uint32_t)NP ? pm[p] : cm[p - NP];
      const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
      const SelScratch& sj = S.sel[q.j];
      Ind& o = S.nxt[q.dst_slot];
      o = src; o.rank = sj.prank[q.pos]; o.crowd = sj.pcrowd[q.pos]; o.err = sj.pe[q.pos]; o.valid = 1;
      nm[q.dst_slot] = mt;
    }
  });
  H33_LAP("copies")
}
