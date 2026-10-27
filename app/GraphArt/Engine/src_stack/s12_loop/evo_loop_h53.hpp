
#pragma once

#include <chrono>

struct SelH53 {
  struct MK { uint16_t rank; uint32_t idx; double crowd; };
  struct Slot { uint64_t h; uint32_t t; uint32_t gen; };
  struct PRec { uint32_t dst_slot, p, pos, j; };
  std::vector<uint64_t> hkey;
  std::vector<uint32_t> hpid, hpos;
  std::vector<Slot> tab;
  uint32_t gen = 0;
  std::vector<int> mstart, cur, fillA;
  std::vector<PRec> plist;
  std::vector<std::vector<MK>> mk;
};
static SelH53 g_sel_h53;

#ifdef H53_PROF
#define H53_T0 auto _h53t = std::chrono::steady_clock::now();
#define H53_LAP(name) { auto _h53n = std::chrono::steady_clock::now(); \
  fprintf(stderr, "h53 %s %.3f ms\n", name, std::chrono::duration<double, std::milli>(_h53n - _h53t).count()); _h53t = _h53n; }
#else
#define H53_T0
#define H53_LAP(name)
#endif

inline void migrants_h53(int d, const evo::Ind* par, uint32_t* out, std::vector<SelH53::MK>& mk) {
  mk.resize((size_t)evo::DECK);
  for (int i = 0; i < evo::DECK; i++) {
    const uint32_t a = (uint32_t)(d * evo::DECK + i);
    mk[(size_t)i] = SelH53::MK{par[a].rank, a, par[a].crowd};
  }
  std::partial_sort(mk.begin(), mk.begin() + evo::MIGR, mk.end(), [](const SelH53::MK& x, const SelH53::MK& y) {
    if (x.rank != y.rank) return x.rank < y.rank;
    if (x.crowd != y.crowd) return x.crowd > y.crowd;
    return x.idx < y.idx;
  });
  for (int k = 0; k < evo::MIGR; k++) out[k] = mk[(size_t)k].idx;
}

inline void choose_survivors_h53(int n, const double* crowd, evo::RankScratch& sc, uint32_t* surv, std::vector<uint32_t>& tmp) {
  int k = 0;
  const int R = (int)sc.off.size() - 2;
  for (int r = 1; r <= R && k < n; r++) {
    const int a = (int)sc.off[r], b = (int)sc.off[r + 1], f = b - a;
    if (k + f <= n) {
      for (int q = a; q < b; q++) surv[k++] = evo::jof(sc.FE[q]);
    } else {
      const int need = n - k;
      tmp.clear();
      for (int q = a; q < b; q++) tmp.push_back(evo::jof(sc.FE[q]));
      std::partial_sort(tmp.begin(), tmp.begin() + need, tmp.end(), [&](uint32_t x, uint32_t y) {
        if (crowd[x] != crowd[y]) return crowd[x] > crowd[y];
        return x < y;
      });
      for (int q = 0; q < need; q++) surv[k++] = tmp[(size_t)q];
    }
  }
}

static void select_h53(evo::Pool& pool, evo::State& S, const std::vector<Meta>& pm, const std::vector<Meta>& cm,
                       const uint8_t* admit, std::vector<Meta>& nm, Times& T) {
  using namespace evo;
  H53_T0
  SelH53& H = g_sel_h53;
  const int NP = NPOP;
  const size_t HMAX = (size_t)2 * DECK + MIGR;
  if (H.hkey.size() != (size_t)NDEME * HMAX) {
    H.hkey.assign((size_t)NDEME * HMAX, 0); H.hpid.assign((size_t)NDEME * HMAX, 0); H.hpos.assign((size_t)NDEME * HMAX, 0);
  }
  if (H.mk.size() != (size_t)NDEME) H.mk.resize(NDEME);
  H.mstart.assign(NDEME, 0); H.cur.assign(NDEME, 0); H.fillA.assign(NDEME, 0);

  pool.run(NDEME, [&](int d) { migrants_h53(d, S.par.data(), S.migtab.data() + d * MIGR, H.mk[(size_t)d]); });
  H53_LAP("migrants")
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
    const int pd = (d + NDEME - 1) % NDEME;
    for (int k = 0; k < MIGR; k++) s.pidx.push_back(S.migtab[pd * MIGR + k]);
    const int m = (int)s.pidx.size();
    s.pe.resize(m); s.pc.resize(m); s.prank.resize(m); s.pcrowd.resize(m); s.surv.resize(m);
    for (int q = 0; q < m; q++) {
      const uint32_t p = s.pidx[q];
      const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
      s.pe[q] = src.err; s.pc[q] = (uint16_t)src.cx;
    }
    rank_crowd(m, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    choose_survivors_h53(m, s.pcrowd.data(), s.rs, s.surv.data(), s.tmp);
    const size_t base = (size_t)d * HMAX;
    for (int q = 0; q < m; q++) {
      const uint32_t pos = s.surv[q], p = s.pidx[pos];
      const Meta& mt = p < (uint32_t)NP ? pm[p] : cm[p - NP];
      H.hkey[base + q] = mt.kh; H.hpid[base + q] = p; H.hpos[base + q] = pos;
    }
  });
  H53_LAP("prep+rank")

  size_t cap = 1024; while (cap < 2 * (size_t)NP + 2) cap <<= 1;
  if (H.tab.size() != cap) { H.tab.assign(cap, SelH53::Slot{0, 0, 0}); H.gen = 0; }
  if (++H.gen == 0) { std::fill(H.tab.begin(), H.tab.end(), SelH53::Slot{0, 0, 0}); H.gen = 1; }
  const uint32_t gen = H.gen;
  const uint64_t mask = (uint64_t)cap - 1;
  auto metaOf = [&](uint32_t p) -> const Meta& { return p < (uint32_t)NP ? pm[p] : cm[p - NP]; };
  auto place_key = [&](uint64_t h, uint32_t p, uint32_t f) -> bool {
    uint64_t i = h & mask;
    while (H.tab[i].gen == gen) {
      if (H.tab[i].h == h && metaOf(H.hpid[H.tab[i].t - 1]).key == metaOf(p).key) return false;
      i = (i + 1) & mask;
    }
    H.tab[i] = SelH53::Slot{h, f + 1, gen}; return true;
  };
  H.plist.clear(); H.plist.reserve((size_t)NP);
  auto record = [&](int dst, int slot, int j, uint32_t p, uint32_t pos) {
    H.plist.push_back(SelH53::PRec{(uint32_t)(dst * DECK + slot), p, pos, (uint32_t)j});
    if (p >= (uint32_t)NP) T.surv_child += 1;
  };
  const size_t PF = 8;
  for (int d = 0; d < NDEME; d++) {
    const size_t base = (size_t)d * HMAX;
    const int m = (int)S.sel[d].pidx.size();
    int k = 0, q = 0;
    for (; q < m && k < DECK; q++) {
      if (q + (int)PF < m) __builtin_prefetch(&H.tab[H.hkey[base + q + PF] & mask]);
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
  H53_LAP("serial placement")
#ifdef H53_PROF
  { size_t pool_n = 0, visits = 0; for (int d = 0; d < NDEME; d++) { pool_n += S.sel[d].pidx.size(); visits += (size_t)H.cur[d]; }
    fprintf(stderr, "h53 select: pool entries %zu, visited %zu, placed %zu (children %.0f), skips %.0f (migrant %.0f), xfill %.0f\n", pool_n, visits, H.plist.size(), T.surv_child, T.sel_skip, T.sel_skip_mig, T.xfill); }
#endif

  const size_t NPL = H.plist.size(), CHK = 1024;
  pool.run((int)((NPL + CHK - 1) / CHK), [&](int c) {
    const size_t r1 = std::min(NPL, (size_t)(c + 1) * CHK);
    for (size_t r = (size_t)c * CHK; r < r1; r++) {
      if (r + 4 < NPL) {
        const uint32_t pn = H.plist[r + 4].p;
        if (pn < (uint32_t)NP) { __builtin_prefetch(&pm[pn]); __builtin_prefetch(&S.par[pn]); }
        else { __builtin_prefetch(&cm[pn - NP]); __builtin_prefetch(&S.chi[pn - NP]); }
      }
      const SelH53::PRec& q = H.plist[r];
      const uint32_t p = q.p;
      const Meta& mt = p < (uint32_t)NP ? pm[p] : cm[p - NP];
      const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
      const SelScratch& sj = S.sel[q.j];
      Ind& o = S.nxt[q.dst_slot];
      o = src; o.rank = sj.prank[q.pos]; o.crowd = sj.pcrowd[q.pos]; o.err = sj.pe[q.pos]; o.valid = 1;
      nm[q.dst_slot] = mt;
    }
  });
  H53_LAP("copies")
}
