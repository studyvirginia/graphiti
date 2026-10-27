#pragma once

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>

enum { H128P_MIG, H128P_POOL, H128P_KEYS, H128P_PASSA, H128P_PASSB, H128P_FREE, H128P_COPY, H128P_NS };
struct H128Prof { double ms[H128P_NS] = {0}; long n = 0; };
static H128Prof g_h128_prof;
static bool g_h128_prof_on = false;
static int g_h128 = 1;
inline std::chrono::steady_clock::time_point h128_now() { return std::chrono::steady_clock::now(); }
inline void h128_add(int k, std::chrono::steady_clock::time_point& t0) {
  auto t1 = h128_now();
  if (g_h128_prof_on) g_h128_prof.ms[k] += std::chrono::duration<double, std::milli>(t1 - t0).count();
  t0 = t1;
}
static double g_h128_sub[64][4];
inline void h128_sub(int d, int k, std::chrono::steady_clock::time_point& t0) {
  auto t1 = h128_now();
  if (g_h128_prof_on) g_h128_sub[d][k] += std::chrono::duration<double, std::milli>(t1 - t0).count();
  t0 = t1;
}
inline void h128_print_sub(double nc) {
  double s[4] = {0, 0, 0, 0};
  for (int d = 0; d < 64; d++) for (int k = 0; k < 4; k++) s[k] += g_h128_sub[d][k];
  printf("H128 POOL SUB (thread-summed per call, ms): fill %.3f | rank_crowd (NDS+crowd) %.3f | survivors %.3f | hot %.3f\n",
         s[0] / nc, s[1] / nc, s[2] / nc, s[3] / nc);
}

struct SelH92 {
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
  std::vector<uint32_t> freeF;
  std::vector<uint8_t> used;
  struct ShT { std::vector<uint64_t> h; std::vector<uint32_t> rep, id; };
  std::vector<uint32_t> kid; std::vector<uint8_t> ksh; std::vector<size_t> kbase; std::vector<uint32_t> shN;
  std::vector<ShT> sh; std::vector<uint8_t> keyseen;
};
static SelH92 g_sel_h92;
static int g_h117_bit4 = 0;

inline evo::Ind& XA(evo::State& S, uint32_t s) {
  const uint32_t NP = (uint32_t)evo::NPOP;
  return s < NP ? S.par[s] : (s < 2 * NP ? S.chi[s - NP] : S.nxt[s - 2 * NP]);
}
inline const evo::Ind& XC(const evo::State& S, uint32_t s) {
  const uint32_t NP = (uint32_t)evo::NPOP;
  return s < NP ? S.par[s] : (s < 2 * NP ? S.chi[s - NP] : S.nxt[s - 2 * NP]);
}

#include <atomic>
static const int g_h740 = getenv("H740") ? atoi(getenv("H740")) : 0;
struct H740E { float err; uint16_t cx; uint8_t valid; uint8_t pad; };
struct H740R { double crowd; uint16_t rank; uint16_t pad[3]; };
static std::vector<H740E> g_h740e;
static std::vector<H740R> g_h740r;
static std::atomic<long> g_h740_bad{0};
static std::atomic<long> g_spd[8];
struct H740Rep { ~H740Rep() { if (g_h740 == 2) printf("H740 verify: %ld flat-vs-Ind mismatches\n", g_h740_bad.load()); } };
static H740Rep g_h740_rep;
inline void h740_rc(const evo::State& S, uint32_t slot, uint16_t& rk, double& cw) {
  if (g_h740 == 1) { const H740R& R = g_h740r[slot]; rk = R.rank; cw = R.crowd; return; }
  const evo::Ind& P = XC(S, slot); rk = P.rank; cw = P.crowd;
  if (g_h740 == 2) { const H740R& R = g_h740r[slot]; if (R.rank != P.rank || std::memcmp(&R.crowd, &P.crowd, sizeof(double))) g_h740_bad.fetch_add(1); }
}
inline void h740_ec(const evo::State& S, uint32_t slot, float& e, uint16_t& c) {
  if (g_h740 == 1) { const H740E& F = g_h740e[slot]; e = F.err; c = (uint16_t)F.cx; return; }
  const evo::Ind& P = XC(S, slot); e = P.err; c = (uint16_t)P.cx;
  if (g_h740 == 2) { const H740E& F = g_h740e[slot]; if (std::memcmp(&F.err, &P.err, sizeof(float)) || F.cx != P.cx) g_h740_bad.fetch_add(1); }
}
inline uint8_t h740_valid(const evo::State& S, uint32_t slot) {
  if (g_h740 == 1) return g_h740e[slot].valid;
  const evo::Ind& P = XC(S, slot);
  if (g_h740 == 2 && g_h740e[slot].valid != P.valid) g_h740_bad.fetch_add(1);
  return P.valid;
}
inline float h740_perr(const evo::State& S, uint32_t slot) {
  if (g_h740 == 1) return g_h740e[slot].err;
  const evo::Ind& P = XC(S, slot);
  if (g_h740 == 2 && std::memcmp(&g_h740e[slot].err, &P.err, sizeof(float))) g_h740_bad.fetch_add(1);
  return P.err;
}

inline void next_children_h92(const uint32_t* pidx, uint32_t* cs) {
  const size_t NS3 = (size_t)3 * evo::NPOP;
  SelH92& H = g_sel_h92;
  H.used.assign(NS3, 0);
  for (int k = 0; k < evo::NPOP; k++) H.used[pidx[k]] = 1;
  size_t c = 0;
  for (size_t s = 0; s < NS3 && c < (size_t)evo::NPOP; s++) if (!H.used[s]) cs[c++] = (uint32_t)s;
  if (c != (size_t)evo::NPOP) { fprintf(stderr, "next_children_h92: only %zu free slots\n", c); exit(1); }
}

inline void breed_h92(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, evo::BreedBuf& BB,
                      uint64_t seed, uint64_t step, uint32_t* rebreeds, uint32_t* failed) {
  using namespace evo;
  pool.run(NDEME, [&](int d) {
    for (int i = d * DECK; i < (d + 1) * DECK; i++) { const Ind& P = XC(S, pidx[i]); BB.G[i] = P.g; if (g_h740) h740_rc(S, pidx[i], BB.prk[i], BB.pcw[i]); else { BB.prk[i] = P.rank; BB.pcw[i] = P.crowd; } }
    Rng r(mix(mix(seed, step), (uint64_t)d));
    uint32_t rb = 0, fl = 0;
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      const bool ok_ = make_child_i12(BB, d * DECK, r, XA(S, cs[i]), rb, g_link320 ? g_link320 + i : nullptr);
      if (!ok_) fl++;
      if (g_h740) g_h740e[cs[i]].valid = ok_ ? 1 : 0;
    }
    rebreeds[d] = rb; failed[d] = fl;
  });
}

inline void migrants_h92(int d, const evo::State& S, const uint32_t* pidx, uint32_t* out, std::vector<SelH92::MK>& mk) {
  using namespace evo;
  mk.resize((size_t)DECK);
  for (int i = 0; i < DECK; i++) {
    const uint32_t a = (uint32_t)(d * DECK + i);
    if (g_h740) { uint16_t rk_; double cw_; h740_rc(S, pidx[a], rk_, cw_); mk[(size_t)i] = SelH92::MK{rk_, a, cw_}; continue; }
    const Ind& P = XC(S, pidx[a]);
    mk[(size_t)i] = SelH92::MK{P.rank, a, P.crowd};
  }
  std::partial_sort(mk.begin(), mk.begin() + MIGR, mk.end(), [](const SelH92::MK& x, const SelH92::MK& y) {
    if (x.rank != y.rank) return x.rank < y.rank;
    if (x.crowd != y.crowd) return x.crowd > y.crowd;
    return x.idx < y.idx;
  });
  for (int k = 0; k < MIGR; k++) out[k] = mk[(size_t)k].idx;
}

inline void h117_keyids(evo::Pool& pool, SelH92& H, const std::vector<Meta>& pm, const std::vector<Meta>& cm, uint32_t NP,
                        const std::vector<int>& mlen, size_t HMAX) {
  using namespace evo;
  const size_t NE = (size_t)NDEME * HMAX;
  const int NSh = 8;
  H.kid.resize(NE); H.ksh.resize(NE);
  std::vector<size_t> cnt((size_t)NSh, 0);
  for (int d = 0; d < NDEME; d++) for (int q = 0; q < mlen[(size_t)d]; q++) {
    const size_t f = (size_t)d * HMAX + (size_t)q; cnt[(size_t)((H.hkey[f] >> 32) % (uint64_t)NSh)]++; }
  if (H.sh.size() != (size_t)NSh) H.sh.resize((size_t)NSh);
  H.shN.assign((size_t)NSh, 0);
  auto metaOf = [&](uint32_t p) -> const Meta& { return p < NP ? pm[p] : cm[p - NP]; };
  pool.run(NSh, [&](int s) {
    SelH92::ShT& A = H.sh[(size_t)s];
    size_t cap = 16; while (cap < 2 * cnt[(size_t)s] + 2) cap <<= 1;
    A.h.assign(cap, 0); A.rep.assign(cap, 0xFFFFFFFFu); A.id.assign(cap, 0);
    const uint64_t mask = (uint64_t)cap - 1;
    uint32_t nid = 0;
    for (int d = 0; d < NDEME; d++) for (int q = 0; q < mlen[(size_t)d]; q++) {
      const size_t f = (size_t)d * HMAX + (size_t)q;
      const uint64_t h = H.hkey[f];
      if ((int)((h >> 32) % (uint64_t)NSh) != s) continue;
      const uint32_t p = H.hpid[f];
      uint64_t i = h & mask; uint32_t id = 0xFFFFFFFFu;
      while (A.rep[i] != 0xFFFFFFFFu) {
        if (A.h[i] == h && metaOf(H.hpid[A.rep[i]]).key == metaOf(p).key) { id = A.id[i]; break; }
        i = (i + 1) & mask;
      }
      if (id == 0xFFFFFFFFu) { id = nid++; A.h[i] = h; A.rep[i] = (uint32_t)f; A.id[i] = id; }
      H.kid[f] = id; H.ksh[f] = (uint8_t)s;
    }
    H.shN[(size_t)s] = nid;
  });
  H.kbase.assign((size_t)NSh, 0);
  size_t acc = 0;
  for (int s = 0; s < NSh; s++) { H.kbase[(size_t)s] = acc; acc += H.shN[(size_t)s]; }
  H.keyseen.assign(acc, 0);
}

static std::vector<uint32_t> g_h128_slist;
static std::vector<size_t> g_h128_sbase;
static std::vector<size_t> g_h128_doff;
inline void h128_keyids2(evo::Pool& pool, SelH92& H, const std::vector<Meta>& pm, const std::vector<Meta>& cm, uint32_t NP,
                         const std::vector<int>& mlen, size_t HMAX) {
  using namespace evo;
  const int NSh = 8;
  const size_t NE = (size_t)NDEME * HMAX;
  H.kid.resize(NE); H.ksh.resize(NE);
  std::vector<size_t> dcnt((size_t)NDEME * NSh, 0);
  pool.run(NDEME, [&](int d) {
    size_t* c = &dcnt[(size_t)d * NSh];
    for (int q = 0; q < mlen[(size_t)d]; q++) {
      const size_t f = (size_t)d * HMAX + (size_t)q;
      const uint8_t s = (uint8_t)((H.hkey[f] >> 32) % (uint64_t)NSh);
      H.ksh[f] = s; c[s]++;
    }
  });
  g_h128_doff.assign((size_t)NDEME * NSh, 0);
  std::vector<size_t> tot((size_t)NSh, 0);
  for (int s = 0; s < NSh; s++) {
    size_t acc = 0;
    for (int d = 0; d < NDEME; d++) { g_h128_doff[(size_t)d * NSh + s] = acc; acc += dcnt[(size_t)d * NSh + s]; }
    tot[(size_t)s] = acc;
  }
  g_h128_sbase.assign((size_t)NSh + 1, 0);
  for (int s = 0; s < NSh; s++) g_h128_sbase[(size_t)s + 1] = g_h128_sbase[(size_t)s] + tot[(size_t)s];
  g_h128_slist.resize(g_h128_sbase[(size_t)NSh]);
  pool.run(NDEME, [&](int d) {
    size_t pos[8];
    for (int s = 0; s < NSh; s++) pos[s] = g_h128_sbase[(size_t)s] + g_h128_doff[(size_t)d * NSh + s];
    for (int q = 0; q < mlen[(size_t)d]; q++) {
      const size_t f = (size_t)d * HMAX + (size_t)q;
      g_h128_slist[pos[H.ksh[f]]++] = (uint32_t)f;
    }
  });
  if (H.sh.size() != (size_t)NSh) H.sh.resize((size_t)NSh);
  H.shN.assign((size_t)NSh, 0);
  auto metaOf = [&](uint32_t p) -> const Meta& { return p < NP ? pm[p] : cm[p - NP]; };
  pool.run(NSh, [&](int s) {
    SelH92::ShT& A = H.sh[(size_t)s];
    size_t cap = 16; while (cap < 2 * tot[(size_t)s] + 2) cap <<= 1;
    A.h.assign(cap, 0); A.rep.assign(cap, 0xFFFFFFFFu); A.id.assign(cap, 0);
    const uint64_t mask = (uint64_t)cap - 1;
    uint32_t nid = 0;
    for (size_t j = g_h128_sbase[(size_t)s]; j < g_h128_sbase[(size_t)s + 1]; j++) {
      const size_t f = g_h128_slist[j];
      const uint64_t h = H.hkey[f];
      const uint32_t p = H.hpid[f];
      uint64_t i = h & mask; uint32_t id = 0xFFFFFFFFu;
      while (A.rep[i] != 0xFFFFFFFFu) {
        if (A.h[i] == h && metaOf(H.hpid[A.rep[i]]).key == metaOf(p).key) { id = A.id[i]; break; }
        i = (i + 1) & mask;
      }
      if (id == 0xFFFFFFFFu) { id = nid++; A.h[i] = h; A.rep[i] = (uint32_t)f; A.id[i] = id; }
      H.kid[f] = id;
    }
    H.shN[(size_t)s] = nid;
  });
  H.kbase.assign((size_t)NSh, 0);
  size_t acc = 0;
  for (int s = 0; s < NSh; s++) { H.kbase[(size_t)s] = acc; acc += H.shN[(size_t)s]; }
  H.keyseen.assign(acc, 0);
}

inline void select_h92(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, uint32_t* pidx_out,
                       std::vector<Meta>& pm, std::vector<Meta>& cm, const uint8_t* admit, std::vector<Meta>& nm, Times& T,
                       bool move_mode) {
  using namespace evo;
  SelH92& H = g_sel_h92;
  const uint32_t NP = (uint32_t)NPOP;
  const size_t HMAX = (size_t)2 * DECK + MIGR;
  if (H.hkey.size() != (size_t)NDEME * HMAX) {
    H.hkey.assign((size_t)NDEME * HMAX, 0); H.hpid.assign((size_t)NDEME * HMAX, 0); H.hpos.assign((size_t)NDEME * HMAX, 0);
  }
  if (H.mk.size() != (size_t)NDEME) H.mk.resize(NDEME);
  H.mstart.assign(NDEME, 0); H.cur.assign(NDEME, 0); H.fillA.assign(NDEME, 0);

  auto tq = h128_now();
  pool.run(NDEME, [&](int d) { migrants_h92(d, S, pidx, S.migtab.data() + d * MIGR, H.mk[(size_t)d]); });
  h128_add(H128P_MIG, tq);
  pool.run(NDEME, [&](int d) {
    auto tA = h128_now();
    SelScratch& s = S.sel[d];
    s.pidx.clear();
    for (int i = d * DECK; i < (d + 1) * DECK; i++) s.pidx.push_back((uint32_t)i);
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      if (!admit[i]) continue;
      if (g_h740 == 1) { s.pidx.push_back((uint32_t)(NP + i)); continue; }
      Ind& c = XA(S, cs[i]);
      float sum; if (g_h419) sum = 0.0f + c.res[0]; else { sum = 0.0f; for (int q = 0; q < NRES; q++) sum += c.res[q]; }
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
      const uint32_t sl_ = p < NP ? pidx[p] : cs[p - NP];
      if (g_h740) { float e_; uint16_t c_; h740_ec(S, sl_, e_, c_); s.pe[q] = e_; s.pc[q] = c_; }
      else { const Ind& src = XC(S, sl_); s.pe[q] = src.err; s.pc[q] = (uint16_t)src.cx; }
    }
    h128_sub(d, 0, tA);
    rank_crowd(m, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    h128_sub(d, 1, tA);
    choose_survivors_h53(m, s.pcrowd.data(), s.rs, s.surv.data(), s.tmp);
    if (g_sp_diag) {
      int bq = -1; double bt = INFINITY;
      for (int q = 0; q < m; q++) { if (!(s.pe[q] < 1e19f)) continue; const uint32_t p = s.pidx[q]; const Meta& mt = p < NP ? pm[p] : cm[p - NP];
        const double t = (double)s.pe[q] + fp::search_pen(mt.kf, 48, mt.nx, mt.ny); if (t < bt) { bt = t; bq = q; } }
      if (bq >= 0) { int pos = 0; while (pos < m && (int)s.surv[pos] != bq) pos++;
        g_spd[0]++; if (s.prank[bq] == 1) g_spd[1]++; else if (s.prank[bq] == 2) g_spd[2]++; else g_spd[3]++; if (pos >= DECK) g_spd[4]++; if (bq >= (int)(H.mstart[d])) g_spd[5]++; } }
    h128_sub(d, 2, tA);
    const size_t base = (size_t)d * HMAX;
    for (int q = 0; q < m; q++) {
      const uint32_t pos = s.surv[q], p = s.pidx[pos];
      const Meta& mt = p < NP ? pm[p] : cm[p - NP];
      H.hkey[base + q] = mt.kh; H.hpid[base + q] = p; H.hpos[base + q] = pos;
    }
    h128_sub(d, 3, tA);
  });

  h128_add(H128P_POOL, tq);
  const bool dk = g_h117_bit4 != 0;
  if (dk) {
    std::vector<int> mlen((size_t)NDEME, 0);
    for (int d = 0; d < NDEME; d++) mlen[(size_t)d] = (int)S.sel[d].pidx.size();
    if (g_h128 & 1) h128_keyids2(pool, H, pm, cm, NP, mlen, HMAX);
    else h117_keyids(pool, H, pm, cm, NP, mlen, HMAX);
  }
  h128_add(H128P_KEYS, tq);
  size_t cap = 1024; while (cap < 2 * (size_t)NP + 2) cap <<= 1;
  if (H.tab.size() != cap) { H.tab.assign(cap, SelH92::Slot{0, 0, 0}); H.gen = 0; }
  if (++H.gen == 0) { std::fill(H.tab.begin(), H.tab.end(), SelH92::Slot{0, 0, 0}); H.gen = 1; }
  const uint32_t gen = H.gen;
  const uint64_t mask = (uint64_t)cap - 1;
  auto metaOf = [&](uint32_t p) -> const Meta& { return p < NP ? pm[p] : cm[p - NP]; };
  auto place_key = [&](uint64_t h, uint32_t p, uint32_t f) -> bool {
    uint64_t i = h & mask;
    while (H.tab[i].gen == gen) {
      if (H.tab[i].h == h && metaOf(H.hpid[H.tab[i].t - 1]).key == metaOf(p).key) return false;
      i = (i + 1) & mask;
    }
    H.tab[i] = SelH92::Slot{h, f + 1, gen}; return true;
  };
  auto placeAt = [&](size_t f, uint32_t p) -> bool {
    if (dk) { const size_t g = H.kbase[H.ksh[f]] + H.kid[f]; if (H.keyseen[g]) return false; H.keyseen[g] = 1; return true; }
    return place_key(H.hkey[f], p, (uint32_t)f);
  };
  H.plist.clear(); H.plist.reserve((size_t)NP);
  auto record = [&](int dst, int slot, int j, uint32_t p, uint32_t pos) {
    H.plist.push_back(SelH92::PRec{(uint32_t)(dst * DECK + slot), p, pos, (uint32_t)j});
    if (p >= NP) T.surv_child += 1;
  };
  const size_t PF = 8;
  for (int d = 0; d < NDEME; d++) {
    const size_t base = (size_t)d * HMAX;
    const int m = (int)S.sel[d].pidx.size();
    int k = 0, q = 0;
    for (; q < m && k < DECK; q++) {
      if (!dk && q + (int)PF < m) __builtin_prefetch(&H.tab[H.hkey[base + q + PF] & mask]);
      const size_t f = base + q;
      const uint32_t p = H.hpid[f], pos = H.hpos[f];
      if (!placeAt(f, p)) { T.sel_skip += 1; if (pos >= (uint32_t)H.mstart[d]) T.sel_skip_mig += 1; continue; }
      record(d, k, d, p, pos); k++;
    }
    H.cur[d] = q; H.fillA[d] = k;
  }
  h128_add(H128P_PASSA, tq);
  for (int d = 0; d < NDEME; d++) {
    int slot = H.fillA[d];
    for (int j = 0; j < NDEME && slot < DECK; j++) {
      const size_t base = (size_t)j * HMAX;
      const int mj = (int)S.sel[j].pidx.size();
      while (H.cur[j] < mj && slot < DECK) {
        const size_t f = base + (size_t)H.cur[j]++;
        const uint32_t p = H.hpid[f], pos = H.hpos[f];
        if (!placeAt(f, p)) { T.sel_skip += 1; if (pos >= (uint32_t)H.mstart[j]) T.sel_skip_mig += 1; continue; }
        record(d, slot, j, p, pos); slot++; T.xfill += 1;
      }
    }
    if (slot < DECK) { fprintf(stderr, "select_h92: deme %d short by %d after both passes\n", d, DECK - slot); exit(1); }
  }
  h128_add(H128P_PASSB, tq);
  const size_t NPL = H.plist.size();
  if (NPL != (size_t)NP) { fprintf(stderr, "select_h92: %zu placements, expected %u\n", NPL, NP); exit(1); }

  if (!move_mode) {
    const size_t NS3 = (size_t)3 * NP;
    H.used.assign(NS3, 0);
    for (uint32_t k = 0; k < NP; k++) { H.used[pidx[k]] = 1; H.used[cs[k]] = 1; }
    H.freeF.clear(); H.freeF.reserve(NP);
    for (size_t s = 0; s < NS3 && H.freeF.size() < (size_t)NP; s++) if (!H.used[s]) H.freeF.push_back((uint32_t)s);
    if (H.freeF.size() != (size_t)NP) { fprintf(stderr, "select_h92: %zu free slots, expected %u\n", H.freeF.size(), NP); exit(1); }
  }
  h128_add(H128P_FREE, tq);

  const size_t CHK = 1024;
  pool.run((int)((NPL + CHK - 1) / CHK), [&](int c) {
    const size_t r1 = std::min(NPL, (size_t)(c + 1) * CHK);
    for (size_t r = (size_t)c * CHK; r < r1; r++) {
      const SelH92::PRec& q = H.plist[r];
      const uint32_t p = q.p;
      const uint32_t src = p < NP ? pidx[p] : cs[p - NP];
      const SelScratch& sj = S.sel[q.j];
      uint32_t dst;
      if (move_mode) {
        dst = src;
      } else {
        dst = H.freeF[q.dst_slot];
        XA(S, dst) = XC(S, src);
      }
      Ind& o = XA(S, dst);
      o.rank = sj.prank[q.pos]; o.crowd = sj.pcrowd[q.pos]; o.err = sj.pe[q.pos]; o.valid = 1;
      if (g_h740) { H740E& fe = g_h740e[dst]; if (dst != src) fe.cx = g_h740e[src].cx; fe.err = sj.pe[q.pos]; fe.valid = 1;
                    H740R& fr = g_h740r[dst]; fr.rank = sj.prank[q.pos]; fr.crowd = sj.pcrowd[q.pos]; }
      pidx_out[q.dst_slot] = dst;
      Meta& mt = p < NP ? pm[p] : cm[p - NP];
      if (move_mode) nm[q.dst_slot] = std::move(mt);
      else nm[q.dst_slot] = mt;
    }
  });
  h128_add(H128P_COPY, tq);
  if (g_h128_prof_on) g_h128_prof.n++;
}
