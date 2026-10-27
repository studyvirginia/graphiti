#pragma once
#include <memory>

static std::vector<const evo::Cand*> g_h777_gp;

__attribute__((always_inline)) inline bool make_child_i12_p(const evo::BreedBuf& BB, const evo::Cand* const* Gp, int base, evo::Rng& r,
                                                            evo::Ind& out, uint32_t& rebreeds, evo::ChildLink320* pv = nullptr) {
  using namespace evo;
  for (int attempt = 0; attempt < MAX_TRY; attempt++) {
    const int ia = tournament_i12(BB, base, r);
    const int ib = tournament_i12(BB, base, r);
    const evo::Cand& A = *Gp[ia];
    const evo::Cand& B = *Gp[ib];
    const bool swap = r.uni() < P_SWAP;
    if (pv) { std::memset(pv, 0, sizeof *pv); pv->ia = ia; pv->ib = ib; pv->swap = swap ? 1 : 0; pv->k[0] = 255; pv->k[1] = 255; }
    evo::Cand g;
    bool ok = true;
    for (int t = 0; t < 2 && ok; t++) {
      if (r.uni() < P_XO) {
        const int tb = t ^ (int)swap;
        const int LA = A.len[t], LB = B.len[tb];
        const int i = (int)r.below(LA), j = (int)r.below(LB);
        const int s = subtree_start(A.op[t], i), szi = i - s + 1;
        const int js = subtree_start(B.op[tb], j), szj = j - js + 1;
        const int nl = LA - szi + szj;
        if (nl > MAXLEN) { ok = false; rebreeds++; break; }
        splice_tree(g, t, A, B, tb, s, i, js, szj, nl);
        if (pv) { pv->xo[t] = 1; pv->s[t] = (uint8_t)s; pv->szi[t] = (uint8_t)szi; pv->js[t] = (uint8_t)js; pv->szj[t] = (uint8_t)szj; }
      } else {
        std::memcpy(g.op[t], A.op[t], MAXLEN);
        std::memcpy(g.c[t], A.c[t], MAXLEN * sizeof(float));
        g.len[t] = A.len[t];
      }
    }
    if (!ok) continue;
    for (int t = 0; t < 2; t++) if (r.uni() < P_MUT) { int kk = 0; mutate(g, t, r, &kk); if (pv) pv->k[t] = (uint8_t)kk; }
    out.g = g;
    out.valid = 1;
    return true;
  }
  out.valid = 0;
  return false;
}

inline void breed_h777(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, evo::BreedBuf& BB,
                       uint64_t seed, uint64_t step, uint32_t* rebreeds, uint32_t* failed) {
  using namespace evo;
  if (g_h777_gp.size() != (size_t)NPOP) g_h777_gp.assign((size_t)NPOP, nullptr);
  const evo::Cand** Gp = g_h777_gp.data();
  pool.run(NDEME, [&](int d) {
    for (int i = d * DECK; i < (d + 1) * DECK; i++) { const Ind& P = XC(S, pidx[i]); Gp[i] = &P.g; if (g_h740) h740_rc(S, pidx[i], BB.prk[i], BB.pcw[i]); else { BB.prk[i] = P.rank; BB.pcw[i] = P.crowd; } }
    Rng r(mix(mix(seed, step), (uint64_t)d));
    uint32_t rb = 0, fl = 0;
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      const bool ok_ = make_child_i12_p(BB, Gp, d * DECK, r, XA(S, cs[i]), rb, g_link320 ? g_link320 + i : nullptr);
      if (!ok_) fl++;
      if (g_h740) g_h740e[cs[i]].valid = ok_ ? 1 : 0;
    }
    rebreeds[d] = rb; failed[d] = fl;
  });
}

struct Sh774 { std::vector<uint64_t> h; std::vector<uint32_t> rep, id; size_t cap = 0, cnt = 0; uint32_t nid = 0; };
struct Sel774 {
  SelH92 H;
  std::vector<Sh774> sh;
  std::vector<uint32_t> shN;
  int next_keys = 0;
  size_t HMAX = 0;
  double ms[4] = {0, 0, 0, 0};
  int passB_demes = 0;
  long passB_total = 0, calls = 0;
};
static Sel774 g_sel774;
static const uint32_t SEL774_EMPTY = 0xFFFFFFFFu;

inline double h774_ms(std::chrono::steady_clock::time_point& t0) {
  auto t1 = std::chrono::steady_clock::now(); const double v = std::chrono::duration<double, std::milli>(t1 - t0).count(); t0 = t1; return v;
}
inline void sh774_reset(Sh774& A) {
  if (A.cap < 1024) A.cap = 1024;
  A.h.assign(A.cap, 0); A.rep.assign(A.cap, SEL774_EMPTY); A.id.assign(A.cap, 0); A.cnt = 0; A.nid = 0;
}
inline void sh774_grow(Sh774& A) {
  const size_t nc = A.cap * 2; const uint64_t mask = (uint64_t)nc - 1;
  std::vector<uint64_t> h(nc, 0); std::vector<uint32_t> rep(nc, SEL774_EMPTY), id(nc, 0);
  for (size_t j = 0; j < A.cap; j++) if (A.rep[j] != SEL774_EMPTY) {
    uint64_t i = A.h[j] & mask; while (rep[i] != SEL774_EMPTY) i = (i + 1) & mask;
    h[i] = A.h[j]; rep[i] = A.rep[j]; id[i] = A.id[j];
  }
  A.h.swap(h); A.rep.swap(rep); A.id.swap(id); A.cap = nc;
}

inline void sel774_begin(evo::Pool& pool, evo::State& S, const uint32_t* pidx) {
  using namespace evo;
  Sel774& G = g_sel774; SelH92& H = G.H;
  auto t0 = std::chrono::steady_clock::now();
  G.HMAX = (size_t)2 * DECK + MIGR;
  if (H.hkey.size() != (size_t)NDEME * G.HMAX) {
    H.hkey.assign((size_t)NDEME * G.HMAX, 0); H.hpid.assign((size_t)NDEME * G.HMAX, 0); H.hpos.assign((size_t)NDEME * G.HMAX, 0);
  }
  if (H.mk.size() != (size_t)NDEME) H.mk.resize(NDEME);
  H.mstart.assign(NDEME, 0); H.cur.assign(NDEME, 0); H.fillA.assign(NDEME, 0);
  if (g_h117_bit4) {
    H.kid.resize((size_t)NDEME * G.HMAX); H.ksh.resize((size_t)NDEME * G.HMAX);
    if (G.sh.size() != 8) G.sh.resize(8);
    for (Sh774& A : G.sh) sh774_reset(A);
  }
  G.next_keys = 0; G.passB_demes = 0;
  for (int k = 0; k < 4; k++) G.ms[k] = 0;
  pool.run(NDEME, [&](int d) { migrants_h92(d, S, pidx, S.migtab.data() + d * MIGR, H.mk[(size_t)d]); });
  G.ms[0] += h774_ms(t0);
}

inline void sel774_demes(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, std::vector<Meta>& pm, std::vector<Meta>& cm,
                         const uint8_t* admit, int lo, int hi) {
  using namespace evo;
  Sel774& G = g_sel774; SelH92& H = G.H;
  const uint32_t NP = (uint32_t)NPOP;
  if (lo != G.next_keys) { fprintf(stderr, "sel774_demes: demes out of order (%d, expected %d)\n", lo, G.next_keys); exit(1); }
  auto t0 = std::chrono::steady_clock::now();
  pool.run(hi - lo, [&](int k) {
    const int d = lo + k;
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
    for (int q = 0; q < MIGR; q++) s.pidx.push_back(S.migtab[pd * MIGR + q]);
    const int m = (int)s.pidx.size();
    s.pe.resize(m); s.pc.resize(m); s.prank.resize(m); s.pcrowd.resize(m); s.surv.resize(m);
    for (int q = 0; q < m; q++) {
      const uint32_t p = s.pidx[q];
      const uint32_t sl_ = p < NP ? pidx[p] : cs[p - NP];
      if (g_h740) { float e_; uint16_t c_; h740_ec(S, sl_, e_, c_); s.pe[q] = e_; s.pc[q] = c_; }
      else { const Ind& src = XC(S, sl_); s.pe[q] = src.err; s.pc[q] = (uint16_t)src.cx; }
    }
    rank_crowd(m, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    choose_survivors_h53(m, s.pcrowd.data(), s.rs, s.surv.data(), s.tmp);
    const size_t base = (size_t)d * G.HMAX;
    for (int q = 0; q < m; q++) {
      const uint32_t pos = s.surv[q], p = s.pidx[pos];
      const Meta& mt = p < NP ? pm[p] : cm[p - NP];
      H.hkey[base + q] = mt.kh; H.hpid[base + q] = p; H.hpos[base + q] = pos;
    }
  });
  G.ms[1] += h774_ms(t0);
  if (g_h117_bit4) {
    const int NSh = 8;
    pool.run(hi - lo, [&](int k) {
      const int d = lo + k; const int m = (int)S.sel[d].pidx.size();
      for (int q = 0; q < m; q++) { const size_t f = (size_t)d * G.HMAX + (size_t)q; H.ksh[f] = (uint8_t)((H.hkey[f] >> 32) % (uint64_t)NSh); }
    });
    auto metaOf = [&](uint32_t p) -> const Meta& { return p < NP ? pm[p] : cm[p - NP]; };
    pool.run(NSh, [&](int s) {
      Sh774& A = G.sh[(size_t)s];
      for (int d = lo; d < hi; d++) {
        const int m = (int)S.sel[d].pidx.size();
        for (int q = 0; q < m; q++) {
          const size_t f = (size_t)d * G.HMAX + (size_t)q;
          if ((int)H.ksh[f] != s) continue;
          if ((A.cnt + 1) * 2 > A.cap) sh774_grow(A);
          const uint64_t mask = (uint64_t)A.cap - 1;
          const uint64_t h = H.hkey[f]; const uint32_t p = H.hpid[f];
          uint64_t i = h & mask; uint32_t id = SEL774_EMPTY;
          while (A.rep[i] != SEL774_EMPTY) {
            if (A.h[i] == h && metaOf(H.hpid[A.rep[i]]).key == metaOf(p).key) { id = A.id[i]; break; }
            i = (i + 1) & mask;
          }
          if (id == SEL774_EMPTY) { id = A.nid++; A.h[i] = h; A.rep[i] = (uint32_t)f; A.id[i] = id; A.cnt++; }
          H.kid[f] = id;
        }
      }
    });
    G.ms[2] += h774_ms(t0);
  }
  G.next_keys = hi;
}

inline void sel774_finish(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, uint32_t* pidx_out,
                          std::vector<Meta>& pm, std::vector<Meta>& cm, std::vector<Meta>& nm, Times& T, bool move_mode) {
  using namespace evo;
  Sel774& G = g_sel774; SelH92& H = G.H;
  if (G.next_keys != NDEME) { fprintf(stderr, "sel774_finish: only %d demes delivered\n", G.next_keys); exit(1); }
  auto t0 = std::chrono::steady_clock::now();
  const uint32_t NP = (uint32_t)NPOP;
  const size_t HMAX = G.HMAX;
  const bool dk = g_h117_bit4 != 0;
  if (dk) {
    G.shN.assign(8, 0); H.kbase.assign(8, 0);
    size_t acc = 0;
    for (int s = 0; s < 8; s++) { G.shN[(size_t)s] = G.sh[(size_t)s].nid; H.kbase[(size_t)s] = acc; acc += G.sh[(size_t)s].nid; }
    H.keyseen.assign(acc, 0);
  }
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
  for (int d = 0; d < NDEME; d++) {
    int slot = H.fillA[d];
    if (slot < DECK) G.passB_demes++;
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
    if (slot < DECK) { fprintf(stderr, "sel774_finish: deme %d short by %d after both passes\n", d, DECK - slot); exit(1); }
  }
  const size_t NPL = H.plist.size();
  if (NPL != (size_t)NP) { fprintf(stderr, "sel774_finish: %zu placements, expected %u\n", NPL, NP); exit(1); }
  if (!move_mode) {
    const size_t NS3 = (size_t)3 * NP;
    H.used.assign(NS3, 0);
    for (uint32_t k = 0; k < NP; k++) { H.used[pidx[k]] = 1; H.used[cs[k]] = 1; }
    H.freeF.clear(); H.freeF.reserve(NP);
    for (size_t s = 0; s < NS3 && H.freeF.size() < (size_t)NP; s++) if (!H.used[s]) H.freeF.push_back((uint32_t)s);
    if (H.freeF.size() != (size_t)NP) { fprintf(stderr, "sel774_finish: %zu free slots, expected %u\n", H.freeF.size(), NP); exit(1); }
  }
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
  G.ms[3] += h774_ms(t0);
  G.passB_total += G.passB_demes; G.calls++;
}

inline void h774_run(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, uint32_t* pidx_out,
                     std::vector<Meta>& pm, std::vector<Meta>& cm, const uint8_t* admit, std::vector<Meta>& nm, Times& T,
                     bool move_mode, int st, const std::vector<int>& cstart, int pre = -1) {
  using namespace evo;
  const size_t NP = (size_t)NPOP;
  std::unique_ptr<State> S2; std::vector<Meta> pm2, cm2, nm2; std::vector<uint32_t> pidx2c, cs2c, out2; std::vector<uint8_t> adm2;
  if (g_h774 == 2) {
    S2.reset(new State(S)); pm2 = pm; cm2 = cm; nm2.resize(nm.size()); pidx2c.assign(pidx, pidx + NP); cs2c.assign(cs, cs + NP);
    adm2.assign(admit, admit + NP); out2.assign(NP, 0);
  }
  Times Tp;
  if (pre < 0) { sel774_begin(pool, S, pidx); sel774_demes(pool, S, pidx, cs, pm, cm, admit, 0, NDEME); }
  else if (pre < NDEME) sel774_demes(pool, S, pidx, cs, pm, cm, admit, pre, NDEME);
  sel774_finish(pool, S, pidx, cs, pidx_out, pm, cm, nm, Tp, move_mode);
  Sel774& G = g_sel774;
  printf("H774 step %d: begin %.3f | demes %.3f | keys %.3f | finish %.3f ms | Pass B demes %d (total %ld over %ld calls)\n",
         st, G.ms[0], G.ms[1], G.ms[2], G.ms[3], G.passB_demes, G.passB_total, G.calls);
  if (g_h774 == 2 || st == 2) {
    printf("H774 ready map step %d (deme:last chunk of %zu; err/admit are final only after post16 fin, i.e. after chunk %zu):", st, cstart.size(), cstart.empty() ? (size_t)0 : cstart.size() - 1);
    for (int d = 0; d < NDEME; d++) {
      const int last = (d + 1) * DECK - 1; int c = 0;
      while (c + 1 < (int)cstart.size() && cstart[(size_t)c + 1] <= last) c++;
      printf(" %d:%d", d, c);
    }
    printf("\n");
  }
  if (g_h774 == 2) {
    bool prof = g_h128_prof_on; g_h128_prof_on = false;
    Times Tr;
    select_h92(pool, *S2, pidx2c.data(), cs2c.data(), out2.data(), pm2, cm2, adm2.data(), nm2, Tr, move_mode);
    g_h128_prof_on = prof;
    long mp = 0, mrk = 0, mcr = 0, mer = 0, mvl = 0, mmt = 0, mmg = 0, mpool = 0, mhot = 0, mkid = 0, mcnt = 0;
    for (size_t k = 0; k < NP; k++) {
      if (pidx_out[k] != out2[k]) mp++;
      const Ind& a = XC(S, pidx_out[k]); const Ind& b = XC(*S2, out2[k]);
      if (a.rank != b.rank) mrk++;
      if (std::memcmp(&a.crowd, &b.crowd, sizeof a.crowd)) mcr++;
      if (std::memcmp(&a.err, &b.err, sizeof a.err)) mer++;
      if (a.valid != b.valid || a.cx != b.cx) mvl++;
      if (nm[k].kh != nm2[k].kh || !(nm[k].key == nm2[k].key) || nm[k].r0 != nm2[k].r0 || nm[k].r1 != nm2[k].r1) mmt++;
    }
    for (size_t k = 0; k < (size_t)NDEME * MIGR; k++) if (S.migtab[k] != S2->migtab[k]) mmg++;
    SelH92& A = G.H; SelH92& R = g_sel_h92;
    for (int d = 0; d < NDEME; d++) {
      const SelScratch& a = S.sel[d]; const SelScratch& b = S2->sel[d];
      if (a.pidx != b.pidx) mpool++;
      const int m = (int)std::min(a.pidx.size(), b.pidx.size());
      if (A.mstart[d] != R.mstart[d]) mpool++;
      for (int q = 0; q < m; q++) {
        const size_t f = (size_t)d * G.HMAX + (size_t)q;
        if (A.hkey[f] != R.hkey[f] || A.hpid[f] != R.hpid[f] || A.hpos[f] != R.hpos[f]) mhot++;
        if (g_h117_bit4 && (A.kid[f] != R.kid[f] || A.ksh[f] != R.ksh[f])) mkid++;
      }
    }
    if (g_h117_bit4) { for (int s = 0; s < 8; s++) { if (G.shN[(size_t)s] != R.shN[(size_t)s]) mkid++; if (A.kbase[(size_t)s] != R.kbase[(size_t)s]) mkid++; } }
    if (Tp.sel_skip != Tr.sel_skip) mcnt++;
    if (Tp.sel_skip_mig != Tr.sel_skip_mig) mcnt++;
    if (Tp.surv_child != Tr.surv_child) mcnt++;
    if (Tp.xfill != Tr.xfill) mcnt++;
    printf("H774 verify step %d: mismatches pidx %ld rank %ld crowd %ld err %ld cx/valid %ld meta %ld migtab %ld pool %ld hot %ld keyids %ld counters %ld | Pass B demes %d (ref xfill %.0f)\n",
           st, mp, mrk, mcr, mer, mvl, mmt, mmg, mpool, mhot, mkid, mcnt, G.passB_demes, Tr.xfill);
  }
  addT(T, Tp);
}
