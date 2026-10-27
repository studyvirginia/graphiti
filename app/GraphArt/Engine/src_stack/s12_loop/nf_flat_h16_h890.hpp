#pragma once
#include "nf_flat_h16_h264.hpp"
#include "../fe_stack/nfc_core.hpp"
#include <mutex>
#include <atomic>
#include <chrono>
#include <algorithm>

namespace nf308 {
using namespace skn;

struct DumpRec { int32_t child; uint8_t side, n; uint8_t op[24]; uint8_t pad[2]; };
static_assert(sizeof(DumpRec) == 32, "DumpRec is 32 bytes");
inline bool g_h320_on = false;
inline std::mutex g_h320_mx;
inline std::vector<DumpRec> g_h320_recs;

struct DagB {
  int32_t lmap[MAXTE]; uint32_t lgen[MAXTE]; uint32_t lg;
  std::vector<PNode> out;
  DagB() : lg(0) { memset(lgen, 0, sizeof lgen); out.reserve(MAXNODES); }
  bool emit_te(const nfc::Ctx& C, int id) {
    if (lgen[id] == lg) return true;
    const TE& T = C.te[id];
    PNode nd{T.kind, T.op, 0, 0, T.leaf};
    if (T.kind == TE_NEG || T.kind == TE_UN) {
      if (!emit_te(C, T.a)) return false;
      nd.a = (uint16_t)lmap[T.a];
    } else if (T.kind == TE_BIN || (T.kind != TE_LIN && T.kind != TE_FIX && T.kind != TE_ONEV && T.kind != TE_LEAF)) {
      if (!emit_te(C, T.a)) return false;
      if (!emit_te(C, T.b)) return false;
      nd.a = (uint16_t)lmap[T.a]; nd.b = (uint16_t)lmap[T.b];
    }
    if (out.size() >= (size_t)MAXNODES) return false;
    lmap[id] = (int32_t)out.size(); lgen[id] = lg;
    out.push_back(nd);
    return true;
  }
  bool build(const nfc::Ctx& C, uint16_t* slot_root) {
    lg++;
    if (lg == 0) { memset(lgen, 0, sizeof lgen); lg = 1; }
    out.clear();
    for (int s = 0; s < C.nth; ++s) {
      if (!emit_te(C, C.slot_te[s])) return false;
      slot_root[s] = (uint16_t)lmap[C.slot_te[s]];
    }
    return true;
  }
};

inline bool g_h315a = false;
inline bool g_h315t = false;
inline std::atomic<uint64_t> g_t_pubcopy_h890{0};
inline std::atomic<uint64_t> g_t_nf{0}, g_t_cmp{0}, g_t_pub{0}, g_t_class{0}, g_t_sort{0}, g_t_eval{0}, g_t_out{0}, g_t_all{0}, g_t_rb{0}, g_n_rb{0}, g_w_obj{0}, g_w_objre{0}, g_w_vp{0}, g_w_vpre{0};
static inline uint64_t tnow() { return std::chrono::steady_clock::now().time_since_epoch().count(); }
inline std::atomic<uint64_t> g_h315_miss{0}, g_h315_inc{0}, g_h315_nodes{0}, g_h315_reused{0}, g_h315_groups{0}, g_h315_stage_hit{0};

struct IncSt {
  bool valid = false; int n = 0, done = 0, monoT = 1;
  uint8_t op[nfc::MAXN];
  struct Snap { int nte, nfac, nmon, nmp, nsum, nsp, nvp; } snap[nfc::MAXN];
  nfc::PR* val[nfc::MAXN]; int vl[nfc::MAXN];
};
static inline void erase_slot(uint64_t* tab, uint32_t gen, uint64_t h, int id) {
  uint32_t p = (uint32_t)(h & (TABN - 1));
  const uint64_t want = ((uint64_t)gen << 32) | (uint32_t)(id + 1);
  for (uint32_t it = 0; it < TABN; ++it) { if (tab[p] == want) { tab[p] = 0; return; } p = (p + 1) & (TABN - 1); }
  abort();
}
static inline void rollback_to(nfc::Ctx& C, const IncSt::Snap& s) {
  for (int j = C.nte - 1; j >= s.nte; --j) if (C.te[j].kind != TE_LEAF) erase_slot(C.tte, C.gen, C.te[j].h, j);
  for (int j = C.nfac - 1; j >= s.nfac; --j) erase_slot(C.tfac, C.gen, C.fac[j].h, j);
  for (int j = C.nmon - 1; j >= s.nmon; --j) erase_slot(C.tmon, C.gen, C.mon[j].h, j);
  for (int j = C.nsum - 1; j >= s.nsum; --j) erase_slot(C.tsum, C.gen, C.sum[j].h, j);
  C.nte = s.nte; C.nfac = s.nfac; C.nmon = s.nmon; C.nmp = s.nmp; C.nsum = s.nsum; C.nsp = s.nsp; C.nvp = s.nvp; C.ovf = false;
}

inline bool nf_tree_inc(nfc::Ctx& C, IncSt& I, const uint8_t* op, int n, nfc::NFOut* out, int* reused) {
  using namespace nfc;
  *reused = 0;
  if (n < 1 || n > MAXN) return false;
  int lc[MAXN], rc[MAXN], sp = 0;
  int stk[MAXN + 2]; int* st = stk + 2;
  bool tf[MAXN];
  static const int8_t ARITY[16] = {-1, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -1};
  for (int i = 0; i < n; ++i) {
    const uint8_t o = op[i];
    const int ar = o < 16 ? ARITY[o] : -1;
    if (ar < 0 || sp < ar) return false;
    const int Rr = st[sp - 1], Ll = st[sp - 2];
    lc[i] = ar == 2 ? Ll : ar == 1 ? Rr : -1;
    rc[i] = ar == 2 ? Rr : -1;
    tf[i] = ar == 0 ? (o == OP_CONST) : ar == 1 ? tf[Rr] : (tf[Ll] && tf[Rr]);
    sp += 1 - ar;
    st[sp - 1] = i;
  }
  if (sp != 1) return false;
  int r = 0;
  if (I.valid) { const int m = std::min(n, I.n); while (r < m && op[r] == I.op[r]) ++r; if (r > I.done) r = I.done; }
  I.valid = false;
  if (r == 0) {
    C.bump_gen();
    C.nte = 3; C.nfac = 1; C.nmon = 0; C.nmp = 0; C.nsum = 0; C.nsp = 0; C.nvp = 0; C.ovf = false;
    C.te[0] = TE{TE_LIN, 0, -1, -1, -1, H_LIN};
    C.te[1] = TE{TE_FIX, 0, -1, -1, -1, H_FIX};
    C.te[2] = TE{TE_ONEV, 0, -1, -1, -1, H_ONEV};
    C.fac[0] = FAC{FAC_T, 0, -1, H_T};
    C.nfac = 1;
    { uint64_t z = 0; C.mon_intern_h(nullptr, &z, 0, H_EMPTY); }
    I.monoT = C.mon_single(0);
  } else {
    if (g_h315t) { const uint64_t a = tnow(); const int before = (C.nte - I.snap[r - 1].nte) + (C.nfac - I.snap[r - 1].nfac) + (C.nmon - I.snap[r - 1].nmon) + (C.nsum - I.snap[r - 1].nsum);
      rollback_to(C, I.snap[r - 1]); g_t_rb += tnow() - a; g_n_rb += (uint64_t)before; }
    else rollback_to(C, I.snap[r - 1]);
  }
  *reused = r;
  const int monoT = I.monoT;
  memcpy(I.op, op, (size_t)n); I.n = n; I.done = r;
  PR** val = I.val; int* vl = I.vl;
  for (int i = r; i < n; ++i) {
    uint8_t o = op[i];
    if (C.nvp + MAXP2 >= MAXPAIR) { I.done = i; I.valid = true; return false; }
    PR* d = C.vp + C.nvp;
    int cnt = 1;
    bool bad = false;
    if (o == OP_CONST) d[0] = PR{C.te_leaf(i), 0};
    else if (o == OP_T) d[0] = PR{C.te_leaf(i), monoT};
    else if (tf[i]) {
      if (o >= OP_SIN) d[0] = PR{C.te_un(o, val[lc[i]][0].coef), 0};
      else d[0] = PR{C.te_bin(o, val[lc[i]][0].coef, val[rc[i]][0].coef), 0};
    } else if (o >= OP_SIN) {
      const PR* V = val[lc[i]]; int nv = vl[lc[i]];
      if (o == OP_SQ || o == OP_CUBE) cnt = op_pow(C, V, nv, o == OP_SQ ? 2 : 3, d);
      else if (o == OP_SIN || o == OP_COS) cnt = op_sincos(C, V, nv, o, d);
      else if (o == OP_EXP) cnt = op_exp(C, V, nv, d);
      else if (o == OP_LOG) cnt = op_log(C, V, nv, d);
      else cnt = op_root(C, V, nv, o, d);
    } else if (o == OP_ADD || o == OP_SUB) {
      cnt = op_addsub(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], o == OP_SUB, d);
    } else if (o == OP_MUL) {
      cnt = op_mul(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], d);
      if (cnt < 0) bad = true;
    } else {
      cnt = op_div(C, val[lc[i]], vl[lc[i]], val[rc[i]], vl[rc[i]], tf[rc[i]], d);
    }
    if (bad || C.ovf) { I.done = i; I.valid = true; return false; }
    val[i] = d; vl[i] = cnt; C.nvp += cnt;
    I.snap[i] = IncSt::Snap{C.nte, C.nfac, C.nmon, C.nmp, C.nsum, C.nsp, C.nvp};
    I.done = i + 1;
  }
  I.valid = true;
  if (g_h315t) { g_w_obj += (uint64_t)(C.nte + C.nfac + C.nmon + C.nsum); g_w_vp += (uint64_t)C.nvp;
    if (r > 0) { const IncSt::Snap& q = I.snap[r - 1]; g_w_objre += (uint64_t)(q.nte + q.nfac + q.nmon + q.nsum); g_w_vpre += (uint64_t)q.nvp; } }
  const PR* root = val[n - 1];
  int q = vl[n - 1];
  if (q > MAXQ) return false;
  nf_add_const(root, q);
  out->cols = root; out->q = q;
  return true;
}

static inline void publish_inplace_h890(nf16::Handle& H, uint64_t h, const NFMemo::Key& k, const nf16::BHdr& X, const void* w, size_t wbytes, const void* nd, size_t ndbytes) {
  nf16::Shared& S = H.S;
  const size_t nbody = wbytes + ndbytes;
  if (S.nent.load(std::memory_order_relaxed) >= nf16::MAXENT_H16) { S.drops.fetch_add(1); return; }
  const size_t bytes = (sizeof(nf16::BHdr) + nbody + 7) & ~(size_t)7;
  if (H.r_next + bytes > H.r_end) { const uint64_t blk = std::max<uint64_t>(nf16::RSV_BYTES, bytes); const uint64_t base = S.bump.fetch_add(blk); H.r_next = base; H.r_end = base + blk; }
  const uint64_t off = H.r_next; H.r_next += bytes;
  if (off + bytes > nf16::ARENA_BYTES) { S.drops.fetch_add(1); return; }
  memcpy(S.ar + off, &X, sizeof X);
  if (wbytes) memcpy(S.ar + off + sizeof X, w, wbytes);
  if (ndbytes) memcpy(S.ar + off + sizeof X + wbytes, nd, ndbytes);
  const uint32_t tag = (uint32_t)(h >> 32);
  const uint64_t packed = ((uint64_t)tag << 32) | (off + 1);
  uint32_t p = (uint32_t)(h & (nf16::TB - 1));
  for (;;) {
    uint64_t expected = 0;
    if (S.tab[p].compare_exchange_strong(expected, packed, std::memory_order_release, std::memory_order_acquire)) {
      S.nent.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    if ((uint32_t)(expected >> 32) == tag) {
      const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + ((uint32_t)expected - 1u));
      if (memcmp(&hd->key, &k, sizeof k) == 0) return;
    }
    p = (p + 1) & (nf16::TB - 1);
  }
}
template <class C> static inline void pf_cand_h890(const C* const* cands, int i, int lim) {
  if ((g_h890 & 2) && !(g_h890 & 32) && i + 16 < lim) { const char* p = (const char*)cands[i + 16]; for (size_t o = 0; o < sizeof(C) && o < 256; o += 64) __builtin_prefetch(p + o, 0, 3); }
}
struct St {
  bool inc = false; IncSt I; uint64_t c_miss = 0, c_inc = 0, c_nodes = 0, c_reused = 0;
  nfc::Ctx* C; nfc::Prog P; DagB DB;
  std::vector<char> body;
  bool chain_ok = false; uint64_t chain = 0;
  std::vector<DumpRec> dump; int cur_child = -1;
  St() : C(new nfc::Ctx()) {}
  ~St() { delete C; }
};

template <int MODE>
static inline bool miss_nfc(nf16::Handle& H, St& S, const NFMemo::Key& k, uint64_t h, const uint8_t* op, const float* c, int n) {
  H.misses++;
  nfc::NFOut out;
  bool good = false;
  uint8_t status = 0;
  S.chain_ok = false;
  bool nfok;
  const uint64_t T0 = g_h315t ? tnow() : 0;
  if (S.inc) { int ru; nfok = nf_tree_inc(*S.C, S.I, op, n, &out, &ru); S.c_miss++; S.c_nodes += (uint64_t)n; S.c_reused += (uint64_t)ru; S.c_inc += ru > 0; }
  else nfok = nfc::nf_tree(*S.C, op, c, n, &out);
  const uint64_t T1 = g_h315t ? tnow() : 0;
  if (nfok) { good = nfc::compile(*S.C, out, c, S.P); status = good ? 2 : 1; }
  const uint64_t T2 = g_h315t ? tnow() : 0;
  nf16::BHdr X; memset(&X, 0, sizeof X);
  X.touch = h156::g_touch;
  memcpy(&X.key, &k, sizeof k);
  X.status = status;
  if (good) {
    X.nabs = (uint16_t)S.P.nabs; X.nS = (uint16_t)S.P.nS; X.nV = (uint16_t)S.P.nV; X.Q = (uint8_t)S.P.Q;
    X.nth = (uint16_t)S.P.nth; X.wide = S.P.wide ? 1 : 0; X.h1 = 0; X.h2 = 0;
    if (S.DB.build(*S.C, X.slot_root)) {
      X.pn_len = (uint16_t)S.DB.out.size();
      X.chain = nf16::chain_prefix(S.P.w, (size_t)S.P.nabs, (uint32_t)(S.P.nabs + S.P.nth));
      S.chain_ok = true; S.chain = X.chain;
      const size_t nb = (size_t)S.P.nabs * 4 + S.DB.out.size() * sizeof(PNode);
      const uint64_t P0h890 = g_h315t ? tnow() : 0;
      if (MODE >= 2 && (g_h890 & 2)) {
        publish_inplace_h890(H, h, k, X, S.P.w, (size_t)S.P.nabs * 4, S.DB.out.data(), S.DB.out.size() * sizeof(PNode));
      } else if (MODE >= 2) {
        if (S.body.size() < nb) S.body.resize(nb + 4096);
        if (S.P.nabs) memcpy(S.body.data(), S.P.w, (size_t)S.P.nabs * 4);
        if (!S.DB.out.empty()) memcpy(S.body.data() + (size_t)S.P.nabs * 4, S.DB.out.data(), S.DB.out.size() * sizeof(PNode));
        nf16::publish(H, h, k, X, S.body.data(), nb);
      } else {
        std::vector<char> body(nb);
        if (S.P.nabs) memcpy(body.data(), S.P.w, (size_t)S.P.nabs * 4);
        if (!S.DB.out.empty()) memcpy(body.data() + (size_t)S.P.nabs * 4, S.DB.out.data(), S.DB.out.size() * sizeof(PNode));
        nf16::publish(H, h, k, X, body.data(), body.size());
      }
      if (g_h315t) g_t_pubcopy_h890 += tnow() - P0h890;
    } else H.bigskip++;
  } else nf16::publish(H, h, k, X, nullptr, 0);
  if (g_h315t) { const uint64_t T3 = tnow(); g_t_nf += T1 - T0; g_t_cmp += T2 - T1; g_t_pub += T3 - T2; }
  return good;
}

template <int MODE, class C>
static inline bool tree_seg_h308(nf16::Handle& H, St& S, const C& cd, int t, const NFMemo::Key& k, uint64_t hk, nf16::Block& B, nf16::Child& ch, uint64_t& r) {
  const size_t t0 = B.w.size();
  const int n = cd.len[t];
  uint64_t boff = 0;
  bool ok = false;
  if (nf16::find2(H.S, k, hk, boff)) {
    H.hits++;
    const nf16::BHdr* hd = (const nf16::BHdr*)(H.S.ar + boff);
    if (h156::g_on) h156::touch_u32(const_cast<uint32_t*>(&hd->touch), h156::g_touch);
    if (hd->status == 2 && hd->Q <= 16) {
      const uint32_t nabs = hd->nabs, nth = hd->nth, nw = nabs + nth;
      B.w.push_back(nw);
      const size_t r0 = B.w.size();
      B.w.resize(r0 + nw);
      const uint32_t* aw = (const uint32_t*)(H.S.ar + boff + sizeof(nf16::BHdr));
      if (nabs) memcpy(B.w.data() + r0, aw, sizeof(uint32_t) * nabs);
      const PNode* nd = (const PNode*)(H.S.ar + boff + sizeof(nf16::BHdr) + sizeof(uint32_t) * nabs);
      nf16::replay_to(nd, hd->pn_len, hd->slot_root, (int)nth, cd.c[t], B.w.data() + r0 + nabs);
      uint64_t h = hd->chain;
      for (uint32_t s = 0; s < nth; s++) h = nf16::sm(h ^ B.w[r0 + nabs + s]);
      ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw; r = h; ok = true;
    }
  } else {
    if (g_h320_on) { DumpRec d; std::memset(&d, 0, sizeof d); d.child = S.cur_child; d.side = (uint8_t)t; d.n = (uint8_t)n; std::memcpy(d.op, cd.op[t], std::min(n, 24)); S.dump.push_back(d); }
    const bool good = miss_nfc<MODE>(H, S, k, hk, cd.op[t], cd.c[t], n);
    if (good) {
      const uint32_t nw = (uint32_t)S.P.nw;
      B.w.push_back(nw);
      const size_t r0 = B.w.size();
      B.w.resize(r0 + nw);
      if (nw) memcpy(B.w.data() + r0, S.P.w, sizeof(uint32_t) * nw);
      ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw;
      if (MODE >= 2 && S.chain_ok) {
        uint64_t h = S.chain;
        for (int s = S.P.nabs; s < S.P.nw; s++) h = nf16::sm(h ^ S.P.w[s]);
        r = h;
      } else r = nf16::hash_words(B.w.data() + r0, nw);
      ok = S.P.Q <= 16;
    }
  }
  if (!ok) {
    H.fails++;
    B.w.resize(t0);
    const size_t s0 = B.w.size();
    B.w.push_back(0xFFFFFFFFu); B.w.push_back((uint32_t)cd.len[t]);
    for (int j = 0; j < cd.len[t]; j++) { uint32_t b; std::memcpy(&b, &cd.c[t][j], 4); B.w.push_back(cd.op[t][j]); B.w.push_back(b); }
    r = nf16::hash_words(B.w.data() + s0, B.w.size() - s0);
    ch.rlen[t] = 0;
  }
  return ok;
}

template <int MODE, class C>
static void chunk_h308(nf16::Handle& H, St& S, const C* const* cands, const uint8_t* valid, int b0, int b1, nf16::Block& B, nf16::Child* out) {
  B.w.clear();
  const int m = b1 - b0;
  struct TA { uint64_t t0; TA() : t0(g_h315t ? tnow() : 0) {} ~TA() { if (g_h315t) g_t_all += tnow() - t0; } } ta_;
  NFMemo::Key K[64][2];
  uint64_t HK[64][2];
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    pf_cand_h890(cands, i, b1);
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      nf16::make_key(cands[i]->op[t], cands[i]->len[t], K[j][t]);
      HK[j][t] = NFMemo::hkey(K[j][t]);
      __builtin_prefetch(&H.S.tab[(uint32_t)(HK[j][t] & (nf16::TB - 1))], 0, 3);
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      const uint64_t v = H.S.tab[(uint32_t)(HK[j][t] & (nf16::TB - 1))].load(std::memory_order_relaxed);
      if (v && (uint32_t)(v >> 32) == (uint32_t)(HK[j][t] >> 32)) {
        const char* bp = H.S.ar + ((uint32_t)v - 1u);
        __builtin_prefetch(bp, 0, 3); __builtin_prefetch(bp + 64, 0, 3); __builtin_prefetch(bp + 128, 0, 3);
      }
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      const uint64_t v = H.S.tab[(uint32_t)(HK[j][t] & (nf16::TB - 1))].load(std::memory_order_acquire);
      if (!v || (uint32_t)(v >> 32) != (uint32_t)(HK[j][t] >> 32)) continue;
      const char* bp = H.S.ar + ((uint32_t)v - 1u);
      const nf16::BHdr* hd = (const nf16::BHdr*)bp;
      const size_t nb = (size_t)hd->nabs * 4 + (size_t)hd->pn_len * sizeof(PNode);
      const size_t lim = std::min<size_t>(nb, 768);
      for (size_t o = 192; o < sizeof(nf16::BHdr) + lim; o += 64) __builtin_prefetch(bp + o, 0, 3);
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    nf16::Child& ch = out[i];
    ch = nf16::Child();
    if (!valid[i]) { ch.koff = (uint32_t)B.w.size(); continue; }
    const C& cd = *cands[i];
    ch.koff = (uint32_t)B.w.size();
    S.cur_child = i;
    uint64_t r[2] = {0, 0};
    ch.ok[0] = tree_seg_h308<MODE>(H, S, cd, 0, K[j][0], HK[j][0], B, ch, r[0]) ? 1 : 0;
    ch.ok[1] = tree_seg_h308<MODE>(H, S, cd, 1, K[j][1], HK[j][1], B, ch, r[1]) ? 1 : 0;
    ch.klen = (uint32_t)(B.w.size() - ch.koff);
    ch.r0 = r[0]; ch.r1 = r[1];
    ch.kh = nf16::sm(r[0] ^ nf16::sm(r[1] + ch.klen));
  }
  if (g_h320_on && !S.dump.empty()) { std::lock_guard<std::mutex> lk(g_h320_mx); g_h320_recs.insert(g_h320_recs.end(), S.dump.begin(), S.dump.end()); S.dump.clear(); }
}

struct MissKey { uint64_t k[3]; int32_t j; uint8_t t, n; };
struct StageRec { uint8_t state; uint32_t off, nw; uint64_t r; };
struct GroupBuf {
  std::vector<MissKey> mk; std::vector<uint32_t> words; StageRec sr[512][2];
  NFMemo::Key K[512][2]; uint64_t HK[512][2];
};

template <int MODE, class C>
static void group_h315(nf16::Handle& H, St& S, const C* const* cands, const uint8_t* valid, int a0, int a1, nf16::Block* blk0, nf16::Child* out) {
  static thread_local GroupBuf* GBp = nullptr;
  if (!GBp) GBp = new GroupBuf();
  GroupBuf& G = *GBp;
  const int m = a1 - a0;
  const uint64_t U0 = g_h315t ? tnow() : 0;
  G.mk.clear(); G.words.clear();
  for (int j = 0; j < m; j++) { G.sr[j][0].state = 0; G.sr[j][1].state = 0; }
  for (int j = 0; j < m; j++) {
    const int i = a0 + j;
    pf_cand_h890(cands, i, a1);
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      nf16::make_key(cands[i]->op[t], cands[i]->len[t], G.K[j][t]);
      G.HK[j][t] = NFMemo::hkey(G.K[j][t]);
      __builtin_prefetch(&H.S.tab[(uint32_t)(G.HK[j][t] & (nf16::TB - 1))], 0, 3);
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = a0 + j;
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      uint64_t boff;
      if (nf16::find2(H.S, G.K[j][t], G.HK[j][t], boff)) {
        const char* bp = H.S.ar + boff;
        __builtin_prefetch(bp, 0, 3); __builtin_prefetch(bp + 64, 0, 3); __builtin_prefetch(bp + 128, 0, 3);
        const nf16::BHdr* hd = (const nf16::BHdr*)bp;
        const size_t nb = (size_t)hd->nabs * 4 + (size_t)hd->pn_len * sizeof(PNode);
        const size_t lim = std::min<size_t>(nb, 768);
        for (size_t o = 192; o < sizeof(nf16::BHdr) + lim; o += 64) __builtin_prefetch(bp + o, 0, 3);
        continue;
      }
      const int n = cands[i]->len[t];
      MissKey mk; mk.j = j; mk.t = (uint8_t)t; mk.n = (uint8_t)n;
      uint8_t b[24]; std::memset(b, 0, 24); std::memcpy(b, cands[i]->op[t], std::min(n, 24));
      for (int w = 0; w < 3; w++) { uint64_t v = 0; for (int q = 0; q < 8; q++) v = (v << 8) | b[w * 8 + q]; mk.k[w] = v; }
      G.mk.push_back(mk);
    }
  }
  const uint64_t U1 = g_h315t ? tnow() : 0;
  std::sort(G.mk.begin(), G.mk.end(), [](const MissKey& x, const MissKey& y) {
    if (x.k[0] != y.k[0]) return x.k[0] < y.k[0];
    if (x.k[1] != y.k[1]) return x.k[1] < y.k[1];
    if (x.k[2] != y.k[2]) return x.k[2] < y.k[2];
    if (x.n != y.n) return x.n < y.n;
    return x.j != y.j ? x.j < y.j : x.t < y.t;
  });
  const uint64_t U2 = g_h315t ? tnow() : 0;
  S.inc = true;
  for (const MissKey& mk : G.mk) {
    const int j = mk.j, t = mk.t, i = a0 + j;
    uint64_t boff;
    if (nf16::find2(H.S, G.K[j][t], G.HK[j][t], boff)) { g_h315_stage_hit.fetch_add(1, std::memory_order_relaxed); continue; }
    const C& cd = *cands[i];
    const int n = cd.len[t];
    if (g_h320_on) { DumpRec d; std::memset(&d, 0, sizeof d); d.child = i; d.side = (uint8_t)t; d.n = (uint8_t)n; std::memcpy(d.op, cd.op[t], std::min(n, 24)); S.dump.push_back(d); }
    const bool good = miss_nfc<MODE>(H, S, G.K[j][t], G.HK[j][t], cd.op[t], cd.c[t], n);
    StageRec& R = G.sr[j][t];
    R.state = 2;
    if (good && S.P.Q <= 16) {
      const uint32_t nw = (uint32_t)S.P.nw;
      R.off = (uint32_t)G.words.size(); R.nw = nw;
      G.words.resize(G.words.size() + nw);
      if (nw) std::memcpy(G.words.data() + R.off, S.P.w, sizeof(uint32_t) * nw);
      if (MODE >= 2 && S.chain_ok) {
        uint64_t h = S.chain;
        for (int s = S.P.nabs; s < S.P.nw; s++) h = nf16::sm(h ^ S.P.w[s]);
        R.r = h;
      } else R.r = nf16::hash_words(S.P.w, nw);
      R.state = 1;
    }
  }
  S.inc = false;
  const uint64_t U3 = g_h315t ? tnow() : 0;
  for (int c0 = 0; c0 < m; c0 += 64) {
    nf16::Block& B = blk0[c0 / 64];
    B.w.clear();
    for (int j = c0; j < std::min(m, c0 + 64); j++) {
      const int i = a0 + j;
      nf16::Child& ch = out[i];
      ch = nf16::Child();
      if (!valid[i]) { ch.koff = (uint32_t)B.w.size(); continue; }
      const C& cd = *cands[i];
      ch.koff = (uint32_t)B.w.size();
      S.cur_child = i;
      uint64_t r[2] = {0, 0};
      for (int t = 0; t < 2; t++) {
        const StageRec& R = G.sr[j][t];
        bool ok;
        if (R.state == 1) {
          const size_t r0 = B.w.size() + 1;
          B.w.push_back(R.nw);
          B.w.resize(r0 + R.nw);
          if (R.nw) std::memcpy(B.w.data() + r0, G.words.data() + R.off, sizeof(uint32_t) * R.nw);
          ch.roff[t] = (uint32_t)r0; ch.rlen[t] = R.nw; r[t] = R.r; ok = true;
        } else if (R.state == 2) {
          H.fails++;
          const size_t s0 = B.w.size();
          B.w.push_back(0xFFFFFFFFu); B.w.push_back((uint32_t)cd.len[t]);
          for (int q = 0; q < cd.len[t]; q++) { uint32_t b; std::memcpy(&b, &cd.c[t][q], 4); B.w.push_back(cd.op[t][q]); B.w.push_back(b); }
          r[t] = nf16::hash_words(B.w.data() + s0, B.w.size() - s0);
          ch.rlen[t] = 0; ok = false;
        } else {
          ok = tree_seg_h308<MODE>(H, S, cd, t, G.K[j][t], G.HK[j][t], B, ch, r[t]);
        }
        ch.ok[t] = ok ? 1 : 0;
      }
      ch.klen = (uint32_t)(B.w.size() - ch.koff);
      ch.r0 = r[0]; ch.r1 = r[1];
      ch.kh = nf16::sm(r[0] ^ nf16::sm(r[1] + ch.klen));
    }
  }
  if (g_h315t) { const uint64_t U4 = tnow(); g_t_class += U1 - U0; g_t_sort += U2 - U1; g_t_eval += U3 - U2; g_t_out += U4 - U3; g_t_all += U4 - U0; }
  g_h315_miss.fetch_add(S.c_miss, std::memory_order_relaxed); g_h315_inc.fetch_add(S.c_inc, std::memory_order_relaxed);
  g_h315_nodes.fetch_add(S.c_nodes, std::memory_order_relaxed); g_h315_reused.fetch_add(S.c_reused, std::memory_order_relaxed);
  g_h315_groups.fetch_add(1, std::memory_order_relaxed);
  S.c_miss = S.c_inc = S.c_nodes = S.c_reused = 0;
  if (g_h320_on && !S.dump.empty()) { std::lock_guard<std::mutex> lk(g_h320_mx); g_h320_recs.insert(g_h320_recs.end(), S.dump.begin(), S.dump.end()); S.dump.clear(); }
}

}
