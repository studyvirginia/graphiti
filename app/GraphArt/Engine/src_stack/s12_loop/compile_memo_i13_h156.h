#pragma once
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "compile_cache_h14_i13_h156.h"

namespace cmi6 {

struct TreeEnc {
  uint8_t s_ok = 0, g_ok = 0;
  uint32_t sg[36], gg[36];
  uint32_t ghw = 0, gper = 0;
  std::vector<uint32_t> sp4, gp4;
  const uint32_t* sp4p = nullptr; const uint32_t* gp4p = nullptr; uint32_t nsp4p = 0, ngp4p = 0;
  size_t nsp() const { return sp4p ? (size_t)nsp4p : sp4.size(); }
  size_t ngp() const { return gp4p ? (size_t)ngp4p : gp4.size(); }
  const uint32_t* psp() const { return sp4p ? sp4p : sp4.data(); }
  const uint32_t* pgp() const { return gp4p ? gp4p : gp4.data(); }
  std::vector<uint32_t> tk; const uint32_t* tkp = nullptr; uint32_t ntkp = 0;
  size_t ntkw() const { return tkp ? (size_t)ntkp : tk.size(); }
  const uint32_t* ptk() const { return tkp ? tkp : tk.data(); }
  const uint32_t* gtp = nullptr;
};

static inline void encode_tree_i6(const Tree& tr, TreeEnc& te) {
  std::memset(te.sg, 0, sizeof te.sg); std::memset(te.gg, 0, sizeof te.gg);
  { std::vector<uint32_t> p4; uint32_t sg[36]; std::memset(sg, 0, sizeof sg);
    if (zconvS(tr, p4, sg, tr.nth == 0 ? 2 : 16)) { te.s_ok = 1; std::memcpy(te.sg, sg, sizeof sg); te.sp4 = std::move(p4); } }
  if (tr.nth >= 1 && tr.nth <= GN_TH_MAX) {
    std::vector<uint32_t> p4; uint32_t sg[36]; uint32_t hw = 0; size_t per = 0; std::memset(sg, 0, sizeof sg);
    if (zgconv(tr, tr.nth <= 4 ? 1 + (tr.nth + 1) / 2 : 3, p4, sg, hw, per)) {
      te.g_ok = 1; std::memcpy(te.gg, sg, sizeof sg); te.ghw = hw; te.gper = (uint32_t)per; te.gp4 = std::move(p4); } }
}

struct Slot {
  uint32_t st;
  uint32_t ghw, gper, sp_off, nsp4, gp_off, ngp4;
  uint8_t s_ok, g_ok;
  uint8_t cmode, Qs;
  uint32_t sg[36], gg[36];
};
static_assert(sizeof(Slot) == 32 + 288, "Slot layout (H313a: cmode/Qs take the former padding)");
inline bool g_a3 = false;
inline bool g_b1 = false;
inline bool g_a1 = false;
static inline bool tab_canon(const uint32_t* t, int Q) {
  for (int q = Q + 1; q < 18; q++) if (t[q] != 0) return false;
  for (int q = Q; q < 18; q++) if (t[18 + q] != t[35]) return false;
  return true;
}
static inline void tab_pack(uint32_t* d, const uint32_t* t, int Q) {
  for (int q = 0; q <= Q; q++) d[q] = t[q];
  for (int q = 0; q < Q; q++) d[Q + 1 + q] = t[18 + q];
  d[2 * Q + 1] = t[35];
}
static inline void tab_unpack(uint32_t* t, const uint32_t* d, int Q) {
  for (int q = 0; q <= Q; q++) t[q] = d[q];
  for (int q = Q + 1; q < 18; q++) t[q] = 0;
  for (int q = 0; q < Q; q++) t[18 + q] = d[Q + 1 + q];
#ifdef H313A_NEG3
  const uint32_t tail = d[2 * Q + 1] + 1;
#else
  const uint32_t tail = d[2 * Q + 1];
#endif
  for (int q = Q; q < 18; q++) t[18 + q] = tail;
}
inline bool g_gt = false;
inline bool g_gt_prof = false;
typedef bool (*GtFn)(const Tree&, std::vector<uint32_t>&, uint32_t*, uint32_t&);
inline GtFn g_gtfn = nullptr;
inline std::atomic<uint64_t> g_gt_ns{0}, g_gt_n{0}, g_gt_w{0}, g_gt_fail{0};
static constexpr uint64_t GT_HDR = 40;
inline bool g_h913 = false, g_h913_chk = false;
inline std::atomic<uint64_t> g_h913_n{0}, g_h913_w{0}, g_h913_cmc[4], g_h913_q[19];
static constexpr size_t H913_SLOT = 32;
static inline size_t slot_hdr_bytes() { return g_h913 ? H913_SLOT : sizeof(Slot); }
static inline uint32_t tabw(uint8_t cm, uint8_t Q) { return ((cm & 1) ? 2u * Q + 2u : 36u) + ((cm & 2) ? 2u * Q + 2u : 36u); }
static inline const uint32_t* tabs_ptr(const Slot& s, const uint32_t* ar) {
  uint64_t o = (uint64_t)s.gp_off + s.ngp4;
  if (s.cmode & 4) o += GT_HDR + ar[o];
  return ar + o;
}
static inline uint8_t h913_cm(const TreeEnc& te, int Q) {
  uint8_t cm = 0;
  if ((unsigned)Q <= 17u) { if (tab_canon(te.sg, Q)) cm |= 1; if (tab_canon(te.gg, Q)) cm |= 2; }
  return cm;
}
static inline void h913_get(const uint32_t* tp, uint8_t cm, uint8_t Q, uint32_t* sg, uint32_t* gg) {
  if (cm & 1) { tab_unpack(sg, tp, Q); tp += 2u * Q + 2u; } else { std::memcpy(sg, tp, 36 * sizeof(uint32_t)); tp += 36; }
  if (cm & 2) tab_unpack(gg, tp, Q); else std::memcpy(gg, tp, 36 * sizeof(uint32_t));
}
static inline void h913_put(uint32_t* tp, const TreeEnc& te, uint8_t cm, int Q) {
  uint32_t* t0 = tp;
  if (cm & 1) { tab_pack(tp, te.sg, Q); tp += 2 * Q + 2; } else { std::memcpy(tp, te.sg, 36 * sizeof(uint32_t)); tp += 36; }
  if (cm & 2) tab_pack(tp, te.gg, Q); else std::memcpy(tp, te.gg, 36 * sizeof(uint32_t));
  g_h913_n.fetch_add(1, std::memory_order_relaxed); g_h913_w.fetch_add(tabw(cm, (uint8_t)Q), std::memory_order_relaxed);
  g_h913_cmc[cm & 3].fetch_add(1, std::memory_order_relaxed); g_h913_q[(unsigned)Q <= 17u ? Q : 18].fetch_add(1, std::memory_order_relaxed);
  if (g_h913_chk) {
    uint32_t a[36], b[36]; h913_get(t0, cm, (uint8_t)Q, a, b);
    if (std::memcmp(a, te.sg, sizeof a) != 0 || std::memcmp(b, te.gg, sizeof b) != 0) { fprintf(stderr, "H913 self-check FAILED (cm %d Q %d)\n", (int)cm, Q); abort(); }
  }
}
struct GtScr { std::vector<uint32_t> p4; uint32_t sgl[36]; uint32_t hw = 0; bool ok = false; };
static inline uint64_t gt_make(const Tree& T, GtScr& g) {
  const auto t0 = g_gt_prof ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
  std::memset(g.sgl, 0, sizeof g.sgl); g.hw = 0; g.p4.clear();
  g.ok = g_gtfn(T, g.p4, g.sgl, g.hw);
  if (!g.ok) g.p4.clear();
  if (g_gt_prof) { g_gt_ns.fetch_add((uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count(), std::memory_order_relaxed);
    g_gt_n.fetch_add(1, std::memory_order_relaxed); g_gt_w.fetch_add(GT_HDR + g.p4.size(), std::memory_order_relaxed); if (!g.ok) g_gt_fail.fetch_add(1, std::memory_order_relaxed); }
  return GT_HDR + g.p4.size();
}
static inline void gt_put(uint32_t* d, const GtScr& g) {
  d[0] = (uint32_t)g.p4.size(); d[1] = g.ok ? g.hw : 0u; d[2] = g.ok ? 1u : 0u; d[3] = 0u;
  std::memcpy(d + 4, g.sgl, sizeof g.sgl);
  if (!g.p4.empty()) std::memcpy(d + GT_HDR, g.p4.data(), g.p4.size() * 4);
}
static constexpr size_t ARENA_WORDS = (size_t)1 << 27;

struct Memo {
  size_t stride;
  Slot* slot;
  uint32_t* arena;
  std::atomic<uint64_t> bump;
  std::atomic<uint32_t> nst{0};
  Slot& at(size_t i) { return *(Slot*)((char*)slot + i * stride); }
  const Slot& at(size_t i) const { return *(const Slot*)((const char*)slot + i * stride); }
  Memo() : stride(g_h913 ? H913_SLOT : sizeof(Slot)), slot((Slot*)cc14::ios_bigalloc(cc14::Shared::MAXENT, stride, "slot")), arena((uint32_t*)cc14::ios_bigalloc(ARENA_WORDS, sizeof(uint32_t), "arena")), bump(0) {}
  ~Memo() {}
  Memo(const Memo&) = delete;
  Memo& operator=(const Memo&) = delete;
};

struct alignas(128) Handle {
  cc14::Handle C;
  Memo& M;
  long enc_hit = 0, enc_miss = 0, enc_local = 0;
  Handle(cc14::Shared& s, Memo& m) : C(s), M(m) {}
};

static inline double since_ms(std::chrono::steady_clock::time_point t) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}
static inline void timed_encode(const Tree& T, TreeEnc& te, double& enc_ms) {
  const auto t0 = std::chrono::steady_clock::now();
  encode_tree_i6(T, te);
  enc_ms += since_ms(t0);
}

static inline bool locate(const cc14::Shared& S, const uint32_t* w, size_t nk, uint64_t h, uint32_t& idx) {
  uint32_t p = (uint32_t)(h & (cc14::Shared::TB - 1));
  for (;;) {
    const uint32_t v = S.tab[p].load(std::memory_order_acquire);
    if (v == 0) return false;
    const cc14::Ent& E = S.ent[v - 1];
    if (E.h == h && E.nkey == nk && std::memcmp(S.key + E.key_off, w, sizeof(uint32_t) * nk) == 0) { idx = v - 1; return true; }
    p = (p + 1) & (cc14::Shared::TB - 1);
  }
}

static inline uint64_t slot_words(const Slot& s, const uint32_t* ar) {
  uint64_t w = (uint64_t)s.nsp4 + s.ngp4;
  if ((s.cmode & 4) && ar) w += GT_HDR + ar[(uint64_t)s.gp_off + s.ngp4];
  if (g_h913) w += tabw(s.cmode, s.Qs);
  return w;
}
static inline void load_enc(const Slot& s, const uint32_t* ar, TreeEnc& te) {
  te.s_ok = s.s_ok; te.g_ok = s.g_ok;
  te.gtp = (g_gt && (s.cmode & 4)) ? ar + s.gp_off + s.ngp4 : nullptr;
  if (g_h913) h913_get(tabs_ptr(s, ar), s.cmode, s.Qs, te.sg, te.gg);
  else {
  if (s.cmode & 1) tab_unpack(te.sg, s.sg, s.Qs); else std::memcpy(te.sg, s.sg, sizeof te.sg);
  if (s.cmode & 2) tab_unpack(te.gg, s.gg, s.Qs); else std::memcpy(te.gg, s.gg, sizeof te.gg);
  }
  te.ghw = s.ghw; te.gper = s.gper;
  if (g_b1) {
    if (g_a1) { te.tkp = ar + s.sp_off; te.ntkp = s.nsp4; } else te.tk.assign(ar + s.sp_off, ar + s.sp_off + s.nsp4);
    return;
  }
#ifdef H313A_NEG1
  if (g_a1) { te.sp4p = ar + s.sp_off; te.nsp4p = s.nsp4 ? s.nsp4 - 1 : 0; te.gp4p = ar + s.gp_off; te.ngp4p = s.ngp4; }
#else
  if (g_a1) { te.sp4p = ar + s.sp_off; te.nsp4p = s.nsp4; te.gp4p = ar + s.gp_off; te.ngp4p = s.ngp4; }
#endif
  else {
    te.sp4.assign(ar + s.sp_off, ar + s.sp_off + s.nsp4);
    te.gp4.assign(ar + s.gp_off, ar + s.gp_off + s.ngp4);
  }
}

static inline bool compile_enc(Handle& H, const View& V, Tree& T, TreeEnc& te, double& enc_ms) {
  if (!cc14::compile_cached(H.C, V, T)) return false;
  const size_t nk = (size_t)(V.theta - V.w);
  const uint64_t h = cc14::hkey(V.w, nk);
  uint32_t idx = 0;
  if (!locate(H.C.S, V.w, nk, h, idx)) { H.enc_local++; timed_encode(T, te, enc_ms); return true; }
  Slot& s = H.M.at(idx);
  const uint32_t st = __atomic_load_n(&s.st, __ATOMIC_ACQUIRE);
  if (st == 2) { load_enc(s, H.M.arena, te); H.enc_hit++; return true; }
  timed_encode(T, te, enc_ms);
  H.enc_miss++;
  uint32_t zero = 0;
  if (st == 0 && __atomic_compare_exchange_n(&s.st, &zero, 1u, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
    const uint64_t n4 = te.sp4.size(), ng = te.gp4.size(), need = n4 + ng;
    const uint64_t off = H.M.bump.fetch_add(need);
    if (off + need <= ARENA_WORDS) {
      if (n4) std::memcpy(H.M.arena + off, te.sp4.data(), n4 * sizeof(uint32_t));
      if (ng) std::memcpy(H.M.arena + off + n4, te.gp4.data(), ng * sizeof(uint32_t));
      if (g_h913) { fprintf(stderr, "H913: store via cmi6::compile_enc is not supported\n"); abort(); }
      s.s_ok = te.s_ok; s.g_ok = te.g_ok; s.ghw = te.ghw; s.gper = te.gper;
      std::memcpy(s.sg, te.sg, sizeof s.sg); std::memcpy(s.gg, te.gg, sizeof s.gg);
      s.sp_off = (uint32_t)off; s.nsp4 = (uint32_t)n4; s.gp_off = (uint32_t)(off + n4); s.ngp4 = (uint32_t)ng;
      __atomic_store_n(&s.st, 2u, __ATOMIC_RELEASE);
    } else {
      __atomic_store_n(&s.st, 3u, __ATOMIC_RELEASE);
    }
  }
  return true;
}

static inline void prefetch_batch(const cc14::Shared& S, const Memo& M, const uint32_t* const* w, const size_t* nk, int n,
                                  const uint64_t* hin = nullptr) {
  uint64_t hh[32]; uint32_t pp[32]; uint32_t ix[32];
  for (int i = 0; i < n; i++) { hh[i] = hin ? hin[i] : cc14::hkey(w[i], nk[i]); pp[i] = (uint32_t)(hh[i] & (cc14::Shared::TB - 1)); __builtin_prefetch(&S.tab[pp[i]], 0, 3); }
  for (int i = 0; i < n; i++) {
    const uint32_t v = S.tab[pp[i]].load(std::memory_order_acquire);
    ix[i] = v;
    if (v) { __builtin_prefetch(&S.ent[v - 1], 0, 3); __builtin_prefetch((const char*)&S.ent[v - 1] + 64, 0, 3); }
  }
  for (int i = 0; i < n; i++) {
    if (!ix[i]) continue;
    const cc14::Ent& E = S.ent[ix[i] - 1];
    if (E.h != hh[i] || E.nkey != nk[i]) continue;
    __builtin_prefetch(S.key + E.key_off, 0, 3); __builtin_prefetch(S.prog + E.prog_off, 0, 3); __builtin_prefetch(S.sp + E.sp_off, 0, 3);
    __builtin_prefetch(&M.at(ix[i] - 1), 0, 3);
    if (!g_h913) { __builtin_prefetch((const char*)&M.at(ix[i] - 1) + 64, 0, 3);
    __builtin_prefetch((const char*)&M.at(ix[i] - 1) + 128, 0, 3); __builtin_prefetch((const char*)&M.at(ix[i] - 1) + 192, 0, 3); }
  }
  for (int i = 0; i < n; i++) {
    if (!ix[i]) continue;
    const cc14::Ent& E = S.ent[ix[i] - 1];
    if (E.h != hh[i] || E.nkey != nk[i]) continue;
    const Slot& sl = M.at(ix[i] - 1);
    if (__atomic_load_n(&sl.st, __ATOMIC_ACQUIRE) == 2) { __builtin_prefetch(M.arena + sl.sp_off, 0, 3); __builtin_prefetch(M.arena + sl.gp_off, 0, 3); }
  }
}

}
