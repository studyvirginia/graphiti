#include <mutex>
#include <string>
#include <map>
#pragma once
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <sys/mman.h>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include "h156_flags.h"
#include "compile_tree_h14.h"

namespace cc14 {

struct Ent {
  uint64_t h;
  uint32_t key_off, nkey, prog_off, nprog, sp_off, nsp;
  int32_t Q, nth, nS, nslot;
  uint8_t ok;  uint32_t touch;
  uint16_t segs[36];
};

static inline uint64_t hkey(const uint32_t* w, size_t n) {
  uint64_t h = 0x243F6A8885A308D3ull ^ (uint64_t)n;
  for (size_t i = 0; i < n; i++) { h = (h ^ w[i]) * 0x9E3779B97F4A7C15ull; h ^= h >> 31; }
  h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 29;
  return h;
}

inline std::mutex g_bigalloc_m; inline std::map<std::string, std::pair<void*, size_t>> g_bigalloc_prev;
static inline void* ios_bigalloc(size_t n, size_t sz, const char* what){
  size_t b=n*sz;
  { std::lock_guard<std::mutex> lk(g_bigalloc_m); auto it = g_bigalloc_prev.find(what);
    if (it != g_bigalloc_prev.end()) { munmap(it->second.first, it->second.second); g_bigalloc_prev.erase(it); } }
  void* p=mmap(nullptr,b,PROT_READ|PROT_WRITE,MAP_ANON|MAP_PRIVATE|MAP_NORESERVE,-1,0);
  if(p==MAP_FAILED){ fprintf(stderr,"IOS_BIGALLOC FAIL %s bytes=%zu errno=%d\n",what,b,errno); fflush(stderr); return nullptr; }
  { std::lock_guard<std::mutex> lk(g_bigalloc_m); g_bigalloc_prev[what] = {p, b}; }
  return p; }
struct Shared {
  static constexpr uint32_t TB = 1u << 24;
  static constexpr uint32_t MAXENT = 1u << 21;
  static constexpr size_t MAXKEY = (size_t)1 << 27;
  static constexpr size_t MAXPROG = (size_t)1 << 26;
  static constexpr size_t MAXSP = (size_t)1 << 26;
  std::atomic<uint32_t>* tab;
  Ent* ent;
  uint32_t* key;
  uint32_t* prog;
  uint16_t* sp;
  std::atomic<uint32_t> nent, nkey, nprog, nsp;
  Shared() : tab(new std::atomic<uint32_t>[TB]()) {
    ent = (Ent*)ios_bigalloc(MAXENT, sizeof(Ent), "ent");
    key = (uint32_t*)ios_bigalloc(MAXKEY, sizeof(uint32_t), "key");
    prog = (uint32_t*)ios_bigalloc(MAXPROG, sizeof(uint32_t), "prog");
    sp = (uint16_t*)ios_bigalloc(MAXSP, sizeof(uint16_t), "sp");
    nent.store(0); nkey.store(0); nprog.store(0); nsp.store(0);
  }
  ~Shared() { delete[] tab;  }
  void prefault(size_t bytes) {
    memset(key, 0, std::min(bytes, sizeof(uint32_t) * MAXKEY));
    memset(prog, 0, std::min(bytes, sizeof(uint32_t) * MAXPROG));
    memset(sp, 0, std::min(bytes, sizeof(uint16_t) * MAXSP));
    memset(ent, 0, std::min(bytes, sizeof(Ent) * (size_t)MAXENT));
  }
  Shared(const Shared&) = delete;
  Shared& operator=(const Shared&) = delete;
};

static constexpr uint32_t RSV_ENT = 256, RSV_WORDS = 1u << 14;
static inline bool take_block(std::atomic<uint32_t>& ctr, uint32_t& nxt, uint32_t& end, uint32_t n, uint32_t blk, size_t cap, uint32_t& off) {
  if (end - nxt < n) { const uint32_t want = std::max(blk, n); const uint32_t base = ctr.fetch_add(want, std::memory_order_relaxed); nxt = base; end = base + want; }
  off = nxt;
  if ((size_t)off + n > cap) return false;
  nxt += n;
  return true;
}

struct alignas(128) Handle {
  Shared& S;
  long hits, misses, neg_hits, drops;
  uint32_t r_ent_next = 0, r_ent_end = 0, r_key_next = 0, r_key_end = 0, r_prog_next = 0, r_prog_end = 0, r_sp_next = 0, r_sp_end = 0;
  explicit Handle(Shared& s) : S(s), hits(0), misses(0), neg_hits(0), drops(0) {}

  void publish(uint64_t h, Ent X, const uint32_t* kw, size_t nk, const uint32_t* pw, size_t np, const uint16_t* spw, size_t ns) {
    uint32_t idx, koff, poff, soff;
    if (!take_block(S.nent, r_ent_next, r_ent_end, 1, RSV_ENT, Shared::MAXENT, idx)) { ++drops; return; }
    if (!take_block(S.nkey, r_key_next, r_key_end, (uint32_t)nk, RSV_WORDS, Shared::MAXKEY, koff)) { ++drops; return; }
    if (!take_block(S.nprog, r_prog_next, r_prog_end, (uint32_t)np, RSV_WORDS, Shared::MAXPROG, poff)) { ++drops; return; }
    if (!take_block(S.nsp, r_sp_next, r_sp_end, (uint32_t)ns, RSV_WORDS, Shared::MAXSP, soff)) { ++drops; return; }
    X.touch = h156::g_touch; X.key_off = koff; X.nkey = (uint32_t)nk; X.prog_off = poff; X.nprog = (uint32_t)np; X.sp_off = soff; X.nsp = (uint32_t)ns;
    Ent& D = S.ent[idx];
    D = X;
    memcpy(S.key + koff, kw, sizeof(uint32_t) * nk);
    if (np) memcpy(S.prog + poff, pw, sizeof(uint32_t) * np);
    if (ns) memcpy(S.sp + soff, spw, sizeof(uint16_t) * ns);
    uint32_t p = (uint32_t)(h & (Shared::TB - 1));
    for (;;) {
      uint32_t expected = 0;
      if (S.tab[p].compare_exchange_strong(expected, idx + 1, std::memory_order_release, std::memory_order_acquire)) return;
      const Ent& E = S.ent[expected - 1];
      if (E.h == h && E.nkey == nk && memcmp(S.key + E.key_off, kw, sizeof(uint32_t) * nk) == 0) return;
      p = (p + 1) & (Shared::TB - 1);
    }
  }
};

static inline bool compile_cached(Handle& H, const View& V, Tree& T) {
  const size_t nk = (size_t)(V.theta - V.w);
  const uint64_t h = hkey(V.w, nk);
  uint32_t p = (uint32_t)(h & (Shared::TB - 1));
  for (;;) {
    const uint32_t v = H.S.tab[p].load(std::memory_order_acquire);
    if (v == 0) break;
    const Ent& E = H.S.ent[v - 1];
    if (E.h == h && E.nkey == nk && memcmp(H.S.key + E.key_off, V.w, sizeof(uint32_t) * nk) == 0) {
      ++H.hits;
      if (!E.ok) { ++H.neg_hits; return false; }
      T.Q = E.Q; T.nth = E.nth; T.nS = E.nS; T.nslot = E.nslot;
      T.prog.assign(H.S.prog + E.prog_off, H.S.prog + E.prog_off + E.nprog);
      T.sp.assign(H.S.sp + E.sp_off, H.S.sp + E.sp_off + E.nsp);
      memcpy(T.segs, E.segs, sizeof T.segs);
      T.th.resize((size_t)V.nth);
      for (int s = 0; s < V.nth; s++) T.th[s] = V.th(s);
      return true;
    }
    p = (p + 1) & (Shared::TB - 1);
  }
  ++H.misses;
  const bool ok = compile_tree_corder_h14(V, T);
  Ent X; memset(&X, 0, sizeof X);
  X.h = h; X.ok = ok ? 1 : 0;
  if (ok) {
    X.Q = T.Q; X.nth = T.nth; X.nS = T.nS; X.nslot = T.nslot;
    memcpy(X.segs, T.segs, sizeof T.segs);
    H.publish(h, X, V.w, nk, T.prog.data(), T.prog.size(), T.sp.data(), T.sp.size());
  } else {
    H.publish(h, X, V.w, nk, nullptr, 0, nullptr, 0);
  }
  return ok;
}

}
