#pragma once
#include "../s02_nf/sk_nf_h2b.hpp"
#include "h156_flags.h"
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>
#include <cstdio>
#include <vector>

namespace nf16 {
using namespace skn;

constexpr uint32_t TB = 1u << 22;
constexpr size_t   ARENA_BYTES = (size_t)1 << 30;
constexpr uint32_t MAXENT_H16 = 1u << 21;
constexpr uint64_t RSV_BYTES = 1u << 20;

struct BHdr {
  NFMemo::Key key;
  uint8_t status, Q, wide;
  uint16_t nabs, nS, nV, nth, pn_len;
  uint64_t h1, h2;
  uint64_t chain;
  uint32_t touch, touch_rsv;
  uint16_t slot_root[MAXTH];
};
static_assert(sizeof(BHdr) % 8 == 0, "block header must keep 8-byte alignment");

static inline uint64_t sm(uint64_t x) {
  x += 0x9e3779b97f4a7c15ull; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}
static inline uint64_t hash_words(const uint32_t* w, size_t n) {
  uint64_t h = 0x243f6a8885a308d3ull ^ n;
  for (size_t i = 0; i < n; i++) h = sm(h ^ w[i]);
  return h;
}
static inline uint64_t chain_prefix(const uint32_t* w, size_t nabs, uint32_t nw) {
  uint64_t h = 0x243f6a8885a308d3ull ^ nw;
  for (size_t i = 0; i < nabs; i++) h = sm(h ^ w[i]);
  return h;
}

struct Shared {
  std::atomic<uint64_t>* tab;
  char* ar;
  std::atomic<uint64_t> bump;
  std::atomic<uint32_t> nent;
  std::atomic<long> drops;
  Shared() : tab(new std::atomic<uint64_t>[TB]()) {
    static char* prev_ar = nullptr; if (prev_ar) munmap(prev_ar, ARENA_BYTES);
    ar = (char*)mmap(nullptr, ARENA_BYTES, PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE|MAP_NORESERVE, -1, 0); if (ar == (char*)MAP_FAILED) { fprintf(stderr, "IOS_BIGALLOC FAIL nf arena\n"); ar = nullptr; }
    prev_ar = ar;
    bump.store(0); nent.store(0); drops.store(0);
  }
  ~Shared() { delete[] tab; }
  Shared(const Shared&) = delete;
  Shared& operator=(const Shared&) = delete;
  void prefault(size_t bytes) { memset(ar, 0, std::min(bytes, ARENA_BYTES)); }
};

static inline void replay_to(const PNode* nd, int pn_len, const uint16_t* slot_root, int nth, const float* c, uint32_t* th) {
  double v[MAXNODES];
  for (int i = 0; i < pn_len; ++i) {
    const PNode& N = nd[i];
    double r;
    switch (N.kind) {
      case TE_LIN: r = 0; break;
      case TE_FIX: r = 1; break;
      case TE_ONEV: r = 1; break;
      case TE_LEAF: r = (double)c[N.leaf]; break;
      case TE_NEG: r = -v[N.a]; break;
      case TE_UN: {
        double x = v[N.a];
        switch (N.op) {
          case OP_SIN: r = std::sin(x); break;
          case OP_COS: r = std::cos(x); break;
          case OP_EXP: r = std::exp(x); break;
          case OP_LOG: r = std::log(std::fabs(x)); break;
          case OP_SQ: r = x * x; break;
          case OP_CUBE: r = x * x * x; break;
          case OP_SQRT: r = std::sqrt(std::fabs(x)); break;
          case OP_CBRT: r = std::cbrt(x); break;
          default: r = x;
        }
        break;
      }
      default: {
        double x = v[N.a], y = v[N.b];
        switch (N.op) {
          case OP_ADD: r = x + y; break;
          case OP_SUB: r = x - y; break;
          case OP_MUL: r = x * y; break;
          default: r = x / y;
        }
      }
    }
    v[i] = r;
  }
  for (int s = 0; s < nth; ++s) {
    float f = (float)v[slot_root[s]];
    memcpy(th + s, &f, 4);
  }
}

struct alignas(128) Handle {
  Shared& S;
  Ctx* C;
  uint64_t r_next = 0, r_end = 0;
  DagBuilder DB;
  Prog P;
  long hits = 0, misses = 0, fails = 0, bigskip = 0;
  Handle(Shared& s, Ctx* c) : S(s), C(c) {}
};

static inline bool find(const Shared& S, const NFMemo::Key& k, uint64_t h, uint64_t& boff) {
  const uint32_t tag = (uint32_t)(h >> 32);
  uint32_t p = (uint32_t)(h & (TB - 1));
  for (;;) {
    const uint64_t v = S.tab[p].load(std::memory_order_acquire);
    if (v == 0) return false;
    if ((uint32_t)(v >> 32) == tag) {
      const uint64_t off = (uint64_t)((uint32_t)v - 1u);
      const BHdr* hd = (const BHdr*)(S.ar + off);
      if (memcmp(&hd->key, &k, sizeof k) == 0) { boff = off; return true; }
    }
    p = (p + 1) & (TB - 1);
  }
}

static inline void publish(Handle& H, uint64_t h, const NFMemo::Key& k, const BHdr& X, const char* body, size_t nbody) {
  Shared& S = H.S;
  if (S.nent.load(std::memory_order_relaxed) >= MAXENT_H16) { S.drops.fetch_add(1); return; }
  const size_t bytes = (sizeof(BHdr) + nbody + 7) & ~(size_t)7;
  if (H.r_next + bytes > H.r_end) { const uint64_t blk = std::max<uint64_t>(RSV_BYTES, bytes); const uint64_t base = S.bump.fetch_add(blk); H.r_next = base; H.r_end = base + blk; }
  const uint64_t off = H.r_next; H.r_next += bytes;
  if (off + bytes > ARENA_BYTES) { S.drops.fetch_add(1); return; }
  memcpy(S.ar + off, &X, sizeof X);
  if (nbody) memcpy(S.ar + off + sizeof X, body, nbody);
  const uint32_t tag = (uint32_t)(h >> 32);
  const uint64_t packed = ((uint64_t)tag << 32) | (off + 1);
  uint32_t p = (uint32_t)(h & (TB - 1));
  for (;;) {
    uint64_t expected = 0;
    if (S.tab[p].compare_exchange_strong(expected, packed, std::memory_order_release, std::memory_order_acquire)) {
      S.nent.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    if ((uint32_t)(expected >> 32) == tag) {
      const BHdr* hd = (const BHdr*)(S.ar + ((uint32_t)expected - 1u));
      if (memcmp(&hd->key, &k, sizeof k) == 0) return;
    }
    p = (p + 1) & (TB - 1);
  }
}

static inline bool miss_compute(Handle& H, const NFMemo::Key& k, uint64_t h, const uint8_t* op, const float* c, int n) {
  H.misses++;
  NFOut out;
  bool good = false;
  uint8_t status = 0;
  if (nf_tree(*H.C, op, c, n, &out)) { good = compile(*H.C, out, c, H.P); status = good ? 2 : 1; }
  BHdr X; memset(&X, 0, sizeof X);
  X.touch = h156::g_touch;
  memcpy(&X.key, &k, sizeof k);
  X.status = status;
  if (good) {
    X.nabs = (uint16_t)H.P.nabs; X.nS = (uint16_t)H.P.nS; X.nV = (uint16_t)H.P.nV; X.Q = (uint8_t)H.P.Q;
    X.nth = (uint16_t)H.P.nth; X.wide = H.P.wide ? 1 : 0; X.h1 = H.P.h1; X.h2 = H.P.h2;
    if (H.DB.build(*H.C, X.slot_root)) {
      X.pn_len = (uint16_t)H.DB.out.size();
      X.chain = chain_prefix(H.P.w, (size_t)H.P.nabs, (uint32_t)(H.P.nabs + H.P.nth));
      static const bool g_arena = getenv("NF_ARENA") && atoi(getenv("NF_ARENA")) == 1;
      if (g_arena) {
        static thread_local std::vector<char> body;
        body.resize((size_t)H.P.nabs * 4 + H.DB.out.size() * sizeof(PNode));
        if (H.P.nabs) memcpy(body.data(), H.P.w, (size_t)H.P.nabs * 4);
        if (!H.DB.out.empty()) memcpy(body.data() + (size_t)H.P.nabs * 4, H.DB.out.data(), H.DB.out.size() * sizeof(PNode));
        publish(H, h, k, X, body.data(), body.size());
      } else {
      std::vector<char> body((size_t)H.P.nabs * 4 + H.DB.out.size() * sizeof(PNode));
      if (H.P.nabs) memcpy(body.data(), H.P.w, (size_t)H.P.nabs * 4);
      if (!H.DB.out.empty()) memcpy(body.data() + (size_t)H.P.nabs * 4, H.DB.out.data(), H.DB.out.size() * sizeof(PNode));
      publish(H, h, k, X, body.data(), body.size());
      }
    } else H.bigskip++;
  } else publish(H, h, k, X, nullptr, 0);
  return good;
}

struct Block { std::vector<uint32_t> w; };
struct Child {
  uint32_t koff = 0, klen = 0;
  uint32_t roff[2] = {0, 0}, rlen[2] = {0, 0};
  uint8_t ok[2] = {0, 0};
  uint64_t r0 = 0, r1 = 0, kh = 0;
};

static const uint64_t KEYMASK[25][3] = {
#define KM(n) {n >= 8 ? ~0ull : (n <= 0 ? 0ull : (~0ull >> (8 * (8 - n)))), \
               n >= 16 ? ~0ull : (n <= 8 ? 0ull : (~0ull >> (8 * (16 - n)))), \
               n >= 24 ? ~0ull : (n <= 16 ? 0ull : (~0ull >> (8 * (24 - n))))}
    KM(0), KM(1), KM(2), KM(3), KM(4), KM(5), KM(6), KM(7), KM(8), KM(9), KM(10), KM(11), KM(12),
    KM(13), KM(14), KM(15), KM(16), KM(17), KM(18), KM(19), KM(20), KM(21), KM(22), KM(23), KM(24)
#undef KM
};
static inline void make_key(const uint8_t* op, int n, NFMemo::Key& k) {
  uint64_t w[3];
  memcpy(w, op, 24);
  const int nn = (n >= 0 && n <= MAXN) ? n : 0;
  for (int j = 0; j < 3; j++) w[j] &= KEYMASK[nn][j];
  k.len = (uint8_t)n;
  memcpy(k.op, w, 24);
}
static inline bool key_eq(const NFMemo::Key& a, const NFMemo::Key& b) {
  uint64_t x[4] = {0, 0, 0, 0}, y[4] = {0, 0, 0, 0};
  memcpy(x, &a, sizeof a); memcpy(y, &b, sizeof b);
  return ((x[0] ^ y[0]) | (x[1] ^ y[1]) | (x[2] ^ y[2]) | ((x[3] ^ y[3]) & 0xFFu)) == 0;
}

static inline bool find2(const Shared& S, const NFMemo::Key& k, uint64_t h, uint64_t& boff) {
  const uint32_t tag = (uint32_t)(h >> 32);
  uint32_t p = (uint32_t)(h & (TB - 1));
  for (;;) {
    const uint64_t v = S.tab[p].load(std::memory_order_acquire);
    if (v == 0) return false;
    if ((uint32_t)(v >> 32) == tag) {
      const uint64_t off = (uint64_t)((uint32_t)v - 1u);
      if (key_eq(((const BHdr*)(S.ar + off))->key, k)) { boff = off; return true; }
    }
    p = (p + 1) & (TB - 1);
  }
}

static inline bool tree_seg(Handle& H, const skn::Cand& cd, int t, const NFMemo::Key& k, uint64_t hk, Block& B, Child& ch, uint64_t& r) {
  const size_t t0 = B.w.size();
  const int n = cd.len[t];
  uint64_t boff = 0;
  bool ok = false;
  if (find2(H.S, k, hk, boff)) {
    H.hits++;
    const BHdr* hd = (const BHdr*)(H.S.ar + boff);
    if (h156::g_on) h156::touch_u32(const_cast<uint32_t*>(&hd->touch), h156::g_touch);
    if (hd->status == 2 && hd->Q <= 16) {
      const uint32_t nabs = hd->nabs, nth = hd->nth, nw = nabs + nth;
      B.w.push_back(nw);
      const size_t r0 = B.w.size();
      B.w.resize(r0 + nw);
      const uint32_t* aw = (const uint32_t*)(H.S.ar + boff + sizeof(BHdr));
      if (nabs) memcpy(B.w.data() + r0, aw, sizeof(uint32_t) * nabs);
      const PNode* nd = (const PNode*)(H.S.ar + boff + sizeof(BHdr) + sizeof(uint32_t) * nabs);
      replay_to(nd, hd->pn_len, hd->slot_root, (int)nth, cd.c[t], B.w.data() + r0 + nabs);
      uint64_t h = hd->chain;
      for (uint32_t s = 0; s < nth; s++) h = sm(h ^ B.w[r0 + nabs + s]);
      ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw; r = h; ok = true;
    }
  } else {
    const bool good = miss_compute(H, k, hk, cd.op[t], cd.c[t], n);
    if (good) {
      const uint32_t nw = (uint32_t)H.P.nw;
      B.w.push_back(nw);
      const size_t r0 = B.w.size();
      B.w.resize(r0 + nw);
      if (nw) memcpy(B.w.data() + r0, H.P.w, sizeof(uint32_t) * nw);
      ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw;
      r = hash_words(B.w.data() + r0, nw);
      ok = H.P.Q <= 16;
    }
  }
  if (!ok) {
    H.fails++;
    B.w.resize(t0);
    const size_t s0 = B.w.size();
    B.w.push_back(0xFFFFFFFFu); B.w.push_back((uint32_t)cd.len[t]);
    for (int j = 0; j < cd.len[t]; j++) { uint32_t b; std::memcpy(&b, &cd.c[t][j], 4); B.w.push_back(cd.op[t][j]); B.w.push_back(b); }
    r = hash_words(B.w.data() + s0, B.w.size() - s0);
    ch.rlen[t] = 0;
  }
  return ok;
}

static void chunk(Handle& H, const skn::Cand* cands, const uint8_t* valid, int b0, int b1, Block& B, Child* out) {
  B.w.clear();
  const int m = b1 - b0;
  NFMemo::Key K[64][2];
  uint64_t HK[64][2];
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      make_key(cands[i].op[t], cands[i].len[t], K[j][t]);
      HK[j][t] = NFMemo::hkey(K[j][t]);
      __builtin_prefetch(&H.S.tab[(uint32_t)(HK[j][t] & (TB - 1))], 0, 3);
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    if (!valid[i]) continue;
    for (int t = 0; t < 2; t++) {
      const uint64_t v = H.S.tab[(uint32_t)(HK[j][t] & (TB - 1))].load(std::memory_order_relaxed);
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
      const uint64_t v = H.S.tab[(uint32_t)(HK[j][t] & (TB - 1))].load(std::memory_order_acquire);
      if (!v || (uint32_t)(v >> 32) != (uint32_t)(HK[j][t] >> 32)) continue;
      const char* bp = H.S.ar + ((uint32_t)v - 1u);
      const BHdr* hd = (const BHdr*)bp;
      const size_t nb = (size_t)hd->nabs * 4 + (size_t)hd->pn_len * sizeof(PNode);
      const size_t lim = std::min<size_t>(nb, 768);
      for (size_t o = 192; o < sizeof(BHdr) + lim; o += 64) __builtin_prefetch(bp + o, 0, 3);
    }
  }
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
    Child& ch = out[i];
    ch = Child();
    if (!valid[i]) { ch.koff = (uint32_t)B.w.size(); continue; }
    const skn::Cand& cd = cands[i];
    ch.koff = (uint32_t)B.w.size();
    uint64_t r[2] = {0, 0};
    ch.ok[0] = tree_seg(H, cd, 0, K[j][0], HK[j][0], B, ch, r[0]) ? 1 : 0;
    ch.ok[1] = tree_seg(H, cd, 1, K[j][1], HK[j][1], B, ch, r[1]) ? 1 : 0;
    ch.klen = (uint32_t)(B.w.size() - ch.koff);
    ch.r0 = r[0]; ch.r1 = r[1];
    ch.kh = sm(r[0] ^ sm(r[1] + ch.klen));
  }
}

}
