#include <chrono>
#pragma once
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <functional>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>
#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#endif
#ifndef RK_BITS
#define RK_BITS 8
#endif

namespace evo {

constexpr int MAXLEN = 24;
constexpr int NCASE = 96;
#ifndef IND_NRES
#define IND_NRES 1
#endif
constexpr int NRES = IND_NRES;
#ifndef LOOP_DECK
#define LOOP_DECK 512
#endif
constexpr int DECK = LOOP_DECK;
constexpr int NDEME = 16;
constexpr int NPOP = DECK * NDEME;
constexpr int TOURN = 7;
constexpr int MIGR = DECK / 100;
constexpr int LEX_CASES = 64;
inline int mem_lean() { static const int v = getenv("MEM_LEAN") ? atoi(getenv("MEM_LEAN")) : 7; return v; }
static_assert(2 * MAXLEN < 64, "cval_of masks complexity to 6 bits; F[64] is indexed by it");
constexpr int MAX_TRY = 8;
constexpr double P_XO = 0.9, P_SWAP = 0.5, P_MUT = 0.25;

enum : uint8_t { NOP = 0, CONST = 1, TVAR = 2, ADD = 3, SUB = 4, MUL = 5, DIV = 6,
                 SIN = 7, COS = 8, EXP = 9, LOG = 10, SQ = 11, CUBE = 12, SQRT = 13, CBRT = 14 };

inline int arity(uint8_t o) { return (o >= ADD && o <= DIV) ? 2 : ((o >= SIN && o <= CBRT) ? 1 : 0); }

struct Cand {
  uint8_t len[2];
  uint8_t op[2][MAXLEN];
  float c[2][MAXLEN];
};
static_assert(sizeof(Cand) == 244, "Cand must be 244 B");

inline int cplx(const Cand& g) { return g.len[0] + g.len[1]; }

struct Ind {
  Cand g;
  uint8_t ext[2][MAXLEN];
  float res[NRES];
  float err;
  double crowd;
  uint16_t rank;
  uint16_t cx;
  uint8_t valid;
};

inline uint64_t splitmix(uint64_t x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}
inline uint64_t mix(uint64_t a, uint64_t b) { return splitmix(a * 0x100000001B3ull ^ splitmix(b)); }

struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(splitmix(seed)) { if (!s) s = 0x9E3779B97F4A7C15ull; }
  uint64_t next() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 0x2545F4914F6CDD1Dull; }
  uint32_t below(uint32_t n) { return (uint32_t)(((next() >> 32) * (uint64_t)n) >> 32); }
  double uni() { return (double)(next() >> 11) * 0x1.0p-53; }
  double gauss() {
    double u1 = uni(); if (u1 < 1e-300) u1 = 1e-300;
    double u2 = uni();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
  }
};

inline void extents(const uint8_t* op, int L, uint8_t* ext) {
  uint8_t st[MAXLEN];
  int sp = 0;
  for (int k = 0; k < L; k++) {
    int a = arity(op[k]);
    if (a == 0) { ext[k] = (uint8_t)k; st[sp++] = (uint8_t)k; }
    else if (a == 1) { ext[k] = st[sp - 1]; }
    else { sp--; ext[k] = st[sp - 1]; }
  }
}

inline bool tree_valid(const uint8_t* op, int L) {
  if (L < 1 || L > MAXLEN) return false;
  int depth = 0;
  for (int k = 0; k < L; k++) {
    uint8_t o = op[k];
    if (o == NOP || o > CBRT) return false;
    int a = arity(o);
    if (a == 0) depth++;
    else if (a == 1) { if (depth < 1) return false; }
    else { if (depth < 2) return false; depth--; }
  }
  return depth == 1;
}
inline bool cand_valid(const Cand& g) { return tree_valid(g.op[0], g.len[0]) && tree_valid(g.op[1], g.len[1]); }

inline void mutate(Cand& g, int t, Rng& r, int* kout = nullptr) {
  int k = (int)r.below(g.len[t]);
  if (kout) *kout = k;
  uint8_t o = g.op[t][k];
  if (o == CONST) {
    float c = g.c[t][k];
    g.c[t][k] = c + (float)(r.gauss() * (0.1 * std::fabs(c) + 0.1));
  } else if (o == TVAR) {
    g.op[t][k] = CONST; g.c[t][k] = (float)r.gauss();
  } else if (arity(o) == 1) {
    static const uint8_t U[8] = {SIN, COS, EXP, LOG, SQ, CUBE, SQRT, CBRT};
    uint8_t cand[7]; int m = 0;
    for (int q = 0; q < 8; q++) if (U[q] != o) cand[m++] = U[q];
    g.op[t][k] = cand[r.below(7)];
  } else {
    static const uint8_t B[4] = {ADD, SUB, MUL, DIV};
    uint8_t cand[3]; int m = 0;
    for (int q = 0; q < 4; q++) if (B[q] != o) cand[m++] = B[q];
    g.op[t][k] = cand[r.below(3)];
  }
}

inline int tournament(const Ind* par, int base, Rng& r) {
  int best = base + (int)r.below(DECK);
  for (int k = 1; k < TOURN; k++) {
    int c = base + (int)r.below(DECK);
    if (par[c].rank < par[best].rank || (par[c].rank == par[best].rank && par[c].crowd > par[best].crowd)) best = c;
  }
  return best;
}

inline bool make_child(const Ind* par, int base, Rng& r, Ind& out, uint32_t& rebreeds) {
  for (int attempt = 0; attempt < MAX_TRY; attempt++) {
    const Ind& A = par[tournament(par, base, r)];
    const Ind& B = par[tournament(par, base, r)];
    bool swap = r.uni() < P_SWAP;
    Cand g = A.g;
    bool ok = true;
    for (int t = 0; t < 2 && ok; t++) {
      if (r.uni() < P_XO) {
        int tb = swap ? 1 - t : t;
        int LA = A.g.len[t], LB = B.g.len[tb];
        int i = (int)r.below(LA), j = (int)r.below(LB);
        int s = A.ext[t][i], szi = i - s + 1;
        int js = B.ext[tb][j], szj = j - js + 1;
        int nl = LA - szi + szj;
        if (nl > MAXLEN) { ok = false; rebreeds++; break; }
        const uint8_t* aop = A.g.op[t]; const float* ac = A.g.c[t];
        const uint8_t* bop = B.g.op[tb]; const float* bc = B.g.c[tb];
        uint8_t* dop = g.op[t]; float* dc = g.c[t];
        int tail = LA - i - 1;
        std::memcpy(dop, aop, s);             std::memcpy(dc, ac, s * sizeof(float));
        std::memcpy(dop + s, bop + js, szj);  std::memcpy(dc + s, bc + js, szj * sizeof(float));
        std::memcpy(dop + s + szj, aop + i + 1, tail);
        std::memcpy(dc + s + szj, ac + i + 1, tail * sizeof(float));
        g.len[t] = (uint8_t)nl;
      }
    }
    if (!ok) continue;
    for (int t = 0; t < 2; t++) if (r.uni() < P_MUT) mutate(g, t, r);
    out.g = g;
    out.valid = 1;
    return true;
  }
  out.valid = 0;
  return false;
}

inline bool (*g_pool_trace_on)() = nullptr; inline double (*g_pool_trace_now)() = nullptr; inline void (*g_pool_trace_rec)(double, double) = nullptr;
inline void (*g_pool_trace_wait)(double, double) = nullptr;
inline int g_pool_h403 = 0;
inline int g_pool_h744 = 0; inline double g_pool_h744s = 20.0;
class Pool {
 public:
  explicit Pool(int nthreads) : nth_(std::max(1, nthreads)) {
    for (int t = 1; t < nth_; t++) workers_.emplace_back([this] { loop(); });
  }
  ~Pool() {
    { std::lock_guard<std::mutex> g(m_); stop_ = true; gen_++; gen_a_.store(gen_, std::memory_order_release); }
    cv_.notify_all();
    for (auto& th : workers_) th.join();
  }
  void run(int ntasks, const std::function<void(int)>& fn) {
    if (nth_ == 1) { for (int i = 0; i < ntasks; i++) fn(i); return; }
    if (g_pool_h403 && g_pool_h744) { run744(ntasks, fn); return; }
    if (g_pool_h403) { run403(ntasks, fn); return; }
    { std::lock_guard<std::mutex> g(m_); fn_ = &fn; ntasks_ = ntasks; next_ = 0; busy_ = nth_ - 1; gen_++; }
    cv_.notify_all();
    drain();
    const bool trw = g_pool_trace_on && g_pool_trace_on(); const double w0 = trw ? g_pool_trace_now() : 0;
    std::unique_lock<std::mutex> lk(m_);
    done_.wait(lk, [&] { return busy_ == 0; });
    if (trw) { lk.unlock(); g_pool_trace_wait(w0, g_pool_trace_now()); }
  }
 private:
  void run744(int ntasks, const std::function<void(int)>& fn) {
    for (;;) { { std::lock_guard<std::mutex> g(m_); if (inflight_a_.load(std::memory_order_acquire) == 0) {
        fn_ = &fn; ntasks_ = ntasks; next_ = 0; left_a_.store(ntasks, std::memory_order_relaxed); gen_++; gen_a_.store(gen_, std::memory_order_release); break; } }
      __builtin_arm_yield(); }
    if (parked_.load(std::memory_order_acquire) > 0) cv_.notify_all();
    const bool tr = g_pool_trace_on && g_pool_trace_on(); const double a = tr ? g_pool_trace_now() : 0;
    int i, k = 0; while ((i = next_.fetch_add(1)) < ntasks) { fn(i); k++; }
    if (tr && k) g_pool_trace_rec(a, g_pool_trace_now());
    if (left_a_.fetch_sub(k, std::memory_order_acq_rel) - k == 0) return;
    const double w0 = tr ? g_pool_trace_now() : 0;
    while (left_a_.load(std::memory_order_acquire) != 0) __builtin_arm_yield();
    if (tr) g_pool_trace_wait(w0, g_pool_trace_now());
  }
  void run403(int ntasks, const std::function<void(int)>& fn) {
    { std::unique_lock<std::mutex> lk(m_); done_.wait(lk, [&] { return inflight_ == 0; });
      fn_ = &fn; ntasks_ = ntasks; next_ = 0; left_ = ntasks; gen_++; }
    cv_.notify_all();
    const bool tr = g_pool_trace_on && g_pool_trace_on(); const double a = tr ? g_pool_trace_now() : 0;
    int i, k = 0; while ((i = next_.fetch_add(1)) < ntasks) { fn(i); k++; }
    if (tr && k) g_pool_trace_rec(a, g_pool_trace_now());
    const double w0 = tr ? g_pool_trace_now() : 0;
    std::unique_lock<std::mutex> lk(m_); left_ -= k;
    done_.wait(lk, [&] { return left_ == 0; });
    if (tr) { lk.unlock(); g_pool_trace_wait(w0, g_pool_trace_now()); }
  }
  void drain() { int i; const bool tr = g_pool_trace_on && g_pool_trace_on(); double a = tr ? g_pool_trace_now() : 0; bool any = false;
    while ((i = next_.fetch_add(1)) < ntasks_) { any = true; (*fn_)(i); }
    if (tr && any) g_pool_trace_rec(a, g_pool_trace_now()); }
  void loop() {
    uint64_t seen = 0;
    for (;;) {
      const std::function<void(int)>* f = nullptr; int n = 0;
      {
        if (g_pool_h744) { const auto t0 = std::chrono::steady_clock::now(); int sp = 0;
          while (gen_a_.load(std::memory_order_acquire) == seen) { __builtin_arm_yield();
            if ((++sp & 255) == 0 && std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() > g_pool_h744s) break; } }
        std::unique_lock<std::mutex> lk(m_);
        if (g_pool_h744) { if (gen_ == seen) { parked_++; cv_.wait(lk, [&] { return gen_ != seen; }); parked_--; } }
        else cv_.wait(lk, [&] { return gen_ != seen; });
        seen = gen_;
        if (stop_) return;
        if (g_pool_h403 && g_pool_h744) { inflight_a_.fetch_add(1, std::memory_order_acq_rel); f = fn_; n = ntasks_; }
        else if (g_pool_h403) { inflight_++; f = fn_; n = ntasks_; }
      }
      if (g_pool_h403 && g_pool_h744) {
        const bool tr = g_pool_trace_on && g_pool_trace_on(); const double a = tr ? g_pool_trace_now() : 0;
        int i, k = 0; while ((i = next_.fetch_add(1)) < n) { (*f)(i); k++; }
        if (tr && k) g_pool_trace_rec(a, g_pool_trace_now());
        left_a_.fetch_sub(k, std::memory_order_acq_rel); inflight_a_.fetch_sub(1, std::memory_order_acq_rel);
        continue;
      }
      if (g_pool_h403) {
        const bool tr = g_pool_trace_on && g_pool_trace_on(); const double a = tr ? g_pool_trace_now() : 0;
        int i, k = 0; while ((i = next_.fetch_add(1)) < n) { (*f)(i); k++; }
        if (tr && k) g_pool_trace_rec(a, g_pool_trace_now());
        std::lock_guard<std::mutex> g(m_); inflight_--; left_ -= k;
        if (left_ == 0 || inflight_ == 0) done_.notify_all();
        continue;
      }
      drain();
      std::lock_guard<std::mutex> g(m_);
      if (--busy_ == 0) done_.notify_one();
    }
  }
  int nth_;
  std::vector<std::thread> workers_;
  std::mutex m_;
  std::condition_variable cv_, done_;
  const std::function<void(int)>* fn_ = nullptr;
  int ntasks_ = 0;
  std::atomic<int> next_{0};
  int busy_ = 0;
  int inflight_ = 0, left_ = 0;
  uint64_t gen_ = 0;
  alignas(128) std::atomic<uint64_t> gen_a_{0}; alignas(128) std::atomic<int> parked_{0};
  alignas(128) std::atomic<int> left_a_{0}; alignas(128) std::atomic<int> inflight_a_{0}; char pad744_[128];
  bool stop_ = false;
};

inline void precompute_ext(Pool& pool, Ind* par) {
  pool.run(NDEME, [&](int d) {
    for (int i = d * DECK; i < (d + 1) * DECK; i++)
      for (int t = 0; t < 2; t++) extents(par[i].g.op[t], par[i].g.len[t], par[i].ext[t]);
  });
}

inline void breed(Pool& pool, const Ind* par, Ind* chi, uint64_t seed, uint64_t step, uint32_t* rebreeds, uint32_t* failed) {
  pool.run(NDEME, [&](int d) {
    Rng r(mix(mix(seed, step), (uint64_t)d));
    uint32_t rb = 0, fl = 0;
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      if (!make_child(par, d * DECK, r, chi[i], rb)) fl++;
    }
    rebreeds[d] = rb; failed[d] = fl;
  });
}

inline int subtree_start(const uint8_t* op, int k) {
  static constexpr int8_t DA[16] = {0, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int need = 1, q = k;
  for (;;) { need += DA[op[q] & 15]; if (need == 0) return q; q--; }
}

struct BreedBuf {
  std::vector<Cand> G;
  std::vector<uint16_t> prk;
  std::vector<double> pcw;
  BreedBuf() : G(((mem_lean() & 2) && getenv("H777") && atoi(getenv("H777")) != 0) ? 0 : NPOP), prk(NPOP), pcw(NPOP) {}
};

inline int tournament_i12(const BreedBuf& BB, int base, Rng& r) {
  const uint16_t* prk = BB.prk.data();
  const double* pcw = BB.pcw.data();
  uint32_t best = (uint32_t)(base + (int)r.below(DECK));
  for (int k = 1; k < TOURN; k++) {
    const uint32_t c = (uint32_t)(base + (int)r.below(DECK));
    const uint32_t better = (uint32_t)(prk[c] < prk[best]) | ((uint32_t)(prk[c] == prk[best]) & (uint32_t)(pcw[c] > pcw[best]));
    const uint32_t m = 0u - better;
    best = (c & m) | (best & ~m);
  }
  return (int)best;
}

__attribute__((always_inline)) inline void splice_tree(Cand& g, int t, const Cand& A, const Cand& B, int tb,
                                                       int s, int i, int js, int szj, int nl) {
  uint8_t U[2 * MAXLEN];
  float Uc[2 * MAXLEN];
  std::memcpy(U, A.op[t], MAXLEN);            std::memcpy(U + MAXLEN, B.op[tb], MAXLEN);
  std::memcpy(Uc, A.c[t], MAXLEN * sizeof(float)); std::memcpy(Uc + MAXLEN, B.c[tb], MAXLEN * sizeof(float));
  const int se = s + szj;
  uint8_t idx[MAXLEN];
#if defined(__ARM_NEON) || defined(__aarch64__)
  const uint8x16_t kA = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  const uint8x16_t kB = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
  const uint8x16_t vs = vdupq_n_u8((uint8_t)s), vse = vdupq_n_u8((uint8_t)se), vnl = vdupq_n_u8((uint8_t)nl);
  const uint8x16_t vb = vdupq_n_u8((uint8_t)(MAXLEN + js - s)), vc = vdupq_n_u8((uint8_t)(i + 1 - se));
  auto map = [&](uint8x16_t k) {
    const uint8x16_t b3 = vaddq_u8(k, vc), b2 = vaddq_u8(k, vb);
    const uint8x16_t v3 = vbslq_u8(vcltq_u8(k, vnl), b3, k);
    const uint8x16_t v2 = vbslq_u8(vcltq_u8(k, vse), b2, v3);
    return vbslq_u8(vcltq_u8(k, vs), k, v2);
  };
  const uint8x16_t ia = map(kA), ib = map(kB);
  const uint8x16x3_t T = {{vld1q_u8(U), vld1q_u8(U + 16), vld1q_u8(U + 32)}};
  vst1q_u8(g.op[t], vqtbl3q_u8(T, ia));
  vst1_u8(g.op[t] + 16, vget_low_u8(vqtbl3q_u8(T, ib)));
  vst1q_u8(idx, ia); vst1_u8(idx + 16, vget_low_u8(ib));
  for (int k = 0; k < MAXLEN; k++) g.c[t][k] = Uc[idx[k]];
#else
  for (int k = 0; k < MAXLEN; k++) {
    const uint32_t a = (uint32_t)k, b = (uint32_t)(MAXLEN + js + k - s), c3 = (uint32_t)(i + 1 + k - se);
    const uint32_t m3 = 0u - (uint32_t)(k < nl), m2 = 0u - (uint32_t)(k < se), m1 = 0u - (uint32_t)(k < s);
    const uint32_t v3 = (c3 & m3) | (a & ~m3);
    const uint32_t v2 = (b & m2) | (v3 & ~m2);
    const uint32_t id = (a & m1) | (v2 & ~m1);
    idx[k] = (uint8_t)id;
    g.op[t][k] = U[id];
    g.c[t][k] = Uc[id];
  }
#endif
  g.len[t] = (uint8_t)nl;
}
struct ChildLink320 { int32_t ia, ib; uint8_t swap, xo[2], pad; uint8_t s[2], szi[2], js[2], szj[2], k[2], pad2[2]; };
static_assert(sizeof(ChildLink320) == 24, "ChildLink320 is 24 bytes");
inline ChildLink320* g_link320 = nullptr;
__attribute__((always_inline)) inline bool make_child_i12(const BreedBuf& BB, int base, Rng& r, Ind& out, uint32_t& rebreeds, ChildLink320* pv = nullptr) {
  for (int attempt = 0; attempt < MAX_TRY; attempt++) {
    const int ia = tournament_i12(BB, base, r);
    const int ib = tournament_i12(BB, base, r);
    const Cand& A = BB.G[ia];
    const Cand& B = BB.G[ib];
    const bool swap = r.uni() < P_SWAP;
    if (pv) { std::memset(pv, 0, sizeof *pv); pv->ia = ia; pv->ib = ib; pv->swap = swap ? 1 : 0; pv->k[0] = 255; pv->k[1] = 255; }
    Cand g;
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

inline void breed_i12(Pool& pool, const Ind* par, Ind* chi, BreedBuf& BB, uint64_t seed, uint64_t step, uint32_t* rebreeds, uint32_t* failed) {
  pool.run(NDEME, [&](int d) {
    for (int i = d * DECK; i < (d + 1) * DECK; i++) { BB.G[i] = par[i].g; BB.prk[i] = par[i].rank; BB.pcw[i] = par[i].crowd; }
    Rng r(mix(mix(seed, step), (uint64_t)d));
    uint32_t rb = 0, fl = 0;
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      if (!make_child_i12(BB, d * DECK, r, chi[i], rb)) fl++;
    }
    rebreeds[d] = rb; failed[d] = fl;
  });
}

struct RankScratch {
  std::vector<uint64_t> A, B, FE, tmp;
  std::vector<uint16_t> rS;
  std::vector<uint32_t> off, pos;
};

inline void f_raise(uint16_t* F, int c, uint16_t r) {
#if defined(__ARM_NEON) || defined(__aarch64__)
  static const uint16_t IDX[64] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                   16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                                   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                                   48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
  const uint16x8_t rv = vdupq_n_u16(r), cv = vdupq_n_u16((uint16_t)c), zv = vdupq_n_u16(0);
  for (int v = 0; v < 8; v++) {
    const uint16x8_t ge = vcgeq_u16(vld1q_u16(IDX + 8 * v), cv);
    const uint16x8_t val = vbslq_u16(ge, rv, zv);
    vst1q_u16(F + 8 * v, vmaxq_u16(vld1q_u16(F + 8 * v), val));
  }
#else
  for (int i = c; i < 64 && F[i] < r; i++) F[i] = r;
#endif
}

inline float fix_err(float e) {
  if (std::isnan(e) || e > std::numeric_limits<float>::max()) return std::numeric_limits<float>::max();
  return e < -std::numeric_limits<float>::max() ? -std::numeric_limits<float>::max() : e;
}
inline uint32_t fkey(float f) { uint32_t b; std::memcpy(&b, &f, 4); return (b & 0x80000000u) ? ~b : (b | 0x80000000u); }
inline float funkey(uint32_t k) { const uint32_t b = (k & 0x80000000u) ? (k & 0x7fffffffu) : ~k; float f; std::memcpy(&f, &b, 4); return f; }

inline double dval_of(uint64_t x) { return (double)funkey((uint32_t)(x >> 32)); }
inline int cval_of(uint64_t x) { return (int)((x >> 16) & 63); }
inline uint32_t jof(uint64_t x) { return (uint32_t)(x & 0xffffu); }

inline const double* g_pen_tab = nullptr;
inline void rank_crowd_pen(int m, const float* err, const uint16_t* cp, uint16_t* rank, double* crowd, RankScratch& sc);
inline void rank_crowd(int m, const float* err, const uint16_t* cp, uint16_t* rank, double* crowd, RankScratch& sc) {
  if (g_pen_tab) { rank_crowd_pen(m, err, cp, rank, crowd, sc); return; }
  sc.A.resize(m); sc.B.resize(m); sc.FE.resize(m); sc.tmp.resize(m); sc.rS.resize(m);
  constexpr int RB = RK_BITS, RNB = 1 << RB, RNP = (32 + RB - 1) / RB;
  static thread_local uint16_t h[RNP][RNB];
  std::memset(h, 0, sizeof h);
  uint64_t* src = sc.A.data();
  uint64_t* dst = sc.B.data();
  for (int j = 0; j < m; j++) {
    float e = fix_err(err[j]);
    uint32_t b = fkey(e);
    uint32_t c = cp[j];
    src[j] = ((uint64_t)b << 32) | ((uint64_t)c << 16) | (uint64_t)j;
    for (int p = 0; p < RNP; p++) h[p][(b >> (RB * p)) & (RNB - 1)]++;
  }
  for (int p = 0; p < RNP; p++) {
    const int sh = 32 + RB * p;
    uint16_t* hp = h[p];
    uint32_t sum = 0;
    bool trivial = false;
    for (int d = 0; d < RNB; d++) { uint32_t t = hp[d]; if (t == (uint32_t)m) trivial = true; hp[d] = (uint16_t)sum; sum += t; }
    if (trivial) continue;
    for (int k = 0; k < m; k++) { uint64_t x = src[k]; dst[hp[(x >> sh) & (RNB - 1)]++] = x; }
    std::swap(src, dst);
  }
  uint16_t F[64];
  std::memset(F, 0, sizeof F);
  uint16_t* rs = sc.rS.data();
  int R = 0;
  for (int k = 0; k < m;) {
    const uint32_t eb = (uint32_t)(src[k] >> 32);
    int e = k + 1;
    while (e < m && (uint32_t)(src[e] >> 32) == eb) e++;
    if (e - k > 1) std::sort(src + k, src + e, [](uint64_t x, uint64_t y) { return (uint32_t)x < (uint32_t)y; });
    for (int g = k; g < e;) {
      const uint64_t key = src[g] >> 16;
      int ge = g + 1;
      while (ge < e && (src[ge] >> 16) == key) ge++;
      const int c = (int)(key & 63);
      uint16_t r = (uint16_t)(F[c] + 1);
      for (int q = g; q < ge; q++) { rs[q] = r; rank[jof(src[q])] = r; }
      if ((int)r > R) R = r;
      f_raise(F, c, r);
      g = ge;
    }
    k = e;
  }
  sc.off.assign(R + 2, 0);
  for (int k = 0; k < m; k++) sc.off[rs[k] + 1]++;
  for (int r = 1; r <= R + 1; r++) sc.off[r] += sc.off[r - 1];
  sc.pos = sc.off;
  for (int k = 0; k < m; k++) sc.FE[sc.pos[rs[k]]++] = src[k];

  const double INF = std::numeric_limits<double>::infinity();
  for (int r = 1; r <= R; r++) {
    const uint32_t a = sc.off[r], b = sc.off[r + 1];
    const int n = (int)(b - a);
    if (n == 0) continue;
    const uint64_t* FE = sc.FE.data() + a;
    uint64_t* T = sc.tmp.data();
    const uint32_t jE0 = jof(FE[0]), jE1 = jof(FE[n - 1]);
    const double rge = dval_of(FE[n - 1]) - dval_of(FE[0]);
    crowd[jE0] = INF; crowd[jE1] = INF;
    for (int q = 1; q < n - 1; q++)
      crowd[jof(FE[q])] = rge > 0 ? (dval_of(FE[q + 1]) - dval_of(FE[q - 1])) / rge : 0.0;
    for (int i = 0; i < n; i++) T[i] = FE[n - 1 - i];
    for (int i = 0; i < n;) {
      uint64_t key = T[i] >> 16;
      int e = i + 1;
      while (e < n && (T[e] >> 16) == key) e++;
      if (e - i > 1) std::reverse(T + i, T + e);
      i = e;
    }
    const double rgc = (double)cval_of(T[n - 1]) - (double)cval_of(T[0]);
    for (int p = 0; p < n; p++) {
      const uint32_t j = jof(T[p]);
      const bool endC = (p == 0 || p == n - 1);
      if (endC || j == jE0 || j == jE1) crowd[j] = INF;
      else if (rgc > 0) crowd[j] += ((double)cval_of(T[p + 1]) - (double)cval_of(T[p - 1])) / rgc;
    }
  }
}

inline void rank_crowd_pen(int m, const float* err, const uint16_t* cp, uint16_t* rank, double* crowd, RankScratch& sc) {
  if (m > 65535) { fprintf(stderr, "rank_crowd_pen: pool of %d members; the 16-bit index/histogram layout (as rank_crowd) holds at most 65535\n", m); std::abort(); }
  sc.A.resize(m); sc.B.resize(m); sc.FE.resize(m); sc.tmp.resize(m); sc.rS.resize(m);
  constexpr int RB = RK_BITS, RNB = 1 << RB, RNP = (32 + RB - 1) / RB;
  static thread_local uint16_t h[RNP][RNB];
  std::memset(h, 0, sizeof h);
  uint64_t* src = sc.A.data();
  uint64_t* dst = sc.B.data();
  uint32_t maxc = 0;
  for (int j = 0; j < m; j++) {
    float e = fix_err(err[j]);
    uint32_t b = fkey(e);
    uint32_t c = cp[j]; if (c > maxc) maxc = c;
    src[j] = ((uint64_t)b << 32) | ((uint64_t)c << 16) | (uint64_t)j;
    for (int p = 0; p < RNP; p++) h[p][(b >> (RB * p)) & (RNB - 1)]++;
  }
  for (int p = 0; p < RNP; p++) {
    const int sh = 32 + RB * p;
    uint16_t* hp = h[p];
    uint32_t sum = 0;
    bool trivial = false;
    for (int d = 0; d < RNB; d++) { uint32_t t = hp[d]; if (t == (uint32_t)m) trivial = true; hp[d] = (uint16_t)sum; sum += t; }
    if (trivial) continue;
    for (int k = 0; k < m; k++) { uint64_t x = src[k]; dst[hp[(x >> sh) & (RNB - 1)]++] = x; }
    std::swap(src, dst);
  }
  static thread_local std::vector<uint16_t> fw; static thread_local std::vector<uint32_t> fs; static thread_local uint32_t fcur = 0;
  const uint32_t L = maxc + 1;
  if (fw.size() < L + 1) { fw.assign(L + 1, 0); fs.assign(L + 1, 0); fcur = 0; }
  if (++fcur == 0) { std::fill(fs.begin(), fs.end(), 0u); fcur = 1; }
  auto fq = [&](uint32_t c) { uint16_t r = 0; for (uint32_t i = c + 1; i > 0; i -= i & (0u - i)) if (fs[i] == fcur && fw[i] > r) r = fw[i]; return r; };
  auto fu = [&](uint32_t c, uint16_t r) { for (uint32_t i = c + 1; i <= L; i += i & (0u - i)) { if (fs[i] != fcur) { fs[i] = fcur; fw[i] = 0; } if (fw[i] < r) fw[i] = r; } };
  uint16_t* rs = sc.rS.data();
  int R = 0;
  for (int k = 0; k < m;) {
    const uint32_t eb = (uint32_t)(src[k] >> 32);
    int e = k + 1;
    while (e < m && (uint32_t)(src[e] >> 32) == eb) e++;
    if (e - k > 1) std::sort(src + k, src + e, [](uint64_t x, uint64_t y) { return (uint32_t)x < (uint32_t)y; });
    for (int g = k; g < e;) {
      const uint64_t key = src[g] >> 16;
      int ge = g + 1;
      while (ge < e && (src[ge] >> 16) == key) ge++;
      const uint32_t c = (uint32_t)(key & 0xffffu);
      const uint16_t r = (uint16_t)(fq(c) + 1);
      for (int q = g; q < ge; q++) { rs[q] = r; rank[jof(src[q])] = r; }
      if ((int)r > R) R = r;
      fu(c, r);
      g = ge;
    }
    k = e;
  }
  sc.off.assign(R + 2, 0);
  for (int k = 0; k < m; k++) sc.off[rs[k] + 1]++;
  for (int r = 1; r <= R + 1; r++) sc.off[r] += sc.off[r - 1];
  sc.pos = sc.off;
  for (int k = 0; k < m; k++) sc.FE[sc.pos[rs[k]]++] = src[k];
  const double* PT = g_pen_tab;
  auto pv = [&](uint64_t x) { return PT[(x >> 16) & 0xffffu]; };
  const double INF = std::numeric_limits<double>::infinity();
  for (int r = 1; r <= R; r++) {
    const uint32_t a = sc.off[r], b = sc.off[r + 1];
    const int n = (int)(b - a);
    if (n == 0) continue;
    const uint64_t* FE = sc.FE.data() + a;
    uint64_t* T = sc.tmp.data();
    const uint32_t jE0 = jof(FE[0]), jE1 = jof(FE[n - 1]);
    const double rge = dval_of(FE[n - 1]) - dval_of(FE[0]);
    crowd[jE0] = INF; crowd[jE1] = INF;
    for (int q = 1; q < n - 1; q++)
      crowd[jof(FE[q])] = rge > 0 ? (dval_of(FE[q + 1]) - dval_of(FE[q - 1])) / rge : 0.0;
    for (int i = 0; i < n; i++) T[i] = FE[n - 1 - i];
    for (int i = 0; i < n;) {
      uint64_t key = T[i] >> 16;
      int e = i + 1;
      while (e < n && (T[e] >> 16) == key) e++;
      if (e - i > 1) std::reverse(T + i, T + e);
      i = e;
    }
    const double rgc = pv(T[n - 1]) - pv(T[0]);
    for (int p = 0; p < n; p++) {
      const uint32_t j = jof(T[p]);
      const bool endC = (p == 0 || p == n - 1);
      if (endC || j == jE0 || j == jE1) crowd[j] = INF;
      else if (rgc > 0) crowd[j] += (pv(T[p + 1]) - pv(T[p - 1])) / rgc;
    }
  }
}

inline void choose_survivors(int n, const double* crowd, RankScratch& sc, uint32_t* surv, std::vector<uint32_t>& tmp) {
  int k = 0;
  const int R = (int)sc.off.size() - 2;
  for (int r = 1; r <= R && k < n; r++) {
    const uint32_t a = sc.off[r], b = sc.off[r + 1];
    const int f = (int)(b - a);
    if (k + f <= n) {
      for (uint32_t q = a; q < b; q++) surv[k++] = jof(sc.FE[q]);
    } else {
      tmp.clear();
      for (uint32_t q = a; q < b; q++) tmp.push_back(jof(sc.FE[q]));
      std::sort(tmp.begin(), tmp.end(), [&](uint32_t x, uint32_t y) {
        if (crowd[x] != crowd[y]) return crowd[x] > crowd[y];
        return x < y;
      });
      for (int q = 0; k < n; q++) surv[k++] = tmp[q];
    }
  }
}

struct SelScratch {
  RankScratch rs;
  std::vector<uint32_t> pidx, surv, tmp, mtmp;
  std::vector<float> pe;
  std::vector<uint16_t> pc;
  std::vector<uint16_t> prank;
  std::vector<double> pcrowd;
};

inline void migrants(int d, const Ind* par, uint32_t* out, std::vector<uint32_t>& tmp) {
  tmp.resize(DECK);
  for (int i = 0; i < DECK; i++) tmp[i] = (uint32_t)(d * DECK + i);
  std::partial_sort(tmp.begin(), tmp.begin() + MIGR, tmp.end(), [&](uint32_t a, uint32_t b) {
    if (par[a].rank != par[b].rank) return par[a].rank < par[b].rank;
    if (par[a].crowd != par[b].crowd) return par[a].crowd > par[b].crowd;
    return a < b;
  });
  for (int k = 0; k < MIGR; k++) out[k] = tmp[k];
}

inline void select_deme(int d, const Ind* par, Ind* chi, Ind* nxt, const uint32_t* migtab, SelScratch& S) {
  S.pidx.clear();
  for (int i = d * DECK; i < (d + 1) * DECK; i++) S.pidx.push_back((uint32_t)i);
  for (int i = d * DECK; i < (d + 1) * DECK; i++) {
    Ind& c = chi[i];
    if (!c.valid) continue;
    float s = 0.0f;
    for (int q = 0; q < NRES; q++) s += c.res[q];
    c.err = s;
    c.cx = (uint16_t)cplx(c.g);
    S.pidx.push_back((uint32_t)(NPOP + i));
  }
  int pd = (d + NDEME - 1) % NDEME;
  for (int k = 0; k < MIGR; k++) S.pidx.push_back(migtab[pd * MIGR + k]);
  int m = (int)S.pidx.size();
  S.pe.resize(m); S.pc.resize(m); S.prank.resize(m); S.pcrowd.resize(m); S.surv.resize(DECK);
  for (int q = 0; q < m; q++) {
    uint32_t p = S.pidx[q];
    const Ind& src = p < (uint32_t)NPOP ? par[p] : chi[p - NPOP];
    S.pe[q] = src.err;
    S.pc[q] = (uint16_t)src.cx;
  }
  rank_crowd(m, S.pe.data(), S.pc.data(), S.prank.data(), S.pcrowd.data(), S.rs);
  choose_survivors(DECK, S.pcrowd.data(), S.rs, S.surv.data(), S.tmp);
  for (int s = 0; s < DECK; s++) {
    uint32_t q = S.surv[s];
    uint32_t p = S.pidx[q];
    const Ind& src = p < (uint32_t)NPOP ? par[p] : chi[p - NPOP];
    Ind& o = nxt[d * DECK + s];
    o = src;
    o.rank = S.prank[q];
    o.crowd = S.pcrowd[q];
    o.err = S.pe[q];
    o.valid = 1;
  }
}

inline void select_all(Pool& pool, const Ind* par, Ind* chi, Ind* nxt, uint32_t* migtab, SelScratch* S) {
  pool.run(NDEME, [&](int d) { migrants(d, par, migtab + d * MIGR, S[d].mtmp); });
  pool.run(NDEME, [&](int d) { select_deme(d, par, chi, nxt, migtab, S[d]); });
}

struct LexState {
  std::vector<float> sv, eps, col;
  std::vector<uint16_t> si;
  struct Scr { std::vector<float> v, dev; std::vector<uint32_t> ka, kb; std::vector<uint16_t> ia, ib; };
  std::vector<Scr> scr;
  LexState() : sv((mem_lean() & 1) ? 0 : (size_t)NDEME * NCASE * DECK), eps((mem_lean() & 1) ? 0 : (size_t)NDEME * NCASE), col((mem_lean() & 1) ? 0 : (size_t)NDEME * NCASE * DECK),
               si((mem_lean() & 1) ? 0 : (size_t)NDEME * NCASE * DECK), scr((mem_lean() & 1) ? 0 : NDEME) {}
};

inline void lex_sort_col(const float* v, int n, float* out_v, uint16_t* out_i, LexState::Scr& s) {
  s.ka.resize(n); s.kb.resize(n); s.ia.resize(n); s.ib.resize(n);
  for (int i = 0; i < n; i++) {
    s.ka[i] = fkey(v[i]); s.ia[i] = (uint16_t)i;
  }
  uint32_t* ka = s.ka.data(); uint32_t* kb = s.kb.data();
  uint16_t* ia = s.ia.data(); uint16_t* ib = s.ib.data();
  static const int SH[3] = {0, 11, 22};
  for (int p = 0; p < 3; p++) {
    uint32_t cnt[2048];
    std::memset(cnt, 0, sizeof cnt);
    for (int i = 0; i < n; i++) cnt[(ka[i] >> SH[p]) & 2047]++;
    uint32_t sum = 0;
    for (int d = 0; d < 2048; d++) { uint32_t c = cnt[d]; cnt[d] = sum; sum += c; }
    for (int i = 0; i < n; i++) {
      uint32_t pos = cnt[(ka[i] >> SH[p]) & 2047]++;
      kb[pos] = ka[i]; ib[pos] = ia[i];
    }
    std::swap(ka, kb); std::swap(ia, ib);
  }
  for (int i = 0; i < n; i++) {
    out_v[i] = funkey(ka[i]); out_i[i] = ia[i];
  }
}

inline void lex_prepare(Pool& pool, const Ind* par, LexState& L) {
  pool.run(NDEME, [&](int d) {
    LexState::Scr& s = L.scr[d];
    s.v.resize(DECK); s.dev.resize(DECK);
    for (int c = 0; c < NCASE; c++) {
      size_t off = ((size_t)d * NCASE + c) * DECK;
      for (int i = 0; i < DECK; i++) s.v[i] = par[d * DECK + i].res[c];
      std::memcpy(&L.col[off], s.v.data(), sizeof(float) * DECK);
      lex_sort_col(s.v.data(), DECK, &L.sv[off], &L.si[off], s);
      float med = L.sv[off + DECK / 2];
      for (int i = 0; i < DECK; i++) s.dev[i] = std::fabs(s.v[i] - med);
      std::nth_element(s.dev.begin(), s.dev.begin() + DECK / 2, s.dev.end());
      L.eps[d * NCASE + c] = s.dev[DECK / 2];
    }
  });
}

inline int lex_one(const LexState& L, int d, Rng& r, uint16_t* alive) {
  int perm[NCASE];
  for (int c = 0; c < NCASE; c++) perm[c] = c;
  for (int t = 0; t < LEX_CASES; t++) {
    int k = t + (int)r.below(NCASE - t);
    int tmp = perm[t]; perm[t] = perm[k]; perm[k] = tmp;
  }
  int c0 = perm[0];
  size_t off0 = ((size_t)d * NCASE + c0) * DECK;
  const float* sv0 = &L.sv[off0];
  const uint16_t* si0 = &L.si[off0];
  float thr = sv0[0] + L.eps[d * NCASE + c0];
  int na = 0;
  for (int k = 0; k < DECK && sv0[k] <= thr; k++) alive[na++] = si0[k];
  for (int t = 1; t < LEX_CASES && na > 1; t++) {
    int c = perm[t];
    const float* col = &L.col[((size_t)d * NCASE + c) * DECK];
    float mn = std::numeric_limits<float>::infinity();
    for (int i = 0; i < na; i++) { float v = col[alive[i]]; if (v < mn) mn = v; }
    thr = mn + L.eps[d * NCASE + c];
    int nk = 0;
    for (int i = 0; i < na; i++) if (col[alive[i]] <= thr) alive[nk++] = alive[i];
    na = nk;
  }
  return alive[na == 1 ? 0 : (int)r.below(na)];
}

inline void lex_select_all(Pool& pool, const LexState& L, uint64_t seed, uint64_t step, uint32_t* sel) {
  pool.run(NDEME, [&](int d) {
    Rng r(mix(mix(seed ^ 0xA5A5A5A5ull, step), (uint64_t)d));
    uint16_t alive[DECK];
    for (int i = 0; i < DECK; i++)
      sel[d * DECK + i] = (uint32_t)(d * DECK + lex_one(L, d, r, alive));
  });
}

inline void init_rank(Pool& pool, Ind* par, SelScratch* S) {
  pool.run(NDEME, [&](int d) {
    SelScratch& s = S[d];
    s.pe.resize(DECK); s.pc.resize(DECK); s.prank.resize(DECK); s.pcrowd.resize(DECK);
    for (int i = 0; i < DECK; i++) { s.pe[i] = par[d * DECK + i].err; s.pc[i] = (uint16_t)par[d * DECK + i].cx; }
    rank_crowd(DECK, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    for (int i = 0; i < DECK; i++) { par[d * DECK + i].rank = s.prank[i]; par[d * DECK + i].crowd = s.pcrowd[i]; }
  });
}

inline void fill_rows(std::vector<float>& tab, int nrows, uint64_t seed) {
  tab.resize((size_t)nrows * NCASE);
  Rng r(seed);
  for (auto& v : tab) { double g = r.gauss(); v = (float)(g * g); }
}


struct State {
  std::vector<Ind> par, chi, nxt;
  std::vector<uint32_t> migtab, rb, fl, lexsel;
  BreedBuf BB;
  std::vector<SelScratch> sel;
  LexState lex;
  State() : par(NPOP), chi(NPOP), nxt(NPOP), migtab((size_t)NDEME * MIGR), rb(NDEME), fl(NDEME),
            lexsel((mem_lean() & 1) ? 0 : NPOP), sel(NDEME) {}
};

inline void init_pop(Pool& pool, State& S, const Cand* seeds, int nseeds, const float* rows, int nrows, uint64_t seed) {
  for (int i = 0; i < NPOP; i++) {
    Ind& I = S.par[i];
    I.g = seeds[i % nseeds];
    I.valid = 1;
    Rng r(mix(seed ^ 0x51ull, (uint64_t)i));
    std::memcpy(I.res, rows + (size_t)r.below(nrows) * NCASE, sizeof I.res);
    float s = 0; for (int q = 0; q < NRES; q++) s += I.res[q];
    I.err = s; I.cx = (uint16_t)cplx(I.g);
  }
  init_rank(pool, S.par.data(), S.sel.data());
}

inline void synth_children(State& S, const float* rows, int nrows, uint64_t seed, uint64_t step) {
  for (int i = 0; i < NPOP; i++) {
    Rng r(mix(mix(seed ^ 0x77ull, step), (uint64_t)i));
    std::memcpy(S.chi[i].res, rows + (size_t)r.below(nrows) * NCASE, sizeof S.chi[i].res);
  }
}

inline void breed_step(Pool& pool, State& S, uint64_t seed, uint64_t step) {
  precompute_ext(pool, S.par.data());
  breed(pool, S.par.data(), S.chi.data(), seed, step, S.rb.data(), S.fl.data());
}

inline void select_step(Pool& pool, State& S) {
  select_all(pool, S.par.data(), S.chi.data(), S.nxt.data(), S.migtab.data(), S.sel.data());
  std::swap(S.par, S.nxt);
}

}
