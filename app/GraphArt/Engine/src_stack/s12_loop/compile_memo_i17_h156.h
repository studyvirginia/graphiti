#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "compile_fused_h42_h156.h"
#ifdef I17_VERIFY
#include "compile_tree_h30.h"
#endif

namespace cmi6i17 {

#ifdef I17_VERIFY
struct VerifyCount {
  std::atomic<long> checked{0}, bad{0};
  ~VerifyCount() { fprintf(stderr, "I17_VERIFY: %ld misses checked vs h30+h34, %ld mismatches\n", checked.load(), bad.load()); }
};
static VerifyCount g_verify;

static inline bool same_tree(const Tree& A, const Tree& B) {
  return A.Q == B.Q && A.nth == B.nth && A.nS == B.nS && A.nslot == B.nslot && A.prog == B.prog && A.sp == B.sp &&
         (A.th.size() == B.th.size() && (A.th.empty() || std::memcmp(A.th.data(), B.th.data(), A.th.size() * sizeof(float)) == 0)) && std::memcmp(A.segs, B.segs, sizeof A.segs) == 0;
}
static inline bool same_enc(const cmi6::TreeEnc& A, const cmi6::TreeEnc& B) {
  return A.s_ok == B.s_ok && A.g_ok == B.g_ok && std::memcmp(A.sg, B.sg, sizeof A.sg) == 0 &&
         std::memcmp(A.gg, B.gg, sizeof A.gg) == 0 && A.ghw == B.ghw && A.gper == B.gper && A.sp4 == B.sp4 && A.gp4 == B.gp4;
}
#endif

static inline bool compile_enc_i17(cmi6::Handle& H, const View& V, Tree& T, cmi6::TreeEnc& te, double& enc_ms,
                                   const uint64_t* hpre = nullptr) {
  const size_t nk = (size_t)(V.theta - V.w);
  const uint64_t h = hpre ? *hpre : cc14::hkey(V.w, nk);

  bool hit = false;
  uint32_t p = (uint32_t)(h & (cc14::Shared::TB - 1));
  for (;;) {
    const uint32_t v = H.C.S.tab[p].load(std::memory_order_acquire);
    if (v == 0) break;
    const cc14::Ent& E = H.C.S.ent[v - 1];
    if (E.h == h && E.nkey == nk && std::memcmp(H.C.S.key + E.key_off, V.w, sizeof(uint32_t) * nk) == 0) {
      ++H.C.hits;
      if (h156::g_on) h156::touch_u32(&H.C.S.ent[v - 1].touch, h156::g_touch);
      if (!E.ok) { ++H.C.neg_hits; return false; }
      T.Q = E.Q; T.nth = E.nth; T.nS = E.nS; T.nslot = E.nslot;
      T.prog.assign(H.C.S.prog + E.prog_off, H.C.S.prog + E.prog_off + E.nprog);
      T.sp.assign(H.C.S.sp + E.sp_off, H.C.S.sp + E.sp_off + E.nsp);
      std::memcpy(T.segs, E.segs, sizeof T.segs);
      T.th.resize((size_t)V.nth);
      for (int s = 0; s < V.nth; s++) T.th[s] = V.th(s);
      hit = true;
      break;
    }
    p = (p + 1) & (cc14::Shared::TB - 1);
  }

  if (!hit) {
    ++H.C.misses;
    f38::Grid g;
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = f42::compile_fused_stage(V, T, te, g, 3, nullptr);
    enc_ms += cmi6::since_ms(t0);
#ifdef I17_VERIFY
    {
      Tree R; cmi6::TreeEnc RE;
      const bool rok = h30::compile_tree_h30(V, R);
      if (rok) enc34::encode_tree_h34(R, RE);
      g_verify.checked++;
      if (rok != ok || (ok && (!same_tree(T, R) || !same_enc(te, RE)))) {
        g_verify.bad++;
        fprintf(stderr, "I17_VERIFY MISMATCH (ok %d ref %d, Q %d)\n", (int)ok, (int)rok, V.Q);
        fprintf(stderr, " Q %d/%d nth %d/%d nS %d/%d nslot %d/%d prog %zu/%zu sp %zu/%zu th %zu/%zu segs %d enc s_ok %d/%d g_ok %d/%d sp4 %zu/%zu gp4 %zu/%zu ghw %u/%u gper %u/%u sg %d gg %d\n",
          T.Q, R.Q, T.nth, R.nth, T.nS, R.nS, T.nslot, R.nslot, T.prog.size(), R.prog.size(), T.sp.size(), R.sp.size(), T.th.size(), R.th.size(),
          (int)std::memcmp(T.segs, R.segs, sizeof T.segs), (int)te.s_ok, (int)RE.s_ok, (int)te.g_ok, (int)RE.g_ok, te.sp4.size(), RE.sp4.size(), te.gp4.size(), RE.gp4.size(),
          te.ghw, RE.ghw, te.gper, RE.gper, (int)std::memcmp(te.sg, RE.sg, sizeof te.sg), (int)std::memcmp(te.gg, RE.gg, sizeof te.gg));
        for (size_t i = 0; i < T.prog.size() && i < R.prog.size(); i++) if (T.prog[i] != R.prog[i]) { fprintf(stderr, " first prog diff at %zu: %08x vs %08x\n", i, T.prog[i], R.prog[i]); break; }
        for (size_t i = 0; i < te.sp4.size() && i < RE.sp4.size(); i++) if (te.sp4[i] != RE.sp4[i]) { fprintf(stderr, " first sp4 diff at %zu: %u vs %u\n", i, te.sp4[i], RE.sp4[i]); break; }
        for (size_t i = 0; i < T.sp.size() && i < R.sp.size(); i++) if (T.sp[i] != R.sp[i]) { fprintf(stderr, " first sp diff at %zu\n", i); break; }
        for (size_t i = 0; i < T.th.size() && i < R.th.size(); i++) if (std::memcmp(&T.th[i], &R.th[i], sizeof(float)) != 0) { uint32_t a, b; std::memcpy(&a, &T.th[i], 4); std::memcpy(&b, &R.th[i], 4); fprintf(stderr, " first th diff at %zu: bits %08x vs %08x (%g vs %g)\n", i, a, b, T.th[i], R.th[i]); break; }
        std::abort();
      }
    }
#endif
    cc14::Ent X; std::memset(&X, 0, sizeof X);
    X.h = h; X.ok = ok ? 1 : 0;
    if (ok) {
      X.Q = T.Q; X.nth = T.nth; X.nS = T.nS; X.nslot = T.nslot;
      std::memcpy(X.segs, T.segs, sizeof T.segs);
      H.C.publish(h, X, V.w, nk, T.prog.data(), T.prog.size(), T.sp.data(), T.sp.size());
    } else {
      H.C.publish(h, X, V.w, nk, nullptr, 0, nullptr, 0);
      return false;
    }
  }

  uint32_t idx = 0;
  if (!cmi6::locate(H.C.S, V.w, nk, h, idx)) {
    if (hit) { const auto t1 = std::chrono::steady_clock::now(); enc34::encode_tree_h34(T, te); enc_ms += cmi6::since_ms(t1); }
    H.enc_local++;
    return true;
  }
  cmi6::Slot& s = H.M.at(idx);
  const uint32_t st = __atomic_load_n(&s.st, __ATOMIC_ACQUIRE);
  if (st == 2) { cmi6::load_enc(s, H.M.arena, te); H.enc_hit++; return true; }
  if (hit) { const auto t1 = std::chrono::steady_clock::now(); enc34::encode_tree_h34(T, te); enc_ms += cmi6::since_ms(t1); }
  H.enc_miss++;
  uint32_t zero = 0;
  if (st == 0 && __atomic_compare_exchange_n(&s.st, &zero, 1u, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
    const uint64_t n4 = te.sp4.size(), ng = te.gp4.size(), need = n4 + ng;
    const uint64_t off = H.M.bump.fetch_add(need);
    if (off + need <= cmi6::ARENA_WORDS) {
      if (n4) std::memcpy(H.M.arena + off, te.sp4.data(), n4 * sizeof(uint32_t));
      if (ng) std::memcpy(H.M.arena + off + n4, te.gp4.data(), ng * sizeof(uint32_t));
      if (cmi6::g_h913) { fprintf(stderr, "H913: store via cmi6i17::compile_enc_i17 is not supported\n"); abort(); }
      s.s_ok = te.s_ok; s.g_ok = te.g_ok; s.ghw = te.ghw; s.gper = te.gper;
      std::memcpy(s.sg, te.sg, sizeof s.sg); std::memcpy(s.gg, te.gg, sizeof s.gg);
      s.sp_off = (uint32_t)off; s.nsp4 = (uint32_t)n4; s.gp_off = (uint32_t)(off + n4); s.ngp4 = (uint32_t)ng;
      H.M.nst.fetch_add(1, std::memory_order_relaxed);
      __atomic_store_n(&s.st, 2u, __ATOMIC_RELEASE);
    } else {
      __atomic_store_n(&s.st, 3u, __ATOMIC_RELEASE);
    }
  }
  return true;
}

}
