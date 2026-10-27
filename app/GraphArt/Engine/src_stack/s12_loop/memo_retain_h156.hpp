#pragma once
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
#include "h156_flags.h"
#include "nf_flat_h16_i13_h156.hpp"
#include "compile_memo_i17_h156.h"

namespace h156 {

struct RetainStats {
  long nf_ev = 0, cc_ev = 0, enc_ev = 0;
  long nf_kept = 0, cc_kept = 0, enc_kept = 0;
  uint32_t thr = 0;
};

static inline uint64_t nf_block_bytes(const nf16::BHdr* hd) {
  const uint64_t nb = (uint64_t)hd->nabs * 4 + (uint64_t)hd->pn_len * sizeof(skn::PNode);
  return ((uint64_t)sizeof(nf16::BHdr) + nb + 7) & ~(uint64_t)7;
}

static inline double cc_entry_bytes(const cc14::Ent& e, const cmi6::Slot& s, const uint32_t* ar = nullptr) {
  double b = (double)sizeof(cc14::Ent) + 4.0 * e.nkey + 4.0 * e.nprog + 2.0 * e.nsp + (double)cmi6::slot_hdr_bytes();
  if (s.st == 2) b += 4.0 * (double)(ar ? cmi6::slot_words(s, ar) : (uint64_t)s.nsp4 + s.ngp4 + (cmi6::g_h913 ? cmi6::tabw(s.cmode, s.Qs) : 0u));
  return b;
}

static inline uint32_t gen_of(uint32_t t, uint32_t nt) { return t < nt ? t : nt - 1; }

static inline void compact_nf(nf16::Shared& S, uint32_t thr, RetainStats& out) {
  std::vector<uint64_t> keep;
  long live = 0;
  for (uint32_t p = 0; p < nf16::TB; p++) {
    const uint64_t v = S.tab[p].load(std::memory_order_relaxed);
    if (!v) continue;
    live++;
    const uint64_t off = (uint64_t)((uint32_t)v - 1u);
    if (((const nf16::BHdr*)(S.ar + off))->touch >= thr) keep.push_back(off);
  }
  std::sort(keep.begin(), keep.end());
  for (uint32_t p = 0; p < nf16::TB; p++) S.tab[p].store(0, std::memory_order_relaxed);
  uint64_t dst = 0;
  for (uint64_t off : keep) {
    const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + off);
    const uint64_t bytes = nf_block_bytes(hd);
    const uint64_t h = skn::NFMemo::hkey(hd->key);
    if (dst != off) memmove(S.ar + dst, S.ar + off, (size_t)bytes);
    uint32_t p = (uint32_t)(h & (nf16::TB - 1));
    while (S.tab[p].load(std::memory_order_relaxed) != 0) p = (p + 1) & (nf16::TB - 1);
    S.tab[p].store(((uint64_t)(uint32_t)(h >> 32) << 32) | (dst + 1), std::memory_order_relaxed);
    dst += bytes;
  }
  S.bump.store(dst, std::memory_order_relaxed);
  S.nent.store((uint32_t)keep.size(), std::memory_order_relaxed);
  out.nf_ev = live - (long)keep.size();
  out.nf_kept = (long)keep.size();
}

static inline void compact_cc(cc14::Shared& C, cmi6::Memo& M, uint32_t thr, RetainStats& out) {
  const uint32_t ne = std::min<uint32_t>(C.nent.load(std::memory_order_relaxed), cc14::Shared::MAXENT);
  std::vector<uint8_t> live(ne, 0);
  for (uint32_t p = 0; p < cc14::Shared::TB; p++) {
    const uint32_t v = C.tab[p].load(std::memory_order_relaxed);
    if (v && v - 1 < ne) live[v - 1] = 1;
  }
  long nlive = 0, nlive_enc = 0;
  std::vector<uint32_t> keep;
  for (uint32_t i = 0; i < ne; i++) {
    if (!live[i]) continue;
    nlive++;
    if (M.at(i).st == 2) nlive_enc++;
    if (C.ent[i].touch >= thr) keep.push_back(i);
  }
  const uint32_t K = (uint32_t)keep.size();

  std::vector<uint32_t> nk_off(ne, 0), np_off(ne, 0), ns_off(ne, 0), nt_off(ne, 0);
  auto slide = [](const std::vector<uint32_t>& lst, auto offf, auto cntf, auto* base, std::vector<uint32_t>& nw) -> uint64_t {
    std::vector<uint32_t> o(lst);
    std::sort(o.begin(), o.end(), [&](uint32_t a, uint32_t b) { return offf(a) < offf(b); });
    uint64_t dst = 0;
    for (uint32_t i : o) {
      const uint64_t src = offf(i), n = cntf(i);
      nw[i] = (uint32_t)dst;
      if (dst != src && n) memmove(base + dst, base + src, (size_t)(n * sizeof(*base)));
      dst += n;
    }
    return dst;
  };
  const uint64_t tot_k = slide(keep, [&](uint32_t i) { return (uint64_t)C.ent[i].key_off; },
                               [&](uint32_t i) { return (uint64_t)C.ent[i].nkey; }, C.key, nk_off);
  const uint64_t tot_p = slide(keep, [&](uint32_t i) { return (uint64_t)C.ent[i].prog_off; },
                               [&](uint32_t i) { return (uint64_t)C.ent[i].nprog; }, C.prog, np_off);
  const uint64_t tot_s = slide(keep, [&](uint32_t i) { return (uint64_t)C.ent[i].sp_off; },
                               [&](uint32_t i) { return (uint64_t)C.ent[i].nsp; }, C.sp, ns_off);
  std::vector<uint32_t> enc_keep;
  for (uint32_t i : keep) if (M.at(i).st == 2) enc_keep.push_back(i);
  const uint64_t tot_t = slide(enc_keep, [&](uint32_t i) { return (uint64_t)M.at(i).sp_off; },
                               [&](uint32_t i) { return cmi6::slot_words(M.at(i), M.arena); }, M.arena, nt_off);

  for (uint32_t r = 0; r < K; r++) {
    const uint32_t i = keep[r];
    if (r != i) { C.ent[r] = C.ent[i]; std::memcpy((void*)&M.at(r), (const void*)&M.at(i), M.stride); }
    cc14::Ent& e = C.ent[r];
    e.key_off = nk_off[i]; e.prog_off = np_off[i]; e.sp_off = ns_off[i];
    cmi6::Slot& s = M.at(r);
    if (s.st == 2) { s.sp_off = nt_off[i]; s.gp_off = nt_off[i] + s.nsp4; }
  }
  if (ne > K) memset((void*)&M.at(K), 0, M.stride * (ne - K));

  for (uint32_t p = 0; p < cc14::Shared::TB; p++) C.tab[p].store(0, std::memory_order_relaxed);
  for (uint32_t r = 0; r < K; r++) {
    uint32_t p = (uint32_t)(C.ent[r].h & (cc14::Shared::TB - 1));
    while (C.tab[p].load(std::memory_order_relaxed) != 0) p = (p + 1) & (cc14::Shared::TB - 1);
    C.tab[p].store(r + 1, std::memory_order_relaxed);
  }
  C.nent.store(K, std::memory_order_relaxed);
  C.nkey.store((uint32_t)tot_k, std::memory_order_relaxed);
  C.nprog.store((uint32_t)tot_p, std::memory_order_relaxed);
  C.nsp.store((uint32_t)tot_s, std::memory_order_relaxed);
  M.bump.store(tot_t, std::memory_order_relaxed);
  long enc_kept = 0;
  for (uint32_t r = 0; r < K; r++) if (M.at(r).st == 2) enc_kept++;
  M.nst.store((uint32_t)enc_kept, std::memory_order_relaxed);
  out.cc_ev = nlive - (long)K;
  out.cc_kept = (long)K;
  out.enc_ev = nlive_enc - enc_kept;
  out.enc_kept = enc_kept;
}

static inline void retain(nf16::Shared& S, cc14::Shared& C, cmi6::Memo& M, uint32_t st, double cap_bytes, RetainStats& out) {
  const uint32_t NT = st + 2;
  std::vector<double> hist(NT, 0.0);
  for (uint32_t p = 0; p < nf16::TB; p++) {
    const uint64_t v = S.tab[p].load(std::memory_order_relaxed);
    if (!v) continue;
    const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + (uint64_t)((uint32_t)v - 1u));
    hist[gen_of(hd->touch, NT)] += (double)nf_block_bytes(hd);
  }
  const uint32_t ne = std::min<uint32_t>(C.nent.load(std::memory_order_relaxed), cc14::Shared::MAXENT);
  std::vector<uint8_t> live(ne, 0);
  for (uint32_t p = 0; p < cc14::Shared::TB; p++) {
    const uint32_t v = C.tab[p].load(std::memory_order_relaxed);
    if (v && v - 1 < ne) live[v - 1] = 1;
  }
  for (uint32_t i = 0; i < ne; i++)
    if (live[i]) hist[gen_of(C.ent[i].touch, NT)] += cc_entry_bytes(C.ent[i], M.at(i), M.arena);

  const double target = 0.5 * cap_bytes;
  uint32_t thr = st > 0 ? st - 1 : 0;
  double cum = 0;
  for (int64_t g = (int64_t)st - 1; g >= 0; g--) {
    cum += hist[(size_t)g];
    if (cum > target) break;
    thr = (uint32_t)g;
  }
  out.thr = thr;
  compact_nf(S, thr, out);
  compact_cc(C, M, thr, out);
}

}
