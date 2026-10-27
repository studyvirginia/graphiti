#pragma once
#include "nf_flat_h16_i13_h156.hpp"
#include "h156_flags.h"

namespace nf16h264 {
using namespace skn;

template <class C>
static inline bool tree_seg_h264(nf16::Handle& H, const C& cd, int t, const NFMemo::Key& k, uint64_t hk, nf16::Block& B, nf16::Child& ch, uint64_t& r) {
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
    const bool good = nf16::miss_compute(H, k, hk, cd.op[t], cd.c[t], n);
    if (good) {
      const uint32_t nw = (uint32_t)H.P.nw;
      B.w.push_back(nw);
      const size_t r0 = B.w.size();
      B.w.resize(r0 + nw);
      if (nw) memcpy(B.w.data() + r0, H.P.w, sizeof(uint32_t) * nw);
      ch.roff[t] = (uint32_t)r0; ch.rlen[t] = nw;
      r = nf16::hash_words(B.w.data() + r0, nw);
      ok = H.P.Q <= 16;
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

template <class C>
static void chunk_h264(nf16::Handle& H, const C* const* cands, const uint8_t* valid, int b0, int b1, nf16::Block& B, nf16::Child* out) {
  B.w.clear();
  const int m = b1 - b0;
  NFMemo::Key K[64][2];
  uint64_t HK[64][2];
  for (int j = 0; j < m; j++) {
    const int i = b0 + j;
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
    uint64_t r[2] = {0, 0};
    ch.ok[0] = tree_seg_h264(H, cd, 0, K[j][0], HK[j][0], B, ch, r[0]) ? 1 : 0;
    ch.ok[1] = tree_seg_h264(H, cd, 1, K[j][1], HK[j][1], B, ch, r[1]) ? 1 : 0;
    ch.klen = (uint32_t)(B.w.size() - ch.koff);
    ch.r0 = r[0]; ch.r1 = r[1];
    ch.kh = nf16::sm(r[0] ^ nf16::sm(r[1] + ch.klen));
  }
}

}
