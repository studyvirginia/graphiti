#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <atomic>

namespace h313b {

inline bool g_on = false;
inline std::atomic<long> g_fix{0};

static inline uint32_t scan_word(uint32_t c, uint32_t SL, uint32_t TT) {
  if (c < 64u) return c * SL;
  if (c == 1022u) return (1u << 31) | (3u << 29) | (TT + 48u);
  if (c == 900u) return 2u << 29;
  if (c == 901u) return (3u << 29) | (TT + 48u);
  if (c == 2000u) return (3u << 29) | (TT + 49u);
  return (1u << 29) | (c - 64u);
}
static inline uint32_t gn_word(uint32_t c, uint32_t CD48, uint32_t SB, uint32_t NS) {
  if (c < 64u) return c * CD48;
  if (c == 1022u) return 1u << 31;
  if (c == 900u) return 1u << 29;
  if (c == 901u) return (1u << 30) | (SB + NS);
  if (c == 2000u) return (1u << 30) | (SB + NS + 1u);
  return (1u << 30) | (SB + (c - 64u));
}

static inline uint32_t map_code(uint32_t c, uint32_t NS) { return c == 1022u ? 0xFFFEu : c == 900u ? 0xFFFDu : c == 901u ? 64u + NS : c == 2000u ? 65u + NS : c; }
static inline uint32_t unmap_code(uint32_t c, uint32_t NS) {
  if (c < 64u) return c;
  if (c == 0xFFFEu) return 1022u;
  if (c == 0xFFFDu) return 900u;
  if (c == 64u + NS) return 901u;
  if (c == 65u + NS) return 2000u;
  return c < 64u + NS ? c : 0xFFFFFFFFu;
}

struct ScanP { uint32_t SL, TT, NS; };
static inline ScanP scan_p(int nth, int nslot, int nS) {
  const uint32_t SS = nth == 0 ? 2u : 16u, SL = SS * 48u;
  const uint32_t nsc = nS <= 4 ? 4u : nS <= 8 ? 8u : nS <= 16 ? 16u : 32u;
  return ScanP{SL, (uint32_t)nslot * SL + SS * nsc, (uint32_t)nS};
}
struct GnP { uint32_t CD48, SB, NS; };
static inline GnP gn_p(uint32_t ghw) {
  const uint32_t nth = ghw >> 24, CD = nth <= 4 ? 1u + (nth + 1u) / 2u : 3u;
  return GnP{CD * 48u, (ghw >> 8) & 0xFFFFu, (ghw & 255u) - 2u};
}

static inline bool expand_scan(const uint32_t* tk, size_t ntok, ScanP p, uint32_t* out4) {
  bool ok = true;
  for (size_t i = 0; i < ntok; i++) {
    const uint32_t x = tk[2 * i], y = tk[2 * i + 1], zo = x & 255u, d = x >> 8;
    if (zo >= 16u || !(d < 40u || d == 63u)) ok = false;
    out4[4 * i] = zo; out4[4 * i + 1] = d == 63u ? 0xFFFFFFFFu : d * p.SL;
    const uint32_t A = unmap_code(y & 0xFFFFu, p.NS), B = unmap_code(y >> 16, p.NS);
    if (A == 0xFFFFFFFFu || B == 0xFFFFFFFFu) ok = false;
    out4[4 * i + 2] = scan_word(A, p.SL, p.TT); out4[4 * i + 3] = scan_word(B, p.SL, p.TT);
  }
  return ok;
}
static inline bool expand_gn(const uint32_t* tk, size_t ntok, GnP p, uint32_t* out4) {
  bool ok = true;
  for (size_t i = 0; i < ntok; i++) {
    const uint32_t x = tk[2 * i], y = tk[2 * i + 1], zo = x & 255u, d = x >> 8;
    if (zo >= 16u || !(d < 40u || d == 63u)) ok = false;
    out4[4 * i] = zo; out4[4 * i + 1] = d == 63u ? 0xFFFFFFFFu : d * p.CD48;
    const uint32_t A = unmap_code(y & 0xFFFFu, p.NS), B = unmap_code(y >> 16, p.NS);
    if (A == 0xFFFFFFFFu || B == 0xFFFFFFFFu) ok = false;
    out4[4 * i + 2] = gn_word(A, p.CD48, p.SB, p.NS); out4[4 * i + 3] = gn_word(B, p.CD48, p.SB, p.NS);
  }
  return ok;
}

static inline uint32_t scan_code(uint32_t w, ScanP p) {
  if (w >> 31) return w == ((1u << 31) | (3u << 29) | (p.TT + 48u)) ? 1022u : 0xFFFFFFFFu;
  switch (w >> 29) {
    case 0: return (p.SL && w % p.SL == 0 && w / p.SL < 64u) ? w / p.SL : 0xFFFFFFFFu;
    case 1: return (w & 0x1FFFFFFFu) + 64u;
    case 2: return w == (2u << 29) ? 900u : 0xFFFFFFFFu;
    default: return w == ((3u << 29) | (p.TT + 48u)) ? 901u : w == ((3u << 29) | (p.TT + 49u)) ? 2000u : 0xFFFFFFFFu;
  }
}
static inline bool contract_scan(const uint32_t* p4, size_t ntok, ScanP p, std::vector<uint32_t>& tk) {
  tk.resize(2 * ntok);
  for (size_t i = 0; i < ntok; i++) {
    const uint32_t zo = p4[4 * i], y = p4[4 * i + 1];
    if (zo >= 16u) return false;
    uint32_t d;
    if (y == 0xFFFFFFFFu) d = 63u; else { if (y % p.SL) return false; d = y / p.SL; if (d >= 40u) return false; }
    const uint32_t A = map_code(scan_code(p4[4 * i + 2], p), p.NS), B = map_code(scan_code(p4[4 * i + 3], p), p.NS);
    if (A > 0xFFFFu || B > 0xFFFFu) return false;
    tk[2 * i] = zo | (d << 8); tk[2 * i + 1] = A | (B << 16);
  }
  std::vector<uint32_t> chk(4 * ntok);
  expand_scan(tk.data(), ntok, p, chk.data());
  for (size_t i = 0; i < 4 * ntok; i++) if (chk[i] != p4[i]) return false;
  return true;
}

static inline uint32_t gn_code(uint32_t w, GnP p) {
  if (w >> 31) return w == (1u << 31) ? 1022u : 0xFFFFFFFFu;
  if ((w >> 30) & 1u) { const uint32_t o = w & 0x1FFFFFFFu; if (w >> 29 & 1u) return 0xFFFFFFFFu; if (o == p.SB + p.NS) return 901u; if (o == p.SB + p.NS + 1u) return 2000u; return o >= p.SB ? o - p.SB + 64u : 0xFFFFFFFFu; }
  if ((w >> 29) & 1u) return w == (1u << 29) ? 900u : 0xFFFFFFFFu;
  return (w % p.CD48 == 0 && w / p.CD48 < 64u) ? w / p.CD48 : 0xFFFFFFFFu;
}
static inline bool contract_gn(const uint32_t* p4, size_t ntok, GnP p, std::vector<uint32_t>& tk) {
  tk.resize(2 * ntok);
  for (size_t i = 0; i < ntok; i++) {
    const uint32_t zo = p4[4 * i], y = p4[4 * i + 1];
    if (zo >= 16u) return false;
    uint32_t d;
    if (y == 0xFFFFFFFFu) d = 63u; else { if (y % p.CD48) return false; d = y / p.CD48; if (d >= 40u) return false; }
    const uint32_t A = map_code(gn_code(p4[4 * i + 2], p), p.NS), B = map_code(gn_code(p4[4 * i + 3], p), p.NS);
    if (A > 0xFFFFu || B > 0xFFFFu) return false;
    tk[2 * i] = zo | (d << 8); tk[2 * i + 1] = A | (B << 16);
  }
  std::vector<uint32_t> chk(4 * ntok);
  expand_gn(tk.data(), ntok, p, chk.data());
  for (size_t i = 0; i < 4 * ntok; i++) if (chk[i] != p4[i]) return false;
  return true;
}

}
