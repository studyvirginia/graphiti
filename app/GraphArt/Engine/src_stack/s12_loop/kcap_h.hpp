#pragma once
#include "kcap_fmt.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
inline int g_kcap_step = 0;
static inline const char* kcap_dir() { static const char* d = getenv("FJN_KCAP"); return (d && d[0]) ? d : nullptr; }
static inline int kcap_want_step() { static const int s = getenv("FJN_KCAP_STEP") ? atoi(getenv("FJN_KCAP_STEP")) : 3; return s; }
struct KcapState { FILE* f = nullptr; uint32_t ncalls = 0; };
static inline KcapState& kcap_state() { static KcapState s; return s; }
static inline bool kcap_active() { return kcap_dir() && g_kcap_step == kcap_want_step(); }
struct KcapDispIn { const char* name; int Q, Tn; size_t ebo, plo, tgb, ntg, tpg, nmem, per; };
static inline void kcap_write_call(const std::vector<KcapDispIn>& ds, const void* const bufs[8], const size_t lens[8], const void* cf, size_t cfBytes, size_t ebLen, size_t outLen,
                                   size_t npr, bool closed, int maxit, float flr, const float* yv96, const float* Smv, bool res48b, bool cfh, size_t r48s) {
  KcapState& S = kcap_state();
  if (!S.f) { std::string p = std::string(kcap_dir()) + "/kcap.bin"; S.f = fopen(p.c_str(), "wb"); if (!S.f) { fprintf(stderr, "KCAP: cannot open %s\n", p.c_str()); return; }
    KcHdr h; std::memset(&h, 0, sizeof h); std::memcpy(h.magic, "KCAP0002", 8); h.ncalls = 0; h.flags = (res48b ? 1u : 0u) | (cfh ? 2u : 0u); fwrite(&h, sizeof h, 1, S.f); }
  KcCall c; std::memset(&c, 0, sizeof c); c.magic = KC_CALL_MAGIC; c.step = g_kcap_step; c.npr = (uint32_t)npr; c.ndisp = (uint32_t)ds.size(); c.closed = closed; c.maxit = (uint32_t)maxit; c.flr = flr; c.r48s = (uint32_t)r48s;
  for (int i = 0; i < 8; i++) c.sz[i] = lens[i];
  c.sz[8] = cfBytes; c.sz[9] = ebLen; c.sz[10] = outLen;
  std::memcpy(c.yv, yv96, sizeof c.yv); std::memcpy(c.Smv, Smv, sizeof c.Smv);
  { unsigned a = 1, b = 1, cc = 0, e2 = 1, f = 0; if (const char* ev = getenv("FJN_REPS")) sscanf(ev, "%u,%u,%u,%u,%u", &a, &b, &cc, &e2, &f); c.rpv[0] = a; c.rpv[1] = b; c.rpv[2] = cc; c.rpv[3] = e2; c.bwr = (float)f; }
  fwrite(&c, sizeof c, 1, S.f);
  for (int i = 0; i < 8; i++) if (lens[i]) fwrite(bufs[i], 1, lens[i], S.f);
  if (cfBytes) fwrite(cf, 1, cfBytes, S.f);
  for (const KcapDispIn& d : ds) { KcDisp k; std::memset(&k, 0, sizeof k); std::strncpy(k.name, d.name, 31); k.Q = (uint32_t)d.Q; k.Tn = (uint32_t)d.Tn; k.ebo = d.ebo; k.plo = d.plo; k.tgb = d.tgb; k.ntg = (uint32_t)d.ntg; k.tpg = (uint32_t)d.tpg; k.nmem = (uint32_t)d.nmem; k.per = (uint32_t)d.per; fwrite(&k, sizeof k, 1, S.f); }
  S.ncalls++; fflush(S.f);
  fseek(S.f, offsetof(KcHdr, ncalls), SEEK_SET); fwrite(&S.ncalls, 4, 1, S.f); fseek(S.f, 0, SEEK_END); fflush(S.f);
}
