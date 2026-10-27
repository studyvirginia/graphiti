#pragma once
#include <cstdint>
#include <cstring>
#include <vector>
namespace h313a {
static inline size_t nprog(const Tree& t) { return t.progp ? (size_t)t.nprogp : t.prog.size(); }
static inline const uint32_t* prog(const Tree& t) { return t.progp ? t.progp : t.prog.data(); }
static inline const uint16_t* segs(const Tree& t) { return t.segsp ? t.segsp : t.segs; }
static inline void own_prog(Tree& t) {
  if (!t.progp) return;
  t.prog.assign(t.progp, t.progp + t.nprogp); std::memcpy(t.segs, t.segsp, sizeof t.segs);
  t.progp = nullptr; t.segsp = nullptr; t.nprogp = 0;
}
}
#define H313_NPROG(t) (h313a::nprog(t))
#define H313_PROG(t) (h313a::prog(t))
#define H313_SEGS(t, out) do { const uint16_t* sg_ = h313a::segs(t); for (int q_ = 0; q_ < 36; q_++) (out)[q_] = sg_[q_]; } while (0)
#define H313_NSP(e) ((e).nsp())
#define H313_PSP(e) ((e).psp())
#define H313_NGP(e) ((e).ngp())
#define H313_PGP(e) ((e).pgp())
