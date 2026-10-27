#pragma once
#include <cstdint>
#include <vector>

struct GNInputs {
  const uint32_t* tok;
  const uint32_t* hdr;
  const uint16_t* sprog;
  const uint32_t* segs;
  const float* yv;
  const float* tlin;
  const float* th0;
  uint32_t per;
};

struct GNResult {
  float th[32];
  float a[16];
  float rss0, rcur;
  uint32_t steps, npp, ndp;
};

inline void gn_tlin(float* t) { for (int i = 0; i < 48; i++) t[i] = (float)i * (1.0f / 47); }

class GNWork {
 public:
  GNWork();
  void polish(const GNInputs& in, uint32_t tree, int Q, float tau, float flr, int itmax, GNResult& out);
  void roots_debug(const GNInputs& in, uint32_t tree, const float* TH, int C, int Q, float* out);
 private:
  struct Tree { uint32_t x, nS, z, nsc2, SB, nth; };
  void eval(const GNInputs& in, const Tree& T, const uint32_t* sg, const float* TH, const float* tvec, int C, int Q);
  void exec(const uint32_t* w4, int C);
  void exec417(const uint32_t* w4, int C);
  void exec421(const uint32_t* w4, int C);
  void exec836(const uint32_t* w4, int C);
  struct R644 { uint32_t op; uint8_t sa, sb; const float* a; const float* b; float* acc; float* cp; };
  struct K644 { const uint32_t* tok = nullptr; const uint32_t* sg = nullptr; int C = -1, Q = -1; uint32_t x = 0, cd = 0, nsc2 = 0, SB = 0;
                const float *R = nullptr, *S = nullptr, *TS = nullptr, *acc = nullptr, *acc2 = nullptr; bool h313 = false;
                std::vector<R644> rec; std::vector<uint32_t> beg; std::vector<const float*> root; std::vector<uint8_t> rsc; float* accp_end = nullptr; };
  K644 k644_[2];
  void eval644(const GNInputs& in, const Tree& T, const uint32_t* sg, int C, int Q);
  std::vector<uint32_t> smask_; std::vector<uint64_t> slotst_; uint32_t rootm_[16]; bool smask_ok_ = false, rootm_valid_ = false; uint32_t accm_ = 0, gen941_ = 0, full941_ = 0;
  inline bool opmask941(uint32_t w, uint32_t& m);
  bool exec940(const uint32_t* w4, int C);
  bool eval_impl(const GNInputs& in, const Tree& T, const uint32_t* sg, const float* TH, const float* tvec, int C, int Q, bool masked);
  const float* load(uint32_t w, float* buf, int C);
  float* rdst_ = nullptr; size_t rstr_ = 0;
  size_t vslot(uint32_t off) const { return (size_t)(cd_ == 1u ? off / 48u : cd_ == 2u ? off / 96u : cd_ == 3u ? off / 144u : off / (cd_ * 48u)) * (size_t)(cstride_ * 48u); }
  std::vector<float> R_, S_, TS_, acc_, acc2_, roots_, a_, b_; float* accp_ = nullptr;
  int sb842_ = 1, tsz842_ = 1; bool ts0_842_ = false;
  uint32_t nsc2_ = 0, SB_ = 0, per_ = 0, cd_ = 1, cstride_ = 1;
};
