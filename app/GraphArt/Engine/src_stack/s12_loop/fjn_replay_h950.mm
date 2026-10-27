#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include "fjn_h950.hpp"
#include "fjn_run_h950.inc"

struct SideStore { int nth = 0, Q = 0, nS = 0, nslot = 0; std::vector<float> th0; std::vector<uint32_t> prog; std::vector<uint16_t> sp; uint16_t segs[36]; };
struct PairStore { int child = 0; SideStore s[2]; double Fx = 0, Fy = 0; double itrec = 0, fix = 0; };

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "usage: %s dump.txt fjn.metallib reps [tpgcap] [maxit] [minQ,minT filter]\n", argv[0]); return 2; }
  const int reps = atoi(argv[3]); const size_t tpgcap = argc > 4 ? (size_t)atoi(argv[4]) : 8; int maxit = argc > 5 ? atoi(argv[5]) : 60;
  std::ifstream f(argv[1]); if (!f) { fprintf(stderr, "cannot open dump\n"); return 1; }
  std::string line; bool closed = false; double floor_ = 0; float yv[96] = {0};
  std::vector<PairStore> ps; int sd = -1;
  while (std::getline(f, line)) {
    std::istringstream is(line); std::string key; is >> key;
    if (key == "CLOSED") { int c; is >> c; closed = c != 0; }
    else if (key == "FLOOR") is >> floor_;
    else if (key == "Y") { for (int i = 0; i < 96; i++) is >> yv[i]; }
    else if (key == "PAIR") { PairStore p; is >> p.child; ps.push_back(p); sd = -1; }
    else if (key == "SIDE") { std::string k; int s_; is >> s_; sd = s_; SideStore& S = ps.back().s[sd];
      while (is >> k) { if (k == "nth") is >> S.nth; else if (k == "Q") is >> S.Q; else if (k == "nS") is >> S.nS; else if (k == "nslot") is >> S.nslot; else if (k == "fixedrfin") { double d; is >> d; ps.back().fix += d; } } }
    else if (key == "TH0") { SideStore& S = ps.back().s[sd]; float v; S.th0.clear(); while (is >> v) S.th0.push_back(v); }
    else if (key == "PROG") { SideStore& S = ps.back().s[sd]; uint32_t v; S.prog.clear(); while (is >> v) S.prog.push_back(v); }
    else if (key == "SP") { SideStore& S = ps.back().s[sd]; unsigned v; S.sp.clear(); while (is >> v) S.sp.push_back((uint16_t)v); }
    else if (key == "SEGS") { SideStore& S = ps.back().s[sd]; for (int i = 0; i < 36; i++) { unsigned v; is >> v; S.segs[i] = (uint16_t)v; } }
    else if (key == "OUT") { std::string k; double v; PairStore& p = ps.back(); while (is >> k >> v) { if (k == "Fx") p.Fx = v; else if (k == "Fy") p.Fy = v; else if (k == "it") p.itrec = v; } }
  }
  if (getenv("FJN_FORCE_OPEN")) closed = false;
  printf("replay: %zu pairs, closed %d, floor %.6g\n", ps.size(), (int)closed, floor_);
  std::vector<FrPairIn> prs(ps.size());
  for (size_t i = 0; i < ps.size(); i++) {
    for (int s = 0; s < 2; s++) { SideStore& S = ps[i].s[s]; FrTreeRef r; r.Q = S.Q; r.nth = S.nth; r.nS = S.nS; r.nslot = S.nslot; r.prog = S.prog.data(); r.nprog = S.prog.size(); r.sp = S.sp.data(); r.nsp = S.sp.size(); r.segs = S.segs;
      if (s == 0) { prs[i].x = r; prs[i].thx = S.th0.data(); } else { prs[i].y = r; prs[i].thy = S.th0.data(); } }
  }
  id<MTLDevice> dev = MTLCreateSystemDefaultDevice(); id<MTLCommandQueue> cq = [dev newCommandQueue];
  NSError* e = nil; id<MTLLibrary> lib = [dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:argv[2]]] error:&e];
  if (!lib) { fprintf(stderr, "lib: %s\n", e.description.UTF8String); return 1; }
  std::map<std::string, id<MTLComputePipelineState>> psc;
  auto pipe = [&](const std::string& nm) -> id<MTLComputePipelineState> {
    auto it = psc.find(nm); if (it != psc.end()) return it->second;
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:nm.c_str()]]; id<MTLComputePipelineState> p = nil;
    if (fn) { NSError* e2 = nil; p = [dev newComputePipelineStateWithFunction:fn error:&e2]; if (!p) { fprintf(stderr, "pso %s: %s\n", nm.c_str(), e2.description.UTF8String); exit(1); } }
    psc[nm] = p; return p; };
  fjn::WarpMat WM; fjn::make_warp(WM);
  std::vector<double> gms;
  for (int r = 0; r < reps; r++) {
    std::vector<FrOut> outs; FrRunStats st;
    fjn_run(dev, cq, pipe, [&]() { return [cq commandBuffer]; }, prs, yv, (float)floor_, closed, maxit, WM, outs, st, tpgcap);
    gms.push_back(st.gpu_ms);
    if (r == 0) {
      size_t same = 0, ran = 0; double itsum = 0, dsum = 0, psum = 0, ssum = 0, tsum = 0; std::vector<double> ratio;
      for (size_t i = 0; i < ps.size(); i++) if (outs[i].status) { ran++; itsum += outs[i].it; dsum += outs[i].dual; psum += outs[i].plain; ssum += outs[i].raw[7]; tsum += outs[i].trials;
        const double Fk = (double)outs[i].Fx + outs[i].Fy, Fr = ps[i].Fx + ps[i].Fy; ratio.push_back(Fk / Fr); if (std::fabs(Fk - Fr) <= 1e-5 * std::fabs(Fr)) same++; }
      if (getenv("FJN_PAIRS")) for (size_t i = 0; i < ps.size(); i++) if (outs[i].status) printf("PAIRF %d %.9g %.9g %.9g %.9g\n", ps[i].child, (double)outs[i].Fx + outs[i].Fy, ps[i].Fx + ps[i].Fy, (double)outs[i].raw[127], ps[i].fix);
      if (const char* dof = getenv("FJN_DUMPOUT")) { FILE* fo = fopen(dof, "wb"); for (size_t i = 0; i < ps.size(); i++) { float rec[256] = {0}; const FrOut& o = outs[i]; rec[0]=o.Fx; rec[1]=o.Fy; rec[2]=o.it; rec[3]=o.dual; rec[4]=o.plain; rec[5]=o.trials; rec[6]=(float)o.status; memcpy(rec+7, o.raw, 128*4); memcpy(rec+135, o.th, 64*4); memcpy(rec+199, o.a, 32*4); memcpy(rec+231, o.u, 10*4); memcpy(rec+241, o.traj, 14*4); fwrite(rec, 4, 256, fo);} fclose(fo); }
      std::sort(ratio.begin(), ratio.end());
      printf("COUNTS per pair: it %.3f dual %.3f plain %.3f solves %.3f trials %.3f\n", ran ? itsum / ran : 0, ran ? dsum / ran : 0, ran ? psum / ran : 0, ran ? ssum / ran : 0, ran ? tsum / ran : 0);
      printf("rep0: ran %zu of %zu (nokernel %ld tg %ld enc %ld) bins %d | F same as dump (1e-5) %zu | F/Frec p10 %.4f p50 %.4f p90 %.4f | mean it %.2f dual %.2f plain %.2f | gpu %.2f ms host %.2f ms\n",
        ran, ps.size(), st.nokernel, st.tg, st.enc, st.bins, same, ratio.empty() ? 0 : ratio[ratio.size() / 10], ratio.empty() ? 0 : ratio[ratio.size() / 2], ratio.empty() ? 0 : ratio[ratio.size() * 9 / 10], ran ? itsum / ran : 0, ran ? dsum / ran : 0, ran ? psum / ran : 0, st.gpu_ms, st.host_ms);
      std::map<std::pair<int, int>, int> cls; for (auto& p : prs) cls[{std::max(p.x.Q, p.y.Q), std::max(p.x.nth, p.y.nth)}]++;
      if (getenv("FJN_CLASSES")) for (auto& kv : cls) printf("  class Q%d T%d: %d\n", kv.first.first, kv.first.second, kv.second);
    }
  }
  std::vector<double> g2 = gms; std::sort(g2.begin(), g2.end());
  printf("GPU ms per rep: min %.2f median %.2f max %.2f over %d reps (%.2f us per pair at median)\n", g2.front(), g2[g2.size() / 2], g2.back(), reps, 1e3 * g2[g2.size() / 2] / (double)ps.size());
  return 0;
}
