#include "gpu_fail_h.hpp"
#pragma once
#include <tuple>
#if __has_include("fjn_tmax_h950.h")
#include "fjn_tmax_h950.h"
#endif
#ifndef FJN_TMAX
#define FJN_TMAX 12
#endif
static inline bool fjn_tpgrule_on() { static const bool v = getenv("FJN_TPGRULE") && atoi(getenv("FJN_TPGRULE")) == 1; return v; }
struct FrTreeRef { int Q = 0, nth = 0, nS = 0, nslot = 0; const uint32_t* prog = nullptr; size_t nprog = 0; const uint16_t* sp = nullptr; size_t nsp = 0; const uint16_t* segs = nullptr; };
struct FrPairIn { FrTreeRef x, y; const float* thx = nullptr; const float* thy = nullptr; };
struct FrOut { float Fx = INFINITY, Fy = INFINITY, stalled = 0, it = 0, trials = 0, dual = 0, plain = 0; float th[2][32]; float a[2][16]; float u[10]; float traj[14]; float raw[128]; float res48[96]; float cf = 0; int status = 0; };
struct FrRunStats { long pairs = 0, nokernel = 0, tg = 0, enc = 0, ran = 0, b_q = 0, b_tn = 0, b_pso = 0, b_slot = 0, b_tpg = 0; int bins = 0; double gpu_ms = 0, host_ms = 0, wait_ms = 0, cb_start = 0, cb_end = 0; };

struct FjnPending {
  struct Disp { id<MTLComputePipelineState> ps; int Q, Tn; std::vector<int> mem; size_t tpg, ntg, tgb, per, ebo, plo; id<MTLBuffer> pl; };
  bool valid = false; id<MTLCommandBuffer> cb = nil; id<MTLBuffer> bOut = nil; std::vector<Disp> ds; std::vector<id<MTLBuffer>> keep;
  std::chrono::steady_clock::time_point tw; std::vector<FrOut>* out = nullptr; FrRunStats* S = nullptr;
};
static const int g_wb_par_on = getenv("WB_PAR") ? atoi(getenv("WB_PAR")) : 0;
static std::function<void(int, const std::function<void(int)>&)> g_wb_par;
static const bool g_cfh = getenv("FJN_CISFUSE_HOST") != nullptr;
static const int g_cfchk = (g_cfh && getenv("FJN_CISFUSE_CHECK")) ? std::max(0, atoi(getenv("FJN_CISFUSE_CHECK"))) : 0;
static const size_t R48S = g_cfh ? 97 : 96;
static std::vector<float> g_cf_data; static id<MTLBuffer> g_cf_buf = nil; static bool g_cf_on = false;
static const bool g_res48b = getenv("FJN_RES48_HOST") != nullptr || g_cfh;
static void fjn_unpack(id<MTLCommandBuffer> cb, id<MTLBuffer> bOut, const std::vector<FjnPending::Disp>& ds, std::vector<FrOut>& out, FrRunStats& S) {
  using clk_ = std::chrono::steady_clock;
  auto msnow = [](clk_::time_point a) { return std::chrono::duration<double, std::milli>(clk_::now() - a).count(); };
  if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "fjn_run command buffer failed: %s\n", cb.error ? cb.error.description.UTF8String : "?"); gpu_fail(); }
  S.gpu_ms += (cb.GPUEndTime - cb.GPUStartTime) * 1e3; S.cb_start = cb.GPUStartTime; S.cb_end = cb.GPUEndTime;
  const auto tr0 = clk_::now();
  const float* po = (const float*)bOut.contents;
  auto unpack1 = [&](int k) {
    const float* o = po + (size_t)k * 128; FrOut& r = out[(size_t)k];
    static const bool full_ = getenv("FJN_UNPACK_FULL") && atoi(getenv("FJN_UNPACK_FULL")) == 1;
    if (full_) std::memcpy(r.raw, o, 128 * 4);
    r.Fx = o[0]; r.Fy = o[1]; r.stalled = o[2]; r.it = o[3]; r.trials = o[4]; r.dual = o[5]; r.plain = o[6];
    std::memcpy(r.th[0], o + 8, 32 * 4); std::memcpy(r.th[1], o + 40, 32 * 4); std::memcpy(r.a[0], o + 72, 16 * 4); std::memcpy(r.a[1], o + 88, 16 * 4); std::memcpy(r.u, o + 104, 10 * 4);
    if (full_) std::memcpy(r.traj, o + 114, 14 * 4);
    if (g_res48b && (!g_cfh || g_cfchk > 0)) std::memcpy(r.res48, po + out.size() * 128 + (size_t)k * R48S, 96 * 4);
    if (g_cfh) r.cf = po[out.size() * 128 + (size_t)k * R48S + 96];
    r.status = 1;
  };
  if (g_wb_par_on && g_wb_par) {
    std::vector<int> ks; for (const FjnPending::Disp& d : ds) ks.insert(ks.end(), d.mem.begin(), d.mem.end());
    const int NR = std::min(8, (int)ks.size() / 64 + 1);
    g_wb_par(NR, [&](int r) { const size_t n = ks.size(); for (size_t i = n * (size_t)r / NR; i < n * (size_t)(r + 1) / NR; i++) unpack1(ks[i]); });
    S.ran += (long)ks.size();
  } else {
    for (const FjnPending::Disp& d : ds) for (int k : d.mem) { unpack1(k); S.ran++; }
  }
  S.host_ms += msnow(tr0);
}

static const int g_fcpw = getenv("FCP_WORKERS") ? atoi(getenv("FCP_WORKERS")) : 0;
static const int g_kgp = (getenv("KGP") ? atoi(getenv("KGP")) : 0) || g_fcpw > 0 || (!(getenv("FJN_HOSTPAR") && atoi(getenv("FJN_HOSTPAR")) == 0));
static std::function<void(int, const std::function<void(int)>&)> g_kgp_par;
static inline bool fr_enc_pair(int nth, int nS, int nslot, int CD, const uint32_t* prg, size_t nprg, const uint16_t* tsg,
                               std::vector<uint32_t>& pP, std::vector<uint32_t>& pD, uint32_t* segsP, uint32_t* segsD,
                               uint32_t& hwP, uint32_t& hwD, uint32_t& SBp_o, uint32_t& SBd_o) {
  const uint32_t SSp = 48u, SSd = 48u * (uint32_t)CD, NSC = (uint32_t)nS + 2;
  const uint32_t SBp = (uint32_t)nslot * SSp, SBd = (uint32_t)nslot * SSd; SBp_o = SBp; SBd_o = SBd;
  auto enc = [&](uint32_t c, uint32_t SS, uint32_t SB) -> uint32_t { if (c == 1022) return 1u << 31; if (c < 64) return c * SS; if (c == 900) return 1u << 29;
    if (c == 901) return (1u << 30) | (SB + (uint32_t)nS); return (1u << 30) | (SB + (c - 64)); };
  const uint32_t TWOp = (1u << 30) | (SBp + (uint32_t)nS + 1), TWOd = (1u << 30) | (SBd + (uint32_t)nS + 1);
  auto form = [&](uint32_t op, uint32_t a, uint32_t b, uint32_t SS, uint32_t SB, uint32_t TWO, uint32_t& zo, uint32_t& A, uint32_t& B) -> bool {
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a, SS, SB); B = enc(901, SS, SB); } else if (op == 3) { zo = 0; A = enc(a, SS, SB); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b, SS, SB); B = enc(a, SS, SB); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a, SS, SB); B = enc(b, SS, SB); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a, SS, SB); B = enc(901, SS, SB); }
    else { zo = 6; A = enc(a, SS, SB); B = enc(901, SS, SB); }
    return true; };
  pP.clear(); pD.clear();
  for (size_t pi = 0; pi < nprg; pi++) { const uint32_t w = prg[pi]; const uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22;
    uint32_t zoP, AP, BP, zoD, AD, BD;
    if (!form(op, a, b, SSp, SBp, TWOp, zoP, AP, BP)) return false;
    if (!form(op, a, b, SSd, SBd, TWOd, zoD, AD, BD)) return false;
    pP.push_back(zoP); pP.push_back(d == 63 ? 0xFFFFFFFFu : d * SSp); pP.push_back(AP); pP.push_back(BP);
    pD.push_back(zoD); pD.push_back(d == 63 ? 0xFFFFFFFFu : d * SSd); pD.push_back(AD); pD.push_back(BD); }
  for (int q = 0; q < 18; q++) { segsP[q] = tsg[q]; segsD[q] = tsg[q]; }
  for (int q = 0; q < 18; q++) { segsP[18 + q] = enc(tsg[18 + q], SSp, SBp); segsD[18 + q] = enc(tsg[18 + q], SSd, SBd); }
  if (SBp >= 65536 || SBd >= 65536) return false;
  hwP = NSC | (SBp << 8) | ((uint32_t)nth << 24); hwD = NSC | (SBd << 8) | ((uint32_t)nth << 24);
  return true;
}

static void fjn_run(id<MTLDevice> dev, id<MTLCommandQueue> cq, const std::function<id<MTLComputePipelineState>(const std::string&)>& pipe,
                    const std::function<id<MTLCommandBuffer>()>& mkcb, const std::vector<FrPairIn>& pairs, const float* yv96, float flr, bool closed, int maxit,
                    const fjn::WarpMat& WM, std::vector<FrOut>& out, FrRunStats& S, size_t tpgcap = 8, FjnPending* pend = nullptr) {
  using clk_ = std::chrono::steady_clock;
  auto msnow = [](clk_::time_point a) { return std::chrono::duration<double, std::milli>(clk_::now() - a).count(); };
  const auto th0 = clk_::now();
  const size_t NPR = pairs.size(); out.assign(NPR, FrOut()); S.pairs += (long)NPR;
  if (NPR == 0) return;
  struct RecP { fjn::Rec r; std::vector<uint16_t> sp; bool ok = false; std::vector<uint32_t> kprog; std::vector<uint16_t> kseg; int kq = 0, knth = 0, knS = 0, knslot = 0, kcd = 0; };
  std::vector<RecP> recs; recs.reserve(NPR * 2);
  std::unordered_multimap<uint64_t, int> rmap; rmap.reserve(NPR * 3);
  auto getrec = [&](const FrTreeRef& tr, int CD) -> int {
    uint64_t key = 1469598103934665603ull;
    auto mix = [&](uint64_t v) { key = (key ^ v) * 1099511628211ull; };
    for (size_t i = 0; i < tr.nprog; i++) mix((uint64_t)tr.prog[i]); mix((uint64_t)tr.nprog); mix((uint64_t)CD); mix((uint64_t)tr.Q); mix((uint64_t)tr.nth); mix((uint64_t)tr.nS); mix((uint64_t)tr.nslot); mix((uint64_t)tr.nsp);
    for (int i = 0; i < 36; i++) mix((uint64_t)tr.segs[i]);
    for (size_t i = 0; i < tr.nsp; i++) mix((uint64_t)tr.sp[i]);
    auto rg = rmap.equal_range(key);
    for (auto it = rg.first; it != rg.second; ++it) { const RecP& c = recs[(size_t)it->second];
      if (c.kq == tr.Q && c.knth == tr.nth && c.knS == tr.nS && c.knslot == tr.nslot && c.kcd == CD && c.kprog.size() == tr.nprog && c.sp.size() == tr.nsp &&
          (tr.nprog == 0 || std::memcmp(c.kprog.data(), tr.prog, tr.nprog * 4) == 0) && std::memcmp(c.kseg.data(), tr.segs, 36 * 2) == 0 && (tr.nsp == 0 || std::memcmp(c.sp.data(), tr.sp, tr.nsp * 2) == 0))
        return it->second; }
    RecP rp; fjn::Rec& r = rp.r; rp.kq = tr.Q; rp.knth = tr.nth; rp.knS = tr.nS; rp.knslot = tr.nslot; rp.kcd = CD; rp.kprog.assign(tr.prog, tr.prog + tr.nprog); rp.kseg.assign(tr.segs, tr.segs + 36); r.nth = tr.nth; r.Q = tr.Q; r.nS = tr.nS;
    if (g_kgp) { rp.sp.assign(tr.sp, tr.sp + tr.nsp); rp.ok = false; recs.push_back(std::move(rp)); const int id = (int)recs.size() - 1; rmap.emplace(key, id); return id; }
    struct TT { int nth, nS, nslot; } tt{tr.nth, tr.nS, tr.nslot};
    uint32_t SBp = 0, SBd = 0; const std::pair<const uint32_t*, size_t> pg{tr.prog, tr.nprog};
    static const bool g_epair = getenv("FJN_EPAIR") && atoi(getenv("FJN_EPAIR")) == 1;
    bool ok = g_epair ? fr_enc_pair(tr.nth, tr.nS, tr.nslot, CD, tr.prog, tr.nprog, tr.segs, r.tokP, r.tokD, r.sgP, r.sgD, r.hP[3], r.hD[3], SBp, SBd)
                      : fjn::enc_tree(tt, 48u, pg, tr.segs, r.tokP, r.sgP, r.hP[3], SBp) && fjn::enc_tree(tt, 48u * (uint32_t)CD, pg, tr.segs, r.tokD, r.sgD, r.hD[3], SBd);
    if (ok) { r.hP[1] = r.hD[1] = (uint32_t)tr.nS | ((uint32_t)tr.Q << 16); r.hP[0] = r.hP[2] = r.hD[0] = r.hD[2] = 0;
      r.perP = SBp + (uint32_t)(tr.nS + 2); r.perD = SBd + (uint32_t)CD * (uint32_t)(tr.nS + 2); rp.sp.assign(tr.sp, tr.sp + tr.nsp); }
    rp.sp.assign(tr.sp, tr.sp + tr.nsp); rp.ok = ok; recs.push_back(std::move(rp)); const int id = (int)recs.size() - 1; rmap.emplace(key, id); return id; };
  auto encrec = [&](size_t id) {
    RecP& rp = recs[id]; fjn::Rec& r = rp.r; const int CD = rp.kcd;
    struct TT { int nth, nS, nslot; } tt{rp.knth, rp.knS, rp.knslot};
    uint32_t SBp = 0, SBd = 0; const std::pair<const uint32_t*, size_t> pg{rp.kprog.data(), rp.kprog.size()};
    bool ok = fjn::enc_tree(tt, 48u, pg, rp.kseg.data(), r.tokP, r.sgP, r.hP[3], SBp) && fjn::enc_tree(tt, 48u * (uint32_t)CD, pg, rp.kseg.data(), r.tokD, r.sgD, r.hD[3], SBd);
    if (ok) { r.hP[1] = r.hD[1] = (uint32_t)rp.knS | ((uint32_t)rp.kq << 16); r.hP[0] = r.hP[2] = r.hD[0] = r.hD[2] = 0;
      r.perP = SBp + (uint32_t)(rp.knS + 2); r.perD = SBd + (uint32_t)CD * (uint32_t)(rp.knS + 2); }
    rp.ok = ok; };
  static const bool g_pq = !(getenv("FJN_PQ") && atoi(getenv("FJN_PQ")) == 0);
  const size_t v5_lpnk = closed ? 10 : 9;
  auto v5_lpneed = [&](const FrPairIn& p) -> size_t { if (!(fjn::g_v5 & 1)) return 0; const size_t n = (size_t)p.x.nth + (size_t)p.y.nth + v5_lpnk; return n > 16 ? n * (n + 1) / 2 : 0; };
  std::map<std::tuple<int, int, int>, std::vector<int>> bins;
  std::vector<int> recx(NPR, -1), recy(NPR, -1);
  std::vector<int> gA, gB;
  if (g_kgp) { gA.resize(NPR); gB.resize(NPR);
    for (size_t k = 0; k < NPR; k++) { const FrPairIn& p = pairs[k]; const int CD = fjn::geom_CD(std::max(p.x.nth, p.y.nth)); gA[k] = getrec(p.x, CD); gB[k] = getrec(p.y, CD); }
    if (g_kgp_par) g_kgp_par((int)recs.size(), [&](int i) { encrec((size_t)i); });
    else { const size_t NR = recs.size(); const int nw = g_fcpw > 1 ? std::min<int>(g_fcpw, (int)std::max<size_t>(1, NR / 64)) : 1;
      if (nw <= 1) { for (size_t i = 0; i < NR; i++) encrec(i); }
      else { std::vector<std::thread> th; for (int w = 0; w < nw; w++) th.emplace_back([&, w] { for (size_t i = (size_t)w; i < NR; i += (size_t)nw) encrec(i); }); for (auto& t : th) t.join(); } }
  }
  for (size_t k = 0; k < NPR; k++) {
    const FrPairIn& p = pairs[k]; const int Qc = std::max(p.x.Q, p.y.Q), Tc = std::max(p.x.nth, p.y.nth); const int CD = fjn::geom_CD(Tc);
    const int a = g_kgp ? gA[k] : getrec(p.x, CD), b = g_kgp ? gB[k] : getrec(p.y, CD);
    if (!recs[(size_t)a].ok || !recs[(size_t)b].ok) { S.enc++; continue; }
    size_t pr = 0; { const fjn::Rec& ra = recs[(size_t)a].r; const fjn::Rec& rb = recs[(size_t)b].r; pr = std::max(std::max(ra.perP, ra.perD), std::max(rb.perP, rb.perD)); }
    pr = std::max(pr, v5_lpneed(p)); pr = (pr + 3) / 4 * 4; int pq = 0;
    if (g_pq) { size_t v = 64; while (v < pr) v += v / 5 + 1; pq = (int)((v + 3) / 4 * 4); }
    recx[k] = a; recy[k] = b; bins[std::make_tuple(Qc, Tc, pq)].push_back((int)k);
  }
  size_t ntokP = 0, ntokD = 0, nsp = 0;
  for (RecP& rp : recs) if (rp.ok) { ntokP += rp.r.tokP.size() / 4; ntokD += rp.r.tokD.size() / 4; nsp += rp.sp.size(); }
  std::vector<uint32_t> tokP(ntokP * 4 + 4), tokD(ntokD * 4 + 4), segv(recs.size() * 72 + 72), hdrv(recs.size() * 8 + 8); std::vector<uint16_t> spv(nsp + 4);
  { size_t oP = 0, oD = 0, oS = 0;
    for (size_t r = 0; r < recs.size(); r++) { RecP& rp = recs[r]; if (!rp.ok) continue; fjn::Rec& q = rp.r;
      std::memcpy(&tokP[oP * 4], q.tokP.data(), q.tokP.size() * 4); std::memcpy(&tokD[oD * 4], q.tokD.data(), q.tokD.size() * 4); if (!rp.sp.empty()) std::memcpy(&spv[oS], rp.sp.data(), rp.sp.size() * 2);
      uint32_t* h = &hdrv[r * 8]; h[0] = (uint32_t)oP; h[1] = q.hP[1]; h[2] = (uint32_t)oS; h[3] = q.hP[3]; h[4] = (uint32_t)oD; h[5] = q.hD[1]; h[6] = (uint32_t)oS; h[7] = q.hD[3];
      std::memcpy(&segv[r * 72], q.sgP, 36 * 4); std::memcpy(&segv[r * 72 + 36], q.sgD, 36 * 4);
      oP += q.tokP.size() / 4; oD += q.tokD.size() / 4; oS += rp.sp.size(); } }
  std::vector<uint32_t> prv(NPR * 2); std::vector<float> th0v(NPR * 64, 0.0f);
  for (size_t k = 0; k < NPR; k++) { if (recx[k] < 0) continue; prv[k * 2] = (uint32_t)recx[k]; prv[k * 2 + 1] = (uint32_t)recy[k];
    if (pairs[k].thx) std::memcpy(&th0v[k * 64], pairs[k].thx, (size_t)pairs[k].x.nth * 4); if (pairs[k].thy) std::memcpy(&th0v[k * 64 + 32], pairs[k].thy, (size_t)pairs[k].y.nth * 4); }
  std::vector<float> Smv(48 * 9 + 9); std::memcpy(Smv.data(), WM.S, 48 * 9 * 4); std::memcpy(Smv.data() + 48 * 9, WM.u0, 9 * 4);
  static const bool fjnUntracked = getenv("FJN_UNTRACKED") != nullptr && getenv("FJN_UNTRACKED")[0] == '1';
  const MTLResourceOptions inOpt = MTLResourceStorageModeShared | (fjnUntracked ? MTLResourceHazardTrackingModeUntracked : MTLResourceHazardTrackingModeTracked);
  id<MTLBuffer> bP = [dev newBufferWithBytes:tokP.data() length:tokP.size() * 4 options:inOpt];
  id<MTLBuffer> bD = [dev newBufferWithBytes:tokD.data() length:tokD.size() * 4 options:inOpt];
  id<MTLBuffer> bH = [dev newBufferWithBytes:hdrv.data() length:hdrv.size() * 4 options:inOpt];
  id<MTLBuffer> bSp = [dev newBufferWithBytes:spv.data() length:spv.size() * 2 options:inOpt];
  id<MTLBuffer> bSg = [dev newBufferWithBytes:segv.data() length:segv.size() * 4 options:inOpt];
  id<MTLBuffer> bPr = [dev newBufferWithBytes:prv.data() length:prv.size() * 4 options:inOpt];
  id<MTLBuffer> bTh = [dev newBufferWithBytes:th0v.data() length:th0v.size() * 4 options:inOpt];
  id<MTLBuffer> bOut = [dev newBufferWithLength:NPR * (128 + (g_res48b ? R48S : 0)) * 4 options:MTLResourceStorageModeShared];
  for (id<MTLBuffer> b_ : {bP, bD, bH, bSp, bSg, bPr, bTh, bOut}) memdiag::metal_note(b_, (size_t)b_.length, false);
  std::memset(bOut.contents, 0xFF, NPR * 128 * 4);
  struct Disp { id<MTLComputePipelineState> ps; int Q, Tn; std::vector<int> mem; size_t tpg, ntg, tgb, per, ebo, plo; id<MTLBuffer> pl; };
  std::vector<Disp> ds; size_t ebTot = 0; std::vector<uint32_t> plAll;
  struct BWork { std::tuple<int, int, int> key; std::vector<int> mem; bool merged; };
  static const bool g_binmerge = getenv("FJN_BINMERGE") && atoi(getenv("FJN_BINMERGE")) != 0;
  std::vector<BWork> work;
  if (g_binmerge) { std::map<std::pair<int, int>, std::vector<int>> g;
    for (auto& kv : bins) { auto& v = g[std::make_pair(std::get<0>(kv.first), std::get<1>(kv.first))]; v.insert(v.end(), kv.second.begin(), kv.second.end()); }
    for (auto& g1 : g) {
      const int Q = g1.first.first, Tn = g1.first.second; size_t per = 0;
      for (int k : g1.second) { const fjn::Rec& a = recs[(size_t)recx[(size_t)k]].r; const fjn::Rec& b = recs[(size_t)recy[(size_t)k]].r;
        per = std::max<size_t>(per, std::max(std::max(a.perP, a.perD), std::max(b.perP, b.perD))); per = std::max<size_t>(per, v5_lpneed(pairs[(size_t)k])); }
      per = (per + 3) / 4 * 4; per = std::max<size_t>(per, (fjn::g_v5 & 1) ? (size_t)0 : (size_t)fjn::per_min(Tn)); per = (per + 3) / 4 * 4;
      const size_t slot = ((size_t)fjn::slotc_floats(Q, Tn) + (fjn::g_lay == 2 ? std::max<size_t>((size_t)fjn::ovl_floats(Q, Tn), 2 * per) : 2 * per)) * 4;
      if (slot <= 32768) work.push_back({std::make_tuple(Q, Tn, 0), g1.second, true});
      else for (auto& kv : bins) if (std::get<0>(kv.first) == Q && std::get<1>(kv.first) == Tn) work.push_back({kv.first, kv.second, false}); } }
  else for (auto& kv : bins) work.push_back({kv.first, kv.second, false});
  static const bool g_lpt = getenv("FJN_LPT") && atoi(getenv("FJN_LPT")) != 0;
  if (g_lpt) std::stable_sort(work.begin(), work.end(), [](const BWork& a, const BWork& b) {
    auto w = [](const BWork& x) { const long Q = std::get<0>(x.key), T = std::get<1>(x.key); return (Q + 2) * (T + 1) * (T + 2); };
    const long wa = w(a), wb = w(b); if (wa != wb) return wa > wb; return std::get<2>(a.key) > std::get<2>(b.key); });
  for (size_t wi = 0; wi < work.size(); wi++) {
    const int Q = std::get<0>(work[wi].key), Tn = std::get<1>(work[wi].key); std::vector<int>& mem = work[wi].mem;
    if (Q < 1 || Q > 16) { S.nokernel += (long)mem.size(); S.b_q += (long)mem.size(); continue; }
    if (Tn > FJN_TMAX) { S.nokernel += (long)mem.size(); S.b_tn += (long)mem.size(); continue; }
    id<MTLComputePipelineState> ps = pipe("fjn_q" + std::to_string(Q) + "_t" + std::to_string(Tn));
    if (!ps) { S.nokernel += (long)mem.size(); S.b_pso += (long)mem.size(); continue; }
    size_t per = 0; for (int k : mem) { const fjn::Rec& a = recs[(size_t)recx[(size_t)k]].r; const fjn::Rec& b = recs[(size_t)recy[(size_t)k]].r;
      per = std::max<size_t>(per, std::max(std::max(a.perP, a.perD), std::max(b.perP, b.perD))); per = std::max<size_t>(per, v5_lpneed(pairs[(size_t)k])); }
    per = (per + 3) / 4 * 4; if (std::get<2>(work[wi].key) > 0 && !(fjn::g_v5 & 2)) per = (size_t)std::get<2>(work[wi].key);
    per = std::max<size_t>(per, (fjn::g_v5 & 1) ? (size_t)4 : (size_t)fjn::per_min(Tn)); per = (per + 3) / 4 * 4;
    const size_t slot = ((size_t)fjn::slotc_floats(Q, Tn) + (fjn::g_lay == 2 ? std::max<size_t>((size_t)fjn::ovl_floats(Q, Tn), 2 * per) : 2 * per)) * 4;
    if (slot > 32768) { S.tg += (long)mem.size(); S.b_slot += (long)mem.size(); continue; }
    size_t tpg = std::min<size_t>(tpgcap, std::max<size_t>(1, 32768 / slot)); tpg = std::min<size_t>(tpg, (size_t)ps.maxTotalThreadsPerThreadgroup / 32);
    if (fjn_tpgrule_on()) {
      const size_t cap = std::min<size_t>(mem.size(), (size_t)ps.maxTotalThreadsPerThreadgroup / 32);
      size_t best = 1; long bestSim = -1;
      for (size_t t = 1; t <= std::max<size_t>(1, cap); t++) {
        const size_t g = ((size_t)slot * t + 15) & ~(size_t)15; if (g > 32768 || g == 0) break;
        const long sim = std::min<long>(96, (long)(65536 / g) * (long)t);
        if (sim >= bestSim) { bestSim = sim; best = t; }
      }
      tpg = best;
    }
    if (tpg < 1) { S.nokernel += (long)mem.size(); S.b_tpg += (long)mem.size(); continue; }
    Disp d; d.ps = ps; d.Q = Q; d.Tn = Tn; d.mem = mem; d.tpg = tpg; d.ntg = (mem.size() + tpg - 1) / tpg; d.tgb = (tpg * slot + 15) / 16 * 16 + (getenv("FJN_TGPAD") ? (size_t)atoi(getenv("FJN_TGPAD")) : 0);    d.per = per; d.ebo = ebTot; d.plo = plAll.size() * 4;
    ebTot += ((mem.size() + tpg) * (size_t)fjn::eb_stride(Q, Tn) * 4 + 255) / 256 * 256;
    if (getenv("FJN_VERBOSE")) printf("fjn bin Q%d T%d: pairs %zu per %zu slot %zu B tpg %zu ntg %zu maxTotal %lu execWidth %lu\n", Q, Tn, mem.size(), per, slot, tpg, d.ntg, (unsigned long)ps.maxTotalThreadsPerThreadgroup, (unsigned long)ps.threadExecutionWidth);
    plAll.insert(plAll.end(), mem.begin(), mem.end());
    ds.push_back(std::move(d));
  }
  if (ds.empty()) { S.host_ms += msnow(th0); return; }
  id<MTLBuffer> bEB = [dev newBufferWithLength:std::max<size_t>(ebTot, 256) options:MTLResourceStorageModePrivate];
  id<MTLBuffer> bPlAll = [dev newBufferWithBytes:(plAll.empty() ? (const void*)"\0\0\0\0" : (const void*)plAll.data()) length:std::max<size_t>(plAll.size() * 4, 4) options:MTLResourceStorageModeShared];
  memdiag::metal_note(bEB, (size_t)bEB.length, true); memdiag::metal_note(bPlAll, (size_t)bPlAll.length, false);
  for (Disp& d : ds) d.pl = bPlAll;
  S.host_ms += msnow(th0);
  id<MTLCommandBuffer> cb = mkcb();
  id<MTLComputeCommandEncoder> en = [cb computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
  if (getenv("TL_LABELS") && getenv("TL_LABELS")[0] == '1') en.label = [NSString stringWithFormat:@"fjn_run chunk Disp=%zu", ds.size()];
  for (Disp& d : ds) {
    [en setComputePipelineState:d.ps];
    id<MTLBuffer> bs[9] = {bP, bD, bH, bSp, bSg, nil, bPr, bTh, bOut};
    for (int k = 0; k < 9; k++) if (bs[k]) [en setBuffer:bs[k] offset:0 atIndex:k];
    if (g_res48b) [en setBuffer:bOut offset:(NSUInteger)(NPR * 128 * 4) atIndex:15];
    if (g_cfh) { if (!g_cf_buf) { std::vector<float> z(400, 0.0f); const std::vector<float>& s0 = (g_cf_data.size() == 400 || g_cf_data.size() == 496) ? g_cf_data : z; g_cf_buf = [dev newBufferWithBytes:s0.data() length:s0.size() * 4 options:MTLResourceStorageModeShared]; } [en setBuffer:g_cf_buf offset:0 atIndex:16]; }
    [en setBytes:yv96 length:96 * 4 atIndex:5];
    static uint32_t rpv[4] = {1, 1, 0, 1}; static float bwr = 0; static bool rpinit = false;
    if (!rpinit) { rpinit = true; if (const char* ev = getenv("FJN_REPS")) { unsigned a = 1, b = 1, c = 0, e2 = 1, f = 0; sscanf(ev, "%u,%u,%u,%u,%u", &a, &b, &c, &e2, &f); rpv[0] = a; rpv[1] = b; rpv[2] = c; rpv[3] = e2; bwr = (float)f; } }
    float tf[4] = {(float)d.Q * 9.0f * ldexpf(1.0f, -24), flr, bwr, 0.0f}; [en setBytes:tf length:16 atIndex:9];
    [en setBytes:rpv length:16 atIndex:14];
    [en setBuffer:bEB offset:d.ebo atIndex:10];
    [en setBytes:Smv.data() length:Smv.size() * 4 atIndex:11];
    [en setBuffer:d.pl offset:d.plo atIndex:12];
    uint32_t wq[4] = {(uint32_t)d.mem.size(), (uint32_t)d.per, closed ? 1u : 0u, (uint32_t)maxit}; [en setBytes:wq length:16 atIndex:13];
    [en setThreadgroupMemoryLength:d.tgb atIndex:0];
    [en dispatchThreadgroups:MTLSizeMake(d.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * d.tpg, 1, 1)];
    S.bins++;
  }
  [en endEncoding]; const auto tw = clk_::now(); [cb commit];
  if (pend) { pend->cb = cb; pend->bOut = bOut; pend->ds.clear(); for (Disp& d : ds) pend->ds.push_back({d.ps, d.Q, d.Tn, d.mem, d.tpg, d.ntg, d.tgb, d.per, d.ebo, d.plo, d.pl}); pend->keep = {bP, bD, bH, bSp, bSg, bPr, bTh, bEB, bOut, bPlAll}; pend->tw = tw; pend->out = &out; pend->S = &S; pend->valid = true; S.host_ms += msnow(th0); return; }
  [cb waitUntilCompleted]; S.wait_ms += msnow(tw);
  std::vector<FjnPending::Disp> dv; for (Disp& d : ds) dv.push_back({d.ps, d.Q, d.Tn, d.mem, d.tpg, d.ntg, d.tgb, d.per, d.ebo, d.plo, d.pl});
  fjn_unpack(cb, bOut, dv, out, S);
  S.host_ms += msnow(th0);
}

static void fjn_collect(FjnPending& pend) {
  if (!pend.valid) return;
  using clk_ = std::chrono::steady_clock;
  const auto tw = clk_::now(); [pend.cb waitUntilCompleted]; pend.S->wait_ms += std::chrono::duration<double, std::milli>(clk_::now() - tw).count();
  fjn_unpack(pend.cb, pend.bOut, pend.ds, *pend.out, *pend.S);
  pend.valid = false; pend.keep.clear(); pend.bOut = nil; pend.cb = nil; pend.ds.clear();
}
