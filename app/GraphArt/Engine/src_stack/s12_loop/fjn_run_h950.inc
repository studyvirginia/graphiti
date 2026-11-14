#include "gpu_fail_h.hpp"
#pragma once
#include <tuple>
#include <array>
#include <deque>
#include <exception>
#include <mutex>
#include "kcap_h.hpp"
#if __has_include("fjn_tmax_h950.h")
#include "fjn_tmax_h950.h"
#endif
#ifndef FJN_TMAX
#define FJN_TMAX 12
#endif
static inline bool fjn_tpgrule_on() { static const bool v = getenv("FJN_TPGRULE") && atoi(getenv("FJN_TPGRULE")) == 1; return v; }
struct FrTreeRef { int Q = 0, nth = 0, nS = 0, nslot = 0; const uint32_t* prog = nullptr; size_t nprog = 0; const uint16_t* sp = nullptr; size_t nsp = 0; const uint16_t* segs = nullptr; };
struct FrPairIn { FrTreeRef x, y; const float* thx = nullptr; const float* thy = nullptr; int id = -1; uint64_t k1 = 0, k2 = 0; };
struct FrOut { float Fx = INFINITY, Fy = INFINITY, stalled = 0, it = 0, trials = 0, dual = 0, plain = 0; float th[2][32]; float a[2][16]; float u[10]; float traj[14]; float raw[128]; float res48[96]; float cf = 0; int status = 0; };
struct FrRunStats { long pairs = 0, nokernel = 0, tg = 0, enc = 0, ran = 0, b_q = 0, b_tn = 0, b_pso = 0, b_slot = 0, b_tpg = 0; int bins = 0; double gpu_ms = 0, host_ms = 0, wait_ms = 0, cb_start = 0, cb_end = 0; };
#include <mach/mach_time.h>
static inline double gap_now() { static mach_timebase_info_data_t tb = {0, 0}; if (!tb.denom) mach_timebase_info(&tb); return (double)mach_absolute_time() * (double)tb.numer / (double)tb.denom * 1e-9; }
static const bool g_cblog = getenv("FJN_CBLOG") && atoi(getenv("FJN_CBLOG")) != 0;
static bool g_cofit_bypass = false;
static bool g_gap_in_replay = false; static int g_gap_stepno = 1;
#include "cofit_lane_h.hpp"
struct FjnPending {
  struct Disp { id<MTLComputePipelineState> ps; int Q, Tn; std::vector<int> mem; size_t tpg, ntg, tgb, per, ebo, plo; id<MTLBuffer> pl; };
  bool valid = false; id<MTLCommandBuffer> cb = nil; id<MTLBuffer> bOut = nil; std::vector<Disp> ds; std::vector<id<MTLBuffer>> keep;
  std::chrono::steady_clock::time_point tw; std::vector<FrOut>* out = nullptr; FrRunStats* S = nullptr; std::shared_ptr<CoJob> cojob;
  int stepno = 0; double commit_t = 0; size_t npg = 0, ncpu = 0; int nbins = 0;
  std::function<id<MTLCommandBuffer>()> rerun;
};
static const int g_wb_par_on = getenv("WB_PAR") ? atoi(getenv("WB_PAR")) : 0;
static std::function<void(int, const std::function<void(int)>&)> g_wb_par;
static const bool g_cfh = getenv("FJN_CISFUSE_HOST") != nullptr;
static const int g_cfchk = (g_cfh && getenv("FJN_CISFUSE_CHECK")) ? std::max(0, atoi(getenv("FJN_CISFUSE_CHECK"))) : 0;
static const size_t R48S = g_cfh ? 97 : 96;
static std::vector<float> g_cf_data; static id<MTLBuffer> g_cf_buf = nil; static bool g_cf_on = false;
static const bool g_res48b = getenv("FJN_RES48_HOST") != nullptr || g_cfh;
static const int g_hc1 = getenv("FJN_HC1") ? atoi(getenv("FJN_HC1")) : 0;
static const int g_aghc = getenv("FJN_AGHC") ? atoi(getenv("FJN_AGHC")) : 0;
static const size_t g_aghc_div = getenv("AGHC_GPUDIV") ? (size_t)std::max(1, atoi(getenv("AGHC_GPUDIV"))) : 1;
static thread_local bool tl_fjn_worker = false;
static std::atomic<unsigned> g_fjn_lib_gen{0};
static void (*g_fjn_notify)() = nullptr;
static void fjn_unpack(id<MTLCommandBuffer> cb, id<MTLBuffer> bOut, const std::vector<FjnPending::Disp>& ds, std::vector<FrOut>& out, FrRunStats& S,
                       const std::function<id<MTLCommandBuffer>()>* rerun = nullptr) {
  using clk_ = std::chrono::steady_clock;
  auto msnow = [](clk_::time_point a) { return std::chrono::duration<double, std::milli>(clk_::now() - a).count(); };
  static const bool rr_test_ = getenv("FJN_RERUN_TEST") && atoi(getenv("FJN_RERUN_TEST")) == 1;
  if (rr_test_ && rerun && *rerun && cb.status == MTLCommandBufferStatusCompleted) { const double rt0_ = gap_now(); cb = (*rerun)(); fprintf(stderr, "RERUN_TEST_DONE ms %.1f\n", (gap_now() - rt0_) * 1e3); }
  for (int a_ = 1; cb.status != MTLCommandBufferStatusCompleted && rerun && *rerun && a_ <= 3; a_++) {
    fprintf(stderr, "fjn_run command buffer failed: %s; re-run %d of the same dispatches\n", cb.error ? cb.error.localizedDescription.UTF8String : "?", a_);
    const double rr0_ = gap_now(); cb = (*rerun)();
    fprintf(stderr, "RERUN_DONE t_start %.4f t_end %.4f ms %.1f status %d\n", rr0_, gap_now(), (gap_now() - rr0_) * 1e3, (int)cb.status);
  }
  if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "fjn_run command buffer failed: %s\n", cb.error ? cb.error.description.UTF8String : "?"); gpu_fail(); }
  S.gpu_ms += (cb.GPUEndTime - cb.GPUStartTime) * 1e3; S.cb_start = cb.GPUStartTime; S.cb_end = cb.GPUEndTime;
  if (g_cofit_auto) { std::lock_guard<std::mutex> lk(g_auto_mx); g_auto_iv.emplace_back(cb.GPUStartTime, cb.GPUEndTime); }
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
  if (!tl_fjn_worker && g_wb_par_on && g_wb_par) {
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
static const int g_hc2 = getenv("FJN_HC2") ? atoi(getenv("FJN_HC2")) : 0;
static const int g_scrh = getenv("FJN_SCR_HOST") ? atoi(getenv("FJN_SCR_HOST")) : 0;
static const bool g_bufpool = (getenv("FJN_BUFPOOL") && atoi(getenv("FJN_BUFPOOL")) != 0) || (g_hc2 & 1);
static const size_t g_bufpool_cap = (size_t)(getenv("FJN_HC2_POOLMB") ? std::max(0, atoi(getenv("FJN_HC2_POOLMB"))) : 256) << 20;
struct BufPool { std::mutex m; size_t bytes = 0; std::map<std::pair<unsigned long, size_t>, std::vector<id<MTLBuffer>>> fr; };
static BufPool& bufpool() { static BufPool* p = new BufPool(); return *p; }
static std::atomic<unsigned long long> g_bpc[4];
static id<MTLBuffer> bp_get(id<MTLDevice> dev, size_t bytes, MTLResourceOptions opt, const void* src, bool zero_rest) {
  { const int pv = (opt & MTLResourceStorageModePrivate) ? 1 : 0; g_bpc[pv] += bytes; if (!g_bufpool) g_bpc[2 + pv] += bytes; }
  if (!g_bufpool) return src ? [dev newBufferWithBytes:src length:bytes options:opt] : [dev newBufferWithLength:bytes options:opt];
  size_t cap = 4096; while (cap < bytes) cap <<= 1; id<MTLBuffer> b = nil;
  { BufPool& P = bufpool(); std::lock_guard<std::mutex> lk(P.m); auto& v = P.fr[{(unsigned long)(opt & 0xFF), cap}]; if (!v.empty()) { b = v.back(); v.pop_back(); P.bytes -= cap; } }
  if (!b) { b = [dev newBufferWithLength:cap options:opt]; g_bpc[2 + ((opt & MTLResourceStorageModePrivate) ? 1 : 0)] += cap; }
  if (src && (opt & MTLResourceStorageModePrivate) == 0) { std::memcpy(b.contents, src, bytes); if (zero_rest && cap > bytes) std::memset((char*)b.contents + bytes, 0, cap - bytes); }
  return b;
}
static void bp_put(const std::vector<id<MTLBuffer>>& bs) {
  if (!g_bufpool) return; BufPool& P = bufpool(); std::lock_guard<std::mutex> lk(P.m);
  for (id<MTLBuffer> b : bs) if (b) { if ((g_hc2 & 1) && P.bytes + (size_t)b.length > g_bufpool_cap) continue;
    P.bytes += (size_t)b.length; P.fr[{(unsigned long)(b.resourceOptions & 0xFF), (size_t)b.length}].push_back(b); }
}
struct GapPair { int piece = 0; FrPairIn p; std::vector<uint32_t> px, py; std::vector<uint16_t> spx, spy; uint16_t sgx[36], sgy[36]; float thx[32], thy[32]; };
static const bool g_gap_collect_on = getenv("FJN_GAP_REPLAY") != nullptr;
static int g_gap_piece = 0; static std::vector<std::unique_ptr<GapPair>> g_gap_pairs; static float g_gap_yv[96]; static float g_gap_flr = 0;
static void gap_collect(const FrPairIn& q) {
  auto g = std::make_unique<GapPair>(); g->p = q; g->piece = g_gap_piece;
  g->px.assign(q.x.prog, q.x.prog + q.x.nprog); g->py.assign(q.y.prog, q.y.prog + q.y.nprog); g->spx.assign(q.x.sp, q.x.sp + q.x.nsp); g->spy.assign(q.y.sp, q.y.sp + q.y.nsp);
  std::memcpy(g->sgx, q.x.segs, 72); std::memcpy(g->sgy, q.y.segs, 72); std::memset(g->thx, 0, sizeof g->thx); std::memset(g->thy, 0, sizeof g->thy);
  if (q.thx) std::memcpy(g->thx, q.thx, (size_t)q.x.nth * 4); if (q.thy) std::memcpy(g->thy, q.thy, (size_t)q.y.nth * 4);
  g->p.x.prog = g->px.data(); g->p.y.prog = g->py.data(); g->p.x.sp = g->spx.data(); g->p.y.sp = g->spy.data(); g->p.x.segs = g->sgx; g->p.y.segs = g->sgy; g->p.thx = q.thx ? g->thx : nullptr; g->p.thy = q.thy ? g->thy : nullptr;
  g_gap_pairs.push_back(std::move(g)); }
struct FjnRecP { fjn::Rec r; std::vector<uint16_t> sp; bool ok = false; std::vector<uint32_t> kprog; std::vector<uint16_t> kseg; int kq = 0, knth = 0, knS = 0, knslot = 0, kcd = 0; };
struct FjnScratch {
  std::vector<FjnRecP> v; size_t n = 0;
  std::vector<int32_t> tab; std::vector<uint64_t> tkey; size_t mask = 0;
  std::vector<uint32_t> tokP, tokD, segv, hdrv, prv; std::vector<uint16_t> spv; std::vector<float> th0v;
};
static inline uint64_t fjn_hbytes(const void* p, size_t n, uint64_t h) {
  const uint8_t* b = (const uint8_t*)p; const uint64_t K = 0x9E3779B97F4A7C15ull, K2 = 0xBF58476D1CE4E5B9ull;
  uint64_t a0 = h ^ K, a1 = h + K2, a2 = ~h, a3 = h * K + 1; size_t i = 0;
  auto st = [&](uint64_t& a, uint64_t w) { a = (a ^ w) * K2; a ^= a >> 31; };
  for (; i + 32 <= n; i += 32) { uint64_t w0, w1, w2, w3; std::memcpy(&w0, b + i, 8); std::memcpy(&w1, b + i + 8, 8); std::memcpy(&w2, b + i + 16, 8); std::memcpy(&w3, b + i + 24, 8); st(a0, w0); st(a1, w1); st(a2, w2); st(a3, w3); }
  if (i < n) { uint8_t t[32] = {0}; std::memcpy(t, b + i, n - i); uint64_t w0, w1, w2, w3; std::memcpy(&w0, t, 8); std::memcpy(&w1, t + 8, 8); std::memcpy(&w2, t + 16, 8); std::memcpy(&w3, t + 24, 8); st(a0, w0); st(a1, w1); st(a2, w2); st(a3, w3); }
  uint64_t r = (a0 * K) ^ (a1 * K2) ^ ((a2 >> 17) * K) ^ (a3 * 0xD6E8FEB86659FD93ull) ^ (uint64_t)n; r ^= r >> 32; return r * K2; }
static const bool g_recc2 = getenv("RECC2") && atoi(getenv("RECC2")) == 1; static std::mutex g_rc2_m;
static long g_rc2_rec = 0, g_rc2_pairs = 0, g_rc2_calls = 0; static std::unordered_set<uint64_t> g_rc2_U, g_rc2_UP, g_rc2_prev, g_rc2_all;
static const bool g_f8h = getenv("FJN_F8H") && atoi(getenv("FJN_F8H")) == 1;
static inline uint32_t f8h_w(uint32_t v, uint32_t cd, uint32_t k) { return (v >> 31) ? v : (((v >> 30) & 1u) ? v + k : (((v >> 29) & 1u) ? v : v * cd)); }
static FjnScratch& fjn_tl_scratch() { static thread_local FjnScratch s; return s; }
static void fjn_run(id<MTLDevice> dev, id<MTLCommandQueue> cq, const std::function<id<MTLComputePipelineState>(const std::string&)>& pipe,
                    const std::function<id<MTLCommandBuffer>()>& mkcb, const std::vector<FrPairIn>& pairs, const float* yv96, float flr, bool closed, int maxit,
                    const fjn::WarpMat& WM, std::vector<FrOut>& out, FrRunStats& S, size_t tpgcap = 8, FjnPending* pend = nullptr) {
  using clk_ = std::chrono::steady_clock;
  auto msnow = [](clk_::time_point a) { return std::chrono::duration<double, std::milli>(clk_::now() - a).count(); };
  TmScope tm_fr_("fjn: cofit_start+setup");
  const auto th0 = clk_::now();
  const size_t NPR = pairs.size(); out.assign(NPR, FrOut()); S.pairs += (long)NPR;
  if (NPR == 0) return;
  std::vector<uint8_t> cosel; std::shared_ptr<CoJob> cojob;
  if (g_cofit_pm > 0 && !g_cofit_bypass) { cosel.assign(NPR, 0); cojob = cofit_start(pairs, cosel, yv96, flr, closed, maxit, WM, out, FJN_TMAX); }
  if (g_cofit_auto) { long nc_ = 0; if (cojob) for (uint8_t v_ : cosel) nc_ += v_ != 0; g_auto_np += (long)NPR - nc_; }
  tm_mark("fjn: rec table");
  FjnScratch lsc; FjnScratch& SC = (g_hc2 & 2) ? fjn_tl_scratch() : lsc; const bool flat = (g_hc2 & 2) != 0;
  std::vector<FjnRecP>& recs = SC.v; SC.n = 0; if (recs.capacity() < NPR * 2) recs.reserve(NPR * 2);
  std::unordered_multimap<uint64_t, int> rmap;
  if (flat) { size_t cap = 16; while (cap < NPR * 4 + 16) cap <<= 1; SC.tab.assign(cap, -1); if (SC.tkey.size() < cap) SC.tkey.resize(cap); SC.mask = cap - 1; } else rmap.reserve(NPR * 3);
  auto fresh = [&]() -> FjnRecP& { if (SC.n == recs.size()) recs.emplace_back(); FjnRecP& q = recs[SC.n++]; q.ok = false; q.r.perP = q.r.perD = 0; return q; };
  auto same = [&](const FjnRecP& c, const FrTreeRef& tr, int CD) -> bool {
    return c.kq == tr.Q && c.knth == tr.nth && c.knS == tr.nS && c.knslot == tr.nslot && c.kcd == CD && c.kprog.size() == tr.nprog && c.sp.size() == tr.nsp &&
           (tr.nprog == 0 || std::memcmp(c.kprog.data(), tr.prog, tr.nprog * 4) == 0) && std::memcmp(c.kseg.data(), tr.segs, 36 * 2) == 0 && (tr.nsp == 0 || std::memcmp(c.sp.data(), tr.sp, tr.nsp * 2) == 0); };
  auto getrec = [&](const FrTreeRef& tr, int CD) -> int {
    uint64_t key = 1469598103934665603ull;
    auto mix = [&](uint64_t v) { key = (key ^ v) * 1099511628211ull; };
    if (flat && (g_hc2 & 8)) {
      key = fjn_hbytes(tr.prog, tr.nprog * 4, key); key = fjn_hbytes(tr.segs, 72, key); key = fjn_hbytes(tr.sp, tr.nsp * 2, key);
      mix((uint64_t)tr.nprog); mix((uint64_t)CD); mix((uint64_t)tr.Q); mix((uint64_t)tr.nth); mix((uint64_t)tr.nS); mix((uint64_t)tr.nslot); mix((uint64_t)tr.nsp);
    } else {
    for (size_t i = 0; i < tr.nprog; i++) mix((uint64_t)tr.prog[i]); mix((uint64_t)tr.nprog); mix((uint64_t)CD); mix((uint64_t)tr.Q); mix((uint64_t)tr.nth); mix((uint64_t)tr.nS); mix((uint64_t)tr.nslot); mix((uint64_t)tr.nsp);
    for (int i = 0; i < 36; i++) mix((uint64_t)tr.segs[i]);
    for (size_t i = 0; i < tr.nsp; i++) mix((uint64_t)tr.sp[i]); }
    size_t slot = 0;
    if (flat) { for (slot = (size_t)(key ^ (key >> 29)) & SC.mask; SC.tab[slot] >= 0; slot = (slot + 1) & SC.mask) if (SC.tkey[slot] == key && same(recs[(size_t)SC.tab[slot]], tr, CD)) return SC.tab[slot]; }
    else { auto rg = rmap.equal_range(key); for (auto it = rg.first; it != rg.second; ++it) if (same(recs[(size_t)it->second], tr, CD)) return it->second; }
    FjnRecP& rp = fresh(); fjn::Rec& r = rp.r; rp.kq = tr.Q; rp.knth = tr.nth; rp.knS = tr.nS; rp.knslot = tr.nslot; rp.kcd = CD; rp.kprog.assign(tr.prog, tr.prog + tr.nprog); rp.kseg.assign(tr.segs, tr.segs + 36); r.nth = tr.nth; r.Q = tr.Q; r.nS = tr.nS;
    const int id = (int)SC.n - 1;
    auto reg = [&]() { if (flat) { SC.tab[slot] = id; SC.tkey[slot] = key; } else rmap.emplace(key, id); };
    if (g_kgp) { rp.sp.assign(tr.sp, tr.sp + tr.nsp); rp.ok = false; reg(); return id; }
    struct TT { int nth, nS, nslot; } tt{tr.nth, tr.nS, tr.nslot};
    uint32_t SBp = 0, SBd = 0; const std::pair<const uint32_t*, size_t> pg{tr.prog, tr.nprog};
    static const bool g_epair = getenv("FJN_EPAIR") && atoi(getenv("FJN_EPAIR")) == 1;
    bool ok = g_epair ? fr_enc_pair(tr.nth, tr.nS, tr.nslot, CD, tr.prog, tr.nprog, tr.segs, r.tokP, r.tokD, r.sgP, r.sgD, r.hP[3], r.hD[3], SBp, SBd)
                      : fjn::enc_tree(tt, 48u, pg, tr.segs, r.tokP, r.sgP, r.hP[3], SBp) && fjn::enc_tree(tt, 48u * (uint32_t)CD, pg, tr.segs, r.tokD, r.sgD, r.hD[3], SBd);
    if (ok) { r.hP[1] = r.hD[1] = (uint32_t)tr.nS | ((uint32_t)tr.Q << 16); r.hP[0] = r.hP[2] = r.hD[0] = r.hD[2] = 0;
      if (fjn::g_wd & 1) { fjn::wd_s3_rewrite(r.tokP, r.sgP, SBp); fjn::wd_s3_rewrite(r.tokD, r.sgD, SBd); }
      r.perP = SBp + fjn::wd_k() * (uint32_t)(tr.nS + 2) + fjn::wd_pad(); r.perD = SBd + fjn::wd_k() * (uint32_t)CD * (uint32_t)(tr.nS + 2) + fjn::wd_pad(); }
    rp.sp.assign(tr.sp, tr.sp + tr.nsp); rp.ok = ok; reg(); return id; };
  auto encrec = [&](size_t id) {
    FjnRecP& rp = recs[id]; fjn::Rec& r = rp.r; const int CD = rp.kcd;
    struct TT { int nth, nS, nslot; } tt{rp.knth, rp.knS, rp.knslot};
    uint32_t SBp = 0, SBd = 0; const std::pair<const uint32_t*, size_t> pg{rp.kprog.data(), rp.kprog.size()};
    bool ok;
    if (g_f8h) { ok = fjn::enc_tree(tt, 48u, pg, rp.kseg.data(), r.tokP, r.sgP, r.hP[3], SBp);
      static const bool g_f8c = getenv("FJN_F8C") && atoi(getenv("FJN_F8C")) == 1;
      if (ok && g_f8c) { const size_t n = r.tokP.size() / 4; std::vector<char> bnd(n + 1, 0); for (int q = 0; q <= rp.kq; q++) if (r.sgP[q] <= n) bnd[r.sgP[q]] = 1;
        std::vector<uint32_t> o; o.reserve(r.tokP.size()); std::vector<uint32_t> ni(n + 1, 0); size_t k = 0; const uint32_t* t4 = r.tokP.data();
        while (k < n) { ni[k] = (uint32_t)(o.size() / 4); const uint32_t* a = t4 + k * 4;
          if (k + 1 < n && !bnd[k + 1] && a[0] == 0u && a[1] == 0xFFFFFFFFu) { const uint32_t* b = t4 + (k + 1) * 4;
            if (b[0] >= 7u && b[0] <= 15u && (b[2] >> 31)) { o.push_back(b[0]); o.push_back(b[1]); o.push_back(a[2]); o.push_back(a[3]); ni[k + 1] = ni[k]; k += 2; continue; } }
          o.insert(o.end(), a, a + 4); k++; }
        ni[n] = (uint32_t)(o.size() / 4); for (int q = 0; q <= rp.kq; q++) if (r.sgP[q] <= n) r.sgP[q] = ni[r.sgP[q]]; r.tokP.swap(o); }
      if (ok) { const uint32_t cd = (uint32_t)CD; SBd = SBp * cd; ok = SBd < 65536u;
        if (ok) { const uint32_t k8 = SBd - SBp; const size_t nw = r.tokP.size(); r.tokD.resize(nw); const uint32_t* s_ = r.tokP.data(); uint32_t* d_ = r.tokD.data();
          for (size_t i = 0; i < nw; i += 4) { d_[i] = s_[i]; d_[i + 1] = s_[i + 1] == 0xFFFFFFFFu ? s_[i + 1] : s_[i + 1] * cd; d_[i + 2] = f8h_w(s_[i + 2], cd, k8); d_[i + 3] = f8h_w(s_[i + 3], cd, k8); }
          for (int q = 0; q < 18; q++) r.sgD[q] = r.sgP[q]; for (int q = 18; q < 36; q++) r.sgD[q] = f8h_w(r.sgP[q], cd, k8);
          r.hD[3] = ((uint32_t)rp.knS + 2u) | (SBd << 8) | ((uint32_t)rp.knth << 24); } } }
    else ok = fjn::enc_tree(tt, 48u, pg, rp.kseg.data(), r.tokP, r.sgP, r.hP[3], SBp) && fjn::enc_tree(tt, 48u * (uint32_t)CD, pg, rp.kseg.data(), r.tokD, r.sgD, r.hD[3], SBd);
    if (ok) { r.hP[1] = r.hD[1] = (uint32_t)rp.knS | ((uint32_t)rp.kq << 16); r.hP[0] = r.hP[2] = r.hD[0] = r.hD[2] = 0;
      if (fjn::g_wd & 1) { fjn::wd_s3_rewrite(r.tokP, r.sgP, SBp); fjn::wd_s3_rewrite(r.tokD, r.sgD, SBd); }
      r.perP = SBp + fjn::wd_k() * (uint32_t)(rp.knS + 2) + fjn::wd_pad(); r.perD = SBd + fjn::wd_k() * (uint32_t)CD * (uint32_t)(rp.knS + 2) + fjn::wd_pad(); }
    rp.ok = ok; };
  static const bool g_pq = !(getenv("FJN_PQ") && atoi(getenv("FJN_PQ")) == 0);
  const size_t v5_lpnk = closed ? 10 : 9;
  auto v5_lpneed = [&](const FrPairIn& p) -> size_t { if (!(fjn::g_v5 & 1)) return 0; const size_t n = (size_t)p.x.nth + (size_t)p.y.nth + v5_lpnk; return n > 16 ? n * (n + 1) / 2 : 0; };
  std::map<std::tuple<int, int, int>, std::vector<int>> bins;
  std::vector<int> recx(NPR, -1), recy(NPR, -1);
  std::vector<int> gA, gB;
  if (g_kgp) { gA.resize(NPR); gB.resize(NPR);
    for (size_t k = 0; k < NPR; k++) { if (!cosel.empty() && cosel[k]) { gA[k] = gB[k] = -1; continue; } const FrPairIn& p = pairs[k]; const int CD = fjn::geom_CD(std::max(p.x.nth, p.y.nth)); gA[k] = getrec(p.x, CD); gB[k] = getrec(p.y, CD); }
    tm_mark("fjn: encrec (pool)");
    if (g_recc2) { std::lock_guard<std::mutex> lk_(g_rc2_m); g_rc2_calls++; g_rc2_rec += (long)SC.n; for (size_t k = 0; k < NPR; k++) if (gA[k] >= 0) g_rc2_pairs++;
      for (size_t i = 0; i < SC.n; i++) { const FjnRecP& q = recs[i]; uint64_t h = fjn_hbytes(q.kprog.data(), q.kprog.size() * 4, 7); h = fjn_hbytes(q.kseg.data(), q.kseg.size() * 2, h); h = fjn_hbytes(q.sp.data(), q.sp.size() * 2, h);
        const uint64_t m[4] = {(uint64_t)q.kq, (uint64_t)q.knth, (uint64_t)q.knS, (uint64_t)q.knslot}; h = fjn_hbytes(m, sizeof m, h); g_rc2_UP.insert(h);
        const uint64_t c = (uint64_t)q.kcd; g_rc2_U.insert(fjn_hbytes(&c, 8, h)); } }
    if (!tl_fjn_worker && g_kgp_par) g_kgp_par((int)SC.n, [&](int i) { encrec((size_t)i); });
    else { const size_t NR = SC.n; const int nw = (g_fcpw > 1 && !tl_fjn_worker) ? std::min<int>(g_fcpw, (int)std::max<size_t>(1, NR / 64)) : 1;
      if (nw <= 1) { for (size_t i = 0; i < NR; i++) encrec(i); }
      else { std::vector<std::thread> th; for (int w = 0; w < nw; w++) th.emplace_back([&, w] { for (size_t i = (size_t)w; i < NR; i += (size_t)nw) encrec(i); }); for (auto& t : th) t.join(); } }
  }
  tm_mark("fjn: bin assign");
  for (size_t k = 0; k < NPR; k++) {
    if (!cosel.empty() && cosel[k]) continue;
    const FrPairIn& p = pairs[k]; const int Qc = std::max(p.x.Q, p.y.Q), Tc = std::max(p.x.nth, p.y.nth); const int CD = fjn::geom_CD(Tc);
    const int a = g_kgp ? gA[k] : getrec(p.x, CD), b = g_kgp ? gB[k] : getrec(p.y, CD);
    if (!recs[(size_t)a].ok || !recs[(size_t)b].ok) { S.enc++; continue; }
    size_t pr = 0; { const fjn::Rec& ra = recs[(size_t)a].r; const fjn::Rec& rb = recs[(size_t)b].r; pr = std::max(std::max(ra.perP, ra.perD), std::max(rb.perP, rb.perD)); }
    pr = std::max(pr, v5_lpneed(p)); pr = (pr + 3) / 4 * 4; int pq = 0;
    if (g_pq) { size_t v = 64; while (v < pr) v += v / 5 + 1; pq = (int)((v + 3) / 4 * 4); }
    recx[k] = a; recy[k] = b; bins[std::make_tuple(Qc, Tc, pq)].push_back((int)k);
  }
  tm_mark("fjn: token pack");
  size_t ntokP = 0, ntokD = 0, nsp = 0;
  for (size_t r = 0; r < SC.n; r++) { const FjnRecP& rp = recs[r]; if (rp.ok) { ntokP += rp.r.tokP.size() / 4; ntokD += rp.r.tokD.size() / 4; nsp += rp.sp.size(); } }
  FjnScratch lsc2; FjnScratch& PK = flat ? SC : lsc2;
  std::vector<uint32_t>& tokP = PK.tokP; std::vector<uint32_t>& tokD = PK.tokD; std::vector<uint32_t>& segv = PK.segv; std::vector<uint32_t>& hdrv = PK.hdrv; std::vector<uint16_t>& spv = PK.spv;
  tokP.assign(ntokP * 4 + 4, 0u); tokD.assign(ntokD * 4 + 4, 0u); segv.assign(SC.n * 72 + 72, 0u); hdrv.assign(SC.n * 8 + 8, 0u); spv.assign(nsp + 4, (uint16_t)0);
  { size_t oP = 0, oD = 0, oS = 0;
    for (size_t r = 0; r < SC.n; r++) { FjnRecP& rp = recs[r]; if (!rp.ok) continue; fjn::Rec& q = rp.r;
      std::memcpy(&tokP[oP * 4], q.tokP.data(), q.tokP.size() * 4); std::memcpy(&tokD[oD * 4], q.tokD.data(), q.tokD.size() * 4); if (!rp.sp.empty()) std::memcpy(&spv[oS], rp.sp.data(), rp.sp.size() * 2);
      uint32_t* h = &hdrv[r * 8]; h[0] = (uint32_t)oP; h[1] = q.hP[1]; h[2] = (uint32_t)oS; h[3] = q.hP[3]; h[4] = (uint32_t)oD; h[5] = q.hD[1]; h[6] = (uint32_t)oS; h[7] = q.hD[3];
      std::memcpy(&segv[r * 72], q.sgP, 36 * 4); std::memcpy(&segv[r * 72 + 36], q.sgD, 36 * 4);
      oP += q.tokP.size() / 4; oD += q.tokD.size() / 4; oS += rp.sp.size(); } }
  std::vector<uint32_t>& prv = PK.prv; prv.assign(NPR * 2, 0u); std::vector<float>& th0v = PK.th0v; th0v.assign(NPR * 64, 0.0f);
  for (size_t k = 0; k < NPR; k++) { if (recx[k] < 0) continue; prv[k * 2] = (uint32_t)recx[k]; prv[k * 2 + 1] = (uint32_t)recy[k];
    if (pairs[k].thx) std::memcpy(&th0v[k * 64], pairs[k].thx, (size_t)pairs[k].x.nth * 4); if (pairs[k].thy) std::memcpy(&th0v[k * 64 + 32], pairs[k].thy, (size_t)pairs[k].y.nth * 4); }
  std::vector<float> Smv(48 * 9 + 9); std::memcpy(Smv.data(), WM.S, 48 * 9 * 4); std::memcpy(Smv.data() + 48 * 9, WM.u0, 9 * 4);
  static const bool fjnUntracked = getenv("FJN_UNTRACKED") != nullptr && getenv("FJN_UNTRACKED")[0] == '1';
  const MTLResourceOptions inOpt = MTLResourceStorageModeShared | (fjnUntracked ? MTLResourceHazardTrackingModeUntracked : MTLResourceHazardTrackingModeTracked);
  tm_mark("fjn: newBuffer x8");
  id<MTLBuffer> bP = bp_get(dev, tokP.size() * 4, inOpt, tokP.data(), false);
  id<MTLBuffer> bD = bp_get(dev, tokD.size() * 4, inOpt, tokD.data(), false);
  id<MTLBuffer> bH = bp_get(dev, hdrv.size() * 4, inOpt, hdrv.data(), false);
  id<MTLBuffer> bSp = bp_get(dev, spv.size() * 2, inOpt, spv.data(), false);
  id<MTLBuffer> bSg = bp_get(dev, segv.size() * 4, inOpt, segv.data(), false);
  id<MTLBuffer> bPr = bp_get(dev, prv.size() * 4, inOpt, prv.data(), false);
  id<MTLBuffer> bTh = bp_get(dev, th0v.size() * 4, inOpt, th0v.data(), false);
  id<MTLBuffer> bOut = bp_get(dev, NPR * (128 + (g_res48b ? R48S : 0)) * 4, MTLResourceStorageModeShared, nullptr, false);
  if (g_bufpool) std::memset(bOut.contents, 0, bOut.length);
  if (!g_bufpool) for (id<MTLBuffer> b_ : {bP, bD, bH, bSp, bSg, bPr, bTh, bOut}) memdiag::metal_note(b_, (size_t)b_.length, false);
  std::memset(bOut.contents, 0xFF, NPR * 128 * 4);
  struct Disp { id<MTLComputePipelineState> ps; int Q, Tn; std::vector<int> mem; size_t tpg, ntg, tgb, per, ebo, plo; id<MTLBuffer> pl; size_t scro = 0; };
  tm_mark("fjn: bin plan");
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
      if (slot <= 32768 || g_scrh) work.push_back({std::make_tuple(Q, Tn, 0), g1.second, true});
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
    id<MTLComputePipelineState> ps = nil;
    if (g_hc2 & 4) { struct PT { id<MTLComputePipelineState> p[17][FJN_TMAX + 1]; unsigned gen = 0; }; static thread_local PT pso_tab;
      if (pso_tab.gen != g_fjn_lib_gen.load()) { for (auto& row : pso_tab.p) for (auto& e : row) e = nil; pso_tab.gen = g_fjn_lib_gen.load(); }
      ps = pso_tab.p[Q][Tn]; if (!ps) { ps = pipe("fjn_q" + std::to_string(Q) + "_t" + std::to_string(Tn)); pso_tab.p[Q][Tn] = ps; } }
    else ps = pipe("fjn_q" + std::to_string(Q) + "_t" + std::to_string(Tn));
    if (!ps) { S.nokernel += (long)mem.size(); S.b_pso += (long)mem.size(); continue; }
    size_t per = 0; for (int k : mem) { const fjn::Rec& a = recs[(size_t)recx[(size_t)k]].r; const fjn::Rec& b = recs[(size_t)recy[(size_t)k]].r;
      per = std::max<size_t>(per, std::max(std::max(a.perP, a.perD), std::max(b.perP, b.perD))); per = std::max<size_t>(per, v5_lpneed(pairs[(size_t)k])); }
    per = (per + 3) / 4 * 4; if (std::get<2>(work[wi].key) > 0 && !(fjn::g_v5 & 2)) per = (size_t)std::get<2>(work[wi].key);
    per = std::max<size_t>(per, (fjn::g_v5 & 1) ? (size_t)4 : (size_t)fjn::per_min(Tn)); per = (per + 3) / 4 * 4;
#ifdef H308_PERPAD
    per += (size_t)H308_PERPAD;
#endif
    const size_t slot = ((size_t)fjn::slotc_floats(Q, Tn) + (fjn::g_lay == 2 ? std::max<size_t>((size_t)fjn::ovl_floats(Q, Tn), 2 * per) : 2 * per)) * 4;
    if (slot > 32768 && !g_scrh) { S.tg += (long)mem.size(); S.b_slot += (long)mem.size(); continue; }
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
    Disp d; d.ps = ps; d.Q = Q; d.Tn = Tn; d.mem = mem; d.tpg = tpg; d.ntg = (mem.size() + tpg - 1) / tpg; d.tgb = (tpg * slot + 15) / 16 * 16 + (getenv("FJN_TGPAD") ? (size_t)atoi(getenv("FJN_TGPAD")) : 0);     d.per = per; d.ebo = ebTot; d.plo = plAll.size() * 4;
    ebTot += ((mem.size() + tpg) * (size_t)fjn::eb_stride(Q, Tn) * 4 + 255) / 256 * 256;
    if (getenv("FJN_VERBOSE")) printf("fjn bin Q%d T%d: pairs %zu per %zu slot %zu B tpg %zu ntg %zu maxTotal %lu execWidth %lu\n", Q, Tn, mem.size(), per, slot, tpg, d.ntg, (unsigned long)ps.maxTotalThreadsPerThreadgroup, (unsigned long)ps.threadExecutionWidth);
    plAll.insert(plAll.end(), mem.begin(), mem.end());
#ifdef H308_TRACE
    for (int k : mem) { const FrPairIn& pp = pairs[(size_t)k]; float f[6] = {(float)std::get<2>(work[wi].key), (float)per, (float)tpg, (float)pp.x.nslot, (float)pp.y.nslot, (float)mem.size()}; H308T('J', pp.k1, pp.k2, f, 6); }
#endif
    ds.push_back(std::move(d));
  }
  if (ds.empty()) { cofit_join(cojob); S.host_ms += msnow(th0); return; }
  if (g_gap_collect_on && !g_cofit_bypass) { g_gap_piece++; for (Disp& d : ds) for (int k : d.mem) gap_collect(pairs[(size_t)k]); }
  tm_mark("fjn: EB/pl buffers");
  id<MTLBuffer> bEB = bp_get(dev, std::max<size_t>(ebTot, 256), MTLResourceStorageModePrivate, nullptr, false);
  id<MTLBuffer> bPlAll = bp_get(dev, std::max<size_t>(plAll.size() * 4, 4), MTLResourceStorageModeShared, (plAll.empty() ? (const void*)"\0\0\0\0" : (const void*)plAll.data()), false);
  if (!g_bufpool) { memdiag::metal_note(bEB, (size_t)bEB.length, true); memdiag::metal_note(bPlAll, (size_t)bPlAll.length, false); }
  for (Disp& d : ds) d.pl = bPlAll;
  id<MTLBuffer> bScr = nil;
  if (g_scrh) { size_t so = 0; for (Disp& d : ds) { d.scro = so; so += d.mem.size() * d.tgb; } bScr = bp_get(dev, std::max<size_t>(so, 256), MTLResourceStorageModePrivate, nullptr, false);
    if (getenv("FJN_SCR_LOG")) { size_t np_ = 0; for (Disp& d : ds) np_ += d.mem.size(); printf("FJN_SCR call: bins %zu pairs %zu scratch %.2f MB (buffer %.2f MB)\n", ds.size(), np_, so / 1048576.0, (double)bScr.length / 1048576.0); } }
  S.host_ms += msnow(th0);
  const std::vector<float> yvc(yv96, yv96 + 96);
  struct EncCtx { std::vector<Disp> ds; id<MTLBuffer> bP, bD, bH, bSp, bSg, bPr, bTh, bOut, bEB, bPl, bYv, bSm, bScr; };
  std::shared_ptr<EncCtx> cx;
  if (g_hc2 & 4) {
    cx = std::make_shared<EncCtx>(); cx->ds = std::move(ds); cx->bP = bP; cx->bD = bD; cx->bH = bH; cx->bSp = bSp; cx->bSg = bSg; cx->bPr = bPr; cx->bTh = bTh; cx->bOut = bOut; cx->bEB = bEB; cx->bPl = bPlAll; cx->bScr = bScr;
    struct CC { id<MTLBuffer> yv = nil, sm = nil; std::vector<float> yvv, smv; }; static thread_local CC cc;
    if (!cc.yv || cc.yvv != yvc) { cc.yvv = yvc; cc.yv = [dev newBufferWithBytes:yvc.data() length:yvc.size() * 4 options:MTLResourceStorageModeShared]; }
    if (!cc.sm || cc.smv != Smv) { cc.smv = Smv; cc.sm = [dev newBufferWithBytes:Smv.data() length:Smv.size() * 4 options:MTLResourceStorageModeShared]; }
    cx->bYv = cc.yv; cx->bSm = cc.sm;
  }
  const std::vector<Disp>& dsv = cx ? cx->ds : ds;
  std::function<void(id<MTLCommandBuffer>, size_t, size_t)> enc;
  if (cx) enc = [cx, dev, NPR, closed, maxit, flr](id<MTLCommandBuffer> cb, size_t lo, size_t hi) {
    id<MTLComputeCommandEncoder> en = [cb computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
    if (getenv("TL_LABELS") && getenv("TL_LABELS")[0] == '1') en.label = [NSString stringWithFormat:@"fjn_run chunk Disp=%zu", cx->ds.size()];
    id<MTLBuffer> bs[9] = {cx->bP, cx->bD, cx->bH, cx->bSp, cx->bSg, nil, cx->bPr, cx->bTh, cx->bOut};
    for (int k = 0; k < 9; k++) if (bs[k]) [en setBuffer:bs[k] offset:0 atIndex:k];
    if (g_res48b) [en setBuffer:cx->bOut offset:(NSUInteger)(NPR * 128 * 4) atIndex:15];
    if (g_cfh) { if (!g_cf_buf) { std::vector<float> z(400, 0.0f); const std::vector<float>& s0 = (g_cf_data.size() == 400 || g_cf_data.size() == 496) ? g_cf_data : z; g_cf_buf = [dev newBufferWithBytes:s0.data() length:s0.size() * 4 options:MTLResourceStorageModeShared]; } [en setBuffer:g_cf_buf offset:0 atIndex:16]; }
    static uint32_t rpv[4] = {1, 1, 0, 1}; static float bwr = 0; static bool rpinit = false;
    if (!rpinit) { rpinit = true; if (const char* ev = getenv("FJN_REPS")) { unsigned a = 1, b = 1, c = 0, e2 = 1, f = 0; sscanf(ev, "%u,%u,%u,%u,%u", &a, &b, &c, &e2, &f); rpv[0] = a; rpv[1] = b; rpv[2] = c; rpv[3] = e2; bwr = (float)f; } }
    [en setBuffer:cx->bYv offset:0 atIndex:5];
    [en setBytes:rpv length:16 atIndex:14];
    [en setBuffer:cx->bSm offset:0 atIndex:11];
    [en setBuffer:cx->bEB offset:0 atIndex:10];
    [en setBuffer:cx->bPl offset:0 atIndex:12];
    for (size_t di = lo; di < hi; di++) { const Disp& d = cx->ds[di];
      [en setComputePipelineState:d.ps];
      float tf[4] = {(float)d.Q * 9.0f * ldexpf(1.0f, -24), flr, bwr, 0.0f}; [en setBytes:tf length:16 atIndex:9];
      [en setBufferOffset:d.ebo atIndex:10];
      [en setBufferOffset:d.plo atIndex:12];
      uint32_t wq[4] = {(uint32_t)d.mem.size(), (uint32_t)d.per, closed ? 1u : 0u, (uint32_t)maxit}; [en setBytes:wq length:16 atIndex:13];
      if (cx->bScr) { [en setBuffer:cx->bScr offset:d.scro atIndex:17]; uint32_t sq[4] = {(uint32_t)(d.tgb / 4), (uint32_t)(d.tgb / 4), (uint32_t)d.Q, (uint32_t)d.Tn}; [en setBytes:sq length:16 atIndex:18]; }
      else [en setThreadgroupMemoryLength:d.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g_aghc_div > 1 ? std::max<size_t>(1, d.ntg / g_aghc_div) : d.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * d.tpg, 1, 1)];
    }
    [en endEncoding];
  };
  else enc = [=](id<MTLCommandBuffer> cb, size_t lo, size_t hi) {
  id<MTLComputeCommandEncoder> en = [cb computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
  if (getenv("TL_LABELS") && getenv("TL_LABELS")[0] == '1') en.label = [NSString stringWithFormat:@"fjn_run chunk Disp=%zu", ds.size()];
  for (size_t di = lo; di < hi; di++) { const Disp& d = ds[di];
    [en setComputePipelineState:d.ps];
    id<MTLBuffer> bs[9] = {bP, bD, bH, bSp, bSg, nil, bPr, bTh, bOut};
    for (int k = 0; k < 9; k++) if (bs[k]) [en setBuffer:bs[k] offset:0 atIndex:k];
    if (g_res48b) [en setBuffer:bOut offset:(NSUInteger)(NPR * 128 * 4) atIndex:15];
    if (g_cfh) { if (!g_cf_buf) { std::vector<float> z(400, 0.0f); const std::vector<float>& s0 = (g_cf_data.size() == 400 || g_cf_data.size() == 496) ? g_cf_data : z; g_cf_buf = [dev newBufferWithBytes:s0.data() length:s0.size() * 4 options:MTLResourceStorageModeShared]; } [en setBuffer:g_cf_buf offset:0 atIndex:16]; }
    [en setBytes:yvc.data() length:96 * 4 atIndex:5];
    static uint32_t rpv[4] = {1, 1, 0, 1}; static float bwr = 0; static bool rpinit = false;
    if (!rpinit) { rpinit = true; if (const char* ev = getenv("FJN_REPS")) { unsigned a = 1, b = 1, c = 0, e2 = 1, f = 0; sscanf(ev, "%u,%u,%u,%u,%u", &a, &b, &c, &e2, &f); rpv[0] = a; rpv[1] = b; rpv[2] = c; rpv[3] = e2; bwr = (float)f; } }
    float tf[4] = {(float)d.Q * 9.0f * ldexpf(1.0f, -24), flr, bwr, 0.0f}; [en setBytes:tf length:16 atIndex:9];
    [en setBytes:rpv length:16 atIndex:14];
    [en setBuffer:bEB offset:d.ebo atIndex:10];
    [en setBytes:Smv.data() length:Smv.size() * 4 atIndex:11];
    [en setBuffer:d.pl offset:d.plo atIndex:12];
    uint32_t wq[4] = {(uint32_t)d.mem.size(), (uint32_t)d.per, closed ? 1u : 0u, (uint32_t)maxit}; [en setBytes:wq length:16 atIndex:13];
    if (bScr) { [en setBuffer:bScr offset:d.scro atIndex:17]; uint32_t sq[4] = {(uint32_t)(d.tgb / 4), (uint32_t)(d.tgb / 4), (uint32_t)d.Q, (uint32_t)d.Tn}; [en setBytes:sq length:16 atIndex:18]; }
      else [en setThreadgroupMemoryLength:d.tgb atIndex:0];
    [en dispatchThreadgroups:MTLSizeMake(g_aghc_div > 1 ? std::max<size_t>(1, d.ntg / g_aghc_div) : d.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * d.tpg, 1, 1)];
  }
  [en endEncoding];
  };
  S.bins += (int)dsv.size();
  const std::function<id<MTLCommandBuffer>()> mk = mkcb;
  static const size_t rrg_ = getenv("FJN_RERUN_GRP") ? (size_t)std::max(1, atoi(getenv("FJN_RERUN_GRP"))) : 1;
  std::function<id<MTLCommandBuffer>()> rerun = [enc, mk, nd = dsv.size()]() -> id<MTLCommandBuffer> {
    const size_t ng = (nd + rrg_ - 1) / rrg_;
    std::vector<id<MTLCommandBuffer>> cs(ng);
    for (size_t gi = 0; gi < ng; gi++) { cs[gi] = mk(); enc(cs[gi], gi * rrg_, std::min(nd, (gi + 1) * rrg_)); [cs[gi] commit]; }
    for (size_t gi = 0; gi < ng; gi++) {
      [cs[gi] waitUntilCompleted];
      for (int a_ = 0; a_ < 3 && cs[gi].status != MTLCommandBufferStatusCompleted; a_++) { cs[gi] = mk(); enc(cs[gi], gi * rrg_, std::min(nd, (gi + 1) * rrg_)); [cs[gi] commit]; [cs[gi] waitUntilCompleted]; }
      if (cs[gi].status != MTLCommandBufferStatusCompleted) { for (size_t gj = gi + 1; gj < ng; gj++) [cs[gj] waitUntilCompleted]; return cs[gi]; }
    }
    return ng ? cs[ng - 1] : mk(); };
  tm_mark("fjn: encode CB+commit");
  if (kcap_active() && !g_gap_in_replay) {
    static std::mutex kcap_m_; std::lock_guard<std::mutex> kcap_lk_(kcap_m_);
    std::vector<KcapDispIn> kd_; std::vector<std::string> kn_; kn_.reserve(dsv.size());
    for (const Disp& d : dsv) { kn_.push_back("fjn_q" + std::to_string(d.Q) + "_t" + std::to_string(d.Tn)); kd_.push_back({kn_.back().c_str(), d.Q, d.Tn, d.ebo, d.plo, d.tgb, d.ntg, d.tpg, d.mem.size(), d.per}); }
    const void* kb_[8] = {bP.contents, bD.contents, bH.contents, bSp.contents, bSg.contents, bPr.contents, bTh.contents, bPlAll.contents};
    const size_t kl_[8] = {(size_t)bP.length, (size_t)bD.length, (size_t)bH.length, (size_t)bSp.length, (size_t)bSg.length, (size_t)bPr.length, (size_t)bTh.length, (size_t)bPlAll.length};
    const void* cfp_ = nullptr; size_t cfn_ = 0; std::vector<float> cfz_;
    if (g_cfh) { if (g_cf_buf) { cfp_ = g_cf_buf.contents; cfn_ = (size_t)g_cf_buf.length; } else if (g_cf_data.size() == 400 || g_cf_data.size() == 496) { cfp_ = g_cf_data.data(); cfn_ = g_cf_data.size() * 4; } else { cfz_.assign(400, 0.0f); cfp_ = cfz_.data(); cfn_ = cfz_.size() * 4; } }
    kcap_write_call(kd_, kb_, kl_, cfp_, cfn_, (size_t)bEB.length, (size_t)bOut.length, NPR, closed, maxit, flr, yvc.data(), Smv.data(), g_res48b, g_cfh, R48S);
  }
  id<MTLCommandBuffer> cb = mkcb(); enc(cb, 0, dsv.size());
  if (pend && g_fjn_notify) [cb addCompletedHandler:^(id<MTLCommandBuffer> c_) { (void)c_; if (g_fjn_notify) g_fjn_notify(); }];
  const auto tw = clk_::now(); const double gap_ct = gap_now(); [cb commit];
  if (pend) { pend->cb = cb; pend->bOut = bOut; pend->ds.clear(); for (const Disp& d : dsv) pend->ds.push_back({d.ps, d.Q, d.Tn, d.mem, d.tpg, d.ntg, d.tgb, d.per, d.ebo, d.plo, d.pl}); pend->keep = {bP, bD, bH, bSp, bSg, bPr, bTh, bEB, bOut, bPlAll}; if (bScr) pend->keep.push_back(bScr); pend->tw = tw; pend->stepno = g_gap_stepno; pend->commit_t = gap_ct; pend->npg = 0; for (const Disp& d : dsv) pend->npg += d.mem.size(); pend->ncpu = cojob ? cojob->prs.size() : 0; pend->nbins = (int)dsv.size(); pend->out = &out; pend->S = &S; pend->cojob = cojob; pend->rerun = rerun; pend->valid = true; S.host_ms += msnow(th0); return; }
  tm_mark("fjn: sync CB wait");
  [cb waitUntilCompleted]; S.wait_ms += msnow(tw);
  tm_mark("fjn: sync cofit_join");
  cofit_join(cojob);
  std::vector<FjnPending::Disp> dv; for (const Disp& d : dsv) dv.push_back({d.ps, d.Q, d.Tn, d.mem, d.tpg, d.ntg, d.tgb, d.per, d.ebo, d.plo, d.pl});
  fjn_unpack(cb, bOut, dv, out, S, &rerun);
  bp_put({bP, bD, bH, bSp, bSg, bPr, bTh, bEB, bOut, bPlAll, bScr});
  S.host_ms += msnow(th0);
}
struct FjnDeleter { std::mutex m; std::condition_variable cv; std::deque<std::shared_ptr<CoJob>> q; bool started = false;
  void push(std::shared_ptr<CoJob>&& j) {
    if (!j) return; std::lock_guard<std::mutex> lk(m);
    if (!started) { started = true; std::thread([this] { for (;;) { std::shared_ptr<CoJob> x; { std::unique_lock<std::mutex> lk2(m); cv.wait(lk2, [&] { return !q.empty(); }); x = std::move(q.front()); q.pop_front(); } x.reset(); } }).detach(); }
    q.push_back(std::move(j)); cv.notify_one(); } };
static FjnDeleter& fjn_deleter() { static FjnDeleter* d = new FjnDeleter(); return *d; }
static void fjn_collect(FjnPending& pend) {
  if (!pend.valid) return;
  TmScope tm_fc_("fc: GPU CB wait");
  using clk_ = std::chrono::steady_clock;
  const double gap_c0 = gap_now();
  const auto tw = clk_::now(); [pend.cb waitUntilCompleted]; pend.S->wait_ms += std::chrono::duration<double, std::milli>(clk_::now() - tw).count();
  const double gap_c1 = gap_now();
  tm_mark("fc: cofit_join wait"); cofit_join(pend.cojob);
  if (g_cblog) { static double base = gap_now(); const double gs = pend.cb.GPUStartTime, ge = pend.cb.GPUEndTime;
    printf("CBLOG%s S%d commit %.4f gpu_start %.4f gpu_end %.4f dur %.2f bins %d gp %zu cp %zu | collect_enter %.4f cb_done_seen %.4f join_done %.4f cofit_done %.4f | abs %.4f %.4f\n", g_gap_in_replay ? "R" : "", pend.stepno, pend.commit_t - base, gs - base, ge - base, (ge - gs) * 1e3, pend.nbins, pend.npg, pend.ncpu,
      gap_c0 - base, gap_c1 - base, gap_now() - base, pend.cojob ? pend.cojob->done_t - base : 0.0, gs, ge); }
  if (g_hc1 & 2) fjn_deleter().push(std::move(pend.cojob)); else pend.cojob.reset();
  tm_mark("fc: fjn_unpack");
  fjn_unpack(pend.cb, pend.bOut, pend.ds, *pend.out, *pend.S, &pend.rerun);
  pend.rerun = nullptr; pend.valid = false; bp_put(pend.keep); pend.keep.clear(); pend.bOut = nil; pend.cb = nil; pend.ds.clear();
}
struct FitTicket { std::mutex m; std::condition_variable cv; bool done = false; std::exception_ptr ex; };
struct FitWorker {
  struct Post { FjnPending* p; std::shared_ptr<FitTicket> tk; };
  std::mutex m; std::condition_variable cv; std::deque<std::function<void()>> runq; std::deque<Post> inflight; bool started = false; qos_class_t qos = QOS_CLASS_UNSPECIFIED;
  std::atomic<long> busy_us{0}, wait_us{0}, njobs{0}, nposts{0}, outstanding{0}; std::mutex qm; std::condition_variable qcv;
  static FitWorker& get() { static FitWorker* w = new FitWorker(); return *w; }
  static bool ready(const Post& po) { return po.p->cb.status >= MTLCommandBufferStatusCompleted && cofit_done(po.p->cojob); }
  static void finish(const std::shared_ptr<FitTicket>& tk) { { std::lock_guard<std::mutex> lk(tk->m); tk->done = true; } tk->cv.notify_all();
    FitWorker& w = get(); { std::lock_guard<std::mutex> lk(w.qm); w.outstanding--; } w.qcv.notify_all(); }
  void quiesce() { std::unique_lock<std::mutex> lk(qm); qcv.wait(lk, [&] { return outstanding.load() == 0; }); }
  void kick() { std::lock_guard<std::mutex> lk(m); cv.notify_all(); }
  void submit(std::function<void()> f) {
    outstanding++; std::lock_guard<std::mutex> lk(m);
    if (!started) { started = true; int rp = 0; pthread_get_qos_class_np(pthread_self(), &qos, &rp);
      g_fjn_notify = [] { FitWorker::get().kick(); }; g_cofit_fin_hook = [] { FitWorker::get().kick(); };
      std::thread([this] { loop(); }).detach(); }
    runq.push_back(std::move(f)); cv.notify_all(); }
  void add_inflight(const Post& po) { std::lock_guard<std::mutex> lk(m); inflight.push_back(po); }
  void post(Post& po) { try { fjn_collect(*po.p); } catch (...) { po.tk->ex = std::current_exception(); } finish(po.tk); nposts++; }
  void loop() {
    if (qos != QOS_CLASS_UNSPECIFIED) pthread_set_qos_class_self_np(qos, 0);
    tl_fjn_worker = true;
    for (;;) {
      std::function<void()> job;
      { std::unique_lock<std::mutex> lk(m);
        for (;;) { if (!runq.empty()) { job = std::move(runq.front()); runq.pop_front(); break; }
          if (!inflight.empty() && ready(inflight.front())) break;
          if (inflight.empty()) cv.wait(lk); else cv.wait_for(lk, std::chrono::milliseconds(2)); } }
      const auto t0 = std::chrono::steady_clock::now();
      if (job) { @autoreleasepool { job(); } njobs++; }
      for (;;) { Post po; { std::lock_guard<std::mutex> lk(m); if (inflight.empty() || !ready(inflight.front())) break; po = inflight.front(); inflight.pop_front(); }
        @autoreleasepool { post(po); } }
      busy_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
    }
  }
};
static void fit_wait(const std::shared_ptr<FitTicket>& tk) {
  const auto t0 = std::chrono::steady_clock::now();
  { std::unique_lock<std::mutex> lk(tk->m); tk->cv.wait(lk, [&] { return tk->done; }); }
  FitWorker::get().wait_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
  if (tk->ex) { FitWorker::get().quiesce(); std::rethrow_exception(tk->ex); }
}
