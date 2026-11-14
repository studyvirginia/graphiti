#include "fjn_h950.hpp"
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <map>
#include "fjn_run_h950.inc"
static int g_fr_maxit = 60, g_fr_dumpat = 0, g_fr_dumpn = 0, g_fr_dumpk = 1, g_fr_verbose = 1;
static bool g_fr_closed = false;
static std::string g_fr_dump;
static id<MTLLibrary> g_fjnL = nil;
static fjn::WarpMat g_fr_warp;
static const evo::Cand* const* g_fit_cands = nullptr;
static float g_fit_yv[96];
static const bool g_fit_chk = getenv("FIT_CHECK") && atoi(getenv("FIT_CHECK")) != 0;
static struct FitChk { long wn = 0, wr3 = 0; double wmax = 0; long n = 0, bad = 0, r3 = 0, r5 = 0, r2 = 0; double maxabs = 0, maxrel = 0, worst_f = 0, worst_r = 0; } g_fchk;
static double fit_recon_rss(const fitx::Side& S, const double* t, const float* Y) {
  std::vector<double> v(S.n.size() + 1); double rss = 0;
  for (int i = 0; i < 48; i++) { const double r = (double)Y[i] - fitx::eval_side(S, t[i], v.data()); rss += r * r; }
  return rss;
}
static void fit_chk_one(double F, double R, double flr) {
  if (!std::isfinite(F) || !std::isfinite(R)) { g_fchk.bad++; return; }
  const double d = std::fabs(F - R), rel = d / std::max(std::max(F, R), flr);
  g_fchk.n++; if (rel > 1e-2) g_fchk.r2++; if (rel > 1e-3) g_fchk.r3++; if (rel > 1e-5) g_fchk.r5++; if (d > g_fchk.maxabs) g_fchk.maxabs = d; if (rel > g_fchk.maxrel) { g_fchk.maxrel = rel; g_fchk.worst_f = F; g_fchk.worst_r = R; }
}
static void fit_chk_print(int tag) {
  if (g_fit_chk && g_fchk.n + g_fchk.bad > 0) { printf("FIT_CHECK well-conditioned (Q<=4, |a|<100): sides %ld rel>1e-3 %ld max rel %.3g\n", g_fchk.wn, g_fchk.wr3, g_fchk.wmax); printf("FIT_CHECK %s: sides %ld nonfinite %ld | rel>1e-5 %ld rel>1e-3 %ld rel>1e-2 %ld | max|rss_kernel - rss_recon| %.6g | max rel (vs max(rss, floor)) %.6g | worst kernel %.9g recon %.9g\n", tag == 1 ? "free-t" : "fixed-t", g_fchk.n, g_fchk.bad, g_fchk.r5, g_fchk.r3, g_fchk.r2, g_fchk.maxabs, g_fchk.maxrel, g_fchk.worst_f, g_fchk.worst_r); }
  g_fchk = FitChk();
}
struct FrStats { long chunks = 0, pairs = 0, disp = 0, fb_fitted = 0, fb_nokernel = 0, fb_tg = 0, fb_enc = 0, fb_fit = 0, bins = 0, iters = 0, dual = 0, plain = 0, stalled = 0, worse = 0; long b_g = 0, b_q = 0, b_tn = 0, b_pso = 0, b_slot = 0, b_tpg = 0, b_nonfin = 0; double gpu_ms = 0, host_ms = 0, wait_ms = 0; };
static size_t g_fjn_tpg = getenv("FJN_TPG") ? (size_t)std::max(1, atoi(getenv("FJN_TPG"))) : 1;
static FrStats g_fr, g_frs;
static const bool g_pprobe = getenv("FJN_PAIRPROBE") && atoi(getenv("FJN_PAIRPROBE")) != 0;
static std::unordered_set<uint64_t> g_pp_seen; static long g_pp_hit = 0, g_pp_all = 0;
static const bool g_glue = getenv("H_GLUE_SKIP") && atoi(getenv("H_GLUE_SKIP")) != 0;
static void glue_fb_pair(const StepBatch& SBr, int a, int b, int c, float* err, float* rx, float* ry) {
  const float fa = SBr.F[0].rfin[(size_t)a], fb = SBr.F[1].rfin[(size_t)b];
  rx[c] = fa; ry[c] = fb;
  if (!isfinite(fa) || !isfinite(fb) || g_fr_score96) { err[c] = DEAD; return; }
  const double s = (double)fa + (double)fb; err[c] = (float)std::min(s, (double)DEAD);
}
static std::mutex g_cfm; static long g_cfn = 0, g_cfdead = 0; static double g_cfmaxabs = 0, g_cfmaxrel = 0, g_cfwg = 0, g_cfwc = 0;
static void cf_record(double gs, double cs) {
  std::lock_guard<std::mutex> lk(g_cfm);
  if (gs >= 1e29 || cs >= 1e29) { if (!(gs >= 1e29 && cs >= 1e29)) g_cfdead++; return; }
  const double d = std::fabs(gs - cs), rel = d / std::max(std::fabs(cs), 1.0);
  g_cfn++; if (d > g_cfmaxabs) { g_cfmaxabs = d; g_cfwg = gs; g_cfwc = cs; } if (rel > g_cfmaxrel) g_cfmaxrel = rel;
}
static bool g_qsplit = getenv("FREET_QSPLIT") && atoi(getenv("FREET_QSPLIT")) != 0;
static id<MTLCommandQueue> g_fjnq = nil;
static id<MTLCommandQueue> fjn_queue(Engine& E) { if (!g_qsplit) return E.cq; if (!g_fjnq) g_fjnq = [E.dev newCommandQueue]; return g_fjnq; }
static const bool g_fjn_q2 = getenv("FJN_Q2") && atoi(getenv("FJN_Q2")) != 0; static id<MTLCommandQueue> g_fjnq2 = nil; static unsigned long g_fjn_qn = 0;
static id<MTLCommandQueue> fjn_queue_piece(Engine& E) { if (!g_fjn_q2 || !g_qsplit) return fjn_queue(E); if (!g_fjnq2) g_fjnq2 = [E.dev newCommandQueue]; return (g_fjn_qn++ & 1u) ? g_fjnq2 : fjn_queue(E); }
static void freet_step_print(int st) {
  if (g_recc2) { size_t h1 = 0, hR = 0; for (uint64_t k : g_rc2_U) { h1 += g_rc2_prev.count(k); hR += g_rc2_all.count(k); }
    printf("RECC2 step %d: gpu pairs %ld | R %ld (%.3f/pair) | U %zu (%.3f/pair) | UP %zu | hit prev step %zu (%.1f %% of U) | hit any earlier %zu (%.1f %%) | calls %ld\n", st, g_rc2_pairs, g_rc2_rec,
           g_rc2_pairs ? (double)g_rc2_rec / g_rc2_pairs : 0.0, g_rc2_U.size(), g_rc2_pairs ? (double)g_rc2_U.size() / g_rc2_pairs : 0.0, g_rc2_UP.size(), h1, g_rc2_U.empty() ? 0.0 : 100.0 * h1 / g_rc2_U.size(), hR, g_rc2_U.empty() ? 0.0 : 100.0 * hR / g_rc2_U.size(), g_rc2_calls);
    for (uint64_t k : g_rc2_U) g_rc2_all.insert(k); g_rc2_prev.swap(g_rc2_U); g_rc2_U.clear(); g_rc2_UP.clear(); g_rc2_rec = g_rc2_pairs = g_rc2_calls = 0; }
  g_gap_stepno = st + 1; g_cbx_step = st + 1;
  if (g_cofit_auto) {
    std::vector<std::pair<double, double>> iv; { std::lock_guard<std::mutex> lk(g_auto_mx); iv.swap(g_auto_iv); }
    std::sort(iv.begin(), iv.end()); double bg = 0, ce = -1; for (auto& x : iv) { if (x.first > ce) { bg += x.second - x.first; ce = x.second; } else if (x.second > ce) { bg += x.second - ce; ce = x.second; } }
    const double bc = g_auto_lane_us.exchange(0) / 1e3; bg *= 1e3; const long ng = g_auto_np.exchange(0), nc = g_cofit_pairs.load();
    int next = g_cofit_pm; double ps = -1;
    if (nc > 0 && ng > 0 && bc > 0 && bg > 0) { const double rc = (double)nc / bc, rg = (double)ng / bg; ps = 1000.0 * rc / (rc + rg);
      static int nsamp = 0; const double w = nsamp++ == 0 ? 1.0 : 0.5;
      const double tgt = (1.0 - w) * g_cofit_pm + w * ps; next = std::max(0, std::min(900, (int)(tgt / 25.0 + 0.5) * 25)); }
    else if (g_cofit_pm == 0 && g_cofit_auto >= 2) next = 100;
    printf("COFIT_AUTO step %d: share %d | lane %ld pairs in %.1f ms wall-busy | gpu %ld pairs in %.1f ms (union of %zu CBs) | balance share %.0f -> next %d%s\n", st, g_cofit_pm, nc, bc, ng, bg, iv.size(), ps, next, g_cofit_auto >= 2 ? " (applied)" : "");
    if (g_cofit_auto >= 2) g_cofit_pm = next; }
  if (g_bcount) { size_t nu; { std::lock_guard<std::mutex> lk(g_bc_m); nu = g_bc_u.size(); g_bc_u.clear(); }
    printf("BCOUNT step %d: lane builds %ld distinct %zu build cpu %.1f ms | cofit_start calls %ld nodes %ld main cpu %.1f ms wall %.1f ms\n", st, g_bc_builds.exchange(0), nu, g_bc_bns.exchange(0) / 1e6, g_bc_calls.exchange(0), g_bc_nodes.exchange(0), g_bc_sns.exchange(0) / 1e6, g_bc_swall.exchange(0) / 1e6); }
  if (g_cofit_pm > 0 && g_cofit_pairs.load() > 0) { printf("COFIT step %d: permille %d cpu pairs %ld jobs %ld cpu_busy %.1f ms (sum over workers) join_wait %.1f ms encfail %ld scorer %.1f ms\n", st, g_cofit_pm, g_cofit_pairs.load(), g_cofit_jobs.load(), g_cofit_busy_us.load() / 1e3, g_cofit_wait_us.load() / 1e3, g_cofit_encfail.load(), g_cofit_score_us.exchange(0) / 1e3); for (int w = 0; w < std::min(64, g_cofit_nwk.load()); w++) { CoWk& W = g_cofit_wk[w]; printf("COWK step %d: w%d %s pairs %ld wall %.1f ms cpu %.1f ms\n", st, w, W.ew.load() ? "E" : "P", W.pairs.exchange(0), W.wall_us.exchange(0) / 1e3, W.cpu_us.exchange(0) / 1e3); } { mach_timebase_info_data_t tb; mach_timebase_info(&tb); const double k = (double)tb.numer / (double)tb.denom / 1e6; printf("COPH step %d:", st); for (int q = 0; q < 8; q++) printf(" ph%d %.1f", q, (double)g_cofit_ph[q].exchange(0) * k); printf(" ms | plain %ld dual %ld solves %ld it %ld | cone0 %ld cone1 %ld\n", g_cofit_ph[8].exchange(0), g_cofit_ph[9].exchange(0), g_cofit_ph[10].exchange(0), g_cofit_ph[11].exchange(0), g_cofit_ph[12].exchange(0), g_cofit_ph[13].exchange(0)); } g_cofit_pairs = 0; g_cofit_jobs = 0; g_cofit_busy_us = 0; g_cofit_wait_us = 0; g_cofit_encfail = 0; }
  if (g_hc1 & 1) { FitWorker& w = FitWorker::get(); printf("FITW step %d: fit worker busy %.1f ms (run jobs %ld, posts %ld; this is the fit host path moved off the main thread) | main waited on piece tickets %.1f ms\n", st, (double)w.busy_us.exchange(0) / 1e3, w.njobs.exchange(0), w.nposts.exchange(0), (double)w.wait_us.exchange(0) / 1e3); }
  if (g_cfchk > 0 && g_cfn + g_cfdead > 0) { printf("CISFUSE_CHECK step %d: checked %ld | max|gpu-cpu| %.6g | max rel (|d|/max(|cpu|,1)) %.6g | dead-mismatch %ld | worst gpu %.9g cpu %.9g\n", st, g_cfn, g_cfmaxabs, g_cfmaxrel, g_cfdead, g_cfwg, g_cfwc); g_cfn = g_cfdead = 0; g_cfmaxabs = g_cfmaxrel = 0; }
  FrStats& s = g_frs;
  if (g_pprobe && g_pp_all) { printf("PAIRPROBE step %d: pairs %ld already fitted earlier in the run %ld (%.1f%%) | seen %zu\n", st, g_pp_all, g_pp_hit, 100.0 * (double)g_pp_hit / (double)g_pp_all, g_pp_seen.size()); g_pp_all = g_pp_hit = 0; }
  if (g_freet && g_fr_verbose && s.chunks) printf("H950 step %d: chunks %ld pairs %ld dispatched %ld bins %ld | fallback nokernel %ld tg %ld enc %ld nonfinite %ld | LM iters mean %.1f dual %.1f plain %.1f stalled %ld | joint worse than start %ld | gpu %.2f ms host %.2f ms wait %.2f ms\n",
    st, s.chunks, s.pairs, s.disp, s.bins, s.fb_nokernel, s.fb_tg, s.fb_enc, s.fb_fit, s.disp ? (double)s.iters / s.disp : 0.0, s.disp ? (double)s.dual / s.disp : 0.0, s.disp ? (double)s.plain / s.disp : 0.0, s.stalled, s.worse, s.gpu_ms, s.host_ms, s.wait_ms);
  if (g_freet && g_fr_verbose && s.chunks) printf("H950 step %d branches: is_g %ld Q %ld Tn>max %ld pso %ld slot>32K %ld tpg<1 %ld nonfinite %ld\n", st, s.b_g, s.b_q, s.b_tn, s.b_pso, s.b_slot, s.b_tpg, s.b_nonfin);
  if (g_freet && s.chunks) {
    const long skipped = s.pairs - s.disp - s.fb_fitted;
    { static double pv[2] = {0, 0}; double cu[2]; cofit_cput(cu); printf("CPUT call %d: process cpu user %.1f sys %.1f ms since the previous fit call\n", st, cu[0] - pv[0], cu[1] - pv[1]); pv[0] = cu[0]; pv[1] = cu[1]; }
    cofit_thrt(st);
    printf("BPC call %d: requested shared %llu private %llu | fresh shared %llu private %llu bytes\n", st, g_bpc[0].load(), g_bpc[1].load(), g_bpc[2].load(), g_bpc[3].load());
    if (g_fr_verbose) printf("H950 step %d pairs check: pairs %ld = dispatched %ld + fallback_fitted %ld (skipped %ld)\n", st, s.pairs, s.disp, s.fb_fitted, skipped);
    if (skipped != 0) { fprintf(stderr, "H950 ASSERT step %d: pairs %ld != dispatched %ld + fallback_fitted %ld\n", st, s.pairs, s.disp, s.fb_fitted); abort(); }
  }
  g_fr.fb_fitted += s.fb_fitted; g_fr.b_g += s.b_g; g_fr.b_q += s.b_q; g_fr.b_tn += s.b_tn; g_fr.b_pso += s.b_pso; g_fr.b_slot += s.b_slot; g_fr.b_tpg += s.b_tpg; g_fr.b_nonfin += s.b_nonfin;
  g_fr.chunks += s.chunks; g_fr.pairs += s.pairs; g_fr.disp += s.disp; g_fr.fb_nokernel += s.fb_nokernel; g_fr.fb_tg += s.fb_tg; g_fr.fb_enc += s.fb_enc; g_fr.fb_fit += s.fb_fit; g_fr.bins += s.bins; g_fr.iters += s.iters; g_fr.dual += s.dual; g_fr.plain += s.plain; g_fr.stalled += s.stalled; g_fr.worse += s.worse; g_fr.gpu_ms += s.gpu_ms; g_fr.host_ms += s.host_ms; g_fr.wait_ms += s.wait_ms;
  s = FrStats();
}
static bool g_freet_async = getenv("FREET_ASYNC") && atoi(getenv("FREET_ASYNC")) != 0;
struct FrOutPool { std::mutex m; std::vector<std::vector<FrOut>> v; size_t bytes = 0; };
static FrOutPool& frout_pool() { static FrOutPool* p = new FrOutPool(); return *p; }
static std::vector<FrOut> frout_take(size_t n) {
  FrOutPool& P = frout_pool(); std::lock_guard<std::mutex> lk(P.m); int best = -1;
  for (int i = 0; i < (int)P.v.size(); i++) if (P.v[(size_t)i].capacity() >= n && (best < 0 || P.v[(size_t)i].capacity() < P.v[(size_t)best].capacity())) best = i;
  if (best < 0) return std::vector<FrOut>();
  std::vector<FrOut> r = std::move(P.v[(size_t)best]); P.v.erase(P.v.begin() + best); P.bytes -= r.capacity() * sizeof(FrOut); r.clear(); return r; }
static void frout_give(std::vector<FrOut>&& o) {
  if (o.capacity() == 0) return; FrOutPool& P = frout_pool(); std::lock_guard<std::mutex> lk(P.m);
  if (P.bytes + o.capacity() * sizeof(FrOut) > ((size_t)256 << 20)) return;
  P.bytes += o.capacity() * sizeof(FrOut); o.clear(); P.v.push_back(std::move(o)); }
struct FreetDeferred { ~FreetDeferred() { if ((g_hc2 & 16) && outs.capacity()) frout_give(std::move(outs)); } FjnPending p; FrRunStats rs; std::vector<FrOut> outs; std::vector<FrPairIn> prs; std::vector<int> child, ua, ub; std::shared_ptr<FitTicket> tk; std::vector<float> thsnap; const StepBatch* sbp = nullptr; float* err = nullptr; float* rx = nullptr; float* ry = nullptr; bool any = false;
  std::shared_ptr<CoJob> cojob; std::vector<std::shared_ptr<struct CbCall>> calls; uint32_t npend = 0; bool cbk = false; };
static void cbk_prepare(FreetDeferred& D);
static FreetDeferred g_fdq[2]; static int g_fdi = 0;
static const int g_hc3 = getenv("FJN_HC3") ? atoi(getenv("FJN_HC3")) : 0;
static void freet_writeback_outs(const std::vector<FrOut>& outs, const std::vector<int>& child, const std::vector<int>& ua, const std::vector<int>& ub, const std::vector<FrPairIn>& prs,
                                 const StepBatch& SBr, float* err, float* rx, float* ry, FrStats& S, const Unique& UQ, const float* yvv, float flrs, FreetDeferred* own = nullptr) {
  TmScope tm_wb_("wb: score+install (pool)");
  auto body = [&](size_t k, FrStats& A) {
    const FrOut& o = outs[k]; const int c = child[k];
    H308T('F', UQ.hash[(size_t)ua[k]], UQ.hash[(size_t)ub[k]], (const float*)&o, sizeof(FrOut) / 4);
#ifdef H308_TRACE
    { const float* fo = (const float*)&o; const size_t nraw = (size_t)((const float*)&o.raw[0] - fo); H308T('E', UQ.hash[(size_t)ua[k]], UQ.hash[(size_t)ub[k]], fo, nraw); H308T('W', UQ.hash[(size_t)ua[k]], UQ.hash[(size_t)ub[k]], (const float*)&o.raw[0], 128); }
#endif
    { const FrPairIn& pp = prs[k]; float f[8] = {(float)pp.x.Q, (float)pp.y.Q, (float)pp.x.nth, (float)pp.y.nth, (float)pp.x.nS, (float)pp.y.nS, (float)pp.x.nslot, (float)pp.y.nslot}; H308T('I', UQ.hash[(size_t)ua[k]], UQ.hash[(size_t)ub[k]], f, 8);
      H308T('P', UQ.hash[(size_t)ua[k]], 0, (const float*)pp.x.prog, pp.x.nprog); H308T('P', UQ.hash[(size_t)ub[k]], 0, (const float*)pp.y.prog, pp.y.nprog); }
    if (!o.status) { A.fb_fitted++; if (g_glue) glue_fb_pair(SBr, ua[k], ub[k], c, err, rx, ry); return; }
    if (!isfinite(o.Fx) || !isfinite(o.Fy)) { A.fb_fit++; A.fb_fitted++; if (g_glue) glue_fb_pair(SBr, ua[k], ub[k], c, err, rx, ry); return; }
    const double s = (double)o.Fx + (double)o.Fy; const double s0 = (double)SBr.F[0].rfin[(size_t)ua[k]] + (double)SBr.F[1].rfin[(size_t)ub[k]];
    if (s > s0 * (1.0 + 1e-4)) A.worse++;
    rx[c] = o.Fx; ry[c] = o.Fy; err[c] = (float)std::min(s, (double)DEAD);
    if (g_cf_on) {
      err[c] = (float)std::min((double)o.cf, (double)DEAD);
      if (g_cfchk > 0 && (int)(k % (size_t)g_cfchk) == 0) cf_record((double)err[c], (double)(float)std::min(g_fr_score96(o.res48), (double)DEAD));
    } else
    if (g_fr_score96) err[c] = (float)std::min(g_fr_score96(o.res48), (double)DEAD);
    A.disp++; A.iters += (long)o.it; A.dual += (long)o.dual; A.plain += (long)o.plain; if (o.stalled != 0.0f) A.stalled++;
  };
  auto addc = [](FrStats& D, const FrStats& A) { D.fb_fitted += A.fb_fitted; D.fb_fit += A.fb_fit; D.worse += A.worse; D.disp += A.disp; D.iters += A.iters; D.dual += A.dual; D.plain += A.plain; D.stalled += A.stalled; };
  if (g_wb_par_on && g_wb_par && prs.size() >= 256) {
    const int NR = 8; std::vector<FrStats> part((size_t)NR);
    g_wb_par(NR, [&](int r) { const size_t n = prs.size(); for (size_t k = n * (size_t)r / NR; k < n * (size_t)(r + 1) / NR; k++) body(k, part[(size_t)r]); });
    for (int r = 0; r < NR; r++) addc(S, part[(size_t)r]);
  } else {
    for (size_t k = 0; k < prs.size(); k++) body(k, S);
  }
  (void)UQ; (void)yvv; (void)flrs;
  tm_mark("wb: fitx capture submit");
  const auto fx_t0 = std::chrono::steady_clock::now();
  struct FxT { std::chrono::steady_clock::time_point t; ~FxT() { fitx::g_ms_wb += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); } } fx_tm{fx_t0};
  static const bool fx_skip_ = getenv("FX_SKIP") && atoi(getenv("FX_SKIP")) != 0;
  if (fitx::on() && g_fit_cands && !fx_skip_) {
    static const bool fx_async = getenv("FX_ASYNC") && atoi(getenv("FX_ASYNC")) == 1;
    if (!fx_async || g_fit_chk) fitx::retain_join();
    tm_mark("wb: fitx retain_join+alloc");
    double tw[48]; fitx::Rec r; const size_t NK = prs.size();
    std::vector<uint64_t> kk(NK), hx(NK), hy(NK);
    tm_mark("wb: fitx batch_gate");
    const float bg = fitx::batch_gate(NK, [&](size_t k) { const FrOut& o = outs[k]; return (o.status && isfinite(o.Fx) && isfinite(o.Fy)) ? err[child[k]] : INFINITY; }, 1);
    tm_mark("wb: fitx job build+submit");
    if (fx_async && !g_fit_chk && (g_hc3 & 2) && own) {
      struct FxOwn { ~FxOwn() { if ((g_hc2 & 16) && outs.capacity()) frout_give(std::move(outs)); } std::vector<FrOut> outs; std::vector<FrPairIn> prs; std::vector<int> child; std::vector<float> errs; float bg = 0; size_t ne = 0; };
      auto ow = std::make_shared<FxOwn>(); ow->errs.resize(NK); ow->bg = bg;
      for (size_t k = 0; k < NK; k++) { const FrOut& o = outs[k]; const float ev = err[child[k]]; ow->errs[k] = ev; if (o.status && isfinite(o.Fx) && isfinite(o.Fy) && isfinite(ev) && !(ev > bg)) ow->ne++; }
      const evo::Cand* const* cands = g_fit_cands;
      if (ow->ne) { ow->outs = std::move(own->outs); ow->prs = std::move(own->prs); ow->child = std::move(own->child);
        fitx::cap_submit([ow, cands]() {
          fitx::retain_join();
          fitx::st().cur.reserve(fitx::st().cur.size() + ow->ne);
          for (size_t k = 0; k < ow->errs.size(); k++) { const FrOut& o = ow->outs[k]; const float ev = ow->errs[k];
            if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(ev) || ev > ow->bg) continue;
            const int c = ow->child[k]; FrTreeRef xr = ow->prs[k].x, yr = ow->prs[k].y;
            const uint64_t key = fitx::gkey(*cands[c]), h0 = fitx::shash(xr), h1 = fitx::shash(yr);
            fitx::capture(c, key, ev, 1, &xr, h0, o.th[0], o.a[0], &yr, h1, o.th[1], o.a[1], o.u); } }); }
    } else
    if (fx_async && !g_fit_chk) {
      struct FxItem { int c; float err; FrTreeRef x, y; float th[2][32]; float a[2][16]; float u[10]; };
      auto job = std::make_shared<std::vector<FxItem>>();
      if (g_hc3 & 1) { size_t ne = 0; for (size_t k = 0; k < NK; k++) { const FrOut& o = outs[k]; if (o.status && isfinite(o.Fx) && isfinite(o.Fy) && isfinite(err[child[k]]) && !(err[child[k]] > bg)) ne++; } job->reserve(ne); }
      else job->reserve(NK);
      for (size_t k = 0; k < NK; k++) { const FrOut& o = outs[k]; const int c = child[k];
        if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(err[c]) || err[c] > bg) continue;
        job->emplace_back(); FxItem& it = job->back(); it.c = c; it.err = err[c]; it.x = prs[k].x; it.y = prs[k].y;
        std::memcpy(it.th, o.th, sizeof it.th); std::memcpy(it.a, o.a, sizeof it.a); std::memcpy(it.u, o.u, sizeof it.u); }
      const evo::Cand* const* cands = g_fit_cands;
      if (!job->empty()) fitx::cap_submit([job, cands]() {
        fitx::retain_join();
        fitx::st().cur.reserve(fitx::st().cur.size() + job->size());
        for (FxItem& it : *job) {
          const uint64_t key = fitx::gkey(*cands[it.c]), h0 = fitx::shash(it.x), h1 = fitx::shash(it.y);
          fitx::capture(it.c, key, it.err, 1, &it.x, h0, it.th[0], it.a[0], &it.y, h1, it.th[1], it.a[1], it.u); } });
    } else {
    auto prep = [&](size_t k) { const FrOut& o = outs[k]; const int c = child[k]; if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(err[c]) || err[c] > bg) return;
      kk[k] = fitx::gkey(*g_fit_cands[c]); hx[k] = fitx::shash(prs[k].x); hy[k] = fitx::shash(prs[k].y); };
    if (g_wb_par_on && g_wb_par && NK >= 256) g_wb_par(8, [&](int rr) { for (size_t k = NK * (size_t)rr / 8; k < NK * (size_t)(rr + 1) / 8; k++) prep(k); });
    else for (size_t k = 0; k < NK; k++) prep(k);
    std::vector<std::shared_ptr<const fitx::Side>> fx_keep;
    { static const bool fx_par = getenv("FX_PAR") && atoi(getenv("FX_PAR")) == 1;
      if (fx_par && g_wb_par_on && g_wb_par && NK >= 256) {
        fitx::Store& FS = fitx::st(); struct PJob { uint64_t h; const FrTreeRef* t; }; std::vector<PJob> jobs; std::unordered_set<uint64_t> seenh; seenh.reserve(NK);
        for (size_t k = 0; k < NK; k++) { const FrOut& o = outs[k]; const int c = child[k];
          if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(err[c]) || err[c] > bg) continue;
          const uint64_t hh[2] = {hx[k], hy[k]}; const FrTreeRef* tt[2] = {&prs[k].x, &prs[k].y};
          for (int sd = 0; sd < 2; sd++) { auto ii = FS.intern.find(hh[sd]); if (ii != FS.intern.end() && !ii->second.expired()) continue; if (seenh.insert(hh[sd]).second) jobs.push_back({hh[sd], tt[sd]}); } }
        if (!jobs.empty()) {
          std::vector<std::shared_ptr<const fitx::Side>> built(jobs.size());
          g_wb_par(8, [&](int rr) { const size_t n = jobs.size(); for (size_t i = n * (size_t)rr / 8; i < n * (size_t)(rr + 1) / 8; i++) { float z[32] = {0}; auto sd = std::make_shared<fitx::Side>(); if (fitx::build_side(*jobs[i].t, z, z, *sd)) built[i] = sd; } });
          for (size_t i = 0; i < jobs.size(); i++) if (built[i]) FS.intern[jobs[i].h] = built[i];
          fx_keep = std::move(built);
        } } }
    fitx::st().cur.reserve(fitx::st().cur.size() + NK);
    for (size_t k = 0; k < prs.size(); k++) {
      const FrOut& o = outs[k]; const int c = child[k];
      if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(err[c]) || err[c] > bg) continue;
      if (!fitx::capture(c, kk[k], err[c], 1, &prs[k].x, hx[k], o.th[0], o.a[0], &prs[k].y, hy[k], o.th[1], o.a[1], o.u)) continue;
      if (g_fit_chk) {
        fitx::materialize(fitx::st().cur.back(), r);
        double U = 0; for (int q = 0; q < 9; q++) U += (double)o.u[q] * o.u[q];
        for (int i = 0; i < 48; i++) { double a = 0; for (int q = 0; q < 9; q++) a += g_fr_warp.Sd[i * 9 + q] * (double)o.u[q] * o.u[q]; tw[i] = (double)o.u[9] + a / U; }
        double yyx = 0, yyy = 0; for (int i = 0; i < 48; i++) { yyx += (double)g_fit_yv[i] * g_fit_yv[i]; yyy += (double)g_fit_yv[48 + i] * g_fit_yv[48 + i]; }
        fit_chk_one((double)o.Fx, fit_recon_rss(r.s[0], tw, g_fit_yv), 1e-5 * yyx); fit_chk_one((double)o.Fy, fit_recon_rss(r.s[1], tw, g_fit_yv + 48), 1e-5 * yyy);
        for (int sd = 0; sd < 2; sd++) { double am = 0; for (int q = 0; q < r.s[sd].Q; q++) am = std::max(am, (double)std::fabs(r.s[sd].a[q]));
          if (r.s[sd].Q <= 4 && am < 100) { const double F = sd ? o.Fy : o.Fx, R = fit_recon_rss(r.s[sd], tw, g_fit_yv + 48 * sd), fl = 1e-5 * (sd ? yyy : yyx); const double rel = std::fabs(F - R) / std::max(std::max(F, R), fl); g_fchk.wn++; if (rel > 1e-3) g_fchk.wr3++; g_fchk.wmax = std::max(g_fchk.wmax, rel); } }
        { static int dbg = 0; const double Fs[2] = {o.Fx, o.Fy}; const double Rs[2] = {fit_recon_rss(r.s[0], tw, g_fit_yv), fit_recon_rss(r.s[1], tw, g_fit_yv + 48)};
          for (int sd = 0; sd < 2 && getenv("FIT_CHECK_DBG") && dbg < 12; sd++) if (std::isfinite(Rs[sd]) && std::fabs(Fs[sd] - Rs[sd]) > 1e-2 * std::max(Fs[sd], 1e-3)) { dbg++;
            printf("FIT_CHECK_DBG side %d Q %d nth %d kernel %.6g recon %.6g a [", sd, r.s[sd].Q, r.s[sd].nth, Fs[sd], Rs[sd]); for (int q = 0; q < r.s[sd].Q; q++) printf("%.4g ", (double)r.s[sd].a[q]); printf("] th ["); for (int q = 0; q < r.s[sd].nth; q++) printf("%.4g ", (double)r.s[sd].th[q]); printf("] u9 %.4g\n", (double)o.u[9]); } }
      }
    }
    }
  }
}
static void freet_writeback(FreetDeferred& D, FrStats& S) {
  const auto tr0 = clk::now();
  static FrStats* unused = nullptr; (void)unused;
  freet_writeback_outs(D.outs, D.child, D.ua, D.ub, D.prs, *D.sbp, D.err, D.rx, D.ry, S, D.sbp->UQ, nullptr, 0.0f, &D);
  S.host_ms += ms_since(tr0);
}
static void freet_drain(FreetDeferred& D) {
  if (!D.any) return;
  TmScope tm_fd_("fd: collect");
  FrStats& S = g_frs; const FrRunStats& r = D.rs;
  if (D.cbk) cbk_prepare(D);
  else if (D.tk) { tm_mark("fd: fit-worker wait"); std::shared_ptr<FitTicket> tk = std::move(D.tk); D.tk.reset(); fit_wait(tk); }
  else fjn_collect(D.p);
  tm_mark("fd: stats");
  S.fb_nokernel += r.nokernel; S.fb_tg += r.tg; S.fb_enc += r.enc; S.bins += r.bins; S.gpu_ms += r.gpu_ms; S.host_ms += r.host_ms; S.wait_ms += r.wait_ms;
  S.b_q += r.b_q; S.b_tn += r.b_tn; S.b_pso += r.b_pso; S.b_slot += r.b_slot; S.b_tpg += r.b_tpg;
  freet_writeback(D, S); D.any = false;
}
static bool g_stream = getenv("STREAM") && atoi(getenv("STREAM")) != 0 && g_freet_async;
static int g_stream_n = getenv("STREAM_N") ? std::max(1, atoi(getenv("STREAM_N"))) : 4;
static const bool g_stream_auto = !getenv("STREAM_N");
static const size_t g_stream_cap = getenv("STREAM_CAP") ? (size_t)std::max(1, atoi(getenv("STREAM_CAP"))) : 1200;
struct StreamPiece { std::unique_ptr<FreetDeferred> d; long chunk; };
static std::vector<StreamPiece> g_stp; static long g_st_chunk = 0;
static void stream_drain(long keep_from) {
  size_t k = 0; while (k < g_stp.size() && g_stp[k].chunk < keep_from) { freet_drain(*g_stp[k].d); k++; }
  g_stp.erase(g_stp.begin(), g_stp.begin() + (long)k);
}
static const int g_holdk = getenv("FJN_HOLDK") ? std::max(0, atoi(getenv("FJN_HOLDK"))) : 0;
static const int g_hold_edge = getenv("FJN_HOLD_EDGE") ? std::max(0, atoi(getenv("FJN_HOLD_EDGE"))) : 0;
static const size_t g_holdk_maxgp = getenv("FJN_HOLDK_MAXGP") ? (size_t)std::max(0, atoi(getenv("FJN_HOLDK_MAXGP"))) : 0;
struct HoldState { std::unique_ptr<FreetDeferred> d; long chunk0 = -1, chunk1 = -1; int npiece = 0; float yv[96]; float flr = 0; Engine* E = nullptr; };
static HoldState g_hold;
static void freet_fjn_submit(Engine& E, FreetDeferred& N, const float* yv96, float flr) {
  id<MTLCommandQueue> pq_ = fjn_queue_piece(E);
  if ((g_hc2 & 16) && N.outs.empty()) N.outs = frout_take(N.prs.size());
  if (!(g_hc1 & 1) || g_gap_collect_on) {
    fjn_run(E.dev, pq_, [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&, pq_]() { return h182_cb(pq_); }, N.prs, yv96, flr, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, N.rs, g_fjn_tpg, &N.p);
    return; }
  const size_t NPR = N.prs.size(); size_t nf = 0;
  for (size_t k = 0; k < NPR; k++) { const FrPairIn& p = N.prs[k]; if (p.thx) nf += (size_t)p.x.nth; if (p.thy) nf += (size_t)p.y.nth; }
  N.thsnap.assign(nf + 1, 0.0f); size_t off = 0;
  for (size_t k = 0; k < NPR; k++) { FrPairIn& p = N.prs[k];
    if (p.thx) { if (p.x.nth > 0) std::memcpy(&N.thsnap[off], p.thx, (size_t)p.x.nth * 4); p.thx = &N.thsnap[off]; off += (size_t)p.x.nth; }
    if (p.thy) { if (p.y.nth > 0) std::memcpy(&N.thsnap[off], p.thy, (size_t)p.y.nth * 4); p.thy = &N.thsnap[off]; off += (size_t)p.y.nth; } }
  auto tk = std::make_shared<FitTicket>(); N.tk = tk;
  std::array<float, 96> yv; std::memcpy(yv.data(), yv96, 96 * 4); const int stepno = g_gap_stepno; (void)stepno;
  FreetDeferred* Np = &N; Engine* Ep = &E;
  FitWorker::get().submit([Np, Ep, pq_, yv, flr, tk]() {
    try {
      FreetDeferred& N = *Np; Engine& E = *Ep;
      fjn_run(E.dev, pq_, [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&, pq_]() { return h182_cb(pq_); }, N.prs, yv.data(), flr, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, N.rs, g_fjn_tpg, &N.p);
      if (N.p.valid) FitWorker::get().add_inflight({&N.p, tk}); else FitWorker::finish(tk);
    } catch (...) { tk->ex = std::current_exception(); FitWorker::finish(tk); } });
}
static void hold_commit() {
  if (!g_hold.d) return;
  Engine& E = *g_hold.E; StreamPiece sp; sp.chunk = g_hold.chunk1; sp.d = std::move(g_hold.d); FreetDeferred& N = *sp.d; N.any = true;
  freet_fjn_submit(E, N, g_hold.yv, g_hold.flr);
  g_stp.push_back(std::move(sp)); g_hold.chunk0 = g_hold.chunk1 = -1; g_hold.npiece = 0;
}
static const size_t g_cbk = getenv("FJN_CBUCKET") ? (size_t)std::max(0, atoi(getenv("FJN_CBUCKET"))) : 0;
static const bool g_cbk_on = g_cbk > 0 && g_stream && g_holdk == 0 && !(g_hc1 & 1);
static const size_t g_cbk_cap = getenv("FJN_CBUCKET_CAP") ? (size_t)std::max(1, atoi(getenv("FJN_CBUCKET_CAP"))) : g_stream_cap;
static const long g_cbk_w = getenv("FJN_CBUCKET_W") ? (long)std::max(1, atoi(getenv("FJN_CBUCKET_W"))) : 1;
struct CbItem { FreetDeferred* P; uint32_t k; };
struct CbCall {
  ~CbCall() { if ((g_hc2 & 16) && outs.capacity()) frout_give(std::move(outs)); }
  FjnPending p; FrRunStats rs; std::vector<FrOut> outs; std::vector<FrPairIn> prs; std::vector<CbItem> items; bool done = false; };
struct CbBucket { std::vector<CbItem> items; long oldest = 0; };
struct CbState { std::map<int, CbBucket> b; float yv[96]; float flr = 0; const StepBatch* sbp = nullptr; Engine* E = nullptr; size_t npend = 0;
  long st_calls = 0, st_pairs = 0, st_classes = 0, st_maxcls = 0, st_early = 0, st_forced = 0; double st_copy_ms = 0; };
static CbState g_cbs;
static inline int cbk_key(const FrPairIn& p) { return std::max(p.x.Q, p.y.Q) * 1024 + std::max(p.x.nth, p.y.nth); }
static void cbk_collect(CbCall& c) {
  if (c.done) return;
  fjn_collect(c.p);
  tm_mark("cb: copy-back");
  const auto t0 = clk::now();
  for (size_t i = 0; i < c.items.size(); i++) c.items[i].P->outs[c.items[i].k] = c.outs[i];
  g_cbs.st_copy_ms += ms_since(t0);
  FrStats& S = g_frs; const FrRunStats& r = c.rs;
  S.fb_nokernel += r.nokernel; S.fb_tg += r.tg; S.fb_enc += r.enc; S.bins += r.bins; S.gpu_ms += r.gpu_ms; S.host_ms += r.host_ms; S.wait_ms += r.wait_ms;
  S.b_q += r.b_q; S.b_tn += r.b_tn; S.b_pso += r.b_pso; S.b_slot += r.b_slot; S.b_tpg += r.b_tpg;
  c.done = true; std::vector<FrOut>().swap(c.outs); std::vector<FrPairIn>().swap(c.prs); std::vector<CbItem>().swap(c.items);
}
static void cbk_emit(const std::vector<int>& keys, const std::function<void()>& after_first = nullptr) {
  if (keys.empty()) { if (after_first) after_first(); return; }
  Engine& E = *g_cbs.E;
  std::vector<CbItem> all; std::vector<int> ak;
  for (int key : keys) { auto it = g_cbs.b.find(key); if (it == g_cbs.b.end()) continue; for (const CbItem& x : it->second.items) { all.push_back(x); ak.push_back(key); } g_cbs.b.erase(it); }
  const size_t n = all.size(); if (!n) { if (after_first) after_first(); return; }
  bool af_done = !after_first;
  const size_t nc = (n + g_cbk_cap - 1) / g_cbk_cap;
  for (size_t m = 0; m < nc; m++) {
    const size_t lo = n * m / nc, hi = n * (m + 1) / nc; if (lo == hi) continue;
    auto c = std::make_shared<CbCall>(); c->prs.reserve(hi - lo); c->items.assign(all.begin() + (long)lo, all.begin() + (long)hi);
    long ncls = 0; int pk = -1;
    for (size_t i = lo; i < hi; i++) { FreetDeferred* P = all[i].P; c->prs.push_back(P->prs[all[i].k]); if (ak[i] != pk) { ncls++; pk = ak[i]; }
      P->npend--; if (P->calls.empty() || P->calls.back().get() != c.get()) P->calls.push_back(c); }
    g_cbs.npend -= hi - lo;
    g_cbs.st_calls++; g_cbs.st_pairs += (long)(hi - lo); g_cbs.st_classes += ncls; g_cbs.st_maxcls = std::max(g_cbs.st_maxcls, ncls);
    tm_mark("cb: fjn_run (enc+commit)");
    id<MTLCommandQueue> pq_ = fjn_queue_piece(E);
    g_cofit_bypass = true;
    fjn_run(E.dev, pq_, [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&, pq_]() { return h182_cb(pq_); }, c->prs, g_cbs.yv, g_cbs.flr, g_fr_closed, g_fr_maxit, g_fr_warp, c->outs, c->rs, g_fjn_tpg, &c->p);
    g_cofit_bypass = false;
    if (!af_done) { af_done = true; after_first(); }
  }
  if (!af_done) after_first();
}
static void cbk_flush_all() { std::vector<int> ks; for (auto& kv : g_cbs.b) ks.push_back(kv.first); cbk_emit(ks); }
static void cbk_prepare(FreetDeferred& D) {
  if (D.npend) cbk_flush_all();
  for (auto& c : D.calls) cbk_collect(*c);
  D.calls.clear();
  if (D.cojob) { tm_mark("fc: cofit_join wait"); cofit_join(D.cojob); D.cojob.reset(); }
}
static void cbk_chunk(Engine& E, const StepBatch& SBr, long ck, const std::vector<FrPairIn>& prs, const std::vector<int>& child, const std::vector<int>& ua, const std::vector<int>& ub,
                      const float* yvv, float flrs, float* err, float* rx, float* ry, size_t M) {
  if (g_cbs.sbp && (g_cbs.sbp != &SBr || std::memcmp(g_cbs.yv, yvv, 96 * 4) != 0 || g_cbs.flr != flrs)) cbk_flush_all();
  std::memcpy(g_cbs.yv, yvv, 96 * 4); g_cbs.flr = flrs; g_cbs.sbp = &SBr; g_cbs.E = &E;
  const size_t NP_ = prs.size();
  std::vector<std::vector<uint8_t>> lane_sel; std::vector<FreetDeferred*> lane_pc;
  for (size_t m = 0; m < M; m++) {
    const size_t lo = NP_ * m / M, hi = NP_ * (m + 1) / M; if (lo == hi) continue;
    StreamPiece sp; sp.chunk = ck; sp.d.reset(new FreetDeferred()); FreetDeferred& N = *sp.d;
    N.any = true; N.cbk = true; N.prs.assign(prs.begin() + (long)lo, prs.begin() + (long)hi); N.child.assign(child.begin() + (long)lo, child.begin() + (long)hi);
    N.ua.assign(ua.begin() + (long)lo, ua.begin() + (long)hi); N.ub.assign(ub.begin() + (long)lo, ub.begin() + (long)hi);
    N.sbp = &SBr; N.err = err; N.rx = rx; N.ry = ry;
    const size_t n = N.prs.size();
    if ((g_hc2 & 16) && N.outs.empty()) N.outs = frout_take(n);
    N.outs.assign(n, FrOut());
    std::vector<uint8_t> sel;
    if (g_cofit_pm > 0 && !g_cofit_bypass) { sel.assign(n, 0);
      if (g_aghc & 6) { tm_mark("cb: cofit_select"); cofit_select(N.prs, sel, FJN_TMAX); lane_sel.push_back(sel); lane_pc.push_back(&N); }
      else { tm_mark("cb: cofit_start"); N.cojob = cofit_start(N.prs, sel, yvv, flrs, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, FJN_TMAX); } }
    tm_mark("cb: route");
    for (size_t k = 0; k < n; k++) { if (!sel.empty() && sel[k]) continue;
      CbBucket& bk = g_cbs.b[cbk_key(N.prs[k])]; if (bk.items.empty()) bk.oldest = ck; bk.items.push_back({&N, (uint32_t)k}); N.npend++; g_cbs.npend++; }
    g_stp.push_back(std::move(sp));
  }
  std::vector<int> forced, early; size_t nearly = 0;
  for (auto& kv : g_cbs.b) { if (kv.second.oldest <= ck - g_cbk_w + 1) forced.push_back(kv.first); else if (kv.second.items.size() >= g_cbk) { early.push_back(kv.first); nearly += kv.second.items.size(); } }
  std::vector<int> ks;
  if (!forced.empty()) { ks = forced; ks.insert(ks.end(), early.begin(), early.end()); std::sort(ks.begin(), ks.end()); g_cbs.st_forced++; }
  else if (nearly * 2 >= g_cbk_cap) { ks = early; g_cbs.st_early++; }
  auto lane_build = [&]() { for (size_t i = 0; i < lane_pc.size(); i++) {
    FreetDeferred& N = *lane_pc[i]; std::vector<uint8_t> sel2(N.prs.size(), 0); tm_mark("cb: cofit_start (deferred)");
    N.cojob = cofit_start(N.prs, sel2, yvv, flrs, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, FJN_TMAX);
    if (sel2 != lane_sel[i]) { fprintf(stderr, "AGHC ASSERT: deferred co-fit selection differs from the routed one\n"); abort(); } } lane_pc.clear(); };
  if ((g_aghc & 6) == 2) cbk_emit(ks, lane_build); else cbk_emit(ks);
  if (g_aghc & 6) lane_build();
}
static void cbk_step_print() {
  if (!g_cbk_on || !g_cbs.st_calls) return;
  printf("CBUCKET step: B %zu W %ld cap %zu | calls %ld pairs %ld | classes per call mean %.1f max %ld | end-of-chunk flushes forced %ld early %ld | copy-back %.1f ms\n", g_cbk, g_cbk_w, g_cbk_cap, g_cbs.st_calls, g_cbs.st_pairs,
    (double)g_cbs.st_classes / (double)g_cbs.st_calls, g_cbs.st_maxcls, g_cbs.st_forced, g_cbs.st_early, g_cbs.st_copy_ms);
  g_cbs.st_calls = g_cbs.st_pairs = g_cbs.st_classes = g_cbs.st_maxcls = g_cbs.st_early = g_cbs.st_forced = 0; g_cbs.st_copy_ms = 0;
}
static bool aghc_piece_ready(FreetDeferred& D) {
  if (!D.any) return true;
  if (D.cbk) { if (D.npend) return false; for (auto& c : D.calls) if (!c->done && c->p.valid && c->p.cb.status < MTLCommandBufferStatusCompleted) return false; return !D.cojob || cofit_done(D.cojob); }
  if (D.tk) { std::lock_guard<std::mutex> lk(D.tk->m); return D.tk->done; }
  return !D.p.valid || (D.p.cb.status >= MTLCommandBufferStatusCompleted && cofit_done(D.p.cojob));
}
static bool aghc_fill() {
  if (!g_stream || g_holdk > 0 || g_stp.empty()) return false;
  const long keep = g_cbk_on ? g_st_chunk - g_cbk_w : g_st_chunk - 1;
  if (g_stp[0].chunk >= keep || !aghc_piece_ready(*g_stp[0].d)) return false;
  freet_drain(*g_stp[0].d); g_stp.erase(g_stp.begin());
  return !g_stp.empty() && g_stp[0].chunk < keep;
}
static struct AghcInit { AghcInit() { if (g_aghc & 8) g_afill_fn = aghc_fill; } } g_aghc_init;
static void freet_drain_prior() {
  if (g_stream) { if (g_hold.d && g_hold.chunk0 < g_st_chunk - 1) hold_commit(); stream_drain(g_st_chunk - 1); return; }
  freet_drain(g_fdq[g_fdi]);
}
static void freet_flush(Times& T) {
  (void)T; freet_drain(g_fdq[g_fdi]); freet_drain(g_fdq[g_fdi ^ 1]);
  hold_commit(); if (g_cbk_on) cbk_flush_all(); stream_drain(1L << 60); cbk_step_print();
}
static void freet_gap_replay(Engine& E) {
  static int nstep = 0; static const char* ev = getenv("FJN_GAP_REPLAY"); if (!ev) return;
  static const int maxst = getenv("FJN_GAP_STEPS") ? atoi(getenv("FJN_GAP_STEPS")) : 1;
  nstep++; if (nstep > maxst || g_gap_pairs.empty()) { g_gap_pairs.clear(); g_gap_piece = 0; return; }
  struct Cfg { std::string name; int K; bool sorted, eng, burn; }; std::vector<Cfg> cf;
  { std::string t(ev); size_t i = 0; while (i < t.size()) { size_t j = t.find(',', i); if (j == std::string::npos) j = t.size(); std::string tok = t.substr(i, j - i); i = j + 1; if (tok.empty()) continue;
      Cfg c{tok, 1, false, false, false}; if (tok.back() == 'b') { c.burn = true; tok.pop_back(); } if (tok == "E") c.eng = true; else { if (tok[0] == 'S') { c.sorted = true; tok.erase(0, 1); } c.K = std::max(1, atoi(tok.c_str())); } cf.push_back(c); } }
  const size_t N = g_gap_pairs.size(); std::vector<FrPairIn> all; std::vector<int> pc; all.reserve(N); for (auto& g : g_gap_pairs) { all.push_back(g->p); pc.push_back(g->piece); }
  std::vector<size_t> ord(N); for (size_t i = 0; i < N; i++) ord[i] = i;
  std::vector<size_t> srt = ord; std::stable_sort(srt.begin(), srt.end(), [&](size_t a, size_t b) { const int qa = std::max(all[a].x.Q, all[a].y.Q), qb = std::max(all[b].x.Q, all[b].y.Q), ta = std::max(all[a].x.nth, all[a].y.nth), tb = std::max(all[b].x.nth, all[b].y.nth); return qa != qb ? qa < qb : ta < tb; });
  float yvv[96]; std::memcpy(yvv, g_gap_yv, sizeof yvv); const float flr = g_gap_flr;
  g_cofit_bypass = true; g_gap_in_replay = true;
  std::map<std::string, std::vector<std::array<double, 4>>> res; std::map<std::string, int> npart;
  const int R = 3;
  for (int r = 0; r < R; r++) for (size_t ki = 0; ki < cf.size(); ki++) {
    const Cfg& c = cf[(r & 1) ? cf.size() - 1 - ki : ki];
    std::vector<std::vector<FrPairIn>> parts;
    if (c.eng) { int np = 0; for (int x : pc) np = std::max(np, x); parts.resize((size_t)np + 1); for (size_t i = 0; i < N; i++) parts[(size_t)pc[i]].push_back(all[i]); }
    else { parts.resize((size_t)c.K); const std::vector<size_t>& o = c.sorted ? srt : ord; for (int m = 0; m < c.K; m++) for (size_t i = N * (size_t)m / (size_t)c.K; i < N * (size_t)(m + 1) / (size_t)c.K; i++) parts[(size_t)m].push_back(all[o[i]]); }
    parts.erase(std::remove_if(parts.begin(), parts.end(), [](auto& v) { return v.empty(); }), parts.end());
    const int K = (int)parts.size(); npart[c.name] = K;
    std::vector<FjnPending> pend((size_t)K); std::vector<std::vector<FrOut>> outs((size_t)K); std::vector<FrRunStats> rs((size_t)K);
    std::shared_ptr<CoJob> burn; std::vector<FrPairIn> bp; std::vector<FrOut> bo; std::vector<uint8_t> bsel;
    if (c.burn && g_cofit_pm > 0) { for (int rep = 0; rep < 4; rep++) bp.insert(bp.end(), all.begin(), all.end()); for (size_t i = 0; i < bp.size(); i++) bp[i].id = (int)i; bo.assign(bp.size(), FrOut()); bsel.assign(bp.size(), 0); g_cofit_bypass = false; burn = cofit_start(bp, bsel, yvv, flr, g_fr_closed, g_fr_maxit, g_fr_warp, bo, FJN_TMAX); g_cofit_bypass = true; }
    const double t0 = gap_now();
    for (int m = 0; m < K; m++) { id<MTLCommandQueue> q = fjn_queue(E); fjn_run(E.dev, q, [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [q]() { return h182_cb(q); }, parts[(size_t)m], yvv, flr, g_fr_closed, g_fr_maxit, g_fr_warp, outs[(size_t)m], rs[(size_t)m], g_fjn_tpg, &pend[(size_t)m]); }
    double sum = 0, lo = 1e30, hi = 0, bins = 0;
    for (int m = 0; m < K; m++) { fjn_collect(pend[(size_t)m]); sum += (rs[(size_t)m].cb_end - rs[(size_t)m].cb_start) * 1e3; lo = std::min(lo, rs[(size_t)m].cb_start); hi = std::max(hi, rs[(size_t)m].cb_end); bins += rs[(size_t)m].bins; }
    const double wall = (gap_now() - t0) * 1e3;
    cofit_join(burn);
    res[c.name].push_back({sum, (hi - lo) * 1e3, wall, bins});
  }
  g_cofit_bypass = false; g_gap_in_replay = false;
  for (const Cfg& c : cf) { auto v = res[c.name]; std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a[1] < b[1]; }); auto& m = v[v.size() / 2];
    printf("GAP_REPLAY step %d cfg %-5s: pairs %zu CBs %d | sumCB %.1f ms span %.1f ms wall(enc+gpu) %.1f ms bins %.0f | span us/pair %.2f | rounds span:", nstep, c.name.c_str(), N, npart[c.name], m[0], m[1], m[2], m[3], m[1] * 1e3 / (double)N);
    for (auto& x : res[c.name]) printf(" %.1f", x[1]); printf("\n"); }
  fflush(stdout); g_gap_pairs.clear(); g_gap_piece = 0;
}
static void freet_chunk(Engine& E, PipeStep& PS, const Batch& Bt, float* err, float* rx, float* ry, Times& T) {
  TmScope tm_ft_("ft: glue_range(fixed-t)");
  const auto th0 = clk::now();
  StepBatch& SBr = *PS.SB; const Unique& UQ = SBr.UQ; const std::vector<uint8_t>& admit = *PS.admit;
  FrStats& S = g_frs; S.chunks++;
  if (!g_glue) glue_range(PS, Bt.b0, Bt.b1, err, rx, ry, T);
  tm_mark("ft: pair gather");
  std::vector<FrPairIn> prs; std::vector<int> child, ua, ub; prs.reserve((size_t)(Bt.b1 - Bt.b0));
  auto ref = [&](const Tree& t) { FrTreeRef r; r.Q = t.Q; r.nth = t.nth; r.nS = t.nS; r.nslot = t.nslot; r.prog = h313a::prog(t); r.nprog = h313a::nprog(t); r.sp = t.sp.data(); r.nsp = t.sp.size(); r.segs = h313a::segs(t); return r; };
  std::vector<uint8_t> gm; if (g_glue) gm.assign((size_t)(Bt.b1 - Bt.b0), 0);
  for (int i = Bt.b0; i < Bt.b1; i++) {
    if (!admit[i]) continue;
    const int a = SBr.ct[(size_t)i * 2], b = SBr.ct[(size_t)i * 2 + 1]; if (a < 0 || b < 0) continue;
    const Tree& ta = UQ.U[(size_t)a]; const Tree& tb = UQ.U[(size_t)b];
    if (!isfinite(SBr.F[0].rfin[(size_t)a]) || !isfinite(SBr.F[1].rfin[(size_t)b])) { S.b_nonfin++; continue; }
    S.pairs++;
    if ((ta.nth >= 1 && !SBr.F[0].is_g[(size_t)a]) || (tb.nth >= 1 && !SBr.F[1].is_g[(size_t)b])) { S.fb_nokernel++; S.b_g++; S.fb_fitted++; continue; }
    FrPairIn p; p.x = ref(ta); p.y = ref(tb); p.thx = &F14.th[0][32 * (size_t)a]; p.thy = &F14.th[1][32 * (size_t)b]; p.id = i; p.k1 = UQ.hash[(size_t)a]; p.k2 = UQ.hash[(size_t)b];
    prs.push_back(p); child.push_back(i); ua.push_back(a); ub.push_back(b);
    if (g_pprobe) { const uint64_t k = (UQ.hash[(size_t)a] * 0x9E3779B97F4A7C15ull) ^ (UQ.hash[(size_t)b] + 0x7F4A7C15ull + (UQ.hash[(size_t)a] << 6)); g_pp_all++; if (!g_pp_seen.insert(k).second) g_pp_hit++; }
    if (g_glue) gm[(size_t)(i - Bt.b0)] = 1;
  }
  if (g_glue) for (int i = Bt.b0; i < Bt.b1;) { if (gm[(size_t)(i - Bt.b0)]) { i++; continue; } int j = i; while (j < Bt.b1 && !gm[(size_t)(j - Bt.b0)]) j++; glue_range(PS, i, j, err, rx, ry, T); i = j; }
  if (prs.empty()) { S.host_ms += ms_since(th0); return; }
  tm_mark("ft: setup");
  g_kgp_par = [&E](int n, const std::function<void(int)>& f) { PL(E).run(n, f); };
  if (g_cofit_pm > 0) {
    if (g_cofit_scorer != g_fr_score96) g_cofit_scorer = g_fr_score96;
    if (!g_cofit_class_ok) g_cofit_class_ok = [](int Q, int T) { static std::unordered_set<std::string> nm; static bool in = false; if (!in) { in = true; for (NSString* n in g_fjnL.functionNames) nm.insert(n.UTF8String); } return nm.count("fjn_q" + std::to_string(Q) + "_t" + std::to_string(T)) > 0; };
  }
  g_wb_par = [&E](int n, const std::function<void(int)>& f) { PL(E).run(n, f); };
  float yvv[96]; std::memcpy(yvv, E.tg[0].yf, 48 * 4); std::memcpy(yvv + 48, E.tg[1].yf, 48 * 4);
  const float flrs = E.tg[0].floor_ + E.tg[1].floor_;
  std::vector<FrOut> outs; FrRunStats rs;
  if (g_gap_collect_on) { std::memcpy(g_gap_yv, yvv, sizeof yvv); g_gap_flr = flrs; }
  g_fit_cands = PS.cands; std::memcpy(g_fit_yv, yvv, sizeof g_fit_yv); fitx::g_yv1 = g_fit_yv;
  S.host_ms += ms_since(th0);
  if (g_stream) {
    const long ck = g_st_chunk++; if (g_hold.d && (g_hold.d->sbp != &SBr || g_hold.chunk0 < ck - 1 || std::memcmp(g_hold.yv, yvv, sizeof yvv) != 0 || g_hold.flr != flrs)) hold_commit();
    const bool defer_drain_ = (g_aghc & 1) && g_holdk == 0; const long keep_ = g_cbk_on ? ck - g_cbk_w : ck - 1;
    tm_mark("ft: stream_drain(ck-1)"); if (!defer_drain_) stream_drain(keep_); tm_mark("ft: stream piece split");
    const size_t NP_ = prs.size(); size_t NL_ = NP_;
    if (const char* tm_ = getenv("FJN_COFIT_TMIN")) { const int t_ = atoi(tm_); if (t_ > 0 && g_cofit_pm > 0) { NL_ = 0; for (const FrPairIn& q_ : prs) NL_ += std::max(q_.x.nth, q_.y.nth) < t_; } }
    const size_t NG_ = g_stream_auto ? (NL_ * (size_t)(1000 - g_cofit_pm) + 999) / 1000 : 0;
    const size_t M = std::min(g_stream_auto ? std::max<size_t>(1, (NG_ + g_stream_cap - 1) / g_stream_cap) : (size_t)g_stream_n, std::max<size_t>(1, NP_));
    if (g_cbk_on) { cbk_chunk(E, SBr, ck, prs, child, ua, ub, yvv, flrs, err, rx, ry, M); if (defer_drain_) { tm_mark("ft: deferred stream_drain"); stream_drain(keep_); } S.host_ms += ms_since(th0); return; }
    for (size_t m = 0; m < M; m++) {
      const size_t lo = NP_ * m / M, hi = NP_ * (m + 1) / M; if (lo == hi) continue;
      StreamPiece sp; sp.chunk = ck; sp.d.reset(new FreetDeferred()); FreetDeferred& N = *sp.d;
      N.any = true; N.prs.assign(prs.begin() + (long)lo, prs.begin() + (long)hi); N.child.assign(child.begin() + (long)lo, child.begin() + (long)hi);
      N.ua.assign(ua.begin() + (long)lo, ua.begin() + (long)hi); N.ub.assign(ub.begin() + (long)lo, ub.begin() + (long)hi);
      N.sbp = &SBr; N.err = err; N.rx = rx; N.ry = ry;
      const bool hold_here_ = g_holdk > 0 && !(g_hold_edge >= 1 && Bt.b1 >= PS.n) && !(g_hold_edge >= 2 && Bt.ci == 0);
      if (g_holdk > 0 && !hold_here_) hold_commit();
      if (hold_here_) {
        if (g_hold.d) { FreetDeferred& H = *g_hold.d; H.prs.insert(H.prs.end(), N.prs.begin(), N.prs.end()); H.child.insert(H.child.end(), N.child.begin(), N.child.end()); H.ua.insert(H.ua.end(), N.ua.begin(), N.ua.end()); H.ub.insert(H.ub.end(), N.ub.begin(), N.ub.end()); g_hold.chunk1 = ck; }
        else { g_hold.d = std::move(sp.d); g_hold.chunk0 = g_hold.chunk1 = ck; std::memcpy(g_hold.yv, yvv, sizeof yvv); g_hold.flr = flrs; g_hold.E = &E; }
        if (++g_hold.npiece >= g_holdk || (g_holdk_maxgp && g_hold.d->prs.size() * (size_t)(1000 - g_cofit_pm) / 1000 >= g_holdk_maxgp)) hold_commit();
        continue;
      }
      tm_mark("ft: fjn_run (enc+commit)");
      freet_fjn_submit(E, N, yvv, flrs);
      g_stp.push_back(std::move(sp)); tm_mark("ft: stream piece split");
    }
    if (defer_drain_) { tm_mark("ft: deferred stream_drain"); stream_drain(keep_); }
    S.host_ms += ms_since(th0);
    return;
  }
  if (g_freet_async) {
    FreetDeferred& N = g_fdq[g_fdi]; freet_drain(N); g_fdi ^= 1;
    N.any = true; N.prs = prs; N.child = child; N.ua = ua; N.ub = ub; N.sbp = &SBr; N.err = err; N.rx = rx; N.ry = ry;
    N.p = FjnPending();
    N.rs = FrRunStats();
    fjn_run(E.dev, fjn_queue(E), [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&]() { return h182_cb(fjn_queue(E)); }, N.prs, yvv, flrs, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, N.rs, g_fjn_tpg, &N.p);
    S.host_ms += ms_since(th0);
    return;
  }
  fjn_run(E.dev, fjn_queue(E), [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&]() { return h182_cb(fjn_queue(E)); }, prs, yvv, flrs, g_fr_closed, g_fr_maxit, g_fr_warp, outs, rs, g_fjn_tpg);
  if (rs.cb_end > 0) PS.iv.push_back({rs.cb_start, rs.cb_end});
  S.b_q += rs.b_q; S.b_tn += rs.b_tn; S.b_pso += rs.b_pso; S.b_slot += rs.b_slot; S.b_tpg += rs.b_tpg;
  S.fb_nokernel += rs.nokernel; S.fb_tg += rs.tg; S.fb_enc += rs.enc; S.bins += rs.bins; S.gpu_ms += rs.gpu_ms; S.host_ms += rs.host_ms; S.wait_ms += rs.wait_ms;
  const auto tr0 = clk::now();
  freet_writeback_outs(outs, child, ua, ub, prs, SBr, err, rx, ry, S, UQ, yvv, flrs);
  S.host_ms += ms_since(tr0);
}
