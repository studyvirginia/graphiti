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
static void freet_step_print(int st) {
  if (g_cfchk > 0 && g_cfn + g_cfdead > 0) { printf("CISFUSE_CHECK step %d: checked %ld | max|gpu-cpu| %.6g | max rel (|d|/max(|cpu|,1)) %.6g | dead-mismatch %ld | worst gpu %.9g cpu %.9g\n", st, g_cfn, g_cfmaxabs, g_cfmaxrel, g_cfdead, g_cfwg, g_cfwc); g_cfn = g_cfdead = 0; g_cfmaxabs = g_cfmaxrel = 0; }
  FrStats& s = g_frs;
  if (g_pprobe && g_pp_all) { printf("PAIRPROBE step %d: pairs %ld already fitted earlier in the run %ld (%.1f%%) | seen %zu\n", st, g_pp_all, g_pp_hit, 100.0 * (double)g_pp_hit / (double)g_pp_all, g_pp_seen.size()); g_pp_all = g_pp_hit = 0; }
  if (g_freet && g_fr_verbose && s.chunks) printf("H950 step %d: chunks %ld pairs %ld dispatched %ld bins %ld | fallback nokernel %ld tg %ld enc %ld nonfinite %ld | LM iters mean %.1f dual %.1f plain %.1f stalled %ld | joint worse than start %ld | gpu %.2f ms host %.2f ms wait %.2f ms\n",
    st, s.chunks, s.pairs, s.disp, s.bins, s.fb_nokernel, s.fb_tg, s.fb_enc, s.fb_fit, s.disp ? (double)s.iters / s.disp : 0.0, s.disp ? (double)s.dual / s.disp : 0.0, s.disp ? (double)s.plain / s.disp : 0.0, s.stalled, s.worse, s.gpu_ms, s.host_ms, s.wait_ms);
  if (g_freet && g_fr_verbose && s.chunks) printf("H950 step %d branches: is_g %ld Q %ld Tn>max %ld pso %ld slot>32K %ld tpg<1 %ld nonfinite %ld\n", st, s.b_g, s.b_q, s.b_tn, s.b_pso, s.b_slot, s.b_tpg, s.b_nonfin);
  if (g_freet && s.chunks) {
    const long skipped = s.pairs - s.disp - s.fb_fitted;
    if (g_fr_verbose) printf("H950 step %d pairs check: pairs %ld = dispatched %ld + fallback_fitted %ld (skipped %ld)\n", st, s.pairs, s.disp, s.fb_fitted, skipped);
    if (skipped != 0) { fprintf(stderr, "H950 ASSERT step %d: pairs %ld != dispatched %ld + fallback_fitted %ld\n", st, s.pairs, s.disp, s.fb_fitted); abort(); }
  }
  g_fr.fb_fitted += s.fb_fitted; g_fr.b_g += s.b_g; g_fr.b_q += s.b_q; g_fr.b_tn += s.b_tn; g_fr.b_pso += s.b_pso; g_fr.b_slot += s.b_slot; g_fr.b_tpg += s.b_tpg; g_fr.b_nonfin += s.b_nonfin;
  g_fr.chunks += s.chunks; g_fr.pairs += s.pairs; g_fr.disp += s.disp; g_fr.fb_nokernel += s.fb_nokernel; g_fr.fb_tg += s.fb_tg; g_fr.fb_enc += s.fb_enc; g_fr.fb_fit += s.fb_fit; g_fr.bins += s.bins; g_fr.iters += s.iters; g_fr.dual += s.dual; g_fr.plain += s.plain; g_fr.stalled += s.stalled; g_fr.worse += s.worse; g_fr.gpu_ms += s.gpu_ms; g_fr.host_ms += s.host_ms; g_fr.wait_ms += s.wait_ms;
  s = FrStats();
}

static bool g_freet_async = getenv("FREET_ASYNC") && atoi(getenv("FREET_ASYNC")) != 0;
struct FreetDeferred { FjnPending p; FrRunStats rs; std::vector<FrOut> outs; std::vector<FrPairIn> prs; std::vector<int> child, ua, ub; const StepBatch* sbp = nullptr; float* err = nullptr; float* rx = nullptr; float* ry = nullptr; bool any = false; };
static FreetDeferred g_fdq[2]; static int g_fdi = 0;
static void freet_writeback_outs(const std::vector<FrOut>& outs, const std::vector<int>& child, const std::vector<int>& ua, const std::vector<int>& ub, const std::vector<FrPairIn>& prs,
                                 const StepBatch& SBr, float* err, float* rx, float* ry, FrStats& S, const Unique& UQ, const float* yvv, float flrs) {
  auto body = [&](size_t k, FrStats& A) {
    const FrOut& o = outs[k]; const int c = child[k];
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
  const auto fx_t0 = std::chrono::steady_clock::now();
  struct FxT { std::chrono::steady_clock::time_point t; ~FxT() { fitx::g_ms_wb += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); } } fx_tm{fx_t0};
  if (fitx::on() && g_fit_cands) {
    double tw[48]; fitx::Rec r; const size_t NK = prs.size();
    std::vector<uint64_t> kk(NK), hx(NK), hy(NK);
    const float bg = fitx::batch_gate(NK, [&](size_t k) { const FrOut& o = outs[k]; return (o.status && isfinite(o.Fx) && isfinite(o.Fy)) ? err[child[k]] : INFINITY; }, 1);
    auto prep = [&](size_t k) { const FrOut& o = outs[k]; const int c = child[k]; if (!o.status || !isfinite(o.Fx) || !isfinite(o.Fy) || !isfinite(err[c]) || err[c] > bg) return;
      kk[k] = fitx::gkey(*g_fit_cands[c]); hx[k] = fitx::shash(prs[k].x); hy[k] = fitx::shash(prs[k].y); };
    if (g_wb_par_on && g_wb_par && NK >= 256) g_wb_par(8, [&](int rr) { for (size_t k = NK * (size_t)rr / 8; k < NK * (size_t)(rr + 1) / 8; k++) prep(k); });
    else for (size_t k = 0; k < NK; k++) prep(k);
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
static void freet_writeback(FreetDeferred& D, FrStats& S) {
  const auto tr0 = clk::now();
  static FrStats* unused = nullptr; (void)unused;
  freet_writeback_outs(D.outs, D.child, D.ua, D.ub, D.prs, *D.sbp, D.err, D.rx, D.ry, S, D.sbp->UQ, nullptr, 0.0f);
  S.host_ms += ms_since(tr0);
}
static void freet_drain(FreetDeferred& D) {
  if (!D.any) return;
  FrStats& S = g_frs; const FrRunStats& r = D.rs;
  fjn_collect(D.p);
  S.fb_nokernel += r.nokernel; S.fb_tg += r.tg; S.fb_enc += r.enc; S.bins += r.bins; S.gpu_ms += r.gpu_ms; S.host_ms += r.host_ms; S.wait_ms += r.wait_ms;
  S.b_q += r.b_q; S.b_tn += r.b_tn; S.b_pso += r.b_pso; S.b_slot += r.b_slot; S.b_tpg += r.b_tpg;
  freet_writeback(D, S); D.any = false;
}
static bool g_stream = getenv("STREAM") && atoi(getenv("STREAM")) != 0 && g_freet_async;
static int g_stream_n = getenv("STREAM_N") ? std::max(1, atoi(getenv("STREAM_N"))) : 4;
struct StreamPiece { std::unique_ptr<FreetDeferred> d; long chunk; };
static std::vector<StreamPiece> g_stp; static long g_st_chunk = 0;
static void stream_drain(long keep_from) {
  size_t k = 0; while (k < g_stp.size() && g_stp[k].chunk < keep_from) { freet_drain(*g_stp[k].d); k++; }
  g_stp.erase(g_stp.begin(), g_stp.begin() + (long)k);
}
static void freet_drain_prior() {
  if (g_stream) { stream_drain(g_st_chunk - 1); return; }
  freet_drain(g_fdq[g_fdi]);
}
static void freet_flush(Times& T) {
  (void)T; freet_drain(g_fdq[g_fdi]); freet_drain(g_fdq[g_fdi ^ 1]);
  stream_drain(1L << 60);
}

static void freet_chunk(Engine& E, PipeStep& PS, const Batch& Bt, float* err, float* rx, float* ry, Times& T) {
  const auto th0 = clk::now();
  StepBatch& SBr = *PS.SB; const Unique& UQ = SBr.UQ; const std::vector<uint8_t>& admit = *PS.admit;
  FrStats& S = g_frs; S.chunks++;
  if (!g_glue) glue_range(PS, Bt.b0, Bt.b1, err, rx, ry, T);
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
    FrPairIn p; p.x = ref(ta); p.y = ref(tb); p.thx = &F14.th[0][32 * (size_t)a]; p.thy = &F14.th[1][32 * (size_t)b];
    prs.push_back(p); child.push_back(i); ua.push_back(a); ub.push_back(b);
    if (g_pprobe) { const uint64_t k = (UQ.hash[(size_t)a] * 0x9E3779B97F4A7C15ull) ^ (UQ.hash[(size_t)b] + 0x7F4A7C15ull + (UQ.hash[(size_t)a] << 6)); g_pp_all++; if (!g_pp_seen.insert(k).second) g_pp_hit++; }
    if (g_glue) gm[(size_t)(i - Bt.b0)] = 1;
  }
  if (g_glue) for (int i = Bt.b0; i < Bt.b1;) { if (gm[(size_t)(i - Bt.b0)]) { i++; continue; } int j = i; while (j < Bt.b1 && !gm[(size_t)(j - Bt.b0)]) j++; glue_range(PS, i, j, err, rx, ry, T); i = j; }
  if (prs.empty()) { S.host_ms += ms_since(th0); return; }
  g_kgp_par = [&E](int n, const std::function<void(int)>& f) { PL(E).run(n, f); };
  g_wb_par = [&E](int n, const std::function<void(int)>& f) { PL(E).run(n, f); };
  float yvv[96]; std::memcpy(yvv, E.tg[0].yf, 48 * 4); std::memcpy(yvv + 48, E.tg[1].yf, 48 * 4);
  const float flrs = E.tg[0].floor_ + E.tg[1].floor_;
  std::vector<FrOut> outs; FrRunStats rs;
  g_fit_cands = PS.cands; std::memcpy(g_fit_yv, yvv, sizeof g_fit_yv); fitx::g_yv1 = g_fit_yv;
  S.host_ms += ms_since(th0);
  if (g_stream) {
    const long ck = g_st_chunk++; stream_drain(ck - 1);
    const size_t NP_ = prs.size(), M = std::min((size_t)g_stream_n, NP_);
    for (size_t m = 0; m < M; m++) {
      const size_t lo = NP_ * m / M, hi = NP_ * (m + 1) / M; if (lo == hi) continue;
      StreamPiece sp; sp.chunk = ck; sp.d.reset(new FreetDeferred()); FreetDeferred& N = *sp.d;
      N.any = true; N.prs.assign(prs.begin() + (long)lo, prs.begin() + (long)hi); N.child.assign(child.begin() + (long)lo, child.begin() + (long)hi);
      N.ua.assign(ua.begin() + (long)lo, ua.begin() + (long)hi); N.ub.assign(ub.begin() + (long)lo, ub.begin() + (long)hi);
      N.sbp = &SBr; N.err = err; N.rx = rx; N.ry = ry;
      fjn_run(E.dev, fjn_queue(E), [&](const std::string& nm) { return E.pipe(g_fjnL, nm); }, [&]() { return h182_cb(fjn_queue(E)); }, N.prs, yvv, flrs, g_fr_closed, g_fr_maxit, g_fr_warp, N.outs, N.rs, g_fjn_tpg, &N.p);
      g_stp.push_back(std::move(sp));
    }
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
