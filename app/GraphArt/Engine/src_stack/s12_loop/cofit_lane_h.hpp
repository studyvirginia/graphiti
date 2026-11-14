#pragma once
#include "cofit_cpu.hpp"
#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <sys/sysctl.h>
#include <pthread.h>
static int g_cofit_pm = getenv("FJN_COFIT") ? std::max(0, std::min(1000, atoi(getenv("FJN_COFIT")))) : 0;
static const int g_cofit_auto = getenv("FJN_COFIT_AUTO") ? atoi(getenv("FJN_COFIT_AUTO")) : 0;
static std::atomic<long> g_auto_lane_us{0}, g_auto_np{0}; static std::mutex g_auto_mx; static std::vector<std::pair<double, double>> g_auto_iv; static double g_auto_lane_t0 = 0;
static inline double auto_now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double (*g_cofit_scorer)(const float*) = nullptr;
static std::function<bool(int, int)> g_cofit_class_ok;
static void (*g_cofit_fin_hook)() = nullptr;
static std::atomic<long> g_cofit_score_us{0};
static std::atomic<long> g_cofit_pairs{0}, g_cofit_busy_us{0}; static std::atomic<long> g_cofit_wait_us{0}, g_cofit_encfail{0}; static std::atomic<long> g_cofit_jobs{0};
struct CoWk { std::atomic<long> pairs{0}, wall_us{0}, cpu_us{0}; std::atomic<int> ew{0}; }; static CoWk g_cofit_wk[64]; static std::atomic<int> g_cofit_nwk{0};
#include <time.h>
#include <sys/resource.h>
#include <mach/mach.h>
static inline void cofit_thrt(int st) {
  thread_act_array_t th; mach_msg_type_number_t n = 0;
  if (task_threads(mach_task_self(), &th, &n) == KERN_SUCCESS) {
    for (unsigned i = 0; i < n; i++) {
      thread_basic_info_data_t bi; mach_msg_type_number_t c = THREAD_BASIC_INFO_COUNT;
      thread_identifier_info_data_t ii; mach_msg_type_number_t c2 = THREAD_IDENTIFIER_INFO_COUNT;
      if (thread_info(th[i], THREAD_BASIC_INFO, (thread_info_t)&bi, &c) == KERN_SUCCESS && thread_info(th[i], THREAD_IDENTIFIER_INFO, (thread_info_t)&ii, &c2) == KERN_SUCCESS) {
        const double u = bi.user_time.seconds * 1e3 + bi.user_time.microseconds / 1e3, y = bi.system_time.seconds * 1e3 + bi.system_time.microseconds / 1e3;
        if (u + y >= 0.5) printf("THRT call %d: tid %llu user %.1f sys %.1f ms\n", st, (unsigned long long)ii.thread_id, u, y);
      }
      mach_port_deallocate(mach_task_self(), th[i]);
    }
    vm_deallocate(mach_task_self(), (vm_address_t)th, n * sizeof(thread_act_t));
  }
  mach_task_basic_info_data_t tb; mach_msg_type_number_t c3 = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&tb, &c3) == KERN_SUCCESS)
    printf("THRT call %d: dead user %.1f sys %.1f ms\n", st, tb.user_time.seconds * 1e3 + tb.user_time.microseconds / 1e3, tb.system_time.seconds * 1e3 + tb.system_time.microseconds / 1e3);
  task_events_info_data_t te; mach_msg_type_number_t c4 = TASK_EVENTS_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_EVENTS_INFO, (task_info_t)&te, &c4) == KERN_SUCCESS)
    printf("TEVT call %d: faults %d cow %d mach_sys %d unix_sys %d csw %d\n", st, te.faults, te.cow_faults, te.syscalls_mach, te.syscalls_unix, te.csw);
}
static inline void cofit_cput(double* us) { struct rusage r; getrusage(RUSAGE_SELF, &r); us[0] = r.ru_utime.tv_sec * 1e3 + r.ru_utime.tv_usec / 1e3; us[1] = r.ru_stime.tv_sec * 1e3 + r.ru_stime.tv_usec / 1e3; }
static const bool g_bcount = getenv("BCOUNT") && atoi(getenv("BCOUNT")) == 1;
static std::atomic<long> g_bc_builds{0}, g_bc_bns{0}, g_bc_sns{0}, g_bc_swall{0}, g_bc_calls{0}, g_bc_nodes{0}; static std::mutex g_bc_m; static std::unordered_set<uint64_t> g_bc_u;
static inline long bc_tcpu() { return (long)clock_gettime_nsec_np(CLOCK_THREAD_CPUTIME_ID); }
static inline uint64_t bc_hash(int Q, int nth, int nS, int nslot, const uint32_t* p, size_t np, const uint16_t* sp, size_t nsp, const uint16_t* sg) {
  uint64_t h = 1469598103934665603ull; auto mix = [&](const void* d, size_t n) { const unsigned char* c = (const unsigned char*)d; for (size_t i = 0; i < n; i++) { h ^= c[i]; h *= 1099511628211ull; } };
  int hd[6] = {Q, nth, nS, nslot, (int)np, (int)nsp}; mix(hd, sizeof hd); if (np) mix(p, np * 4); if (nsp) mix(sp, nsp * 2); mix(sg, 72); return h; }
static std::atomic<long> g_cofit_ph[16];
struct CoNode { std::vector<uint32_t> prog; std::vector<uint16_t> sp; uint16_t segs[36]; int Q = 0, nth = 0, nS = 0, nslot = 0; cofit::Tree T; std::once_flag once; bool ok = false; };
struct CoPair { uint32_t k; int nx, ny; int a, b; int Qc; float thx[32], thy[32]; };
struct CoJob {
  std::vector<std::unique_ptr<CoNode>> nodes; std::vector<CoPair> prs; std::vector<FrOut>* out = nullptr;
  float yv[96]; float flr = 0; bool closed = false; int maxit = 60; const fjn::WarpMat* WM = nullptr; double (*scorer)(const float*) = nullptr;
  std::atomic<size_t> next{0}, done{0}; std::mutex m; std::condition_variable cv; bool fin = false; double done_t = 0;
};
struct CoPool {
  std::mutex m; std::condition_variable cv; std::deque<std::shared_ptr<CoJob>> q; std::vector<std::thread> th; bool started = false;
  static int nthreads() {
    if (const char* e = getenv("FJN_COFIT_THREADS")) return std::max(1, atoi(e));
    int v = 0; size_t sz = sizeof v; if (sysctlbyname("hw.perflevel0.physicalcpu", &v, &sz, nullptr, 0) != 0 || v < 1) { sz = sizeof v; if (sysctlbyname("hw.physicalcpu", &v, &sz, nullptr, 0) != 0 || v < 2) v = 2; else v = std::max(1, v / 2); }
    return v; }
  static int ethreads() {
    int v = 0; if (const char* e = getenv("FJN_COFIT_ETHREADS")) v = std::max(0, atoi(e)); else return 0;
    int ne = 0; size_t sz = sizeof ne; if (sysctlbyname("hw.perflevel1.physicalcpu", &ne, &sz, nullptr, 0) != 0) ne = 0; return std::min(v, ne); }
  static void run_one(CoJob& J, size_t i, cofit::Scratch& S) {
    const CoPair& P = J.prs[i]; CoNode& X = *J.nodes[(size_t)P.a]; CoNode& Y = *J.nodes[(size_t)P.b];
    for (CoNode* N : {&X, &Y}) std::call_once(N->once, [N] { const long c0 = g_bcount ? bc_tcpu() : 0; cofit::TreeIn in{N->Q, N->nth, N->nS, N->nslot, N->prog.data(), N->prog.size(), N->sp.data(), N->sp.size(), N->segs}; N->ok = cofit::build_tree(N->T, in);
      if (g_bcount) { g_bc_bns += bc_tcpu() - c0; g_bc_builds++; const uint64_t h = bc_hash(N->Q, N->nth, N->nS, N->nslot, N->prog.data(), N->prog.size(), N->sp.data(), N->sp.size(), N->segs); std::lock_guard<std::mutex> lk(g_bc_m); g_bc_u.insert(h); } });
    FrOut& o = (*J.out)[P.k];
    if (!X.ok || !Y.ok) { g_cofit_encfail++; return; }
    static thread_local cofit::Out r; cofit::Ctx cx{J.yv, J.flr, J.closed, J.maxit, J.WM};
    cofit::fit_pair(cx, X.T, Y.T, P.thx, P.thy, P.Qc, S, r);
    o.Fx = r.Fx; o.Fy = r.Fy; o.stalled = r.stalled; o.it = r.it; o.trials = r.trials; o.dual = r.dual; o.plain = r.plain;
    std::memcpy(o.th[0], r.th[0], sizeof o.th); std::memcpy(o.a[0], r.a[0], sizeof o.a); std::memcpy(o.u, r.u, sizeof o.u); std::memcpy(o.traj, r.traj, sizeof o.traj); std::memcpy(o.res48, r.res48, sizeof o.res48);
    o.raw[0] = r.Fx; o.raw[1] = r.Fy; o.raw[2] = r.stalled; o.raw[3] = r.it; o.raw[4] = r.trials; o.raw[5] = r.dual; o.raw[6] = r.plain; o.raw[7] = r.solves;
    g_cofit_ph[8] += (long)r.plain; g_cofit_ph[9] += (long)r.dual; g_cofit_ph[10] += (long)r.solves; g_cofit_ph[11] += (long)r.it;
    if (std::isfinite(r.Fx) && std::isfinite(r.Fy) && J.scorer) { const auto ts0 = std::chrono::steady_clock::now(); const double s = J.scorer(r.res48); g_cofit_score_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - ts0).count(); o.cf = std::isfinite(s) ? std::min((float)s, 1.0e30f) : 1.0e30f; } else o.cf = 1.0e30f;
    o.status = 1;
  }
  static int g_cofit_eq() { static const int v = getenv("FJN_COFIT_EQOS") ? atoi(getenv("FJN_COFIT_EQOS")) : 0; return v; }
  void worker(bool eworker) {
    const int wid = g_cofit_nwk.fetch_add(1) & 63; g_cofit_wk[wid].ew = eworker ? 1 : 0;
    { static const int pq = getenv("FJN_COFIT_PQOS") ? atoi(getenv("FJN_COFIT_PQOS")) : 0;
      pthread_set_qos_class_self_np(eworker ? (g_cofit_eq() == 3 ? QOS_CLASS_USER_INITIATED : g_cofit_eq() == 2 ? QOS_CLASS_DEFAULT : g_cofit_eq() == 1 ? QOS_CLASS_UTILITY : QOS_CLASS_BACKGROUND) : (pq == 2 ? QOS_CLASS_UTILITY : pq == 1 ? QOS_CLASS_DEFAULT : QOS_CLASS_USER_INITIATED), 0); }
    cofit::set_ftz(); cofit::Scratch S;
    for (;;) {
      std::shared_ptr<CoJob> J;
      { std::unique_lock<std::mutex> lk(m); cv.wait(lk, [&] { return !q.empty(); }); J = q.front(); }
      const size_t n = J->prs.size();
      for (;;) {
        const size_t i0 = J->next.fetch_add(4); if (i0 >= n) break; const size_t i1 = std::min(n, i0 + 4);
        const uint64_t c0 = clock_gettime_nsec_np(CLOCK_THREAD_CPUTIME_ID);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = i0; i < i1; i++) run_one(*J, i, S);
#ifdef COFIT_PROF
        for (int q = 0; q < 8; q++) { g_cofit_ph[q] += (long)g_ptn[q]; g_ptn[q] = 0; } for (int q = 0; q < 2; q++) { g_cofit_ph[12 + q] += g_cnt[q]; g_cnt[q] = 0; }
#endif
        { CoWk& W = g_cofit_wk[wid]; W.pairs += (long)(i1 - i0); W.cpu_us += (long)((clock_gettime_nsec_np(CLOCK_THREAD_CPUTIME_ID) - c0) / 1000); W.wall_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count(); }
        g_cofit_busy_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
        if (J->done.fetch_add(i1 - i0) + (i1 - i0) == n) { { std::lock_guard<std::mutex> lk(J->m); J->done_t = gap_now(); J->fin = true; J->cv.notify_all(); } if (g_cofit_fin_hook) g_cofit_fin_hook(); }
      }
      { std::lock_guard<std::mutex> lk(m); if (!q.empty() && q.front() == J) { q.pop_front(); if (g_cofit_auto && q.empty() && g_auto_lane_t0 > 0) { g_auto_lane_us += (long)((auto_now() - g_auto_lane_t0) * 1e6); g_auto_lane_t0 = 0; } } }
    }
  }
  void submit(std::shared_ptr<CoJob> J) {
    std::lock_guard<std::mutex> lk(m);
    if (!started) { started = true; const int n = nthreads(), ne = ethreads(); for (int i = 0; i < n + ne; i++) { const bool ew = i >= n; th.emplace_back([this, ew] { worker(ew); }); th.back().detach(); } }
    if (g_cofit_auto && q.empty()) g_auto_lane_t0 = auto_now(); q.push_back(std::move(J)); cv.notify_all();
  }
};
static CoPool& cofit_pool() { static CoPool* p = new CoPool(); return *p; }
static inline bool cofit_pick(size_t k, int pm) { return ((k + 1) * (size_t)pm) / 1000 != (k * (size_t)pm) / 1000; }
static void cofit_join(const std::shared_ptr<CoJob>& J) {
  if (!J) return; const auto t0 = std::chrono::steady_clock::now();
  { std::unique_lock<std::mutex> lk(J->m); J->cv.wait(lk, [&] { return J->fin; }); }
  g_cofit_wait_us += (long)std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
}
static bool cofit_done(const std::shared_ptr<CoJob>& J) { if (!J) return true; std::lock_guard<std::mutex> lk(J->m); return J->fin; }
static bool cofit_eligible(const FrPairIn& p, size_t k, int fjn_tmax) {
  static const int tmin_ = getenv("FJN_COFIT_TMIN") ? atoi(getenv("FJN_COFIT_TMIN")) : 0; const int Qc = std::max(p.x.Q, p.y.Q), Tc = std::max(p.x.nth, p.y.nth);
  if (!((tmin_ > 0 && Tc >= tmin_) || cofit_pick(p.id >= 0 ? (size_t)p.id : k, g_cofit_pm))) return false;
  if (Qc < 1 || Qc > 16 || Tc > fjn_tmax || p.x.nth > 14 || p.y.nth > 14 || (p.x.nth > 0 && !p.thx) || (p.y.nth > 0 && !p.thy)) return false;
  if (g_cofit_class_ok && !g_cofit_class_ok(Qc, Tc)) return false;
  return true; }
static void cofit_select(const std::vector<FrPairIn>& pairs, std::vector<uint8_t>& sel, int fjn_tmax) { if (g_cofit_pm <= 0) return; for (size_t k = 0; k < pairs.size(); k++) sel[k] = cofit_eligible(pairs[k], k, fjn_tmax) ? 1 : 0; }
static std::shared_ptr<CoJob> cofit_start(const std::vector<FrPairIn>& pairs, std::vector<uint8_t>& sel, const float* yv96, float flr, bool closed, int maxit, const fjn::WarpMat& WM, std::vector<FrOut>& out, int fjn_tmax) {
  if (g_cofit_pm <= 0) return nullptr;
  struct BcT { long c0 = 0; std::chrono::steady_clock::time_point w0; std::shared_ptr<CoJob>* J = nullptr; BcT() { if (g_bcount) { c0 = bc_tcpu(); w0 = std::chrono::steady_clock::now(); } }
    ~BcT() { if (g_bcount) { g_bc_sns += bc_tcpu() - c0; g_bc_swall += (long)std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - w0).count(); g_bc_calls++; if (J && *J) g_bc_nodes += (long)(*J)->nodes.size(); } } } bct_;
  auto J = std::make_shared<CoJob>(); bct_.J = &J; std::unordered_map<uint64_t, int> nmap; nmap.reserve(256);
  auto node = [&](const FrTreeRef& t) -> int {
    const uint64_t key = ((uint64_t)(uintptr_t)t.prog * 1000003ull) ^ ((uint64_t)(uintptr_t)t.sp * 0x9E3779B97F4A7C15ull) ^ ((uint64_t)(uintptr_t)t.segs * 31ull) ^ ((uint64_t)t.nprog << 48) ^ ((uint64_t)t.nsp << 40);
    auto it = nmap.find(key);
    if (it != nmap.end()) { const CoNode& c = *J->nodes[(size_t)it->second]; if (c.Q == t.Q && c.nth == t.nth && c.nS == t.nS && c.nslot == t.nslot && c.prog.size() == t.nprog && c.sp.size() == t.nsp && (t.nprog == 0 || !std::memcmp(c.prog.data(), t.prog, t.nprog * 4)) && !std::memcmp(c.segs, t.segs, 72) && (t.nsp == 0 || !std::memcmp(c.sp.data(), t.sp, t.nsp * 2))) return it->second; }
    auto N = std::make_unique<CoNode>(); N->prog.assign(t.prog, t.prog + t.nprog); N->sp.assign(t.sp, t.sp + t.nsp); std::memcpy(N->segs, t.segs, 72); N->Q = t.Q; N->nth = t.nth; N->nS = t.nS; N->nslot = t.nslot;
    J->nodes.push_back(std::move(N)); const int id = (int)J->nodes.size() - 1; nmap[key] = id; return id; };
  for (size_t k = 0; k < pairs.size(); k++) {
    static const int tmin_ = getenv("FJN_COFIT_TMIN") ? atoi(getenv("FJN_COFIT_TMIN")) : 0;
    sel[k] = 0; const FrPairIn& p = pairs[k]; const int Qc = std::max(p.x.Q, p.y.Q), Tc = std::max(p.x.nth, p.y.nth);
    if (!((tmin_ > 0 && Tc >= tmin_) || cofit_pick(pairs[k].id >= 0 ? (size_t)pairs[k].id : k, g_cofit_pm))) continue;
    if (Qc < 1 || Qc > 16 || Tc > fjn_tmax || p.x.nth > 14 || p.y.nth > 14 || (p.x.nth > 0 && !p.thx) || (p.y.nth > 0 && !p.thy)) continue;
    if (g_cofit_class_ok && !g_cofit_class_ok(Qc, Tc)) continue;
    CoPair c; c.k = (uint32_t)k; c.nx = p.x.nth; c.ny = p.y.nth; c.Qc = Qc; std::memset(c.thx, 0, sizeof c.thx); std::memset(c.thy, 0, sizeof c.thy);
    std::memcpy(c.thx, p.thx, (size_t)p.x.nth * 4); std::memcpy(c.thy, p.thy, (size_t)p.y.nth * 4);
    c.a = node(p.x); c.b = node(p.y); J->prs.push_back(c); sel[k] = 1;
  }
  if (J->prs.empty()) return nullptr;
  std::memcpy(J->yv, yv96, 96 * 4); J->flr = flr; J->closed = closed; J->maxit = maxit; J->WM = &WM; J->out = &out; J->scorer = g_cofit_scorer;
  g_cofit_pairs += (long)J->prs.size(); g_cofit_jobs++;
  cofit_pool().submit(J);
  return J;
}
