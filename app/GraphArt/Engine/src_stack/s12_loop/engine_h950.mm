


#include "app_env_early.hpp"
#include "s1_twin_h950.hpp"
#import <Metal/Metal.h>
#include <chrono>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <functional>
#include <memory>
static thread_local bool tl_in_fitw = false;
static int g_fw_persist = getenv("H242D_PERSIST") ? atoi(getenv("H242D_PERSIST")) : 0;
static bool fw_persist() { return g_fw_persist != 0; }
#include <pthread/qos.h>
#include <sys/mman.h>
#include "memdiag.hpp"
#include "memdiag_mtl.hpp"
#define mmap(...) memdiag::mmap_w(__VA_ARGS__)
#define munmap(...) memdiag::munmap_w(__VA_ARGS__)
static const int g_h298q = getenv("H298Q") ? atoi(getenv("H298Q")) : 0;
static inline void h298_util(int bit) { if (g_h298q & bit) pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0); }
static const int g_kqos = getenv("K_QOS") ? atoi(getenv("K_QOS")) : 0;
static inline void k_qos_set(int bit) { if (g_kqos & bit) pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0); }
struct FitJob242 { std::mutex m; std::condition_variable cv; bool done = false; };
struct FitWorkers242 {
  std::mutex m; std::condition_variable cv;
  std::deque<std::pair<std::function<void()>, std::shared_ptr<FitJob242>>> q; bool started = false;
  void start() {
    const int W = getenv("H242D_W") ? std::max(1, atoi(getenv("H242D_W"))) : 4;
    for (int i = 0; i < W; i++) std::thread([this] { loop(); }).detach();
    started = true;
  }
  void loop() {
    tl_in_fitw = true;
    k_qos_set(1);
    for (;;) {
      std::pair<std::function<void()>, std::shared_ptr<FitJob242>> it;
      { std::unique_lock<std::mutex> lk(m); cv.wait(lk, [&] { return !q.empty(); }); it = std::move(q.front()); q.pop_front(); }
      it.first();
      { std::lock_guard<std::mutex> lk(it.second->m); it.second->done = true; }
      it.second->cv.notify_all();
    }
  }
};
static FitWorkers242* g_fw = new FitWorkers242();
struct FitHandle {
  std::shared_ptr<FitJob242> j; std::shared_ptr<std::thread> t;
  bool joinable() const { return j != nullptr || t != nullptr; }
  void join() {
    if (j) { std::unique_lock<std::mutex> lk(j->m); j->cv.wait(lk, [&] { return j->done; }); lk.unlock(); j.reset(); }
    if (t) { t->join(); t.reset(); }
  }
};

struct CpuSpan265 { std::mutex m; std::chrono::steady_clock::time_point s, e; bool any = false; };
static CpuSpan265* g_fit_span = nullptr;
static void h401_fitjob(const std::function<void()>& f);

#include <atomic>
static const int g_h873 = getenv("H873") ? atoi(getenv("H873")) : 0;
struct H873Team {
  struct W { std::atomic<int> st{0}; std::mutex m; std::condition_variable cv; std::function<void()> fn, fin; };
  struct Grp { std::mutex m; std::condition_variable cv; int left = 0; };
  std::mutex gm; std::vector<W*> idle[4];
  W* take(int cls) {
    { std::lock_guard<std::mutex> g(gm); if (!idle[cls].empty()) { W* w = idle[cls].back(); idle[cls].pop_back(); return w; } }
    W* w = new W; std::thread([this, w, cls] { loop(w, cls); }).detach(); return w;
  }
  void kick(W* w) { { std::lock_guard<std::mutex> lk(w->m); w->st.store(1, std::memory_order_release); } w->cv.notify_one(); }
  void loop(W* w, int cls) {
    for (;;) {
      const auto t0 = std::chrono::steady_clock::now(); int sp = 0;
      while (w->st.load(std::memory_order_acquire) != 1) { __builtin_arm_yield();
        if ((++sp & 255) == 0 && std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() > 20.0) break; }
      if (w->st.load(std::memory_order_acquire) != 1) { std::unique_lock<std::mutex> lk(w->m); w->cv.wait(lk, [&] { return w->st.load(std::memory_order_acquire) == 1; }); }
      std::function<void()> f = std::move(w->fn), fin = std::move(w->fin);
      f(); f = nullptr;
      w->st.store(0, std::memory_order_release);
      { std::lock_guard<std::mutex> g(gm); idle[cls].push_back(w); }
      fin();
    }
  }
  template <class F> void run(int cls, int nt, F& body) {
    Grp g; g.left = nt;
    for (int k = 0; k < nt; k++) { W* w = take(cls); w->fn = [&body] { body(); };
      w->fin = [&g] { std::lock_guard<std::mutex> lk(g.m); if (--g.left == 0) g.cv.notify_all(); }; kick(w); }
    std::unique_lock<std::mutex> lk(g.m); g.cv.wait(lk, [&] { return g.left == 0; });
  }
  void submit(int cls, std::function<void()> f, std::function<void()> fin) { W* w = take(cls); w->fn = std::move(f); w->fin = std::move(fin); kick(w); }
};
static H873Team g_h873t;
static FitHandle fw_submit(std::function<void()> f0) {
  std::function<void()> f = f0;
  if (CpuSpan265* R = g_fit_span) f = [f0, R]() { const auto a = std::chrono::steady_clock::now(); f0(); const auto b = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lk(R->m); if (!R->any || a < R->s) R->s = a; if (!R->any || b > R->e) R->e = b; R->any = true; };
  FitHandle h;
  if (!fw_persist() && (g_h873 & 2)) { h.j = std::make_shared<FitJob242>(); auto jj = h.j;
    g_h873t.submit(3, [f = std::move(f)] { h401_fitjob(f); }, [jj] { { std::lock_guard<std::mutex> lk(jj->m); jj->done = true; } jj->cv.notify_all(); }); return h; }
  if (!fw_persist()) { h.t = std::make_shared<std::thread>([f = std::move(f)] { h401_fitjob(f); }); return h; }
  h.j = std::make_shared<FitJob242>();
  { std::lock_guard<std::mutex> lk(g_fw->m); if (!g_fw->started) g_fw->start(); g_fw->q.push_back({std::move(f), h.j}); }
  g_fw->cv.notify_one();
  return h;
}
#include <atomic>
#include <mach/mach_time.h>
#include <unordered_set>
static const int g_h890 = (getenv("H890") ? atoi(getenv("H890")) : 0) | ((getenv("H917") ? atoi(getenv("H917")) : 0) & 4 ? 1 : 0);
#include "single_fast_h1000.hpp"
#include "compile_memo_i17_h156.h"
#include "nf_flat_h16_i13_h156.hpp"
#include "nf_flat_h16_h264.hpp"
#include "compile_h308c.hpp"
#include "compile_h720.hpp"
#include "final_pick_h.hpp"
#include "search_pen_h.hpp"
#include "stroke_input_f.hpp"
#include "nf_flat_h16_h890.hpp"
#include "evo_loop_i12_h744.hpp"
#include <unordered_map>
#include "sicis_host.hpp"
#include "ui_bridge_h.hpp"
#include "geo_pick_h.hpp"
#include "fit_export_h.hpp"
#include <sstream>
#include <mach/mach.h>
#include <sys/sysctl.h>
#include <memory>
#include "../s10_band/band_lib.h"
#include "memo_retain_h156.hpp"
#include <thread>
#include "scan_cpu_h290.hpp"
using cmi6::TreeEnc;

using clk = std::chrono::steady_clock;
static double ms_since(clk::time_point t) { return std::chrono::duration<double, std::milli>(clk::now() - t).count(); }
static uint64_t splitmix64(uint64_t x) { x += 0x9e3779b97f4a7c15ull; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull; return x ^ (x >> 31); }
static uint64_t hash_rec(const std::vector<uint32_t>& r) { uint64_t h = 0x243f6a8885a308d3ull ^ r.size(); for (uint32_t w : r) h = splitmix64(h ^ w); return h; }
static const float DEAD = 1.0e30f;

#define TFIELDS(X) X(uq) X(breed) X(nf_only) X(cmp_only) X(nf) X(dedup) X(comp) X(scan_enc) X(scan_wall) X(scan_gpu) X(gn_enc) X(gn_wall) X(gn_gpu) \
  X(glue) X(select) X(dump) X(check) X(nchild) X(admitted) X(rej_within) X(rej_pop) X(shortfall) X(nuniq) X(nfit0) X(nfit1) \
  X(nf_fail) X(cmp_fail) X(scan_fail) X(gn_drop) X(gn_kmiss) X(s_alive) X(g_alive) X(nf_hit) X(nf_miss) X(nf_drop) \
  X(steps_sum) X(steps_n) X(gn_cand) X(gn_floor) X(sel_skip) X(sel_skip_mig) X(surv_child) X(xfill) X(enc_tree) X(pipe_wall) X(w_s) X(w_g) X(rebuild) X(head) X(glue_e) X(gpu_union) \
  X(p_enc) X(p_grid) X(p_band) X(p_pa) X(p_pb) X(p_pa_wall) X(p_pb_wall) X(p_cpu) X(p_cand) X(p_acc0) X(p_acc1) X(p_viol) X(p_fb) \
  X(p_chk_n) X(p_chk_objmax) X(p_chk_mis) X(p_chk_mono) X(p_chk_corr) X(p_ndrop) X(p_fin_n) X(p_fin_bad) X(p_fin_max) X(p_fin_lo) X(p_fin_lobad) X(p_wall) X(p_ebmb) X(p_side_up) \
  X(h_s_alloc) X(h_s_pre) X(h_s_sort) X(h_s_bins) X(h_s_fill) X(h_s_cmd) X(h_s_pipe) X(h_g_alloc) X(h_g_pre) X(h_g_sort) X(h_g_plan) X(h_g_fill) X(h_g_cmd) X(h_g_pipe) \
  X(h131_iter) X(h131_pre) X(h131_alloc) X(h131_bind) X(h131_enc) X(h131_inst) X(h131_h4wait) X(h131_h4) \
  X(u_p1) X(u_find) X(u_pref) X(u_ctor) X(u_comp) X(u_ids) X(u_ct) X(c_runs) X(c_recs) X(c_wall) X(c_busy) X(c_cpu) X(c_lat) X(c_tail) X(c_ord)


static const int g_h257 = getenv("H257") ? atoi(getenv("H257")) : 0;
static const int g_h930_memfix = getenv("H930_MEMFIX") ? atoi(getenv("H930_MEMFIX")) : 0;
static const int g_h899 = getenv("H899") ? atoi(getenv("H899")) : 0;
static const int g_h898 = getenv("H898") ? atoi(getenv("H898")) : 0;
static double g_h898_loop_ms = 0, g_h898_fin_ms = 0, g_h898_arr_ms = 0, g_h898_bin_ms = 0, g_h898_grp_ms = 0;
static int g_h897 = getenv("H897") ? atoi(getenv("H897")) : 0;
static const int g_h894 = getenv("H894") ? atoi(getenv("H894")) : 0;

static const int g_h900 = getenv("H900") ? atoi(getenv("H900")) : 0;

#include <malloc/malloc.h>
typedef uint64_t mach_vm_address_t; typedef uint64_t mach_vm_size_t;
static inline kern_return_t mach_vm_region(...) { return KERN_FAILURE; }
#include <mach/mach.h>
extern "C" { void* objc_autoreleasePoolPush(void); void objc_autoreleasePoolPop(void*); }
static const int g_h917 = getenv("H917") ? atoi(getenv("H917")) : 0;
static const int g_h917_every = getenv("H917_EVERY") ? std::max(1, atoi(getenv("H917_EVERY"))) : 1;
static double g_h917_relief_ms = 0; static int g_h917_relief_n = 0;
static double now_ms_h917() { return (double)clock_gettime_nsec_np(CLOCK_UPTIME_RAW) / 1e6; }
static bool h917_at(int st) {
  const char* e = getenv("H917_AT"); std::string s = e ? e : "1,25,50"; s = "," + s + ",";
  return s.find("," + std::to_string(st) + ",") != std::string::npos;
}
static void h917_report(int st) {
  const double MB = 1048576.0, pg = (double)vm_page_size; task_vm_info_data_t vi; mach_msg_type_number_t cnt = TASK_VM_INFO_COUNT; task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vi, &cnt);
  malloc_statistics_t ms; malloc_zone_statistics(nullptr, &ms);
  double res[4] = {0, 0, 0, 0}, swp[4] = {0, 0, 0, 0}; mach_vm_address_t a = 0;
  for (;;) { mach_vm_size_t sz = 0; vm_region_extended_info_data_t ei; mach_msg_type_number_t c = VM_REGION_EXTENDED_INFO_COUNT; mach_port_t on = 0;
    if (mach_vm_region(mach_task_self(), &a, &sz, VM_REGION_EXTENDED_INFO, (vm_region_info_t)&ei, &c, &on) != KERN_SUCCESS) break;
    const int t = ei.user_tag & 255; const int k = t == 2 ? 0 : t == 3 ? 1 : t == 7 ? 2 : t == 11 ? 3 : -1; if (k >= 0) { res[k] += ei.pages_resident * pg; swp[k] += ei.pages_swapped_out * pg; } a += sz; }
  fprintf(stderr, "H917 st=%d footprint %.1f MB | compressed %.1f | resident %.1f | malloc blocks %u in-use %.1f MB | SMALL res %.1f swp %.1f | LARGE res %.1f swp %.1f | TINY res %.1f swp %.1f | NANO res %.1f swp %.1f | relief %d calls %.2f ms\n",
    st, vi.phys_footprint / MB, vi.compressed / MB, vi.resident_size / MB, ms.blocks_in_use, ms.size_in_use / MB, res[0] / MB, swp[0] / MB, res[1] / MB, swp[1] / MB, res[2] / MB, swp[2] / MB, res[3] / MB, swp[3] / MB, g_h917_relief_n, g_h917_relief_ms);
  g_h917_relief_n = 0; g_h917_relief_ms = 0;
}
static int g_h900_st = 0;
#include <sys/mman.h>
#include <pthread.h>
#include <execinfo.h>
#include <mach-o/dyld.h>
static double h900_now_ms() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
struct H900Reg { const char* nm; int line; int st; const void* p; size_t n; bool shared; bool main; };
static std::vector<H900Reg> g_h900_reg; static std::mutex g_h900_mu;
static const char* h900_thr() { return pthread_main_np() ? "main" : "other"; }
static void h900_note(const char* nm, int line, size_t n, bool shared, id<MTLBuffer> b) {
  if (n < ((size_t)1 << 20)) return;
  std::lock_guard<std::mutex> lk(g_h900_mu);
  g_h900_reg.push_back({nm, line, g_h900_st, shared ? b.contents : nullptr, n, shared, (bool)pthread_main_np()});
  void* fr[8]; const int nf = backtrace(fr, 8); const intptr_t slide = _dyld_get_image_vmaddr_slide(0);
  fprintf(stderr, "H900 alloc st=%d %s line=%d %.3f MB %s thr=%s @", g_h900_st, nm, line, (double)n / 1048576.0, shared ? "shared" : "private", h900_thr());
  for (int k = 1; k < nf; k++) fprintf(stderr, " 0x%lx", (unsigned long)((intptr_t)fr[k] - slide)); fprintf(stderr, "\n");
}
static void h900_resident(int st) {
  std::lock_guard<std::mutex> lk(g_h900_mu);
  for (auto& r : g_h900_reg) { if (!r.p) continue; const size_t pg = 16384; const size_t np = (r.n + pg - 1) / pg; std::vector<char> v(np);
    if (mincore((void*)r.p, r.n, v.data()) != 0) continue; size_t c = 0; for (char x : v) c += (x & 1);
    fprintf(stderr, "H900 resident end_st=%d %s line=%d allocated_st=%d alloc %.3f MB touched %.3f MB\n", st, r.nm, r.line, r.st, (double)r.n / 1048576.0, (double)c * pg / 1048576.0); }
}
static id<MTLBuffer> h900_nbl(id<MTLDevice> dev, int line, size_t n, MTLResourceOptions o) {
  id<MTLBuffer> b = [dev newBufferWithLength:n options:o]; memdiag::metal_note(b, n, (o & MTLResourceStorageModePrivate) != 0);
  if (g_h900 & 1) h900_note("nbl", line, n, !(o & MTLResourceStorageModePrivate), b);
  return b;
}
static id<MTLBuffer> h900_nbb(id<MTLDevice> dev, int line, const void* p, size_t n, MTLResourceOptions o) {
  id<MTLBuffer> b = [dev newBufferWithBytes:p length:n options:o]; memdiag::metal_note(b, n, (o & MTLResourceStorageModePrivate) != 0);
  if (g_h900 & 1) h900_note("nbb", line, n, !(o & MTLResourceStorageModePrivate), b);
  return b;
}

static const int g_h262 = getenv("H262") ? atoi(getenv("H262")) : 0;


static const int g_h263_raw = getenv("H263") ? atoi(getenv("H263")) : 0;
static const int g_h263 = g_h263_raw == 1 ? 7 : (g_h263_raw & 7);


static const int g_h264 = (getenv("H264") ? atoi(getenv("H264")) : 0) || (getenv("H308") ? atoi(getenv("H308")) != 0 : 0);
static const int g_h308 = getenv("H308") ? atoi(getenv("H308")) : 0;

static const int g_h843 = getenv("H843") ? atoi(getenv("H843")) : 0;
static inline bool compile_enc_sel(cmi6::Handle& H, const View& V, Tree& T, cmi6::TreeEnc& te, double& em, const uint64_t* hp) {
  if (g_h843 && g_h308 >= 3) return c720::compile_enc(H, V, T, te, em, hp);
  return g_h308 >= 3 ? c308::compile_enc(H, V, T, te, em, hp) : cmi6i17::compile_enc_i17(H, V, T, te, em, hp);
}
static const double g_h257f = getenv("H257F") ? atof(getenv("H257F")) : -1.0;

struct Times {
#define DECL(f) double f = 0;
  TFIELDS(DECL)
#undef DECL
  double timed() const { return breed + pipe_wall + glue + select + rebuild; }
};
static void addT(Times& a, const Times& b) {
#define ADD(f) a.f += b.f;
  TFIELDS(ADD)
#undef ADD
}
static void scaleT(Times& a, double s) {
#define SC(f) a.f *= s;
  TFIELDS(SC)
#undef SC
}


struct Stroke { std::string id, label; std::vector<double> xs, ys; };
static bool load_strokes(const char* path, std::vector<Stroke>& out) {
  std::ifstream f(path); if (!f) return false;
  int cnt = 0; f >> cnt;
  for (int k = 0; k < cnt; k++) {
    Stroke s; int n = 0; f >> s.id >> s.label >> n;
    if (!f || n < 2) return false;
    s.xs.resize(n); s.ys.resize(n);
    for (int i = 0; i < n; i++) f >> s.xs[i] >> s.ys[i];
    if (!f) return false;
    out.push_back(s);
  }
  return true;
}


static void to_npt(const std::vector<double>& X, const std::vector<double>& Y, std::vector<double>& ox, std::vector<double>& oy) {
  if ((int)X.size() == NPT) { ox = X; oy = Y; return; }
  std::vector<double> cx{X[0]}, cy{Y[0]};
  for (size_t i = 1; i < X.size(); i++) if (X[i] != cx.back() || Y[i] != cy.back()) { cx.push_back(X[i]); cy.push_back(Y[i]); }
  std::vector<double> s{0.0};
  for (size_t i = 1; i < cx.size(); i++) s.push_back(s.back() + std::hypot(cx[i] - cx[i - 1], cy[i] - cy[i - 1]));
  const double L = s.back();
  ox.assign(NPT, 0); oy.assign(NPT, 0);
  size_t j = 0;
  for (int i = 0; i < NPT; i++) {
    const double t = L * i / (NPT - 1);
    while (j + 1 < s.size() - 1 && s[j + 1] < t) j++;
    const double seg = s[j + 1] - s[j];
    double a = seg == 0 ? 0.0 : (t - s[j]) / seg;
    a = std::min(std::max(a, 0.0), 1.0);
    ox[i] = cx[j] + a * (cx[j + 1] - cx[j]);
    oy[i] = cy[j] + a * (cy[j + 1] - cy[j]);
  }
}

struct Target { double y[NPT]; float yf[NPT]; double yy; float floor_; };


struct alignas(128) Worker { Ctx* C = nullptr; nf16::Handle* H = nullptr; nf308::St* S308 = nullptr; cmi6::Handle* CH = nullptr; double nf_only = 0, nf_fail = 0, cmp_only = 0, enc_only = 0; };


static int h178_on() { static const int v = getenv("H178") ? atoi(getenv("H178")) : 1; return v; }
static int h178_pre() { static const int v = getenv("H178P") ? atoi(getenv("H178P")) : h178_on(); return v; }
static std::string h178_path() { const char* p = getenv("H178F"); return p && *p ? std::string(p) : std::string("./psos_i22.bin"); }
static id<MTLBinaryArchive> g_h178_arch = nil;


static int h876_mask() { static const int v = getenv("H876") ? atoi(getenv("H876")) : 0; return v; }
static bool g_h178_exists = false, g_h178_dirty = false;
static long g_h178_created = 0, g_h178_absent = 0, g_h178_hits = 0, g_h178_miss = 0, g_h178_addfail = 0;
static double g_h178_ms = 0;
static bool g_h178_released = false;
static id<MTLComputePipelineState> h178_make(id<MTLDevice> dev, id<MTLFunction> fn) {
  if (!g_h178_arch) { NSError* e0 = nil; id<MTLComputePipelineState> p0 = [dev newComputePipelineStateWithFunction:fn error:&e0]; if (!p0) { fprintf(stderr, "pso h876: %s\n", e0.description.UTF8String); exit(1); } return p0; }
  MTLComputePipelineDescriptor* d = [MTLComputePipelineDescriptor new];
  d.computeFunction = fn;
  if (g_h178_exists) {
    d.binaryArchives = @[ g_h178_arch ]; NSError* e = nil;
    id<MTLComputePipelineState> ps = [dev newComputePipelineStateWithDescriptor:d options:MTLPipelineOptionFailOnBinaryArchiveMiss reflection:nil error:&e];
    if (ps) { g_h178_hits++; return ps; }
    g_h178_miss++;
  }

  NSError* ea = nil; d.binaryArchives = nil;
  if ([g_h178_arch addComputePipelineFunctionsWithDescriptor:d error:&ea]) g_h178_dirty = true;
  else { g_h178_addfail++; if (g_h178_addfail <= 3) fprintf(stderr, "H178: addComputePipelineFunctions: %s\n", ea.description.UTF8String); }
  d.binaryArchives = @[ g_h178_arch ]; NSError* e2 = nil;
  id<MTLComputePipelineState> ps = [dev newComputePipelineStateWithDescriptor:d options:MTLPipelineOptionNone reflection:nil error:&e2];
  if (!ps) { fprintf(stderr, "H178: pso from archive failed (%s); plain compile\n", e2.description.UTF8String);
    NSError* e3 = nil; ps = [dev newComputePipelineStateWithFunction:fn error:&e3];
    if (!ps) { fprintf(stderr, "pso %s: %s\n", "h178", e3.description.UTF8String); exit(1); } }
  return ps; }

static void h178_flush() {
  if (!g_h178_arch || !g_h178_dirty) return;
  if ((h876_mask() & 2) && g_h178_exists) return;
  NSFileManager* fm = [NSFileManager defaultManager]; const char* hw = getenv("H178W");
  NSString* ps = [NSString stringWithUTF8String:(hw && *hw ? std::string(hw) : h178_path()).c_str()];
  NSString* tmp = [ps stringByAppendingString:@".tmp"];
  [fm removeItemAtPath:tmp error:nil];
  NSError* e = nil;
  if (![g_h178_arch serializeToURL:[NSURL fileURLWithPath:tmp] error:&e]) { fprintf(stderr, "H178: serialize failed: %s\n", e.description.UTF8String); [fm removeItemAtPath:tmp error:nil]; return; }
  [fm removeItemAtPath:ps error:nil];
  if (![fm moveItemAtPath:tmp toPath:ps error:&e]) { fprintf(stderr, "H178: archive move failed: %s\n", e.description.UTF8String); [fm removeItemAtPath:tmp error:nil]; return; }
  NSDictionary* at = [fm attributesOfItemAtPath:ps error:nil];
  printf("H178 archive written: %s | %llu bytes | pipelines created this run %ld | archive hits %ld misses %ld add failures %ld\n",
         h178_path().c_str(), (unsigned long long)[at fileSize], g_h178_created, g_h178_hits, g_h178_miss, g_h178_addfail);
  fflush(stdout); }
static void h178_open(id<MTLDevice> dev) {
  NSString* ps = [NSString stringWithUTF8String:h178_path().c_str()];
  g_h178_exists = [[NSFileManager defaultManager] fileExistsAtPath:ps];
  MTLBinaryArchiveDescriptor* ad = [MTLBinaryArchiveDescriptor new];
  if (g_h178_exists) ad.url = [NSURL fileURLWithPath:ps];
  NSError* e = nil; g_h178_arch = [dev newBinaryArchiveWithDescriptor:ad error:&e];
  if (!g_h178_arch && g_h178_exists) {
    fprintf(stderr, "H178: archive %s unreadable (%s); starting a new one\n", h178_path().c_str(), e.description.UTF8String);
    g_h178_exists = false; e = nil; g_h178_arch = [dev newBinaryArchiveWithDescriptor:[MTLBinaryArchiveDescriptor new] error:&e]; }
  if (!g_h178_arch) { fprintf(stderr, "H178: no binary archive (%s); compiling without it\n", e ? e.description.UTF8String : "?"); g_h178_exists = false; return; }
  g_h178_dirty = !g_h178_exists;
  printf("H178=1: binary archive %s (%s)\n", h178_path().c_str(), g_h178_exists ? "loaded" : "new, written at exit");
  atexit(h178_flush); }

struct Engine {
  id<MTLDevice> dev; id<MTLCommandQueue> cq; id<MTLLibrary> Ls, Lg, Lp, Lgr, Lb, Lgp, Lpp;
  std::map<std::string, id<MTLComputePipelineState>> pso;
  std::mutex pso_m;
  nf16::Shared* SH; cc14::Shared* CS; cmi6::Memo* MM; std::vector<Worker> W; evo::Pool* pool = nullptr; Target tg[2];
  void set_target(int s, const std::vector<double>& v) {
    Target& T = tg[s]; T.yy = 0;
    for (int t = 0; t < NPT; t++) { T.y[t] = v[t]; T.yf[t] = (float)v[t]; T.yy += v[t] * v[t]; }
    T.floor_ = (float)(1e-5 * T.yy);
  }
  id<MTLComputePipelineState> pipe(id<MTLLibrary> L, const std::string& name) {
    std::unique_lock<std::mutex> pso_lk(pso_m);
    auto it = pso.find(name); if (it != pso.end()) return it->second;
    if (g_h178_arch || g_h178_released) {
      id<MTLFunction> fn = [L newFunctionWithName:[NSString stringWithUTF8String:name.c_str()]];
      id<MTLComputePipelineState> ps = nil;
      if (fn) { auto tc = clk::now(); ps = h178_make(dev, fn); g_h178_created++; g_h178_ms += ms_since(tc); } else g_h178_absent++;
      pso[name] = ps; return ps; }
    pso_lk.unlock();
    id<MTLFunction> fn = [L newFunctionWithName:[NSString stringWithUTF8String:name.c_str()]];
    id<MTLComputePipelineState> ps = nil; double cms = 0;
    if (fn) { auto tc = clk::now(); NSError* e = nil; ps = [dev newComputePipelineStateWithFunction:fn error:&e];
      if (!ps) { fprintf(stderr, "pso %s: %s\n", name.c_str(), e.description.UTF8String); exit(1); } cms = ms_since(tc); }
    pso_lk.lock();
    auto it2 = pso.find(name); if (it2 != pso.end()) return it2->second;
    if (fn) { g_h178_created++; g_h178_ms += cms; } else g_h178_absent++;
    pso[name] = ps; return ps; }
  id<MTLBuffer> B(const void* p, size_t n) { if (!p || n == 0) return h900_nbl(dev, __LINE__, 16, MTLResourceStorageModeShared); return h900_nbb(dev, __LINE__, p, n, MTLResourceStorageModeShared); }
};


static thread_local evo::Pool* tl_pool = nullptr;
static thread_local std::vector<Worker>* tl_W = nullptr;
static evo::Pool* g_pool2 = nullptr;
static std::vector<Worker> g_W2;
static int g_h217 = 0;
static evo::Pool& PL(Engine& E) { return tl_pool ? *tl_pool : *E.pool; }
static std::vector<Worker>& WK(Engine& E) { return tl_W ? *tl_W : E.W; }
static const std::vector<Worker>& WK(const Engine& E) { return tl_W ? *tl_W : E.W; }
static int g_h870 = 0;
struct Unique { std::vector<Tree> U; std::vector<uint64_t> hash; std::vector<TreeEnc> enc; };
#include "h308d_take.hpp"


static const int g_h301r = getenv("H301R") ? atoi(getenv("H301R")) : 1;
static std::vector<Tree> g_h301r_t; static std::vector<TreeEnc> g_h301r_e; static std::atomic<size_t> g_h301r_n{0};


static thread_local h308d::Walk tl_h308d;

static const int g_h308x = getenv("H308X") ? atoi(getenv("H308X")) : 0;


static const int g_h313a_mask = g_h308 >= 5 ? (getenv("H313A") ? atoi(getenv("H313A")) : 7) : 0;
static const int g_h313x = (g_h308 >= 5 && getenv("H313X")) ? atoi(getenv("H313X")) : 0;
static const bool g_h313_defer = (g_h313a_mask & 5) != 0;
static const int g_h313dbg = getenv("H313DBG") ? atoi(getenv("H313DBG")) : 0;
static const bool g_h313_init = (c308::g_a2 = c720::g_a2 = (g_h313a_mask & 1) != 0, cmi6::g_a3 = (g_h313a_mask & 2) != 0, cmi6::g_a1 = (g_h313a_mask & 4) != 0, true);
static bool g_t4 = g_h308 >= 4;
static const int g_h308db = getenv("H308DB") ? atoi(getenv("H308DB")) : 1;
static thread_local unsigned tl_h308d_nb = 0;
static inline void h308d_claim(size_t nb) { h308d::begin(tl_h308d, nb, g_h308db == 100 ? (int)(tl_h308d_nb++ & 1) : g_h308db); }
static inline void h301r_take(Tree& t, TreeEnc& e) {
  if (g_t4) { h308d::take(tl_h308d, g_h301r_n, g_h301r_t, g_h301r_e, t, e); return; }
  const size_t r = g_h301r_n.fetch_add(1, std::memory_order_relaxed);
  if (r < g_h301r_t.size()) { std::swap(t, g_h301r_t[r]); std::swap(e, g_h301r_e[r]); }
  t.Q = 0; t.nth = 0; t.nS = 0; t.nslot = 0; t.prog.clear(); t.sp.clear(); t.th.clear(); t.rec.clear();
  std::memset(t.segs, 0, sizeof t.segs);
  e.s_ok = 0; e.g_ok = 0; std::memset(e.sg, 0, sizeof e.sg); std::memset(e.gg, 0, sizeof e.gg); e.ghw = 0; e.gper = 0;
  e.sp4.clear(); e.gp4.clear();
}
struct FitOut {
  std::vector<float> rfin, smin; std::vector<int> sarg; std::vector<char> is_g;
  void init(int M) { rfin.assign(M, INFINITY); smin.assign(M, INFINITY); sarg.assign(M, -1); is_g.assign(M, 0); }
};
struct StepBatch { std::vector<int> ct; Unique UQ; FitOut F[2]; std::vector<int> ids[2]; };


static const int g_h874 = getenv("H874") ? atoi(getenv("H874")) : 0;


class Key {
  std::vector<uint32_t> v_; const uint32_t* vp_ = nullptr; size_t vn_ = 0; bool view_ = false;
public:
  Key() = default;
  Key(const Key& o) { v_.assign(o.data(), o.data() + o.size()); }
  Key(Key&& o) noexcept { if (o.view_) v_.assign(o.vp_, o.vp_ + o.vn_); else v_ = std::move(o.v_); }
  Key& operator=(const Key& o) { if (this != &o) { view_ = false; v_.assign(o.data(), o.data() + o.size()); } return *this; }
  Key& operator=(Key&& o) noexcept { if (this != &o) { view_ = false; if (o.view_) v_.assign(o.vp_, o.vp_ + o.vn_); else v_ = std::move(o.v_); } return *this; }
  const uint32_t* data() const { return view_ ? vp_ : v_.data(); }
  size_t size() const { return view_ ? vn_ : v_.size(); }
  bool empty() const { return size() == 0; }
  const uint32_t* begin() const { return data(); }
  const uint32_t* end() const { return data() + size(); }
  uint32_t operator[](size_t i) const { return data()[i]; }
  void clear() { view_ = false; v_.clear(); }
  void assign(const uint32_t* a, const uint32_t* b) { view_ = false; v_.assign(a, b); }
  void set_view(const uint32_t* p, size_t n) { view_ = true; vp_ = p; vn_ = n; }
  bool is_view() const { return view_; }
  friend bool operator==(const Key& a, const Key& b) { return a.size() == b.size() && (a.size() == 0 || memcmp(a.data(), b.data(), a.size() * sizeof(uint32_t)) == 0); }
};
static uint64_t hash_rec(const Key& r) { uint64_t h = 0x243f6a8885a308d3ull ^ r.size(); for (size_t i = 0, n = r.size(); i < n; i++) h = splitmix64(h ^ r.data()[i]); return h; }
struct KeyHash { size_t operator()(const Key& k) const noexcept { return (size_t)hash_rec(k); } };


struct FlatKeySet {
  std::vector<FlatKeySet> sh;
  static unsigned shard_of(uint64_t h) { return (unsigned)(h >> 60); }
  std::vector<uint32_t> slot;
  std::vector<uint64_t> eh, eoff;
  std::vector<uint32_t> elen;
  std::vector<uint32_t> arena;
  bool inc = false; std::vector<uint8_t> dead; size_t ndead = 0;
  bool byref = false; std::vector<const uint32_t*> eptr;
  size_t size() const { if (!sh.empty()) { size_t s = 0; for (const auto& x : sh) s += x.size(); return s; } return eh.size() - ndead; }
  void clear() { for (auto& x : sh) x.clear(); inc = false; dead.clear(); ndead = 0; slot.clear(); eh.clear(); eoff.clear(); elen.clear(); arena.clear(); eptr.clear(); }
  void rehash(size_t cap) {
    slot.assign(cap, 0);
    const uint64_t m = cap - 1;
    for (size_t e = 0; e < eh.size(); e++) { uint64_t p = eh[e] & m; while (slot[p]) p = (p + 1) & m; slot[p] = (uint32_t)e + 1; }
  }
  void reserve(size_t n) { if (g_h874 & 1) { eh.reserve(n); eoff.reserve(n); elen.reserve(n); }
    size_t cap = 1024; while (cap < 2 * n + 2) cap <<= 1; if (cap > slot.size()) rehash(cap); }
  bool same(size_t e, const uint32_t* k, uint32_t n) const { return elen[e] == n && (n == 0 || memcmp(byref ? eptr[e] : arena.data() + eoff[e], k, sizeof(uint32_t) * n) == 0); }
  uint32_t find_or_add(uint64_t h, const uint32_t* k, uint32_t n, bool& fresh) {
    if ((eh.size() + 1) * 2 > slot.size()) rehash(slot.empty() ? 1024 : slot.size() * 2);
    const uint64_t m = slot.size() - 1;
    uint64_t p = h & m;
    while (slot[p]) {
      const size_t e = slot[p] - 1;
      if (eh[e] == h && same(e, k, n)) { fresh = false; return (uint32_t)e; }
      p = (p + 1) & m;
    }
    const uint32_t e = (uint32_t)eh.size();
    eh.push_back(h); elen.push_back(n);
    if (byref) eptr.push_back(k); else { eoff.push_back(arena.size()); arena.insert(arena.end(), k, k + n); }
    slot[p] = e + 1; fresh = true;
    return e;
  }
  bool contains(uint64_t h, const uint32_t* k, uint32_t n) const {
    if (!sh.empty()) return sh[shard_of(h)].contains(h, k, n);
    if (slot.empty()) return false;
    const uint64_t m = slot.size() - 1;
    uint64_t p = h & m;
    while (slot[p]) { const size_t e = slot[p] - 1; if (eh[e] == h && same(e, k, n)) return !(inc && dead[e]); p = (p + 1) & m; }
    return false;
  }
  bool insert(uint64_t h, const uint32_t* k, uint32_t n) { if (!sh.empty()) return sh[shard_of(h)].insert(h, k, n); bool f = false; find_or_add(h, k, n, f); return f; }

  inline void pf_slot(uint64_t h) const { if (!sh.empty()) { sh[shard_of(h)].pf_slot(h); return; } if (!slot.empty()) __builtin_prefetch(&slot[h & (slot.size() - 1)], 0, 3); }
  inline void pf_entry(uint64_t h) const {
    if (!sh.empty()) { sh[shard_of(h)].pf_entry(h); return; }
    if (slot.empty()) return;
    const uint32_t s = slot[h & (slot.size() - 1)]; if (s) __builtin_prefetch(&eh[s - 1], 0, 3);
  }
  inline void pf_arena(uint64_t h) const {
    if (!sh.empty()) { sh[shard_of(h)].pf_arena(h); return; }
    if (slot.empty()) return;
    const uint32_t s = slot[h & (slot.size() - 1)];
    if (s && eh[s - 1] == h) { __builtin_prefetch(&elen[s - 1], 0, 3); if (byref) __builtin_prefetch(&eptr[s - 1], 0, 3); else { __builtin_prefetch(&eoff[s - 1], 0, 3); __builtin_prefetch(arena.data() + eoff[s - 1], 0, 3); } }
  }
};

struct Meta { Key key; uint64_t r0 = 0, r1 = 0, kh = 0; uint8_t kf = 0, nx = 0, ny = 0; };


static const int g_sp_env = getenv("SEARCH_PEN") ? atoi(getenv("SEARCH_PEN")) : 0;
static const int g_sp_diag = getenv("SEARCH_PEN_DIAG") ? atoi(getenv("SEARCH_PEN_DIAG")) : 0;
static bool g_sp_on = false;
static std::atomic<long> g_sp_kunk{0};
static inline int sp_show_cx(const evo::Ind& I) { return g_sp_on ? (int)evo::cplx(I.g) : (int)I.cx; }
static inline uint16_t sp_cx(const evo::Cand& g, const Meta& m) { return g_sp_on ? sp::level(m.nx, m.ny, m.kf) : (uint16_t)evo::cplx(g); }

static std::atomic<uint64_t> g_h890_key_ns{0}; static std::atomic<uint64_t> g_h890_key_bytes{0};
struct H890Tm { bool on; clk::time_point t; uint64_t nb = 0; H890Tm() : on(g_h890 != 0) { if (on) t = clk::now(); }
  ~H890Tm() { if (on) { g_h890_key_ns.fetch_add((uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(clk::now() - t).count(), std::memory_order_relaxed); g_h890_key_bytes.fetch_add(nb, std::memory_order_relaxed); } } };
static inline void h890_setkey(Meta& m, const nf16::Block& B, const nf16::Child& ch, H890Tm* tm = nullptr) {
  if (g_h890 & 1) m.key.set_view(B.w.data() + ch.koff, ch.klen); else m.key.assign(B.w.data() + ch.koff, B.w.data() + ch.koff + ch.klen);
  if (tm) tm->nb += (uint64_t)ch.klen * 4;
}


struct NFData { std::vector<nf16::Child> ch; std::vector<nf16::Block> blk; std::vector<skn::Cand> sc; };
static NFData NFD;


static inline size_t al5g(size_t x, size_t a) { return (x + a - 1) / a * a; }


static const int g_mem_grow = getenv("MEM_GROW") ? atoi(getenv("MEM_GROW")) : 1;
static inline size_t mem_grow_cap(size_t bytes, size_t cap) { if (g_mem_grow && bytes >= ((size_t)1 << 20)) return al5g(bytes + bytes / 16, 16384); return std::max(bytes, cap + cap / 2); }
struct GBuf {
  id<MTLBuffer> b = nil; size_t cap = 0;
  void ensure(id<MTLDevice> dev, size_t bytes, MTLResourceOptions opt) {
    bytes = std::max<size_t>(bytes, 256);
    if (!b || cap < bytes) { const size_t nc = mem_grow_cap(bytes, cap); b = h900_nbl(dev, __LINE__, nc, opt); cap = nc; }
  }
};
static const MTLResourceOptions MSH = MTLResourceStorageModeShared, MPV = MTLResourceStorageModePrivate;
static void h900_ens(GBuf& g, const char* nm, id<MTLDevice> dev, size_t bytes, MTLResourceOptions opt) {
  if (!(g_h900 & 1)) { g.ensure(dev, bytes, opt); return; }
  const size_t oc = g.cap; g.ensure(dev, bytes, opt);
  if (g.cap != oc) h900_note(nm, 0, g.cap, !(opt & MTLResourceStorageModePrivate), g.b);
}
static const size_t BIN_ALIGN = 64;
static size_t al5(size_t x, size_t a) { return (x + a - 1) / a * a; }

static void sort_by_key5(const std::vector<uint32_t>& seq, const std::vector<uint32_t>& key, size_t DN,
                         std::vector<uint32_t>& start, std::vector<uint32_t>& fill, std::vector<uint32_t>& out) {
  start.assign(DN + 1, 0);
  for (uint32_t x : seq) start[key[x] + 1]++;
  for (size_t k = 0; k < DN; k++) start[k + 1] += start[k];
  fill.assign(start.begin(), start.begin() + (std::ptrdiff_t)DN);
  out.resize(seq.size());
  for (uint32_t x : seq) out[fill[key[x]]++] = x;
}


struct Tick { clk::time_point t = clk::now(); double lap() { const auto n = clk::now(); const double d = std::chrono::duration<double, std::milli>(n - t).count(); t = n; return d; } };
static int g_h916 = 0;
static int g_h215 = 0;
static bool g_h93 = true, g_h93_print = true;
static thread_local std::vector<uint32_t> h93_hist;

static std::mutex h93_pso_m;

static id<MTLComputePipelineState> h93_scan_ps(Engine& E, int s2, int Q) {
  std::lock_guard<std::mutex> lk(h93_pso_m);
  static id<MTLComputePipelineState> cache[2][128]; static bool have[2][128];
  if (Q < 0 || Q >= 128) return E.pipe(E.Ls, s2 ? "rcscan_s2_q" + std::to_string(Q) : "rcscan_q" + std::to_string(Q));
  if (!have[s2][Q]) { cache[s2][Q] = E.pipe(E.Ls, s2 ? "rcscan_s2_q" + std::to_string(Q) : "rcscan_q" + std::to_string(Q)); have[s2][Q] = true; }
  return cache[s2][Q];
}
static id<MTLComputePipelineState> h93_gn_ps(Engine& E, int Q, int Tk) {
  std::lock_guard<std::mutex> lk(h93_pso_m);
  static id<MTLComputePipelineState> cache[128][25]; static bool have[128][25];
  if (Q < 0 || Q >= 128 || Tk < 0 || Tk >= 25) return E.pipe(E.Lg, "gnw_q" + std::to_string(Q) + "_t" + std::to_string(Tk));
  if (!have[Q][Tk]) { cache[Q][Tk] = E.pipe(E.Lg, "gnw_q" + std::to_string(Q) + "_t" + std::to_string(Tk)); have[Q][Tk] = true; }
  return cache[Q][Tk];
}


static void sort_by_key5_p(Engine& E, const std::vector<uint32_t>& seq, const std::vector<uint32_t>& key, size_t DN,
                           std::vector<uint32_t>& start, std::vector<uint32_t>& out) {
  const size_t n = seq.size();
  out.resize(n);
  const int NB = (int)std::min<size_t>((size_t)(4 * WK(E).size()), n / 1024);
  if (NB < 2) {
    start.assign(DN + 1, 0);
    for (uint32_t x : seq) start[key[x] + 1]++;
    for (size_t k = 0; k < DN; k++) start[k + 1] += start[k];
    std::vector<uint32_t> fill(start.begin(), start.begin() + (std::ptrdiff_t)DN);
    for (uint32_t x : seq) out[fill[key[x]]++] = x;
    return;
  }
  std::vector<uint32_t>& H = h93_hist; H.assign((size_t)NB * DN, 0);
  PL(E).run(NB, [&](int b) {
    const size_t a = n * b / NB, z = n * (b + 1) / NB; uint32_t* h = H.data() + (size_t)b * DN;
    for (size_t i = a; i < z; i++) h[key[seq[i]]]++;
  });
  start.assign(DN + 1, 0);
  PL(E).run(NB, [&](int kb) {
    const size_t k0 = DN * kb / NB, k1 = DN * (kb + 1) / NB;
    for (size_t k = k0; k < k1; k++) { uint32_t t = 0; for (int b = 0; b < NB; b++) t += H[(size_t)b * DN + k]; start[k + 1] = t; }
  });
  for (size_t k = 0; k < DN; k++) start[k + 1] += start[k];
  PL(E).run(NB, [&](int kb) {
    const size_t k0 = DN * kb / NB, k1 = DN * (kb + 1) / NB;
    for (size_t k = k0; k < k1; k++) { uint32_t run = start[k];
      for (int b = 0; b < NB; b++) { uint32_t& h = H[(size_t)b * DN + k]; const uint32_t t = h; h = run; run += t; } }
  });
  PL(E).run(NB, [&](int b) {
    const size_t a = n * b / NB, z = n * (b + 1) / NB; uint32_t* h = H.data() + (size_t)b * DN;
    for (size_t i = a; i < z; i++) { const uint32_t x = seq[i]; out[h[key[x]]++] = x; }
  });
}


static int sort_count_h263(Engine& E, const std::vector<uint32_t>& seq, const std::vector<uint32_t>& key, size_t DN, std::vector<uint32_t>& start) {
  const size_t n = seq.size();
  const int NB = std::max(1, (int)std::min<size_t>((size_t)(4 * E.W.size()), n / 1024));
  std::vector<uint32_t>& H = h93_hist; H.assign((size_t)NB * DN, 0);
  PL(E).run(NB, [&](int b) {
    const size_t a = n * b / NB, z = n * (b + 1) / NB; uint32_t* h = H.data() + (size_t)b * DN;
    for (size_t i = a; i < z; i++) h[key[seq[i]]]++;
  });
  start.assign(DN + 1, 0);
  PL(E).run(NB, [&](int kb) {
    const size_t k0 = DN * kb / NB, k1 = DN * (kb + 1) / NB;
    for (size_t k = k0; k < k1; k++) { uint32_t t = 0; for (int b = 0; b < NB; b++) t += H[(size_t)b * DN + k]; start[k + 1] = t; }
  });
  for (size_t k = 0; k < DN; k++) start[k + 1] += start[k];
  PL(E).run(NB, [&](int kb) {
    const size_t k0 = DN * kb / NB, k1 = DN * (kb + 1) / NB;
    for (size_t k = k0; k < k1; k++) { uint32_t run = start[k];
      for (int b = 0; b < NB; b++) { uint32_t& h = H[(size_t)b * DN + k]; const uint32_t t = h; h = run; run += t; } }
  });
  return NB;
}
static void sort_scatter_h263(Engine& E, const std::vector<uint32_t>& seq, const std::vector<uint32_t>& key, int NB, size_t DN,
                              const std::vector<int64_t>& shift, uint32_t* dst) {
  const size_t n = seq.size();
  std::vector<uint32_t>& H = h93_hist;
  PL(E).run(NB, [&](int b) {
    const size_t a = n * b / NB, z = n * (b + 1) / NB; uint32_t* h = H.data() + (size_t)b * DN;
    for (size_t i = a; i < z; i++) { const uint32_t x = seq[i]; const uint32_t k = key[x]; dst[(size_t)((int64_t)(h[k]++) + shift[k])] = x; }
  });
}

struct H263ScanScr { std::vector<uint32_t> bti, bsi, bhi, bns, oti, osi, ohi, osq; std::vector<int> bq; std::vector<double> bfl; std::vector<int64_t> shv; };
struct H263GnScr { std::vector<uint32_t> bgi, bgs, bns, ogt, ogs, osq; std::vector<int> bqg, bcg, blm; std::vector<double> bal, bkm;
                   std::vector<size_t> ne, perb, perb2; };
static H263ScanScr h263_keep_scan;
static H263GnScr h263_keep_gn;


struct H263Bind {
  bool on = false; const void* ps = nullptr; const void* lb[16] = {}; size_t lo[16] = {}; bool hv[16] = {};
  void buf(id<MTLComputeCommandEncoder> en, id<MTLBuffer> b, size_t off, int k) {
    if (on && hv[k] && lb[k] == (__bridge const void*)b && lo[k] == off) return;
    [en setBuffer:b offset:off atIndex:k]; lb[k] = (__bridge const void*)b; lo[k] = off; hv[k] = true;
  }
  void pipe(id<MTLComputeCommandEncoder> en, id<MTLComputePipelineState> p) {
    if (on && ps == (__bridge const void*)p) return;
    [en setComputePipelineState:p]; ps = (__bridge const void*)p;
  }
  void forget(int k) { hv[k] = false; }
  void forget_pipe() { ps = nullptr; }
};

struct SBin5 { int Q, ns, nc; bool s2; size_t off, cnt, sst; };


static int h150p_on() { static const int v = getenv("H150P") ? atoi(getenv("H150P")) : 1; return v; }
static int h150s_on() { static const int v = getenv("H150S") ? atoi(getenv("H150S")) : 0; return v; }
static int h150c_on() { static const int v = getenv("H150C") ? atoi(getenv("H150C")) : 0; return v; }
struct H150Acc { long gnDisp = 0, gnLaunch = 0, gnOcc = 0, polDisp = 0, polLaunch = 0, polOcc = 0, polBins = 0, polChanged = 0,
                      scDisp = 0, scTG = 0, scThr = 0, scSimd = 0; };
static H150Acc g_h150;
static void h150_step_reset() { g_h150 = H150Acc(); }
static size_t h150_count_tl(const std::vector<uint32_t>& tl) { size_t c = 0; for (size_t k = 1; k < tl.size(); k++) if (tl[k] != 0xFFFFFFFFu) c++; return c; }

static long h150_sim_pol(size_t pers, int t) {
  const size_t g = ((size_t)t * pers * 4 + 15) & ~(size_t)15;
  return std::min<long>(96, (long)(65536 / std::max<size_t>(g, 1)) * (long)t);
}
static int h150_pick_pol(size_t pers, size_t cnt, int cap, int t0) {
  const int tmax = std::max(1, std::min<int>(cap, (int)std::min<size_t>(cnt, 1024)));
  int best = t0; long bs = -1; int bd = 1 << 30;
  for (int t = 1; t <= tmax; t++) {
    if ((size_t)t * pers * 4 > 32768) break;
    const long s = h150_sim_pol(pers, t); const int d = std::abs(t - t0);
    if (s > bs || (s == bs && (d < bd || (d == bd && t > best)))) { best = t; bs = s; bd = d; }
  }
  return best;
}
static void h150_step_print(int st) {
  if (!h150c_on() && !h150s_on()) return;
  const H150Acc& a = g_h150;
  const double gE = a.gnLaunch ? 100.0 * (double)(a.gnLaunch - a.gnOcc) / (double)a.gnLaunch : 0.0;
  const double pE = a.polLaunch ? 100.0 * (double)(a.polLaunch - a.polOcc) / (double)a.polLaunch : 0.0;
  printf("H150 step %2d | P %d S %d C %d | GN disp %ld launched SIMD %ld occupied %ld empty %.1f%% | polish disp %ld launched SIMD %ld occupied %ld empty %.1f%% | polish bins %ld tpg changed %ld | scan disp %ld TGs %ld launched threads %ld SIMDs %ld occupied n/a\n",
         st, h150p_on(), h150s_on(), h150c_on(), a.gnDisp, a.gnLaunch, a.gnOcc, gE, a.polDisp, a.polLaunch, a.polOcc, pE,
         a.polBins, a.polChanged, a.scDisp, a.scTG, a.scThr, a.scSimd);
}

struct GBin5 { int side, Q, nth, tpg; size_t cnt, nt, ntg, tgb, per, per32 = 0, sst, tl_off, eb_off, eb_bytes; uint32_t wix; id<MTLComputePipelineState> ps; };


static int h182_on(void) { static int v = -1; if (v < 0) { const char* e = getenv("H182"); v = (e && e[0] == '0') ? 0 : 1; } return v; }


static const int g_h301 = getenv("H301") ? atoi(getenv("H301")) : 0;
static std::mutex g_h301_m;
static std::vector<std::pair<double, double>> g_h301_cb;
static std::vector<std::pair<const char*, double>> g_h301_mk;
static std::thread::id g_h301_main;
static double h301_now() { static const mach_timebase_info_data_t tb = [] { mach_timebase_info_data_t t; mach_timebase_info(&t); return t; }();
  return (double)mach_absolute_time() * (double)tb.numer / (double)tb.denom * 1e-9; }

static const int g_h401 = getenv("H401") ? atoi(getenv("H401")) : 0;
static std::atomic<int> g_h401_on{0};
static std::mutex g_h401_m;
struct H401Rec { int kind, tid; const char* lab; double a, b; };
static std::vector<H401Rec> g_h401_r;
static std::atomic<int> g_h401_tn{0};
static int h401_tid() { static thread_local int t = -1; if (t < 0) t = g_h401_tn.fetch_add(1); return t; }
static std::atomic<const char*> g_h401_cur{"(none)"};

static const int g_h892 = getenv("H892") ? atoi(getenv("H892")) : 0;
static const size_t g_h892_cap = (size_t)1 << 20;
static double g_h892_t0 = 0; static H401Rec* g_h892_ring = nullptr; static std::atomic<size_t> g_h892_n{0};
static void h401_push(int kind, const char* lab, double a, double b) { const int t = h401_tid();
  if (g_h892) { const size_t i = g_h892_n.fetch_add(1, std::memory_order_relaxed); if (i < g_h892_cap) g_h892_ring[i] = {kind, t, lab, a, b}; return; }
  std::lock_guard<std::mutex> lk(g_h401_m); g_h401_r.push_back({kind, t, lab, a, b}); }
struct H401S { const char* l; double a = 0; bool on; explicit H401S(const char* l_) : l(l_), on(g_h401_on.load(std::memory_order_relaxed) != 0) { if (on) a = h301_now(); }
  ~H401S() { if (on) h401_push(2, l, a, h301_now()); } };
#include <sys/resource.h>

static const int g_h657 = getenv("H657") ? atoi(getenv("H657")) : 0;
static struct rusage g_h657_ru0; static bool g_h657_init = false;
static void h657_step(int st) {
  if (!g_h657) return;
  struct rusage r; getrusage(RUSAGE_SELF, &r);
  if (g_h657_init) {
    auto ms = [](const timeval& a, const timeval& b) { return (a.tv_sec - b.tv_sec) * 1e3 + (a.tv_usec - b.tv_usec) * 1e-3; };
    printf("H657 step %2d: minflt %ld majflt %ld nvcsw %ld nivcsw %ld | cpu user %.1f sys %.1f ms\n", st, r.ru_minflt - g_h657_ru0.ru_minflt,
           r.ru_majflt - g_h657_ru0.ru_majflt, r.ru_nvcsw - g_h657_ru0.ru_nvcsw, r.ru_nivcsw - g_h657_ru0.ru_nivcsw,
           ms(r.ru_utime, g_h657_ru0.ru_utime), ms(r.ru_stime, g_h657_ru0.ru_stime));
  }
  g_h657_ru0 = r; g_h657_init = true;
}

static const int g_h651 = getenv("H651") ? atoi(getenv("H651")) : 0;
static const int g_h652 = getenv("H652") ? atoi(getenv("H652")) : 0;

static const char* g_h650d = getenv("H650D");
static const double g_h650ms = getenv("H650MS") ? atof(getenv("H650MS")) : 10.0;
static const double g_h650ns = getenv("H650NS") ? atof(getenv("H650NS")) : 1000.0;
static bool h650_is(const char* s) { return g_h650d && !strcmp(g_h650d, s); }
static void h650_spin_ns(double ns) {
  if (ns <= 0) return;
  const auto t0 = std::chrono::steady_clock::now();
  while (std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() < ns) { }
}
static std::atomic<uint64_t> g_h650_items{0};

static const int g_h891 = getenv("H891") ? atoi(getenv("H891")) : 0;
static const double g_h891_cpu_us = getenv("H891_CPU_US") ? atof(getenv("H891_CPU_US")) : 0;
static const int g_h891_site = getenv("H891_CPU_SITE") ? atoi(getenv("H891_CPU_SITE")) : 1;
static const double g_h891_gpu_us = getenv("H891_GPU_US") ? atof(getenv("H891_GPU_US")) : 0;
static const double g_h891_head_us = getenv("H891_HEAD_US") ? atof(getenv("H891_HEAD_US")) : 0;
static const double g_h891_tail_us = getenv("H891_TAIL_US") ? atof(getenv("H891_TAIL_US")) : 0;
static std::atomic<uint64_t> g_h891_cpu_ns{0}, g_h891_cpu_n{0};
static double g_h891_head_ms = 0, g_h891_tail_ms = 0;
static double g_h891_gpu_n = 0, g_h891_gpu_req_ms = 0, g_h891_gpu_cal_ms = 0;
static double h891_spin_ns_actual(double ns) {
  const auto t0 = std::chrono::steady_clock::now(); double el = 0;
  do { el = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count(); } while (el < ns);
  return el;
}
static inline void h891_cpu_batch(int site) {
  if (g_h891 && g_h891_cpu_us > 0 && (g_h891_site & site)) {
    const double el = h891_spin_ns_actual(g_h891_cpu_us * 1e3);
    g_h891_cpu_ns.fetch_add((uint64_t)el, std::memory_order_relaxed); g_h891_cpu_n.fetch_add(1, std::memory_order_relaxed); }
}
static inline double h891_main_spin_ms(double us) { return us > 0 ? h891_spin_ns_actual(us * 1e3) * 1e-6 : 0.0; }
static void h401_fitjob(const std::function<void()>& f) { H401S h401s_("fit job thread"); f(); }
static void h401_sub(const char* lab) { if (std::this_thread::get_id() != g_h301_main || !g_h401_on.load(std::memory_order_relaxed)) return; const double t = h301_now(); h401_push(5, lab, t, t); }
static bool h401_pool_on() { return g_h401_on.load(std::memory_order_relaxed) != 0; }
static void h401_pool_rec(double a, double b) { h401_push(1, g_h401_cur.load(std::memory_order_relaxed), a, b); }
static void h401_pool_wait(double a, double b) { h401_push(4, g_h401_cur.load(std::memory_order_relaxed), a, b); }
static void h401_dump(double t0) {
  const char* fn = getenv("H401_OUT") ? getenv("H401_OUT") : "h401_trace.csv";
  FILE* f = fopen(fn, "w"); if (!f) return;
  std::lock_guard<std::mutex> lk(g_h401_m);
  fprintf(f, "kind,tid,label,a_ms,b_ms\n");
  for (const H401Rec& r : g_h401_r) fprintf(f, "%d,%d,\"%s\",%.4f,%.4f\n", r.kind, r.tid, r.lab, (r.a - t0) * 1e3, (r.b - t0) * 1e3);
  fclose(f); g_h401_r.clear(); printf("H401: trace written to %s\n", fn);
}
static void h892_dump(double t0) {
  const char* fn = getenv("H892_OUT") ? getenv("H892_OUT") : "h892_trace.csv";
  FILE* f = fopen(fn, "w"); if (!f) return;
  const size_t n = std::min(g_h892_n.load(), g_h892_cap);
  fprintf(f, "kind,tid,label,a_ms,b_ms\n");
  for (size_t i = 0; i < n; i++) { const H401Rec& r = g_h892_ring[i]; fprintf(f, "%d,%d,\"%s\",%.4f,%.4f\n", r.kind, r.tid, r.lab, (r.a - t0) * 1e3, (r.b - t0) * 1e3); }
  fclose(f); printf("H892: %zu records%s written to %s\n", n, g_h892_n.load() > g_h892_cap ? " (RING OVERFLOW)" : "", fn);
}
static void h301_mark(const char* lab) {
  if (std::this_thread::get_id() == g_h301_main) { g_h401_cur.store(lab, std::memory_order_relaxed); if (g_h401_on.load(std::memory_order_relaxed)) { const double t = h301_now(); h401_push(0, lab, t, t); } }
  if (!g_h301 || std::this_thread::get_id() != g_h301_main) return;
  std::lock_guard<std::mutex> lk(g_h301_m); g_h301_mk.push_back({lab, h301_now()}); }


static const int g_h301d = getenv("H301D") ? atoi(getenv("H301D")) : 0;


static const int g_h301t = getenv("H301T") ? atoi(getenv("H301T")) : 0;

static const int g_h402 = getenv("H402") ? atoi(getenv("H402")) : 0;
static const double g_h402r = getenv("H402R") ? atof(getenv("H402R")) : 0.0;
static const int g_h402v = getenv("H402V") ? atoi(getenv("H402V")) : 0;
static double g_h402_rho = 1.4;


static const int g_h302 = getenv("H302") ? atoi(getenv("H302")) : 0;
static const int g_h301mb = getenv("H301M") ? atoi(getenv("H301M")) : 0;


static std::function<void()> g_h301_after_front;
static void h301_track(id<MTLCommandBuffer> cb) {
  if (!g_h301) return;
  [cb addCompletedHandler:^(id<MTLCommandBuffer> c) { std::lock_guard<std::mutex> lk(g_h301_m); g_h301_cb.push_back({c.GPUStartTime, c.GPUEndTime}); }];
}

static std::map<std::string, double> g_h301_all;
static double h301_window(double t0, double t1, std::map<std::string, double>& acc, double& idle) {
  std::vector<std::pair<double, double>> iv; std::vector<std::pair<const char*, double>> mk;
  { std::lock_guard<std::mutex> lk(g_h301_m);
    for (auto& q : g_h301_cb) if (q.second > t0 && q.first < t1) iv.push_back({std::max(q.first, t0), std::min(q.second, t1)});
    for (auto& m : g_h301_mk) if (m.second >= t0 && m.second < t1) mk.push_back(m);
    std::vector<std::pair<double, double>> kc; for (auto& q : g_h301_cb) if (q.second >= t1) kc.push_back(q); g_h301_cb.swap(kc);
    std::vector<std::pair<const char*, double>> km; for (auto& m : g_h301_mk) if (m.second >= t1) km.push_back(m); g_h301_mk.swap(km); }
  for (size_t i = 0; i < mk.size(); i++) g_h301_all[mk[i].first] += ((i + 1 < mk.size() ? mk[i + 1].second : t1) - mk[i].second) * 1e3;
  std::sort(iv.begin(), iv.end());
  std::vector<std::pair<double, double>> gaps; double busy = 0, cur = t0;
  for (auto& q : iv) { if (q.first > cur) gaps.push_back({cur, q.first}); if (q.second > cur) { busy += q.second - std::max(q.first, cur); cur = q.second; } }
  if (cur < t1) gaps.push_back({cur, t1});
  for (auto& g : gaps) {
    idle += (g.second - g.first) * 1e3;
    for (size_t i = 0; i < mk.size(); i++) {
      const double a = mk[i].second, b = i + 1 < mk.size() ? mk[i + 1].second : t1;
      const double o = std::min(b, g.second) - std::max(a, g.first);
      if (o > 0) acc[mk[i].first] += o * 1e3;
    }
    if (!mk.empty() && g.first < mk[0].second) acc["(before first mark)"] += (std::min(g.second, mk[0].second) - g.first) * 1e3;
  }
  return busy * 1e3;
}

static bool tl_on() { static const bool on = getenv("TL_LABELS") && getenv("TL_LABELS")[0] == '1'; return on; }
static id<MTLComputeCommandEncoder> tl_enc(id<MTLCommandBuffer> cb, MTLDispatchType t = MTLDispatchTypeSerial,
                                           const char* f = __builtin_FUNCTION(), int l = __builtin_LINE()) {
  id<MTLComputeCommandEncoder> en = [cb computeCommandEncoderWithDispatchType:t];
  if (tl_on()) en.label = [NSString stringWithFormat:@"%s@%d", f, l];
  return en;
}
static id<MTLCommandBuffer> h182_cb(id<MTLCommandQueue> q, const char* tlf = __builtin_FUNCTION(), int tll = __builtin_LINE()) {
  id<MTLCommandBuffer> cb;
  if (!h182_on()) cb = [q commandBuffer];
  else { MTLCommandBufferDescriptor* d = [[MTLCommandBufferDescriptor alloc] init];
    d.errorOptions = MTLCommandBufferErrorOptionEncoderExecutionStatus;
    cb = [q commandBufferWithDescriptor:d]; }
  if (tl_on()) cb.label = [NSString stringWithFormat:@"%s@%d", tlf, tll];
  h301_track(cb);
  if (g_h401_on.load(std::memory_order_relaxed)) { const char* lab = g_h401_cur.load(std::memory_order_relaxed); const int tid = h401_tid();
    [cb addCompletedHandler:^(id<MTLCommandBuffer> c) { if (g_h892) { const size_t i = g_h892_n.fetch_add(1, std::memory_order_relaxed); if (i < g_h892_cap) g_h892_ring[i] = {3, tid, lab, c.GPUStartTime, c.GPUEndTime}; } else { std::lock_guard<std::mutex> lk(g_h401_m); g_h401_r.push_back({3, tid, lab, c.GPUStartTime, c.GPUEndTime}); } }]; }
  return cb;
}
static NSString* h182_gname(const GBin5& g) {
  return [NSString stringWithFormat:@"GN side=%d Q=%d tpg=%d ntg=%lu tgb=%lu cnt=%lu", g.side, g.Q, g.tpg,
          (unsigned long)g.ntg, (unsigned long)g.tgb, (unsigned long)g.cnt];
}
static const char* h182_status(MTLCommandBufferStatus s) {
  switch (s) {
    case MTLCommandBufferStatusNotEnqueued: return "NotEnqueued"; case MTLCommandBufferStatusEnqueued: return "Enqueued";
    case MTLCommandBufferStatusCommitted: return "Committed"; case MTLCommandBufferStatusScheduled: return "Scheduled";
    case MTLCommandBufferStatusCompleted: return "Completed"; case MTLCommandBufferStatusError: return "Error";
    default: return "?";
  }
}
static const char* h182_errcode(NSInteger c) {
  switch (c) {
    case 0: return "None"; case 1: return "Internal"; case 2: return "Timeout (Caused GPU Timeout Error)"; case 3: return "PageFault";
    case 4: return "Blacklisted"; case 7: return "NotPermitted"; case 8: return "OutOfMemory"; case 9: return "InvalidResource";
    case 10: return "Memoryless"; case 11: return "DeviceRemoved"; case 12: return "StackOverflow"; default: return "unknown";
  }
}
static const char* h182_estate(NSInteger s) {
  switch (s) {
    case 0: return "Unknown"; case 1: return "Completed"; case 2: return "Affected"; case 3: return "Pending"; case 4: return "Faulted";
    default: return "?";
  }
}
#include "gpu_fail_h.hpp"
static void h182_dump(id<MTLCommandBuffer> cb, const char* site) {
  fprintf(stderr, "H182 dump site=%s status=%ld(%s) label=%s GPUStart=%.6f GPUEnd=%.6f kernelStart=%.6f kernelEnd=%.6f H182=%d\n",
          site, (long)cb.status, h182_status(cb.status), cb.label ? cb.label.UTF8String : "(nil)",
          cb.GPUStartTime, cb.GPUEndTime, cb.kernelStartTime, cb.kernelEndTime, h182_on());
  NSError* e = cb.error;
  if (!e) { fprintf(stderr, "H182 error: nil\n"); return; }
  fprintf(stderr, "H182 error domain=%s code=%ld(%s)\n", e.domain.UTF8String, (long)e.code, h182_errcode(e.code));
  fprintf(stderr, "H182 error localizedDescription=%s\n", e.localizedDescription.UTF8String);
  fprintf(stderr, "H182 error userInfo=%s\n", e.userInfo.description.UTF8String);
  NSArray* infos = e.userInfo[MTLCommandBufferEncoderInfoErrorKey];
  if (!infos) { fprintf(stderr, "H182 encoder infos: none (H182 unset, or no encoder info reported)\n"); return; }
  for (id<MTLCommandBufferEncoderInfo> info in infos) {
    fprintf(stderr, "H182 encoder label=%s errorState=%ld(%s) signposts=%s\n", info.label ? info.label.UTF8String : "(nil)",
            (long)info.errorState, h182_estate(info.errorState), [info.debugSignposts componentsJoinedByString:@" | "].UTF8String);
  }
}
struct Side5 {
  std::vector<uint32_t> cpuJ; FitHandle cpuT; double cpu_ms = 0;
  GBuf tok, hdr, sp, th, seg, bl, out;
  GBuf gtok, ghdr, gsp, gseg, gth, grs, gst, gai, gbl;
  id<MTLBuffer> yb = nil;
  std::vector<char> listed, sok, gok;
  std::vector<int> ls;
  std::vector<uint32_t> ninst, toff, spoff, thoff, gninst, gtoff, gspoff, gper;
  std::vector<uint32_t> seq, seq2, dkey, keyl, start, fill, sorted;
  std::vector<SBin5> sbins;
  uint32_t ti = 0, si = 0, hi = 0, gti = 0, gsi = 0;
  size_t bl_words = 0, gbl_words = 0;
  GBuf h139o; size_t h139n = 0, h139t = 0;
  std::vector<std::pair<size_t, size_t>> h335b;
  std::vector<uint32_t> lid, glid; uint32_t nl = 0, gnl = 0;
};
static inline size_t lid916(const Side5& X, size_t j) { return g_h916 ? (size_t)X.lid[j] : j; }
static inline size_t glid916(const Side5& X, size_t j) { return g_h916 ? (size_t)X.glid[j] : j; }
struct Persist5 {
  Side5 S[2];
  GBuf scr, eb, ctr;
  id<MTLBuffer> bct = nil, bsob = nil;
  GBuf h139scr;
  std::vector<GBin5> gbins;
};


struct UqRep { const uint32_t* rec; uint32_t nrec; uint64_t h; };
struct PipeStep {
  int n = 0, NCH = 0;
  const evo::Cand* const* cands = nullptr; const uint8_t* valid = nullptr; const FlatKeySet* popset = nullptr;
  std::vector<uint8_t>* admit = nullptr; std::vector<Meta>* M = nullptr; StepBatch* SB = nullptr;
  FlatKeySet step, uidx;
  std::vector<FlatKeySet> stepSh, uidxSh; std::vector<std::vector<uint32_t>> shMap;
  std::vector<UqRep> rep;
  std::vector<int> fid;
  std::vector<uint8_t> fitted[2];
  std::vector<int> allids[2];
  std::vector<std::pair<double, double>> iv;
};


static int g_h139 = 2;
static id<MTLLibrary> g_h139L = nil;
static size_t g_h139_viol = 0;

static size_t h139_tgb(const SBin5& b) {
  size_t t = (size_t)b.ns * (b.s2 ? 96 : 768) * 4 + (b.s2 ? 2 : 16) * (size_t)b.nc * 4 + 56 * 4;
  return (t + 15) & ~(size_t)15;
}
static id<MTLComputePipelineState> h139_ps(Engine& E, int s2, int Q) {
  return E.pipe(g_h139L, s2 ? "rcscand_s2_q" + std::to_string(Q) : "rcscand_q" + std::to_string(Q));
}

static int g_h335 = 0; static id<MTLLibrary> g_h335L = nil;
static double g_h335_n = 0, g_h335_mis = 0, g_h335_kept = 0, g_h335_lr = 0, g_h335_bins = 0, g_h335_trees = 0;
static const int g_h802 = getenv("H802") ? atoi(getenv("H802")) : 0;
static const int g_h803 = getenv("H803") ? atoi(getenv("H803")) : 0;

static inline int h335_thr(int four, int Q) { const int G = (Q <= 5 || (four ? g_h803 : g_h802)) ? 6 : 3; return (four ? 192 : 384) / G; }
static id<MTLComputePipelineState> h335_ps(Engine& E, int four, int Q) {
  static id<MTLComputePipelineState> c[2][17]; if (Q < 1 || Q > 16) return nil;
  if (!c[four][Q]) { c[four][Q] = E.pipe(g_h335L, (four ? "rcscan4x_q" : "rcscan24_q") + std::to_string(Q));
    if ((g_h802 || g_h803) && c[four][Q]) { const int mt = (int)c[four][Q].maxTotalThreadsPerThreadgroup, th = h335_thr(four, Q);
      printf("H802/H803: %s q%d maxTotalThreadsPerThreadgroup %d, dispatch %d threads%s\n", four ? "rcscan4x" : "rcscan24", Q, mt, th, mt < th ? "  ** TOO SMALL **" : "");
      if (mt < th) { fprintf(stderr, "H802/H803: %s_q%d maxTotalThreadsPerThreadgroup %d < %d\n", four ? "rcscan4x" : "rcscan24", Q, mt, th); exit(1); } } }
  return c[four][Q];
}

static int g_h333 = 0; static GBuf g_h333hand[3][2], g_h333hoff[3][2]; static double g_h333_trees = 0, g_h333_mb = 0;
static inline size_t h333_hstr(int Q) { return (size_t)(Q + 1) * 48 + 20; }
static size_t h335_tgb(const SBin5& b, int four) { size_t t = (size_t)b.ns * (four ? 192 : 384) * 4 + (four ? 4 : 16) * (size_t)b.nc * 4 + 56 * 4; return (t + 15) & ~(size_t)15; }

static void h139_cmp(const Side5& X, int s, const Unique& UQ) {
  (void)UQ;
  const float* SO = (const float*)X.out.b.contents; const float* O2 = (const float*)X.h139o.b.contents;
  const uint32_t* bl = (const uint32_t*)X.bl.b.contents;
  long nt = 0, nbit = 0, nfin = 0, narg = 0; double mabs = 0, mrel = 0; int shown = 0;
  for (const SBin5& b : X.sbins) {
    if (h139_tgb(b) <= 32768) continue;
    for (size_t k = 0; k < b.cnt; k++) {
      const uint32_t tr = bl[b.off + k]; const float* a = SO + (size_t)tr * 16; const float* c = O2 + (size_t)tr * 16;
      nt++; bool bd = false; int ia = 0, ic = 0;
      for (int q = 0; q < 16; q++) {
        if (memcmp(&a[q], &c[q], 4) != 0) bd = true;
        const bool fa = std::isfinite(a[q]), fc = std::isfinite(c[q]);
        if (fa != fc) nfin++;
        else if (fa) { const double d = fabs((double)a[q] - (double)c[q]); if (d > mabs) mabs = d;
                       const double r = d / std::max(fabs((double)a[q]), 1e-30); if (r > mrel) mrel = r; }
        if (a[q] < a[ia]) ia = q;
        if (c[q] < c[ic]) ic = q;
      }
      if (ia != ic) narg++;
      if (bd) { nbit++; if (shown < 4) { shown++;
        printf("H139 cmp side %d tree %u: bit-differs; orig argmin %d = %.9g, dev argmin %d = %.9g\n", s, tr, ia, a[ia], ic, c[ic]); } }
    }
  }
  printf("H139 cmp side %d: over-limit TGs %ld, bit-differing trees %ld, finite-mismatch values %ld, max abs %.6g, max rel %.6g, argmin mismatches %ld\n",
         s, nt, nbit, nfin, mabs, mrel, narg);
}


#define private public
#include "gn_cpu_h940.hpp"
#undef private

static const int g_h779 = getenv("H779") ? atoi(getenv("H779")) : 0;
static const int g_h771 = getenv("H771") ? atoi(getenv("H771")) : 0;
static const int g_h770_raw = g_h771 ? 7 : (getenv("H770") ? atoi(getenv("H770")) : 0);
static const int g_h770 = (g_h770_raw & 5) ? (g_h770_raw | 5) : g_h770_raw;
static const int g_h770v = (getenv("H770V") && g_h770) ? atoi(getenv("H770V")) : 0;
static const int g_h777 = getenv("H777") ? atoi(getenv("H777")) : 0;
static const int g_ovl = getenv("STEP_OVL") ? atoi(getenv("STEP_OVL")) : 0;
static std::function<void(int)> g_ovl_hook;
static const int g_kesel = getenv("K_ESEL") ? atoi(getenv("K_ESEL")) : 0;
static int g_h774 = getenv("H774") ? atoi(getenv("H774")) : (g_ovl == 2 ? 1 : 0);
static std::vector<int> g_h774_cstart;
static const int g_h826 = getenv("H826") ? atoi(getenv("H826")) : 0;
static std::function<void(int)> g_h826_head;
static std::function<void()> g_h826_rest;
static const int g_h712 = getenv("H712") ? atoi(getenv("H712")) : 0;
static const int g_h713 = getenv("H713") ? atoi(getenv("H713")) : 0;
static const int g_h711 = getenv("H711") ? atoi(getenv("H711")) : 0;
static size_t ebt_bytes(int Q, int nth, size_t ntg, int tpg) {
  if (!g_h711) return 0;
  const size_t CH = nth <= 4 ? (size_t)(nth + 1) / 2 : 2;
  return ntg * (size_t)tpg * 2 * (size_t)Q * CH * 48 * 4; }
static const int g_h790 = getenv("H790") ? atoi(getenv("H790")) : 0;
#include "gn_cpu_h940.cpp"
static const int g_h250a = getenv("H250A") ? atoi(getenv("H250A")) : 0;
static int g_freet = 0;
static const int g_postlin = getenv("POST_LIN") ? atoi(getenv("POST_LIN")) : 0;
static int g_sf_axis = 0;


static const int g_sf_pin = getenv("SINGLE_FAST_PIN") ? atoi(getenv("SINGLE_FAST_PIN")) : 1;
static void sf_pin_driver(evo::Cand& c) {
  if (!g_sf_axis || !g_sf_pin) return;
  const int t = g_sf_axis == 1 ? 0 : 1;
  c.len[t] = 1; c.op[t][0] = evo::TVAR; c.c[t][0] = 0.0f;
}
static const int g_pinfix = getenv("PIN_FIX") ? atoi(getenv("PIN_FIX")) : 1;
static const int g_pincensus = getenv("PIN_CENSUS") ? atoi(getenv("PIN_CENSUS")) : 0;
static inline bool sf_is_driver(const evo::Cand& c) { const int t = g_sf_axis == 1 ? 0 : 1; return c.len[t] == 1 && c.op[t][0] == evo::TVAR; }

static inline bool postlin_on(int s) { return g_postlin || (g_pinfix && g_sf_axis && g_sf_pin && s == (g_sf_axis == 1 ? 0 : 1)); }
static const int g_h250at = getenv("H250AT") ? std::max(1, atoi(getenv("H250AT"))) : 6;
static double g_h250a_rg = 0, g_h250a_rc = 0;


static const int g_h300b = getenv("H300B") ? atoi(getenv("H300B")) : 0;
static double g_h300b_rg[16] = {0}, g_h300b_rc[16] = {0};
struct H250aItem { int side, Q, nth; float floor_; uint32_t tree, per; };
struct H770B;
struct Batch {
  int b0 = 0, b1 = 0, M = 0; Persist5* P = nullptr;
  std::vector<int> need[2];
  id<MTLCommandBuffer> cbS = nil, cbG = nil;
  unsigned long long t882S = 0, t882G = 0;
  double scan_gpu = 0, gn_gpu = 0;
  bool sdone = false, gdone = false, gbegun = false;
  bool egn = false;
  bool fdone = false;
  double fms = 0;
  int ci = 0;
  bool pdone = false;
  std::vector<size_t> h334_gc;
  std::shared_ptr<H770B> h770; std::vector<char> h770cpu[2];
  std::vector<H250aItem> h250_items; size_t h250_ngpu = 0; FitHandle h250_thr; double h250_cpu_ms = 0, h250_f = 0; bool h250_on = false;
};
static Persist5 PP[3];
static int g_chunks = 8;


static void nf_range(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  auto t0 = clk::now();
  for (Worker& Wk : WK(E)) { Wk.nf_only = 0; Wk.nf_fail = 0; }
  const int c1 = (b1 + 63) / 64;
  std::atomic<int> next(b0 / 64);
  PL(E).run((int)WK(E).size(), [&](int w) {
    Worker& Wk = WK(E)[w];
    for (;;) {
      const int c = next.fetch_add(1);
      if (c >= c1) break;
      const int a0 = c * 64, a1 = std::min(PS.n, a0 + 64);
      for (int i = a0; i < a1; i++) std::memcpy(&NFD.sc[i], PS.cands[i], sizeof(skn::Cand));
      auto a = clk::now();
      nf16::chunk(*Wk.H, NFD.sc.data(), PS.valid, a0, a1, NFD.blk[c], NFD.ch.data());
      Wk.nf_only += ms_since(a);
      const nf16::Block& B = NFD.blk[c];
      for (int i = a0; i < a1; i++) {
        if (!PS.valid[i]) continue;
        const nf16::Child& ch = NFD.ch[i];
        Meta& m = (*PS.M)[i];
        h890_setkey(m, B, ch);
        m.r0 = ch.r0; m.r1 = ch.r1; m.kh = ch.kh;
        Wk.nf_fail += (ch.ok[0] ? 0 : 1) + (ch.ok[1] ? 0 : 1);
      }
    }
  });
  for (Worker& Wk : WK(E)) { T.nf_only += Wk.nf_only; T.nf_fail += Wk.nf_fail; }
  T.nf += ms_since(t0);
}


static const char* g_h660d = getenv("H660D");
static const double g_h660ns = getenv("H660NS") ? atof(getenv("H660NS")) : 500.0;
static const int g_h660 = getenv("H660") ? atoi(getenv("H660")) : 0;
static const int g_h670 = getenv("H670") ? atoi(getenv("H670")) : 0;
static const int g_h673 = getenv("H673") ? atoi(getenv("H673")) : 0;
static const int g_h675 = getenv("H675") ? atoi(getenv("H675")) : 0;
static const int g_h675t = getenv("H675T") ? atoi(getenv("H675T")) : 4;
static thread_local bool tl_h675_job = false;
static void h670_join();
static std::vector<std::pair<int, float>> g_h670_res; static bool g_h670_have = false;
static std::vector<uint32_t> g_h670_cs;
static int g_h670_st = 0, g_h670_jst = -1;
static double g_h670_wait = 0; static long g_h670_n = 0;
static const double g_h660cpu = getenv("H660CPU") ? atof(getenv("H660CPU")) : 58.0;
static void h660_spin_ns(double ns) {
  if (ns <= 0) return;
  const auto t0 = std::chrono::steady_clock::now();
  while (std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() < ns) { }
}
static void h660_burn() {
  for (int i = 0; i < 8; i++) std::thread([] { pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0); h660_spin_ns(g_h660cpu * 1e6); }).detach();
}

static const int g_h680 = getenv("H680") ? atoi(getenv("H680")) : 0;
static const double g_h680cpu = getenv("H680CPU") ? atof(getenv("H680CPU")) : 466.0;
static const double g_h680g = getenv("H680G") ? atof(getenv("H680G")) : 55.0;
static id<MTLCommandQueue> g_h680q = nil; static id<MTLComputePipelineState> g_h680ps = nil; static id<MTLBuffer> g_h680o = nil;
static double g_h680_ms_per_iter = 0;
static std::vector<std::thread> g_h680_th; static std::vector<id<MTLCommandBuffer>> g_h680_cb;
static void h680_init(id<MTLDevice> dev) {
  if (g_h680ps) return;
  const char* src = "#include <metal_stdlib>\nusing namespace metal;\n"
    "kernel void h680_spin(device float* o [[buffer(0)]], constant uint& n [[buffer(1)]], uint gid [[thread_position_in_grid]]) {\n"
    "  float a = float(gid) * 1e-7f, b = 1.0001f;\n"
    "  for (uint i = 0; i < n; i++) { a = fma(a, b, 0.5f); b = fma(b, 0.9999f, 1e-7f); }\n"
    "  if (a == 12345.0f) o[0] = a + b; }\n";
  NSError* e = nil; id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:src] options:nil error:&e];
  if (!lib) { fprintf(stderr, "H680: spin kernel compile failed\n"); exit(1); }
  g_h680ps = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"h680_spin"] error:&e];
  g_h680q = [dev newCommandQueue]; g_h680o = h900_nbl(dev, __LINE__, 16, MTLResourceStorageModeShared);
  double best = 1e30; const uint32_t n0 = 1u << 14;
  for (int r = 0; r < 5; r++) {
    id<MTLCommandBuffer> cb = [g_h680q commandBuffer]; id<MTLComputeCommandEncoder> en = tl_enc(cb);
    [en setComputePipelineState:g_h680ps]; [en setBuffer:g_h680o offset:0 atIndex:0]; [en setBytes:&n0 length:4 atIndex:1];
    [en dispatchThreadgroups:MTLSizeMake(256, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)]; [en endEncoding];
    [cb commit]; [cb waitUntilCompleted]; best = std::min(best, (cb.GPUEndTime - cb.GPUStartTime) * 1e3); }
  g_h680_ms_per_iter = best / n0;
  printf("H680: spin calibration %.3f GPU-ms for %u iters (65536 threads)\n", best, n0);
}
static void h680_spin_cb(id<MTLCommandQueue> q, double gpu_ms) {
  const uint32_t n = (uint32_t)std::max(1.0, gpu_ms / g_h680_ms_per_iter);
  id<MTLCommandBuffer> cb = h182_cb(q); id<MTLComputeCommandEncoder> en = tl_enc(cb);
  [en setComputePipelineState:g_h680ps]; [en setBuffer:g_h680o offset:0 atIndex:0]; [en setBytes:&n length:4 atIndex:1];
  [en dispatchThreadgroups:MTLSizeMake(256, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)]; [en endEncoding]; [cb commit];
}
static void h860_launch(int Kn);
static void h680_inject(id<MTLDevice> dev, int Kn, bool last) {
  h680_init(dev);
  for (auto& t : g_h680_th) t.join(); g_h680_th.clear();
  const double core_ms = g_h680cpu / 8.0 / Kn;
  if (g_h680 == 2) h860_launch(Kn);
  else for (int i = 0; i < 8; i++) g_h680_th.emplace_back([core_ms] { pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0); h660_spin_ns(core_ms * 1e6); });
  const uint32_t n = (uint32_t)std::max(1.0, g_h680g / Kn / g_h680_ms_per_iter);
  id<MTLCommandBuffer> cb = h182_cb(g_h680q); id<MTLComputeCommandEncoder> en = tl_enc(cb);
  [en setComputePipelineState:g_h680ps]; [en setBuffer:g_h680o offset:0 atIndex:0]; [en setBytes:&n length:4 atIndex:1];
  [en dispatchThreadgroups:MTLSizeMake(256, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)]; [en endEncoding]; [cb commit];
  g_h680_cb.push_back(cb);
  if (last) { h301_mark("H680 residual join"); for (auto& t : g_h680_th) t.join(); g_h680_th.clear();
    for (id<MTLCommandBuffer> c : g_h680_cb) [c waitUntilCompleted]; g_h680_cb.clear(); }
}
static const int g_h743 = getenv("H743") ? atoi(getenv("H743")) : 0;
static void nf_range_h264(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  auto t0 = clk::now();
  for (Worker& Wk : WK(E)) { Wk.nf_only = 0; Wk.nf_fail = 0; }
  const int c1 = (b1 + 63) / 64;
  std::atomic<int> next(b0 / 64);
  PL(E).run((int)WK(E).size(), [&](int w) {
    Worker& Wk = WK(E)[w];
    for (;;) {
      if (nf308::g_h315a) {
        const int g = next.fetch_add(8);
        if (g >= c1) break;
        const int gc1 = std::min(c1, g + 8);
        const int a0 = g * 64, a1 = std::min(PS.n, gc1 * 64);
        auto a = clk::now();
        if (!Wk.S308) Wk.S308 = new nf308::St();
        if (g_h308 >= 2) nf308::group_h315<2>(*Wk.H, *Wk.S308, PS.cands, PS.valid, a0, a1, &NFD.blk[g], NFD.ch.data());
        else nf308::group_h315<1>(*Wk.H, *Wk.S308, PS.cands, PS.valid, a0, a1, &NFD.blk[g], NFD.ch.data());
        Wk.nf_only += ms_since(a);
        if (g_h660d && !strcmp(g_h660d, "frontpool")) h660_spin_ns(g_h660ns * (double)(a1 - a0));
        H890Tm h890tm;
        for (int c = g; c < gc1; c++) {
          const nf16::Block& B = NFD.blk[c];
          for (int i = c * 64, ie = std::min(PS.n, c * 64 + 64); i < ie; i++) {
            Meta& m = (*PS.M)[i];
            if (!PS.valid[i]) { m.key.clear(); m.r0 = 0; m.r1 = 0; m.kh = 0; continue; }
            const nf16::Child& ch = NFD.ch[i];
            h890_setkey(m, B, ch, &h890tm);
            m.r0 = ch.r0; m.r1 = ch.r1; m.kh = ch.kh;
            Wk.nf_fail += (ch.ok[0] ? 0 : 1) + (ch.ok[1] ? 0 : 1);
          }
        }
        continue;
      }
      const int c = next.fetch_add(1);
      if (c >= c1) break;
      const int a0 = c * 64, a1 = std::min(PS.n, a0 + 64);
      auto a = clk::now();
      if (g_h308) {
        if (!Wk.S308) Wk.S308 = new nf308::St();
        if (g_h308 >= 2) nf308::chunk_h308<2>(*Wk.H, *Wk.S308, PS.cands, PS.valid, a0, a1, NFD.blk[c], NFD.ch.data());
        else nf308::chunk_h308<1>(*Wk.H, *Wk.S308, PS.cands, PS.valid, a0, a1, NFD.blk[c], NFD.ch.data());
      } else
      nf16h264::chunk_h264(*Wk.H, PS.cands, PS.valid, a0, a1, NFD.blk[c], NFD.ch.data());
      Wk.nf_only += ms_since(a);
      const nf16::Block& B = NFD.blk[c];
      H890Tm h890tm;
      for (int i = a0; i < a1; i++) {
        Meta& m = (*PS.M)[i];
        if (!PS.valid[i]) { m.key.clear(); m.r0 = 0; m.r1 = 0; m.kh = 0; continue; }
        const nf16::Child& ch = NFD.ch[i];
        h890_setkey(m, B, ch, &h890tm);
        m.r0 = ch.r0; m.r1 = ch.r1; m.kh = ch.kh;
        Wk.nf_fail += (ch.ok[0] ? 0 : 1) + (ch.ok[1] ? 0 : 1);
      }
    }
  });
  for (Worker& Wk : WK(E)) { T.nf_only += Wk.nf_only; T.nf_fail += Wk.nf_fail; }
  T.nf += ms_since(t0);
}


static const int g_h896 = getenv("H896") ? atoi(getenv("H896")) : 0;
static const int g_h921 = getenv("H921") ? atoi(getenv("H921")) : 0;
static bool g_h921_have = false;
static std::vector<uint8_t> g_h921_keep; static std::vector<uint32_t> g_h921_code, g_h921_code2; static std::vector<std::vector<uint32_t>> g_h921_bk, g_h921_la; static std::vector<std::vector<std::pair<uint64_t, uint32_t>>> g_h921_ba; static std::vector<size_t> g_h921_bw;
struct H896Pop {
  std::vector<std::thread> th;
  void join() { for (auto& t : th) if (t.joinable()) t.join(); th.clear(); }
  ~H896Pop() { join(); }
};
static H896Pop* g_h896p = nullptr;
static inline void h896_join() { if (g_h896p && !g_h896p->th.empty()) { const auto tj = clk::now(); g_h896p->join(); if (g_h896 & 2) fprintf(stderr, "H896 join wait %.3f ms\n", ms_since(tj)); } }

static void dedup_range(PipeStep& PS, int b0, int b1, Times& T) {
  h896_join();
  auto t0 = clk::now();
  const std::vector<Meta>& M = *PS.M; std::vector<uint8_t>& admit = *PS.admit;
  for (int i = b0; i < b1; i++) {
    admit[i] = 0;
    if (!PS.valid[i]) continue;
    T.nchild += 1;
    const uint32_t kn = (uint32_t)M[i].key.size();
    if (PS.popset->contains(M[i].kh, M[i].key.data(), kn)) { T.rej_pop += 1; continue; }
    if (!PS.step.insert(M[i].kh, M[i].key.data(), kn)) { T.rej_within += 1; continue; }
    admit[i] = 1; T.admitted += 1;
  }
  T.dedup += ms_since(t0);
}


static int g_h98 = 1;


static void build_unique_range(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  auto t0 = clk::now();
  Unique& UQ = PS.SB->UQ; const std::vector<Meta>& M = *PS.M; const std::vector<uint8_t>& admit = *PS.admit;
  std::vector<int>& ct = PS.SB->ct;
  const uint32_t NONE = 0xFFFFFFFFu;
  std::vector<uint32_t> occ((size_t)(b1 - b0) * 2, NONE);
  const size_t rep_b = PS.rep.size();
  for (int i = b0; i < b1; i++) {
    if (!admit[i]) continue;
    const nf16::Child& ch = NFD.ch[i];
    const uint32_t* base = NFD.blk[i / 64].w.data();
    for (int t = 0; t < 2; t++) {
      if (!ch.ok[t]) continue;
      const uint32_t* rec = base + ch.roff[t]; const uint32_t nrec = ch.rlen[t];
      const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
      bool fresh = false;
      const uint32_t k = PS.uidx.find_or_add(h, rec, nrec, fresh);
      if (fresh) PS.rep.push_back(UqRep{rec, nrec, h});
      occ[(size_t)(i - b0) * 2 + t] = k;
    }
  }
  const size_t K = PS.rep.size() - rep_b;
  std::vector<Tree> ctree(K); std::vector<char> cok(K, 0); std::vector<TreeEnc> cenc(K);
  for (Worker& Wk : WK(E)) { Wk.cmp_only = 0; Wk.enc_only = 0; }
  std::atomic<size_t> nxt(0);
  PL(E).run((int)WK(E).size(), [&](int w) {
    Worker& Wk = WK(E)[w];
    for (;;) {
      const size_t b0k = nxt.fetch_add(32);
      if (b0k >= K) break;
      const size_t b1k = std::min(K, b0k + 32);

      uint64_t hb[32];
      {
        const uint32_t* wv[32]; size_t nkv[32]; int nbt = 0;
        for (size_t k = b0k; k < b1k; k++, nbt++) {
          View Vt(PS.rep[rep_b + k].rec); wv[nbt] = Vt.w; nkv[nbt] = (size_t)(Vt.theta - Vt.w);
          if (g_h98 & 1) hb[nbt] = cc14::hkey(Vt.w, nkv[nbt]);
        }
        cmi6::prefetch_batch(Wk.CH->C.S, Wk.CH->M, wv, nkv, nbt, (g_h98 & 1) ? hb : nullptr);
      }
      for (size_t k = b0k; k < b1k; k++) {
        View V(PS.rep[rep_b + k].rec);
        auto c0 = clk::now();
        double em = 0;
        cok[k] = compile_enc_sel(*Wk.CH, V, ctree[k], cenc[k], em, (g_h98 & 1) ? &hb[k - b0k] : nullptr) ? 1 : 0;
        Wk.cmp_only += ms_since(c0);
        Wk.enc_only += em;
      }
    }
  });
  for (size_t k = 0; k < K; k++) {
    if (!cok[k]) { T.cmp_fail += 1; PS.fid.push_back(-1); continue; }
    const UqRep& r = PS.rep[rep_b + k];
    Tree& tr = ctree[k]; tr.rec.assign(r.rec, r.rec + r.nrec);
    PS.fid.push_back((int)UQ.U.size());
    UQ.U.push_back(std::move(tr)); UQ.hash.push_back(r.h); UQ.enc.push_back(std::move(cenc[k]));
  }
  for (int i = b0; i < b1; i++) for (int t = 0; t < 2; t++) {
    const uint32_t o = occ[(size_t)(i - b0) * 2 + t];
    if (o != NONE) ct[(size_t)i * 2 + t] = PS.fid[o];
  }
  for (Worker& Wk : WK(E)) { T.cmp_only += Wk.cmp_only; T.enc_tree += Wk.enc_only; }
  T.comp += ms_since(t0);
}


static const int g_h842 = getenv("H842") ? atoi(getenv("H842")) : 0;
static inline uint32_t shard_h92(uint64_t h, int NS) { return (uint32_t)((h >> 32) % (uint64_t)NS); }

static void dedup_h92(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  h896_join();
  auto t0 = clk::now();
  const std::vector<Meta>& M = *PS.M; std::vector<uint8_t>& admit = *PS.admit;
  for (int i = b0; i < b1; i++) admit[i] = 0;
  const int NS = (int)PS.stepSh.size();
  struct C4 { double nchild = 0, rej_pop = 0, rej_within = 0, admitted = 0; };
  std::vector<C4> cnt((size_t)NS);

  const int n = b1 - b0;
  std::vector<uint32_t> shard_count((size_t)NS, 0);
  std::vector<uint8_t> shard_id((size_t)n, 0xFF);
  for (int i = b0; i < b1; i++) {
    if (!PS.valid[i]) continue;
    uint32_t s = shard_h92(M[i].kh, NS);
    shard_id[(size_t)(i - b0)] = (uint8_t)s;
    shard_count[(size_t)s]++;
  }

  std::vector<uint32_t> shard_offset((size_t)NS + 1, 0);
  for (int s = 0; s < NS; s++) {
    shard_offset[(size_t)s + 1] = shard_offset[(size_t)s] + shard_count[(size_t)s];
  }
  std::vector<uint32_t> shard_pos = shard_offset;
  std::vector<int> shard_cand((size_t)shard_offset[(size_t)NS]);
  for (int i = b0; i < b1; i++) {
    uint8_t s = shard_id[(size_t)(i - b0)];
    if (s != 0xFF) {
      shard_cand[shard_pos[(size_t)s]++] = i;
    }
  }

  PL(E).run(NS, [&](int s) {
    C4& C = cnt[(size_t)s]; FlatKeySet& SS = PS.stepSh[(size_t)s];
    const uint32_t start_idx = shard_offset[(size_t)s];
    const uint32_t end_idx = shard_offset[(size_t)s + 1];
    const bool pf842 = (g_h842 & 4) != 0;
    constexpr uint32_t D1 = 16;
    auto st1 = [&](uint32_t j) { const Meta& m = M[shard_cand[(size_t)j]]; __builtin_prefetch(&m, 0, 3); const char* kp = (const char*)m.key.data(); const size_t kb = m.key.size() * 4;
      for (size_t o = 0; o < kb; o += 64) __builtin_prefetch(kp + o, 0, 3); PS.popset->pf_slot(m.kh); };
    auto st2 = [&](uint32_t j) { PS.popset->pf_entry(M[shard_cand[(size_t)j]].kh); };
    auto st3 = [&](uint32_t j) { const uint64_t hh = M[shard_cand[(size_t)j]].kh; PS.popset->pf_arena(hh); SS.pf_slot(hh); };
    if (pf842) {
      for (uint32_t j = start_idx; j < end_idx && j < start_idx + D1; j++) st1(j);
      for (uint32_t j = start_idx; j < end_idx && j < start_idx + D1 / 2; j++) st2(j);
      for (uint32_t j = start_idx; j < end_idx && j < start_idx + D1 / 4; j++) st3(j);
    }
    for (uint32_t j = start_idx; j < end_idx; j++) {
      if (pf842) { if (j + D1 < end_idx) st1(j + D1); if (j + D1 / 2 < end_idx) st2(j + D1 / 2); if (j + D1 / 4 < end_idx) st3(j + D1 / 4); }
      const int i = shard_cand[(size_t)j];
      const uint64_t h = M[i].kh;
      C.nchild += 1;
      const uint32_t kn = (uint32_t)M[i].key.size();
      if (PS.popset->contains(h, M[i].key.data(), kn)) { C.rej_pop += 1; continue; }
      if (!SS.insert(h, M[i].key.data(), kn)) { C.rej_within += 1; continue; }
      admit[i] = 1; C.admitted += 1;
    }
  });
  for (const C4& C : cnt) { T.nchild += C.nchild; T.rej_pop += C.rej_pop; T.rej_within += C.rej_within; T.admitted += C.admitted; }
  T.dedup += ms_since(t0);
}


static const int g_h864 = getenv("H864") ? atoi(getenv("H864")) : 0;
static inline void h864_lap(Tick& tk, double& f) { if (g_h864 & 1) f += tk.lap(); }

static const int g_h879 = getenv("H879") ? atoi(getenv("H879")) : 0;
static int g_h879_lanes = 8; static std::mutex g_h879_m; static double g_h879_wb[64]; static double g_h879_wc[64]; static long g_h879_wn[64];
static inline double h879_cpu_ms() { return (double)clock_gettime_nsec_np(CLOCK_THREAD_CPUTIME_ID) * 1e-6; }
static inline uint64_t h864_mix(uint64_t h, uint64_t v) { return splitmix64(h ^ v); }
template <class V> static inline uint64_t h864_vec(uint64_t h, const V& v) {
  h = h864_mix(h, (uint64_t)v.size());
  const unsigned char* p = (const unsigned char*)v.data(); const size_t nb = v.size() * sizeof(v[0]);
  for (size_t i = 0; i + 8 <= nb; i += 8) { uint64_t w; std::memcpy(&w, p + i, 8); h = h864_mix(h, w); }
  for (size_t i = nb & ~(size_t)7; i < nb; i++) h = h864_mix(h, p[i]);
  return h;
}
static uint64_t h864_fp(const Tree& t, const cmi6::TreeEnc& e) {
  uint64_t h = 0x864;
  h = h864_mix(h, (uint64_t)(uint32_t)t.Q); h = h864_mix(h, (uint64_t)(uint32_t)t.nth); h = h864_mix(h, (uint64_t)(uint32_t)t.nS); h = h864_mix(h, (uint64_t)(uint32_t)t.nslot);
  h = h864_vec(h, t.prog); h = h864_vec(h, t.sp); h = h864_vec(h, t.th); h = h864_vec(h, t.rec);
  for (int i = 0; i < 36; i++) h = h864_mix(h, t.segs[i]);
  h = h864_mix(h, (uint64_t)(uintptr_t)t.progp); h = h864_mix(h, (uint64_t)(uintptr_t)t.segsp); h = h864_mix(h, t.nprogp);
  h = h864_mix(h, e.s_ok); h = h864_mix(h, e.g_ok);
  for (int i = 0; i < 36; i++) { h = h864_mix(h, e.sg[i]); h = h864_mix(h, e.gg[i]); }
  h = h864_mix(h, e.ghw); h = h864_mix(h, e.gper);
  h = h864_vec(h, e.sp4); h = h864_vec(h, e.gp4); h = h864_vec(h, e.tk);
  h = h864_mix(h, (uint64_t)(uintptr_t)e.sp4p); h = h864_mix(h, (uint64_t)(uintptr_t)e.gp4p); h = h864_mix(h, e.nsp4p); h = h864_mix(h, e.ngp4p);
  h = h864_mix(h, (uint64_t)(uintptr_t)e.tkp); h = h864_mix(h, e.ntkp); h = h864_mix(h, (uint64_t)(uintptr_t)e.gtp);
  return h;
}
static void build_unique_h92(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  auto t0 = clk::now();
  Tick tk864;
  Unique& UQ = PS.SB->UQ; const std::vector<Meta>& M = *PS.M; const std::vector<uint8_t>& admit = *PS.admit;
  std::vector<int>& ct = PS.SB->ct;
  const uint32_t NONE = 0xFFFFFFFFu;
  const int NS = (int)PS.uidxSh.size();
  const size_t NE = (size_t)(b1 - b0) * 2;
  std::vector<uint32_t> evLe(NE, NONE), evSh(NE, 0), occ(NE, NONE);
  std::vector<uint8_t> evFresh(NE, 0);
  const size_t rep_b = PS.rep.size();
  auto tu = clk::now();

  std::vector<uint32_t> ev_shard_count((size_t)NS, 0);
  std::vector<uint8_t> ev_shard_id(NE, 0xFF);
  for (int i = b0; i < b1; i++) {
    if (!admit[i]) continue;
    const nf16::Child& ch = NFD.ch[i];
    for (int t = 0; t < 2; t++) {
      if (!ch.ok[t]) continue;
      const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
      const uint32_t s = shard_h92(h, NS);
      const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
      ev_shard_id[e] = (uint8_t)s;
      ev_shard_count[(size_t)s]++;
    }
  }

  std::vector<uint32_t> ev_shard_offset((size_t)NS + 1, 0);
  for (int s = 0; s < NS; s++) {
    ev_shard_offset[(size_t)s + 1] = ev_shard_offset[(size_t)s] + ev_shard_count[(size_t)s];
  }
  std::vector<uint32_t> ev_shard_pos = ev_shard_offset;
  std::vector<uint32_t> ev_shard_events((size_t)ev_shard_offset[(size_t)NS]);
  for (size_t e = 0; e < NE; e++) {
    uint8_t s = ev_shard_id[e];
    if (s != 0xFF) {
      ev_shard_events[ev_shard_pos[(size_t)s]++] = (uint32_t)e;
    }
  }

  h864_lap(tk864, T.u_p1);
  PL(E).run(NS, [&](int s) {
    FlatKeySet& U = PS.uidxSh[(size_t)s];
    const uint32_t start_idx = ev_shard_offset[(size_t)s];
    const uint32_t end_idx = ev_shard_offset[(size_t)s + 1];
    for (uint32_t j = start_idx; j < end_idx; j++) {
      const size_t e = (size_t)ev_shard_events[(size_t)j];
      const int i = b0 + (int)(e / 2);
      const int t = (int)(e % 2);
      const nf16::Child& ch = NFD.ch[i];
      const uint32_t* base = NFD.blk[i / 64].w.data();
      const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
      bool fresh = false;
      evLe[e] = U.find_or_add(h, base + ch.roff[t], ch.rlen[t], fresh);
      evFresh[e] = fresh ? 1 : 0; evSh[e] = (uint32_t)s;
    }
  });
  h864_lap(tk864, T.u_find);

  for (int i = b0; i < b1; i++) {
    if (!admit[i]) continue;
    const nf16::Child& ch = NFD.ch[i];
    const uint32_t* base = NFD.blk[i / 64].w.data();
    for (int t = 0; t < 2; t++) {
      if (!ch.ok[t]) continue;
      const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
      const uint32_t s = evSh[e], le = evLe[e];
      std::vector<uint32_t>& sm = PS.shMap[s];
      if (evFresh[e]) {
        if (sm.size() != (size_t)le) { fprintf(stderr, "H92: shard %u local entry %u out of order\n", s, le); exit(1); }
        const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
        sm.push_back((uint32_t)PS.rep.size());
        PS.rep.push_back(UqRep{base + ch.roff[t], ch.rlen[t], h});
      }
      occ[e] = sm[le];
    }
  }
  h864_lap(tk864, T.u_pref);
  T.uq += ms_since(tu);
  const size_t K = PS.rep.size() - rep_b;
  const bool b2 = (g_h864 & 2) != 0, b4 = (g_h864 & 4) != 0, b8 = (g_h864 & 8) != 0;

  const size_t U0 = UQ.U.size();
  std::vector<Tree> ctree; std::vector<TreeEnc> cenc; std::vector<char> cok(K, 0);
  std::vector<uint64_t> vfp;
  if (b8) vfp.assign((K + 63) / 64, 0);
  if (b4) { UQ.U.resize(U0 + K); UQ.enc.resize(U0 + K); UQ.hash.resize(U0 + K); }
  else { ctree.resize(K); cenc.resize(K); }
  h864_lap(tk864, T.u_ctor);
  for (Worker& Wk : WK(E)) { Wk.cmp_only = 0; Wk.enc_only = 0; }
  const int h879 = g_h879;
  const size_t BS = (h879 & 12) == 4 ? 16 : (h879 & 12) == 8 ? 64 : (h879 & 12) == 12 ? 8 : 32;
  const size_t NBL = (K + BS - 1) / BS;
  std::vector<uint32_t> bord, rord;
  auto to0 = clk::now();
  if (h879 & 2) {
    std::vector<uint64_t> bc(NBL, 0);
    for (size_t k = 0; k < K; k++) bc[k / BS] += PS.rep[rep_b + k].nrec;
    bord.resize(NBL); for (size_t q = 0; q < NBL; q++) bord[q] = (uint32_t)q;
    std::stable_sort(bord.begin(), bord.end(), [&](uint32_t x, uint32_t y) { return bc[x] > bc[y]; });
  }
  if (h879 & 16) {
    uint32_t mx = 0; for (size_t k = 0; k < K; k++) mx = std::max(mx, PS.rep[rep_b + k].nrec);
    std::vector<uint32_t> cnt((size_t)mx + 2, 0);
    for (size_t k = 0; k < K; k++) cnt[(size_t)(mx - PS.rep[rep_b + k].nrec) + 1]++;
    for (size_t q = 1; q < cnt.size(); q++) cnt[q] += cnt[q - 1];
    rord.resize(K);
    for (size_t k = 0; k < K; k++) rord[cnt[(size_t)(mx - PS.rep[rep_b + k].nrec)]++] = (uint32_t)k;
  }
  const double ord_ms = ms_since(to0);
  const size_t NW = WK(E).size();
  std::vector<double> d_st(NW, 0), d_en(NW, 0), d_cpu(NW, 0), d_busy(NW, 0); std::vector<long> d_n(NW, 0);
  const bool dg = (h879 & 1) != 0; g_h879_lanes = (int)NW;
  auto tl0 = clk::now();
  std::atomic<size_t> nxt(0);
  PL(E).run((int)NW, [&](int w) {
    Worker& Wk = WK(E)[w];
    const double cpu0 = dg ? h879_cpu_ms() : 0, busy0 = Wk.cmp_only, ts = dg ? ms_since(tl0) : 0; long nclaim = 0;
    for (;;) {
      size_t kk[64]; size_t nb = 0;
      if (bord.empty()) {
        const size_t b0k = nxt.fetch_add(BS);
        if (b0k >= K) break;
        const size_t b1k = std::min(K, b0k + BS);
        for (size_t k = b0k; k < b1k; k++) kk[nb++] = rord.empty() ? k : (size_t)rord[k];
      } else {
        const size_t q = nxt.fetch_add(1);
        if (q >= NBL) break;
        const size_t b0k = (size_t)bord[q] * BS, b1k = std::min(K, b0k + BS);
        for (size_t k = b0k; k < b1k; k++) kk[nb++] = k;
      }
      nclaim++;

      uint64_t hb[64];
      {
        const uint32_t* wv[64]; size_t nkv[64];
        for (size_t j = 0; j < nb; j++) {
          View Vt(PS.rep[rep_b + kk[j]].rec); wv[j] = Vt.w; nkv[j] = (size_t)(Vt.theta - Vt.w);
          if (g_h98 & 1) hb[j] = cc14::hkey(Vt.w, nkv[j]);
        }
        for (size_t j0 = 0; j0 < nb; j0 += 32) cmi6::prefetch_batch(Wk.CH->C.S, Wk.CH->M, wv + j0, nkv + j0, (int)std::min<size_t>(32, nb - j0), (g_h98 & 1) ? hb + j0 : nullptr);
      }
      if (g_h301r && g_t4) h308d_claim(nb);
      for (size_t j = 0; j < nb; j++) {
        const size_t k = kk[j];
        View V(PS.rep[rep_b + k].rec);
        auto c0 = clk::now();
        double em = 0;
        Tree& trk = b4 ? UQ.U[U0 + k] : ctree[k]; TreeEnc& enk = b4 ? UQ.enc[U0 + k] : cenc[k];
        if (g_h301r) h301r_take(trk, enk);
        cok[k] = compile_enc_sel(*Wk.CH, V, trk, enk, em, (g_h98 & 1) ? &hb[j] : nullptr) ? 1 : 0;
        Wk.cmp_only += ms_since(c0);
        Wk.enc_only += em;
        if (b4) {
          const UqRep& r = PS.rep[rep_b + k];
          UQ.hash[U0 + k] = r.h;
          if (cok[k]) trk.rec.assign(r.rec, r.rec + r.nrec);
        } else if (b8 && cok[k]) { const UqRep& r = PS.rep[rep_b + k]; trk.rec.assign(r.rec, r.rec + r.nrec); }
        if (b8 && (k & 63) == 0) vfp[k >> 6] = cok[k] ? h864_fp(trk, enk) : 0;
      }
    }
    if (dg) { d_st[w] = ts; d_en[w] = ms_since(tl0); d_cpu[w] = h879_cpu_ms() - cpu0; d_busy[w] = Wk.cmp_only - busy0; d_n[w] = nclaim; }
  });
  if (dg) {
    const double wall = ms_since(tl0); double sc = 0, sb = 0, lat = 0, tail = 0;
    for (size_t w = 0; w < NW; w++) { sc += d_cpu[w]; sb += d_busy[w]; if (d_n[w]) lat += d_st[w]; tail += wall - d_en[w]; }
    T.c_runs += 1; T.c_recs += (double)K; T.c_wall += wall; T.c_busy += sb; T.c_cpu += sc; T.c_lat += lat; T.c_tail += tail; T.c_ord += ord_ms;
    std::lock_guard<std::mutex> lg(g_h879_m);
    for (size_t w = 0; w < NW && w < 64; w++) { g_h879_wb[w] += d_busy[w]; g_h879_wc[w] += d_cpu[w]; g_h879_wn[w] += d_n[w]; }
  }
  h864_lap(tk864, T.u_comp);
  const size_t fidB = PS.fid.size();
  size_t nokK = 0;
  if (!b4 && !b2) {
    for (size_t k = 0; k < K; k++) {
      if (!cok[k]) { T.cmp_fail += 1; PS.fid.push_back(-1); continue; }
      const UqRep& r = PS.rep[rep_b + k];
      Tree& tr = ctree[k]; tr.rec.assign(r.rec, r.rec + r.nrec);
      PS.fid.push_back((int)UQ.U.size());
      UQ.U.push_back(std::move(tr)); UQ.hash.push_back(r.h); UQ.enc.push_back(std::move(cenc[k]));
    }
  } else if (!b4) {

    PS.fid.resize(fidB + K);
    const size_t BKr = 256, NBk = (K + BKr - 1) / BKr;
    std::vector<size_t> okc(NBk, 0);
    PL(E).run((int)NBk, [&](int bk) {
      const size_t k1 = std::min(K, ((size_t)bk + 1) * BKr);
      size_t c = 0;
      for (size_t k = (size_t)bk * BKr; k < k1; k++) c += cok[k] ? 1 : 0;
      okc[(size_t)bk] = c;
    });
    std::vector<size_t> okb(NBk + 1, 0);
    for (size_t bk = 0; bk < NBk; bk++) okb[bk + 1] = okb[bk] + okc[bk];
    const size_t nok = okb[NBk];
    UQ.U.resize(U0 + nok); UQ.hash.resize(U0 + nok); UQ.enc.resize(U0 + nok);
    PL(E).run((int)NBk, [&](int bk) {
      const size_t k1 = std::min(K, ((size_t)bk + 1) * BKr);
      size_t id = U0 + okb[(size_t)bk];
      for (size_t k = (size_t)bk * BKr; k < k1; k++) {
        if (!cok[k]) { PS.fid[fidB + k] = -1; continue; }
        const UqRep& r = PS.rep[rep_b + k];
        ctree[k].rec.assign(r.rec, r.rec + r.nrec);
        PS.fid[fidB + k] = (int)id;
        UQ.U[id] = std::move(ctree[k]); UQ.hash[id] = r.h; UQ.enc[id] = std::move(cenc[k]);
        id++;
      }
    });
    T.cmp_fail += (double)(K - nok);
  } else {

    PS.fid.resize(fidB + K);
    for (size_t k = 0; k < K; k++) nokK += cok[k] ? 1 : 0;
    if (nokK == K) {
      for (size_t k = 0; k < K; k++) PS.fid[fidB + k] = (int)(U0 + k);
    } else {
      size_t wr = U0;
      for (size_t k = 0; k < K; k++) {
        if (!cok[k]) { T.cmp_fail += 1; PS.fid[fidB + k] = -1; continue; }
        if (wr != U0 + k) { UQ.U[wr] = std::move(UQ.U[U0 + k]); UQ.hash[wr] = UQ.hash[U0 + k]; UQ.enc[wr] = std::move(UQ.enc[U0 + k]); }
        PS.fid[fidB + k] = (int)wr; wr++;
      }
      UQ.U.resize(wr); UQ.hash.resize(wr); UQ.enc.resize(wr);
    }
  }
  h864_lap(tk864, T.u_ids);
  if (b2) {
    const int BK = 256, NB = (b1 - b0 + BK - 1) / BK;
    PL(E).run(NB, [&](int bi) {
      const int i1 = std::min(b1, b0 + (bi + 1) * BK);
      for (int i = b0 + bi * BK; i < i1; i++) for (int t = 0; t < 2; t++) {
        const uint32_t o = occ[(size_t)(i - b0) * 2 + (size_t)t];
        if (o != NONE) ct[(size_t)i * 2 + (size_t)t] = PS.fid[o];
      }
    });
  } else
  for (int i = b0; i < b1; i++) for (int t = 0; t < 2; t++) {
    const uint32_t o = occ[(size_t)(i - b0) * 2 + (size_t)t];
    if (o != NONE) ct[(size_t)i * 2 + (size_t)t] = PS.fid[o];
  }
  h864_lap(tk864, T.u_ct);
  if (b8) {
    size_t id = U0, bad = 0;
    for (size_t k = 0; k < K; k++) {
      const int want = cok[k] ? (int)id : -1;
      if (PS.fid[fidB + k] != want) { bad++; if (bad < 4) fprintf(stderr, "H864 verify: fid[%zu] %d want %d\n", k, PS.fid[fidB + k], want); }
      if (cok[k]) {
        const UqRep& r = PS.rep[rep_b + k];
        if (UQ.hash[id] != r.h || UQ.U[id].rec.size() != r.nrec || std::memcmp(UQ.U[id].rec.data(), r.rec, r.nrec * sizeof(uint32_t)) != 0) { bad++; if (bad < 4) fprintf(stderr, "H864 verify: hash/rec at id %zu (k %zu)\n", id, k); }
        if ((k & 63) == 0 && h864_fp(UQ.U[id], UQ.enc[id]) != vfp[k >> 6]) { bad++; if (bad < 4) fprintf(stderr, "H864 verify: Tree/TreeEnc fingerprint at id %zu (k %zu)\n", id, k); }
        id++;
      }
    }
    if (UQ.U.size() != id || UQ.hash.size() != id || UQ.enc.size() != id) { bad++; fprintf(stderr, "H864 verify: sizes %zu %zu %zu want %zu\n", UQ.U.size(), UQ.hash.size(), UQ.enc.size(), id); }
    for (int i = b0; i < b1; i++) for (int t = 0; t < 2; t++) {
      const uint32_t o = occ[(size_t)(i - b0) * 2 + (size_t)t];
      if (o != NONE && ct[(size_t)i * 2 + (size_t)t] != PS.fid[o]) { bad++; if (bad < 4) fprintf(stderr, "H864 verify: ct[%d,%d]\n", i, t); }
    }
    if (bad) { fprintf(stderr, "H864 verify: %zu mismatches (chunk %d..%d)\n", bad, b0, b1); exit(1); }
  }
  for (Worker& Wk : WK(E)) { T.cmp_only += Wk.cmp_only; T.enc_tree += Wk.enc_only; }
  T.comp += ms_since(t0);
}


static void build_unique_h257(Engine& E, PipeStep& PS, int b0, int b1, Times& T) {
  auto t0 = clk::now();
  Tick tk864;
  Unique& UQ = PS.SB->UQ; const std::vector<Meta>& M = *PS.M; const std::vector<uint8_t>& admit = *PS.admit;
  std::vector<int>& ct = PS.SB->ct;
  const uint32_t NONE = 0xFFFFFFFFu;
  const int NS = (int)PS.uidxSh.size();
  const size_t NE = (size_t)(b1 - b0) * 2;
  std::vector<uint32_t> evLe(NE, NONE), evSh(NE, 0), occ(NE, NONE);
  std::vector<uint8_t> evFresh(NE, 0);
  const size_t rep_b = PS.rep.size();
  const int BK = 256, NB = (b1 - b0 + BK - 1) / BK;
  auto tu = clk::now();
  if (g_h264) {


    std::vector<uint8_t> evS1(NE, 0);
    std::vector<uint64_t> evH(NE, 0);
    PL(E).run(NB, [&](int bi) {
      const int i1 = std::min(b1, b0 + (bi + 1) * BK);
      for (int i = b0 + bi * BK; i < i1; i++) {
        if (!admit[i]) continue;
        const nf16::Child& ch = NFD.ch[i];
        for (int t = 0; t < 2; t++) {
          if (!ch.ok[t]) continue;
          const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
          const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
          evH[e] = h; evS1[e] = (uint8_t)(shard_h92(h, NS) + 1);
        }
      }
    });
    PL(E).run(NS, [&](int s) {
      FlatKeySet& U = PS.uidxSh[(size_t)s];
      for (size_t e = 0; e < NE; e++) {
        if (evS1[e] != (uint8_t)(s + 1)) continue;
        const int i = b0 + (int)(e >> 1); const int t = (int)(e & 1);
        const nf16::Child& ch = NFD.ch[i];
        const uint32_t* base = NFD.blk[i / 64].w.data();
        bool fresh = false;
        evLe[e] = U.find_or_add(evH[e], base + ch.roff[t], ch.rlen[t], fresh);
        evFresh[e] = fresh ? 1 : 0; evSh[e] = (uint32_t)s;
      }
    });
  } else

  PL(E).run(NS, [&](int s) {
    FlatKeySet& U = PS.uidxSh[(size_t)s];
    for (int i = b0; i < b1; i++) {
      if (!admit[i]) continue;
      const nf16::Child& ch = NFD.ch[i];
      const uint32_t* base = NFD.blk[i / 64].w.data();
      for (int t = 0; t < 2; t++) {
        if (!ch.ok[t]) continue;
        const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
        if (shard_h92(h, NS) != (uint32_t)s) continue;
        const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
        bool fresh = false;
        evLe[e] = U.find_or_add(h, base + ch.roff[t], ch.rlen[t], fresh);
        evFresh[e] = fresh ? 1 : 0; evSh[e] = (uint32_t)s;
      }
    }
  });
  h864_lap(tk864, T.u_find);

  std::vector<size_t> cntF((size_t)NB, 0), cntS((size_t)NB * NS, 0);
  PL(E).run(NB, [&](int bi) {
    const int i1 = std::min(b1, b0 + (bi + 1) * BK);
    for (int i = b0 + bi * BK; i < i1; i++) {
      if (!admit[i]) continue;
      for (int t = 0; t < 2; t++) {
        const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
        if (!evFresh[e]) continue;
        cntF[(size_t)bi]++; cntS[(size_t)bi * NS + evSh[e]]++;
      }
    }
  });

  std::vector<size_t> offF((size_t)NB + 1, 0), shBeg((size_t)NS, 0), shBase((size_t)NB * NS, 0), shTot((size_t)NS, 0);
  for (int bi = 0; bi < NB; bi++) offF[(size_t)bi + 1] = offF[(size_t)bi] + cntF[(size_t)bi];
  for (int s = 0; s < NS; s++) {
    shBeg[(size_t)s] = PS.shMap[(size_t)s].size();
    size_t acc = shBeg[(size_t)s];
    for (int bi = 0; bi < NB; bi++) { shBase[(size_t)bi * NS + s] = acc; acc += cntS[(size_t)bi * NS + s]; }
    shTot[(size_t)s] = acc - shBeg[(size_t)s];
  }
  PS.rep.resize(rep_b + offF[(size_t)NB]);
  for (int s = 0; s < NS; s++) PS.shMap[(size_t)s].resize(shBeg[(size_t)s] + shTot[(size_t)s]);

  PL(E).run(NB, [&](int bi) {
    size_t rf = rep_b + offF[(size_t)bi];
    std::vector<size_t> nx((size_t)NS);
    for (int s = 0; s < NS; s++) nx[(size_t)s] = shBase[(size_t)bi * NS + s];
    const int i1 = std::min(b1, b0 + (bi + 1) * BK);
    for (int i = b0 + bi * BK; i < i1; i++) {
      if (!admit[i]) continue;
      const nf16::Child& ch = NFD.ch[i];
      const uint32_t* base = NFD.blk[i / 64].w.data();
      for (int t = 0; t < 2; t++) {
        if (!ch.ok[t]) continue;
        const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
        if (!evFresh[e]) continue;
        const uint32_t s = evSh[e];
        if ((size_t)evLe[e] != nx[s]) { fprintf(stderr, "H257: shard %u local entry %u out of order\n", s, evLe[e]); exit(1); }
        const uint64_t h = t == 0 ? M[i].r0 : M[i].r1;
        PS.shMap[s][evLe[e]] = (uint32_t)rf;
        PS.rep[rf] = UqRep{base + ch.roff[t], ch.rlen[t], h};
        rf++; nx[s]++;
      }
    }
  });

  PL(E).run(NB, [&](int bi) {
    const int i1 = std::min(b1, b0 + (bi + 1) * BK);
    for (int i = b0 + bi * BK; i < i1; i++) {
      if (!admit[i]) continue;
      const nf16::Child& ch = NFD.ch[i];
      for (int t = 0; t < 2; t++) {
        if (!ch.ok[t]) continue;
        const size_t e = (size_t)(i - b0) * 2 + (size_t)t;
        occ[e] = PS.shMap[evSh[e]][evLe[e]];
      }
    }
  });
  h864_lap(tk864, T.u_pref);
  T.uq += ms_since(tu);
  const size_t K = PS.rep.size() - rep_b;

  std::vector<Tree> ctree(K); std::vector<char> cok(K, 0); std::vector<TreeEnc> cenc(K);
  h864_lap(tk864, T.u_ctor);
  for (Worker& Wk : WK(E)) { Wk.cmp_only = 0; Wk.enc_only = 0; }
  std::atomic<size_t> nxt(0);
  PL(E).run((int)WK(E).size(), [&](int w) {
    Worker& Wk = WK(E)[w];
    for (;;) {
      const size_t b0k = nxt.fetch_add(32);
      if (b0k >= K) break;
      const size_t b1k = std::min(K, b0k + 32);

      uint64_t hb[32];
      {
        const uint32_t* wv[32]; size_t nkv[32]; int nbt = 0;
        for (size_t k = b0k; k < b1k; k++, nbt++) {
          View Vt(PS.rep[rep_b + k].rec); wv[nbt] = Vt.w; nkv[nbt] = (size_t)(Vt.theta - Vt.w);
          if (g_h98 & 1) hb[nbt] = cc14::hkey(Vt.w, nkv[nbt]);
        }
        cmi6::prefetch_batch(Wk.CH->C.S, Wk.CH->M, wv, nkv, nbt, (g_h98 & 1) ? hb : nullptr);
      }
      if (g_h301r && g_t4) h308d_claim(b1k - b0k);
      for (size_t k = b0k; k < b1k; k++) {
        View V(PS.rep[rep_b + k].rec);
        auto c0 = clk::now();
        double em = 0;
        if (g_h301r) h301r_take(ctree[k], cenc[k]);
        cok[k] = compile_enc_sel(*Wk.CH, V, ctree[k], cenc[k], em, (g_h98 & 1) ? &hb[k - b0k] : nullptr) ? 1 : 0;
        Wk.cmp_only += ms_since(c0);
        Wk.enc_only += em;
      }
    }
  });
  h864_lap(tk864, T.u_comp);

  const size_t fidB = PS.fid.size();
  PS.fid.resize(fidB + K);
  const size_t BKr = 256, NBk = (K + BKr - 1) / BKr;
  std::vector<size_t> okc(NBk, 0);
  PL(E).run((int)NBk, [&](int bk) {
    const size_t k1 = std::min(K, ((size_t)bk + 1) * BKr);
    size_t c = 0;
    for (size_t k = (size_t)bk * BKr; k < k1; k++) c += cok[k] ? 1 : 0;
    okc[(size_t)bk] = c;
  });
  std::vector<size_t> okb(NBk + 1, 0);
  for (size_t bk = 0; bk < NBk; bk++) okb[bk + 1] = okb[bk] + okc[bk];
  const size_t U0 = UQ.U.size(), nok = okb[NBk];
  UQ.U.resize(U0 + nok); UQ.hash.resize(U0 + nok); UQ.enc.resize(U0 + nok);
  PL(E).run((int)NBk, [&](int bk) {
    const size_t k1 = std::min(K, ((size_t)bk + 1) * BKr);
    size_t id = U0 + okb[(size_t)bk];
    for (size_t k = (size_t)bk * BKr; k < k1; k++) {
      if (!cok[k]) { PS.fid[fidB + k] = -1; continue; }
      const UqRep& r = PS.rep[rep_b + k];
      ctree[k].rec.assign(r.rec, r.rec + r.nrec);
      PS.fid[fidB + k] = (int)id;
      UQ.U[id] = std::move(ctree[k]); UQ.hash[id] = r.h; UQ.enc[id] = std::move(cenc[k]);
      id++;
    }
  });
  T.cmp_fail += (double)(K - nok);
  h864_lap(tk864, T.u_ids);

  PL(E).run(NB, [&](int bi) {
    const int i1 = std::min(b1, b0 + (bi + 1) * BK);
    for (int i = b0 + bi * BK; i < i1; i++) for (int t = 0; t < 2; t++) {
      const uint32_t o = occ[(size_t)(i - b0) * 2 + (size_t)t];
      if (o != NONE) ct[(size_t)i * 2 + (size_t)t] = PS.fid[o];
    }
  });
  h864_lap(tk864, T.u_ct);
  for (Worker& Wk : WK(E)) { T.cmp_only += Wk.cmp_only; T.enc_tree += Wk.enc_only; }
  T.comp += ms_since(t0);
}

static constexpr uint32_t H99_NONE = 0xFFFFFFFFu;
static uint32_t h99_find(uint64_t hr, uint32_t t, const uint32_t* k, uint32_t n);
static void h99_lookup_count();
static void h99_hit_count();
static int g_h99 = 0;
static void h99_apply(PipeStep& PS, int t, int id, uint32_t e, Times& T);
static void h99_put(PipeStep& PS, int t, int j, const FitOut& o, uint32_t steps);


static int g_h242 = 0;
static int g_h413 = 0;
static const bool g_h413dbg = getenv("H413DBG") != nullptr;
static uint64_t g_h413_hk = 0, g_h413_ik = 0; static long g_h413_h0 = 0; static uint64_t g_h413_f[6] = {0, 0, 0, 0, 0, 0}; static std::atomic<long> g_h413_tail{0}, g_h413_la{0}, g_h413_lt{0};

static int g_h242d_bal = getenv("H242D_BAL") ? atoi(getenv("H242D_BAL")) : 0;
static const int g_h410ns = getenv("H410NS") ? atoi(getenv("H410NS")) : 0;
struct StepBatch; static StepBatch* g_h415_spare = nullptr; static double g_h415_join_ms = 0;
static const int g_h415 = getenv("H415") ? atoi(getenv("H415")) : 0;
static const int g_h405 = getenv("H405") ? atoi(getenv("H405")) : 0;
static const double g_h404gn = getenv("H404GN") ? atof(getenv("H404GN")) : -1.0;
static double g_bal_s = 1.0, g_h242d_union = 0; static int g_h242d_steps = 0;

static const int g_h648 = getenv("H648") ? atoi(getenv("H648")) : 0;
static const double g_h648t = getenv("H648T") ? atof(getenv("H648T")) : 10.0, g_h648k = getenv("H648K") ? atof(getenv("H648K")) : 0.5;
static double g_h648_s = getenv("H648S0") ? atof(getenv("H648S0")) : 1.0, g_h648_idle = 0, g_h648_mid = 0;
static double g_h242d_fix = getenv("H242D_FIX") ? atof(getenv("H242D_FIX")) : -1.0;
static inline double h242d_f(double f) { if (g_h242d_fix >= 0) return std::min(1.0, g_h242d_fix); if (!g_h242d_bal) return std::min(1.0, std::max(0.0, f)); return std::min(0.5, std::max(0.0, f * g_bal_s)); }
static void h242_schedule(Engine& E, PipeStep& PS, Batch& B, Times& T, const std::vector<std::pair<int, int>>& pr);
static void h771_schedule(Engine& E, PipeStep& PS, Batch& B, Times& T, const std::vector<std::pair<int, int>>& pr);


static void schedule_range(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const std::vector<uint8_t>& admit = *PS.admit; const std::vector<int>& ct = PS.SB->ct;
  if (!g_h99 && !g_h242 && !g_h771) {
    for (int i = B.b0; i < B.b1; i++) {
      if (!admit[i]) continue;
      for (int t = 0; t < 2; t++) {
        const int id = ct[(size_t)i * 2 + t];
        if (id < 0 || PS.fitted[t][id]) continue;
        PS.fitted[t][id] = 1; B.need[t].push_back(id); PS.allids[t].push_back(id);
      }
    }
    T.nfit0 += (double)B.need[0].size(); T.nfit1 += (double)B.need[1].size();
    return;
  }
  static const int g_kh4 = getenv("KH4") ? atoi(getenv("KH4")) : 0;
  std::vector<std::pair<int, int>> pr;
  if (g_kh4) pr.reserve((size_t)(B.b1 - B.b0) * 2);
  for (int i = B.b0; i < B.b1; i++) {
    if (!admit[i]) continue;
    for (int t = 0; t < 2; t++) {
      const int id = ct[(size_t)i * 2 + t];
      if (id < 0 || PS.fitted[t][id]) continue;
      PS.fitted[t][id] = 1; PS.allids[t].push_back(id); pr.push_back({t, id});
    }
  }
  const size_t np = pr.size();
  if (g_h771) { h771_schedule(E, PS, B, T, pr); return; }
  if (g_h242) { h242_schedule(E, PS, B, T, pr); return; }
  std::vector<uint32_t> fe(np, H99_NONE);
  if (np) {
    const int nw = (int)WK(E).size();
    PL(E).run(nw, [&](int w) {
      for (size_t q = (size_t)w; q < np; q += (size_t)nw) {
        const Tree& tr = PS.SB->UQ.U[(size_t)pr[q].second];
        fe[q] = h99_find(PS.SB->UQ.hash[(size_t)pr[q].second], (uint32_t)pr[q].first, tr.rec.data(), (uint32_t)tr.rec.size());
      }
    });
  }
  size_t nsch[2] = {0, 0};
  for (size_t q = 0; q < np; q++) {
    const int t = pr[q].first, id = pr[q].second;
    nsch[t]++; h99_lookup_count();
    if (fe[q] != H99_NONE) { h99_hit_count(); h99_apply(PS, t, id, fe[q], T); }
    else B.need[t].push_back(id);
  }
  T.nfit0 += (double)nsch[0]; T.nfit1 += (double)nsch[1];
}


static const int g_h250b = getenv("H250B") ? atoi(getenv("H250B")) : 0;
static const int g_h334 = getenv("H334") ? atoi(getenv("H334")) : 0;
static const int g_h250bt = getenv("H250BT") ? std::max(1, atoi(getenv("H250BT"))) : 6;
static double g_h250b_rg = 0, g_h250b_rc = 0;
struct H250bTree { uint32_t j; int Q; bool s2; int ns, nc; float tau; };
struct H250bIn { const uint32_t* hdr; const uint32_t* tok; const uint16_t* sp; const float* th; const uint32_t* seg; float* out; const float* sob; const float* yb; const float* bct; };
static void h250b_cpu(H250bIn in, std::vector<H250bTree> v, int nthr, double* ms_out) {
  auto t0 = clk::now();
  std::sort(v.begin(), v.end(), [](const H250bTree& a, const H250bTree& b) {
    if (a.Q != b.Q) return a.Q < b.Q; if (a.s2 != b.s2) return a.s2 < b.s2; if (a.ns != b.ns) return a.ns < b.ns; return a.nc < b.nc; });
  s250b::Inputs I{in.hdr, in.tok, in.sp, in.th, in.sob, in.yb, in.seg, in.bct};
  auto work = [&](size_t a, size_t b) { H401S h401s_("H250B cpu");
    if (g_kqos & 1) k_qos_set(1); else h298_util(1);
    s250b::Work<float> W; int pQ = -1, ps = -1, pn = -1, pc = -1;
    for (size_t i = a; i < b; i++) {
      const H250bTree& t = v[i];
      if (t.Q != pQ || (int)t.s2 != ps || t.ns != pn || t.nc != pc) { W.prepare(t.Q, t.s2, t.ns, t.nc); pQ = t.Q; ps = (int)t.s2; pn = t.ns; pc = t.nc; }
      s250b::run_tree<float>(W, I, t.j, t.tau, in.out + 16 * (size_t)t.j);
    }
  };
  const int nt = std::max(1, std::min(nthr, (int)v.size()));
  std::vector<std::thread> th;
  if (tl_in_fitw) work(0, v.size());
  else {
  for (int k = 0; k < nt; k++) th.emplace_back(work, v.size() * k / nt, v.size() * (k + 1) / nt);
  for (auto& x : th) x.join();
  }
  *ms_out = ms_since(t0);
}


static void scan_begin(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;
  const int M = (int)UQ.U.size();
  B.M = M;
  Persist5& P = *B.P;
  id<MTLDevice> dev = E.dev;
  if (!P.bct) { const float ctv[16] = {1,0,0,0, 0,1,1,0, 0,1,-1,0, 0,-1,-1,0}; P.bct = E.B(ctv, sizeof ctv); P.bsob = E.B(SOB_U, sizeof(SOB_U)); }
  for (int s = 0; s < 2; s++) if (!P.S[s].yb) P.S[s].yb = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  const int NCH = 4 * (int)WK(E).size();
  auto span = [M, NCH](int c, int& j0, int& j1) { j0 = (int)((long)M * c / NCH); j1 = (int)((long)M * (c + 1) / NCH); };

  auto t0 = clk::now();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s];
    Tick tk;
    X.listed.assign(M, 0); X.sok.assign(M, 0); X.ninst.assign(M, 0);
    X.toff.assign(M, 0); X.spoff.assign(M, 0); X.thoff.assign(M, 0);
    for (int j : B.need[s]) X.listed[j] = 1;
    T.h_s_alloc += tk.lap();
    X.ti = X.si = X.hi = 0;
    for (int j = 0; j < M; j++) {
      X.toff[j] = X.ti; X.spoff[j] = X.si; X.thoff[j] = X.hi;
      if (!X.listed[j]) continue;
      const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
      if (te.s_ok) { X.sok[j] = 1; X.ninst[j] = (uint32_t)(te.nsp() / 4); X.ti += X.ninst[j]; } else T.scan_fail += 1;
      X.si += (uint32_t)tr.sp.size(); X.hi += (uint32_t)tr.th.size();
    }
    int Qmax = 0; for (int j = 0; j < M; j++) if (X.sok[j]) Qmax = std::max(Qmax, UQ.U[j].Q);
    T.h_s_pre += tk.lap();
    const size_t DN = (size_t)(Qmax + 1) * 41 * 8;
    X.seq.clear(); X.dkey.assign(M, 0); X.cpuJ.clear();
    for (int j = 0; j < M; j++) {
      if (!X.sok[j]) continue;
      const Tree& tr = UQ.U[j];
      const uint32_t cls = (tr.nS <= 4 ? 0u : tr.nS <= 8 ? 1u : tr.nS <= 16 ? 2u : 3u) + (tr.nth == 0 ? 4u : 0u);
      X.dkey[j] = (uint32_t)(((size_t)tr.Q * 41 + (size_t)tr.nslot) * 8 + cls);
      X.seq.push_back((uint32_t)j);
    }
    T.h_s_pre += tk.lap();
    if (g_h250b) {
      static bool once_s = true; if (once_s) { once_s = false; fprintf(stderr, "H250B active path=s\n"); }
      const double f = h242d_f((g_h250b_rg > 0 && g_h250b_rc > 0) ? g_h250b_rc / (g_h250b_rc + g_h250b_rg) : 0.5);
      const uint32_t thr = (uint32_t)(f * 1000.0 + 0.5);
      std::vector<uint32_t> keep; keep.reserve(X.seq.size());
      for (uint32_t j : X.seq) { if ((j * 2654435761u) % 1000u < thr) X.cpuJ.push_back(j); else keep.push_back(j); }
      X.seq.swap(keep);
    }
    sort_by_key5(X.seq, X.dkey, DN, X.start, X.fill, X.sorted);
    T.h_s_sort += tk.lap();
    X.sbins.clear(); size_t cur = 0;
    for (size_t dk = 0; dk < DN; dk++) {
      const size_t c = X.start[dk + 1] - X.start[dk]; if (!c) continue;
      SBin5 b; b.Q = (int)(dk / (41 * 8)); b.ns = (int)((dk / 8) % 41); const int cls = (int)(dk % 8);
      b.s2 = cls >= 4; b.nc = 4 << (cls & 3); b.off = cur; b.cnt = c; b.sst = X.start[dk];
      X.sbins.push_back(b); cur = al5(cur + c, BIN_ALIGN);
    }
    X.bl_words = cur;
    T.h_s_bins += tk.lap();
    h900_ens(X.tok, "X.tok", dev, ((size_t)X.ti * 4 + 8) * 4, MSH);
    h900_ens(X.hdr, "X.hdr", dev, (size_t)M * 16, MSH);
    h900_ens(X.sp, "X.sp", dev, (size_t)X.si * 2, MSH);
    h900_ens(X.th, "X.th", dev, (size_t)X.hi * 4, MSH);
    h900_ens(X.seg, "X.seg", dev, (size_t)M * 144, MSH);
    h900_ens(X.out, "X.out", dev, (size_t)M * 64, MSH);
    h900_ens(X.bl, "X.bl", dev, X.bl_words * 4, MSH);
    T.h_s_alloc += tk.lap();
    uint32_t* tok = (uint32_t*)X.tok.b.contents; uint32_t* hdr = (uint32_t*)X.hdr.b.contents; uint16_t* sp = (uint16_t*)X.sp.b.contents;
    float* th = (float*)X.th.b.contents; uint32_t* seg = (uint32_t*)X.seg.b.contents; float* out = (float*)X.out.b.contents;
    uint32_t* bl = (uint32_t*)X.bl.b.contents;
    T.h_s_alloc += tk.lap();
    for (const SBin5& b : X.sbins) std::memcpy(bl + b.off, X.sorted.data() + b.sst, b.cnt * 4);
    T.h_s_bins += tk.lap();
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      for (int j = j0; j < j1; j++) {
        uint32_t* h = hdr + 4 * (size_t)j; float* ob = out + 16 * (size_t)j; uint32_t* sg = seg + 36 * (size_t)j;
        h[0] = X.toff[j]; h[2] = X.spoff[j];
        for (int q = 0; q < 16; q++) ob[q] = INFINITY;
        if (!X.listed[j]) { h[1] = 0; h[3] = 0; std::memset(sg, 0, 144); continue; }
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.ninst[j] | ((uint32_t)tr.nS << 16);
        h[3] = (uint32_t)tr.nth | (X.thoff[j] << 8);
        if (te.s_ok) { std::memcpy(sg, te.sg, 144); std::memcpy(tok + 4 * (size_t)X.toff[j], te.psp(), te.nsp() * 4); }
        else std::memset(sg, 0, 144);
        if (!tr.sp.empty()) std::memcpy(sp + X.spoff[j], tr.sp.data(), tr.sp.size() * 2);
        if (!tr.th.empty()) std::memcpy(th + X.thoff[j], tr.th.data(), tr.th.size() * 4);
      }
    });
    for (int z = 0; z < 8; z++) tok[4 * (size_t)X.ti + z] = 0;
    if (g_h250b && !X.cpuJ.empty()) {
      std::vector<H250bTree> v; v.reserve(X.cpuJ.size());
      for (uint32_t j : X.cpuJ) {
        const Tree& tr = UQ.U[j];
        const uint32_t cls = (tr.nS <= 4 ? 0u : tr.nS <= 8 ? 1u : tr.nS <= 16 ? 2u : 3u) + (tr.nth == 0 ? 4u : 0u);
        v.push_back({j, (int)tr.Q, cls >= 4u, (int)tr.nslot, (int)(4u << (cls & 3u)), (float)tr.Q * 9.0f * ldexpf(1.0f, -24)});
      }
      H250bIn in{ (const uint32_t*)X.hdr.b.contents, (const uint32_t*)X.tok.b.contents, (const uint16_t*)X.sp.b.contents,
                  (const float*)X.th.b.contents, (const uint32_t*)X.seg.b.contents, (float*)X.out.b.contents,
                  (const float*)P.bsob.contents, (const float*)X.yb.contents, (const float*)P.bct.contents };
      X.cpu_ms = 0; X.cpuT = fw_submit([in, v, pms = &X.cpu_ms]() mutable { h250b_cpu(in, v, fw_persist() ? 1 : g_h250bt, pms); });
    }
    T.h_s_fill += tk.lap();
  }
  T.scan_enc += ms_since(t0);


  Tick tc;
  size_t nsb = 0; for (int s = 0; s < 2; s++) nsb += P.S[s].sbins.size();
  if (nsb) {
    h900_ens(P.scr, "P.scr", dev, (size_t)M * 16 * 16 * 4, MPV);


    id<MTLCommandBuffer> cb = h182_cb(E.cq); id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (int s = 0; s < 2; s++) {
      Side5& X = P.S[s];
      for (const SBin5& b : X.sbins) {
        Tick tp;
        id<MTLComputePipelineState> ps = E.pipe(E.Ls, b.s2 ? "rcscan_s2_q" + std::to_string(b.Q) : "rcscan_q" + std::to_string(b.Q));
        T.h_s_pipe += tp.lap();
        if (!ps) { fprintf(stderr, "missing scan kernel Q%d\n", b.Q); exit(1); }
        [en setComputePipelineState:ps];
        id<MTLBuffer> bs[9] = {X.tok.b, X.hdr.b, X.sp.b, X.th.b, P.bsob, X.yb, X.out.b, X.bl.b, X.seg.b};
        for (int k = 0; k < 9; k++) [en setBuffer:bs[k] offset:(k == 7 ? b.off * 4 : 0) atIndex:k];
        uint32_t ns = (uint32_t)b.ns, nc = (uint32_t)b.nc; float tau = b.Q * 9.0f * ldexpf(1.0f, -24);
        [en setBuffer:P.bct offset:0 atIndex:11]; [en setBuffer:P.scr.b offset:0 atIndex:12];
        [en setBytes:&tau length:4 atIndex:13]; [en setBytes:&ns length:4 atIndex:9]; [en setBytes:&nc length:4 atIndex:10];
        size_t tgb = (size_t)ns * (b.s2 ? 96 : 768) * 4 + (b.s2 ? 2 : 16) * nc * 4 + 56 * 4; tgb = (tgb + 15) & ~(size_t)15;
        uint32_t sstr = (uint32_t)(((size_t)ns * 96 + 2 * nc + 59) & ~(size_t)3); [en setBytes:&sstr length:4 atIndex:14];
        [en setThreadgroupMemoryLength:tgb atIndex:0];
        const int thr = b.s2 ? 96 / 3 : (b.Q <= 5 ? 768 / 6 : 768 / 3);
        if (h150s_on() || h150c_on()) { g_h150.scDisp++; g_h150.scTG += (long)b.cnt; g_h150.scThr += (long)thr * (long)b.cnt; g_h150.scSimd += (long)b.cnt * (long)(thr / 32); }
        [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)];
      }
    }
    [en endEncoding]; [cb commit];
    B.cbS = cb;
  }
  T.h_s_cmd += tc.lap();
}


#include <deque>
static const int g_h882 = getenv("H882") ? atoi(getenv("H882")) : 0;
struct H882Q {
  std::mutex m; std::condition_variable cv, cvd; std::deque<std::function<void()>> q;
  unsigned long long sub = 0, done = 0; bool started = false;
  double busy[2] = {0, 0}, wait = 0;
  void loop() {
    for (;;) {
      std::function<void()> f;
      { std::unique_lock<std::mutex> lk(m); cv.wait(lk, [&] { return !q.empty(); }); f = std::move(q.front()); q.pop_front(); }
      if (g_h917 & 1) { @autoreleasepool { f(); } } else f();
      { std::lock_guard<std::mutex> lk(m); ++done; }
      cvd.notify_all();
    }
  }
  unsigned long long submit(int kind, std::function<void()> body) {
    std::lock_guard<std::mutex> lk(m);
    if (!started) { started = true; (void)h182_on(); (void)h150c_on(); std::thread([this] { this->loop(); }).detach(); }
    q.push_back([this, kind, body = std::move(body)]() { const auto a = clk::now(); body(); const double d = ms_since(a); std::lock_guard<std::mutex> lk2(m); busy[kind] += d; });
    cv.notify_one();
    return ++sub;
  }
  bool is_done(unsigned long long t) { if (!t) return true; std::lock_guard<std::mutex> lk(m); return done >= t; }
  void wait_for(unsigned long long t) {
    if (!t) return;
    std::unique_lock<std::mutex> lk(m);
    if (done >= t) return;
    const auto a = clk::now();
    cvd.wait(lk, [&] { return done >= t; });
    wait += ms_since(a);
  }
};
static H882Q* g_h882q = new H882Q();
static void h882_drain() { g_h882q->wait_for(g_h882q->sub); }
static void h882_report() {
  if (!(g_h882 & 4)) return;
  std::lock_guard<std::mutex> lk(g_h882q->m);
  printf("H882: helper busy scan %.2f GN %.2f ms | main wait-for-helper %.2f ms\n", g_h882q->busy[0], g_h882q->busy[1], g_h882q->wait);
  g_h882q->busy[0] = g_h882q->busy[1] = g_h882q->wait = 0;
}


static const int g_kahead = getenv("K_AHEAD") ? atoi(getenv("K_AHEAD")) : 1;
static const int g_kpopkey = getenv("K_POPKEY") ? atoi(getenv("K_POPKEY")) : 1;
static const int g_kcap = getenv("K_CAP") ? atoi(getenv("K_CAP")) : 1;
static const int g_kverb = getenv("K_VERB") ? atoi(getenv("K_VERB")) : 0;
static std::function<bool()> g_kfill; static std::function<void()> g_kfin; 
static std::vector<uint64_t> g_kskey; static std::vector<uint8_t> g_kvalid; static std::vector<uint32_t> g_kcs_prev;
static inline void k_flush() { if (g_kfill) { while (g_kfill()) {} g_kfill = nullptr; } if (g_kfin) { auto f = std::move(g_kfin); g_kfin = nullptr; f(); } }
static std::function<bool()> g_tfill;
static inline void k_wait(id<MTLCommandBuffer> cb) {
  if (g_tfill) { while (cb.status < MTLCommandBufferStatusCompleted) { if (!g_tfill()) break; } }
  if (g_kfill) { while (cb.status < MTLCommandBufferStatusCompleted) { if (!g_kfill()) break; } }
  [cb waitUntilCompleted];
}
static void scan_finish(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  auto t0 = clk::now();
  g_h882q->wait_for(B.t882S);
  if (B.cbS) {
    k_wait(B.cbS);
    if (B.cbS.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "scan command buffer failed\n"); h182_dump(B.cbS, "scan"); gpu_fail(); }
    B.scan_gpu = (B.cbS.GPUEndTime - B.cbS.GPUStartTime) * 1e3;
    PS.iv.push_back({B.cbS.GPUStartTime, B.cbS.GPUEndTime});
  }
  if (g_h250b) {
    size_t ng = 0, nc = 0; double cms = 0;
    for (int s = 0; s < 2; s++) {
      Side5& X = B.P->S[s];
      if (X.cpuT.joinable()) X.cpuT.join();
      ng += X.seq.size(); nc += X.cpuJ.size(); cms = std::max(cms, X.cpu_ms);
    }
    if (B.scan_gpu > 0 && ng > 0) { const double r = (double)ng / B.scan_gpu; g_h250b_rg = g_h250b_rg > 0 ? 0.5 * g_h250b_rg + 0.5 * r : r; }
    if (cms > 0 && nc > 0) { const double r = (double)nc / cms; g_h250b_rc = g_h250b_rc > 0 ? 0.5 * g_h250b_rc + 0.5 * r : r; }
    printf("H250B: f_cpu %.3f cpu trees %zu cpu ms %.3f gpu trees %zu gpu ms %.3f\n", nc + ng ? (double)nc / (double)(nc + ng) : 0.0, nc, cms, ng, B.scan_gpu);
  }
  { const double w = ms_since(t0); T.w_s += w; T.scan_wall += w; }
  T.scan_gpu += B.scan_gpu;
  const Unique& UQ = PS.SB->UQ; Persist5& P = *B.P;
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s]; const Target& tg = E.tg[s];
    const float* SO = (const float*)X.out.b.contents;
    for (int j : B.need[s]) {
      if (!X.sok[j]) continue;
      const float* r = SO + lid916(X, (size_t)j) * 16;
      if (UQ.U[j].nth == 0) { o.smin[j] = r[0]; o.sarg[j] = 0; }
      else for (int q = 0; q < 16; q++) if (r[q] < o.smin[j]) { o.smin[j] = r[q]; o.sarg[j] = q; }
      T.s_alive += isfinite(o.smin[j]) ? 1 : 0;
      if (UQ.U[j].nth >= 1 && isfinite(o.smin[j])) { T.gn_cand += 1; if (o.smin[j] < tg.floor_) T.gn_floor += 1; }
    }
    for (int j : B.need[s]) if (UQ.U[j].nth == 0) o.rfin[j] = o.smin[j];
  }
  if (g_h139 == 1) for (int s = 0; s < 2; s++) if (P.S[s].h139n) h139_cmp(P.S[s], s, UQ);
  if (g_h335 == 2) {
    double n = 0, mis = 0, kept = 0, lr = 0;
    for (int s = 0; s < 2; s++) { const Side5& X = P.S[s]; const float* SO = (const float*)X.out.b.contents; const float* PO = (const float*)X.h139o.b.contents;
      const uint32_t* BL = (const uint32_t*)X.bl.b.contents;
      for (const auto& bb : X.h335b) for (size_t k = 0; k < bb.second; k++) { const uint32_t j = BL[bb.first + k]; if (j == 0xFFFFFFFFu) continue;
        const float* c = SO + (size_t)j * 16; const float* p = PO + (size_t)j * 16; float cm = INFINITY, pm = INFINITY; int pa = -1;
        for (int q = 0; q < 16; q++) { if (isfinite(c[q]) && memcmp(&c[q], &p[q], 4) != 0) mis += 1; if (c[q] < cm) cm = c[q]; if (p[q] < pm) { pm = p[q]; pa = q; } }
        n += 1; if (pa >= 0 && c[pa] == p[pa] && isfinite(c[pa])) kept += 1; if (isfinite(cm) && isfinite(pm) && pm > 0) lr += log((double)cm / (double)pm); } }
    g_h335_n += n; g_h335_mis += mis; g_h335_kept += kept; g_h335_lr += lr;
    printf("H335 CHK: trees %.0f | chosen rows != production %.0f | production winner kept %.4f | geo cut/prod min %.6f\n", n, mis, n ? kept / n : 0.0, n ? exp(lr / n) : 1.0);
  }
  B.sdone = true;
}


static int h94_mode() { static const int m = getenv("H94") ? atoi(getenv("H94")) : 0; return m; }
static void h94_cap(std::vector<size_t>& ntg, const std::vector<size_t>& nat, const std::vector<size_t>& tgb, const std::vector<size_t>& cnt) {
  static const double cores = getenv("H94_CORES") ? atof(getenv("H94_CORES")) : 7.0;
  double N = 0; for (size_t c : cnt) N += (double)c;
  for (size_t b = 0; b < ntg.size(); b++) {
    const double rpc = (double)std::max<size_t>(1, 65536 / std::max<size_t>(1, tgb[b]));
    const double cap = std::max(1.0, std::floor(cores * rpc * (double)cnt[b] / std::max(1.0, N)));
    ntg[b] = std::min<size_t>(nat[b], (size_t)cap);
  }
}


static int h141a_on() { static const int v = getenv("H141A") ? atoi(getenv("H141A")) : 1; return v; }
static int h141b_on() { static const int v = getenv("H141B") ? atoi(getenv("H141B")) : 1; return v; }
static std::atomic<long> g_h141_badL{0}, g_h141_badO{0};
struct H141Acc { long bins = 0, changed = 0, shrunkTrees = 0; double trees = 0, simB = 0, simA = 0, wB = 0, wA = 0; };
static H141Acc g_h141;
static void h141_step_reset() { g_h141 = H141Acc(); g_h141_badL = 0; g_h141_badO = 0; }
static uint32_t h141_pe(const TreeEnc& te, uint32_t gper) {
  if (!h141a_on()) return gper;
  const uint32_t nsc = te.ghw & 255u, SB = (te.ghw >> 8) & 0xFFFFu, nth = te.ghw >> 24;
  const uint32_t CD = nth <= 4 ? 1u + (nth + 1u) / 2u : 3u;
  if (gper != SB + 2u * CD * nsc) { g_h141_badL++; return gper; }
  const uint32_t full = SB + CD * nsc;
  auto ok = [&](uint32_t w) -> bool {
    if ((w >> 31) || !((w >> 30) & 1u)) return true;
    return (w & 0x1FFFFFFFu) + (CD - 1u) * nsc < full;
  };
  for (size_t k = 0; k + 3 < te.ngp(); k += 4) if (!ok(te.pgp()[k + 2]) || !ok(te.pgp()[k + 3])) { g_h141_badO++; return gper; }
  for (int q = 18; q < 36; q++) if (!ok(te.gg[q])) { g_h141_badO++; return gper; }
  return full;
}

static int g_h332 = getenv("H332") ? atoi(getenv("H332")) : 0;
static uint32_t h332_pe(const TreeEnc& te, uint32_t pe) {
  if (!g_h332) return pe;
  const uint32_t nsc = te.ghw & 255u, SB = (te.ghw >> 8) & 0xFFFFu, nth = te.ghw >> 24;
  const uint32_t CD = nth <= 4 ? 1u + (nth + 1u) / 2u : 3u, Ns = SB / (CD * 48u), SBf = Ns * (48u + (CD - 1u) * 24u);
  return pe == SB + CD * nsc ? SBf + CD * nsc : SBf + 2u * CD * nsc;
}

static int h141_fit_tpg(size_t shb, size_t pers, int cap) {
  int t = 8; while (t > 1 && (shb + pers * (size_t)t) * 4 > 32768) t--;
  return std::min<int>(t, cap);
}

static long h141_sim(size_t shb, size_t pers, int t) {
  const size_t g = ((shb + pers * (size_t)t) * 4 + 15) & ~(size_t)15;
  return std::min<long>(96, (long)(65536 / std::max<size_t>(g, 1)) * (long)t);
}
static int h141_pick(size_t shb, size_t persN, size_t cnt, int cap) {
  const int t0 = h141_fit_tpg(shb, persN, cap);
  if (!h141b_on()) return t0;
  const int tmax = std::max(1, std::min<int>(cap, (int)std::min<size_t>(cnt, 1024)));
  int best = t0; long bs = -1; int bd = 1 << 30;
  for (int t = 1; t <= tmax; t++) {
    if ((shb + persN * (size_t)t) * 4 > 32768) break;
    const long s = h141_sim(shb, persN, t); const int d = std::abs(t - t0);
    if (s > bs || (s == bs && (d < bd || (d == bd && t > best)))) { best = t; bs = s; bd = d; }
  }
  return best;
}
static void h141_note(size_t shb, size_t pers, size_t persN, size_t cnt, int cap, int tpg) {
  if (!h141a_on() && !h141b_on()) return;
  const int tU = h141_fit_tpg(shb, pers, cap);
  g_h141.bins++; g_h141.changed += (tpg != tU);
  const double sB = (double)h141_sim(shb, pers, tU), sA = (double)h141_sim(shb, persN, tpg);
  g_h141.simB += sB; g_h141.simA += sA; g_h141.wB += sB * (double)cnt; g_h141.wA += sA * (double)cnt; g_h141.trees += (double)cnt;
  if (persN < pers) g_h141.shrunkTrees += (long)cnt;
}

static size_t h713_xf(int Tk) { return 4 + 2 * (size_t)Tk * 48; }
static size_t h713_gn(Engine& E, int Qk, int Tk, size_t shb, size_t persN, __strong id<MTLComputePipelineState>& ps, int& tpg) {
  if (!g_h713 || Tk < 9) return 0;
  id<MTLComputePipelineState> pp = E.pipe(E.Lgp, "gnp_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
  if (!pp || (int)pp.maxTotalThreadsPerThreadgroup < 64) return 0;
  const size_t xf = h713_xf(Tk);
  if ((shb + 2 * persN + xf) * 4 > 32768) return 0;
  ps = pp; tpg = 2; return xf; }
static void h141_step_print(int st) {
  if (!h141a_on() && !h141b_on()) return;
  const double n = (double)std::max<long>(g_h141.bins, 1), tr = std::max(g_h141.trees, 1.0);
  printf("H141 step %2d | A %d B %d | GN bins %ld | tpg changed %ld | resident SIMD/core mean per bin before %.2f after %.2f | tree-weighted %.2f -> %.2f | trees in smaller-TG bins %ld | layout fallback %ld | operand fallback %ld\n",
         st, h141a_on(), h141b_on(), g_h141.bins, g_h141.changed, g_h141.simB / n, g_h141.simA / n, g_h141.wB / tr, g_h141.wA / tr,
         g_h141.shrunkTrees, g_h141_badL.load(), g_h141_badO.load());
}


struct H874GNPool { std::mutex m; std::vector<GNWork*> v;
  GNWork* get() { { std::lock_guard<std::mutex> l(m); if (!v.empty()) { GNWork* p = v.back(); v.pop_back(); p->k644_[0] = GNWork::K644(); p->k644_[1] = GNWork::K644(); return p; } } return new GNWork(); }
  void put(GNWork* p) { std::lock_guard<std::mutex> l(m); v.push_back(p); } };
static H874GNPool g_h874_gn;
static void h250a_body(Batch* B, Persist5* P, double* ms_out) {
  const auto t0 = clk::now();
  const size_t n = B->h250_items.size();
  std::atomic<size_t> next{0};
  std::vector<std::thread> th;
  auto h250a_wb = [&]() { H401S h401s_("GN cpu (H250A)");
    if (g_kqos & 1) k_qos_set(1); else h298_util(1);
    GNWork* W874 = (g_h874 & 4) ? g_h874_gn.get() : nullptr; GNWork Wl; GNWork& W = W874 ? *W874 : Wl;
    GNResult r; float tl[48]; gn_tlin(tl);
    const bool pf836 = (g836_bits() & 8) != 0;
    auto mkin836 = [&](const H250aItem& jt) { Side5& Y = P->S[jt.side];
      return GNInputs{(const uint32_t*)Y.gtok.b.contents, (const uint32_t*)Y.ghdr.b.contents, (const uint16_t*)Y.gsp.b.contents,
                      (const uint32_t*)Y.gseg.b.contents, (const float*)[Y.yb contents], tl, (const float*)Y.gth.b.contents, jt.per}; };
    auto pf1grp836 = [&](size_t g0) { for (size_t k = g0; k < std::min(n, g0 + 4); k++) { const H250aItem& jt = B->h250_items[k]; g836_pf1(mkin836(jt), jt.tree); } };
    size_t cur836 = next.fetch_add(4), nxt836 = (size_t)-1;
    if (pf836 && cur836 < n) { pf1grp836(cur836); nxt836 = next.fetch_add(4); if (nxt836 < n) pf1grp836(nxt836); }
    for (;;) {
      const size_t i0 = cur836;
      if (i0 >= n) break;
      h891_cpu_batch(1);
      for (size_t i = i0; i < std::min(n, i0 + 4); i++) {
        if (pf836) { const size_t jn = (i + 1 < std::min(n, i0 + 4)) ? i + 1 : nxt836; if (jn < n) { const H250aItem& jt = B->h250_items[jn]; g836_pf2(mkin836(jt), jt.tree, jt.Q); } }
        if (h650_is("gncpu")) { h650_spin_ns(g_h650ns); g_h650_items.fetch_add(1, std::memory_order_relaxed); }
        const H250aItem& it = B->h250_items[i];
        Side5& X = P->S[it.side];
        GNInputs in{(const uint32_t*)X.gtok.b.contents, (const uint32_t*)X.ghdr.b.contents, (const uint16_t*)X.gsp.b.contents,
                    (const uint32_t*)X.gseg.b.contents, (const float*)[X.yb contents], tl, (const float*)X.gth.b.contents, it.per};
        W.polish(in, it.tree, it.Q, (float)it.Q * 9.0f * ldexpf(1.0f, -24), it.floor_, 3, r);
        float* thp = (float*)X.gth.b.contents + 32 * (size_t)it.tree;
        for (int j = 0; j < it.nth; j++) thp[j] = r.th[j];
        float* ap = (float*)X.gai.b.contents + 16 * (size_t)it.tree;
        for (int q = 0; q < 16; q++) ap[q] = r.a[q];
        float* rsp = (float*)X.grs.b.contents + 2 * (size_t)it.tree;
        rsp[0] = r.rss0; rsp[1] = r.rcur;
        ((uint32_t*)X.gst.b.contents)[it.tree] = r.steps | (r.npp << 8) | (r.ndp << 16);
      }
      if (pf836) { cur836 = nxt836; nxt836 = (size_t)-1; if (cur836 < n) { nxt836 = next.fetch_add(4); if (nxt836 < n) pf1grp836(nxt836); } }
      else cur836 = next.fetch_add(4);
    }
    if (W874) g_h874_gn.put(W874);
  };
  if (tl_in_fitw) h250a_wb();
  else if (g_h873 & 1) g_h873t.run(0, g_h250at, h250a_wb);
  else { for (int t = 0; t < g_h250at; t++) th.emplace_back(h250a_wb); for (auto& x : th) x.join(); }
  *ms_out = ms_since(t0);
}
static Batch* g_h250_pending = nullptr;
static const int g_h722 = getenv("H722") ? atoi(getenv("H722")) : 0;
static void h250a_join(Batch& B) {
  if (B.h250_thr.joinable()) B.h250_thr.join();
  if (g_h250_pending == &B) g_h250_pending = nullptr;
}
static void gn_begin(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;


  const int M = B.M;
  if ((!g_h217 && (int)PS.SB->UQ.U.size() < M) || (int)B.P->S[0].listed.size() != M || (int)B.P->S[1].listed.size() != M) { fprintf(stderr, "pipeline: slot arrays do not match chunk %d\n", B.b0); exit(1); }
  Persist5& P = *B.P; id<MTLDevice> dev = E.dev;
  const int NCH = 4 * (int)WK(E).size();
  auto span = [M, NCH](int c, int& j0, int& j1) { j0 = (int)((long)M * c / NCH); j1 = (int)((long)M * (c + 1) / NCH); };
  auto t0 = clk::now();
  Tick tk;
  P.gbins.clear();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    X.gok.assign(M, 0); X.gninst.assign(M, 0); X.gtoff.assign(M, 0); X.gspoff.assign(M, 0); X.gper.assign(M, 0);
    X.gti = X.gsi = 0;
    T.h_g_alloc += tk.lap();
    for (int j = 0; j < M; j++) {
      X.gtoff[j] = X.gti; X.gspoff[j] = X.gsi;
      if (!X.listed[j]) continue;
      const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
      if (tr.nth >= 1 && tr.nth > GN_TH_MAX) T.gn_kmiss++;
      if (te.g_ok && isfinite(o.smin[j])) {
        X.gok[j] = 1; o.is_g[j] = 1; T.g_alive++;
        X.gninst[j] = (uint32_t)(te.ngp() / 4); X.gper[j] = te.gper; X.gti += X.gninst[j];
      }
      X.gsi += (uint32_t)tr.sp.size();
    }
    int Qg = 0, Cg = 0, LM = 0;
    X.seq.clear();
    for (int j = 0; j < M; j++) if (X.gok[j]) { X.seq.push_back((uint32_t)j); Qg = std::max(Qg, UQ.U[j].Q); Cg = std::max(Cg, (int)((X.gper[j] + 63) / 64)); }
    if (Qg >= 100) { fprintf(stderr, "GN: Q %d >= 100 breaks the bin key\n", Qg); exit(1); }
    for (uint32_t j : X.seq) LM = std::max(LM, (int)h313a::nprog(UQ.U[j]));
    X.keyl.assign(M, 0);
    for (uint32_t j : X.seq) X.keyl[j] = (uint32_t)(LM - (int)h313a::nprog(UQ.U[j]));
    T.h_g_pre += tk.lap();
    sort_by_key5(X.seq, X.keyl, (size_t)LM + 1, X.start, X.fill, X.seq2);
    T.h_g_sort += tk.lap();
    X.dkey.assign(M, 0);
    for (uint32_t j : X.seq2) { const Tree& tr = UQ.U[j]; const size_t cls = (X.gper[j] + 63) / 64;
      X.dkey[j] = (uint32_t)((cls * (size_t)(Qg + 1) + (size_t)tr.Q) * 25 + (size_t)tr.nth); }
    const size_t DNg = (size_t)(Cg + 1) * (size_t)(Qg + 1) * 25;
    T.h_g_plan += tk.lap();
    sort_by_key5(X.seq2, X.dkey, DNg, X.start, X.fill, X.sorted);
    T.h_g_sort += tk.lap();
    for (size_t dk = 0; dk < DNg; dk++) {
      const size_t c = X.start[dk + 1] - X.start[dk]; if (!c) continue;
      const size_t sst = X.start[dk];
      const int Tk = (int)(dk % 25); const size_t rem = dk / 25; const int Qk = (int)(rem % (size_t)(Qg + 1));
      Tick tp;
      id<MTLComputePipelineState> ps = E.pipe(E.Lg, "gnw_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
      T.h_g_pipe += tp.lap();
      if (!ps) { T.gn_kmiss += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      size_t per = 0; for (size_t k = 0; k < c; k++) per = std::max<size_t>(per, X.gper[X.sorted[sst + k]]);
      size_t pen = 0, pen32 = 0; for (size_t k = 0; k < c; k++) { const size_t p32_ = h141_pe(UQ.enc[X.sorted[sst + k]], X.gper[X.sorted[sst + k]]); pen32 = std::max<size_t>(pen32, p32_); pen = std::max<size_t>(pen, h332_pe(UQ.enc[X.sorted[sst + k]], (uint32_t)p32_)); }
      const int NPk = std::max(4, (Tk + 3) / 4 * 4), QPk = std::max(4, (Qk + 3) / 4 * 4); const size_t WSk = (4 * NPk + 3 * QPk + 6) / 4 * 4;
      const size_t pers = WSk + 2 * per; const size_t shb = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
      if (pers * 4 > 32768) { T.gn_drop += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      const size_t persN = WSk + 2 * pen;
      const int cap141 = (int)ps.maxTotalThreadsPerThreadgroup / 32;
      int tpg = h141_pick(shb, persN, c, cap141);
      h141_note(shb, pers, persN, c, cap141, tpg);
      const size_t xf713 = h713_gn(E, Qk, Tk, shb, persN, ps, tpg);
      GBin5 g{}; g.side = s; g.Q = Qk; g.nth = Tk; g.tpg = tpg; g.cnt = c;
      g.nt = (c + (size_t)tpg - 1) / (size_t)tpg * (size_t)tpg;
      g.ntg = xf713 ? g.nt : g.nt / (size_t)tpg;
      g.tgb = ((shb + persN * tpg + xf713) * 4 + 15) & ~(size_t)15; g.per = pen; g.per32 = pen32; g.sst = sst; g.ps = ps;
      g.eb_bytes = g.ntg * (size_t)tpg * 2 * (size_t)(Qk + 1) * 64 * 4;
      P.gbins.push_back(g);
    }
    T.h_g_plan += tk.lap();
  }
  if (h94_mode() >= 1) {
    std::vector<size_t> nat, tgb, cnt, cap; size_t before = 0, after = 0, ntr = 0;
    for (const GBin5& g : P.gbins) { nat.push_back(g.ntg); tgb.push_back(g.tgb); cnt.push_back(g.cnt); before += g.ntg; ntr += g.cnt; }
    cap.assign(nat.size(), 0); h94_cap(cap, nat, tgb, cnt);
    for (size_t b = 0; b < P.gbins.size(); b++) { GBin5& g = P.gbins[b]; g.ntg = cap[b]; g.eb_bytes = g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4; after += g.ntg; }
    fprintf(stderr, "H94 GN dispatch: bins %zu trees %zu TGs before %zu after %zu\n", P.gbins.size(), ntr, before, after);
  }
  size_t tlc[2] = {0, 0}, ebc = 0; uint32_t NB = 0;
  for (GBin5& g : P.gbins) {
    g.wix = NB++;
    g.tl_off = tlc[g.side]; tlc[g.side] = al5(tlc[g.side] + 1 + g.nt, BIN_ALIGN);
    g.eb_off = ebc; ebc = al5(ebc + g.eb_bytes + ebt_bytes(g.Q, g.nth, g.ntg, g.tpg), 256);
  }
  T.h_g_plan += tk.lap();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    h900_ens(X.gtok, "X.gtok", dev, ((size_t)X.gti * 4 + 8) * 4, MSH);
    h900_ens(X.ghdr, "X.ghdr", dev, (size_t)M * 16, MSH);
    h900_ens(X.gsp, "X.gsp", dev, (size_t)X.gsi * 2, MSH);
    h900_ens(X.gseg, "X.gseg", dev, (size_t)M * 144, MSH);
    h900_ens(X.gth, "X.gth", dev, (size_t)M * 128, MSH);
    h900_ens(X.grs, "X.grs", dev, (size_t)M * 8, MSH);
    h900_ens(X.gst, "X.gst", dev, (size_t)M * 4, MSH);
    h900_ens(X.gai, "X.gai", dev, (size_t)M * 64, MSH);
    h900_ens(X.gbl, "X.gbl", dev, tlc[s] * 4, MSH);
    T.h_g_alloc += tk.lap();
    uint32_t* gtok = (uint32_t*)X.gtok.b.contents; uint32_t* ghdr = (uint32_t*)X.ghdr.b.contents; uint16_t* gsp = (uint16_t*)X.gsp.b.contents;
    uint32_t* gseg = (uint32_t*)X.gseg.b.contents; float* gth = (float*)X.gth.b.contents;
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      for (int j = j0; j < j1; j++) {
        uint32_t* h = ghdr + 4 * (size_t)j; float* th = gth + 32 * (size_t)j; uint32_t* sg = gseg + 36 * (size_t)j;
        for (int q = 0; q < 32; q++) th[q] = 0.0f;
        h[0] = X.gtoff[j]; h[2] = X.gspoff[j];
        if (!X.listed[j]) { h[1] = 0; h[3] = 0; std::memset(sg, 0, 144); continue; }
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.gninst[j] | ((uint32_t)tr.nS << 16);
        if (X.gok[j]) {
          h[3] = te.ghw; std::memcpy(sg, te.gg, 144);
          std::memcpy(gtok + 4 * (size_t)X.gtoff[j], te.pgp(), te.ngp() * 4);
          double th2[32]; start_theta(tr.th, tr.nth, PS.SB->F[s].sarg[j], th2);
          for (int q = 0; q < tr.nth; q++) th[q] = (float)th2[q];
        } else { h[3] = 0; std::memset(sg, 0, 144); }
        if (!tr.sp.empty()) std::memcpy(gsp + X.gspoff[j], tr.sp.data(), tr.sp.size() * 2);
      }
    });
    for (int z = 0; z < 8; z++) gtok[4 * (size_t)X.gti + z] = 0;
    T.h_g_fill += tk.lap();
    uint32_t* gbl = (uint32_t*)X.gbl.b.contents;
    for (const GBin5& g : P.gbins) {
      if (g.side != s) continue;
      uint32_t* tl = gbl + g.tl_off; tl[0] = (uint32_t)g.per;
      for (size_t k = 0; k < g.cnt; k++) tl[1 + k] = X.sorted[g.sst + k];
      for (size_t k = g.cnt; k < g.nt; k++) tl[1 + k] = 0xFFFFFFFFu;
    }
    T.h_g_fill += tk.lap();
  }
  h900_ens(P.eb, "P.eb", dev, ebc, MPV);
  h900_ens(P.ctr, "P.ctr", dev, 4 * (size_t)NB + 16, MSH);
  T.h_g_alloc += tk.lap();
  T.gn_enc += ms_since(t0);
  Tick tc;


  if (NB) {
    std::memset(P.ctr.b.contents, 0, P.ctr.cap);
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (size_t gj779 = 0; gj779 < P.gbins.size(); gj779++) { const GBin5& g = P.gbins[g_h779 ? P.gbins.size() - 1 - gj779 : gj779];
      Side5& X = P.S[g.side];
      if (h182_on()) [en pushDebugGroup:h182_gname(g)];
      [en setComputePipelineState:g.ps];
      id<MTLBuffer> bs[9] = {X.gtok.b, X.ghdr.b, X.gsp.b, X.gseg.b, X.yb, X.gth.b, X.grs.b, X.gst.b, X.gbl.b};
      for (int k = 0; k < 9; k++) [en setBuffer:bs[k] offset:(k == 8 ? g.tl_off * 4 : 0) atIndex:k];
      float tf[4] = {g.Q * 9.0f * ldexpf(1.0f, -24), E.tg[g.side].floor_, (g_freet ? 0.0f : 3.0f), 0.0f};
      [en setBytes:tf length:16 atIndex:9];
      [en setBuffer:X.gai.b offset:0 atIndex:10];
      [en setBuffer:P.eb.b offset:g.eb_off atIndex:11];
      if (g_h711) [en setBuffer:P.eb.b offset:g.eb_off + g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4 atIndex:15];
      uint32_t wq[2] = {g.wix, (uint32_t)g.nt}; [en setBytes:wq length:8 atIndex:12];
      [en setBuffer:P.ctr.b offset:0 atIndex:13];
      [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      if (h150c_on()) { g_h150.gnDisp++; g_h150.gnLaunch += (long)(g.ntg * (size_t)g.tpg); g_h150.gnOcc += (long)g.cnt; }
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
      if (h182_on()) [en popDebugGroup];
    }
    [en endEncoding]; [cb commit];
    B.cbG = cb;
  }
  T.h_g_cmd += tc.lap();
}


static int g_kh1 = getenv("KH1") ? atoi(getenv("KH1")) : 0;
static void scan_begin_p(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;
  const int M = (int)UQ.U.size();
  B.M = M;
  Persist5& P = *B.P;
  id<MTLDevice> dev = E.dev;
  auto t0 = clk::now();
  Tick tk;
  if (!P.bct) { const float ctv[16] = {1,0,0,0, 0,1,1,0, 0,1,-1,0, 0,-1,-1,0}; P.bct = E.B(ctv, sizeof ctv); P.bsob = E.B(SOB_U, sizeof(SOB_U)); }
  for (int s = 0; s < 2; s++) if (!P.S[s].yb) P.S[s].yb = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  const int NCH = 4 * (int)WK(E).size();
  auto span = [M, NCH](int c, int& j0, int& j1) { j0 = (int)((long)M * c / NCH); j1 = (int)((long)M * (c + 1) / NCH); };
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s];
    const bool l215 = (g_h215 & 1) != 0;
    if (l215) {
      if ((int)X.listed.size() < M) X.listed.resize(M, 0);
      for (int j : X.ls) if (j < (int)X.listed.size()) X.listed[j] = 0;
      X.ls.assign(B.need[s].begin(), B.need[s].end()); std::sort(X.ls.begin(), X.ls.end());
      if (!g_kh1) X.ls.erase(std::unique(X.ls.begin(), X.ls.end()), X.ls.end());
      for (int j : X.ls) X.listed[j] = 1;
      X.cpuJ.clear();
    } else {
      X.listed.assign(M, 0); X.cpuJ.clear();
      for (int j : B.need[s]) X.listed[j] = 1;
    }
    const int N215 = l215 ? (int)X.ls.size() : M;
    if (g_h916) { X.nl = (uint32_t)X.ls.size(); X.lid.resize((size_t)M); for (size_t k = 0; k < X.ls.size(); k++) X.lid[(size_t)X.ls[k]] = (uint32_t)k; }
    const size_t TBN = g_h916 ? (size_t)X.nl : (size_t)M;
    auto span = [N215, NCH](int c, int& j0, int& j1) { j0 = (int)((long)N215 * c / NCH); j1 = (int)((long)N215 * (c + 1) / NCH); };
    X.sok.resize(M); X.ninst.resize(M); X.toff.resize(M); X.spoff.resize(M); X.thoff.resize(M); X.dkey.resize(M);
    T.h_s_alloc += tk.lap();

    H263ScanScr ownS; H263ScanScr& sc263 = (g_h263 & 1) ? h263_keep_scan : ownS;
    std::vector<uint32_t>& bti = sc263.bti; std::vector<uint32_t>& bsi = sc263.bsi; std::vector<uint32_t>& bhi = sc263.bhi; std::vector<uint32_t>& bns = sc263.bns;
    std::vector<int>& bq = sc263.bq; std::vector<double>& bfl = sc263.bfl;
    bti.assign(NCH, 0); bsi.assign(NCH, 0); bhi.assign(NCH, 0); bns.assign(NCH, 0);
    bq.assign(NCH, 0); bfl.assign(NCH, 0.0);
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      uint32_t ti = 0, si = 0, hi = 0, ns = 0; int q = 0; double fl = 0;
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        if (!X.listed[j]) continue;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (te.s_ok) { ti += (uint32_t)(te.nsp() / 4); ns += 1; q = std::max(q, UQ.U[j].Q); } else fl += 1;
        si += (uint32_t)tr.sp.size(); hi += (uint32_t)tr.th.size();
      }
      bti[c] = ti; bsi[c] = si; bhi[c] = hi; bns[c] = ns; bq[c] = q; bfl[c] = fl;
    });
    std::vector<uint32_t>& oti = sc263.oti; std::vector<uint32_t>& osi = sc263.osi; std::vector<uint32_t>& ohi = sc263.ohi; std::vector<uint32_t>& osq = sc263.osq;
    oti.resize(NCH); osi.resize(NCH); ohi.resize(NCH); osq.resize(NCH);
    uint32_t TI = 0, SI = 0, HI = 0, NS = 0; int QM = 0;
    for (int c = 0; c < NCH; c++) {
      oti[c] = TI; osi[c] = SI; ohi[c] = HI; osq[c] = NS;
      TI += bti[c]; SI += bsi[c]; HI += bhi[c]; NS += bns[c]; QM = std::max(QM, bq[c]); T.scan_fail += bfl[c];
    }
    X.ti = TI; X.si = SI; X.hi = HI;
    X.seq.resize(NS);

    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      uint32_t ti = oti[c], si = osi[c], hi = ohi[c], sq = osq[c];
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        X.toff[j] = ti; X.spoff[j] = si; X.thoff[j] = hi;
        X.sok[j] = 0; X.ninst[j] = 0; X.dkey[j] = 0;
        if (!X.listed[j]) continue;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (te.s_ok) {
          X.sok[j] = 1; X.ninst[j] = (uint32_t)(te.nsp() / 4); ti += X.ninst[j];
          const uint32_t cls = (tr.nS <= 4 ? 0u : tr.nS <= 8 ? 1u : tr.nS <= 16 ? 2u : 3u) + (tr.nth == 0 ? 4u : 0u);
          X.dkey[j] = (uint32_t)(((size_t)tr.Q * 41 + (size_t)tr.nslot) * 8 + cls);
          X.seq[sq++] = (uint32_t)j;
        }
        si += (uint32_t)tr.sp.size(); hi += (uint32_t)tr.th.size();
      }
    });
    const size_t DN = (size_t)(QM + 1) * 41 * 8;
    T.h_s_pre += tk.lap();
    if (g_h250b) {
      static bool once_p = true; if (once_p) { once_p = false; fprintf(stderr, "H250B active path=p\n"); }
      const double f = h242d_f((g_h250b_rg > 0 && g_h250b_rc > 0) ? g_h250b_rc / (g_h250b_rc + g_h250b_rg) : 0.5);
      const uint32_t thr = (uint32_t)(f * 1000.0 + 0.5);
      std::vector<uint32_t> keep; keep.reserve(X.seq.size());
      for (uint32_t j : X.seq) { if ((j * 2654435761u) % 1000u < thr) X.cpuJ.push_back(j); else keep.push_back(j); }
      X.seq.swap(keep);
    }
    int nbf263 = 0;
    if (g_h263 & 2) nbf263 = sort_count_h263(E, X.seq, X.dkey, DN, X.start);
    else sort_by_key5_p(E, X.seq, X.dkey, DN, X.start, X.sorted);
    T.h_s_sort += tk.lap();
    X.sbins.clear(); size_t cur = 0;
    if (g_h263 & 2) sc263.shv.assign(DN, 0);
    for (size_t dk = 0; dk < DN; dk++) {
      const size_t c = X.start[dk + 1] - X.start[dk]; if (!c) continue;
      SBin5 b; b.Q = (int)(dk / (41 * 8)); b.ns = (int)((dk / 8) % 41); const int cls = (int)(dk % 8);
      b.s2 = cls >= 4; b.nc = 4 << (cls & 3); b.off = cur; b.cnt = c; b.sst = X.start[dk];
      X.sbins.push_back(b); if (g_h263 & 2) sc263.shv[dk] = (int64_t)cur - (int64_t)X.start[dk]; cur = al5(cur + c, BIN_ALIGN);
    }
    X.bl_words = cur;
    T.h_s_bins += tk.lap();
    h900_ens(X.tok, "X.tok", dev, ((size_t)X.ti * 4 + 8) * 4, MSH);
    h900_ens(X.hdr, "X.hdr", dev, TBN * 16, MSH);
    h900_ens(X.sp, "X.sp", dev, (size_t)X.si * 2, MSH);
    h900_ens(X.th, "X.th", dev, (size_t)X.hi * 4, MSH);
    h900_ens(X.seg, "X.seg", dev, TBN * 144, MSH);
    h900_ens(X.out, "X.out", dev, TBN * 64, MSH);
    h900_ens(X.bl, "X.bl", dev, X.bl_words * 4, MSH);
    T.h_s_alloc += tk.lap();
    uint32_t* tok = (uint32_t*)X.tok.b.contents; uint32_t* hdr = (uint32_t*)X.hdr.b.contents; uint16_t* sp = (uint16_t*)X.sp.b.contents;
    float* th = (float*)X.th.b.contents; uint32_t* seg = (uint32_t*)X.seg.b.contents; float* out = (float*)X.out.b.contents;
    uint32_t* bl = (uint32_t*)X.bl.b.contents;
    if (g_h263 & 2) sort_scatter_h263(E, X.seq, X.dkey, nbf263, DN, sc263.shv, bl);
    else PL(E).run((int)X.sbins.size(), [&](int bi) {
      const SBin5& b = X.sbins[bi]; std::memcpy(bl + b.off, X.sorted.data() + b.sst, b.cnt * 4);
    });
    if (g_h916) PL(E).run((int)X.sbins.size(), [&](int bi) {
      const SBin5& b = X.sbins[bi]; for (size_t k = 0; k < b.cnt; k++) bl[b.off + k] = X.lid[bl[b.off + k]];
    });
    T.h_s_bins += tk.lap();
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        const size_t lj = lid916(X, (size_t)j); uint32_t* h = hdr + 4 * lj; float* ob = out + 16 * lj; uint32_t* sg = seg + 36 * lj;
        h[0] = X.toff[j]; h[2] = X.spoff[j];
        for (int q = 0; q < 16; q++) ob[q] = INFINITY;
        if (!X.listed[j]) { h[1] = 0; h[3] = 0; std::memset(sg, 0, 144); continue; }
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.ninst[j] | ((uint32_t)tr.nS << 16);
        h[3] = (uint32_t)tr.nth | (X.thoff[j] << 8);
        if (te.s_ok) { std::memcpy(sg, te.sg, 144); std::memcpy(tok + 4 * (size_t)X.toff[j], te.psp(), te.nsp() * 4); }
        else std::memset(sg, 0, 144);
        if (!tr.sp.empty()) std::memcpy(sp + X.spoff[j], tr.sp.data(), tr.sp.size() * 2);
        if (!tr.th.empty()) std::memcpy(th + X.thoff[j], tr.th.data(), tr.th.size() * 4);
      }
    });
    for (int z = 0; z < 8; z++) tok[4 * (size_t)X.ti + z] = 0;
    if (g_h250b && !X.cpuJ.empty()) {
      std::vector<H250bTree> v; v.reserve(X.cpuJ.size());
      for (uint32_t j : X.cpuJ) {
        const Tree& tr = UQ.U[j];
        const uint32_t cls = (tr.nS <= 4 ? 0u : tr.nS <= 8 ? 1u : tr.nS <= 16 ? 2u : 3u) + (tr.nth == 0 ? 4u : 0u);
        v.push_back({j, (int)tr.Q, cls >= 4u, (int)tr.nslot, (int)(4u << (cls & 3u)), (float)tr.Q * 9.0f * ldexpf(1.0f, -24)});
      }
      H250bIn in{ (const uint32_t*)X.hdr.b.contents, (const uint32_t*)X.tok.b.contents, (const uint16_t*)X.sp.b.contents,
                  (const float*)X.th.b.contents, (const uint32_t*)X.seg.b.contents, (float*)X.out.b.contents,
                  (const float*)P.bsob.contents, (const float*)X.yb.contents, (const float*)P.bct.contents };
      X.cpu_ms = 0; X.cpuT = fw_submit([in, v, pms = &X.cpu_ms]() mutable { h250b_cpu(in, v, fw_persist() ? 1 : g_h250bt, pms); });
    }
    T.h_s_fill += tk.lap();
  }
  T.scan_enc += ms_since(t0);


  Tick tc;
  size_t nsb = 0; for (int s = 0; s < 2; s++) nsb += P.S[s].sbins.size();
  auto h882_scan = [&E, &B, &P, &T, dev, M, nsb]() {
  if (nsb) {
    const size_t NLM = g_h916 ? (size_t)std::max(P.S[0].nl, P.S[1].nl) : (size_t)M;
    h900_ens(P.scr, "P.scr", dev, NLM * 16 * 16 * 4, MPV);
    size_t h139_tot = 0;
    for (int s = 0; s < 2; s++) { Side5& X = P.S[s]; X.h139n = 0; X.h139t = 0;
      if (g_h139) for (const SBin5& b : X.sbins) { const size_t t = h139_tgb(b);
        if (t > 32768) { h139_tot += ((t * b.cnt) + 255) & ~(size_t)255; X.h139n++; X.h139t += b.cnt; } } }
    if (h139_tot) h900_ens(P.h139scr, "P.h139scr", dev, h139_tot, MPV);
    if (g_h139 == 1) for (int s = 0; s < 2; s++) if (P.S[s].h139n) h900_ens(P.S[s].h139o, "P.S[s].h139o", dev, (size_t)M * 64, MSH);
    if (g_h335 == 2) for (int s = 0; s < 2; s++) h900_ens(P.S[s].h139o, "P.S[s].h139o", dev, (size_t)M * 64, MSH);
    for (int s = 0; s < 2; s++) P.S[s].h335b.clear();
    if (g_h333) { const int hsl = B.ci % 3; for (int s = 0; s < 2; s++) { Side5& X = P.S[s]; size_t tot = 0;
      for (const SBin5& b : X.sbins) if (!(b.s2 || h139_tgb(b) > 32768 || b.Q < 1 || b.Q > 16)) tot += b.cnt * h333_hstr(b.Q);
      h900_ens(g_h333hand[hsl][s], "g_h333hand[hsl][s]", dev, tot * 4, MPV); h900_ens(g_h333hoff[hsl][s], "g_h333hoff[hsl][s]", dev, (size_t)M * 4, MSH);
      std::memset(g_h333hoff[hsl][s].b.contents, 0xFF, (size_t)M * 4); g_h333_mb += (double)tot * 4 / 1048576.0; } }
    size_t h139_off = 0;
    id<MTLCommandBuffer> cb = h182_cb(E.cq); id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    H263Bind bc263; bc263.on = (g_h263 & 4) != 0;
    for (int s = 0; s < 2; s++) {
      Side5& X = P.S[s];
      const size_t h779n0 = X.h335b.size();
      for (size_t bj779 = 0; bj779 < X.sbins.size(); bj779++) { const SBin5& b = X.sbins[g_h779 ? X.sbins.size() - 1 - bj779 : bj779];
        Tick tp;
        const size_t tgb = h139_tgb(b); const bool over = tgb > 32768; const bool dpath = over && g_h139 == 2;
        id<MTLComputePipelineState> ps = dpath ? h139_ps(E, b.s2 ? 1 : 0, b.Q) : h93_scan_ps(E, b.s2 ? 1 : 0, b.Q);
        T.h_s_pipe += tp.lap();
        if (!ps) { fprintf(stderr, "missing scan kernel Q%d\n", b.Q); exit(1); }
        bc263.pipe(en, ps);
        id<MTLBuffer> bs[9] = {X.tok.b, X.hdr.b, X.sp.b, X.th.b, P.bsob, X.yb, X.out.b, X.bl.b, X.seg.b};
        for (int k = 0; k < 9; k++) bc263.buf(en, bs[k], (k == 7 ? b.off * 4 : 0), k);
        uint32_t ns = (uint32_t)b.ns, nc = (uint32_t)b.nc; float tau = b.Q * 9.0f * ldexpf(1.0f, -24);
        bc263.buf(en, P.bct, 0, 11); bc263.buf(en, P.scr.b, 0, 12);
        [en setBytes:&tau length:4 atIndex:13]; [en setBytes:&ns length:4 atIndex:9]; [en setBytes:&nc length:4 atIndex:10];
        uint32_t sstr = (uint32_t)(((size_t)ns * 96 + 2 * nc + 59) & ~(size_t)3); [en setBytes:&sstr length:4 atIndex:14];
        const int thr = b.s2 ? 96 / 3 : (b.Q <= 5 ? 768 / 6 : 768 / 3);
        if (g_h335 && !b.s2 && !over && b.Q >= 1 && b.Q <= 16) {
          bc263.pipe(en, h335_ps(E, 0, b.Q)); bc263.buf(en, P.scr.b, (size_t)s * NLM * 16 * 4, 12);
          [en setThreadgroupMemoryLength:h335_tgb(b, 0) atIndex:0];
          { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(h335_thr(0, b.Q), 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(h335_thr(0, b.Q), 1, 1)]; }
          X.h335b.push_back({b.off, b.cnt}); g_h335_bins += 1; g_h335_trees += (double)b.cnt;
          continue;
        }
        if (dpath) {
          const uint32_t stride = (uint32_t)(tgb / 4); const size_t reg = h139_off; h139_off += ((tgb * b.cnt) + 255) & ~(size_t)255;
          [en setBuffer:P.h139scr.b offset:reg atIndex:15]; [en setBytes:&stride length:4 atIndex:16];
          [en setThreadgroupMemoryLength:0 atIndex:0];
          { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; }
        } else {
          if (tgb > 32768) g_h139_viol += b.cnt;
          [en setThreadgroupMemoryLength:tgb atIndex:0];
          { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; }
          if (over && g_h139 == 1) {
            const uint32_t stride = (uint32_t)(tgb / 4); const size_t reg = h139_off; h139_off += ((tgb * b.cnt) + 255) & ~(size_t)255;
            [en setComputePipelineState:h139_ps(E, b.s2 ? 1 : 0, b.Q)]; bc263.forget_pipe();
            [en setBuffer:X.h139o.b offset:0 atIndex:6]; bc263.forget(6);
            [en setBuffer:P.h139scr.b offset:reg atIndex:15]; [en setBytes:&stride length:4 atIndex:16];
            [en setThreadgroupMemoryLength:0 atIndex:0];
            { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)]; }
          }
        }
      }
      if (g_h779) std::reverse(X.h335b.begin() + (long)h779n0, X.h335b.end());
    }
    if (g_h335 && (P.S[0].h335b.size() || P.S[1].h335b.size())) {
      [en memoryBarrierWithScope:MTLBarrierScopeBuffers];
      const int hsl = B.ci % 3;
      for (int s = 0; s < 2; s++) { Side5& X = P.S[s]; uint32_t hb333 = 0;
        for (const SBin5& b : X.sbins) {
          if (b.s2 || h139_tgb(b) > 32768 || b.Q < 1 || b.Q > 16) continue;
          id<MTLBuffer> bs[9] = {X.tok.b, X.hdr.b, X.sp.b, X.th.b, P.bsob, X.yb, X.out.b, X.bl.b, X.seg.b};
          bc263.pipe(en, h335_ps(E, 1, b.Q));
          for (int k = 0; k < 9; k++) bc263.buf(en, bs[k], (k == 7 ? b.off * 4 : 0), k);
          uint32_t ns = (uint32_t)b.ns, nc = (uint32_t)b.nc; float tau = b.Q * 9.0f * ldexpf(1.0f, -24);
          bc263.buf(en, P.bct, 0, 11); bc263.buf(en, P.scr.b, (size_t)s * NLM * 16 * 4, 12);
          [en setBytes:&tau length:4 atIndex:13]; [en setBytes:&ns length:4 atIndex:9]; [en setBytes:&nc length:4 atIndex:10];
          if (g_h333) {
            [en setBuffer:g_h333hand[hsl][s].b offset:0 atIndex:16];
 [en setBytes:&hb333 length:4 atIndex:17];
            uint32_t* HO = (uint32_t*)g_h333hoff[hsl][s].b.contents; const uint32_t* BL = (const uint32_t*)X.bl.b.contents; const uint32_t hs = (uint32_t)h333_hstr(b.Q);
            for (size_t k = 0; k < b.cnt; k++) { const uint32_t j = BL[b.off + k]; if (j < (uint32_t)M) HO[j] = hb333 + (uint32_t)k * hs; }
            hb333 += (uint32_t)b.cnt * hs; g_h333_trees += (double)b.cnt; }
          [en setThreadgroupMemoryLength:h335_tgb(b, 1) atIndex:0];
          { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(h335_thr(1, b.Q), 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(h335_thr(1, b.Q), 1, 1)]; }
          if (g_h335 == 2) {
            bc263.pipe(en, h93_scan_ps(E, 0, b.Q)); bc263.buf(en, X.h139o.b, 0, 6); bc263.buf(en, P.scr.b, 0, 12);
            [en setThreadgroupMemoryLength:h139_tgb(b) atIndex:0];
            { [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(b.Q <= 5 ? 768 / 6 : 768 / 3, 1, 1)]; if (g_h660d && !strcmp(g_h660d, "scan2")) [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(b.Q <= 5 ? 768 / 6 : 768 / 3, 1, 1)]; }
          }
        }
      }
    }
    [en endEncoding]; [cb commit];
    B.cbS = cb;
    if (g_h139 && (P.S[0].h139n || P.S[1].h139n))
      printf("H139 commit mode %d: over-limit bins %zu (TGs %zu), TG-memory violations so far %llu\n", g_h139,
             P.S[0].h139n + P.S[1].h139n, P.S[0].h139t + P.S[1].h139t, (unsigned long long)g_h139_viol);
  }
  };
  if (g_h882 & 1) B.t882S = g_h882q->submit(0, h882_scan); else h882_scan();
  T.h_s_cmd += tc.lap();
}


struct H770Hit { int t, id; uint8_t v; int8_t vs; uint8_t st; unsigned char raw[2][128]; };
struct H770B { std::vector<H770Hit> hits; std::vector<int32_t> hix[2]; };
struct H770Stat {
  std::atomic<long> pairs{0}, hits{0}, full{0}, partial{0}, dropg{0}, dropc{0}, nonbin{0}, rowdis{0}, ver{0}, veq{0}, vdiff{0}, vtail{0}, vnos{0}, put_cpu{0}, put_gpu{0};
  std::atomic<long> dw[9]; std::atomic<long> refq[17];
};
static H770Stat g_h770s;
static void h770_apply_raw(PipeStep& PS, int t, int id, const unsigned char* raw, Times& T);

static void h770_list(Persist5& P, H770B& H) {
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; size_t n0 = X.ls.size(); bool any = false;
    for (const H770Hit& h : H.hits) if (h.t == s && (size_t)h.id < X.listed.size()) { X.listed[(size_t)h.id] = 1; if (g_h215 & 1) X.ls.push_back(h.id); any = true; }
    if (any && (g_h215 & 1)) { std::sort(X.ls.begin() + (long)n0, X.ls.end()); std::inplace_merge(X.ls.begin(), X.ls.begin() + (long)n0, X.ls.end()); }
  }
}

static bool h770_cpu_hit(H770B& H, int s, uint32_t j, bool verify) {
  if ((size_t)j >= H.hix[s].size() || H.hix[s][j] < 0) return false;
  H770Hit& h = H.hits[(size_t)H.hix[s][j]]; h.vs = 1;
  if (verify) return false;
  if (h.v & 2) { h.st = 1; g_h770s.dropc++; return true; }
  h.st = 2; return false;
}

static void h770_gpu_part(Persist5& P, H770B& H, const std::vector<size_t>& lo, bool verify, size_t& ndrop) {
  for (size_t b = 0; b < P.gbins.size(); b++) {
    const GBin5& g = P.gbins[b]; Side5& X = P.S[g.side];
    if (H.hix[g.side].empty()) continue;
    uint32_t* gbl = (uint32_t*)X.gbl.b.contents;
    const size_t kend = (b < lo.size() && lo[b] != SIZE_MAX) ? lo[b] : g.cnt;
    for (size_t k = 0; k < kend; k++) {
      const uint32_t j = X.sorted[g.sst + k];
      if ((size_t)j >= H.hix[g.side].size() || H.hix[g.side][j] < 0) continue;
      H770Hit& h = H.hits[(size_t)H.hix[g.side][j]]; h.vs = 0;
      if (verify) continue;
      if (h.v & 1) { h.st = 1; gbl[g.tl_off + 1 + k] = 0xFFFFFFFFu; ndrop++; g_h770s.dropg++; } else h.st = 2;
    }
  }
}

static void h770_finish(PipeStep& PS, Batch& B, H770B& H, Times& T) {
  if (g_h770v) return;
  for (H770Hit& h : H.hits) {
    if (h.vs < 0) { h.st = 1; g_h770s.nonbin++; }
    if (h.st == 1) { h770_apply_raw(PS, h.t, h.id, h.raw[h.vs >= 0 ? h.vs : ((h.v & 1) ? 0 : 1)], T); g_h770s.full++; }
    else if (h.st == 2) { B.need[h.t].push_back(h.id); g_h770s.partial++; }
  }
}

static void gn_begin_p(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;
  const int M = B.M;
  H770B* H7 = (g_h770 && B.h770 && !B.h770->hits.empty()) ? B.h770.get() : nullptr;
  if (H7 && !g_h770v) h770_list(*B.P, *H7);
  if ((!g_h217 && (int)PS.SB->UQ.U.size() < M) || ((g_h215 & 1) ? ((int)B.P->S[0].listed.size() < M || (int)B.P->S[1].listed.size() < M) : ((int)B.P->S[0].listed.size() != M || (int)B.P->S[1].listed.size() != M))) { fprintf(stderr, "pipeline: slot arrays do not match chunk %d\n", B.b0); exit(1); }
  Persist5& P = *B.P; id<MTLDevice> dev = E.dev;
  const int NCH = 4 * (int)WK(E).size();
  auto span = [M, NCH](int c, int& j0, int& j1) { j0 = (int)((long)M * c / NCH); j1 = (int)((long)M * (c + 1) / NCH); };
  auto t0 = clk::now();
  Tick tk;
  P.gbins.clear();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    const bool l215 = (g_h215 & 1) != 0; const int N215 = l215 ? (int)X.ls.size() : M;
    auto span = [N215, NCH](int c, int& j0, int& j1) { j0 = (int)((long)N215 * c / NCH); j1 = (int)((long)N215 * (c + 1) / NCH); };
    X.gok.resize(M); X.gninst.resize(M); X.gper.resize(M); X.gtoff.resize(M); X.gspoff.resize(M); X.keyl.resize(M); X.dkey.resize(M);
    T.h_g_alloc += tk.lap();

    H263GnScr ownG; H263GnScr& sg = (g_h263 & 1) ? h263_keep_gn : ownG;
    std::vector<uint32_t>& bgi = sg.bgi; std::vector<uint32_t>& bgs = sg.bgs; std::vector<uint32_t>& bns = sg.bns;
    std::vector<int>& bqg = sg.bqg; std::vector<int>& bcg = sg.bcg; std::vector<int>& blm = sg.blm;
    std::vector<double>& bal = sg.bal; std::vector<double>& bkm = sg.bkm;
    bgi.assign(NCH, 0); bgs.assign(NCH, 0); bns.assign(NCH, 0);
    bqg.assign(NCH, 0); bcg.assign(NCH, 0); blm.assign(NCH, 0);
    bal.assign(NCH, 0.0); bkm.assign(NCH, 0.0);
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      uint32_t gi = 0, gs = 0, ns = 0; int qg = 0, cg = 0, lm = 0; double al = 0, km = 0;
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        X.gok[j] = 0; X.gninst[j] = 0; X.gper[j] = 0;
        if (!X.listed[j]) continue;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (tr.nth >= 1 && tr.nth > GN_TH_MAX) km += 1;
        if (te.g_ok && isfinite(o.smin[j])) {
          X.gok[j] = 1; o.is_g[j] = 1; al += 1;
          X.gninst[j] = (uint32_t)(te.ngp() / 4); X.gper[j] = te.gper; gi += X.gninst[j];
          ns += 1; qg = std::max(qg, tr.Q); cg = std::max(cg, (int)((X.gper[j] + 63) / 64)); lm = std::max(lm, (int)h313a::nprog(tr));
        }
        gs += (uint32_t)tr.sp.size();
      }
      bgi[c] = gi; bgs[c] = gs; bns[c] = ns; bqg[c] = qg; bcg[c] = cg; blm[c] = lm; bal[c] = al; bkm[c] = km;
    });
    std::vector<uint32_t>& ogt = sg.ogt; std::vector<uint32_t>& ogs = sg.ogs; std::vector<uint32_t>& osq = sg.osq;
    ogt.resize(NCH); ogs.resize(NCH); osq.resize(NCH);
    uint32_t GTI = 0, GSI = 0, NS = 0; int Qg = 0, Cg = 0, LM = 0;
    for (int c = 0; c < NCH; c++) {
      ogt[c] = GTI; ogs[c] = GSI; osq[c] = NS;
      GTI += bgi[c]; GSI += bgs[c]; NS += bns[c];
      Qg = std::max(Qg, bqg[c]); Cg = std::max(Cg, bcg[c]); LM = std::max(LM, blm[c]);
      T.g_alive += bal[c]; T.gn_kmiss += bkm[c];
    }
    X.gti = GTI; X.gsi = GSI;
    X.seq.resize(NS);
    if (Qg >= 100) { fprintf(stderr, "GN: Q %d >= 100 breaks the bin key\n", Qg); exit(1); }

    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      uint32_t gt = ogt[c], gs = ogs[c], sq = osq[c];
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        X.gtoff[j] = gt; X.gspoff[j] = gs; X.keyl[j] = 0; X.dkey[j] = 0;
        if (!X.listed[j]) continue;
        const Tree& tr = UQ.U[j];
        if (X.gok[j]) { X.seq[sq++] = (uint32_t)j; X.keyl[j] = (uint32_t)(LM - (int)h313a::nprog(tr)); gt += X.gninst[j]; }
        gs += (uint32_t)tr.sp.size();
      }
    });
    T.h_g_pre += tk.lap();
    sort_by_key5_p(E, X.seq, X.keyl, (size_t)LM + 1, X.start, X.seq2);
    T.h_g_sort += tk.lap();
    const size_t DNg = (size_t)(Cg + 1) * (size_t)(Qg + 1) * 25;
    PL(E).run(NCH, [&](int c) {
      const size_t n2 = X.seq2.size(), a = n2 * c / NCH, z = n2 * (c + 1) / NCH;
      for (size_t i = a; i < z; i++) {
        const uint32_t j = X.seq2[i]; const Tree& tr = UQ.U[j]; const size_t cls = (X.gper[j] + 63) / 64;
        X.dkey[j] = (uint32_t)((cls * (size_t)(Qg + 1) + (size_t)tr.Q) * 25 + (size_t)tr.nth);
      }
    });
    T.h_g_plan += tk.lap();
    sort_by_key5_p(E, X.seq2, X.dkey, DNg, X.start, X.sorted);
    T.h_g_sort += tk.lap();

    std::vector<size_t>& ne = sg.ne; ne.clear(); for (size_t dk = 0; dk < DNg; dk++) if (X.start[dk + 1] > X.start[dk]) ne.push_back(dk);
    std::vector<size_t>& perb = sg.perb; std::vector<size_t>& perb2 = sg.perb2; perb.assign(ne.size(), 0); perb2.assign(ne.size(), 0); std::vector<size_t> perb3(ne.size(), 0);
    PL(E).run(NCH, [&](int c) {
      const size_t a = ne.size() * c / NCH, z = ne.size() * (c + 1) / NCH;
      for (size_t bi = a; bi < z; bi++) {
        size_t per = 0;
        for (size_t k = X.start[ne[bi]]; k < X.start[ne[bi] + 1]; k++) per = std::max<size_t>(per, X.gper[X.sorted[k]]);
        size_t pen = 0, pen32 = 0; for (size_t k = X.start[ne[bi]]; k < X.start[ne[bi] + 1]; k++) { const size_t p32_ = h141_pe(UQ.enc[X.sorted[k]], X.gper[X.sorted[k]]); pen32 = std::max<size_t>(pen32, p32_); pen = std::max<size_t>(pen, h332_pe(UQ.enc[X.sorted[k]], (uint32_t)p32_)); }
        perb[bi] = per; perb2[bi] = pen; perb3[bi] = pen32;
      }
    });
    for (size_t bi = 0; bi < ne.size(); bi++) {
      const size_t dk = ne[bi];
      const size_t c = X.start[dk + 1] - X.start[dk];
      const size_t sst = X.start[dk];
      const int Tk = (int)(dk % 25); const size_t rem = dk / 25; const int Qk = (int)(rem % (size_t)(Qg + 1));
      Tick tp;
      id<MTLComputePipelineState> ps = h93_gn_ps(E, Qk, Tk);
      T.h_g_pipe += tp.lap();
      if (!ps) { T.gn_kmiss += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      const size_t per = perb[bi]; const size_t pen = perb2[bi]; const size_t pen32 = perb3[bi];
      const int NPk = std::max(4, (Tk + 3) / 4 * 4), QPk = std::max(4, (Qk + 3) / 4 * 4); const size_t WSk = (4 * NPk + 3 * QPk + 6) / 4 * 4;
      const size_t pers = WSk + 2 * per; const size_t shb = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
      if (pers * 4 > 32768) { T.gn_drop += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      const size_t persN = WSk + 2 * pen;
      const int cap141 = (int)ps.maxTotalThreadsPerThreadgroup / 32;
      int tpg = h141_pick(shb, persN, c, cap141);
      h141_note(shb, pers, persN, c, cap141, tpg);
      const size_t xf713 = h713_gn(E, Qk, Tk, shb, persN, ps, tpg);
      GBin5 g{}; g.side = s; g.Q = Qk; g.nth = Tk; g.tpg = tpg; g.cnt = c;
      g.nt = (c + (size_t)tpg - 1) / (size_t)tpg * (size_t)tpg;
      g.ntg = xf713 ? g.nt : g.nt / (size_t)tpg;
      g.tgb = ((shb + persN * tpg + xf713) * 4 + 15) & ~(size_t)15; g.per = pen; g.per32 = pen32; g.sst = sst; g.ps = ps;
      g.eb_bytes = g.ntg * (size_t)tpg * 2 * (size_t)(Qk + 1) * 64 * 4;
      P.gbins.push_back(g);
    }
    T.h_g_plan += tk.lap();
  }

  if (h94_mode() >= 1) {
    std::vector<size_t> nat, tgb, cnt, cap; size_t before = 0, after = 0, ntr = 0;
    for (const GBin5& g : P.gbins) { nat.push_back(g.ntg); tgb.push_back(g.tgb); cnt.push_back(g.cnt); before += g.ntg; ntr += g.cnt; }
    cap.assign(nat.size(), 0); h94_cap(cap, nat, tgb, cnt);
    for (size_t b = 0; b < P.gbins.size(); b++) { GBin5& g = P.gbins[b]; g.ntg = cap[b]; g.eb_bytes = g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4; after += g.ntg; }
    fprintf(stderr, "H94 GN dispatch: bins %zu trees %zu TGs before %zu after %zu\n", P.gbins.size(), ntr, before, after);
  }
  size_t tlc[2] = {0, 0}, ebc = 0; uint32_t NB = 0;
  for (GBin5& g : P.gbins) {
    g.wix = NB++;
    g.tl_off = tlc[g.side]; tlc[g.side] = al5(tlc[g.side] + 1 + g.nt, BIN_ALIGN);
    g.eb_off = ebc; ebc = al5(ebc + g.eb_bytes + ebt_bytes(g.Q, g.nth, g.ntg, g.tpg), 256);
  }
  T.h_g_plan += tk.lap();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    const bool l215 = (g_h215 & 1) != 0; const int N215 = l215 ? (int)X.ls.size() : M;
    auto span = [N215, NCH](int c, int& j0, int& j1) { j0 = (int)((long)N215 * c / NCH); j1 = (int)((long)N215 * (c + 1) / NCH); };
    if (g_h916) { X.gnl = (uint32_t)X.ls.size(); X.glid.resize((size_t)M); for (size_t k = 0; k < X.ls.size(); k++) X.glid[(size_t)X.ls[k]] = (uint32_t)k; }
    const size_t TBN = g_h916 ? (size_t)X.gnl : (size_t)M;
    h900_ens(X.gtok, "X.gtok", dev, ((size_t)X.gti * 4 + 8) * 4, MSH);
    h900_ens(X.ghdr, "X.ghdr", dev, TBN * 16, MSH);
    h900_ens(X.gsp, "X.gsp", dev, (size_t)X.gsi * 2, MSH);
    h900_ens(X.gseg, "X.gseg", dev, TBN * 144, MSH);
    h900_ens(X.gth, "X.gth", dev, TBN * 128, MSH);
    h900_ens(X.grs, "X.grs", dev, TBN * 8, MSH);
    h900_ens(X.gst, "X.gst", dev, TBN * 4, MSH);
    h900_ens(X.gai, "X.gai", dev, TBN * 64, MSH);
    h900_ens(X.gbl, "X.gbl", dev, tlc[s] * 4, MSH);
    T.h_g_alloc += tk.lap();
    uint32_t* gtok = (uint32_t*)X.gtok.b.contents; uint32_t* ghdr = (uint32_t*)X.ghdr.b.contents; uint16_t* gsp = (uint16_t*)X.gsp.b.contents;
    uint32_t* gseg = (uint32_t*)X.gseg.b.contents; float* gth = (float*)X.gth.b.contents;
    PL(E).run(NCH, [&](int c) {
      int j0, j1; span(c, j0, j1);
      for (int jj = j0; jj < j1; jj++) { const int j = l215 ? X.ls[jj] : jj;
        const size_t lj = glid916(X, (size_t)j); uint32_t* h = ghdr + 4 * lj; float* th = gth + 32 * lj; uint32_t* sg = gseg + 36 * lj;
        for (int q = 0; q < 32; q++) th[q] = 0.0f;
        h[0] = X.gtoff[j]; h[2] = X.gspoff[j];
        if (!X.listed[j]) { h[1] = 0; h[3] = 0; std::memset(sg, 0, 144); continue; }
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.gninst[j] | ((uint32_t)tr.nS << 16);
        if (X.gok[j]) {
          h[3] = te.ghw; std::memcpy(sg, te.gg, 144);
          std::memcpy(gtok + 4 * (size_t)X.gtoff[j], te.pgp(), te.ngp() * 4);
          double th2[32]; start_theta(tr.th, tr.nth, PS.SB->F[s].sarg[j], th2);
          for (int q = 0; q < tr.nth; q++) th[q] = (float)th2[q];
        } else { h[3] = 0; std::memset(sg, 0, 144); }
        if (!tr.sp.empty()) std::memcpy(gsp + X.gspoff[j], tr.sp.data(), tr.sp.size() * 2);
      }
    });
    for (int z = 0; z < 8; z++) gtok[4 * (size_t)X.gti + z] = 0;
    T.h_g_fill += tk.lap();
    uint32_t* gbl = (uint32_t*)X.gbl.b.contents;
    std::vector<const GBin5*> sb; for (const GBin5& g : P.gbins) if (g.side == s) sb.push_back(&g);
    if (!sb.empty()) PL(E).run((int)sb.size(), [&](int i) {
      const GBin5& g = *sb[i]; uint32_t* tl = gbl + g.tl_off; tl[0] = (uint32_t)g.per;
      for (size_t k = 0; k < g.cnt; k++) tl[1 + k] = (uint32_t)glid916(X, X.sorted[g.sst + k]);
      for (size_t k = g.cnt; k < g.nt; k++) tl[1 + k] = 0xFFFFFFFFu;
    });
    T.h_g_fill += tk.lap();
  }
  h900_ens(P.eb, "P.eb", dev, ebc, MPV);
  h900_ens(P.ctr, "P.ctr", dev, 4 * (size_t)NB + 16, MSH);
  T.h_g_alloc += tk.lap();
  T.gn_enc += ms_since(t0);


  if (g_h250_pending && !g_h722) h250a_join(*g_h250_pending);
  B.h250_items.clear(); B.h250_ngpu = 0; B.h250_on = false; B.h250_cpu_ms = 0; B.h334_gc.clear();
  std::vector<size_t> h770lo; size_t h770ngd = 0;
  if (g_h770) for (int s_ = 0; s_ < 2; s_++) B.h770cpu[s_].assign(UQ.U.size(), 0);
  if (H7) h770lo.assign(P.gbins.size(), SIZE_MAX);
  if (g_h250a && NB) {
    size_t ntot = 0; for (const GBin5& g : P.gbins) ntot += g.cnt;
    const int pc = std::min(B.ci, 15);
    const bool own = g_h300b && g_h300b_rg[pc] > 0 && g_h300b_rc[pc] > 0;
    const double f = h242d_f(own ? g_h300b_rc[pc] / (g_h300b_rc[pc] + g_h300b_rg[pc])
                             : (g_h250a_rg > 0 && g_h250a_rc > 0) ? g_h250a_rc / (g_h250a_rc + g_h250a_rg) : 0.5);
    const double f404 = g_h648 ? std::min(1.0, std::max(0.0, f * g_h648_s)) : g_h404gn >= 0 ? std::min(1.0, std::max(0.0, f * g_h404gn)) : f;
    const size_t ncpu = (size_t)llround(f404 * (double)ntot);
    B.h250_f = f404;
    if (g_h334) { B.h334_gc.resize(P.gbins.size()); for (size_t b = 0; b < P.gbins.size(); b++) B.h334_gc[b] = P.gbins[b].cnt; }
    size_t h770n = 0;
    for (int b = (int)P.gbins.size() - 1; b >= 0 && h770n < ncpu; b--) {
      const GBin5& g = P.gbins[(size_t)b]; Side5& X = P.S[g.side];
      uint32_t* gbl = (uint32_t*)X.gbl.b.contents;
      for (size_t k = g.cnt; k-- > 0 && h770n < ncpu;) {
        gbl[g.tl_off + 1 + k] = 0xFFFFFFFFu;
        if (g_h334) B.h334_gc[(size_t)b] = k;
        h770n++;
        if (g_h770) { const uint32_t j770 = X.sorted[g.sst + k]; if ((size_t)j770 < B.h770cpu[g.side].size()) B.h770cpu[g.side][j770] = 1; }
        if (H7) { h770lo[(size_t)b] = k; if (h770_cpu_hit(*H7, g.side, X.sorted[g.sst + k], g_h770v != 0)) continue; }
        B.h250_items.push_back(H250aItem{g.side, g.Q, g.nth, E.tg[g.side].floor_, (uint32_t)glid916(X, X.sorted[g.sst + k]), (uint32_t)g.per32});
      }
    }
    B.h250_ngpu = ntot - h770n;
  }
  if (H7) {
    h770_gpu_part(P, *H7, h770lo, g_h770v != 0, h770ngd);
    if (g_h250a && NB) B.h250_ngpu = B.h250_ngpu >= h770ngd ? B.h250_ngpu - h770ngd : 0;
    h770_finish(PS, B, *H7, T);
  }

  Tick tc;
  auto h882_gn = [&E, &B, &P, NB]() {
  if (NB) {
    std::memset(P.ctr.b.contents, 0, P.ctr.cap);
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    H263Bind bg263; bg263.on = (g_h263 & 4) != 0;
    for (size_t gj779 = 0; gj779 < P.gbins.size(); gj779++) { const size_t gi334 = g_h779 ? P.gbins.size() - 1 - gj779 : gj779; const GBin5& g = P.gbins[gi334];
      Side5& X = P.S[g.side];
      const bool t334 = g_h334 && gi334 < B.h334_gc.size();
      const size_t gc334 = t334 ? B.h334_gc[gi334] : g.nt;
      if (t334 && gc334 == 0) continue;
      const size_t ntg334 = t334 ? std::min(g.ntg, (gc334 + (size_t)g.tpg - 1) / (size_t)g.tpg) : g.ntg;
      if (h182_on()) [en pushDebugGroup:h182_gname(g)];
      bg263.pipe(en, g.ps);
      id<MTLBuffer> bs[9] = {X.gtok.b, X.ghdr.b, X.gsp.b, X.gseg.b, X.yb, X.gth.b, X.grs.b, X.gst.b, X.gbl.b};
      for (int k = 0; k < 9; k++) bg263.buf(en, bs[k], (k == 8 ? g.tl_off * 4 : 0), k);
      float tf[4] = {g.Q * 9.0f * ldexpf(1.0f, -24), E.tg[g.side].floor_, (g_freet ? 0.0f : 3.0f), g_h333 ? 1.0f : 0.0f};
      if (g_h333) { const int hsl = B.ci % 3; bg263.buf(en, g_h333hand[hsl][g.side].b, 0, 14); bg263.buf(en, g_h333hoff[hsl][g.side].b, 0, 15); }
      [en setBytes:tf length:16 atIndex:9];
      bg263.buf(en, X.gai.b, 0, 10);
      bg263.buf(en, P.eb.b, g.eb_off, 11);
      if (g_h711) bg263.buf(en, P.eb.b, g.eb_off + g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4, 15);
      uint32_t wq[2] = {g.wix, (uint32_t)gc334}; [en setBytes:wq length:8 atIndex:12];
      bg263.buf(en, P.ctr.b, 0, 13);
      [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      if (h150c_on()) { g_h150.gnDisp++; g_h150.gnLaunch += (long)(g.ntg * (size_t)g.tpg); g_h150.gnOcc += (long)g.cnt; }
      [en dispatchThreadgroups:MTLSizeMake(ntg334, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
      if (h182_on()) [en popDebugGroup];
    }
    [en endEncoding]; [cb commit];
    B.cbG = cb;
  }
  };
  if ((g_h882 & 2) && NB) B.t882G = g_h882q->submit(1, h882_gn); else h882_gn();
  T.h_g_cmd += tc.lap();
  if (g_h722 && g_h250_pending && g_h250_pending != &B) h250a_join(*g_h250_pending);
  if (g_h250a && NB && !B.h250_items.empty()) {
    B.h250_on = true; g_h250_pending = &B;
    B.h250_thr = fw_submit([&B, &P]() { h250a_body(&B, &P, &B.h250_cpu_ms); });
  }
}


struct Fin14 {
  std::vector<float> th[2], a[2];
  void init(size_t nu) {
    static const int h742 = getenv("H742") ? atoi(getenv("H742")) : 0;
    if (h742 == 1 && th[0].size() == nu * 32) return;
    const float v = h742 == 2 ? __builtin_bit_cast(float, 0xFFFFFFFFu) : 0.0f;
    for (int s = 0; s < 2; s++) { th[s].assign(nu * 32, v); a[s].assign(nu * 16, v); } }
};
static Fin14 F14;


struct FitMemo99 {
  static constexpr uint32_t NONE = 0xFFFFFFFFu;
  size_t cap = 0, arena_w = 0;
  size_t aptr = 0, aold = 0;
  uint32_t hd = 0, tl = 0, nlive = 0;
  std::vector<uint32_t> bucket, nxt, eoff, elen, ett, arena;
  std::vector<uint64_t> eh;
  std::vector<float> rf, sm, th, ap;
  std::vector<int32_t> sa; std::vector<uint8_t> ig; std::vector<uint32_t> st;
  long lookups = 0, hits = 0, inserts = 0, evicts = 0, refused = 0, dups = 0; long keyw = 0;
  void init(size_t entries, size_t words) {
    cap = entries; arena_w = words;
    bucket.assign(cap, NONE); nxt.assign(cap, NONE); eoff.assign(cap, 0); elen.assign(cap, 0); ett.assign(cap, 0); eh.assign(cap, 0);
    arena.assign(arena_w, 0);
    rf.assign(cap, 0.0f); sm.assign(cap, 0.0f); th.assign(cap * 32, 0.0f); ap.assign(cap * 16, 0.0f);
    sa.assign(cap, -1); ig.assign(cap, 0); st.assign(cap, 0);
  }
  size_t bytes() const {
    return cap * (4 * 5 + 8 + 4 + 4 + 4 + 1 + 4 + 4 * 32 + 4 * 16) + arena_w * 4;
  }
  uint32_t find(uint64_t h, uint32_t t, const uint32_t* k, uint32_t n) const {
    for (uint32_t e = bucket[(size_t)(h & (cap - 1))]; e != NONE; e = nxt[e])
      if (eh[e] == h && ett[e] == t && elen[e] == n && std::memcmp(arena.data() + eoff[e], k, 4 * (size_t)n) == 0) return e;
    return NONE;
  }
  void evict_head() {
    const uint32_t e = hd, b = (uint32_t)(eh[e] & (cap - 1));
    if (bucket[b] == e) bucket[b] = nxt[e];
    else { uint32_t p = bucket[b]; while (nxt[p] != e) p = nxt[p]; nxt[p] = nxt[e]; }
    hd = (uint32_t)((hd + 1) % cap); nlive--; evicts++;
    aold = nlive ? eoff[hd] : aptr;
  }

  bool overlaps(size_t s0, size_t n) const {
    if (nlive == 0) return false;
    if (aold < aptr) return s0 < aptr && s0 + n > aold;
    return s0 + n > aold || s0 < aptr;
  }

  uint32_t insert(uint64_t h, uint32_t t, const uint32_t* k, uint32_t n) {
    if (n == 0 || (size_t)n > arena_w / 8) { refused++; return NONE; }
    if (nlive == cap) evict_head();
    size_t s0 = 0;
    for (;;) {
      if (nlive == 0) { aptr = 0; aold = 0; }
      s0 = (aptr + n <= arena_w) ? aptr : 0;
      if (!overlaps(s0, n)) break;
      evict_head();
    }
    const uint32_t slot = tl; tl = (uint32_t)((tl + 1) % cap); nlive++; inserts++; keyw += n;
    eh[slot] = h; ett[slot] = t; eoff[slot] = (uint32_t)s0; elen[slot] = n;
    std::memcpy(arena.data() + s0, k, 4 * (size_t)n); aptr = s0 + n;
    const uint32_t b = (uint32_t)(h & (cap - 1)); nxt[slot] = bucket[b]; bucket[b] = slot;
    return slot;
  }
};
static FitMemo99 G99;
static uint64_t h99_key(uint64_t hr, uint32_t t) { return splitmix64(hr ^ (0x9e3779b97f4a7c15ull * (uint64_t)(t + 1))); }

static inline uint64_t h770_key(uint64_t hr, uint32_t t, int cpu) { return splitmix64(hr ^ (0x9e3779b97f4a7c15ull * (uint64_t)(t + 1)) ^ (cpu ? 0xD1B54A32D192ED03ull : 0ull)); }

static uint32_t h99_find(uint64_t hr, uint32_t t, const uint32_t* k, uint32_t n) { return G99.find(h99_key(hr, t), t, k, n); }
static void h99_lookup_count() { G99.lookups++; }
static void h99_hit_count() { G99.hits++; }

static void h99_apply(PipeStep& PS, int t, int id, uint32_t e, Times& T) {
  FitOut& o = PS.SB->F[t];
  o.rfin[id] = G99.rf[e]; o.smin[id] = G99.sm[e]; o.sarg[id] = G99.sa[e]; o.is_g[id] = G99.ig[e];
  std::memcpy(&F14.th[t][32 * (size_t)id], &G99.th[32 * (size_t)e], 32 * sizeof(float));
  std::memcpy(&F14.a[t][16 * (size_t)id], &G99.ap[16 * (size_t)e], 16 * sizeof(float));
  if (G99.ig[e]) { T.steps_sum += (double)G99.st[e]; T.steps_n += 1; }
}

static void h99_put(PipeStep& PS, int t, int j, const FitOut& o, uint32_t steps) {
  const Tree& tr = PS.SB->UQ.U[(size_t)j];
  const uint32_t n = (uint32_t)tr.rec.size();
  const uint64_t h = h99_key(PS.SB->UQ.hash[(size_t)j], (uint32_t)t);
  if (G99.find(h, (uint32_t)t, tr.rec.data(), n) != FitMemo99::NONE) { G99.dups++; return; }
  const uint32_t e = G99.insert(h, (uint32_t)t, tr.rec.data(), n);
  if (e == FitMemo99::NONE) return;
  G99.rf[e] = o.rfin[j]; G99.sm[e] = o.smin[j]; G99.sa[e] = o.sarg[j]; G99.ig[e] = o.is_g[j]; G99.st[e] = steps;
  std::memcpy(&G99.th[32 * (size_t)e], &F14.th[t][32 * (size_t)j], 32 * sizeof(float));
  std::memcpy(&G99.ap[16 * (size_t)e], &F14.a[t][16 * (size_t)j], 16 * sizeof(float));
}


#include <deque>
#include <array>
static constexpr uint32_t H242_ROW = 26;
struct alignas(128) Slot242 { uint64_t key; uint32_t meta; float rfin, smin; int32_t sarg; float row[H242_ROW]; };
static_assert(sizeof(Slot242) == 128, "H242 slot must be one 128-byte line");
struct Stage242 { uint64_t key; uint32_t meta, t, nrec_off, nrec; float rfin, smin; int32_t sarg; float row[H242_ROW]; };
struct PutJob242 { std::vector<Stage242> st; std::vector<uint32_t> recs;
  PipeStep* ps = nullptr; int side = 0; std::vector<int> ids; std::vector<uint32_t> steps; std::vector<uint8_t> cpu; };
struct Memo242 {
  Slot242* tab = nullptr; size_t cap = 0, mask = 0, maxlive = 0, nlive = 0; bool verify = false, sync = false;
  std::vector<uint32_t> vdat; std::vector<uint64_t> voff; std::vector<uint32_t> vlen;
  long lookups = 0, hits = 0, nputs = 0, inserts = 0, dups = 0, refused_wide = 0, clears = 0, vmis = 0, vcoll = 0, njobs = 0;
  double lookup_ms = 0, apply_ms = 0, put_ms = 0, worker_ms = 0;
  std::mutex tab_m;
  std::mutex qm; std::condition_variable qcv, qdone;
  std::deque<PutJob242> q; bool busy = false, stop = false, started = false; std::thread th;
  void init(size_t entries, bool ver) {
    cap = entries; mask = cap - 1; maxlive = cap / 4 * 3; verify = ver;
    void* p = nullptr;
    if (posix_memalign(&p, 128, cap * sizeof(Slot242)) != 0 || !p) { fprintf(stderr, "H242: memo allocation failed\n"); std::exit(1); }
    tab = (Slot242*)p; std::memset((void*)tab, 0, cap * sizeof(Slot242));
    if (verify) { voff.assign(cap, 0); vlen.assign(cap, 0); }
  }
  void clear() { std::memset((void*)tab, 0, cap * sizeof(Slot242)); nlive = 0; vdat.clear(); clears++; }

  uint32_t probe(uint64_t key, uint32_t t) const {
    for (size_t i = (size_t)key & mask;; i = (i + 1) & mask) {
      const Slot242& s = tab[i];
      if (!(s.meta & 1u)) return H99_NONE;
      if (s.key == key && ((s.meta >> 21) & 1u) == t) return (uint32_t)i;
    }
  }
};
static Memo242 G242;
static inline uint32_t h242_bits(float f) { uint32_t u; std::memcpy(&u, &f, 4); return u; }
static inline float h242_flt(uint32_t u) { float f; std::memcpy(&f, &u, 4); return f; }
static bool h242_recmatch(uint32_t e, const uint32_t* k, uint32_t n) {
  return G242.vlen[e] == n && std::memcmp(G242.vdat.data() + G242.voff[e], k, 4 * (size_t)n) == 0;
}

static void h242_insert(const Stage242& g, const std::vector<uint32_t>& recs) {
  const uint32_t f = G242.probe(g.key, g.t);
  if (f != H99_NONE) {
    if (G242.verify && !(G242.vlen[f] == g.nrec && std::memcmp(G242.vdat.data() + G242.voff[f], recs.data() + g.nrec_off, 4 * (size_t)g.nrec) == 0)) G242.vcoll++;
    else G242.dups++;
    return;
  }
  if (G242.nlive >= G242.maxlive) G242.clear();
  size_t i = (size_t)g.key & G242.mask;
  while (G242.tab[i].meta & 1u) i = (i + 1) & G242.mask;
  Slot242& s = G242.tab[i];
  s.key = g.key; s.meta = g.meta; s.rfin = g.rfin; s.smin = g.smin; s.sarg = g.sarg;
  std::memcpy(s.row, g.row, sizeof(g.row));
  if (G242.verify) { G242.voff[i] = G242.vdat.size(); G242.vlen[i] = g.nrec; G242.vdat.insert(G242.vdat.end(), recs.begin() + g.nrec_off, recs.begin() + g.nrec_off + g.nrec); }
  G242.nlive++; G242.inserts++;
  if (g_h413dbg) { auto mx = [](uint64_t& a, uint64_t v) { uint64_t z = (v + 0x9E3779B97F4A7C15ull) * 0xBF58476D1CE4E5B9ull; a += z ^ (z >> 31); }; uint32_t w; mx(g_h413_f[0], g.key); mx(g_h413_f[1], g.meta & ~(255u << 13)); mx(g_h413_f[5], (g.meta >> 13) & 255u); std::memcpy(&w, &g.rfin, 4); mx(g_h413_f[2], w); std::memcpy(&w, &g.smin, 4); mx(g_h413_f[3], w ^ ((uint64_t)(uint32_t)g.sarg << 32)); { uint64_t rh = g.key; for (uint32_t q = 0; q < ((g.meta >> 2) & 63u) + ((g.meta >> 8) & 31u) && q < H242_ROW; q++) { std::memcpy(&w, &g.row[q], 4); rh = (rh ^ w) * 0x9E3779B97F4A7C15ull; } mx(g_h413_f[4], rh); } }
  if (g_h413dbg) { uint64_t h = g_h413_ik ^ g.key; h = h * 0x9E3779B97F4A7C15ull ^ g.meta; uint32_t u; std::memcpy(&u, &g.rfin, 4); h = h * 0x9E3779B97F4A7C15ull ^ u; std::memcpy(&u, &g.smin, 4); h = h * 0x9E3779B97F4A7C15ull ^ u ^ ((uint64_t)(uint32_t)g.sarg << 32); for (uint32_t q = 0; q < ((g.meta >> 2) & 63u) + ((g.meta >> 8) & 31u) && q < H242_ROW; q++) { std::memcpy(&u, &g.row[q], 4); h = h * 0x9E3779B97F4A7C15ull ^ u; } g_h413_ik = h; }
}
static void h242_stage(PipeStep& PS, int t, int j, const FitOut& o, uint32_t steps, PutJob242& J, int cpu = 0);
static void h242_worker_main() {
  for (;;) {
    PutJob242 J;
    {
      std::unique_lock<std::mutex> lk(G242.qm);
      G242.qcv.wait(lk, [] { return G242.stop || !G242.q.empty(); });
      if (G242.q.empty()) return;
      J = std::move(G242.q.front()); G242.q.pop_front(); G242.busy = true;
    }
    auto tw = clk::now();
    if (J.ps) { const FitOut& o = J.ps->SB->F[J.side]; for (size_t q = 0; q < J.ids.size(); q++) h242_stage(*J.ps, J.side, J.ids[q], o, J.steps[q], J, J.cpu.empty() ? 0 : (int)J.cpu[q]); }
    {
      std::lock_guard<std::mutex> tl(G242.tab_m);
      for (const Stage242& g : J.st) h242_insert(g, J.recs);
    }
    G242.worker_ms += std::chrono::duration<double, std::milli>(clk::now() - tw).count(); G242.njobs++;
    { std::lock_guard<std::mutex> lk(G242.qm); G242.busy = false; }
    G242.qdone.notify_all();
  }
}
static void h242_enqueue(PutJob242&& J) {
  if (J.st.empty() && J.ids.empty()) return;
  {
    std::lock_guard<std::mutex> lk(G242.qm);
    if (!G242.started) { G242.th = std::thread(h242_worker_main); G242.started = true; }
    G242.q.push_back(std::move(J));
  }
  G242.qcv.notify_one();
}
static void h242_drain() {
  std::unique_lock<std::mutex> lk(G242.qm);
  G242.qdone.wait(lk, [] { return G242.q.empty() && !G242.busy; });
}
static void h242_shutdown() {
  h242_drain();
  { std::lock_guard<std::mutex> lk(G242.qm); G242.stop = true; }
  G242.qcv.notify_all();
  if (G242.started) G242.th.join();
  G242.started = false;
}

static void h242_stage(PipeStep& PS, int t, int j, const FitOut& o, uint32_t steps, PutJob242& J, int cpu) {
  const Tree& tr = PS.SB->UQ.U[(size_t)j];
  const uint32_t n = (uint32_t)tr.rec.size();
  const float* th = &F14.th[t][32 * (size_t)j]; const float* ap = &F14.a[t][16 * (size_t)j];
  uint32_t lt = 0, la = 0;
  for (uint32_t q = 0; q < 32; q++) if (h242_bits(th[q])) lt = q + 1;
  const uint32_t qd = (uint32_t)std::max(0, std::min(tr.Q, (g_h770 & 2) ? 16 : 4));
  for (uint32_t q = 0; q < qd; q++) if (h242_bits(ap[q])) la = q + 1;
  const uint32_t nth = (uint32_t)std::max(0, std::min(tr.nth, 32));
  if (g_h413dbg) { if (lt > nth) g_h413_tail++; g_h413_la += la; g_h413_lt += lt; }
  if (nth > lt) lt = nth;
  if (lt + la > H242_ROW) { G242.refused_wide++; if (g_h770) g_h770s.refq[std::max(0, std::min(tr.Q, 16))]++; return; }
  if (g_h770) (cpu ? g_h770s.put_cpu : g_h770s.put_gpu)++;
  const uint32_t ig = o.is_g[j] ? 1u : 0u;
  Stage242 g;
  g.key = (g_h770 & 1) ? h770_key(PS.SB->UQ.hash[(size_t)j], (uint32_t)t, cpu) : h99_key(PS.SB->UQ.hash[(size_t)j], (uint32_t)t); g.t = (uint32_t)t;
  g.meta = 1u | (ig << 1) | (lt << 2) | (la << 8) | ((steps & 255u) << 13) | ((uint32_t)t << 21) | ((uint32_t)(cpu ? 1 : 0) << 22);
  g.rfin = o.rfin[j]; g.smin = o.smin[j]; g.sarg = o.sarg[j];
  for (uint32_t q = 0; q < lt; q++) g.row[q] = th[q];
  for (uint32_t q = 0; q < la; q++) g.row[lt + q] = ap[q];
  g.nrec_off = 0; g.nrec = 0;
  if (G242.verify) { g.nrec_off = (uint32_t)J.recs.size(); g.nrec = n; J.recs.insert(J.recs.end(), tr.rec.begin(), tr.rec.end()); }
  J.st.push_back(g);
}


static void h242_apply(PipeStep& PS, int t, int id, uint32_t e, Times& T);
static void h242_lookup(Engine& E, PipeStep& PS, const std::vector<std::pair<int, int>>& pr, std::vector<uint32_t>& fe, std::vector<uint64_t>& kv, Times* TA = nullptr) {
  auto t0 = clk::now();
  std::vector<Times> tw;
  const size_t np = pr.size();
  kv.assign(np, 0);
  for (size_t q = 0; q < np; q++) kv[q] = h99_key(PS.SB->UQ.hash[(size_t)pr[q].second], (uint32_t)pr[q].first);
  auto block = [&](size_t a, size_t b) {
    for (size_t q = a; q < b; q++) {
      if (q + 8 < b) __builtin_prefetch(&G242.tab[(size_t)kv[q + 8] & G242.mask], 0, 3);
      fe[q] = G242.probe(kv[q], (uint32_t)pr[q].first);
    }
  };
  auto blockA = [&](size_t a, size_t b, Times& Tw) {
    block(a, b);
    for (size_t q = a; q < b; q++) if (fe[q] != H99_NONE) h242_apply(PS, pr[q].first, pr[q].second, fe[q], Tw);
  };
  if (TA) {
    const size_t nw = (g_h242 == 1 && WK(E).size() > 1 && np >= 256) ? WK(E).size() : 1;
    tw.assign(nw, Times{});
    for (Times& x : tw) { x.steps_sum = 0; x.steps_n = 0; }
    if (nw > 1) PL(E).run((int)nw, [&](int w) { blockA(np * (size_t)w / nw, np * (size_t)(w + 1) / nw, tw[(size_t)w]); });
    else blockA(0, np, tw[0]);
    for (const Times& x : tw) { TA->steps_sum += x.steps_sum; TA->steps_n += x.steps_n; }
    G242.lookups += (long)np;
    G242.lookup_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
    return;
  }
  if (g_h242 == 1 && WK(E).size() > 1 && np >= 256) {
    const size_t nw = WK(E).size();
    PL(E).run((int)nw, [&](int w) { block(np * (size_t)w / nw, np * (size_t)(w + 1) / nw); });
  } else block(0, np);
  G242.lookups += (long)np;
  G242.lookup_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

static void h242_apply(PipeStep& PS, int t, int id, uint32_t e, Times& T) {
  FitOut& o = PS.SB->F[t];
  const Slot242& s = G242.tab[e];
  const uint32_t lt = (s.meta >> 2) & 63u, la = (s.meta >> 8) & 31u, ig = (s.meta >> 1) & 1u, steps = (s.meta >> 13) & 255u;
  o.rfin[id] = s.rfin; o.smin[id] = s.smin; o.sarg[id] = s.sarg; o.is_g[id] = (char)ig;
  float* th = &F14.th[t][32 * (size_t)id];
  for (uint32_t q = 0; q < 32; q++) th[q] = q < lt ? s.row[q] : 0.0f;
  float* ap = &F14.a[t][16 * (size_t)id];
  for (uint32_t q = 0; q < 16; q++) ap[q] = q < la ? s.row[lt + q] : 0.0f;
  if (ig) { T.steps_sum += (double)steps; T.steps_n += 1; }
}

static void h770_apply_raw(PipeStep& PS, int t, int id, const unsigned char* raw, Times& T) {
  Slot242 s; std::memcpy(&s, raw, sizeof s);
  FitOut& o = PS.SB->F[t];
  const uint32_t lt = (s.meta >> 2) & 63u, la = (s.meta >> 8) & 31u, ig = (s.meta >> 1) & 1u, steps = (s.meta >> 13) & 255u;
  o.rfin[id] = s.rfin; o.smin[id] = s.smin; o.sarg[id] = s.sarg; o.is_g[id] = (char)ig;
  float* th = &F14.th[t][32 * (size_t)id];
  for (uint32_t q = 0; q < 32; q++) th[q] = q < lt ? s.row[q] : 0.0f;
  float* ap = &F14.a[t][16 * (size_t)id];
  for (uint32_t q = 0; q < 16; q++) ap[q] = q < la ? s.row[lt + q] : 0.0f;
  if (ig) { T.steps_sum += (double)steps; T.steps_n += 1; }
}


static void h770_schedule(Engine& E, PipeStep& PS, Batch& B, Times& T, const std::vector<std::pair<int, int>>& pr) {
  const size_t np = pr.size();
  h242_drain();
  std::lock_guard<std::mutex> tl(G242.tab_m);
  auto t0 = clk::now();
  auto Hp = std::make_shared<H770B>(); B.h770 = Hp; H770B& H = *Hp;
  std::vector<uint32_t> e0(np, H99_NONE), e1(np, H99_NONE);
  std::vector<uint64_t> k0(np), k1(np);
  const Unique& UQ = PS.SB->UQ;
  for (size_t q = 0; q < np; q++) { const uint64_t hr = UQ.hash[(size_t)pr[q].second]; k0[q] = h770_key(hr, (uint32_t)pr[q].first, 0); k1[q] = h770_key(hr, (uint32_t)pr[q].first, 1); }
  auto block = [&](size_t a, size_t b) {
    for (size_t q = a; q < b; q++) {
      if (q + 8 < b) { __builtin_prefetch(&G242.tab[(size_t)k0[q + 8] & G242.mask], 0, 3); __builtin_prefetch(&G242.tab[(size_t)k1[q + 8] & G242.mask], 0, 3); }
      e0[q] = G242.probe(k0[q], (uint32_t)pr[q].first); e1[q] = G242.probe(k1[q], (uint32_t)pr[q].first);
    }
  };
  if (g_h242 == 1 && WK(E).size() > 1 && np >= 256) {
    const size_t nw = WK(E).size();
    PL(E).run((int)nw, [&](int w) { block(np * (size_t)w / nw, np * (size_t)(w + 1) / nw); });
  } else block(0, np);
  size_t nsch[2] = {0, 0};
  for (size_t q = 0; q < np; q++) {
    const int t = pr[q].first, id = pr[q].second;
    nsch[t]++;
    if (e0[q] == H99_NONE && e1[q] == H99_NONE) { B.need[t].push_back(id); continue; }
    G242.hits++; g_h770s.hits++;
    H770Hit h; std::memset(&h, 0, sizeof h); h.t = t; h.id = id; h.vs = -1; h.st = 0;
    h.v = (uint8_t)((e0[q] != H99_NONE ? 1 : 0) | (e1[q] != H99_NONE ? 2 : 0));
    if (e0[q] != H99_NONE) std::memcpy(h.raw[0], &G242.tab[e0[q]], 128);
    if (e1[q] != H99_NONE) std::memcpy(h.raw[1], &G242.tab[e1[q]], 128);
    const Slot242& a = G242.tab[e0[q] != H99_NONE ? e0[q] : e1[q]];
    if (e0[q] != H99_NONE && e1[q] != H99_NONE) {
      const Slot242& b = G242.tab[e1[q]];
      if (h242_bits(a.smin) != h242_bits(b.smin) || a.sarg != b.sarg) g_h770s.rowdis++;
    }
    if (g_h770v) B.need[t].push_back(id);
    else { FitOut& o = PS.SB->F[t]; o.smin[id] = a.smin; o.sarg[id] = a.sarg; }
    if (H.hix[t].empty()) H.hix[t].assign(UQ.U.size(), -1);
    if ((size_t)id < H.hix[t].size()) H.hix[t][(size_t)id] = (int32_t)H.hits.size();
    H.hits.push_back(h);
  }
  g_h770s.pairs += (long)np;
  G242.lookups += (long)np;
  G242.lookup_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
  T.nfit0 += (double)nsch[0]; T.nfit1 += (double)nsch[1];
}

static void h770_verify(PipeStep& PS, Batch& B, int s, const uint32_t* ST) {
  FitOut& o = PS.SB->F[s]; H770B& H = *B.h770;
  for (H770Hit& h : H.hits) {
    if (h.t != s) continue;
    const int j = h.id; const Tree& tr = PS.SB->UQ.U[(size_t)j];
    const int c = (h.vs >= 0) ? h.vs : ((h.v & 1) ? 0 : 1);
    if (!((h.v >> c) & 1)) { g_h770s.vnos++; continue; }
    Slot242 sl; std::memcpy(&sl, h.raw[c], sizeof sl);
    const uint32_t lt = (sl.meta >> 2) & 63u, la = (sl.meta >> 8) & 31u, ig = (sl.meta >> 1) & 1u, stp = (sl.meta >> 13) & 255u;
    const float* fth = &F14.th[s][32 * (size_t)j]; const float* fa = &F14.a[s][16 * (size_t)j];
    const uint32_t nth = (uint32_t)std::max(0, std::min(tr.nth, 32));
    const uint32_t wdef = (uint32_t)std::max(0, std::min(tr.Q, (g_h770 & 2) ? 16 : 4));
    long d[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    d[0] = h242_bits(sl.rfin) != h242_bits(o.rfin[j]); d[1] = h242_bits(sl.smin) != h242_bits(o.smin[j]); d[2] = sl.sarg != o.sarg[j];
    d[3] = (int)ig != (o.is_g[j] ? 1 : 0); d[4] = stp != (o.is_g[j] ? (ST[j] & 255u) : 0u);
    for (uint32_t q = 0; q < 32; q++) { const uint32_t sv = h242_bits(q < lt ? sl.row[q] : 0.0f); if (sv != h242_bits(fth[q])) d[q < nth ? 5 : 6]++; }
    for (uint32_t q = 0; q < 16; q++) { const uint32_t sv = h242_bits(q < la ? sl.row[lt + q] : 0.0f); if (sv != h242_bits(fa[q])) d[q < wdef ? 7 : 8]++; }
    const long bad = d[0] + d[1] + d[2] + d[3] + d[4] + d[5] + d[7], tail = d[6] + d[8];
    g_h770s.ver++; if (bad) g_h770s.vdiff++; else g_h770s.veq++; if (tail) g_h770s.vtail++;
    for (int k = 0; k < 9; k++) g_h770s.dw[k] += d[k];
  }
}
static long g_h770_cum[16];
static void h770_step_print(int st) {
  H770Stat& S = g_h770s;
  auto X = [](std::atomic<long>& a) { return a.exchange(0); };
  const long pairs = X(S.pairs), hits = X(S.hits), full = X(S.full), part = X(S.partial), dg = X(S.dropg), dc = X(S.dropc), nb = X(S.nonbin), rd = X(S.rowdis), pc = X(S.put_cpu), pg = X(S.put_gpu);
  printf("H770 step %2d: pairs %ld hits %ld (%.4f) | full hits applied %ld (GPU-part slot skips %ld, CPU-suffix drops %ld, no GN bin %ld) | partial (no row for the virtual side, refit) %ld | pairs whose two side rows disagree on smin/sarg %ld | puts cpu-rows %ld gpu-rows %ld | memo clears so far %ld\n",
         st, pairs, hits, pairs ? (double)hits / (double)pairs : 0.0, full, dg, dc, nb, part, rd, pc, pg, G242.clears);
  long rq[17], rt = 0; for (int q = 0; q <= 16; q++) { rq[q] = S.refq[q].exchange(0); rt += rq[q]; }
  printf("H770 step %2d: refused (row wider than the slot) %ld by Q:", st, rt); for (int q = 0; q <= 16; q++) if (rq[q]) printf(" Q%d=%ld", q, rq[q]); printf("\n");
  const long c0[10] = {pairs, hits, full, part, dg, dc, nb, rd, rt, 0}; for (int k = 0; k < 9; k++) g_h770_cum[k] += c0[k];
  if (g_h770v) {
    const long ver = X(S.ver), veq = X(S.veq), vd = X(S.vdiff), vt = X(S.vtail), vn = X(S.vnos);
    long dw[9]; for (int k = 0; k < 9; k++) dw[k] = S.dw[k].exchange(0);
    printf("H770V step %2d: hits compared %ld | rows equal on every read word %ld | rows differing %ld | rows with differing unread-tail words %ld | no row for the fit side %ld\n", st, ver, veq, vd, vt, vn);
    printf("H770V step %2d: differing words by kind: rfin %ld smin %ld sarg %ld is_g %ld steps %ld | th[0..nth) %ld th[nth..32) %ld | a[0..wdef) %ld a[wdef..16) %ld   (wdef = Q with bit 1, else min(Q,4))\n",
           st, dw[0], dw[1], dw[2], dw[3], dw[4], dw[5], dw[6], dw[7], dw[8]);
    g_h770_cum[9] += ver; g_h770_cum[10] += veq; g_h770_cum[11] += vd; g_h770_cum[12] += vt; g_h770_cum[13] += vn;
  }
}
static void h770_total_print() {
  printf("H770 TOTAL: pairs %ld hits %ld full %ld partial %ld GPU-slot skips %ld CPU drops %ld no-bin %ld side-rows-disagree %ld refused %ld | VERIFY compared %ld equal %ld differing %ld tail-differing %ld noside %ld\n",
         g_h770_cum[0], g_h770_cum[1], g_h770_cum[2], g_h770_cum[3], g_h770_cum[4], g_h770_cum[5], g_h770_cum[6], g_h770_cum[7], g_h770_cum[8], g_h770_cum[9], g_h770_cum[10], g_h770_cum[11], g_h770_cum[12], g_h770_cum[13]);
}


static constexpr int H771_NSUB = 16;
struct H771State {
  Slot242* tab[4][2] = {{nullptr, nullptr}, {nullptr, nullptr}, {nullptr, nullptr}, {nullptr, nullptr}};
  size_t subcap = 0, submask = 0, submax = 0; int ng = 2; uint64_t gen = 0;
  uint32_t nlive[2][H771_NSUB];
  std::vector<uint64_t> key[2];
  std::vector<uint64_t> kv; std::vector<const Slot242*> res;
  std::atomic<long> lookups{0}, found{0}, put_fit{0}, put_carry{0}, ref_cap{0}, ref_wide{0}, dup{0};
  double lookup_ms = 0, put_ms = 0, put_wall_n = 0;
};
static H771State G771;
static inline uint32_t h771_stamp(uint64_t g) { return 1u + (uint32_t)(g % 511u); }
static void h771_init(size_t np, int ng, int log2cap) {
  G771.ng = std::max(2, std::min(4, ng));
  size_t cap = (size_t)1 << log2cap;
  G771.subcap = cap / H771_NSUB; G771.submask = G771.subcap - 1; G771.submax = G771.subcap / 4 * 3;
  for (int g = 0; g < G771.ng; g++) for (int t = 0; t < 2; t++) {
    void* p = nullptr;
    if (posix_memalign(&p, 128, cap * sizeof(Slot242)) != 0 || !p) { fprintf(stderr, "H771: table allocation failed\n"); std::exit(1); }
    std::memset(p, 0, cap * sizeof(Slot242)); G771.tab[g][t] = (Slot242*)p;
  }
  (void)np;
}

static void h771_begin(int n) {
  G771.gen++;
  if (G771.gen % 511u == 0) { const size_t bytes = G771.subcap * H771_NSUB * sizeof(Slot242); for (int g = 0; g < G771.ng; g++) for (int t = 0; t < 2; t++) std::memset((void*)G771.tab[g][t], 0, bytes); }
  for (int t = 0; t < 2; t++) for (int u = 0; u < H771_NSUB; u++) G771.nlive[t][u] = 0;
  for (int t = 0; t < 2; t++) if (G771.key[t].size() < (size_t)2 * n + 2) G771.key[t].resize((size_t)2 * n + 2);
}
static inline const Slot242* h771_base(int g, int t, uint64_t key) { return G771.tab[g][t] + (size_t)(key >> 60) * G771.subcap; }
static inline bool h771_live(const Slot242& s, uint32_t stamp) { return (s.meta & 1u) && (s.meta >> 23) == stamp; }

static const Slot242* h771_probe(int g, int t, uint64_t key, uint32_t stamp) {
  const Slot242* b = h771_base(g, t, key);
  size_t i = (size_t)key & G771.submask;
  for (size_t n = 0; n < G771.subcap; n++, i = (i + 1) & G771.submask) {
    const Slot242& s = b[i];
    if (!h771_live(s, stamp)) return nullptr;
    if (s.key == key) return &s;
  }
  return nullptr;
}

static int h771_insert(int g, int t, const Slot242& src, uint64_t key, uint32_t stamp) {
  const int u = (int)(key >> 60);
  Slot242* b = G771.tab[g][t] + (size_t)u * G771.subcap;
  size_t i = (size_t)key & G771.submask;
  for (size_t n = 0; n < G771.subcap; n++, i = (i + 1) & G771.submask) {
    Slot242& s = b[i];
    if (h771_live(s, stamp)) { if (s.key == key) return 0; continue; }
    if (G771.nlive[t][u] >= G771.submax) return -1;
    std::memcpy((void*)&s, &src, sizeof s);
    s.key = key; s.meta = (src.meta & 0x007FFFFFu) | 1u | (stamp << 23);
    G771.nlive[t][u]++;
    return 1;
  }
  return -1;
}

static void h771_schedule(Engine& E, PipeStep& PS, Batch& B, Times& T, const std::vector<std::pair<int, int>>& pr) {
  const size_t np = pr.size();
  auto t0 = clk::now();
  auto Hp = std::make_shared<H770B>(); B.h770 = Hp; H770B& H = *Hp;
  const Unique& UQ = PS.SB->UQ;
  const int npg = G771.ng - 1;
  uint32_t stp[3] = {0, 0, 0}; int tg[3] = {0, 0, 0}; bool ok[3] = {false, false, false};
  for (int k = 1; k <= npg; k++) { ok[k - 1] = G771.gen > (uint64_t)k; stp[k - 1] = h771_stamp(G771.gen - (uint64_t)k); tg[k - 1] = (int)((G771.gen - (uint64_t)k) % (uint64_t)G771.ng); }
  if (G771.kv.size() < np) G771.kv.resize(np);
  if (G771.res.size() < np * 3) G771.res.resize(np * 3);
  uint64_t* kv = G771.kv.data(); const Slot242** res = G771.res.data();
  auto block = [&](size_t a, size_t b) {
    for (size_t q = a; q < b; q++) { const int t = pr[q].first, id = pr[q].second; const uint64_t kk = h99_key(UQ.hash[(size_t)id], (uint32_t)t); kv[q] = kk; G771.key[t][(size_t)id] = kk; }
    for (size_t q = a; q < b; q++) {
      if (q + 8 < b) for (int k = 0; k < npg; k++) if (ok[k]) __builtin_prefetch(h771_base(tg[k], pr[q + 8].first, kv[q + 8]) + ((size_t)kv[q + 8] & G771.submask), 0, 3);
      for (int k = 0; k < 3; k++) res[q * 3 + (size_t)k] = (k < npg && ok[k]) ? h771_probe(tg[k], pr[q].first, kv[q], stp[k]) : nullptr;
    }
  };
  if (WK(E).size() > 1 && np >= 256) {
    const size_t nw = WK(E).size();
    PL(E).run((int)nw, [&](int w) { block(np * (size_t)w / nw, np * (size_t)(w + 1) / nw); });
  } else block(0, np);
  size_t nsch[2] = {0, 0}; long nhit = 0;
  for (size_t q = 0; q < np; q++) {
    const int t = pr[q].first, id = pr[q].second;
    nsch[t]++;
    const Slot242* rs[2] = {nullptr, nullptr};
    for (int k = 0; k < npg; k++) { const Slot242* s = res[q * 3 + (size_t)k]; if (!s) continue; const int c = (int)((s->meta >> 22) & 1u); if (!rs[c]) rs[c] = s; }
    if (!rs[0] && !rs[1]) { B.need[t].push_back(id); continue; }
    nhit++; g_h770s.hits++;
    H770Hit h; std::memset(&h, 0, sizeof h); h.t = t; h.id = id; h.vs = -1; h.st = 0;
    h.v = (uint8_t)((rs[0] ? 1 : 0) | (rs[1] ? 2 : 0));
    if (rs[0]) std::memcpy(h.raw[0], rs[0], 128);
    if (rs[1]) std::memcpy(h.raw[1], rs[1], 128);
    const Slot242& a = rs[0] ? *rs[0] : *rs[1];
    if (rs[0] && rs[1] && (h242_bits(rs[0]->smin) != h242_bits(rs[1]->smin) || rs[0]->sarg != rs[1]->sarg)) g_h770s.rowdis++;
    if (g_h770v) B.need[t].push_back(id);
    else { FitOut& o = PS.SB->F[t]; o.smin[id] = a.smin; o.sarg[id] = a.sarg; }
    if (H.hix[t].empty()) H.hix[t].assign(UQ.U.size(), -1);
    if ((size_t)id < H.hix[t].size()) H.hix[t][(size_t)id] = (int32_t)H.hits.size();
    H.hits.push_back(h);
  }
  g_h770s.pairs += (long)np; G771.lookups += (long)np; G771.found += nhit;
  G771.lookup_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
  T.nfit0 += (double)nsch[0]; T.nfit1 += (double)nsch[1];
}

static bool h771_fit_row(PipeStep& PS, int t, int j, const FitOut& o, uint32_t steps, int cpu, Slot242& g) {
  const Tree& tr = PS.SB->UQ.U[(size_t)j];
  const float* th = &F14.th[t][32 * (size_t)j]; const float* ap = &F14.a[t][16 * (size_t)j];
  uint32_t lt = 0, la = 0;
  for (uint32_t q = 0; q < 32; q++) if (h242_bits(th[q])) lt = q + 1;
  const uint32_t qd = (uint32_t)std::max(0, std::min(tr.Q, 16));
  for (uint32_t q = 0; q < qd; q++) if (h242_bits(ap[q])) la = q + 1;
  const uint32_t nth = (uint32_t)std::max(0, std::min(tr.nth, 32));
  if (nth > lt) lt = nth;
  if (lt + la > H242_ROW) { g_h770s.refq[std::max(0, std::min(tr.Q, 16))]++; G771.ref_wide++; return false; }
  const uint32_t ig = o.is_g[j] ? 1u : 0u;
  g.key = 0; g.meta = 1u | (ig << 1) | (lt << 2) | (la << 8) | ((steps & 255u) << 13) | ((uint32_t)t << 21) | ((uint32_t)(cpu ? 1 : 0) << 22);
  g.rfin = o.rfin[j]; g.smin = o.smin[j]; g.sarg = o.sarg[j];
  for (uint32_t q = 0; q < lt; q++) g.row[q] = th[q];
  for (uint32_t q = 0; q < la; q++) g.row[lt + q] = ap[q];
  return true;
}


static void h771_put(Engine& E, PipeStep& PS, Batch& B) {
  auto t0 = clk::now();
  Persist5& P = *B.P;
  const uint32_t* ST[2] = {(const uint32_t*)P.S[0].gst.b.contents, (const uint32_t*)P.S[1].gst.b.contents};
  const H770B* H = (!g_h770v && B.h770) ? B.h770.get() : nullptr;
  const int cur = (int)(G771.gen % (uint64_t)G771.ng); const uint32_t stamp = h771_stamp(G771.gen);
  size_t items = B.need[0].size() + B.need[1].size() + (H ? H->hits.size() : 0);
  auto work = [&](int w, int nw) {
    long nf = 0, nc = 0, rc = 0, dp = 0;
    for (int t = 0; t < 2; t++) {
      const FitOut& o = PS.SB->F[t];
      for (int j : B.need[t]) {
        const uint64_t key = G771.key[t][(size_t)j];
        if ((t * H771_NSUB + (int)(key >> 60)) % nw != w) continue;
        const int cpu = ((size_t)j < B.h770cpu[t].size()) ? (int)B.h770cpu[t][(size_t)j] : 0;
        Slot242 g;
        if (!h771_fit_row(PS, t, j, o, o.is_g[j] ? (ST[t][j] & 255u) : 0u, cpu, g)) continue;
        const int r = h771_insert(cur, t, g, key, stamp);
        if (r > 0) nf++; else if (r == 0) dp++; else rc++;
      }
    }
    if (H) for (const H770Hit& h : H->hits) {
      if (h.st != 1) continue;
      const uint64_t key = G771.key[h.t][(size_t)h.id];
      if ((h.t * H771_NSUB + (int)(key >> 60)) % nw != w) continue;
      Slot242 g; std::memcpy(&g, h.raw[h.vs >= 0 ? h.vs : ((h.v & 1) ? 0 : 1)], sizeof g);
      const int r = h771_insert(cur, h.t, g, key, stamp);
      if (r > 0) nc++; else if (r == 0) dp++; else rc++;
    }
    G771.put_fit += nf; G771.put_carry += nc; G771.ref_cap += rc; G771.dup += dp;
  };
  if (WK(E).size() > 1 && items >= 256) { const int nw = (int)WK(E).size(); PL(E).run(nw, [&](int w) { work(w, nw); }); }
  else work(0, 1);
  G771.put_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}
static void h771_step_print(int st) {
  printf("H771 step %2d: lookups %ld found %ld | puts: fits %ld carried hits %ld | refused: sub-table full %ld, row wider than the slot %ld | duplicate keys %ld | lookup %.3f ms (producer wall) put %.3f ms (main wall)\n",
         st, G771.lookups.exchange(0), G771.found.exchange(0), G771.put_fit.exchange(0), G771.put_carry.exchange(0), G771.ref_cap.exchange(0), G771.ref_wide.exchange(0), G771.dup.exchange(0), G771.lookup_ms, G771.put_ms);
  G771.lookup_ms = G771.put_ms = 0;
}


static void h242_schedule(Engine& E, PipeStep& PS, Batch& B, Times& T, const std::vector<std::pair<int, int>>& pr) {
  if (g_h770) { h770_schedule(E, PS, B, T, pr); return; }
  const size_t np = pr.size();
  std::vector<uint32_t> fe(np, H99_NONE); std::vector<uint64_t> kv;
  h242_drain();
  std::lock_guard<std::mutex> tl(G242.tab_m);
  const bool h413a = g_h413 && !G242.verify;
  if (np) h242_lookup(E, PS, pr, fe, kv, h413a ? &T : nullptr);
  size_t nsch[2] = {0, 0};
  auto ta = clk::now();
  if (h413a) {
    for (size_t q = 0; q < np; q++) { nsch[pr[q].first]++; if (fe[q] != H99_NONE) { G242.hits++; if (g_h413dbg) g_h413_hk += kv[q] * 0x9E3779B97F4A7C15ull; } else B.need[pr[q].first].push_back(pr[q].second); }
    G242.apply_ms += std::chrono::duration<double, std::milli>(clk::now() - ta).count();
    T.nfit0 += (double)nsch[0]; T.nfit1 += (double)nsch[1];
    return;
  }
  for (size_t q = 0; q < np; q++) {
    const int t = pr[q].first, id = pr[q].second;
    nsch[t]++;
    uint32_t e = fe[q];
    if (e != H99_NONE && G242.verify) {
      const Tree& tr = PS.SB->UQ.U[(size_t)id];
      if (!h242_recmatch(e, tr.rec.data(), (uint32_t)tr.rec.size())) { G242.vmis++; e = H99_NONE; }
    }
    if (e != H99_NONE) { G242.hits++; if (g_h413dbg) g_h413_hk += kv[q] * 0x9E3779B97F4A7C15ull; h242_apply(PS, t, id, e, T); }
    else B.need[t].push_back(id);
  }
  G242.apply_ms += std::chrono::duration<double, std::milli>(clk::now() - ta).count();
  T.nfit0 += (double)nsch[0]; T.nfit1 += (double)nsch[1];
}

static void gn_finish(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  auto t0 = clk::now();
  g_h882q->wait_for(B.t882G);
  if (B.cbG) {
    k_wait(B.cbG);
    if (B.cbG.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "GN command buffer failed\n"); h182_dump(B.cbG, "GN"); gpu_fail(); }
    B.gn_gpu = (B.cbG.GPUEndTime - B.cbG.GPUStartTime) * 1e3;
    PS.iv.push_back({B.cbG.GPUStartTime, B.cbG.GPUEndTime});
  }
  { const double w = ms_since(t0); T.w_g += w; T.gn_wall += w; }
  if (B.h250_on) {
    h250a_join(B);
    const size_t nc = B.h250_items.size(), ng = B.h250_ngpu;
    if (B.gn_gpu > 0 && ng > 0) g_h250a_rg = g_h250a_rg > 0 ? 0.5 * g_h250a_rg + 0.5 * ((double)ng / B.gn_gpu) : (double)ng / B.gn_gpu;
    if (B.h250_cpu_ms > 0 && nc > 0) g_h250a_rc = g_h250a_rc > 0 ? 0.5 * g_h250a_rc + 0.5 * ((double)nc / B.h250_cpu_ms) : (double)nc / B.h250_cpu_ms;
    { const int pc = std::min(B.ci, 15); double& rg = g_h300b_rg[pc]; double& rc = g_h300b_rc[pc];
      if (B.gn_gpu > 0 && ng > 0) rg = rg > 0 ? 0.5 * rg + 0.5 * ((double)ng / B.gn_gpu) : (double)ng / B.gn_gpu;
      if (B.h250_cpu_ms > 0 && nc > 0) rc = rc > 0 ? 0.5 * rc + 0.5 * ((double)nc / B.h250_cpu_ms) : (double)nc / B.h250_cpu_ms; }
    printf("H250A: f_cpu %.3f cpu ms %.3f gpu ms %.3f trees cpu %zu gpu %zu\n", B.h250_f, B.h250_cpu_ms, B.gn_gpu, nc, ng);
    B.h250_on = false;
  }

  T.gn_gpu += B.gn_gpu;
  Persist5& P = *B.P;
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    const float* RS = (const float*)X.grs.b.contents; const uint32_t* STl = (const uint32_t*)X.gst.b.contents;
    const float* TH8 = (const float*)X.gth.b.contents; const float* A8 = (const float*)X.gai.b.contents;
    static thread_local std::vector<uint32_t> STg916;
    if (g_h916 && STg916.size() < PS.SB->UQ.U.size()) STg916.resize(PS.SB->UQ.U.size());
    for (int j : B.need[s]) if (o.is_g[j]) {
      const size_t lj = glid916(X, (size_t)j);
      o.rfin[j] = RS[lj * 2 + 1];
      if (!isfinite(o.rfin[j])) o.rfin[j] = INFINITY;
      std::memcpy(&F14.th[s][32 * (size_t)j], TH8 + 32 * lj, 32 * sizeof(float)); std::memcpy(&F14.a[s][16 * (size_t)j], A8 + 16 * lj, 16 * sizeof(float));
      T.steps_sum += STl[lj] & 255; T.steps_n += 1; if (g_h916) STg916[(size_t)j] = STl[lj]; }
    const uint32_t* ST = g_h916 ? STg916.data() : STl;
    if (g_h242) {
      auto tp = clk::now();
      PutJob242 J;
      if (g_h770) { if (g_h770v && B.h770) h770_verify(PS, B, s, ST); }
      if (g_h413) { J.ps = &PS; J.side = s; J.ids = B.need[s]; J.steps.resize(J.ids.size()); if (g_h770) { J.cpu.resize(J.ids.size()); for (size_t q = 0; q < J.ids.size(); q++) J.cpu[q] = ((size_t)J.ids[q] < B.h770cpu[s].size()) ? (uint8_t)B.h770cpu[s][(size_t)J.ids[q]] : (uint8_t)0; } for (size_t q = 0; q < J.ids.size(); q++) { const int j = J.ids[q]; J.steps[q] = o.is_g[j] ? (ST[j] & 255u) : 0u; } }
      else for (int j : B.need[s]) h242_stage(PS, s, j, o, o.is_g[j] ? (ST[j] & 255u) : 0u, J, (g_h770 && (size_t)j < B.h770cpu[s].size()) ? (int)B.h770cpu[s][(size_t)j] : 0);
      G242.nputs += (long)B.need[s].size();
      if (G242.sync) {
        std::lock_guard<std::mutex> tl(G242.tab_m);
        for (const Stage242& g : J.st) h242_insert(g, J.recs);
      } else h242_enqueue(std::move(J));
      G242.put_ms += std::chrono::duration<double, std::milli>(clk::now() - tp).count();
    } else if (g_h771) { if (g_h770v && B.h770) h770_verify(PS, B, s, ST); }
    else if (g_h99) for (int j : B.need[s]) h99_put(PS, s, j, o, o.is_g[j] ? (ST[j] & 255u) : 0u);
  }
  if (g_h771) h771_put(E, PS, B);
  B.gdone = true;
}


static double (*g_fr_score96)(const float* r96) = nullptr;


static const bool g_nres = getenv("NORMAL_RES") && atoi(getenv("NORMAL_RES")) != 0;
static float g_nres_n[96];


static const int g_yfx_flag = getenv("YFX_CELLS") ? atoi(getenv("YFX_CELLS")) : 1;
static std::vector<double> g_yfx_txy[2];
static uint64_t g_yfx_epoch = 0; static bool g_yfx = false; static double g_yfx_t[48], g_yfx_y[2][48], g_yfx_n[96]; static float g_yfx_f[2][48], g_yfx_nt[96];
static void nres_init(const double* X, const double* Y) {
  double tx[NPT], ty[NPT];
  for (int i = 0; i < NPT; i++) { const int a = i > 0 ? i - 1 : i, b = i < NPT - 1 ? i + 1 : i; tx[i] = X[b] - X[a]; ty[i] = Y[b] - Y[a]; }
  double px = 0, py = 1;
  for (int i = 0; i < NPT; i++) {
    const double h = std::hypot(tx[i], ty[i]);
    if (h > 0) { px = -ty[i] / h; py = tx[i] / h; }
    g_nres_n[i] = (float)px; g_nres_n[NPT + i] = (float)py;
  }
}
static void glue_range(PipeStep& PS, int b0, int b1, float* err, float* rx, float* ry, Times& T) {
  auto t0 = clk::now();
  const FitOut& F0 = PS.SB->F[0]; const FitOut& F1 = PS.SB->F[1];
  const std::vector<int>& ct = PS.SB->ct; const std::vector<uint8_t>& admit = *PS.admit;
  for (int i = b0; i < b1; i++) {
    if (!admit[i]) { err[i] = DEAD; rx[i] = INFINITY; ry[i] = INFINITY; continue; }
    int a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1];
    float fa = (a >= 0) ? F0.rfin[a] : INFINITY, fb = (b >= 0) ? F1.rfin[b] : INFINITY;
    rx[i] = fa; ry[i] = fb;
    if (!isfinite(fa) || !isfinite(fb) || g_fr_score96) { err[i] = DEAD; continue; }
    double s = (double)fa + (double)fb; err[i] = (float)std::min(s, (double)DEAD);
  }
  T.glue_e += ms_since(t0);
}


#include "freet_h950.inc"
static const int g_h837 = getenv("H837") ? atoi(getenv("H837")) : 0;
static int g_stage8_only = 0, g_src14 = 0, g_chk_band = 0, g_post_call = 0, g_time14 = 0, g_chk_fin = 0;
static std::string g_dump14;
static FILE* g_dumpf = nullptr;


static bool zgconv14(const Tree& T, std::vector<uint32_t>& p4, uint32_t* segs, uint32_t& hw) {
  const uint32_t NSC = (uint32_t)T.nS + 2, SS = 48u, SB = (uint32_t)T.nslot * SS;
  auto enc = [&](uint32_t c) -> uint32_t { if (c == 1022) return 1u << 31; if (c < 64) return c * SS; if (c == 900) return 1u << 29;
    if (c == 901) return (1u << 30) | (SB + (uint32_t)T.nS); return (1u << 30) | (SB + (c - 64)); };
  const uint32_t TWO = (1u << 30) | (SB + (uint32_t)T.nS + 1);
  p4.clear();
  const uint32_t* prg = h313a::prog(T); const size_t nprg = h313a::nprog(T);
  for (size_t pi = 0; pi < nprg; pi++) { const uint32_t w = prg[pi]; uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22, zo, A, B;
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a); B = enc(901); } else if (op == 3) { zo = 0; A = enc(a); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b); B = enc(a); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a); B = enc(b); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a); B = enc(901); }
    else { zo = 6; A = enc(a); B = enc(901); }
    p4.push_back(zo); p4.push_back(d == 63 ? 0xFFFFFFFFu : d * SS); p4.push_back(A); p4.push_back(B); }
  const uint16_t* tsg = h313a::segs(T);
  for (int q = 0; q < 18; q++) segs[q] = tsg[q];
  for (int q = 0; q < 18; q++) segs[18 + q] = enc(tsg[18 + q]);
  if (SB >= 65536) return false;
  hw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
  return true;
}

static bool zgconv14d(const Tree& T, std::vector<uint32_t>& p4, uint32_t* segs, uint32_t& hw) {
  const uint32_t NSC = (uint32_t)T.nS + 2, SS = 48u, SB = (uint32_t)T.nslot * SS;
  auto enc = [&](uint32_t c) -> uint32_t { if (c == 1022) return 1u << 31; if (c < 64) return c * SS; if (c == 900) return 1u << 29;
    if (c == 901) return (1u << 30) | (SB + (uint32_t)T.nS); return (1u << 30) | (SB + (c - 64)); };
  const uint32_t TWO = (1u << 30) | (SB + (uint32_t)T.nS + 1);
  const uint32_t* prg = h313a::prog(T); const size_t nprg = h313a::nprog(T);
  for (size_t pi = 0; pi < nprg; pi++) { const uint32_t w = prg[pi]; uint32_t op = w & 63, d = (w >> 6) & 63, a = (w >> 12) & 1023, b = w >> 22, zo, A, B;
    if (op >= 3 && op <= 6) {
      if (a == 1022 && b == 1022) { if (op == 5) { zo = 11; A = enc(a); B = enc(901); } else if (op == 3) { zo = 0; A = enc(a); B = TWO; } else return false; }
      else if (b == 1022) { A = enc(b); B = enc(a); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 3 : 5; }
      else { A = enc(a); B = enc(b); zo = op == 5 ? 0 : op == 3 ? 1 : op == 4 ? 2 : 4; } }
    else if (op >= 7 && op <= 15) { zo = op; A = enc(a); B = enc(901); }
    else { zo = 6; A = enc(a); B = enc(901); }
    p4.push_back(zo); p4.push_back(d == 63 ? 0xFFFFFFFFu : d * SS); p4.push_back(A); p4.push_back(B); }
  const uint16_t* tsg = h313a::segs(T);
  for (int q = 0; q < 18; q++) segs[q] = tsg[q];
  for (int q = 0; q < 18; q++) segs[18 + q] = enc(tsg[18 + q]);
  if (SB >= 65536) return false;
  hw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
  return true;
}


static inline double segd14(const float* aos, const float S[band::NROW][2], int i, int s, double& u) {
  double Ax = aos[2 * s], Ay = aos[2 * s + 1];
  double Ex = (double)aos[2 * (s + 1)] - Ax, Ey = (double)aos[2 * (s + 1) + 1] - Ay;
  double L2 = std::max(Ex * Ex + Ey * Ey, 1e-30);
  double uu = ((S[i][0] - Ax) * Ex + (S[i][1] - Ay) * Ey) / L2;
  uu = std::min(1.0, std::max(0.0, uu)); u = uu;
  double qx = S[i][0] - (Ax + uu * Ex), qy = S[i][1] - (Ay + uu * Ey);
  return qx * qx + qy * qy;
}
static void ref_cand14(const float* aos, const float S[band::NROW][2], int lev[band::NROW], double& obj) {
  static double D[band::NROW][band::NSEG]; static int arg[band::NROW][band::NSEG];
  const double INF = 1e300;
  for (int i = 0; i < band::NROW; i++) for (int s = 0; s < band::NSEG; s++) {
    if (!band::corridor_ok(s, i)) { D[i][s] = INF; arg[i][s] = -1; continue; }
    double u; double c = segd14(aos, S, i, s, u);
    if (i == 0) { D[0][s] = c; arg[0][s] = -1; continue; }
    double best = INF; int bi = -1;
    for (int s2 = 0; s2 <= s; s2++) if (D[i - 1][s2] < best) { best = D[i - 1][s2]; bi = s2; }
    if (bi < 0) { D[i][s] = INF; arg[i][s] = -1; } else { D[i][s] = c + best; arg[i][s] = bi; }
  }
  double best = INF; int s47 = -1;
  for (int s = 0; s < band::NSEG; s++) if (D[band::NROW - 1][s] < best) { best = D[band::NROW - 1][s]; s47 = s; }
  lev[band::NROW - 1] = s47;
  for (int i = band::NROW - 1; i >= 1; i--) lev[i - 1] = arg[i][lev[i]];
  obj = best;
}
static double path_cost14(const float* aos, const float S[band::NROW][2], const int* lev) {
  double tot = 0; for (int i = 0; i < band::NROW; i++) { double u; tot += segd14(aos, S, i, lev[i], u); } return tot;
}


struct P14Bin { int side = 0, Q = 0, nth = 0, tpg = 1; uint32_t per = 0, wix = 0; size_t nt = 0, ntg = 0, tgb = 0, eb_bytes = 0, eb_off = 0; std::vector<uint32_t> tl; id<MTLComputePipelineState> ps = nil; };
static void p14_bins(Engine& E, int s, const std::vector<int>& qs, const std::vector<uint32_t>& pq, const std::vector<uint32_t>& pn,
                     const std::vector<uint32_t>& pg, std::vector<P14Bin>& out, long& ndrop) {
  std::map<int, std::vector<int>> grp;
  for (int q : qs) { if (postlin_on(s) && pn[q] == 0) continue; grp[(int)pq[q] * 25 + (int)pn[q]].push_back(q); }
  for (auto& kv : grp) {
    const int Qk = kv.first / 25, Tk = kv.first % 25; size_t cnt = kv.second.size();
    id<MTLComputePipelineState> ps = E.pipe(E.Lp, "pol_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
    if (!ps) { ndrop += (long)cnt; continue; }
    {
      const int NPq = std::max(4, (Tk + 3) / 4 * 4), QPq = std::max(4, (Qk + 3) / 4 * 4);
      const size_t WSq = (4 * (size_t)NPq + 3 * (size_t)QPq + 6) / 4 * 4;
      const size_t shq = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
      std::vector<int> kp; for (int q : kv.second) { if ((shq + WSq + 2 * (size_t)pg[q]) * 4 > 32768) ndrop += 1; else kp.push_back(q); }
      kv.second.swap(kp); cnt = kv.second.size(); if (cnt == 0) continue; }
    uint32_t per = 0; for (int q : kv.second) per = std::max(per, pg[q]);
    const int NPk = std::max(4, (Tk + 3) / 4 * 4), QPk = std::max(4, (Qk + 3) / 4 * 4);
    const size_t WSk = (4 * (size_t)NPk + 3 * (size_t)QPk + 6) / 4 * 4;
    const size_t shb = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
    const size_t pers = shb + WSk + 2 * (size_t)per;
    if (pers * 4 > 32768) { ndrop += (long)cnt; continue; }
    int tpg = std::min<int>(8, (int)(32768 / (pers * 4)));
    tpg = std::min<int>(tpg, (int)ps.maxTotalThreadsPerThreadgroup / 32);
    if (tpg < 1) { ndrop += (long)cnt; continue; }
    if (h150p_on()) { const int t0 = tpg; tpg = h150_pick_pol(pers, cnt, (int)ps.maxTotalThreadsPerThreadgroup / 32, tpg); g_h150.polBins++; g_h150.polChanged += (tpg != t0); }
    size_t xf713 = 0;
    if (g_h713 && Tk >= 9) { id<MTLComputePipelineState> pp = E.pipe(E.Lpp, "polp_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
      if (pp && (int)pp.maxTotalThreadsPerThreadgroup >= 64 && (2 * pers + h713_xf(Tk)) * 4 <= 32768) { ps = pp; tpg = 2; xf713 = h713_xf(Tk); } }
    P14Bin g; g.side = s; g.Q = Qk; g.nth = Tk; g.per = per; g.tpg = tpg; g.ps = ps;
    g.nt = (cnt + (size_t)tpg - 1) / (size_t)tpg * (size_t)tpg; g.ntg = xf713 ? g.nt : g.nt / (size_t)tpg;
    g.tgb = ((size_t)tpg * pers + xf713) * 4; g.eb_bytes = g.ntg * (size_t)tpg * 2 * (size_t)(Qk + 1) * 64 * 4;
    g.tl.assign(g.nt + 1, 0xFFFFFFFFu); g.tl[0] = per;
    for (size_t k = 0; k < cnt; k++) g.tl[1 + k] = (uint32_t)kv.second[k];
    out.push_back(g);
  }
}


struct P14Buf { id<MTLBuffer> hdr, segs, th, rs, a, stio, yv; };


static const int g_rec_pol = getenv("REC_POLISH") ? atoi(getenv("REC_POLISH")) :1;
static std::vector<float> g_rp_th[2]; static std::vector<uint8_t> g_rp_ok; static std::atomic<long> g_rp_n{0}; static double g_rp_dsum = 0; static std::atomic<long> g_rp_dn{0};
static inline void rp_begin(int n) { if (!g_rec_pol) return; g_rp_ok.assign((size_t)n, 0); for (int s = 0; s < 2; s++) g_rp_th[s].resize((size_t)n * 32); }
static double g_ebmb_call = 0;
static double p14_run(Engine& E, std::vector<P14Bin>& bins, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp, id<MTLBuffer> tvec, float itmax, float skipf) {
  if (bins.empty()) return 0.0;
  if (h94_mode() >= 1) {
    std::vector<size_t> nat(bins.size()), tgb(bins.size()), cnt(bins.size()), cap(bins.size()); size_t before = 0, after = 0, ntr = 0;
    for (size_t b = 0; b < bins.size(); b++) {
      const P14Bin& g = bins[b]; nat[b] = g.ntg; tgb[b] = g.tgb; size_t c = 0;
      for (size_t k = 1; k < g.tl.size(); k++) if (g.tl[k] != 0xFFFFFFFFu) c++;
      cnt[b] = c; ntr += c; before += g.ntg;
    }
    h94_cap(cap, nat, tgb, cnt);
    for (size_t b = 0; b < bins.size(); b++) { P14Bin& g = bins[b]; g.ntg = cap[b]; g.eb_bytes = g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4; after += g.ntg; }
    fprintf(stderr, "H94 polish dispatch: bins %zu trees %zu TGs before %zu after %zu\n", bins.size(), ntr, before, after);
  }
  size_t eo = 0; for (size_t b = 0; b < bins.size(); b++) { bins[b].eb_off = eo; eo = al5(eo + bins[b].eb_bytes + ebt_bytes(bins[b].Q, bins[b].nth, bins[b].ntg, bins[b].tpg), 256); bins[b].wix = (uint32_t)b; }
  g_ebmb_call += (double)eo / 1048576.0;


  static const bool h82 = !getenv("H82") || atoi(getenv("H82")) != 0;
  static GBuf p82_eb, p82_ctr, p82_tl;
  id<MTLBuffer> ebB, ctr; std::vector<size_t> tlo(bins.size(), 0);
  if (h82) {
    size_t tlb = 0; for (size_t b = 0; b < bins.size(); b++) { tlo[b] = tlb; tlb = al5(tlb + bins[b].tl.size() * 4, 256); }
    h900_ens(p82_eb, "p82_eb", E.dev, eo, MPV); h900_ens(p82_ctr, "p82_ctr", E.dev, 4 * bins.size() + 16, MSH); h900_ens(p82_tl, "p82_tl", E.dev, tlb, MSH);
    ebB = p82_eb.b; ctr = p82_ctr.b;
    for (size_t b = 0; b < bins.size(); b++) std::memcpy((char*)p82_tl.b.contents + tlo[b], bins[b].tl.data(), bins[b].tl.size() * 4);
  } else {
    ebB = h900_nbl(E.dev, __LINE__, std::max<size_t>(eo, 256), MTLResourceStorageModePrivate);
    ctr = h900_nbl(E.dev, __LINE__, 4 * bins.size() + 16, MTLResourceStorageModeShared);
  }
  std::memset(ctr.contents, 0, 4 * bins.size() + 16);
  id<MTLCommandBuffer> cb = h182_cb(E.cq);
  id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
  for (size_t bi = 0; bi < bins.size(); bi++) { auto& g = bins[bi];
    const P14Buf& S = SBf[g.side];
    id<MTLBuffer> tl = h82 ? p82_tl.b : E.B(g.tl.data(), g.tl.size() * 4);
    [en setComputePipelineState:g.ps];
    id<MTLBuffer> bs[15] = {tok, S.hdr, psp, S.segs, S.yv, S.th, S.rs, S.stio, tl, nil, S.a, ebB, nil, ctr, tvec};
    for (int k = 0; k < 15; k++) if (bs[k]) [en setBuffer:bs[k] offset:(k == 11 ? g.eb_off : (k == 8 && h82 ? tlo[bi] : 0)) atIndex:k];
    if (g_h711) [en setBuffer:ebB offset:g.eb_off + g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4 atIndex:15];
    float tf[4] = {(float)g.Q * 9.0f * ldexpf(1.0f, -24), E.tg[g.side].floor_, itmax, skipf};
    [en setBytes:tf length:16 atIndex:9];
    uint32_t wq[2] = {g.wix, (uint32_t)g.nt}; [en setBytes:wq length:8 atIndex:12];
    [en setThreadgroupMemoryLength:g.tgb atIndex:0];
    if (h150c_on()) { g_h150.polDisp++; g_h150.polLaunch += (long)(g.ntg * (size_t)g.tpg); g_h150.polOcc += (long)h150_count_tl(g.tl); }
    [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
  }
  [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
  if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I14 polish command buffer failed\n"); h182_dump(cb, "I14 polish"); gpu_fail(); }
  return (cb.GPUEndTime - cb.GPUStartTime) * 1e3;
}


static bool eval_cols_t14(const View& V, const double* th, const double* tt, std::vector<double>& E) {
  const int n = V.n(), Q = V.Q; std::vector<double> v(n); E.assign((size_t)Q * NPT, 0);
  for (int p = 0; p < NPT; p++) { const double t = tt[p];
    for (int i = 0; i < n; i++) { int o = V.op(i); double A = (o >= 3 && o <= 15) ? v[V.a(i)] : 0, B = (o >= 3 && o <= 6) ? v[V.b(i)] : 0, x = 0;
      switch (o) { case 1: x = th[V.a(i)]; break; case 2: x = t; break; case 16: x = 1; break; case 3: x = A + B; break; case 4: x = A - B; break; case 5: x = A * B; break; case 6: x = A / B; break;
        case 7: x = sin(A); break; case 8: x = cos(A); break; case 9: x = exp(A); break; case 10: x = log(fabs(A)); break; case 11: x = A * A; break; case 12: x = A * A * A; break; case 13: x = sqrt(fabs(A)); break; case 14: x = cbrt(A); break; case 15: x = -A; break; }
      v[i] = x; }
    for (int q = 0; q < Q; q++) { double x = v[V.root(q)]; if (!std::isfinite(x)) return false; E[(size_t)q * NPT + p] = x; } }
  return true; }
static double g_minratio = 0;
static double vp14(const std::vector<double>& E, int Q, const double* y, double tau) {
  std::vector<double> Qm(E); std::vector<double> z(y, y + NPT); std::vector<char> keep(Q, 0); g_minratio = INFINITY;
  for (int k = 0; k < Q; k++) { double on = 0; for (int p = 0; p < NPT; p++) on += E[(size_t)k * NPT + p] * E[(size_t)k * NPT + p]; on = sqrt(on);
    for (int i = 0; i < k; i++) if (keep[i]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)i * NPT + p] * Qm[(size_t)k * NPT + p]; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] -= s * Qm[(size_t)i * NPT + p]; }
    double nn = 0; for (int p = 0; p < NPT; p++) nn += Qm[(size_t)k * NPT + p] * Qm[(size_t)k * NPT + p]; nn = sqrt(nn);
    if (on > 0) g_minratio = std::min(g_minratio, nn / on);
    keep[k] = nn > tau * on; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] = keep[k] ? Qm[(size_t)k * NPT + p] / nn : 0; }
  for (int k = 0; k < Q; k++) if (keep[k]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)k * NPT + p] * z[p]; for (int p = 0; p < NPT; p++) z[p] -= s * Qm[(size_t)k * NPT + p]; }
  double rss = 0; for (int p = 0; p < NPT; p++) rss += z[p] * z[p]; return rss; }

static std::vector<double> vp14z(const std::vector<double>& E, int Q, const double* y, double tau) {
  std::vector<double> Qm(E); std::vector<double> z(y, y + NPT); std::vector<char> keep(Q, 0);
  for (int k = 0; k < Q; k++) { double on = 0; for (int p = 0; p < NPT; p++) on += E[(size_t)k * NPT + p] * E[(size_t)k * NPT + p]; on = sqrt(on);
    for (int i = 0; i < k; i++) if (keep[i]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)i * NPT + p] * Qm[(size_t)k * NPT + p]; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] -= s * Qm[(size_t)i * NPT + p]; }
    double nn = 0; for (int p = 0; p < NPT; p++) nn += Qm[(size_t)k * NPT + p] * Qm[(size_t)k * NPT + p]; nn = sqrt(nn);
    keep[k] = nn > tau * on; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] = keep[k] ? Qm[(size_t)k * NPT + p] / nn : 0; }
  for (int k = 0; k < Q; k++) if (keep[k]) { double s = 0; for (int p = 0; p < NPT; p++) s += Qm[(size_t)k * NPT + p] * z[p]; for (int p = 0; p < NPT; p++) z[p] -= s * Qm[(size_t)k * NPT + p]; }
  return z; }
static const int g_jprof = getenv("J_PROF") ? atoi(getenv("J_PROF")) : 0;
static const int g_jchk = getenv("J_CHK") ? atoi(getenv("J_CHK")) : 0;
static std::atomic<uint64_t> g_jchk_bad{0}, g_jchk_n{0}, g_jchk_qbad{0}, g_jchk_qn{0};
static std::atomic<uint64_t> g_jp_ev{0}, g_jp_qr{0}, g_jp_nn{0}, g_jp_qq{0}, g_jp_sd{0}, g_jp_res{0}, g_jp_sc{0}, g_jp_n{0};
static const Sketch::CisPerCellEnergyInput* g_sicis_en = nullptr;


static inline void jq_dots(const double* const* a, const double* const* b, double* out, int m) {
  for (int q = 0; q < m; q++) { double sd = 0; const double* x = a[q]; const double* y = b[q]; for (int p = 0; p < NPT; p++) sd += x[p] * y[p]; out[q] = sd; }
}
static void qr48_pair(const double* const Ev[2], const int Qs[2], const double* const y[2], double* const zz[2]) {
  double Qm[2][16 * NPT]; char keep[2][16]; double on[2][16];
  { const double* a[32]; double r[32]; int m = 0;
    for (int s = 0; s < 2; s++) { std::memcpy(Qm[s], Ev[s], sizeof(double) * (size_t)Qs[s] * NPT); for (int p = 0; p < NPT; p++) zz[s][p] = y[s][p]; for (int c = 0; c < 16; c++) keep[s][c] = 0;
      for (int c = 0; c < Qs[s]; c++) a[m++] = Ev[s] + (size_t)c * NPT; }
    jq_dots(a, a, r, m); m = 0;
    for (int s = 0; s < 2; s++) for (int c = 0; c < Qs[s]; c++) on[s][c] = sqrt(r[m++]); }
  double tau[2] = {(double)Qs[0] * 9.0 * ldexp(1.0, -24), (double)Qs[1] * 9.0 * ldexp(1.0, -24)};
  const int QM = std::max(Qs[0], Qs[1]);
  for (int i = 0; i < QM; i++) {
    { const double* a[2]; double r[2]; int sd_[2], m = 0;
      for (int s = 0; s < 2; s++) if (i < Qs[s]) { a[m] = Qm[s] + (size_t)i * NPT; sd_[m] = s; m++; }
      jq_dots(a, a, r, m);
      for (int q = 0; q < m; q++) { const int s = sd_[q]; const double nn = sqrt(r[q]); double* c = Qm[s] + (size_t)i * NPT;
        keep[s][i] = nn > tau[s] * on[s][i]; for (int p = 0; p < NPT; p++) c[p] = keep[s][i] ? c[p] / nn : 0; } }
    const double* a[32]; const double* b[32]; double r[32]; int ss[32], jj[32], m = 0;
    for (int s = 0; s < 2; s++) if (i < Qs[s] && keep[s][i]) for (int j = i + 1; j < Qs[s]; j++) { a[m] = Qm[s] + (size_t)i * NPT; b[m] = Qm[s] + (size_t)j * NPT; ss[m] = s; jj[m] = j; m++; }
    jq_dots(a, b, r, m);
    for (int q = 0; q < m; q++) { double* c = Qm[ss[q]] + (size_t)jj[q] * NPT; const double* qi = a[q]; const double sd = r[q]; for (int p = 0; p < NPT; p++) c[p] -= sd * qi[p]; }
  }
  { int idx[2][16], nk[2] = {0, 0}; for (int s = 0; s < 2; s++) for (int c = 0; c < Qs[s]; c++) if (keep[s][c]) idx[s][nk[s]++] = c;
    for (int t = 0; t < std::max(nk[0], nk[1]); t++) {
      const double* a[2]; const double* b[2]; double r[2]; int sd_[2], m = 0;
      for (int s = 0; s < 2; s++) if (t < nk[s]) { a[m] = Qm[s] + (size_t)idx[s][t] * NPT; b[m] = zz[s]; sd_[m] = s; m++; }
      jq_dots(a, b, r, m);
      for (int q = 0; q < m; q++) { const int s = sd_[q]; double* z = zz[s]; const double* qk = a[q]; const double sd = r[q]; for (int p = 0; p < NPT; p++) z[p] -= sd * qk[p]; } } }
}
static bool sicis_res48_fast(Engine& E, const Unique& UQ, const std::vector<int>& ct, const std::vector<int>& K, size_t k,
                             const P14Buf* SBf, id<MTLBuffer> bTV, double z[2][NPT], const float* const* jth = nullptr, const float* jtv = nullptr, int only = -1) {
  thread_local std::vector<double> vb;

  double tb[NPT]; if (g_yfx) { for (int p = 0; p < NPT; p++) tb[p] = g_yfx_t[p]; } else if (jtv) { for (int p = 0; p < NPT; p++) tb[p] = (double)jtv[48 * k + p]; } else for (int p = 0; p < NPT; p++) tb[p] = (double)((const float*)bTV.contents)[48 * k + p];
  static const int jqr = getenv("J_QR") ? atoi(getenv("J_QR")) : 1;
  double EvS[2][16 * NPT]; int QS[2] = {0, 0}; double zref[2][NPT];
  for (int s = 0; s < 2; s++) {
    if (only >= 0 && s != only) continue;
    const int j = ct[(size_t)K[k] * 2 + s];
    const Tree& tr = UQ.U[j]; View V(tr.rec.data());
    const int Q = tr.Q; if (Q > 16 || tr.nth > 32) return false;
    double th[32]; if (jth) { for (int q = 0; q < tr.nth; q++) th[q] = (double)jth[s][32 * k + q]; } else for (int q = 0; q < tr.nth; q++) th[q] = (double)((const float*)SBf[s].th.contents)[32 * k + q];
    const int n = V.n(); if ((int)vb.size() < n) vb.resize(n); double* v = vb.data(); for (int i = 0; i < n; i++) v[i] = 0;
    double* Ev = EvS[s]; QS[s] = Q;
    const uint64_t jq0 = g_jprof ? mach_absolute_time() : 0;


    static const int jev = getenv("J_EVAL") ? atoi(getenv("J_EVAL")) : 1;
    bool done2 = false;
    if (jev) {
      thread_local std::vector<double> vm; if ((int)vm.size() < n * NPT) vm.resize((size_t)n * NPT);
      double* __restrict M = vm.data(); bool okn = true;
      for (int i = 0; i < n && okn; i++) {
        const int o = V.op(i); const int ia = V.a(i), ib = V.b(i); double* __restrict d = M + (size_t)i * NPT;
        if ((o >= 3 && o <= 15 && ia >= i) || (o >= 3 && o <= 6 && ib >= i)) { okn = false; break; }
        const double* A = M + (size_t)ia * NPT; const double* B = M + (size_t)ib * NPT;
        switch (o) {
          case 1: { const double c = th[ia]; for (int p = 0; p < NPT; p++) d[p] = c; } break;
          case 2: for (int p = 0; p < NPT; p++) d[p] = tb[p]; break;
          case 16: for (int p = 0; p < NPT; p++) d[p] = 1; break;
          case 3: for (int p = 0; p < NPT; p++) d[p] = A[p] + B[p]; break;
          case 4: for (int p = 0; p < NPT; p++) d[p] = A[p] - B[p]; break;
          case 5: for (int p = 0; p < NPT; p++) d[p] = A[p] * B[p]; break;
          case 6: for (int p = 0; p < NPT; p++) d[p] = A[p] / B[p]; break;
          case 7: for (int p = 0; p < NPT; p++) d[p] = sin(A[p]); break;
          case 8: for (int p = 0; p < NPT; p++) d[p] = cos(A[p]); break;
          case 9: for (int p = 0; p < NPT; p++) d[p] = exp(A[p]); break;
          case 10: for (int p = 0; p < NPT; p++) d[p] = log(fabs(A[p])); break;
          case 11: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p]; break;
          case 12: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p] * A[p]; break;
          case 13: for (int p = 0; p < NPT; p++) d[p] = sqrt(fabs(A[p])); break;
          case 14: for (int p = 0; p < NPT; p++) d[p] = cbrt(A[p]); break;
          case 15: for (int p = 0; p < NPT; p++) d[p] = -A[p]; break;
          default: for (int p = 0; p < NPT; p++) d[p] = 0; break;
        }
      }
      if (okn) {
        for (int q = 0; q < Q; q++) { const double* R = M + (size_t)V.root(q) * NPT; for (int p = 0; p < NPT; p++) { const double x = R[p]; if (!std::isfinite(x)) return false; Ev[(size_t)q * NPT + p] = x; } }
        done2 = true;
        if (g_jchk) {
          static thread_local std::vector<double> vo; if ((int)vo.size() < n) vo.resize(n); double* w = vo.data(); for (int i = 0; i < n; i++) w[i] = 0; uint64_t bad = 0;
          for (int p = 0; p < NPT; p++) { const double t = tb[p];
            for (int i = 0; i < n; i++) { int o = V.op(i); double A = (o >= 3 && o <= 15) ? w[V.a(i)] : 0, B = (o >= 3 && o <= 6) ? w[V.b(i)] : 0, x = 0;
              switch (o) { case 1: x = th[V.a(i)]; break; case 2: x = t; break; case 16: x = 1; break; case 3: x = A + B; break; case 4: x = A - B; break; case 5: x = A * B; break; case 6: x = A / B; break;
                case 7: x = sin(A); break; case 8: x = cos(A); break; case 9: x = exp(A); break; case 10: x = log(fabs(A)); break; case 11: x = A * A; break; case 12: x = A * A * A; break; case 13: x = sqrt(fabs(A)); break; case 14: x = cbrt(A); break; case 15: x = -A; break; }
              w[i] = x; }
            for (int q = 0; q < Q; q++) { const double x = w[V.root(q)]; if (std::memcmp(&x, &Ev[(size_t)q * NPT + p], 8) != 0) bad++; } }
          g_jchk_bad.fetch_add(bad, std::memory_order_relaxed); g_jchk_n.fetch_add(1, std::memory_order_relaxed);
        }
      }
    }
    if (!done2)
    for (int p = 0; p < NPT; p++) { const double t = tb[p];
      for (int i = 0; i < n; i++) { int o = V.op(i); double A = (o >= 3 && o <= 15) ? v[V.a(i)] : 0, B = (o >= 3 && o <= 6) ? v[V.b(i)] : 0, x = 0;
        switch (o) { case 1: x = th[V.a(i)]; break; case 2: x = t; break; case 16: x = 1; break; case 3: x = A + B; break; case 4: x = A - B; break; case 5: x = A * B; break; case 6: x = A / B; break;
          case 7: x = sin(A); break; case 8: x = cos(A); break; case 9: x = exp(A); break; case 10: x = log(fabs(A)); break; case 11: x = A * A; break; case 12: x = A * A * A; break; case 13: x = sqrt(fabs(A)); break; case 14: x = cbrt(A); break; case 15: x = -A; break; }
        v[i] = x; }
      for (int q = 0; q < Q; q++) { double x = v[V.root(q)]; if (!std::isfinite(x)) return false; Ev[(size_t)q * NPT + p] = x; } }
    const uint64_t jq1 = g_jprof ? mach_absolute_time() : 0;
    if (!jqr || g_jchk) {
    const double tau = (double)Q * 9.0 * ldexp(1.0, -24); const double* y = g_yfx ? g_yfx_y[s] : E.tg[s].y;
    double Qm[16 * NPT]; std::memcpy(Qm, Ev, sizeof(double) * (size_t)Q * NPT);
    double* zz = jqr ? zref[s] : z[s]; for (int p = 0; p < NPT; p++) zz[p] = y[p]; char keep[16] = {0};
    for (int kk = 0; kk < Q; kk++) { double on = 0; for (int p = 0; p < NPT; p++) on += Ev[(size_t)kk * NPT + p] * Ev[(size_t)kk * NPT + p]; on = sqrt(on);
      for (int i = 0; i < kk; i++) if (keep[i]) { double sd = 0; for (int p = 0; p < NPT; p++) sd += Qm[(size_t)i * NPT + p] * Qm[(size_t)kk * NPT + p]; for (int p = 0; p < NPT; p++) Qm[(size_t)kk * NPT + p] -= sd * Qm[(size_t)i * NPT + p]; }
      double nn = 0; for (int p = 0; p < NPT; p++) nn += Qm[(size_t)kk * NPT + p] * Qm[(size_t)kk * NPT + p]; nn = sqrt(nn);
      keep[kk] = nn > tau * on; for (int p = 0; p < NPT; p++) Qm[(size_t)kk * NPT + p] = keep[kk] ? Qm[(size_t)kk * NPT + p] / nn : 0; }
    for (int kk = 0; kk < Q; kk++) if (keep[kk]) { double sd = 0; for (int p = 0; p < NPT; p++) sd += Qm[(size_t)kk * NPT + p] * zz[p]; for (int p = 0; p < NPT; p++) zz[p] -= sd * Qm[(size_t)kk * NPT + p]; }
    }
    if (g_jprof) { const uint64_t jq2 = mach_absolute_time(); g_jp_ev.fetch_add(jq1 - jq0, std::memory_order_relaxed); g_jp_qr.fetch_add(jq2 - jq1, std::memory_order_relaxed); g_jp_nn.fetch_add(n, std::memory_order_relaxed); g_jp_qq.fetch_add(Q, std::memory_order_relaxed); g_jp_sd.fetch_add(1, std::memory_order_relaxed); }
  }
  if (jqr) {
    const uint64_t jq3 = g_jprof ? mach_absolute_time() : 0;
    const double* evp[2] = {EvS[0], EvS[1]}; const double* yp[2] = {g_yfx ? g_yfx_y[0] : E.tg[0].y, g_yfx ? g_yfx_y[1] : E.tg[1].y}; double* zp[2] = {z[0], z[1]};
    qr48_pair(evp, QS, yp, zp);
    if (g_jprof) g_jp_qr.fetch_add(mach_absolute_time() - jq3, std::memory_order_relaxed);
    if (g_jchk) { uint64_t bad = 0; for (int s = 0; s < 2; s++) bad += std::memcmp(z[s], zref[s], sizeof(double) * NPT) != 0; g_jchk_qbad.fetch_add(bad, std::memory_order_relaxed); g_jchk_qn.fetch_add(2, std::memory_order_relaxed); }
  }
  return true;
}


static const int g_yfxfast = (getenv("YFX_FAST") && !(getenv("POST_LIN") && atoi(getenv("POST_LIN")) != 0)) ? atoi(getenv("YFX_FAST")) : 0;
static std::atomic<uint64_t> g_yf_hit{0}, g_yf_miss{0}, g_yf_bad{0}, g_yf_chk{0}, g_yf_schk{0}, g_yf_sbad{0}, g_yf_zchk{0}, g_yf_zbad{0};
static bool sicis_res48_yfxf(Engine& E, const Unique& UQ, const std::vector<int>& ct, const std::vector<int>& K, size_t k,
                             const P14Buf* SBf, id<MTLBuffer> bTV, float* out, const float* const* jth, const float* jtv) {
  struct CEnt { uint64_t ep = ~0ull; std::vector<uint32_t> rec; double z[NPT]; };
  thread_local CEnt cache[8]; thread_local unsigned rr = 0;
  double zf[2][NPT];
  const Tree* tr[2]; int cs = -1; bool hit[2] = {false, false};
  for (int s = 0; s < 2; s++) tr[s] = &UQ.U[ct[(size_t)K[k] * 2 + s]];
  for (int s = 0; s < 2; s++) if (tr[s]->nth == 0 && tr[s]->Q <= 16) {
    for (int c = 0; c < 8; c++) { const CEnt& e = cache[c]; if (e.ep == g_yfx_epoch && e.rec.size() == tr[s]->rec.size() && !std::memcmp(e.rec.data(), tr[s]->rec.data(), 4 * e.rec.size())) { std::memcpy(zf[s], e.z, sizeof(double) * NPT); hit[s] = true; break; } }
  }
  if (hit[0] && hit[1]) { g_yf_hit.fetch_add(2, std::memory_order_relaxed); }
  else if (hit[0] || hit[1]) { const int o = hit[0] ? 1 : 0, h = 1 - o; g_yf_hit.fetch_add(1, std::memory_order_relaxed);
    double sv[NPT]; std::memcpy(sv, zf[h], sizeof sv);
    if (!sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, zf, jth, jtv, o)) return false;
    std::memcpy(zf[h], sv, sizeof sv);
    if (tr[o]->nth == 0 && tr[o]->Q <= 16) { CEnt& e = cache[rr++ & 7]; e.ep = g_yfx_epoch; e.rec = tr[o]->rec; std::memcpy(e.z, zf[o], sizeof(double) * NPT); g_yf_miss.fetch_add(1, std::memory_order_relaxed); }
  } else {
    if (!sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, zf, jth, jtv)) return false;
    for (int s = 0; s < 2; s++) if (tr[s]->nth == 0 && tr[s]->Q <= 16) { CEnt& e = cache[rr++ & 7]; e.ep = g_yfx_epoch; e.rec = tr[s]->rec; std::memcpy(e.z, zf[s], sizeof(double) * NPT); g_yf_miss.fetch_add(1, std::memory_order_relaxed); }
  }
  for (int p = 0; p < NPT; p++) out[p] = (float)(zf[0][p] * g_yfx_n[p] + zf[1][p] * g_yfx_n[NPT + p]);
  return true;
}


#if defined(__aarch64__)
#include <arm_neon.h>


static const int g_yfx_l16 = getenv("YFX_L16") ? atoi(getenv("YFX_L16")) : 0;
template <int L, bool FUSE> static inline void dotL(const double* A, const double* B, double* o) {
  float64x2_t a[L / 2]; for (int v = 0; v < L / 2; v++) a[v] = vdupq_n_f64(0.0);
  for (int p = 0; p < NPT; p++) { const double* x = A + p * L; const double* y = B + p * L;
    for (int v = 0; v < L / 2; v++) { if (FUSE) a[v] = vfmaq_f64(a[v], vld1q_f64(x + 2 * v), vld1q_f64(y + 2 * v)); else a[v] = vaddq_f64(a[v], vmulq_f64(vld1q_f64(x + 2 * v), vld1q_f64(y + 2 * v))); } }
  for (int v = 0; v < L / 2; v++) vst1q_f64(o + 2 * v, a[v]);
}
static const int g_yfx_dm = getenv("YFX_DM") ? atoi(getenv("YFX_DM")) : 3, g_yfx_um = getenv("YFX_UM") ? atoi(getenv("YFX_UM")) : 3;
template <int L> static inline void dotk(int kind, const double* A, const double* B, double* o) { if ((g_yfx_dm >> kind) & 1) dotL<L, true>(A, B, o); else dotL<L, false>(A, B, o); }
static inline float64x2_t updk(bool fz, float64x2_t c, float64x2_t sd, float64x2_t q) { return fz ? vfmsq_f64(c, sd, q) : vsubq_f64(c, vmulq_f64(sd, q)); }

template <int L> static void qr48_bL(const double* const* Ev, int nb, int Q, const double* y, double* const* zo) {
  constexpr int NV = L / 2;
  thread_local std::vector<double> Qv; if (Qv.size() < (size_t)16 * NPT * L) Qv.resize((size_t)16 * NPT * L);
  double* Qm = Qv.data(); alignas(16) double zz[NPT * L];
  for (int q = 0; q < Q; q++) for (int b = 0; b < L; b++) { const double* src = Ev[b < nb ? b : 0] + (size_t)q * NPT; double* dst = Qm + (size_t)q * NPT * L + b; for (int p = 0; p < NPT; p++) dst[p * L] = src[p]; }
  for (int p = 0; p < NPT; p++) for (int b = 0; b < L; b++) zz[p * L + b] = y[p];
  const double tau = (double)Q * 9.0 * ldexp(1.0, -24);
  alignas(16) double on[16][L], r[L];
  for (int q = 0; q < Q; q++) { dotk<L>(0, Qm + (size_t)q * NPT * L, Qm + (size_t)q * NPT * L, r); for (int b = 0; b < L; b++) on[q][b] = sqrt(r[b]); }
  for (int i = 0; i < Q; i++) {
    double* ci = Qm + (size_t)i * NPT * L;
    dotk<L>(1, ci, ci, r);
    uint64x2_t kp[NV]; float64x2_t nv[NV];
    for (int v = 0; v < NV; v++) { nv[v] = vsqrtq_f64(vld1q_f64(r + 2 * v)); kp[v] = vcgtq_f64(nv[v], vmulq_f64(vdupq_n_f64(tau), vld1q_f64(on[i] + 2 * v))); }
    for (int p = 0; p < NPT; p++) for (int v = 0; v < NV; v++) { const float64x2_t x = vld1q_f64(ci + p * L + 2 * v);
      vst1q_f64(ci + p * L + 2 * v, vreinterpretq_f64_u64(vandq_u64(kp[v], vreinterpretq_u64_f64(vdivq_f64(x, nv[v]))))); }
    for (int j = i + 1; j < Q; j++) {
      double* cj = Qm + (size_t)j * NPT * L; dotk<L>(2, ci, cj, r);
      float64x2_t sd[NV]; for (int v = 0; v < NV; v++) sd[v] = vld1q_f64(r + 2 * v);
      for (int p = 0; p < NPT; p++) for (int v = 0; v < NV; v++) vst1q_f64(cj + p * L + 2 * v, updk(g_yfx_um & 1, vld1q_f64(cj + p * L + 2 * v), sd[v], vld1q_f64(ci + p * L + 2 * v)));
    }
  }
  for (int t = 0; t < Q; t++) {
    const double* ct_ = Qm + (size_t)t * NPT * L; dotk<L>(3, ct_, zz, r);
    float64x2_t sd[NV]; for (int v = 0; v < NV; v++) sd[v] = vld1q_f64(r + 2 * v);
    for (int p = 0; p < NPT; p++) for (int v = 0; v < NV; v++) vst1q_f64(zz + p * L + 2 * v, updk((g_yfx_um >> 1) & 1, vld1q_f64(zz + p * L + 2 * v), sd[v], vld1q_f64(ct_ + p * L + 2 * v)));
  }
  for (int b = 0; b < nb; b++) for (int p = 0; p < NPT; p++) zo[b][p] = zz[p * L + b];
}

static bool yfx_eval_side(const Tree& tr, const double* th, double* Ev) {
  View V(tr.rec.data()); const int n = V.n(), Q = tr.Q;
  thread_local std::vector<double> vm; if ((int)vm.size() < n * NPT) vm.resize((size_t)n * NPT);
  double* __restrict M = vm.data();
  for (int i = 0; i < n; i++) {
    const int o = V.op(i); const int ia = V.a(i), ib = V.b(i); double* __restrict d = M + (size_t)i * NPT;
    if ((o >= 3 && o <= 15 && ia >= i) || (o >= 3 && o <= 6 && ib >= i)) return false;
    const double* A = M + (size_t)ia * NPT; const double* B = M + (size_t)ib * NPT; const double* tb = g_yfx_t;
    switch (o) {
      case 1: { const double c = th[ia]; for (int p = 0; p < NPT; p++) d[p] = c; } break;
      case 2: for (int p = 0; p < NPT; p++) d[p] = tb[p]; break;
      case 16: for (int p = 0; p < NPT; p++) d[p] = 1; break;
      case 3: for (int p = 0; p < NPT; p++) d[p] = A[p] + B[p]; break;
      case 4: for (int p = 0; p < NPT; p++) d[p] = A[p] - B[p]; break;
      case 5: for (int p = 0; p < NPT; p++) d[p] = A[p] * B[p]; break;
      case 6: for (int p = 0; p < NPT; p++) d[p] = A[p] / B[p]; break;
      case 7: for (int p = 0; p < NPT; p++) d[p] = sin(A[p]); break;
      case 8: for (int p = 0; p < NPT; p++) d[p] = cos(A[p]); break;
      case 9: for (int p = 0; p < NPT; p++) d[p] = exp(A[p]); break;
      case 10: for (int p = 0; p < NPT; p++) d[p] = log(fabs(A[p])); break;
      case 11: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p]; break;
      case 12: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p] * A[p]; break;
      case 13: for (int p = 0; p < NPT; p++) d[p] = sqrt(fabs(A[p])); break;
      case 14: for (int p = 0; p < NPT; p++) d[p] = cbrt(A[p]); break;
      case 15: for (int p = 0; p < NPT; p++) d[p] = -A[p]; break;
      default: for (int p = 0; p < NPT; p++) d[p] = 0; break;
    }
  }
  for (int q = 0; q < Q; q++) { const double* R = M + (size_t)V.root(q) * NPT; for (int p = 0; p < NPT; p++) { const double x = R[p]; if (!std::isfinite(x)) return false; Ev[(size_t)q * NPT + p] = x; } }
  return true;
}

static void yfx_chunk(Engine& E, const Unique& UQ, const std::vector<int>& ct, const std::vector<int>& K, size_t lo, size_t hi, const uint8_t* need,
                      const P14Buf* SBf, id<MTLBuffer> bTV, const float* const* jth, const float* jtv, float* resOut, float* ceOut, uint8_t* have) {
  const int n = (int)(hi - lo);
  struct CEnt { uint64_t ep = ~0ull; std::vector<uint32_t> rec; double z[NPT]; };
  thread_local CEnt cache[8]; thread_local unsigned rr = 0;
  struct Item { int cand, s, Q; size_t off; };
  thread_local std::vector<Item> items; thread_local std::vector<double> Evb, zb; thread_local std::vector<uint8_t> st;
  items.clear(); st.assign(n, 0); zb.resize((size_t)n * 2 * NPT);
  const uint64_t t0 = g_jprof ? mach_absolute_time() : 0;
  size_t evTot = 0;
  for (int c = 0; c < n; c++) {
    have[c] = 0; if (!need[c]) continue;
    const size_t k = lo + c; const Tree* tr[2] = {&UQ.U[ct[(size_t)K[k] * 2]], &UQ.U[ct[(size_t)K[k] * 2 + 1]]};
    if (tr[0]->Q > 16 || tr[1]->Q > 16 || tr[0]->nth > 32 || tr[1]->nth > 32) { st[c] = 2; continue; }
    st[c] = 1;
    for (int s = 0; s < 2 && st[c] == 1; s++) {
      double* zs = zb.data() + ((size_t)c * 2 + s) * NPT;
      if (tr[s]->nth == 0) {
        const CEnt* ce = nullptr; for (int q = 0; q < 8; q++) { const CEnt& e = cache[q]; if (e.ep == g_yfx_epoch && e.rec.size() == tr[s]->rec.size() && !std::memcmp(e.rec.data(), tr[s]->rec.data(), 4 * e.rec.size())) { ce = &e; break; } }
        if (ce) { std::memcpy(zs, ce->z, sizeof(double) * NPT); g_yf_hit.fetch_add(1, std::memory_order_relaxed); }
        else { double tz[2][NPT]; if (!sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, tz, jth, jtv, s)) { st[c] = 2; break; }
          std::memcpy(zs, tz[s], sizeof(double) * NPT); CEnt& e = cache[rr++ & 7]; e.ep = g_yfx_epoch; e.rec = tr[s]->rec; std::memcpy(e.z, tz[s], sizeof(double) * NPT); g_yf_miss.fetch_add(1, std::memory_order_relaxed); }
      } else { items.push_back({c, s, tr[s]->Q, evTot}); evTot += (size_t)tr[s]->Q * NPT; }
    }
  }
  Evb.resize(evTot + 8);
  const uint64_t t1 = g_jprof ? mach_absolute_time() : 0;
  for (const Item& it : items) { if (st[it.cand] != 1) continue; const size_t k = lo + it.cand; const Tree& tr = UQ.U[ct[(size_t)K[k] * 2 + it.s]];
    double th[32]; for (int q = 0; q < tr.nth; q++) th[q] = (double)jth[it.s][32 * k + q];
    if (!yfx_eval_side(tr, th, Evb.data() + it.off)) st[it.cand] = 2; }
  const uint64_t t2 = g_jprof ? mach_absolute_time() : 0;
  { thread_local std::vector<int> ord; ord.assign(items.size(), 0); int cnt[2 * 32 + 1] = {0};
    for (const Item& it : items) if (st[it.cand] == 1) cnt[it.s * 32 + it.Q + 1]++;
    for (int i = 0; i < 64; i++) cnt[i + 1] += cnt[i];
    const int nord = cnt[64]; int pos[64]; for (int i = 0; i < 64; i++) pos[i] = cnt[i];
    for (int i = 0; i < (int)items.size(); i++) if (st[items[i].cand] == 1) ord[pos[items[i].s * 32 + items[i].Q]++] = i;
    for (int g = 0; g < nord;) {
      const Item& f = items[ord[g]]; int e = g; while (e < nord && items[ord[e]].s == f.s && items[ord[e]].Q == f.Q) e++;
      for (int a = g; a < e;) { const int L = (e - a > 8 && g_yfx_l16) ? 16 : 8; const int nb = std::min(L, e - a); const double* evp[16]; double* zp[16];
        for (int b = 0; b < nb; b++) { const Item& it = items[ord[a + b]]; evp[b] = Evb.data() + it.off; zp[b] = zb.data() + ((size_t)it.cand * 2 + it.s) * NPT; }
        if (L == 16) qr48_bL<16>(evp, nb, f.Q, g_yfx_y[f.s], zp); else qr48_bL<8>(evp, nb, f.Q, g_yfx_y[f.s], zp);
        a += nb; }
      g = e; } }
  const uint64_t t3 = g_jprof ? mach_absolute_time() : 0;
  if (g_yfxfast == 2) for (int c = 0; c < n; c++) if (st[c] == 1) { const size_t k = lo + c; double zr[2][NPT];
      if (sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, zr, jth, jtv)) { const double* z0 = zb.data() + (size_t)c * 2 * NPT;
        for (int s = 0; s < 2; s++) { g_yf_zchk.fetch_add(1, std::memory_order_relaxed); if (std::memcmp(zr[s], z0 + s * NPT, sizeof(double) * NPT)) g_yf_zbad.fetch_add(1, std::memory_order_relaxed); } } }
  thread_local std::vector<float> rc; thread_local std::vector<int> ci; thread_local std::vector<double> sc; thread_local std::vector<uint8_t> fbv;
  ci.clear(); rc.resize((size_t)n * NPT);
  for (int c = 0; c < n; c++) if (st[c] == 1) { const double* z0 = zb.data() + (size_t)c * 2 * NPT; const double* z1 = z0 + NPT; float* o = rc.data() + ci.size() * NPT;
      for (int p = 0; p < NPT; p++) o[p] = (float)(z0[p] * g_yfx_n[p] + z1[p] * g_yfx_n[NPT + p]);
      ci.push_back(c); }
  const uint64_t t4 = g_jprof ? mach_absolute_time() : 0;
  const int m = (int)ci.size(); sc.resize(m + 4); fbv.assign(m + 4, 0);
  bool fbs[1024]; std::vector<uint8_t>* dummy = nullptr; (void)dummy;
  thread_local std::vector<char> fbc; fbc.assign(m + 4, 0);
  bool okb = m > 0 && sicis::scoreBatch48(rc.data(), m, *g_sicis_en, sc.data(), reinterpret_cast<bool*>(fbc.data())); (void)fbs;
  for (int i = 0; i < m; i++) { const int c = ci[i]; const float* r = rc.data() + (size_t)i * NPT; std::memcpy(resOut + (size_t)c * NPT, r, sizeof(float) * NPT);
    const double sv = (okb && !fbc[i]) ? sc[i] : sicis::score48(r, *g_sicis_en); ceOut[c] = (float)std::min(sv, 1.0e30); have[c] = 1; }
  if (g_jprof) { const uint64_t t5 = mach_absolute_time(); g_jp_res.fetch_add(t4 - t0, std::memory_order_relaxed); g_jp_sc.fetch_add(t5 - t4, std::memory_order_relaxed); g_jp_n.fetch_add(m, std::memory_order_relaxed);
    g_jp_ev.fetch_add(t2 - t1, std::memory_order_relaxed); g_jp_qr.fetch_add(t3 - t2, std::memory_order_relaxed); g_jp_sd.fetch_add(m, std::memory_order_relaxed); }
}
#endif
static bool sicis_res48(Engine& E, const Unique& UQ, const std::vector<int>& ct, const std::vector<int>& K, size_t k,
                        const P14Buf* SBf, id<MTLBuffer> bTV, float* out, const float* const* jth = nullptr, const float* jtv = nullptr) {
  double z[2][NPT];
  static const bool fastOn = sicis::fastEnabled();
  static const bool noS1 = !getenv("FJN_S1");
  if (g_yfx && g_yfxfast && fastOn && NPT == 48 && noS1) {
    if (sicis_res48_yfxf(E, UQ, ct, K, k, SBf, bTV, out, jth, jtv)) {
      if (g_yfxfast == 2) { float ref[NPT]; double zf[2][NPT]; if (sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, zf, jth, jtv)) { for (int p = 0; p < NPT; p++) ref[p] = (float)(zf[0][p] * g_yfx_n[p] + zf[1][p] * g_yfx_n[NPT + p]);
          g_yf_chk.fetch_add(1, std::memory_order_relaxed); if (std::memcmp(ref, out, sizeof ref)) g_yf_bad.fetch_add(1, std::memory_order_relaxed); } }
      return true; } }
  if (fastOn && NPT == 48 && noS1) { double zf[2][NPT]; if (sicis_res48_fast(E, UQ, ct, K, k, SBf, bTV, zf, jth, jtv)) {
      if (g_yfx) { for (int p = 0; p < NPT; p++) out[p] = (float)(zf[0][p] * g_yfx_n[p] + zf[1][p] * g_yfx_n[NPT + p]); return true; }
      if (g_nres) { for (int p = 0; p < NPT; p++) out[p] = (float)(zf[0][p] * (double)g_nres_n[p] + zf[1][p] * (double)g_nres_n[NPT + p]); return true; }
      for (int p = 0; p < NPT; p++) { const double h = std::hypot(zf[0][p], zf[1][p]); out[p] = (float)((zf[0][p] + zf[1][p] < 0) ? -h : h); } return true; } }
  if (g_yfx) return false;
  for (int s = 0; s < 2; s++) {
    const int j = ct[(size_t)K[k] * 2 + s];
    const Tree& tr = UQ.U[j]; View V(tr.rec.data());
    double th[32]; for (int q = 0; q < tr.nth; q++) th[q] = (double)((const float*)SBf[s].th.contents)[32 * k + q];
    double tb[NPT]; for (int p = 0; p < NPT; p++) tb[p] = (double)((const float*)bTV.contents)[48 * k + p];
    std::vector<double> Ecol;
    if (!eval_cols_t14(V, th, tb, Ecol)) return false;
    std::vector<double> zs = vp14z(Ecol, tr.Q, E.tg[s].y, (double)tr.Q * 9.0 * ldexp(1.0, -24));
    for (int p = 0; p < NPT; p++) z[s][p] = zs[p];
  }
  static const bool s1on = getenv("FJN_S1") != nullptr;
  if (s1on) {
    double sx[NPT], sy[NPT], fx[NPT], fy[NPT], r[NPT];
    for (int p = 0; p < NPT; p++) {
      sx[p] = E.tg[0].y[p]; sy[p] = E.tg[1].y[p];
      fx[p] = sx[p] - z[0][p]; fy[p] = sy[p] - z[1][p];
    }
    S1::s1Residual(sx, sy, fx, fy, NPT, r);
    for (int p = 0; p < NPT; p++) out[p] = (float)r[p];
    return true;
  }
  if (g_nres) { for (int p = 0; p < NPT; p++) out[p] = (float)(z[0][p] * (double)g_nres_n[p] + z[1][p] * (double)g_nres_n[NPT + p]); return true; }
  for (int p = 0; p < NPT; p++) {
    const double h = std::hypot(z[0][p], z[1][p]);
    out[p] = (float)((z[0][p] + z[1][p] < 0) ? -h : h);
  }
  return true;
}

struct PostWall14 { Times& T; clk::time_point t0; ~PostWall14() { T.p_wall += ms_since(t0); } };
static void post14(Engine& E, PipeStep& PS, StepBatch& SB, int n, float* err, float* rx, float* ry, Times& T) {
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  const std::vector<int>& ct = SB.ct; const std::vector<uint8_t>& admit = *PS.admit;
  const Unique& UQ = SB.UQ; const size_t Mu = UQ.U.size();
  const int call = ++g_post_call;

  std::vector<int> K0;
  for (int i = 0; i < n; i++) {
    if (!admit[i]) continue;
    const int a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1];
    if (a < 0 || b < 0 || !SB.F[0].is_g[a] || !SB.F[1].is_g[b]) continue;
    if (!std::isfinite(SB.F[0].rfin[a]) || !std::isfinite(SB.F[1].rfin[b])) continue;
    K0.push_back(i);
  }
  if (K0.empty()) return;

  std::vector<char> need(Mu, 0), okg(Mu, 0), okp(Mu, 0);
  for (int i : K0) { need[ct[(size_t)i * 2]] = 1; need[ct[(size_t)i * 2 + 1]] = 1; }
  std::vector<uint32_t> gtok, ptok, gh(4 * Mu, 0), ph(4 * Mu, 0), gsg(36 * Mu, 0), psg(36 * Mu, 0), gpj(Mu, 0);
  std::vector<uint16_t> gsp, psp;
  for (size_t j = 0; j < Mu; j++) {
    if (!need[j]) continue;
    const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
    std::vector<uint32_t> p4; uint32_t hw = 0;
    if (zgconv14(tr, p4, gsg.data() + 36 * j, hw)) {
      okg[j] = 1;
      gh[4 * j] = (uint32_t)(gtok.size() / 4); gh[4 * j + 1] = (uint32_t)(p4.size() / 4) | ((uint32_t)tr.nS << 16);
      gh[4 * j + 2] = (uint32_t)gsp.size(); gh[4 * j + 3] = hw;
      gtok.insert(gtok.end(), p4.begin(), p4.end()); gsp.insert(gsp.end(), tr.sp.begin(), tr.sp.end());
    }
    if (te.g_ok) {
      okp[j] = 1; gpj[j] = te.gper;
      ph[4 * j] = (uint32_t)(ptok.size() / 4); ph[4 * j + 1] = (uint32_t)(te.ngp() / 4) | ((uint32_t)tr.nS << 16);
      ph[4 * j + 2] = (uint32_t)psp.size(); ph[4 * j + 3] = te.ghw;
      ptok.insert(ptok.end(), te.pgp(), te.pgp() + te.ngp()); psp.insert(psp.end(), tr.sp.begin(), tr.sp.end());
      for (int q = 0; q < 36; q++) psg[36 * j + q] = te.gg[q];
    }
  }
  for (int z = 0; z < 8; z++) { gtok.push_back(0); ptok.push_back(0); }
  std::vector<int> K;
  for (int i : K0) { const size_t a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1]; if (okg[a] && okg[b] && okp[a] && okp[b]) K.push_back(i); }
  const int nK = (int)K.size();
  T.p_cand += nK;
  if (nK == 0) return;

  std::vector<int> rowid(Mu, -1), ulist;
  for (int i : K) for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; if (rowid[j] < 0) { rowid[j] = (int)ulist.size(); ulist.push_back(j); } }
  const size_t NU = ulist.size();


  id<MTLBuffer> bTok = E.B(gtok.data(), gtok.size() * 4), bHdr = E.B(gh.data(), gh.size() * 4), bSp = E.B(gsp.data(), gsp.size() * 2 + 2), bSg = E.B(gsg.data(), gsg.size() * 4);
  id<MTLBuffer> bPtok = E.B(ptok.data(), ptok.size() * 4), bPsp = E.B(psp.data(), psp.size() * 2 + 2);
  std::vector<uint32_t> qv(Mu, 0); for (size_t j = 0; j < Mu; j++) if (okg[j]) qv[j] = (uint32_t)UQ.U[j].Q;
  id<MTLBuffer> bQ = E.B(qv.data(), Mu * 4);
  id<MTLBuffer> bGout = h900_nbl(E.dev, __LINE__, std::max<size_t>(NU * 128 * 4, 256), MTLResourceStorageModeShared);
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  if (!gps) { fprintf(stderr, "I14: grid_h26 pipeline missing\n"); exit(1); }
  struct G14 { int side; uint32_t pd[2]; int tpg; size_t ntg, tgb; id<MTLBuffer> tl; };
  std::vector<G14> GB;
  std::vector<uint32_t> oidx[2];
  auto tg0 = clk::now();
  for (int s = 0; s < 2; s++) {
    oidx[s].assign(Mu, 0);
    for (size_t r = 0; r < NU; r++) oidx[s][ulist[r]] = (uint32_t)(r * 2 + s);
    std::map<int, std::vector<uint32_t>> bins;
    std::vector<char> seen(Mu, 0);
    for (int i : K) { const int j = ct[(size_t)i * 2 + s]; if (!seen[j]) { seen[j] = 1; const int nsl = std::max(4, (UQ.U[j].nslot + 3) / 4 * 4); bins[nsl].push_back((uint32_t)j); } }
    for (auto& kv : bins) {
      const int nslb = kv.first; const size_t per = (size_t)nslb * 64 + 80;
      int tpg = std::min<int>(8, (int)(32768 / (per * 4)));
      tpg = std::min<int>(tpg, (int)gps.maxTotalThreadsPerThreadgroup / 32);
      if (tpg < 1) { fprintf(stderr, "I14: grid bin nslb %d too large\n", nslb); exit(1); }
      std::vector<uint32_t> tl; tl.push_back((uint32_t)kv.second.size()); tl.insert(tl.end(), kv.second.begin(), kv.second.end());
      G14 g; g.side = s; g.pd[0] = (uint32_t)nslb; g.pd[1] = (uint32_t)per; g.tpg = tpg;
      g.ntg = (kv.second.size() + (size_t)tpg - 1) / (size_t)tpg; g.tgb = (size_t)tpg * per * 4; g.tl = E.B(tl.data(), tl.size() * 4);
      GB.push_back(g);
    }
  }
  id<MTLBuffer> bOid[2] = {E.B(oidx[0].data(), Mu * 4), E.B(oidx[1].data(), Mu * 4)};
  id<MTLBuffer> bTh[2] = {E.B(F14.th[0].data(), F14.th[0].size() * 4), E.B(F14.th[1].data(), F14.th[1].size() * 4)};
  id<MTLBuffer> bAa[2] = {E.B(F14.a[0].data(), F14.a[0].size() * 4), E.B(F14.a[1].data(), F14.a[1].size() * 4)};
  auto grid_pre = ms_since(tg0);
  double gpu_grid = 0; double w_grid = 0;
  {
    auto t0 = clk::now();
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (auto& g : GB) {
      [en setComputePipelineState:gps];

      id<MTLBuffer> bs[10] = {bTok, bHdr, bSp, bSg, bTh[g.side], bAa[g.side], bQ, g.tl, bOid[g.side], bGout};
      for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
      [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
    }
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I14 grid command buffer failed\n"); h182_dump(cb, "I14 grid"); gpu_fail(); }
    gpu_grid = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_grid = ms_since(t0);
  }
  T.p_enc += ms_since(tw) - w_grid;
  T.p_grid += gpu_grid;


  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nTile = ((size_t)nK + 31) / 32;
  id<MTLBuffer> bRows = h900_nbl(E.dev, __LINE__, (size_t)nK * 8, MTLResourceStorageModeShared);
  { uint32_t* rw = (uint32_t*)bRows.contents;
    for (int k = 0; k < nK; k++) { const int i = K[k]; rw[2 * k] = (uint32_t)rowid[ct[(size_t)i * 2]] * 2u; rw[2 * k + 1] = (uint32_t)rowid[ct[(size_t)i * 2 + 1]] * 2u + 1u; } }
  id<MTLBuffer> bCurve = h900_nbl(E.dev, __LINE__, nTile * 4096 * 4, MTLResourceStorageModeShared);
  std::memset(bCurve.contents, 0, nTile * 4096 * 4);
  id<MTLBuffer> bStr = E.B(STR, sizeof STR), bTb = E.B(BT.base, sizeof BT.base), bTd = E.B(BT.del, sizeof BT.del), bTv = E.B(BT.valid, sizeof BT.valid);
  id<MTLBuffer> bOT = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bOS = h900_nbl(E.dev, __LINE__, (size_t)nK * 48, MTLResourceStorageModeShared);
  id<MTLBuffer> bOO = h900_nbl(E.dev, __LINE__, (size_t)nK * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bTV = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4 + 64, MTLResourceStorageModeShared);
  uint32_t Pk = (uint32_t)nK;
  double gpu_band = 0, w_band = 0;
  {
    auto t0 = clk::now();
    id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
    if (!pg || !pb || !pr) { fprintf(stderr, "I14: band pipeline missing\n"); exit(1); }
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pg];
    [en setBuffer:bGout offset:0 atIndex:0]; [en setBuffer:bRows offset:0 atIndex:1]; [en setBuffer:bCurve offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeSerial);
    [en setComputePipelineState:pb];
    [en setBuffer:bCurve offset:0 atIndex:0]; [en setBuffer:bStr offset:0 atIndex:1]; [en setBuffer:bTb offset:0 atIndex:2]; [en setBuffer:bTd offset:0 atIndex:3];
    [en setBuffer:bTv offset:0 atIndex:4]; [en setBuffer:bOT offset:0 atIndex:5]; [en setBuffer:bOS offset:0 atIndex:6]; [en setBuffer:bOO offset:0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
    [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pr];
    [en setBuffer:bOT offset:0 atIndex:0]; [en setBuffer:bOS offset:0 atIndex:1]; [en setBuffer:bTV offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I14 band command buffer failed\n"); h182_dump(cb, "I14 band"); gpu_fail(); }
    gpu_band = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_band = ms_since(t0);
  }
  T.p_band += gpu_band;


  if (g_chk_band > 0) {
    const float* G = (const float*)bGout.contents;
    const uint32_t* rw = (const uint32_t*)bRows.contents; const float* OO = (const float*)bOO.contents;
    const uint8_t* OS = (const uint8_t*)bOS.contents;
    const float* TV = (const float*)bTV.contents;
    int nchk = std::min(g_chk_band, nK);
    for (int k = 0; k < nchk; k++) {
      float aos[128];
      for (int l = 0; l < 64; l++) { aos[2 * l] = G[(size_t)rw[2 * k] * 64 + l]; aos[2 * l + 1] = G[(size_t)rw[2 * k + 1] * 64 + l]; }
      int lev[band::NROW]; double obj;
      ref_cand14(aos, (const float (*)[2])STR, lev, obj);
      int levg[band::NROW]; for (int i = 0; i < band::NROW; i++) levg[i] = OS[(size_t)i * nK + k];
      bool gl_ok = true; for (int i = 0; i < band::NROW; i++) { if (!band::corridor_ok(levg[i], i)) gl_ok = false; }
      const double cg = path_cost14(aos, (const float (*)[2])STR, levg);
      const double rel = std::fabs((double)OO[k] - obj) / std::max(std::fabs(obj), 1e-30);
      const bool tie = std::fabs(cg - obj) <= 1e-5 * std::max(std::fabs(obj), 1e-30);
      bool same = true; for (int i = 0; i < band::NROW; i++) if (levg[i] != lev[i]) same = false;
      T.p_chk_n += 1;
      T.p_chk_objmax = std::max(T.p_chk_objmax, rel);
      if (!same && !tie) T.p_chk_mis += 1;
      if (!gl_ok) T.p_chk_corr += 1;
      for (int i = 1; i < band::NROW; i++) if (TV[k * 48 + i] < TV[k * 48 + i - 1]) { T.p_chk_mono += 1; break; }
    }
  }


  std::vector<uint32_t> PQ[2], PN[2], PG[2], PH[2], PSG[2], PST[2];
  std::vector<float> PTH[2], PRS[2], PA[2];
  for (int s = 0; s < 2; s++) {
    PQ[s].resize(nK); PN[s].resize(nK); PG[s].resize(nK); PH[s].assign(4 * (size_t)nK, 0); PSG[s].assign(36 * (size_t)nK, 0);
    PST[s].assign(nK, 0); PTH[s].assign(32 * (size_t)nK, 0.0f); PRS[s].assign(2 * (size_t)nK, INFINITY); PA[s].assign(16 * (size_t)nK, 0.0f);
    for (int k = 0; k < nK; k++) {
      const int j = ct[(size_t)K[k] * 2 + s];
      PQ[s][k] = (uint32_t)UQ.U[j].Q; PN[s][k] = (uint32_t)UQ.U[j].nth; PG[s][k] = gpj[j];
      std::memcpy(&PH[s][4 * k], &ph[4 * j], 16); std::memcpy(&PSG[s][36 * k], &psg[36 * j], 36 * 4);
      std::memcpy(&PTH[s][32 * k], &F14.th[s][32 * j], 32 * 4);
      PRS[s][2 * k] = PRS[s][2 * k + 1] = SB.F[s].rfin[j];
    }
  }
  P14Buf SBf[2];
  for (int s = 0; s < 2; s++) {
    SBf[s].hdr = E.B(PH[s].data(), PH[s].size() * 4); SBf[s].segs = E.B(PSG[s].data(), PSG[s].size() * 4);
    SBf[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].rs = h900_nbb(E.dev, __LINE__, PRS[s].data(), PRS[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].a = h900_nbb(E.dev, __LINE__, PA[s].data(), PA[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].stio = E.B(PST[s].data(), PST[s].size() * 4);
    SBf[s].yv = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  }
  id<MTLBuffer> bPHtok = bPtok, bPHsp = bPsp;


  std::vector<P14Bin> binsP; long ndrop = 0;
  for (int s = 0; s < 2; s++) { std::vector<int> qs(nK); for (int k = 0; k < nK; k++) qs[k] = k; p14_bins(E, s, qs, PQ[s], PN[s], PG[s], binsP, ndrop); }
  T.p_ndrop += (double)ndrop;
  double gpu_pb = 0;
  { auto t0 = clk::now(); gpu_pb = p14_run(E, binsP, SBf, bPHtok, bPHsp, bTV, 3.0f, 0.0f); T.p_pb += gpu_pb; T.p_pb_wall += ms_since(t0); }


  std::vector<char> acc[2];
  acc[0].assign(nK, 0); acc[1].assign(nK, 0);
  {
    const float* rs0d = (const float*)SBf[0].rs.contents; const float* rs1d = (const float*)SBf[1].rs.contents;
    for (int k = 0; k < nK; k++) {
      const int jx = ct[(size_t)K[k] * 2], jy = ct[(size_t)K[k] * 2 + 1];
      const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
      const double r8x = SB.F[0].rfin[jx], r8y = SB.F[1].rfin[jy];
      if (std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8x + r8y) { acc[0][k] = 1; acc[1][k] = 1; }
    }
  }

  if (g_time14) {
    P14Buf SBt[2];
    for (int s = 0; s < 2; s++) {
      std::vector<float> rs0(2 * (size_t)nK, INFINITY), a0(16 * (size_t)nK, 0.0f);
      SBt[s].hdr = SBf[s].hdr; SBt[s].segs = SBf[s].segs; SBt[s].yv = SBf[s].yv; SBt[s].stio = SBf[s].stio;
      SBt[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
      SBt[s].rs = E.B(rs0.data(), rs0.size() * 4); SBt[s].a = E.B(a0.data(), a0.size() * 4);
    }
    auto t0 = clk::now(); const double g = p14_run(E, binsP, SBt, bPHtok, bPHsp, bTV, 0.0f, 0.0f);
    T.p_pa += g; T.p_pa_wall += ms_since(t0);
  }

  if (g_chk_fin > 0) {
    int done = 0;
    for (int s = 0; s < 2 && done < g_chk_fin; s++) {
      const float* thd = (const float*)SBf[s].th.contents; const float* rsd = (const float*)SBf[s].rs.contents;
      for (int k = 0; k < nK && done < g_chk_fin; k++) {
        if (!acc[s][k]) continue;
        const int j = ct[(size_t)K[k] * 2 + s];
        const Tree& tr = UQ.U[j]; View V(tr.rec.data());
        double th[32]; for (int q = 0; q < tr.nth; q++) th[q] = (double)thd[32 * k + q];
        double tb[NPT]; for (int p = 0; p < NPT; p++) tb[p] = (double)((const float*)bTV.contents)[48 * k + p];
        std::vector<double> Ecol;
        if (!eval_cols_t14(V, th, tb, Ecol)) { T.p_fin_bad += 1; done++; continue; }
        const double ref = vp14(Ecol, tr.Q, E.tg[s].y, (double)tr.Q * 9.0 * ldexp(1.0, -24));
        const double dev = (double)rsd[2 * k + 1];
        const double rel = std::fabs(dev - ref) / std::max(ref, 1e-30);
        T.p_fin_n += 1; T.p_fin_max = std::max(T.p_fin_max, rel);
        const bool bad = !(rel <= 1e-3 || std::fabs(dev - ref) <= 1e-5 * E.tg[s].yy);
        if (bad) T.p_fin_bad += 1;
        const double tau_v = (double)tr.Q * 9.0 * ldexp(1.0, -24);
        if (g_minratio < 10 * tau_v) { T.p_fin_lo += 1; if (bad) T.p_fin_lobad += 1; }

        done++;
      }
    }
  }

  auto to = clk::now();
  if (g_dump14.size() && !g_dumpf) g_dumpf = fopen(g_dump14.c_str(), "w");
  if (g_dumpf && call == 1) fprintf(g_dumpf, "call,child,rx8,ry8,rxf,ryf,acc0,acc1\n");
  for (int k = 0; k < nK; k++) {
    const int i = K[k];
    float fin[2], r8[2]; int ac[2];
    const float rf0 = ((const float*)SBf[0].rs.contents)[2 * k + 1], rf1 = ((const float*)SBf[1].rs.contents)[2 * k + 1];
    for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; r8[s] = SB.F[s].rfin[j]; }
    ac[0] = ac[1] = acc[0][k] ? 1 : 0;
    if (ac[0] && (!std::isfinite(rf0) || !std::isfinite(rf1))) { T.p_fb += 1; ac[0] = ac[1] = 0; }
    fin[0] = ac[0] ? rf0 : r8[0]; fin[1] = ac[1] ? rf1 : r8[1];

    if ((double)fin[0] + fin[1] > (double)r8[0] + r8[1]) T.p_viol += 1;
    if (fin[0] > r8[0]) T.p_side_up += 1;
    if (fin[1] > r8[1]) T.p_side_up += 1;
    if (ac[0]) T.p_acc0 += 1; if (ac[1]) T.p_acc1 += 1;
    rx[i] = fin[0]; ry[i] = fin[1];
    const double sum = (double)fin[0] + (double)fin[1];
    err[i] = (float)std::min(sum, (double)DEAD);
    if (g_dumpf) fprintf(g_dumpf, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%d\n", call, i, r8[0], r8[1], fin[0], fin[1], ac[0], ac[1]);
  }
  T.p_cpu += ms_since(to);
  T.p_ebmb += g_ebmb_call; g_ebmb_call = 0;
  (void)w_grid; (void)gpu_pb; (void)grid_pre;
}


static int g_post16 = 2;
static std::vector<int32_t> g_loc16;
struct Loc16Reset { std::vector<int32_t>& v; const std::vector<int>& u; ~Loc16Reset() { for (int j : u) v[(size_t)j] = -1; } };

static void post16(Engine& E, PipeStep& PS, StepBatch& SB, const std::vector<int>& cand, float* err, float* rx, float* ry, Times& T) {
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  const std::vector<int>& ct = SB.ct; const std::vector<uint8_t>& admit = *PS.admit;
  const Unique& UQ = SB.UQ;
  const int call = ++g_post_call;

  std::vector<int> K0;
  for (int i : cand) {
    if (!admit[i]) continue;
    const int a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1];
    if (a < 0 || b < 0 || !SB.F[0].is_g[a] || !SB.F[1].is_g[b]) continue;
    if (!std::isfinite(SB.F[0].rfin[a]) || !std::isfinite(SB.F[1].rfin[b])) continue;
    K0.push_back(i);
  }
  if (K0.empty()) return;

  if (g_loc16.size() < 2 * (size_t)PS.n + 2) g_loc16.assign(2 * (size_t)PS.n + 2, -1);
  std::vector<int> ulist;
  for (int i : K0) for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; if (g_loc16[(size_t)j] < 0) { g_loc16[(size_t)j] = (int)ulist.size(); ulist.push_back(j); } }
  Loc16Reset loc_reset{g_loc16, ulist};
  const size_t NL = ulist.size();

  std::vector<char> okg(NL, 0), okp(NL, 0);
  std::vector<uint32_t> gtok, ptok, gh(4 * NL, 0), ph(4 * NL, 0), gsg(36 * NL, 0), psg(36 * NL, 0), gpj(NL, 0);
  std::vector<uint16_t> gsp, psp;
  for (size_t l = 0; l < NL; l++) {
    const Tree& tr = UQ.U[(size_t)ulist[l]]; const TreeEnc& te = UQ.enc[(size_t)ulist[l]];
    std::vector<uint32_t> p4; uint32_t hw = 0;
    if (zgconv14(tr, p4, gsg.data() + 36 * l, hw)) {
      okg[l] = 1;
      gh[4 * l] = (uint32_t)(gtok.size() / 4); gh[4 * l + 1] = (uint32_t)(p4.size() / 4) | ((uint32_t)tr.nS << 16);
      gh[4 * l + 2] = (uint32_t)gsp.size(); gh[4 * l + 3] = hw;
      gtok.insert(gtok.end(), p4.begin(), p4.end()); gsp.insert(gsp.end(), tr.sp.begin(), tr.sp.end());
    }
    if (te.g_ok) {
      okp[l] = 1; gpj[l] = te.gper;
      ph[4 * l] = (uint32_t)(ptok.size() / 4); ph[4 * l + 1] = (uint32_t)(te.ngp() / 4) | ((uint32_t)tr.nS << 16);
      ph[4 * l + 2] = (uint32_t)psp.size(); ph[4 * l + 3] = te.ghw;
      ptok.insert(ptok.end(), te.pgp(), te.pgp() + te.ngp()); psp.insert(psp.end(), tr.sp.begin(), tr.sp.end());
      for (int q = 0; q < 36; q++) psg[36 * l + q] = te.gg[q];
    }
  }
  for (int z = 0; z < 8; z++) { gtok.push_back(0); ptok.push_back(0); }

  std::vector<int> K;
  for (int i : K0) {
    const size_t la = (size_t)g_loc16[(size_t)ct[(size_t)i * 2]], lb = (size_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]];
    if (okg[la] && okg[lb] && okp[la] && okp[lb]) K.push_back(i);
  }
  const int nK = (int)K.size();
  T.p_cand += nK;
  if (nK == 0) return;

  std::vector<float> thc[2], aL[2];
  for (int s = 0; s < 2; s++) {
    thc[s].assign(NL * 32, 0.0f); aL[s].assign(NL * 16, 0.0f);
    for (size_t l = 0; l < NL; l++) {
      std::memcpy(&thc[s][32 * l], &F14.th[s][32 * (size_t)ulist[l]], 32 * sizeof(float));
      std::memcpy(&aL[s][16 * l], &F14.a[s][16 * (size_t)ulist[l]], 16 * sizeof(float));
    }
  }

  id<MTLBuffer> bTok = E.B(gtok.data(), gtok.size() * 4), bHdr = E.B(gh.data(), gh.size() * 4), bSp = E.B(gsp.data(), gsp.size() * 2 + 2), bSg = E.B(gsg.data(), gsg.size() * 4);
  id<MTLBuffer> bPtok = E.B(ptok.data(), ptok.size() * 4), bPsp = E.B(psp.data(), psp.size() * 2 + 2);
  std::vector<uint32_t> qv(NL, 0); for (size_t l = 0; l < NL; l++) if (okg[l]) qv[l] = (uint32_t)UQ.U[(size_t)ulist[l]].Q;
  id<MTLBuffer> bQ = E.B(qv.data(), NL * 4);
  id<MTLBuffer> bGout = h900_nbl(E.dev, __LINE__, std::max<size_t>(NL * 128 * 4, 256), MTLResourceStorageModeShared);
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  if (!gps) { fprintf(stderr, "I16: grid_h26 pipeline missing\n"); exit(1); }
  struct G16 { int side; uint32_t pd[2]; int tpg; size_t ntg, tgb; id<MTLBuffer> tl; };
  std::vector<G16> GB;
  std::vector<uint32_t> oidx[2];
  auto tg0 = clk::now();
  for (int s = 0; s < 2; s++) {
    oidx[s].assign(NL, 0);
    for (size_t l = 0; l < NL; l++) oidx[s][l] = (uint32_t)(l * 2 + s);
    std::map<int, std::vector<uint32_t>> bins;
    std::vector<char> seen(NL, 0);
    for (int i : K) {
      const uint32_t l = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + s]];
      if (seen[l]) continue;
      seen[l] = 1;
      const int nsl = std::max(4, (UQ.U[(size_t)ulist[l]].nslot + 3) / 4 * 4);
      bins[nsl].push_back(l);
    }
    for (auto& kv : bins) {
      const int nslb = kv.first; const size_t per = (size_t)nslb * 64 + 80;
      int tpg = std::min<int>(8, (int)(32768 / (per * 4)));
      tpg = std::min<int>(tpg, (int)gps.maxTotalThreadsPerThreadgroup / 32);
      if (tpg < 1) { fprintf(stderr, "I16: grid bin nslb %d too large\n", nslb); exit(1); }
      std::vector<uint32_t> tl; tl.push_back((uint32_t)kv.second.size()); tl.insert(tl.end(), kv.second.begin(), kv.second.end());
      G16 g; g.side = s; g.pd[0] = (uint32_t)nslb; g.pd[1] = (uint32_t)per; g.tpg = tpg;
      g.ntg = (kv.second.size() + (size_t)tpg - 1) / (size_t)tpg; g.tgb = (size_t)tpg * per * 4; g.tl = E.B(tl.data(), tl.size() * 4);
      GB.push_back(g);
    }
  }
  id<MTLBuffer> bOid[2] = {E.B(oidx[0].data(), NL * 4), E.B(oidx[1].data(), NL * 4)};
  id<MTLBuffer> bTh[2] = {E.B(thc[0].data(), thc[0].size() * 4), E.B(thc[1].data(), thc[1].size() * 4)};
  id<MTLBuffer> bAa[2] = {E.B(aL[0].data(), aL[0].size() * 4), E.B(aL[1].data(), aL[1].size() * 4)};
  double gpu_grid = 0, w_grid = 0;
  {
    auto t0 = clk::now();
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (auto& g : GB) {
      [en setComputePipelineState:gps];

      id<MTLBuffer> bs[10] = {bTok, bHdr, bSp, bSg, bTh[g.side], bAa[g.side], bQ, g.tl, bOid[g.side], bGout};
      for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
      [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
    }
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I16 grid command buffer failed\n"); h182_dump(cb, "I16 grid"); gpu_fail(); }
    gpu_grid = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_grid = ms_since(t0);
  }
  T.p_enc += ms_since(tw) - w_grid;
  T.p_grid += gpu_grid;


  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nTile = ((size_t)nK + 31) / 32;
  id<MTLBuffer> bRows = h900_nbl(E.dev, __LINE__, (size_t)nK * 8, MTLResourceStorageModeShared);
  { uint32_t* rw = (uint32_t*)bRows.contents;
    for (int k = 0; k < nK; k++) { const int i = K[k]; rw[2 * k] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2]] * 2u; rw[2 * k + 1] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]] * 2u + 1u; } }
  id<MTLBuffer> bCurve = h900_nbl(E.dev, __LINE__, nTile * 4096 * 4, MTLResourceStorageModeShared);
  std::memset(bCurve.contents, 0, nTile * 4096 * 4);
  id<MTLBuffer> bStr = E.B(STR, sizeof STR), bTb = E.B(BT.base, sizeof BT.base), bTd = E.B(BT.del, sizeof BT.del), bTv = E.B(BT.valid, sizeof BT.valid);
  id<MTLBuffer> bOT = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bOS = h900_nbl(E.dev, __LINE__, (size_t)nK * 48, MTLResourceStorageModeShared);
  id<MTLBuffer> bOO = h900_nbl(E.dev, __LINE__, (size_t)nK * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bTV = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4 + 64, MTLResourceStorageModeShared);
  uint32_t Pk = (uint32_t)nK;
  double gpu_band = 0, w_band = 0;
  {
    auto t0 = clk::now();
    id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
    if (!pg || !pb || !pr) { fprintf(stderr, "I14: band pipeline missing\n"); exit(1); }
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pg];
    [en setBuffer:bGout offset:0 atIndex:0]; [en setBuffer:bRows offset:0 atIndex:1]; [en setBuffer:bCurve offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeSerial);
    [en setComputePipelineState:pb];
    [en setBuffer:bCurve offset:0 atIndex:0]; [en setBuffer:bStr offset:0 atIndex:1]; [en setBuffer:bTb offset:0 atIndex:2]; [en setBuffer:bTd offset:0 atIndex:3];
    [en setBuffer:bTv offset:0 atIndex:4]; [en setBuffer:bOT offset:0 atIndex:5]; [en setBuffer:bOS offset:0 atIndex:6]; [en setBuffer:bOO offset:0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
    [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pr];
    [en setBuffer:bOT offset:0 atIndex:0]; [en setBuffer:bOS offset:0 atIndex:1]; [en setBuffer:bTV offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I14 band command buffer failed\n"); h182_dump(cb, "I14 band"); gpu_fail(); }
    gpu_band = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_band = ms_since(t0);
  }
  T.p_band += gpu_band;


  if (g_chk_band > 0) {
    const float* G = (const float*)bGout.contents;
    const uint32_t* rw = (const uint32_t*)bRows.contents; const float* OO = (const float*)bOO.contents;
    const uint8_t* OS = (const uint8_t*)bOS.contents;
    const float* TV = (const float*)bTV.contents;
    int nchk = std::min(g_chk_band, nK);
    for (int k = 0; k < nchk; k++) {
      float aos[128];
      for (int l = 0; l < 64; l++) { aos[2 * l] = G[(size_t)rw[2 * k] * 64 + l]; aos[2 * l + 1] = G[(size_t)rw[2 * k + 1] * 64 + l]; }
      int lev[band::NROW]; double obj;
      ref_cand14(aos, (const float (*)[2])STR, lev, obj);
      int levg[band::NROW]; for (int i = 0; i < band::NROW; i++) levg[i] = OS[(size_t)i * nK + k];
      bool gl_ok = true; for (int i = 0; i < band::NROW; i++) { if (!band::corridor_ok(levg[i], i)) gl_ok = false; }
      const double cg = path_cost14(aos, (const float (*)[2])STR, levg);
      const double rel = std::fabs((double)OO[k] - obj) / std::max(std::fabs(obj), 1e-30);
      const bool tie = std::fabs(cg - obj) <= 1e-5 * std::max(std::fabs(obj), 1e-30);
      bool same = true; for (int i = 0; i < band::NROW; i++) if (levg[i] != lev[i]) same = false;
      T.p_chk_n += 1;
      T.p_chk_objmax = std::max(T.p_chk_objmax, rel);
      if (!same && !tie) T.p_chk_mis += 1;
      if (!gl_ok) T.p_chk_corr += 1;
      for (int i = 1; i < band::NROW; i++) if (TV[k * 48 + i] < TV[k * 48 + i - 1]) { T.p_chk_mono += 1; break; }
    }
  }


  std::vector<uint32_t> PQ[2], PN[2], PG[2], PH[2], PSG[2], PST[2];
  std::vector<float> PTH[2], PRS[2], PA[2];
  for (int s = 0; s < 2; s++) {
    PQ[s].resize(nK); PN[s].resize(nK); PG[s].resize(nK); PH[s].assign(4 * (size_t)nK, 0); PSG[s].assign(36 * (size_t)nK, 0);
    PST[s].assign(nK, 0); PTH[s].assign(32 * (size_t)nK, 0.0f); PRS[s].assign(2 * (size_t)nK, INFINITY); PA[s].assign(16 * (size_t)nK, 0.0f);
    for (int k = 0; k < nK; k++) {
      const int j = ct[(size_t)K[k] * 2 + s];
      const size_t l = (size_t)g_loc16[(size_t)j];
      PQ[s][k] = (uint32_t)UQ.U[(size_t)j].Q; PN[s][k] = (uint32_t)UQ.U[(size_t)j].nth; PG[s][k] = gpj[l];
      std::memcpy(&PH[s][4 * k], &ph[4 * l], 16); std::memcpy(&PSG[s][36 * k], &psg[36 * l], 36 * 4);
      std::memcpy(&PTH[s][32 * k], &thc[s][32 * l], 32 * 4);
      PRS[s][2 * k] = PRS[s][2 * k + 1] = SB.F[s].rfin[j];
    }
  }
  P14Buf SBf[2];
  for (int s = 0; s < 2; s++) {
    SBf[s].hdr = E.B(PH[s].data(), PH[s].size() * 4); SBf[s].segs = E.B(PSG[s].data(), PSG[s].size() * 4);
    SBf[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].rs = h900_nbb(E.dev, __LINE__, PRS[s].data(), PRS[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].a = h900_nbb(E.dev, __LINE__, PA[s].data(), PA[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].stio = E.B(PST[s].data(), PST[s].size() * 4);
    SBf[s].yv = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  }
  id<MTLBuffer> bPHtok = bPtok, bPHsp = bPsp;


  std::vector<P14Bin> binsP; long ndrop = 0;
  for (int s = 0; s < 2; s++) { std::vector<int> qs(nK); for (int k = 0; k < nK; k++) qs[k] = k; p14_bins(E, s, qs, PQ[s], PN[s], PG[s], binsP, ndrop); }
  T.p_ndrop += (double)ndrop;
  double gpu_pb = 0;
  { auto t0 = clk::now(); gpu_pb = p14_run(E, binsP, SBf, bPHtok, bPHsp, bTV, 3.0f, 0.0f); T.p_pb += gpu_pb; T.p_pb_wall += ms_since(t0); }


  std::vector<char> acc[2];
  acc[0].assign(nK, 0); acc[1].assign(nK, 0);
  {
    const float* rs0d = (const float*)SBf[0].rs.contents; const float* rs1d = (const float*)SBf[1].rs.contents;
    for (int k = 0; k < nK; k++) {
      const int jx = ct[(size_t)K[k] * 2], jy = ct[(size_t)K[k] * 2 + 1];
      const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
      const double r8x = SB.F[0].rfin[jx], r8y = SB.F[1].rfin[jy];
      if (std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8x + r8y) { acc[0][k] = 1; acc[1][k] = 1; }
    }
  }

  if (g_time14) {
    P14Buf SBt[2];
    for (int s = 0; s < 2; s++) {
      std::vector<float> rs0(2 * (size_t)nK, INFINITY), a0(16 * (size_t)nK, 0.0f);
      SBt[s].hdr = SBf[s].hdr; SBt[s].segs = SBf[s].segs; SBt[s].yv = SBf[s].yv; SBt[s].stio = SBf[s].stio;
      SBt[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
      SBt[s].rs = E.B(rs0.data(), rs0.size() * 4); SBt[s].a = E.B(a0.data(), a0.size() * 4);
    }
    auto t0 = clk::now(); const double g = p14_run(E, binsP, SBt, bPHtok, bPHsp, bTV, 0.0f, 0.0f);
    T.p_pa += g; T.p_pa_wall += ms_since(t0);
  }

  if (g_chk_fin > 0) {
    int done = 0;
    for (int s = 0; s < 2 && done < g_chk_fin; s++) {
      const float* thd = (const float*)SBf[s].th.contents; const float* rsd = (const float*)SBf[s].rs.contents;
      for (int k = 0; k < nK && done < g_chk_fin; k++) {
        if (!acc[s][k]) continue;
        const int j = ct[(size_t)K[k] * 2 + s];
        const Tree& tr = UQ.U[j]; View V(tr.rec.data());
        double th[32]; for (int q = 0; q < tr.nth; q++) th[q] = (double)thd[32 * k + q];
        double tb[NPT]; for (int p = 0; p < NPT; p++) tb[p] = (double)((const float*)bTV.contents)[48 * k + p];
        std::vector<double> Ecol;
        if (!eval_cols_t14(V, th, tb, Ecol)) { T.p_fin_bad += 1; done++; continue; }
        const double ref = vp14(Ecol, tr.Q, E.tg[s].y, (double)tr.Q * 9.0 * ldexp(1.0, -24));
        const double dev = (double)rsd[2 * k + 1];
        const double rel = std::fabs(dev - ref) / std::max(ref, 1e-30);
        T.p_fin_n += 1; T.p_fin_max = std::max(T.p_fin_max, rel);
        const bool bad = !(rel <= 1e-3 || std::fabs(dev - ref) <= 1e-5 * E.tg[s].yy);
        if (bad) T.p_fin_bad += 1;
        const double tau_v = (double)tr.Q * 9.0 * ldexp(1.0, -24);
        if (g_minratio < 10 * tau_v) { T.p_fin_lo += 1; if (bad) T.p_fin_lobad += 1; }

        done++;
      }
    }
  }

  auto to = clk::now();
  if (g_dump14.size() && !g_dumpf) g_dumpf = fopen(g_dump14.c_str(), "w");
  if (g_dumpf && call == 1) fprintf(g_dumpf, "call,child,rx8,ry8,rxf,ryf,acc0,acc1\n");
  for (int k = 0; k < nK; k++) {
    const int i = K[k];
    float fin[2], r8[2]; int ac[2];
    const float rf0 = ((const float*)SBf[0].rs.contents)[2 * k + 1], rf1 = ((const float*)SBf[1].rs.contents)[2 * k + 1];
    for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; r8[s] = SB.F[s].rfin[j]; }
    ac[0] = ac[1] = acc[0][k] ? 1 : 0;
    if (ac[0] && (!std::isfinite(rf0) || !std::isfinite(rf1))) { T.p_fb += 1; ac[0] = ac[1] = 0; }
    fin[0] = ac[0] ? rf0 : r8[0]; fin[1] = ac[1] ? rf1 : r8[1];

    if ((double)fin[0] + fin[1] > (double)r8[0] + r8[1]) T.p_viol += 1;
    if (fin[0] > r8[0]) T.p_side_up += 1;
    if (fin[1] > r8[1]) T.p_side_up += 1;
    if (ac[0]) T.p_acc0 += 1; if (ac[1]) T.p_acc1 += 1;
    rx[i] = fin[0]; ry[i] = fin[1];
    const double sum = (double)fin[0] + (double)fin[1];
    err[i] = (float)std::min(sum, (double)DEAD);
    if (g_dumpf) fprintf(g_dumpf, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%d\n", call, i, r8[0], r8[1], fin[0], fin[1], ac[0], ac[1]);
  }
  T.p_cpu += ms_since(to);
  T.p_ebmb += g_ebmb_call; g_ebmb_call = 0;
  (void)gpu_pb;
}


#define H53_CH 1024

static bool h900_inplace(id<MTLDevice> dev, __strong id<MTLBuffer>& b, size_t& cap, char*& ip, size_t& ipres, unsigned char& ipst, size_t bytes, const char* nm) {
  if (ipst == 2) return false;
  size_t len = al5(std::max<size_t>(bytes, 256), 16384);
  if (g_h900 & 8) len = std::max(len, cap);
  if (ipst == 0) {
    size_t R = std::min<size_t>((size_t)dev.maxBufferLength, (size_t)2 << 30); R = R / 16384 * 16384;
    void* p = mmap(nullptr, R, PROT_READ | PROT_WRITE, MAP_ANON | MAP_PRIVATE | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { ipst = 2; return false; }
    ip = (char*)p; ipres = R; ipst = 1;
    if (b && cap) std::memcpy(ip, b.contents, std::min(cap, R));
  }
  if (len > ipres) { ipst = 2; return false; }
  id<MTLBuffer> nb = [dev newBufferWithBytesNoCopy:ip length:len options:MTLResourceStorageModeShared deallocator:nil]; memdiag::metal_note(nb, len, false, true);
  if (!nb) { ipst = 2; return false; }
  if (g_h900 & 1) fprintf(stderr, "H900 inplace st=%d %s %.3f -> %.3f MB (need %.3f)\n", g_h900_st, nm, (double)cap / 1048576.0, (double)len / 1048576.0, (double)bytes / 1048576.0);
  b = nb; cap = len; return true;
}
struct H53Buf {
  id<MTLBuffer> b = nil; size_t cap = 0; char* ip = nullptr; size_t ipres = 0; unsigned char ipst = 0;
  id<MTLBuffer> get(id<MTLDevice> dev, size_t bytes) {
    bytes = std::max<size_t>(bytes, 256);
    if (!b || cap < bytes) { const size_t nc = mem_grow_cap(bytes, cap); b = h900_nbl(dev, __LINE__, nc, MSH); cap = nc; }
    return b;
  }
};
struct H53Task { std::vector<uint32_t> p4, gt, pt; std::vector<uint16_t> gs, ps; size_t og = 0, os = 0, op = 0, ops = 0; size_t l0 = 0, l1 = 0; };

struct H117G { int side = 0; uint32_t pd[2] = {0, 0}; int tpg = 1; size_t ntg = 0, tgb = 0; std::vector<uint32_t> tl; };
struct H53Scr {
  size_t h731_gb = 0, h731_kb = 0; std::vector<id<MTLCommandBuffer>> h731_cbs;
  H53Buf tok, sp, ptok, psp, hdr, phd, sg, qb, gout, oid[2], th[2], aa[2];
  id<MTLBuffer> bHdr = nil, bPH = nil, bSg = nil, bQ = nil, bOid[2] = {nil, nil}, bTh[2] = {nil, nil}, bAa[2] = {nil, nil};
  std::vector<H53Task> task; size_t ntask = 0;
  std::vector<int> K0, K, ulist;
  std::vector<char> okg, okp, seen;
  std::vector<uint32_t> psg, gpj;
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::vector<int>> kb;

  size_t nt_done = 0, k_done = 0, gw_run = 0, gs_run = 0, pw_run = 0, ps_run = 0;
  H53Buf rows, phb[2], psgb[2], pthb[2], prsb[2], pab[2], stiob[2];
  std::vector<uint32_t> PQ[2], PN[2], PG[2];
  std::vector<char> thd[2], gsd[2];
  std::vector<P14Bin> pbins;
  std::vector<H117G> gbins;
  std::vector<std::pair<size_t, size_t>> h898_rng;
  bool open = false; int call = 0;
};
static H53Scr g_h53sb[2]; static H53Scr* g_h53cur = &g_h53sb[0];
#define g_h53s (*g_h53cur)
static H53Scr* g_h53p = &g_h53sb[0];
static int g_h91 = 0;
static int g_h92 = 3;
static int g_h53 = 1;


static int g_h84 = 2;
static int h84_nb(size_t n, size_t gr) { return (g_h84 && n >= 2 * gr) ? (int)((n + gr - 1) / gr) : 1; }
static void h84_run(Engine& E, size_t n, size_t gr, const std::function<void(size_t, size_t, int)>& f) {
  const int nb = h84_nb(n, gr);
  if (nb == 1) { f(0, n, 0); return; }
  PL(E).run(nb, [&](int b) { f((size_t)b * gr, std::min(n, ((size_t)b + 1) * gr), b); });
}
static void kb_fit(H53Scr& HS, int nb) { if ((int)HS.kb.size() < nb) HS.kb.resize((size_t)nb); }


static inline bool postlin_side(const StepBatch& SB, int s, int uid) { return SB.F[s].is_g[(size_t)uid] || (postlin_on(s) && SB.UQ.U[(size_t)uid].nth == 0); }

static void post_lin_a(Engine& E, const Unique& UQ, int s, int uid) {
  const Tree& tr = UQ.U[(size_t)uid]; float* ap = &F14.a[s][16 * (size_t)uid]; for (int q = 0; q < 16; q++) ap[q] = 0.0f;
  const int Q = tr.Q; if (Q < 1 || Q > 16) return;
  static float tgf[48]; static bool ti = false; if (!ti) { for (int p = 0; p < 48; p++) tgf[p] = (float)((double)p / 47.0); ti = true; }
  thread_local std::vector<double> vb, vm; View V(tr.rec.data()); const int n = V.n(); if ((int)vm.size() < n * NPT) vm.resize((size_t)n * NPT);
  double Ev[16 * NPT];
  for (int i = 0; i < n; i++) { const int o = V.op(i), ia = V.a(i), ib = V.b(i); double* d = &vm[(size_t)i * NPT]; const double* A = &vm[(size_t)ia * NPT]; const double* B = &vm[(size_t)ib * NPT];
    switch (o) { case 1: for (int p = 0; p < NPT; p++) d[p] = 0; break; case 2: for (int p = 0; p < NPT; p++) d[p] = (double)tgf[p]; break; case 16: for (int p = 0; p < NPT; p++) d[p] = 1; break;
      case 3: for (int p = 0; p < NPT; p++) d[p] = A[p] + B[p]; break; case 4: for (int p = 0; p < NPT; p++) d[p] = A[p] - B[p]; break; case 5: for (int p = 0; p < NPT; p++) d[p] = A[p] * B[p]; break; case 6: for (int p = 0; p < NPT; p++) d[p] = A[p] / B[p]; break;
      case 7: for (int p = 0; p < NPT; p++) d[p] = sin(A[p]); break; case 8: for (int p = 0; p < NPT; p++) d[p] = cos(A[p]); break; case 9: for (int p = 0; p < NPT; p++) d[p] = exp(A[p]); break; case 10: for (int p = 0; p < NPT; p++) d[p] = log(fabs(A[p])); break;
      case 11: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p]; break; case 12: for (int p = 0; p < NPT; p++) d[p] = A[p] * A[p] * A[p]; break; case 13: for (int p = 0; p < NPT; p++) d[p] = sqrt(fabs(A[p])); break; case 14: for (int p = 0; p < NPT; p++) d[p] = cbrt(A[p]); break; case 15: for (int p = 0; p < NPT; p++) d[p] = -A[p]; break;
      default: for (int p = 0; p < NPT; p++) d[p] = 0; break; } }
  for (int q = 0; q < Q; q++) { const double* R = &vm[(size_t)V.root(q) * NPT]; for (int p = 0; p < NPT; p++) { if (!std::isfinite(R[p])) return; Ev[(size_t)q * NPT + p] = R[p]; } }
  const double* y = E.tg[s].y; double Qm[16 * NPT], Rm[16][16]; bool keep[16]; double tau = (double)Q * 9.0 * ldexp(1.0, -24); std::memcpy(Qm, Ev, sizeof(double) * (size_t)Q * NPT); for (auto& r : Rm) for (double& v : r) v = 0;
  for (int k = 0; k < Q; k++) { double on = 0; for (int p = 0; p < NPT; p++) on += Ev[(size_t)k * NPT + p] * Ev[(size_t)k * NPT + p]; on = sqrt(on);
    for (int i = 0; i < k; i++) if (keep[i]) { double sd = 0; for (int p = 0; p < NPT; p++) sd += Qm[(size_t)i * NPT + p] * Qm[(size_t)k * NPT + p]; Rm[i][k] = sd; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] -= sd * Qm[(size_t)i * NPT + p]; }
    double nn = 0; for (int p = 0; p < NPT; p++) nn += Qm[(size_t)k * NPT + p] * Qm[(size_t)k * NPT + p]; nn = sqrt(nn); keep[k] = nn > tau * on; Rm[k][k] = nn; for (int p = 0; p < NPT; p++) Qm[(size_t)k * NPT + p] = keep[k] ? Qm[(size_t)k * NPT + p] / nn : 0; }
  double c[16]; for (int k = 0; k < Q; k++) { c[k] = 0; if (keep[k]) for (int p = 0; p < NPT; p++) c[k] += Qm[(size_t)k * NPT + p] * y[p]; }
  double a[16]; for (int k = Q - 1; k >= 0; k--) { if (!keep[k]) { a[k] = 0; continue; } double v = c[k]; for (int j = k + 1; j < Q; j++) if (keep[j]) v -= Rm[k][j] * a[j]; a[k] = v / Rm[k][k]; }
  for (int k = 0; k < Q; k++) ap[k] = (float)a[k];
}
static void post16_h53_prep(Engine& E, PipeStep& PS, StepBatch& SB, const std::vector<int>& cand, Times& T) {
  if (g_h670 && !g_h675) h670_join();
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  const std::vector<int>& ct = SB.ct; const std::vector<uint8_t>& admit = *PS.admit;
  const Unique& UQ = SB.UQ;
  H53Scr& HS = *g_h53p;
  id<MTLDevice> dev = E.dev;
  const size_t CAP = g_h91 ? 2 * cand.size() + 2 : 2 * (size_t)PS.n + 2;
  if (!HS.open) {
    HS.open = true; HS.call = ++g_post_call;
    HS.K0.clear(); HS.K.clear(); HS.ulist.clear(); HS.ntask = 0;
    HS.nt_done = 0; HS.k_done = 0; HS.h898_rng.clear();
    if (g_loc16.size() < 2 * (size_t)PS.n + 2) g_loc16.assign(2 * (size_t)PS.n + 2, -1);
    HS.bHdr = HS.hdr.get(E.dev, CAP * 16); HS.bPH = HS.phd.get(dev, CAP * 16); HS.bSg = HS.sg.get(dev, CAP * 144); HS.bQ = HS.qb.get(dev, CAP * 4);
    HS.bOid[0] = HS.oid[0].get(dev, CAP * 4); HS.bOid[1] = HS.oid[1].get(dev, CAP * 4);
    HS.bTh[0] = HS.th[0].get(dev, CAP * 128); HS.bTh[1] = HS.th[1].get(dev, CAP * 128);
    HS.bAa[0] = HS.aa[0].get(dev, CAP * 64); HS.bAa[1] = HS.aa[1].get(dev, CAP * 64);
    HS.okg.resize(CAP); HS.okp.resize(CAP); HS.psg.resize(36 * CAP); HS.gpj.resize(CAP);
  }

  const size_t NC = cand.size(), GR = 4096;
  kb_fit(HS, h84_nb(NC, GR));
  h84_run(E, NC, GR, [&](size_t lo, size_t hi, int bi) {
    std::vector<int>& v = HS.kb[(size_t)bi]; v.clear();
    for (size_t q = lo; q < hi; q++) {
      const int i = cand[q];
      if (!admit[i]) continue;
      const int ua = ct[(size_t)i * 2], ub = ct[(size_t)i * 2 + 1];
      if (ua < 0 || ub < 0 || !postlin_side(SB, 0, ua) || !postlin_side(SB, 1, ub)) continue;
      if (!std::isfinite(SB.F[0].rfin[ua]) || !std::isfinite(SB.F[1].rfin[ub])) continue;
      v.push_back(i);
    }
  });
  const size_t k0b = HS.K0.size();
  for (size_t bi = 0; bi < (size_t)h84_nb(NC, GR); bi++) HS.K0.insert(HS.K0.end(), HS.kb[bi].begin(), HS.kb[bi].end());

  std::vector<int>& ulist = HS.ulist;
  const size_t L0 = ulist.size();
  for (size_t q = k0b; q < HS.K0.size(); q++) {
    const int i = HS.K0[q];
    for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; if (g_loc16[(size_t)j] < 0) { g_loc16[(size_t)j] = (int)ulist.size(); ulist.push_back(j); } }
  }
  const size_t NL = ulist.size();


  const int NCH = (int)((NL - L0 + H53_CH - 1) / H53_CH);
  const size_t T0 = HS.ntask;
  if (HS.task.size() < T0 + (size_t)NCH) HS.task.resize(T0 + (size_t)NCH);
  uint32_t* hdrd = (uint32_t*)HS.bHdr.contents; uint32_t* phd = (uint32_t*)HS.bPH.contents; uint32_t* sgd = (uint32_t*)HS.bSg.contents;
  uint32_t* qd = (uint32_t*)HS.bQ.contents;
  uint32_t* oid0 = (uint32_t*)HS.bOid[0].contents; uint32_t* oid1 = (uint32_t*)HS.bOid[1].contents;
  char* okgd = HS.okg.data(); char* okpd = HS.okp.data();
  uint32_t* psgd = HS.psg.data(); uint32_t* gpjd = HS.gpj.data();
  const int* ul = ulist.data();
  if (NCH > 0) PL(E).run(NCH, [&](int c) {
    H53Task& tk = HS.task[T0 + (size_t)c];
    tk.gt.clear(); tk.gs.clear(); tk.pt.clear(); tk.ps.clear();
    const size_t l0 = L0 + (size_t)c * H53_CH, l1 = std::min(NL, l0 + (size_t)H53_CH);
    tk.l0 = l0; tk.l1 = l1;
    if (g_h837) {
      const bool slim = g_h837 == 1 || g_h837 == 2, pf = g_h837 == 1 || g_h837 == 3;
      auto pf_lines = [](const void* p, size_t bytes) { const char* c = (const char*)p; if (!c) return; for (size_t k = 0; k < bytes; k += 64) __builtin_prefetch(c + k); };
      auto pf_struct = [&](size_t l) { const int jj = ul[l]; pf_lines(&UQ.U[(size_t)jj], sizeof(Tree)); pf_lines(&UQ.enc[(size_t)jj], sizeof(TreeEnc)); };
      auto pf_payload = [&](size_t l) { const int jj = ul[l]; const Tree& t2 = UQ.U[(size_t)jj]; const TreeEnc& e2 = UQ.enc[(size_t)jj];
        pf_lines(h313a::prog(t2), std::min<size_t>(h313a::nprog(t2) * 4, 512)); pf_lines(h313a::segs(t2), 72); pf_lines(t2.sp.data(), std::min<size_t>(t2.sp.size() * 2, 128));
        if (e2.g_ok) pf_lines(e2.pgp(), std::min<size_t>((size_t)e2.ngp() * 4, 512)); };

      if (pf) for (size_t k = 0; k < 4 && l0 + k < l1; k++) { pf_struct(l0 + k); if (k < 2) pf_payload(l0 + k); }
      for (size_t l = l0; l < l1; l++) {
        if (pf) { if (l + 4 < l1) pf_struct(l + 4); if (l + 2 < l1) pf_payload(l + 2); }
        const int j = ul[l];
        const Tree& tr = UQ.U[(size_t)j]; const TreeEnc& te = UQ.enc[(size_t)j];
        uint32_t* hd = hdrd + 4 * l; uint32_t* sgl = sgd + 36 * l; uint32_t* pl = phd + 4 * l;
        okgd[l] = 0; okpd[l] = 0; gpjd[l] = 0;
        uint32_t hw = 0;
        if (slim) {
          const size_t g0 = tk.gt.size();
          if (zgconv14d(tr, tk.gt, sgl, hw)) {
            okgd[l] = 1;
            hd[0] = (uint32_t)(g0 / 4); hd[1] = (uint32_t)((tk.gt.size() - g0) / 4) | ((uint32_t)tr.nS << 16);
            hd[2] = (uint32_t)tk.gs.size(); hd[3] = hw;
            tk.gs.insert(tk.gs.end(), tr.sp.begin(), tr.sp.end());
          } else {
            tk.gt.resize(g0); std::memset(hd, 0, 16); std::memset(sgl, 0, 36 * 4); tk.p4.clear(); zgconv14(tr, tk.p4, sgl, hw);
          }
        } else {
          std::memset(hd, 0, 16); std::memset(sgl, 0, 36 * 4);
          tk.p4.clear();
          if (zgconv14(tr, tk.p4, sgl, hw)) {
            okgd[l] = 1;
            hd[0] = (uint32_t)(tk.gt.size() / 4); hd[1] = (uint32_t)(tk.p4.size() / 4) | ((uint32_t)tr.nS << 16);
            hd[2] = (uint32_t)tk.gs.size(); hd[3] = hw;
            tk.gt.insert(tk.gt.end(), tk.p4.begin(), tk.p4.end()); tk.gs.insert(tk.gs.end(), tr.sp.begin(), tr.sp.end());
          }
        }
        if (te.g_ok) {
          okpd[l] = 1; gpjd[l] = te.gper;
          pl[0] = (uint32_t)(tk.pt.size() / 4); pl[1] = (uint32_t)(te.ngp() / 4) | ((uint32_t)tr.nS << 16);
          pl[2] = (uint32_t)tk.ps.size(); pl[3] = te.ghw;
          tk.pt.insert(tk.pt.end(), te.pgp(), te.pgp() + te.ngp()); tk.ps.insert(tk.ps.end(), tr.sp.begin(), tr.sp.end());
          for (int q = 0; q < 36; q++) psgd[36 * l + q] = te.gg[q];
        } else { std::memset(pl, 0, 16); std::memset(psgd + 36 * l, 0, 36 * 4); }
        qd[l] = okgd[l] ? (uint32_t)UQ.U[(size_t)j].Q : 0u;
        oid0[l] = (uint32_t)(l * 2); oid1[l] = (uint32_t)(l * 2 + 1);
      }
      return;
    }
    for (size_t l = l0; l < l1; l++) {
      const int j = ul[l];
      const Tree& tr = UQ.U[(size_t)j]; const TreeEnc& te = UQ.enc[(size_t)j];
      uint32_t* hd = hdrd + 4 * l; uint32_t* sgl = sgd + 36 * l; uint32_t* pl = phd + 4 * l;
      std::memset(hd, 0, 16); std::memset(sgl, 0, 36 * 4); std::memset(pl, 0, 16);
      okgd[l] = 0; okpd[l] = 0; gpjd[l] = 0; std::memset(psgd + 36 * l, 0, 36 * 4);

      uint32_t hw = 0;
      tk.p4.clear();
      if (zgconv14(tr, tk.p4, sgl, hw)) {
        okgd[l] = 1;
        hd[0] = (uint32_t)(tk.gt.size() / 4); hd[1] = (uint32_t)(tk.p4.size() / 4) | ((uint32_t)tr.nS << 16);
        hd[2] = (uint32_t)tk.gs.size(); hd[3] = hw;
        tk.gt.insert(tk.gt.end(), tk.p4.begin(), tk.p4.end()); tk.gs.insert(tk.gs.end(), tr.sp.begin(), tr.sp.end());
      }

      if (te.g_ok) {
        okpd[l] = 1; gpjd[l] = te.gper;
        pl[0] = (uint32_t)(tk.pt.size() / 4); pl[1] = (uint32_t)(te.ngp() / 4) | ((uint32_t)tr.nS << 16);
        pl[2] = (uint32_t)tk.ps.size(); pl[3] = te.ghw;
        tk.pt.insert(tk.pt.end(), te.pgp(), te.pgp() + te.ngp()); tk.ps.insert(tk.ps.end(), tr.sp.begin(), tr.sp.end());
        for (int q = 0; q < 36; q++) psgd[36 * l + q] = te.gg[q];
      }


      qd[l] = okgd[l] ? (uint32_t)UQ.U[(size_t)j].Q : 0u;
      oid0[l] = (uint32_t)(l * 2); oid1[l] = (uint32_t)(l * 2 + 1);
    }
  });
  HS.ntask = T0 + (size_t)NCH;

  const size_t NK0 = HS.K0.size() - k0b;
  kb_fit(HS, h84_nb(NK0, GR));
  h84_run(E, NK0, GR, [&](size_t lo, size_t hi, int bi) {
    std::vector<int>& v = HS.kb[(size_t)bi]; v.clear();
    for (size_t q = lo; q < hi; q++) {
      const int i = HS.K0[k0b + q];
      const size_t la = (size_t)g_loc16[(size_t)ct[(size_t)i * 2]], lb = (size_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]];
      if (okgd[la] && okgd[lb] && (okpd[la] || (postlin_on(0) && UQ.U[(size_t)ct[(size_t)i * 2]].nth == 0)) && (okpd[lb] || (postlin_on(1) && UQ.U[(size_t)ct[(size_t)i * 2 + 1]].nth == 0))) v.push_back(i);
    }
  });
  if (postlin_on(0) || postlin_on(1)) {
    static thread_local int dummy = 0; (void)dummy;
    std::vector<std::pair<int, int>> need; for (size_t q = k0b; q < HS.K0.size(); q++) { const int i = HS.K0[q]; for (int s2 = 0; s2 < 2; s2++) { const int j = ct[(size_t)i * 2 + s2]; if (UQ.U[(size_t)j].nth == 0 && postlin_on(s2)) need.push_back({s2, j}); } }
    std::sort(need.begin(), need.end()); need.erase(std::unique(need.begin(), need.end()), need.end());
    PL(E).run((int)((need.size() + 63) / 64), [&](int c) { for (size_t q = (size_t)c * 64; q < std::min(need.size(), (size_t)c * 64 + 64); q++) post_lin_a(E, UQ, need[q].first, need[q].second); });
  }
  for (size_t bi = 0; bi < (size_t)h84_nb(NK0, GR); bi++) HS.K.insert(HS.K.end(), HS.kb[bi].begin(), HS.kb[bi].end());
  T.p_enc += ms_since(tw);
}


static void post16_h53_fin(Engine& E, PipeStep& PS, StepBatch& SB, float* err, float* rx, float* ry, Times& T) {
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  H53Scr& HS = *g_h53p;
  id<MTLDevice> dev = E.dev;
  const std::vector<int>& ct = SB.ct;
  const Unique& UQ = SB.UQ;
  struct OpenGuard { bool& o; ~OpenGuard() { o = false; } } open_guard{HS.open};
  Loc16Reset loc_reset{g_loc16, HS.ulist};
  if (!HS.open || HS.K0.empty()) return;
  const size_t NL = HS.ulist.size();


  size_t gw = 0, gsn = 0, pwd = 0, psn = 0;
  for (size_t c = 0; c < HS.ntask; c++) {
    H53Task& tk = HS.task[c];
    tk.og = gw; tk.os = gsn; tk.op = pwd; tk.ops = psn;
    gw += tk.gt.size(); gsn += tk.gs.size(); pwd += tk.pt.size(); psn += tk.ps.size();
  }
  id<MTLBuffer> bTok = HS.tok.get(dev, (gw + 8) * 4), bSp = HS.sp.get(dev, gsn * 2 + 2);
  id<MTLBuffer> bPtok = HS.ptok.get(dev, (pwd + 8) * 4), bPsp = HS.psp.get(dev, psn * 2 + 2);
  uint32_t* tokd = (uint32_t*)bTok.contents; uint16_t* spd = (uint16_t*)bSp.contents;
  uint32_t* ptd = (uint32_t*)bPtok.contents; uint16_t* psd = (uint16_t*)bPsp.contents;
  uint32_t* hdrd = (uint32_t*)HS.bHdr.contents; uint32_t* phd = (uint32_t*)HS.bPH.contents;
  const char* okgd = HS.okg.data(); const char* okpd = HS.okp.data();
  if (HS.ntask > 0) PL(E).run((int)HS.ntask, [&](int c) {
    const H53Task& tk = HS.task[(size_t)c];
    if (!tk.gt.empty()) std::memcpy(tokd + tk.og, tk.gt.data(), tk.gt.size() * 4);
    if (!tk.gs.empty()) std::memcpy(spd + tk.os, tk.gs.data(), tk.gs.size() * 2);
    if (!tk.pt.empty()) std::memcpy(ptd + tk.op, tk.pt.data(), tk.pt.size() * 4);
    if (!tk.ps.empty()) std::memcpy(psd + tk.ops, tk.ps.data(), tk.ps.size() * 2);
    for (size_t l = tk.l0; l < tk.l1; l++) {
      if (okgd[l]) { hdrd[4 * l] += (uint32_t)(tk.og / 4); hdrd[4 * l + 2] += (uint32_t)tk.os; }
      if (okpd[l]) { phd[4 * l] += (uint32_t)(tk.op / 4); phd[4 * l + 2] += (uint32_t)tk.ops; }
    }
  });
  for (int z = 0; z < 8; z++) { tokd[gw + z] = 0; ptd[pwd + z] = 0; }
  { uint8_t* sb = (uint8_t*)bSp.contents; sb[gsn * 2] = 0; sb[gsn * 2 + 1] = 0;
    uint8_t* pb = (uint8_t*)bPsp.contents; pb[psn * 2] = 0; pb[psn * 2 + 1] = 0; }


  {
    float* thw[2] = {(float*)HS.bTh[0].contents, (float*)HS.bTh[1].contents};
    float* aaw[2] = {(float*)HS.bAa[0].contents, (float*)HS.bAa[1].contents};
    const int* ulv = HS.ulist.data();
    h84_run(E, NL, 4096, [&](size_t lo, size_t hi, int) {
      for (size_t l = lo; l < hi; l++) {
        const size_t j = (size_t)ulv[l];
        for (int s = 0; s < 2; s++) {
          std::memcpy(thw[s] + 32 * l, &F14.th[s][32 * j], 32 * sizeof(float));
          std::memcpy(aaw[s] + 16 * l, &F14.a[s][16 * j], 16 * sizeof(float));
        }
      }
    });
  }
  const std::vector<int>& K = HS.K; const int nK = (int)K.size();
  T.p_cand += nK;
  if (nK == 0) return;

  id<MTLBuffer> bHdr = HS.bHdr, bPH = HS.bPH, bSg = HS.bSg, bQ = HS.bQ;
  id<MTLBuffer> bOid[2] = {HS.bOid[0], HS.bOid[1]}, bTh[2] = {HS.bTh[0], HS.bTh[1]}, bAa[2] = {HS.bAa[0], HS.bAa[1]};
  id<MTLBuffer> bPHtok = bPtok, bPHsp = bPsp;
  const uint32_t* ph = (const uint32_t*)bPH.contents;
  const float* thc[2] = {(const float*)bTh[0].contents, (const float*)bTh[1].contents};
  std::vector<uint32_t>& psg = HS.psg; std::vector<uint32_t>& gpj = HS.gpj;

  id<MTLBuffer> bGout = HS.gout.get(dev, std::max<size_t>(NL * 128 * 4, 256));
  struct G16 { int side; uint32_t pd[2]; int tpg; size_t ntg, tgb; id<MTLBuffer> tl; };
  std::vector<G16> GB;
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  if (!gps) { fprintf(stderr, "I16: grid_h26 pipeline missing\n"); exit(1); }
  double gpu_grid = 0, w_grid = 0;
  for (int s = 0; s < 2; s++) {

    for (auto& v : HS.bins) v.clear();
    HS.seen.assign(NL, 0);
    for (int i : K) {
      const uint32_t l = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + s]];
      if (HS.seen[l]) continue;
      HS.seen[l] = 1;
      const int nsl = std::max(4, (UQ.U[(size_t)HS.ulist[l]].nslot + 3) / 4 * 4);
      const size_t bi = (size_t)nsl / 4;
      if (HS.bins.size() <= bi) HS.bins.resize(bi + 1);
      HS.bins[bi].push_back(l);
    }
    for (size_t bi = 0; bi < HS.bins.size(); bi++) {
      const std::vector<uint32_t>& lst = HS.bins[bi];
      if (lst.empty()) continue;
      const int nslb = (int)bi * 4; const size_t per = (size_t)nslb * 64 + 80;
      int tpg = std::min<int>(8, (int)(32768 / (per * 4)));
      tpg = std::min<int>(tpg, (int)gps.maxTotalThreadsPerThreadgroup / 32);
      if (tpg < 1) { fprintf(stderr, "I16: grid bin nslb %d too large\n", nslb); exit(1); }
      std::vector<uint32_t> tl; tl.push_back((uint32_t)lst.size()); tl.insert(tl.end(), lst.begin(), lst.end());
      G16 g; g.side = s; g.pd[0] = (uint32_t)nslb; g.pd[1] = (uint32_t)per; g.tpg = tpg;
      g.ntg = (lst.size() + (size_t)tpg - 1) / (size_t)tpg; g.tgb = (size_t)tpg * per * 4; g.tl = E.B(tl.data(), tl.size() * 4);
      GB.push_back(g);
    }
  }
  {
    auto t0 = clk::now();
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (auto& g : GB) {
      [en setComputePipelineState:gps];

      id<MTLBuffer> bs[10] = {bTok, bHdr, bSp, bSg, bTh[g.side], bAa[g.side], bQ, g.tl, bOid[g.side], bGout};
      for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
      [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
    }
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I16 grid command buffer failed\n"); h182_dump(cb, "I16 grid"); gpu_fail(); }
    gpu_grid = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_grid = ms_since(t0);
  }
  T.p_enc += ms_since(tw) - w_grid;
  T.p_grid += gpu_grid;

  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nTile = ((size_t)nK + 31) / 32;
  id<MTLBuffer> bRows = h900_nbl(E.dev, __LINE__, (size_t)nK * 8, MTLResourceStorageModeShared);
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int) {
    uint32_t* rw = (uint32_t*)bRows.contents;
    for (size_t k = lo; k < hi; k++) { const int i = K[k]; rw[2 * k] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2]] * 2u; rw[2 * k + 1] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]] * 2u + 1u; }
  });
  id<MTLBuffer> bCurve = h900_nbl(E.dev, __LINE__, nTile * 4096 * 4, MTLResourceStorageModeShared);
  std::memset(bCurve.contents, 0, nTile * 4096 * 4);
  id<MTLBuffer> bStr = E.B(STR, sizeof STR), bTb = E.B(BT.base, sizeof BT.base), bTd = E.B(BT.del, sizeof BT.del), bTv = E.B(BT.valid, sizeof BT.valid);
  id<MTLBuffer> bOT = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bOS = h900_nbl(E.dev, __LINE__, (size_t)nK * 48, MTLResourceStorageModeShared);
  id<MTLBuffer> bOO = h900_nbl(E.dev, __LINE__, (size_t)nK * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bTV = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4 + 64, MTLResourceStorageModeShared);
  uint32_t Pk = (uint32_t)nK;
  double gpu_band = 0, w_band = 0;
  {
    auto t0 = clk::now();
    id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
    if (!pg || !pb || !pr) { fprintf(stderr, "I14: band pipeline missing\n"); exit(1); }
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pg];
    [en setBuffer:bGout offset:0 atIndex:0]; [en setBuffer:bRows offset:0 atIndex:1]; [en setBuffer:bCurve offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeSerial);
    [en setComputePipelineState:pb];
    [en setBuffer:bCurve offset:0 atIndex:0]; [en setBuffer:bStr offset:0 atIndex:1]; [en setBuffer:bTb offset:0 atIndex:2]; [en setBuffer:bTd offset:0 atIndex:3];
    [en setBuffer:bTv offset:0 atIndex:4]; [en setBuffer:bOT offset:0 atIndex:5]; [en setBuffer:bOS offset:0 atIndex:6]; [en setBuffer:bOO offset:0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
    [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pr];
    [en setBuffer:bOT offset:0 atIndex:0]; [en setBuffer:bOS offset:0 atIndex:1]; [en setBuffer:bTV offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "I14 band command buffer failed\n"); h182_dump(cb, "I14 band"); gpu_fail(); }
    gpu_band = (cb.GPUEndTime - cb.GPUStartTime) * 1e3; w_band = ms_since(t0);
  }
  T.p_band += gpu_band;


  if (g_chk_band > 0) {
    const float* G = (const float*)bGout.contents;
    const uint32_t* rw = (const uint32_t*)bRows.contents; const float* OO = (const float*)bOO.contents;
    const uint8_t* OS = (const uint8_t*)bOS.contents;
    const float* TV = (const float*)bTV.contents;
    int nchk = std::min(g_chk_band, nK);
    for (int k = 0; k < nchk; k++) {
      float aos[128];
      for (int l = 0; l < 64; l++) { aos[2 * l] = G[(size_t)rw[2 * k] * 64 + l]; aos[2 * l + 1] = G[(size_t)rw[2 * k + 1] * 64 + l]; }
      int lev[band::NROW]; double obj;
      ref_cand14(aos, (const float (*)[2])STR, lev, obj);
      int levg[band::NROW]; for (int i = 0; i < band::NROW; i++) levg[i] = OS[(size_t)i * nK + k];
      bool gl_ok = true; for (int i = 0; i < band::NROW; i++) { if (!band::corridor_ok(levg[i], i)) gl_ok = false; }
      const double cg = path_cost14(aos, (const float (*)[2])STR, levg);
      const double rel = std::fabs((double)OO[k] - obj) / std::max(std::fabs(obj), 1e-30);
      const bool tie = std::fabs(cg - obj) <= 1e-5 * std::max(std::fabs(obj), 1e-30);
      bool same = true; for (int i = 0; i < band::NROW; i++) if (levg[i] != lev[i]) same = false;
      T.p_chk_n += 1;
      T.p_chk_objmax = std::max(T.p_chk_objmax, rel);
      if (!same && !tie) T.p_chk_mis += 1;
      if (!gl_ok) T.p_chk_corr += 1;
      for (int i = 1; i < band::NROW; i++) if (TV[k * 48 + i] < TV[k * 48 + i - 1]) { T.p_chk_mono += 1; break; }
    }
  }


  std::vector<uint32_t> PQ[2], PN[2], PG[2], PH[2], PSG[2], PST[2];
  std::vector<float> PTH[2], PRS[2], PA[2];
  for (int s = 0; s < 2; s++) {
    PQ[s].resize(nK); PN[s].resize(nK); PG[s].resize(nK); PH[s].assign(4 * (size_t)nK, 0); PSG[s].assign(36 * (size_t)nK, 0);
    PST[s].assign(nK, 0); PTH[s].assign(32 * (size_t)nK, 0.0f); PRS[s].assign(2 * (size_t)nK, INFINITY); PA[s].assign(16 * (size_t)nK, 0.0f);
  }
  h84_run(E, (size_t)nK, 2048, [&](size_t lo, size_t hi, int) {
    for (int s = 0; s < 2; s++) for (size_t k = lo; k < hi; k++) {
      const int j = ct[(size_t)K[k] * 2 + s];
      const size_t l = (size_t)g_loc16[(size_t)j];
      PQ[s][k] = (uint32_t)UQ.U[(size_t)j].Q; PN[s][k] = (uint32_t)UQ.U[(size_t)j].nth; PG[s][k] = gpj[l];
      std::memcpy(&PH[s][4 * k], &ph[4 * l], 16); std::memcpy(&PSG[s][36 * k], &psg[36 * l], 36 * 4);
      std::memcpy(&PTH[s][32 * k], &thc[s][32 * l], 32 * 4);
      PRS[s][2 * k] = PRS[s][2 * k + 1] = SB.F[s].rfin[j];
    }
  });
  P14Buf SBf[2];
  for (int s = 0; s < 2; s++) {
    SBf[s].hdr = E.B(PH[s].data(), PH[s].size() * 4); SBf[s].segs = E.B(PSG[s].data(), PSG[s].size() * 4);
    SBf[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].rs = h900_nbb(E.dev, __LINE__, PRS[s].data(), PRS[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].a = h900_nbb(E.dev, __LINE__, PA[s].data(), PA[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].stio = E.B(PST[s].data(), PST[s].size() * 4);
    SBf[s].yv = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  }


  std::vector<P14Bin> binsP; long ndrop = 0;
  for (int s = 0; s < 2; s++) { std::vector<int> qs(nK); for (int k = 0; k < nK; k++) qs[k] = k; p14_bins(E, s, qs, PQ[s], PN[s], PG[s], binsP, ndrop); }
  T.p_ndrop += (double)ndrop;
  double gpu_pb = 0;
  { auto t0 = clk::now(); gpu_pb = p14_run(E, binsP, SBf, bPHtok, bPHsp, bTV, 3.0f, 0.0f); T.p_pb += gpu_pb; T.p_pb_wall += ms_since(t0); }
  const float* rs0d = (const float*)SBf[0].rs.contents; const float* rs1d = (const float*)SBf[1].rs.contents;


  std::vector<char> acc[2];
  acc[0].assign(nK, 0); acc[1].assign(nK, 0);
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int) {
    for (size_t k = lo; k < hi; k++) {
      const int jx = ct[(size_t)K[k] * 2], jy = ct[(size_t)K[k] * 2 + 1];
      const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
      const double r8x = SB.F[0].rfin[jx], r8y = SB.F[1].rfin[jy];
      if (std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8x + r8y) { acc[0][k] = 1; acc[1][k] = 1; }
    }
  });

  if (g_time14) {
    P14Buf SBt[2];
    for (int s = 0; s < 2; s++) {
      std::vector<float> rs0(2 * (size_t)nK, INFINITY), a0(16 * (size_t)nK, 0.0f);
      SBt[s].hdr = SBf[s].hdr; SBt[s].segs = SBf[s].segs; SBt[s].yv = SBf[s].yv; SBt[s].stio = SBf[s].stio;
      SBt[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
      SBt[s].rs = E.B(rs0.data(), rs0.size() * 4); SBt[s].a = E.B(a0.data(), a0.size() * 4);
    }
    auto t0 = clk::now(); const double g = p14_run(E, binsP, SBt, bPHtok, bPHsp, bTV, 0.0f, 0.0f);
    T.p_pa += g; T.p_pa_wall += ms_since(t0);
  }

  if (g_chk_fin > 0) {
    int done = 0;
    for (int s = 0; s < 2 && done < g_chk_fin; s++) {
      const float* thd = (const float*)SBf[s].th.contents; const float* rsd = (const float*)SBf[s].rs.contents;
      for (int k = 0; k < nK && done < g_chk_fin; k++) {
        if (!acc[s][k]) continue;
        const int j = ct[(size_t)K[k] * 2 + s];
        const Tree& tr = UQ.U[j]; View V(tr.rec.data());
        double th[32]; for (int q = 0; q < tr.nth; q++) th[q] = (double)thd[32 * k + q];
        double tb[NPT]; for (int p = 0; p < NPT; p++) tb[p] = (double)((const float*)bTV.contents)[48 * k + p];
        std::vector<double> Ecol;
        if (!eval_cols_t14(V, th, tb, Ecol)) { T.p_fin_bad += 1; done++; continue; }
        const double ref = vp14(Ecol, tr.Q, E.tg[s].y, (double)tr.Q * 9.0 * ldexp(1.0, -24));
        const double dev = (double)rsd[2 * k + 1];
        const double rel = std::fabs(dev - ref) / std::max(ref, 1e-30);
        T.p_fin_n += 1; T.p_fin_max = std::max(T.p_fin_max, rel);
        const bool bad = !(rel <= 1e-3 || std::fabs(dev - ref) <= 1e-5 * E.tg[s].yy);
        if (bad) T.p_fin_bad += 1;
        const double tau_v = (double)tr.Q * 9.0 * ldexp(1.0, -24);
        if (g_minratio < 10 * tau_v) { T.p_fin_lo += 1; if (bad) T.p_fin_lobad += 1; }
        done++;
      }
    }
  }


  auto to = clk::now();
  if (g_dump14.size() && !g_dumpf) g_dumpf = fopen(g_dump14.c_str(), "w");
  if (g_dumpf && HS.call == 1) fprintf(g_dumpf, "call,child,rx8,ry8,rxf,ryf,acc0,acc1\n");
  struct OutCnt { double viol = 0, up = 0, a0 = 0, a1 = 0, fb = 0; };
  std::vector<OutCnt> oc((size_t)h84_nb((size_t)nK, 4096));
  const bool want_dump = g_dumpf != nullptr;
  std::vector<float> ORW(want_dump ? (size_t)nK * 4 : 0);
  std::vector<uint8_t> OAC(want_dump ? (size_t)nK * 2 : 0);
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int bi) {
    OutCnt& cn = oc[(size_t)bi];
    for (size_t k = lo; k < hi; k++) {
      const int i = K[k];
      float fin[2], r8[2]; int ac[2];
      const float rf0 = rs0d[2 * k + 1], rf1 = rs1d[2 * k + 1];
      for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; r8[s] = SB.F[s].rfin[j]; }
      ac[0] = ac[1] = acc[0][k] ? 1 : 0;
      if (ac[0] && (!std::isfinite(rf0) || !std::isfinite(rf1))) { cn.fb += 1; ac[0] = ac[1] = 0; }
      fin[0] = ac[0] ? rf0 : r8[0]; fin[1] = ac[1] ? rf1 : r8[1];

      if ((double)fin[0] + fin[1] > (double)r8[0] + r8[1]) cn.viol += 1;
      if (fin[0] > r8[0]) cn.up += 1;
      if (fin[1] > r8[1]) cn.up += 1;
      if (ac[0]) cn.a0 += 1; if (ac[1]) cn.a1 += 1;
      rx[i] = fin[0]; ry[i] = fin[1];
      const double sum = (double)fin[0] + (double)fin[1];
      err[i] = (float)std::min(sum, (double)DEAD);
      if (want_dump) { ORW[4 * k] = r8[0]; ORW[4 * k + 1] = r8[1]; ORW[4 * k + 2] = fin[0]; ORW[4 * k + 3] = fin[1]; OAC[2 * k] = (uint8_t)ac[0]; OAC[2 * k + 1] = (uint8_t)ac[1]; }
    }
  });
  for (const OutCnt& cn : oc) { T.p_viol += cn.viol; T.p_side_up += cn.up; T.p_acc0 += cn.a0; T.p_acc1 += cn.a1; T.p_fb += cn.fb; }
  if (g_dumpf) for (int k = 0; k < nK; k++)
    fprintf(g_dumpf, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%d\n", HS.call, K[k], ORW[4 * k], ORW[4 * k + 1], ORW[4 * k + 2], ORW[4 * k + 3], OAC[2 * k], OAC[2 * k + 1]);
  T.p_cpu += ms_since(to);
  T.p_ebmb += g_ebmb_call; g_ebmb_call = 0;
  (void)gpu_pb;
}


struct H91Slot {
  H53Scr hs; GBuf pe, pc;
  id<MTLCommandBuffer> cb = nil; bool pend = false;
  std::shared_ptr<FitJob242> h299j; double h299_ms = 0; size_t h299_n = 0;
  id<MTLCommandBuffer> cb406 = nil; size_t n406g = 0;
  std::vector<int> K; P14Buf SBf[2]; id<MTLBuffer> bTV = nil;
};
static std::vector<H91Slot> g_h91s;
static const int g_h414r = getenv("H414R") ? std::max(1, atoi(getenv("H414R"))) : 0;
static const int g_h406 = getenv("H406") ? atoi(getenv("H406")) : 0;
static const double g_h406f = getenv("H406F") ? atof(getenv("H406F")) : -1.0;
static double g_h406_rg = 0, g_h406_rc = 0;
struct H240bItem { int side, Q; float floor_; uint32_t tree, per; };
static void h406_launch(Engine& E, std::shared_ptr<std::vector<H240bItem>> items, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp,
                        id<MTLBuffer> tvec, id<MTLCommandBuffer> cb, H91Slot& SL);

static const int g_h299 = getenv("H299") ? atoi(getenv("H299")) : 0;
struct H299Queue { std::mutex m; std::condition_variable cv; std::deque<std::pair<std::function<void()>, std::shared_ptr<FitJob242>>> q; bool started = false; };
static H299Queue g_h299q;
static std::shared_ptr<FitJob242> h299_submit(std::function<void()> f) {
  auto j = std::make_shared<FitJob242>();
  {
    std::lock_guard<std::mutex> lk(g_h299q.m);
    if (!g_h299q.started) {
      g_h299q.started = true;
      std::thread([] {
        for (;;) {
          std::pair<std::function<void()>, std::shared_ptr<FitJob242>> it;
          { std::unique_lock<std::mutex> lk2(g_h299q.m); g_h299q.cv.wait(lk2, [] { return !g_h299q.q.empty(); }); it = std::move(g_h299q.q.front()); g_h299q.q.pop_front(); }
          it.first();
          { std::lock_guard<std::mutex> l3(it.second->m); it.second->done = true; }
          it.second->cv.notify_all();
        }
      }).detach();
    }
    g_h299q.q.push_back({std::move(f), j});
  }
  g_h299q.cv.notify_one();
  return j;
}
static void h299_join(std::shared_ptr<FitJob242>& j) {
  if (!j) return;
  std::unique_lock<std::mutex> lk(j->m); j->cv.wait(lk, [&] { return j->done; }); lk.unlock(); j.reset();
}

static void h299_launch(Engine& E, const std::vector<P14Bin>& bins, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp,
                        id<MTLBuffer> tvec, id<MTLCommandBuffer> cb, H91Slot& SL);
static id<MTLCommandQueue> g_q91 = nil;


static size_t g_h894_need = 0, g_h894_last_eo = 0; static id<MTLCommandQueue> g_h894q = nil; static id<MTLCommandBuffer> g_h894cb = nil;
static void h894_join() { if (g_h894cb) { [g_h894cb waitUntilCompleted]; g_h894cb = nil; } }
static void h894_pregrow(id<MTLDevice> dev, GBuf& pe) {
  if (!(g_h894 & 1) || g_h894_need == 0) return;
  h894_join();
  const size_t want = g_h894_need + g_h894_need / 2;
  if (pe.b && pe.cap >= want) return;
  pe.b = h900_nbl(dev, __LINE__, want, MPV); pe.cap = want;
  if (!g_h894q) g_h894q = [dev newCommandQueue];
  id<MTLCommandBuffer> cb = [g_h894q commandBuffer];
  id<MTLBlitCommandEncoder> bl = [cb blitCommandEncoder];
  [bl fillBuffer:pe.b range:NSMakeRange(0, pe.cap) value:0];
  [bl endEncoding]; [cb commit]; g_h894cb = cb;
  if (g_h894 & 16) fprintf(stderr, "H894 pregrow: need %.1f MB -> cap %.1f MB\n", (double)g_h894_need / 1048576.0, (double)pe.cap / 1048576.0);
}

static const H91Slot* g_h899_fin = nullptr;
struct H899Pe { id<MTLHeap> heap = nil; size_t hsize = 0; id<MTLCommandBuffer> last = nil; };
static H899Pe g_h899_pe;
static void h899_pe_ensure(Engine& E, GBuf& pe, size_t eo, id<MTLCommandBuffer> cb) {
  eo = al5(std::max<size_t>(eo, 256), 16384);
  H899Pe& H = g_h899_pe;
  if (pe.b && pe.cap >= eo && !(g_h899 & 4)) { H.last = cb; return; }
  if (H.last) { const MTLCommandBufferStatus st = H.last.status;
    if (st >= MTLCommandBufferStatusCommitted && st < MTLCommandBufferStatusCompleted) [H.last waitUntilCompleted]; }
  pe.b = nil; pe.cap = 0;
  const char* how = "fits";
  if (!H.heap || H.hsize < eo) {
    H.heap = nil; H.last = nil;
    MTLHeapDescriptor* hd = [MTLHeapDescriptor new]; hd.size = al5((g_h930_memfix & 1) ? eo : 2 * eo, 16384); hd.storageMode = MTLStorageModePrivate; hd.hazardTrackingMode = MTLHazardTrackingModeTracked;
    H.heap = [E.dev newHeapWithDescriptor:hd]; H.hsize = hd.size; how = "new heap";
    if (!H.heap) { fprintf(stderr, "H899: heap of %zu B failed, using GBuf\n", H.hsize); H.hsize = 0; h900_ens(pe, "pe", E.dev, eo, MPV); H.last = cb; return; }
  } else how = "reuse";
  pe.b = [H.heap newBufferWithLength:eo options:MTLResourceStorageModePrivate | MTLResourceHazardTrackingModeTracked]; pe.cap = eo; memdiag::metal_note(pe.b, eo, true);
  if (!pe.b) { fprintf(stderr, "H899: heap buffer of %zu B failed (heap %zu), using GBuf\n", eo, H.hsize); pe.cap = 0; h900_ens(pe, "pe", E.dev, eo, MPV); }
  if (g_h899 & 1) fprintf(stderr, "H899 pe: need %.3f MB buf %.3f MB heap %.3f MB (%s)\n", (double)eo / 1048576.0, (double)pe.cap / 1048576.0, (double)H.hsize / 1048576.0, how);
  H.last = cb;
}
static void p14_encode_h91(Engine& E, id<MTLCommandBuffer> cb, std::vector<P14Bin>& bins, const P14Buf* SBf, id<MTLBuffer> tok,
                           id<MTLBuffer> psp, id<MTLBuffer> tvec, float itmax, float skipf, H91Slot& SL) {
  if (bins.empty()) return;

  if (h94_mode() >= 1) {
    std::vector<size_t> nat(bins.size()), tgb(bins.size()), cnt(bins.size()), cap(bins.size()); size_t before = 0, after = 0, ntr = 0;
    for (size_t b = 0; b < bins.size(); b++) {
      const P14Bin& g = bins[b]; nat[b] = g.ntg; tgb[b] = g.tgb; size_t c = 0;
      for (size_t k = 1; k < g.tl.size(); k++) if (g.tl[k] != 0xFFFFFFFFu) c++;
      cnt[b] = c; ntr += c; before += g.ntg;
    }
    h94_cap(cap, nat, tgb, cnt);
    for (size_t b = 0; b < bins.size(); b++) { P14Bin& g = bins[b]; g.ntg = cap[b]; g.eb_bytes = g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4; after += g.ntg; }
    fprintf(stderr, "H94 polish dispatch: bins %zu trees %zu TGs before %zu after %zu\n", bins.size(), ntr, before, after);
  }
  size_t eo = 0; for (size_t b = 0; b < bins.size(); b++) { bins[b].eb_off = eo; eo = al5(eo + bins[b].eb_bytes + ebt_bytes(bins[b].Q, bins[b].nth, bins[b].ntg, bins[b].tpg), 256); bins[b].wix = (uint32_t)b; }
  g_ebmb_call += (double)eo / 1048576.0;
  if (g_h894 & 1) { h894_join(); g_h894_last_eo = eo; }
  if (g_h894 & 16) { h301_mark("pol: enc ensure begin"); fprintf(stderr, "H894 pol eb: need %.1f MB cap %.1f MB bins %zu\n", (double)eo / 1048576.0, (double)SL.pe.cap / 1048576.0, bins.size()); }
  if ((g_h899 & 2) && &SL == g_h899_fin) h899_pe_ensure(E, SL.pe, eo, cb); else h900_ens(SL.pe, "SL.pe", E.dev, eo, MPV);
  h900_ens(SL.pc, "SL.pc", E.dev, 4 * bins.size() + 16, MSH);
  if (g_h894 & 16) h301_mark("pol: enc ensure end");
  id<MTLBuffer> ebB = SL.pe.b, ctr = SL.pc.b;
  std::memset(ctr.contents, 0, 4 * bins.size() + 16);
  id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
  for (size_t bi = 0; bi < bins.size(); bi++) { auto& g = bins[bi];
    const P14Buf& S = SBf[g.side];
    std::vector<uint32_t> tl334;
    if (g_h334) { tl334.reserve(g.tl.size()); tl334.push_back(g.tl[0]); for (size_t k = 1; k < g.tl.size(); k++) if (g.tl[k] != 0xFFFFFFFFu) tl334.push_back(g.tl[k]); if (tl334.size() == 1) continue; }
    id<MTLBuffer> tl = g_h334 ? E.B(tl334.data(), tl334.size() * 4) : E.B(g.tl.data(), g.tl.size() * 4);
    const size_t nt334 = g_h334 ? tl334.size() - 1 : g.nt;
    const size_t ntg334 = g_h334 ? std::min(g.ntg, (nt334 + (size_t)g.tpg - 1) / (size_t)g.tpg) : g.ntg;
    [en setComputePipelineState:g.ps];
    id<MTLBuffer> bs[15] = {tok, S.hdr, psp, S.segs, S.yv, S.th, S.rs, S.stio, tl, nil, S.a, ebB, nil, ctr, tvec};
    for (int k = 0; k < 15; k++) if (bs[k]) [en setBuffer:bs[k] offset:(k == 11 ? g.eb_off : 0) atIndex:k];
    if (g_h711) [en setBuffer:ebB offset:g.eb_off + g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4 atIndex:15];
    float tf[4] = {(float)g.Q * 9.0f * ldexpf(1.0f, -24), E.tg[g.side].floor_, itmax, skipf};
    [en setBytes:tf length:16 atIndex:9];
    uint32_t wq[2] = {g.wix, (uint32_t)nt334}; [en setBytes:wq length:8 atIndex:12];
    [en setThreadgroupMemoryLength:g.tgb atIndex:0];
    if (h150c_on()) { g_h150.polDisp++; g_h150.polLaunch += (long)(g.ntg * (size_t)g.tpg); g_h150.polOcc += (long)h150_count_tl(g.tl); }
    [en dispatchThreadgroups:MTLSizeMake(ntg334, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
  }
  [en endEncoding];
}


static void post16_h91_enc(Engine& E, PipeStep& PS, StepBatch& SB, H91Slot& SL, Times& T) {
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  H53Scr& HS = SL.hs; id<MTLDevice> dev = E.dev;
  const std::vector<int>& ct = SB.ct; const Unique& UQ = SB.UQ;
  struct OpenGuard { bool& o; ~OpenGuard() { o = false; } } open_guard{HS.open};
  Loc16Reset loc_reset{g_loc16, HS.ulist};
  SL.pend = false; SL.K.clear();
  if (!HS.open || HS.K0.empty()) return;
  const size_t NL = HS.ulist.size();
  size_t gw = 0, gsn = 0, pwd = 0, psn = 0;
  for (size_t c = 0; c < HS.ntask; c++) { H53Task& tk = HS.task[c]; tk.og = gw; tk.os = gsn; tk.op = pwd; tk.ops = psn;
    gw += tk.gt.size(); gsn += tk.gs.size(); pwd += tk.pt.size(); psn += tk.ps.size(); }
  id<MTLBuffer> bTok = HS.tok.get(dev, (gw + 8) * 4), bSp = HS.sp.get(dev, gsn * 2 + 2);
  id<MTLBuffer> bPtok = HS.ptok.get(dev, (pwd + 8) * 4), bPsp = HS.psp.get(dev, psn * 2 + 2);
  uint32_t* tokd = (uint32_t*)bTok.contents; uint16_t* spd = (uint16_t*)bSp.contents;
  uint32_t* ptd = (uint32_t*)bPtok.contents; uint16_t* psd = (uint16_t*)bPsp.contents;
  uint32_t* hdrd = (uint32_t*)HS.bHdr.contents; uint32_t* phd = (uint32_t*)HS.bPH.contents;
  const char* okgd = HS.okg.data(); const char* okpd = HS.okp.data();
  if (HS.ntask > 0) PL(E).run((int)HS.ntask, [&](int c) {
    const H53Task& tk = HS.task[(size_t)c];
    if (!tk.gt.empty()) std::memcpy(tokd + tk.og, tk.gt.data(), tk.gt.size() * 4);
    if (!tk.gs.empty()) std::memcpy(spd + tk.os, tk.gs.data(), tk.gs.size() * 2);
    if (!tk.pt.empty()) std::memcpy(ptd + tk.op, tk.pt.data(), tk.pt.size() * 4);
    if (!tk.ps.empty()) std::memcpy(psd + tk.ops, tk.ps.data(), tk.ps.size() * 2);
    for (size_t l = tk.l0; l < tk.l1; l++) {
      if (okgd[l]) { hdrd[4 * l] += (uint32_t)(tk.og / 4); hdrd[4 * l + 2] += (uint32_t)tk.os; }
      if (okpd[l]) { phd[4 * l] += (uint32_t)(tk.op / 4); phd[4 * l + 2] += (uint32_t)tk.ops; }
    }
  });
  for (int z = 0; z < 8; z++) { tokd[gw + z] = 0; ptd[pwd + z] = 0; }
  { uint8_t* sb = (uint8_t*)bSp.contents; sb[gsn * 2] = 0; sb[gsn * 2 + 1] = 0;
    uint8_t* pb = (uint8_t*)bPsp.contents; pb[psn * 2] = 0; pb[psn * 2 + 1] = 0; }


  {
    float* thw[2] = {(float*)HS.bTh[0].contents, (float*)HS.bTh[1].contents};
    float* aaw[2] = {(float*)HS.bAa[0].contents, (float*)HS.bAa[1].contents};
    const int* ulv = HS.ulist.data();
    h84_run(E, NL, 4096, [&](size_t lo, size_t hi, int) {
      for (size_t l = lo; l < hi; l++) {
        const size_t j = (size_t)ulv[l];
        for (int s = 0; s < 2; s++) {
          std::memcpy(thw[s] + 32 * l, &F14.th[s][32 * j], 32 * sizeof(float));
          std::memcpy(aaw[s] + 16 * l, &F14.a[s][16 * j], 16 * sizeof(float));
        }
      }
    });
  }
  const std::vector<int>& K = HS.K; const int nK = (int)K.size();
  T.p_cand += nK;
  if (nK == 0) return;
  SL.K = K;
  id<MTLBuffer> bHdr = HS.bHdr, bPH = HS.bPH, bSg = HS.bSg, bQ = HS.bQ;
  id<MTLBuffer> bOid[2] = {HS.bOid[0], HS.bOid[1]}, bTh[2] = {HS.bTh[0], HS.bTh[1]}, bAa[2] = {HS.bAa[0], HS.bAa[1]};
  id<MTLBuffer> bPHtok = bPtok, bPHsp = bPsp;
  const uint32_t* ph = (const uint32_t*)bPH.contents;
  const float* thc[2] = {(const float*)bTh[0].contents, (const float*)bTh[1].contents};
  std::vector<uint32_t>& psg = HS.psg; std::vector<uint32_t>& gpj = HS.gpj;
  id<MTLBuffer> bGout = HS.gout.get(dev, std::max<size_t>(NL * 128 * 4, 256));
  if (!g_q91) g_q91 = [dev newCommandQueue];
  id<MTLCommandBuffer> cb = h182_cb(g_q91);

  struct G91 { int side; uint32_t pd[2]; int tpg; size_t ntg, tgb; id<MTLBuffer> tl; };
  std::vector<G91> GB;
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  if (!gps) { fprintf(stderr, "H91: grid_h26 pipeline missing\n"); exit(1); }
  for (int s = 0; s < 2; s++) {
    for (auto& v : HS.bins) v.clear();
    HS.seen.assign(NL, 0);
    for (int i : K) {
      const uint32_t l = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + s]];
      if (HS.seen[l]) continue;
      HS.seen[l] = 1;
      const int nsl = std::max(4, (UQ.U[(size_t)HS.ulist[l]].nslot + 3) / 4 * 4);
      const size_t bi = (size_t)nsl / 4;
      if (HS.bins.size() <= bi) HS.bins.resize(bi + 1);
      HS.bins[bi].push_back(l);
    }
    for (size_t bi = 0; bi < HS.bins.size(); bi++) {
      const std::vector<uint32_t>& lst = HS.bins[bi];
      if (lst.empty()) continue;
      const int nslb = (int)bi * 4; const size_t per = (size_t)nslb * 64 + 80;
      int tpg = std::min<int>(8, (int)(32768 / (per * 4)));
      tpg = std::min<int>(tpg, (int)gps.maxTotalThreadsPerThreadgroup / 32);
      if (tpg < 1) { fprintf(stderr, "H91: grid bin nslb %d too large\n", nslb); exit(1); }
      std::vector<uint32_t> tl; tl.push_back((uint32_t)lst.size()); tl.insert(tl.end(), lst.begin(), lst.end());
      G91 g; g.side = s; g.pd[0] = (uint32_t)nslb; g.pd[1] = (uint32_t)per; g.tpg = tpg;
      g.ntg = (lst.size() + (size_t)tpg - 1) / (size_t)tpg; g.tgb = (size_t)tpg * per * 4; g.tl = E.B(tl.data(), tl.size() * 4);
      GB.push_back(g);
    }
  }
  {
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (auto& g : GB) {
      [en setComputePipelineState:gps];
      id<MTLBuffer> bs[10] = {bTok, bHdr, bSp, bSg, bTh[g.side], bAa[g.side], bQ, g.tl, bOid[g.side], bGout};
      for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
      [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
    }
    [en endEncoding];
  }

  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nTile = ((size_t)nK + 31) / 32;
  id<MTLBuffer> bRows = h900_nbl(E.dev, __LINE__, (size_t)nK * 8, MTLResourceStorageModeShared);
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int) {
    uint32_t* rw = (uint32_t*)bRows.contents;
    for (size_t k = lo; k < hi; k++) { const int i = K[k]; rw[2 * k] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2]] * 2u; rw[2 * k + 1] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]] * 2u + 1u; }
  });
  id<MTLBuffer> bCurve = h900_nbl(E.dev, __LINE__, nTile * 4096 * 4, MTLResourceStorageModeShared);
  std::memset(bCurve.contents, 0, nTile * 4096 * 4);
  id<MTLBuffer> bStr = E.B(STR, sizeof STR), bTb = E.B(BT.base, sizeof BT.base), bTd = E.B(BT.del, sizeof BT.del), bTv = E.B(BT.valid, sizeof BT.valid);
  id<MTLBuffer> bOT = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bOS = h900_nbl(E.dev, __LINE__, (size_t)nK * 48, MTLResourceStorageModeShared);
  id<MTLBuffer> bOO = h900_nbl(E.dev, __LINE__, (size_t)nK * 4, MTLResourceStorageModeShared);
  id<MTLBuffer> bTV = h900_nbl(E.dev, __LINE__, (size_t)nK * 48 * 4 + 64, MTLResourceStorageModeShared);
  uint32_t Pk = (uint32_t)nK;
  {
    id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
    if (!pg || !pb || !pr) { fprintf(stderr, "H91: band pipeline missing\n"); exit(1); }
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pg];
    [en setBuffer:bGout offset:0 atIndex:0]; [en setBuffer:bRows offset:0 atIndex:1]; [en setBuffer:bCurve offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeSerial);
    [en setComputePipelineState:pb];
    [en setBuffer:bCurve offset:0 atIndex:0]; [en setBuffer:bStr offset:0 atIndex:1]; [en setBuffer:bTb offset:0 atIndex:2]; [en setBuffer:bTd offset:0 atIndex:3];
    [en setBuffer:bTv offset:0 atIndex:4]; [en setBuffer:bOT offset:0 atIndex:5]; [en setBuffer:bOS offset:0 atIndex:6]; [en setBuffer:bOO offset:0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
    [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pr];
    [en setBuffer:bOT offset:0 atIndex:0]; [en setBuffer:bOS offset:0 atIndex:1]; [en setBuffer:bTV offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [en endEncoding];
  }

  std::vector<uint32_t> PQ[2], PN[2], PG[2], PH[2], PSG[2], PST[2];
  std::vector<float> PTH[2], PRS[2], PA[2];
  for (int s = 0; s < 2; s++) {
    PQ[s].resize(nK); PN[s].resize(nK); PG[s].resize(nK); PH[s].assign(4 * (size_t)nK, 0); PSG[s].assign(36 * (size_t)nK, 0);
    PST[s].assign(nK, 0); PTH[s].assign(32 * (size_t)nK, 0.0f); PRS[s].assign(2 * (size_t)nK, INFINITY); PA[s].assign(16 * (size_t)nK, 0.0f);
  }
  h84_run(E, (size_t)nK, 2048, [&](size_t lo, size_t hi, int) {
    for (int s = 0; s < 2; s++) for (size_t k = lo; k < hi; k++) {
      const int j = ct[(size_t)K[k] * 2 + s];
      const size_t l = (size_t)g_loc16[(size_t)j];
      PQ[s][k] = (uint32_t)UQ.U[(size_t)j].Q; PN[s][k] = (uint32_t)UQ.U[(size_t)j].nth; PG[s][k] = gpj[l];
      std::memcpy(&PH[s][4 * k], &ph[4 * l], 16); std::memcpy(&PSG[s][36 * k], &psg[36 * l], 36 * 4);
      std::memcpy(&PTH[s][32 * k], &thc[s][32 * l], 32 * 4);
      PRS[s][2 * k] = PRS[s][2 * k + 1] = SB.F[s].rfin[j];
    }
  });
  P14Buf SBf[2];
  for (int s = 0; s < 2; s++) {
    SBf[s].hdr = E.B(PH[s].data(), PH[s].size() * 4); SBf[s].segs = E.B(PSG[s].data(), PSG[s].size() * 4);
    SBf[s].th = h900_nbb(E.dev, __LINE__, PTH[s].data(), PTH[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].rs = h900_nbb(E.dev, __LINE__, PRS[s].data(), PRS[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].a = h900_nbb(E.dev, __LINE__, PA[s].data(), PA[s].size() * 4, MTLResourceStorageModeShared);
    SBf[s].stio = E.B(PST[s].data(), PST[s].size() * 4);
    SBf[s].yv = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  }

  std::vector<P14Bin> binsP; long ndrop = 0;
  for (int s = 0; s < 2; s++) { std::vector<int> qs(nK); for (int k = 0; k < nK; k++) qs[k] = k; p14_bins(E, s, qs, PQ[s], PN[s], PG[s], binsP, ndrop); }
  T.p_ndrop += (double)ndrop;
  SL.cb406 = nil; SL.n406g = 0;
  if (g_h299 && g_h406) {
    size_t ntot = 0; for (const P14Bin& g : binsP) for (size_t k = 1; k < g.tl.size(); k++) if (g.tl[k] != 0xFFFFFFFFu) ntot++;
    double f = h242d_f((g_h406_rg > 0 && g_h406_rc > 0) ? g_h406_rc / (g_h406_rc + g_h406_rg) : 0.5);
    if (g_h257f >= 0) f = g_h257f;
    if (g_h406f >= 0) f = g_h406f;
    const size_t ncpu = (size_t)llround(f * (double)ntot);
    auto items = std::make_shared<std::vector<H240bItem>>();
    for (int b = (int)binsP.size() - 1; b >= 0 && items->size() < ncpu; b--) {
      P14Bin& g = binsP[(size_t)b];
      for (size_t k = g.tl.size(); k-- > 1 && items->size() < ncpu;) {
        if (g.tl[k] == 0xFFFFFFFFu) continue;
        items->push_back(H240bItem{g.side, g.Q, E.tg[g.side].floor_, g.tl[k], g.tl[0]});
        g.tl[k] = 0xFFFFFFFFu;
      }
    }
    [cb commit];
    h406_launch(E, items, SBf, bPHtok, bPHsp, bTV, cb, SL);
    SL.n406g = ntot - items->size();
    if (SL.n406g) { id<MTLCommandBuffer> cb2 = h182_cb(g_q91); p14_encode_h91(E, cb2, binsP, SBf, bPHtok, bPHsp, bTV, 3.0f, 0.0f, SL); [cb2 commit]; SL.cb406 = cb2; }
  } else if (g_h299) {
    [cb commit];
    h299_launch(E, binsP, SBf, bPHtok, bPHsp, bTV, cb, SL);
  } else {
    p14_encode_h91(E, cb, binsP, SBf, bPHtok, bPHsp, bTV, 3.0f, 0.0f, SL);
    [cb commit];
  }
  SL.cb = cb; SL.pend = true; SL.bTV = bTV;
  for (int s = 0; s < 2; s++) SL.SBf[s] = SBf[s];
  T.p_ebmb += g_ebmb_call; g_ebmb_call = 0;
}


static void post16_h91_col(Engine& E, H91Slot& SL, StepBatch& SB, float* err, float* rx, float* ry, Times& T) {
  if (!SL.pend) return;
  auto t0 = clk::now();
  [SL.cb waitUntilCompleted];
  if (SL.cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H91 post command buffer failed\n"); h182_dump(SL.cb, "H91 post"); gpu_fail(); }
  if (SL.cb406) { [SL.cb406 waitUntilCompleted];
    if (SL.cb406.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H406 polish command buffer failed\n"); h182_dump(SL.cb406, "H406 polish"); gpu_fail(); } }
  if (SL.h299j) { h299_join(SL.h299j); T.p_pb += SL.h299_ms; }
  if (g_h406 && g_h299) {
    const double gms = SL.cb406 ? (SL.cb406.GPUEndTime - SL.cb406.GPUStartTime) * 1e3 : 0.0;
    if (gms > 0 && SL.n406g) g_h406_rg = g_h406_rg > 0 ? 0.5 * g_h406_rg + 0.5 * ((double)SL.n406g / gms) : (double)SL.n406g / gms;
    if (SL.h299_ms > 0 && SL.h299_n) g_h406_rc = g_h406_rc > 0 ? 0.5 * g_h406_rc + 0.5 * ((double)SL.h299_n / SL.h299_ms) : (double)SL.h299_n / SL.h299_ms;
    T.p_pa += gms; SL.cb406 = nil; }
  T.p_wall += ms_since(t0);
  SL.cb = nil; SL.pend = false;
  const std::vector<int>& K = SL.K; const int nK = (int)K.size();
  const std::vector<int>& ct = SB.ct;
  const float* rs0d = (const float*)SL.SBf[0].rs.contents; const float* rs1d = (const float*)SL.SBf[1].rs.contents;
  std::vector<char> acc[2];
  acc[0].assign(nK, 0); acc[1].assign(nK, 0);
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int) {
    for (size_t k = lo; k < hi; k++) {
      const int jx = ct[(size_t)K[k] * 2], jy = ct[(size_t)K[k] * 2 + 1];
      const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
      const double r8x = SB.F[0].rfin[jx], r8y = SB.F[1].rfin[jy];
      if (std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8x + r8y) { acc[0][k] = 1; acc[1][k] = 1; }
    }
  });
  struct OutCnt { double viol = 0, up = 0, a0 = 0, a1 = 0, fb = 0; };
  std::vector<OutCnt> oc((size_t)h84_nb((size_t)nK, 4096));
  h84_run(E, (size_t)nK, 4096, [&](size_t lo, size_t hi, int bi) {
    OutCnt& cn = oc[(size_t)bi];
    for (size_t k = lo; k < hi; k++) {
      const int i = K[k];
      float fin[2], r8[2]; int ac[2];
      const float rf0 = rs0d[2 * k + 1], rf1 = rs1d[2 * k + 1];
      for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; r8[s] = SB.F[s].rfin[j]; }
      ac[0] = ac[1] = acc[0][k] ? 1 : 0;
      if (ac[0] && (!std::isfinite(rf0) || !std::isfinite(rf1))) { cn.fb += 1; ac[0] = ac[1] = 0; }
      fin[0] = ac[0] ? rf0 : r8[0]; fin[1] = ac[1] ? rf1 : r8[1];
      if ((double)fin[0] + fin[1] > (double)r8[0] + r8[1]) cn.viol += 1;
      if (fin[0] > r8[0]) cn.up += 1;
      if (fin[1] > r8[1]) cn.up += 1;
      if (ac[0]) cn.a0 += 1; if (ac[1]) cn.a1 += 1;
      rx[i] = fin[0]; ry[i] = fin[1];
      const double sum = (double)fin[0] + (double)fin[1];
      err[i] = (float)std::min(sum, (double)DEAD);
    }
  });
  for (const OutCnt& cn : oc) { T.p_viol += cn.viol; T.p_side_up += cn.up; T.p_acc0 += cn.a0; T.p_acc1 += cn.a1; T.p_fb += cn.fb; }
}


static int g_h117 = 15;

static int g_h131 = 7;
struct H131Buf { id<MTLBuffer> b = nil; size_t cap = 0; };
struct H131Pool { H131Buf curve, ot, os, oo, tv, str, tb, td, tvl, yv[2]; std::vector<H131Buf> tl; };
static H131Pool g_h131p;
static id<MTLBuffer> h131_get(id<MTLDevice> dev, H131Buf& P, size_t bytes) {
  if (!P.b || P.cap < bytes) { const size_t nb = std::max<size_t>(bytes, 256); P.b = h900_nbl(dev, __LINE__, nb, MSH); P.cap = nb; }
  return P.b;
}
static id<MTLBuffer> h131_put(id<MTLDevice> dev, H131Buf& P, const void* src, size_t bytes) {
  id<MTLBuffer> b = h131_get(dev, P, bytes);
  if (bytes) std::memcpy(b.contents, src, bytes);
  return b;
}

static void h131_prealloc(id<MTLDevice> dev, size_t NP) {
  const size_t nT = (NP + 31) / 32;
  h131_get(dev, g_h131p.curve, nT * 4096 * 4); h131_get(dev, g_h131p.ot, NP * 48 * 4); h131_get(dev, g_h131p.os, NP * 48);
  h131_get(dev, g_h131p.oo, NP * 4); h131_get(dev, g_h131p.tv, NP * 48 * 4 + 64);
}

struct H4Job { double dump_ms = 0, sets_ms = 0; size_t npair = 0, ntree = 0; };
static H91Slot g_h117_sl;
static const int g_h117b = getenv("H117B") ? atoi(getenv("H117B")) : 0;
static std::vector<H91Slot> g_h117b_sl;


static void h117_reserve(H53Buf& hb, id<MTLDevice> dev, size_t bytes, size_t used) {
  bytes = std::max<size_t>(bytes, 256);
  if (hb.b && hb.cap >= bytes && !(g_h900 & 8)) return;
  if ((g_h900 & 2) && h900_inplace(dev, hb.b, hb.cap, hb.ip, hb.ipres, hb.ipst, bytes, "h117_reserve")) return;
  if (hb.b && hb.cap >= bytes) return;
  const size_t nc = mem_grow_cap(bytes, hb.cap);
  id<MTLBuffer> nb = h900_nbl(dev, __LINE__, nc, MSH);
  if (hb.b && used) std::memcpy(nb.contents, hb.b.contents, std::min(used, hb.cap));
  hb.b = nb; hb.cap = nc;
}

static void p14_bins_h898(Engine& E, int s, const std::vector<int>& qs, const std::vector<uint32_t>& pq, const std::vector<uint32_t>& pn,
                     const std::vector<uint32_t>& pg, std::vector<P14Bin>& out, long& ndrop) {
  const auto tg0 = clk::now();
  std::map<int, std::vector<int>> grp;
  { const size_t nq = qs.size(); const int nb = h84_nb(nq, 4096);
    if (nb == 1) { for (int q : qs) { if (postlin_on(s) && pn[q] == 0) continue; grp[(int)pq[q] * 25 + (int)pn[q]].push_back(q); } }
    else { std::vector<std::map<int, std::vector<int>>> part((size_t)nb);
      h84_run(E, nq, 4096, [&](size_t lo, size_t hi, int b) { auto& m = part[(size_t)b]; for (size_t i = lo; i < hi; i++) { const int q = qs[i]; if (postlin_on(s) && pn[q] == 0) continue; m[(int)pq[q] * 25 + (int)pn[q]].push_back(q); } });
      for (int b = 0; b < nb; b++) for (auto& kv : part[(size_t)b]) { auto& d = grp[kv.first]; d.insert(d.end(), kv.second.begin(), kv.second.end()); } } }
  g_h898_grp_ms += ms_since(tg0);
  for (auto& kv : grp) {
    const int Qk = kv.first / 25, Tk = kv.first % 25; size_t cnt = kv.second.size();
    id<MTLComputePipelineState> ps = E.pipe(E.Lp, "pol_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
    if (!ps) { ndrop += (long)cnt; continue; }
    {
      const int NPq = std::max(4, (Tk + 3) / 4 * 4), QPq = std::max(4, (Qk + 3) / 4 * 4);
      const size_t WSq = (4 * (size_t)NPq + 3 * (size_t)QPq + 6) / 4 * 4;
      const size_t shq = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
      std::vector<int> kp; for (int q : kv.second) { if ((shq + WSq + 2 * (size_t)pg[q]) * 4 > 32768) ndrop += 1; else kp.push_back(q); }
      kv.second.swap(kp); cnt = kv.second.size(); if (cnt == 0) continue; }
    uint32_t per = 0; for (int q : kv.second) per = std::max(per, pg[q]);
    const int NPk = std::max(4, (Tk + 3) / 4 * 4), QPk = std::max(4, (Qk + 3) / 4 * 4);
    const size_t WSk = (4 * (size_t)NPk + 3 * (size_t)QPk + 6) / 4 * 4;
    const size_t shb = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
    const size_t pers = shb + WSk + 2 * (size_t)per;
    if (pers * 4 > 32768) { ndrop += (long)cnt; continue; }
    int tpg = std::min<int>(8, (int)(32768 / (pers * 4)));
    tpg = std::min<int>(tpg, (int)ps.maxTotalThreadsPerThreadgroup / 32);
    if (tpg < 1) { ndrop += (long)cnt; continue; }
    if (h150p_on()) { const int t0 = tpg; tpg = h150_pick_pol(pers, cnt, (int)ps.maxTotalThreadsPerThreadgroup / 32, tpg); g_h150.polBins++; g_h150.polChanged += (tpg != t0); }
    size_t xf713 = 0;
    if (g_h713 && Tk >= 9) { id<MTLComputePipelineState> pp = E.pipe(E.Lpp, "polp_q" + std::to_string(Qk) + "_t" + std::to_string(Tk));
      if (pp && (int)pp.maxTotalThreadsPerThreadgroup >= 64 && (2 * pers + h713_xf(Tk)) * 4 <= 32768) { ps = pp; tpg = 2; xf713 = h713_xf(Tk); } }
    P14Bin g; g.side = s; g.Q = Qk; g.nth = Tk; g.per = per; g.tpg = tpg; g.ps = ps;
    g.nt = (cnt + (size_t)tpg - 1) / (size_t)tpg * (size_t)tpg; g.ntg = xf713 ? g.nt : g.nt / (size_t)tpg;
    g.tgb = ((size_t)tpg * pers + xf713) * 4; g.eb_bytes = g.ntg * (size_t)tpg * 2 * (size_t)(Qk + 1) * 64 * 4;
    g.tl.assign(g.nt + 1, 0xFFFFFFFFu); g.tl[0] = per;
    for (size_t k = 0; k < cnt; k++) g.tl[1 + k] = (uint32_t)kv.second[k];
    out.push_back(g);
  }
}


static void h898_polprep(Engine& E, StepBatch& SB, H53Scr& HS, Times& T, size_t k0, size_t k1, int mode) {
  const auto tp0 = clk::now();
  const std::vector<int>& ct = SB.ct; const Unique& UQ = SB.UQ; const std::vector<int>& K = HS.K;
  const uint32_t* ph = (const uint32_t*)HS.bPH.contents;
  const float* thc[2] = {(const float*)HS.bTh[0].contents, (const float*)HS.bTh[1].contents};
  if (mode != 2) h84_run(E, k1 - k0, 2048, [&](size_t lo, size_t hi, int) {
    for (int s = 0; s < 2; s++) for (size_t q = lo; q < hi; q++) {
      const size_t k = k0 + q;
      const int j = ct[(size_t)K[k] * 2 + s];
      const size_t l = (size_t)g_loc16[(size_t)j];
      HS.PQ[s][k] = (uint32_t)UQ.U[(size_t)j].Q; HS.PN[s][k] = (uint32_t)UQ.U[(size_t)j].nth; HS.PG[s][k] = HS.gpj[l];
      if (mode == 1) continue;
      uint32_t* phw = (uint32_t*)HS.phb[s].b.contents; uint32_t* psw = (uint32_t*)HS.psgb[s].b.contents;
      float* thw = (float*)HS.pthb[s].b.contents; float* rsw = (float*)HS.prsb[s].b.contents;
      float* paw = (float*)HS.pab[s].b.contents; uint32_t* stw = (uint32_t*)HS.stiob[s].b.contents;
      std::memcpy(&phw[4 * k], &ph[4 * l], 16); std::memcpy(&psw[36 * k], &HS.psg[36 * l], 36 * 4);
      std::memcpy(&thw[32 * k], &thc[s][32 * l], 32 * 4);
      rsw[2 * k] = rsw[2 * k + 1] = SB.F[s].rfin[j];
      std::memset(&paw[16 * k], 0, 16 * 4); stw[k] = 0;
    }
  });
  else h84_run(E, k1 - k0, 2048, [&](size_t lo, size_t hi, int) {
    for (int s = 0; s < 2; s++) for (size_t q = lo; q < hi; q++) {
      const size_t k = k0 + q;
      const int j = ct[(size_t)K[k] * 2 + s];
      const size_t l = (size_t)g_loc16[(size_t)j];
      uint32_t* phw = (uint32_t*)HS.phb[s].b.contents; uint32_t* psw = (uint32_t*)HS.psgb[s].b.contents;
      float* thw = (float*)HS.pthb[s].b.contents; float* rsw = (float*)HS.prsb[s].b.contents;
      float* paw = (float*)HS.pab[s].b.contents; uint32_t* stw = (uint32_t*)HS.stiob[s].b.contents;
      std::memcpy(&phw[4 * k], &ph[4 * l], 16); std::memcpy(&psw[36 * k], &HS.psg[36 * l], 36 * 4);
      std::memcpy(&thw[32 * k], &thc[s][32 * l], 32 * 4);
      rsw[2 * k] = rsw[2 * k + 1] = SB.F[s].rfin[j];
      std::memset(&paw[16 * k], 0, 16 * 4); stw[k] = 0;
    }
  });
  g_h898_arr_ms += ms_since(tp0);
  if (mode == 2) return;
  const auto tp1 = clk::now();

  long ndrop = 0;
  for (int s = 0; s < 2; s++) {
    std::vector<int> qs(k1 - k0);
    for (size_t q = 0; q < k1 - k0; q++) qs[q] = (int)(k0 + q);
    if (g_h898 & 8) p14_bins_h898(E, s, qs, HS.PQ[s], HS.PN[s], HS.PG[s], HS.pbins, ndrop); else p14_bins(E, s, qs, HS.PQ[s], HS.PN[s], HS.PG[s], HS.pbins, ndrop);
  }
  T.p_ndrop += (double)ndrop;
  g_h898_bin_ms += ms_since(tp1);
}
static void h898_flush(Engine& E, StepBatch& SB, H53Scr& HS, Times& T) {
  const auto tp = clk::now();
  for (const auto& r : HS.h898_rng) h898_polprep(E, SB, HS, T, r.first, r.second, (g_h898 & 1) ? 0 : 2);
  HS.h898_rng.clear();
  g_h898_fin_ms += ms_since(tp);
}


static void h117_cpu(Engine& E, PipeStep& PS, StepBatch& SB, H53Scr& HS, Times& T) {
  auto tw = clk::now();
  id<MTLDevice> dev = E.dev;
  const size_t n = (size_t)PS.n, CAP = HS.okg.size();
  const std::vector<int>& ct = SB.ct; const Unique& UQ = SB.UQ;
  const std::vector<int>& K = HS.K;
  if (HS.nt_done == 0 && HS.k_done == 0) {
    HS.gw_run = HS.gs_run = HS.pw_run = HS.ps_run = 0;
    for (int s = 0; s < 2; s++) { HS.thd[s].assign(CAP, 0); HS.gsd[s].assign(CAP, 0); }
    HS.pbins.clear(); HS.gbins.clear(); HS.h898_rng.clear();
    HS.h731_gb = HS.h731_kb = 0; HS.h731_cbs.clear();
  }

  const size_t t0 = HS.nt_done, t1 = HS.ntask;
  if (t1 > t0) {
    size_t gw = HS.gw_run, gsn = HS.gs_run, pwd = HS.pw_run, psn = HS.ps_run;
    for (size_t c = t0; c < t1; c++) {
      H53Task& tk = HS.task[c];
      tk.og = gw; tk.os = gsn; tk.op = pwd; tk.ops = psn;
      gw += tk.gt.size(); gsn += tk.gs.size(); pwd += tk.pt.size(); psn += tk.ps.size();
    }
    h117_reserve(HS.tok, dev, (gw + 8) * 4, HS.gw_run * 4);
    h117_reserve(HS.sp, dev, gsn * 2 + 2, HS.gs_run * 2);
    h117_reserve(HS.ptok, dev, (pwd + 8) * 4, HS.pw_run * 4);
    h117_reserve(HS.psp, dev, psn * 2 + 2, HS.ps_run * 2);
    uint32_t* tokd = (uint32_t*)HS.tok.b.contents; uint16_t* spd = (uint16_t*)HS.sp.b.contents;
    uint32_t* ptd = (uint32_t*)HS.ptok.b.contents; uint16_t* psd = (uint16_t*)HS.psp.b.contents;
    uint32_t* hdrd = (uint32_t*)HS.bHdr.contents; uint32_t* phd = (uint32_t*)HS.bPH.contents;
    const char* okgd = HS.okg.data(); const char* okpd = HS.okp.data();
    PL(E).run((int)(t1 - t0), [&](int c) {
      const H53Task& tk = HS.task[t0 + (size_t)c];
      if (!tk.gt.empty()) std::memcpy(tokd + tk.og, tk.gt.data(), tk.gt.size() * 4);
      if (!tk.gs.empty()) std::memcpy(spd + tk.os, tk.gs.data(), tk.gs.size() * 2);
      if (!tk.pt.empty()) std::memcpy(ptd + tk.op, tk.pt.data(), tk.pt.size() * 4);
      if (!tk.ps.empty()) std::memcpy(psd + tk.ops, tk.ps.data(), tk.ps.size() * 2);
      for (size_t l = tk.l0; l < tk.l1; l++) {
        if (okgd[l]) { hdrd[4 * l] += (uint32_t)(tk.og / 4); hdrd[4 * l + 2] += (uint32_t)tk.os; }
        if (okpd[l]) { phd[4 * l] += (uint32_t)(tk.op / 4); phd[4 * l + 2] += (uint32_t)tk.ops; }
      }
    });
    for (int z = 0; z < 8; z++) { tokd[gw + z] = 0; ptd[pwd + z] = 0; }
    { uint8_t* sb = (uint8_t*)HS.sp.b.contents; sb[gsn * 2] = 0; sb[gsn * 2 + 1] = 0;
      uint8_t* pb = (uint8_t*)HS.psp.b.contents; pb[psn * 2] = 0; pb[psn * 2 + 1] = 0; }
    HS.gw_run = gw; HS.gs_run = gsn; HS.pw_run = pwd; HS.ps_run = psn; HS.nt_done = t1;
  }

  const size_t k0 = HS.k_done, k1 = K.size();
  if (k1 > k0) {
    id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
    if (!gps) { fprintf(stderr, "H117: grid_h26 pipeline missing\n"); exit(1); }
    if (HS.PQ[0].size() < n) for (int s = 0; s < 2; s++) { HS.PQ[s].resize(n); HS.PN[s].resize(n); HS.PG[s].resize(n); }
    HS.rows.get(dev, n * 8);
    for (int s = 0; s < 2; s++) {
      HS.phb[s].get(dev, n * 16); HS.psgb[s].get(dev, n * 144); HS.pthb[s].get(dev, n * 128);
      HS.prsb[s].get(dev, n * 8); HS.pab[s].get(dev, n * 64); HS.stiob[s].get(dev, n * 4);
    }

    std::vector<std::pair<uint32_t, int>> cp[2];
    for (size_t k = k0; k < k1; k++) {
      const int i = K[k];
      for (int s = 0; s < 2; s++) {
        const int j = ct[(size_t)i * 2 + s]; const size_t l = (size_t)g_loc16[(size_t)j];
        if (!HS.thd[s][l]) { HS.thd[s][l] = 1; cp[s].push_back({(uint32_t)l, j}); }
      }
    }
    for (int s = 0; s < 2; s++) {
      float* thw = (float*)HS.bTh[s].contents; float* aaw = (float*)HS.bAa[s].contents;
      const std::vector<std::pair<uint32_t, int>>& v = cp[s];
      h84_run(E, v.size(), 4096, [&](size_t lo, size_t hi, int) {
        for (size_t q = lo; q < hi; q++) {
          const size_t l = v[q].first, j = (size_t)v[q].second;
          std::memcpy(thw + 32 * l, &F14.th[s][32 * j], 32 * sizeof(float));
          std::memcpy(aaw + 16 * l, &F14.a[s][16 * j], 16 * sizeof(float));
        }
      });
    }

    for (int s = 0; s < 2; s++) {
      std::vector<std::vector<uint32_t>> bl;
      for (size_t k = k0; k < k1; k++) {
        const int j = ct[(size_t)K[k] * 2 + s]; const size_t l = (size_t)g_loc16[(size_t)j];
        if (HS.gsd[s][l]) continue;
        HS.gsd[s][l] = 1;
        const int nsl = std::max(4, (UQ.U[(size_t)j].nslot + 3) / 4 * 4);
        const size_t bi = (size_t)nsl / 4;
        if (bl.size() <= bi) bl.resize(bi + 1);
        bl[bi].push_back((uint32_t)l);
      }
      for (size_t bi = 0; bi < bl.size(); bi++) {
        const std::vector<uint32_t>& lst = bl[bi];
        if (lst.empty()) continue;
        const int nslb = (int)bi * 4; const size_t per = (size_t)nslb * 64 + 80;
        int tpg = std::min<int>(8, (int)(32768 / (per * 4)));
        tpg = std::min<int>(tpg, (int)gps.maxTotalThreadsPerThreadgroup / 32);
        if (tpg < 1) { fprintf(stderr, "H117: grid bin nslb %d too large\n", nslb); exit(1); }
        H117G g; g.side = s; g.pd[0] = (uint32_t)nslb; g.pd[1] = (uint32_t)per; g.tpg = tpg;
        g.ntg = (lst.size() + (size_t)tpg - 1) / (size_t)tpg; g.tgb = (size_t)tpg * per * 4;
        g.tl.reserve(lst.size() + 1); g.tl.push_back((uint32_t)lst.size()); g.tl.insert(g.tl.end(), lst.begin(), lst.end());
        HS.gbins.push_back(std::move(g));
      }
    }

    uint32_t* rw = (uint32_t*)HS.rows.b.contents;
    h84_run(E, k1 - k0, 4096, [&](size_t lo, size_t hi, int) {
      for (size_t q = lo; q < hi; q++) {
        const size_t k = k0 + q; const int i = K[k];
        rw[2 * k] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2]] * 2u; rw[2 * k + 1] = (uint32_t)g_loc16[(size_t)ct[(size_t)i * 2 + 1]] * 2u + 1u;
      }
    });
    if (g_h898 & 1) HS.h898_rng.push_back({k0, k1});
    else if (g_h898 & 4) { HS.h898_rng.push_back({k0, k1}); h898_polprep(E, SB, HS, T, k0, k1, 1); }
    else { const auto tq = clk::now(); h898_polprep(E, SB, HS, T, k0, k1, 0); g_h898_loop_ms += ms_since(tq); }
    HS.k_done = k1;
  }
  T.p_enc += ms_since(tw);
}


#include "polish_cpu_h719.hpp"
#include "polish_cpu_h719.cpp"

struct H860Dump { std::vector<uint8_t> tok, sp, tv, hdr[2], segs[2], yv[2], th0[2]; struct It { int side, Q; uint32_t tree, per; float fl; }; std::vector<It> items; double ns = 0; };
static H860Dump* g_h860 = nullptr; static std::atomic<uint64_t> g_h860_cur{0};
static std::vector<uint8_t> h860_slurp(const std::string& p) {
  FILE* f = fopen(p.c_str(), "rb"); if (!f) { fprintf(stderr, "H860: open %s\n", p.c_str()); exit(1); }
  fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET); std::vector<uint8_t> v((size_t)n);
  if (n && fread(v.data(), 1, (size_t)n, f) != (size_t)n) { fprintf(stderr, "H860: short %s\n", p.c_str()); exit(1); } fclose(f); return v; }
static void h860_item(H860Dump& Dm, PolishWork& W, PolResult& r, size_t i) {
  const H860Dump::It& it = Dm.items[i % Dm.items.size()];
  PolInputs in{(const uint32_t*)Dm.tok.data(), (const uint32_t*)Dm.hdr[it.side].data(), (const uint16_t*)Dm.sp.data(), (const uint32_t*)Dm.segs[it.side].data(),
               (const float*)Dm.yv[it.side].data(), (const float*)Dm.tv.data(), (const float*)Dm.th0[it.side].data(), it.per};
  W.polish(in, it.tree, it.Q, (float)it.Q * 9.0f * ldexpf(1.0f, -24), it.fl, 3, r);
}
static void h860_load() {
  if (g_h860) return;
  const std::string pre = getenv("H860D") ? getenv("H860D") : "dump_pol98k/pol_s11_c1";
  H860Dump* Dm = new H860Dump();
  FILE* mf = fopen((pre + "_meta.txt").c_str(), "r"); if (!mf) { fprintf(stderr, "H860: no meta %s\n", pre.c_str()); exit(1); }
  char line[512]; if (!fgets(line, sizeof line, mf)) exit(1);
  Dm->tok = h860_slurp(pre + "_tok.bin"); Dm->sp = h860_slurp(pre + "_sp.bin"); Dm->tv = h860_slurp(pre + "_tvec.bin"); const std::vector<uint8_t> tlB = h860_slurp(pre + "_tl.bin");
  for (int s = 0; s < 2; s++) { const std::string b = pre + "_s" + std::to_string(s);
    Dm->hdr[s] = h860_slurp(b + "_hdr.bin"); Dm->segs[s] = h860_slurp(b + "_segs.bin"); Dm->yv[s] = h860_slurp(b + "_yv.bin"); Dm->th0[s] = h860_slurp(b + "_th0.bin"); }
  const uint32_t* tl = (const uint32_t*)tlB.data(); size_t off = 0;
  while (fgets(line, sizeof line, mf)) {
    int side = 0, Q = 0, tpg = 1; size_t ntg = 0, tgb = 0, nt = 0, tllen = 0; unsigned wix = 0; float fl = 0;
    if (sscanf(line, "bin side %d Q %d tpg %d ntg %zu tgb %zu nt %zu wix %u floor %f tllen %zu", &side, &Q, &tpg, &ntg, &tgb, &nt, &wix, &fl, &tllen) != 9) continue;
    const uint32_t per = tl[off];
    for (size_t k = 1; k < tllen; k++) { const uint32_t t = tl[off + k]; if (t != 0xFFFFFFFFu) Dm->items.push_back(H860Dump::It{side, Q, t, per, fl}); }
    off += tllen; }
  fclose(mf);
  { PolishWork W; PolResult r; const size_t n = std::min<size_t>(4000, Dm->items.size()); const auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < n; i++) h860_item(*Dm, W, r, i * 26);
    Dm->ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / (double)n; }
  printf("H860: real polish load from %s: %zu items, %.0f ns/item (one thread, calibration)\n", pre.c_str(), Dm->items.size(), Dm->ns);
  g_h860 = Dm;
}
static void h860_launch(int Kn) {
  h860_load();
  const uint64_t quota = (uint64_t)(g_h680cpu * 1e6 / g_h860->ns / Kn);
  const uint64_t base = g_h860_cur.fetch_add(quota);
  auto cur = std::make_shared<std::atomic<uint64_t>>(base);
  for (int t = 0; t < 8; t++) g_h680_th.emplace_back([cur, base, quota] { pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0);
    PolishWork W; PolResult r; uint64_t i;
    while ((i = cur->fetch_add(1)) < base + quota) h860_item(*g_h860, W, r, (size_t)i); });
}
static const int g_h240b = getenv("H240B") ? atoi(getenv("H240B")) : 0;
static const int g_h240bt = getenv("H240BT") ? std::max(1, atoi(getenv("H240BT"))) : 4;
static double g_h240b_rg = 0, g_h240b_rc = 0;


static const int g_h880 = getenv("H880") ? atoi(getenv("H880")) : 0;
static size_t g_h880_nmax = 0;
static void h240b_observe(size_t ntot, size_t ngpu, double gpu_ms, size_t ncpu, double cpu_ms) {
  if (g_h880 & 1) {
    const bool small = 2 * ntot < g_h880_nmax;
    if (ntot > g_h880_nmax) g_h880_nmax = ntot;
    if (small) { if (g_h880 & 2) printf("H880: skip ntot %zu nmax %zu\n", ntot, g_h880_nmax); return; }
  }
  if (gpu_ms > 0 && ngpu > 0) g_h240b_rg = g_h240b_rg > 0 ? 0.5 * g_h240b_rg + 0.5 * ((double)ngpu / gpu_ms) : (double)ngpu / gpu_ms;
  if (cpu_ms > 0 && ncpu > 0) g_h240b_rc = g_h240b_rc > 0 ? 0.5 * g_h240b_rc + 0.5 * ((double)ncpu / cpu_ms) : (double)ncpu / cpu_ms;
}
struct H240bSplit { std::vector<H240bItem> items; std::vector<std::pair<size_t, std::vector<uint32_t>>> saved; std::vector<uint64_t> mod; size_t ntot = 0; double f = 0; };

static H240bSplit h240b_mark(Engine& E, H53Scr& HS) {
  H240bSplit sp;
  for (auto& g : HS.pbins) for (size_t k = 1; k < g.tl.size(); k++) if (g.tl[k] != 0xFFFFFFFFu) sp.ntot++;
  sp.f = h242d_f((g_h240b_rg > 0 && g_h240b_rc > 0) ? g_h240b_rc / (g_h240b_rc + g_h240b_rg) : 0.5);
  if (g_h257f >= 0) sp.f = g_h257f;
  const size_t ncpu = (size_t)llround(sp.f * (double)sp.ntot);
  if (!(g_h900 & 4)) for (size_t b = 0; b < HS.pbins.size(); b++) sp.saved.push_back({b, HS.pbins[b].tl});
  for (int b = (int)HS.pbins.size() - 1; b >= 0 && sp.items.size() < ncpu; b--) {
    P14Bin& g = HS.pbins[(size_t)b];
    for (size_t k = g.tl.size(); k-- > 1 && sp.items.size() < ncpu;) {
      if (g.tl[k] == 0xFFFFFFFFu) continue;
      sp.items.push_back(H240bItem{g.side, g.Q, E.tg[g.side].floor_, g.tl[k], g.tl[0]});
      if (g_h900 & 4) sp.mod.push_back(((uint64_t)b << 32) | (uint64_t)k);
      g.tl[k] = 0xFFFFFFFFu;
    }
  }
  return sp;
}
static void h240b_restore(H53Scr& HS, const H240bSplit& sp) {
  for (const auto& p : sp.saved) HS.pbins[p.first].tl = p.second;
  for (size_t i = 0; i < sp.mod.size(); i++) HS.pbins[(size_t)(sp.mod[i] >> 32)].tl[(size_t)(sp.mod[i] & 0xFFFFFFFFu)] = sp.items[i].tree;
}


static FitHandle fw_submit_nospan265(std::function<void()> f) {
  CpuSpan265* s = g_fit_span; g_fit_span = nullptr; FitHandle h = fw_submit(std::move(f)); g_fit_span = s; return h;
}
static void h240b_cpu(const H240bSplit& sp, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp, id<MTLBuffer> tvec, double& ms_out) {
  const size_t n = sp.items.size();
  auto t0 = clk::now();
  if (n == 0) { ms_out = 0; return; }
  std::atomic<size_t> next{0};
  std::vector<std::thread> th;
  const bool lag675 = tl_h675_job;
  auto h240_wb = [&]() { H401S h401s_("polish cpu (H240B)");
    h298_util(4);
    if (lag675) pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0);
    PolishWork W;
    for (;;) {
      const size_t i0 = next.fetch_add(4);
      if (i0 >= n) break;
      h891_cpu_batch(2);
      for (size_t i = i0; i < std::min(n, i0 + 4); i++) {
        if (h650_is("polcpu")) { h650_spin_ns(g_h650ns); g_h650_items.fetch_add(1, std::memory_order_relaxed); }
        const H240bItem& it = sp.items[i];
        const P14Buf& S = SBf[it.side];
        PolInputs in{(const uint32_t*)tok.contents, (const uint32_t*)S.hdr.contents, (const uint16_t*)psp.contents,
                     (const uint32_t*)S.segs.contents, (const float*)S.yv.contents, (const float*)tvec.contents,
                     (const float*)S.th.contents, it.per};
        PolResult r;
        W.polish(in, it.tree, it.Q, (float)it.Q * 9.0f * ldexpf(1.0f, -24), it.floor_, 3, r);
        const uint32_t* hd = (const uint32_t*)S.hdr.contents + 4 * it.tree;
        const uint32_t nth = hd[3] >> 24;
        float* thp = (float*)S.th.contents + 32 * it.tree;
        for (uint32_t j = 0; j < nth; j++) thp[j] = r.th[j];
        float* rsp = (float*)S.rs.contents + 2 * it.tree;
        rsp[0] = r.rss0; rsp[1] = r.rcur;
        ((uint32_t*)S.stio.contents)[it.tree] = r.steps | (r.npp << 8) | (r.ndp << 16);
      }
    }
  };
  if (tl_in_fitw) h240_wb();
  else if (g_h873 & 1) g_h873t.run(lag675 ? 2 : 1, lag675 ? g_h675t : g_h240bt, h240_wb);
  else { const int nt = lag675 ? g_h675t : g_h240bt; for (int t = 0; t < nt; t++) th.emplace_back(h240_wb); for (auto& x : th) x.join(); }
  ms_out = ms_since(t0);
}
static void h299_launch(Engine& E, const std::vector<P14Bin>& bins, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp,
                        id<MTLBuffer> tvec, id<MTLCommandBuffer> cb, H91Slot& SL) {
  auto sp = std::make_shared<H240bSplit>();
  for (const P14Bin& g : bins) for (size_t k = 1; k < g.tl.size(); k++) {
    if (g.tl[k] == 0xFFFFFFFFu) continue;
    sp->items.push_back(H240bItem{g.side, g.Q, E.tg[g.side].floor_, g.tl[k], g.tl[0]});
  }
  sp->ntot = sp->items.size(); sp->f = 1.0;
  auto bufs = std::make_shared<std::array<P14Buf, 2>>(); (*bufs)[0] = SBf[0]; (*bufs)[1] = SBf[1];
  SL.h299_n = sp->items.size(); SL.h299_ms = 0;
  double* msp = &SL.h299_ms;
  SL.h299j = h299_submit([sp, bufs, tok, psp, tvec, cb, msp]() {
    [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H299 grid/band command buffer failed\n"); gpu_fail(); }
    double ms = 0;
    h240b_cpu(*sp, bufs->data(), tok, psp, tvec, ms);
    *msp = ms;
  });
}

static void h406_launch(Engine& E, std::shared_ptr<std::vector<H240bItem>> items, const P14Buf* SBf, id<MTLBuffer> tok, id<MTLBuffer> psp,
                        id<MTLBuffer> tvec, id<MTLCommandBuffer> cb, H91Slot& SL) {
  auto sp = std::make_shared<H240bSplit>(); sp->items = std::move(*items); sp->ntot = sp->items.size(); sp->f = 1.0;
  auto bufs = std::make_shared<std::array<P14Buf, 2>>(); (*bufs)[0] = SBf[0]; (*bufs)[1] = SBf[1];
  SL.h299_n = sp->items.size(); SL.h299_ms = 0;
  if (sp->items.empty()) return;
  double* msp = &SL.h299_ms;
  SL.h299j = h299_submit([sp, bufs, tok, psp, tvec, cb, msp]() {
    [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H406 grid/band command buffer failed\n"); gpu_fail(); }
    double ms = 0; h240b_cpu(*sp, bufs->data(), tok, psp, tvec, ms); *msp = ms;
  });
}


struct H670Job {
  H240bSplit sp; P14Buf SBf[2]; id<MTLBuffer> tok, psp, tv; id<MTLCommandBuffer> cbG, cbB, cbP;
  std::thread th; double cpu_ms = 0; bool pending = false;
  std::vector<int> K; std::vector<float> r8;
  std::vector<uint32_t> cs; bool cs_ok = false; int step = 0;
};
static H670Job g_h670j;
static id<MTLCommandQueue> g_h670q = nil;
static void h670_join() {
  H670Job& J = g_h670j;
  if (!J.pending) return;
  const auto t0 = clk::now();
  if (J.th.joinable()) J.th.join();
  for (id<MTLCommandBuffer> cb : {J.cbG, J.cbB, J.cbP})
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H670 command buffer failed\n"); h182_dump(cb, "H670"); gpu_fail(); }
  const double gpu_ms = (J.cbP.GPUEndTime - J.cbP.GPUStartTime) * 1e3;
  const size_t ngpu = J.sp.ntot - J.sp.items.size();
  h240b_observe(J.sp.ntot, ngpu, gpu_ms, J.sp.items.size(), J.cpu_ms);

  const float* rs0d = (const float*)J.SBf[0].rs.contents; const float* rs1d = (const float*)J.SBf[1].rs.contents;
  g_h670_res.clear();
  const size_t nK = J.K.size();
  for (size_t k = 0; k < nK; k++) {
    const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
    const double r8x = J.r8[2 * k], r8y = J.r8[2 * k + 1];
    if (!(std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8x + r8y)) continue;
    const float rf0 = rs0d[2 * k + 1], rf1 = rs1d[2 * k + 1];
    if (!std::isfinite(rf0) || !std::isfinite(rf1)) continue;
    g_h670_res.push_back({J.K[k], (float)std::min((double)rf0 + (double)rf1, (double)DEAD)});
  }
  if (!J.cs_ok) { fprintf(stderr, "H670: polish joined before its step's select\n"); exit(1); }
  g_h670_cs.swap(J.cs); J.cs_ok = false; g_h670_jst = J.step;
  g_h670_have = true; J.pending = false;
  J.cbG = nil; J.cbB = nil; J.cbP = nil;
  const double w = ms_since(t0); g_h670_wait += w; g_h670_n++;
  printf("H670 join: wait %.2f ms | cpu %.1f gpu %.1f ms (trees cpu %zu gpu %zu) | accepted %zu of %zu\n", w, J.cpu_ms, gpu_ms, J.sp.items.size(), ngpu, g_h670_res.size(), nK);
}
static const int g_h731 = getenv("H731") ? atoi(getenv("H731")) : 0;
static std::atomic<double> g_h731_gpu[2];
static void h731_wait(H53Scr& HS, double* gpu) {
  for (id<MTLCommandBuffer> cb : HS.h731_cbs) { [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H731 command buffer failed\n"); h182_dump(cb, "H731"); gpu_fail(); }
    if (gpu) gpu[1] += (cb.GPUEndTime - cb.GPUStartTime) * 1e3; }
  HS.h731_cbs.clear();
}

static id<MTLBuffer> h731_gout(id<MTLDevice> dev, H53Scr& HS, size_t bytes) {
  bytes = std::max<size_t>(bytes, 256);
  if (HS.gout.b && HS.gout.cap >= bytes) return HS.gout.b;
  h731_wait(HS, nullptr);
  const size_t nc = mem_grow_cap(bytes, HS.gout.cap);
  id<MTLBuffer> nb = h900_nbl(dev, __LINE__, nc, MSH);
  if (HS.gout.b) std::memcpy(nb.contents, HS.gout.b.contents, HS.gout.cap);
  HS.gout.b = nb; HS.gout.cap = nc; return nb;
}
static void h731_enc_grid(Engine& E, H53Scr& HS, id<MTLCommandBuffer> cb, size_t gb0, id<MTLBuffer> bGout, id<MTLComputePipelineState> gps) {
  if (HS.gbins.size() <= gb0) return;
  id<MTLDevice> dev = E.dev;
  id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
  size_t gi = gb0;
  for (size_t b = gb0; b < HS.gbins.size(); b++) { const H117G& g = HS.gbins[b];
    id<MTLBuffer> tl = (g_h131 & 1) ? h131_put(dev, g_h131p.tl.size() > gi ? g_h131p.tl[gi] : (g_h131p.tl.resize(gi + 1), g_h131p.tl[gi]), g.tl.data(), g.tl.size() * 4)
                                    : E.B(g.tl.data(), g.tl.size() * 4);
    gi++;
    [en setComputePipelineState:gps];
    id<MTLBuffer> bs[10] = {HS.tok.b, HS.bHdr, HS.sp.b, HS.bSg, HS.bTh[g.side], HS.bAa[g.side], HS.bQ, tl, HS.bOid[g.side], bGout};
    for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
    [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
    [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
  }
  [en endEncoding];
}
struct H731BandBufs { id<MTLBuffer> gout, rows, curve, str, tb, td, tv_tab, ot, os, oo, tv; };
static void h731_enc_band(id<MTLCommandBuffer> cb, const H731BandBufs& B, size_t k0, size_t k1,
                          id<MTLComputePipelineState> pg, id<MTLComputePipelineState> pb, id<MTLComputePipelineState> pr) {
  if (k1 <= k0) return;
  const uint32_t Pk = (uint32_t)(k1 - k0); const size_t nTile = ((size_t)Pk + 31) / 32;
  id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
  [en setComputePipelineState:pg];
  [en setBuffer:B.gout offset:0 atIndex:0]; [en setBuffer:B.rows offset:8 * k0 atIndex:1]; [en setBuffer:B.curve offset:(k0 / 32) * 16384 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
  [en dispatchThreadgroups:MTLSizeMake(((size_t)Pk * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
  [en endEncoding];
  en = tl_enc(cb, MTLDispatchTypeSerial);
  [en setComputePipelineState:pb];
  [en setBuffer:B.curve offset:(k0 / 32) * 16384 atIndex:0]; [en setBuffer:B.str offset:0 atIndex:1]; [en setBuffer:B.tb offset:0 atIndex:2]; [en setBuffer:B.td offset:0 atIndex:3];
  [en setBuffer:B.tv_tab offset:0 atIndex:4]; [en setBuffer:B.ot offset:192 * k0 atIndex:5]; [en setBuffer:B.os offset:48 * k0 atIndex:6]; [en setBuffer:B.oo offset:4 * k0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
  [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
  [en endEncoding];
  en = tl_enc(cb, MTLDispatchTypeConcurrent);
  [en setComputePipelineState:pr];
  [en setBuffer:B.ot offset:192 * k0 atIndex:0]; [en setBuffer:B.os offset:48 * k0 atIndex:1]; [en setBuffer:B.tv offset:192 * k0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
  [en dispatchThreadgroups:MTLSizeMake(((size_t)Pk + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
  [en endEncoding];
}

static H731BandBufs h731_bufs(Engine& E, H53Scr& HS, size_t n, size_t NL) {
  id<MTLDevice> dev = E.dev;
  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nT = (n + 31) / 32;
  H731BandBufs B;
  B.gout = h731_gout(dev, HS, NL * 128 * 4); B.rows = HS.rows.b;
  B.curve = h131_get(dev, g_h131p.curve, nT * 4096 * 4); B.ot = h131_get(dev, g_h131p.ot, n * 48 * 4); B.os = h131_get(dev, g_h131p.os, n * 48);
  B.oo = h131_get(dev, g_h131p.oo, n * 4); B.tv = h131_get(dev, g_h131p.tv, n * 48 * 4 + 64);
  B.str = h131_put(dev, g_h131p.str, STR, sizeof STR); B.tb = h131_put(dev, g_h131p.tb, BT.base, sizeof BT.base);
  B.td = h131_put(dev, g_h131p.td, BT.del, sizeof BT.del); B.tv_tab = h131_put(dev, g_h131p.tvl, BT.valid, sizeof BT.valid);
  return B;
}

static void h731_chunk(Engine& E, PipeStep& PS, H53Scr& HS, Times& T) {
  if (!HS.open || HS.K.empty()) return;
  const auto t0 = clk::now();
  const size_t kb1 = HS.k_done / 32 * 32;
  if (HS.gbins.size() <= HS.h731_gb && kb1 <= HS.h731_kb) return;
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
  if (!gps || !pg || !pb || !pr) { fprintf(stderr, "H731: grid/band pipeline missing\n"); exit(1); }
  const H731BandBufs B = h731_bufs(E, HS, (size_t)PS.n, HS.ulist.size());
  id<MTLCommandBuffer> cb = h182_cb(E.cq);
  h731_enc_grid(E, HS, cb, HS.h731_gb, B.gout, gps);
  h731_enc_band(cb, B, HS.h731_kb, kb1, pg, pb, pr);
  [cb commit];
  HS.h731_cbs.push_back(cb);
  HS.h731_gb = HS.gbins.size(); HS.h731_kb = std::max(HS.h731_kb, kb1);
  T.h131_enc += ms_since(t0);
}

static id<MTLBuffer> h900_gout(id<MTLDevice> dev, H53Scr& HS, size_t bytes) {
  H53Buf& g = HS.gout; bytes = std::max<size_t>(bytes, 256);
  if (g.b && g.cap >= bytes && !(g_h900 & 8)) return g.b;
  if ((g_h900 & 2) && h900_inplace(dev, g.b, g.cap, g.ip, g.ipres, g.ipst, bytes, "gout")) return g.b;
  return g.get(dev, bytes);
}


static std::vector<std::vector<P14Bin>> h117b_split(const std::vector<P14Bin>& src, size_t gr, size_t nblk) {
  std::vector<std::vector<P14Bin>> out(nblk);
  for (const P14Bin& g : src) {
    std::vector<std::vector<uint32_t>> ks(nblk);
    for (size_t i = 1; i < g.tl.size(); i++) if (g.tl[i] != 0xFFFFFFFFu) ks[g.tl[i] / gr].push_back(g.tl[i]);
    for (size_t b = 0; b < nblk; b++) {
      if (ks[b].empty()) continue;
      const size_t tp = (size_t)g.tpg;
      P14Bin s = g;
      s.nt = (ks[b].size() + tp - 1) / tp * tp;
      const bool xf = g.tpg > 1 && g.ntg == g.nt;
      s.ntg = xf ? s.nt : s.nt / tp;
      s.eb_bytes = s.ntg * tp * 2 * (size_t)(g.Q + 1) * 64 * 4;
      s.tl.assign(s.nt + 1, 0xFFFFFFFFu); s.tl[0] = g.tl[0];
      for (size_t i = 0; i < ks[b].size(); i++) s.tl[1 + i] = ks[b][i];
      out[b].push_back(s);
    }
  }
  return out;
}


static const int g_post_sicis_gpu = getenv("POST_SICIS_GPU") ? atoi(getenv("POST_SICIS_GPU")) : 0;


static const int g_post_sicis_yfx = getenv("POST_SICIS_YFX") ? atoi(getenv("POST_SICIS_YFX")) : 0;
static const int g_post_sicis_check = getenv("POST_SICIS_CHECK") ? atoi(getenv("POST_SICIS_CHECK")) : 256;
static id<MTLComputePipelineState> g_psRes = nil, g_psScore = nil;
static std::string g_psKdir;

static bool post_sicis_cfc(const Sketch::CisPerCellEnergyInput& se, std::vector<float>& cf) {
  const auto& me = se.memberEnergy; if (me.empty() || me.size() > 8 || Sketch::cisEnsembleModeFromEnv() == 0) return false;
  cf.assign(400, 0.0f); double K[8] = {0}, Kmin = 1e300;
  for (size_t m = 0; m < me.size(); m++) {
    if (me[m].size() != 48) return false;
    for (size_t c = 0; c < 48; c++) { const double p = std::max(me[m][c], 0.0); if (!(p > 0.0) || !std::isfinite(p)) return false;
      const double h = 0.5 / p; if (!std::isfinite((float)h)) return false; cf[8 + m * 48 + c] = (float)h; K[m] += 0.5 * (std::log(2.0 * 3.14159265358979323846) + std::log(p)); }
    Kmin = std::min(Kmin, K[m]);
  }
  const bool joint = Sketch::cisJointCoarseFromEnv() != 0;
  if (joint) {
    const int blkN[7] = {24, 12, 6, 3, 1, 1, 1};
    Kmin = 1e300;
    for (size_t m = 0; m < me.size(); m++) {
      double Km = 0.0;
      for (int b = 0; b < 7; b++) { const double D = blkN[b]; Km += 0.5 * D * std::log(2.0 * 3.14159265358979323846) - std::lgamma(0.5 * D); }
      for (size_t c = 0; c < 48; c++) Km += 0.5 * std::log(std::max(me[m][c], 0.0));
      K[m] = Km; Kmin = std::min(Kmin, Km);
    }
  }
  const float Khi = (float)Kmin; cf[0] = (float)me.size(); cf[1] = 0.0f; cf[2] = Khi; cf[3] = (float)(Kmin - (double)Khi);
  cf[4] = joint ? 1.0f : 0.0f;
  for (size_t m = 0; m < me.size(); m++) cf[392 + m] = (float)(K[m] - Kmin);
  return true;
}
static bool post_sicis_gpu_init(id<MTLDevice> dev) {
  if (g_psRes && g_psScore) return true;
  NSError* e = nil;
  id<MTLLibrary> L = [dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(g_psKdir + "/post_sicis_gpu.metallib").c_str()]] error:&e];
  if (!L) { fprintf(stderr, "post_sicis_gpu.metallib: %s\n", e.description.UTF8String); return false; }
  g_psRes = [dev newComputePipelineStateWithFunction:[L newFunctionWithName:((g_nres || g_yfx) ? @"post_sicis_residual_n" : @"post_sicis_residual")] error:&e];
  if (!g_psRes) { fprintf(stderr, "pso post_sicis_residual: %s\n", e.description.UTF8String); return false; }
  g_psScore = [dev newComputePipelineStateWithFunction:[L newFunctionWithName:@"post_sicis_score"] error:&e];
  if (!g_psScore) { fprintf(stderr, "pso post_sicis_score: %s\n", e.description.UTF8String); return false; }
  printf("POST_SICIS_GPU=1: residual+score kernels from post_sicis_gpu.metallib (check every %d)\n", g_post_sicis_check);
  return true;
}

static id<MTLBuffer> post_sicis_new(id<MTLDevice> dev, size_t bytes) { id<MTLBuffer> b = [dev newBufferWithLength:std::max<size_t>(bytes, 256) options:MTLResourceStorageModeShared]; memdiag::metal_note(b, bytes, false); return b; }

struct PostSicisOut { id<MTLBuffer> err = nil, bad = nil; size_t lo = 0, n = 0; };
static bool post_sicis_enc(Engine& E, id<MTLCommandBuffer> cb, const std::vector<float>& cf, const StepBatch& SB, const std::vector<int>& ct,
                           const std::vector<int>& K, size_t lo, size_t hi, const P14Buf* SBf, id<MTLBuffer> bTV, PostSicisOut& o) {
  id<MTLDevice> dev = E.dev;
  const size_t n = hi - lo;
  std::vector<uint32_t> recpool, offv(2 * n);
  for (size_t q = 0; q < n; q++) for (int s = 0; s < 2; s++) {
    const int j = ct[(size_t)K[lo + q] * 2 + s];
    offv[2 * q + s] = (uint32_t)recpool.size();
    const Tree& tr = SB.UQ.U[j];
    recpool.insert(recpool.end(), tr.rec.begin(), tr.rec.end());
  }
  recpool.push_back(0u);
  id<MTLBuffer> bRec = post_sicis_new(dev, recpool.size() * 4);
  memcpy(bRec.contents, recpool.data(), recpool.size() * 4);
  id<MTLBuffer> bOff = post_sicis_new(dev, offv.size() * 4);
  memcpy(bOff.contents, offv.data(), offv.size() * 4);
  id<MTLBuffer> bY = post_sicis_new(dev, 2 * 48 * 4), bC = post_sicis_new(dev, cf.size() * 4);
  float* yp = (float*)bY.contents;
  for (int s = 0; s < 2; s++) for (int p = 0; p < 48; p++) yp[s * 48 + p] = (float)(g_yfx ? g_yfx_y[s][p] : E.tg[s].y[p]);
  memcpy(bC.contents, cf.data(), cf.size() * 4);
  id<MTLBuffer> bR = post_sicis_new(dev, n * 48 * 4), bB = post_sicis_new(dev, n * 4), bF = post_sicis_new(dev, n * 4);
  o.err = post_sicis_new(dev, n * 4); o.bad = bB; o.lo = lo; o.n = n;
  uint32_t nn = (uint32_t)n;

  id<MTLComputeCommandEncoder> en = [cb computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
  [en setComputePipelineState:g_psRes];
  [en setBuffer:SBf[0].th offset:lo * 32 * 4 atIndex:0]; [en setBuffer:SBf[1].th offset:lo * 32 * 4 atIndex:1]; if (g_yfx) { static id<MTLBuffer> bTfx = nil; if (!bTfx) { bTfx = post_sicis_new(dev, 4096 * 48 * 4); float* tp = (float*)bTfx.contents; for (int q = 0; q < 4096; q++) for (int p = 0; p < 48; p++) tp[q * 48 + p] = (float)g_yfx_t[p]; }
    [en setBuffer:bTfx offset:0 atIndex:2]; } else [en setBuffer:bTV offset:lo * 48 * 4 atIndex:2];
  [en setBuffer:bRec offset:0 atIndex:3]; [en setBuffer:bOff offset:0 atIndex:4]; [en setBuffer:bY offset:0 atIndex:5];
  [en setBuffer:bR offset:0 atIndex:6]; [en setBuffer:bB offset:0 atIndex:7]; [en setBytes:&nn length:4 atIndex:8];
  if (g_yfx) { float nf[96]; for (int p = 0; p < 96; p++) nf[p] = (float)g_yfx_n[p]; [en setBytes:nf length:sizeof nf atIndex:9]; }
  else if (g_nres) [en setBytes:g_nres_n length:sizeof g_nres_n atIndex:9];
  [en dispatchThreadgroups:MTLSizeMake((n + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
  [en endEncoding];
  en = [cb computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
  [en setComputePipelineState:g_psScore];
  [en setBuffer:bR offset:0 atIndex:0]; [en setBuffer:bC offset:0 atIndex:1]; [en setBuffer:o.err offset:0 atIndex:2]; [en setBuffer:bF offset:0 atIndex:3]; [en setBytes:&nn length:4 atIndex:4];
  [en dispatchThreadgroups:MTLSizeMake((n + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
  [en endEncoding];
  return true;
}


static id<MTLComputePipelineState> g_pgGeo = nil, g_pgScore = nil;
static bool pick_gpu_init(id<MTLDevice> dev, double& ms) {
  if (g_pgGeo && g_pgScore) return true;
  const auto t0 = std::chrono::steady_clock::now();
  NSError* e = nil;
  id<MTLLibrary> L1 = [dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(g_psKdir + "/geo_res48.metallib").c_str()]] error:&e];
  if (!L1) { fprintf(stderr, "geo_res48.metallib: %s\n", e.description.UTF8String); return false; }
  id<MTLLibrary> L2 = [dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(g_psKdir + "/post_sicis_gpu.metallib").c_str()]] error:&e];
  if (!L2) { fprintf(stderr, "post_sicis_gpu.metallib: %s\n", e.description.UTF8String); return false; }
  g_pgGeo = [dev newComputePipelineStateWithFunction:[L1 newFunctionWithName:@"geo_res48"] error:&e];
  if (!g_pgGeo) { fprintf(stderr, "pso geo_res48: %s\n", e.description.UTF8String); return false; }
  g_pgScore = [dev newComputePipelineStateWithFunction:[L2 newFunctionWithName:@"post_sicis_score"] error:&e];
  if (!g_pgScore) { fprintf(stderr, "pso post_sicis_score: %s\n", e.description.UTF8String); g_pgGeo = nil; return false; }
  ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return true;
}

static bool pick_gpu_run(Engine& E, const fitx::GeoPack& GP, const float* Yv, const std::vector<float>& cf, int chunk, float gsx, float gsy, int sgn_tan,
                         std::vector<float>& nll, std::vector<uint8_t>& bad, std::vector<float>& cond, double& ms_init, double& ms_gpu) {
  if (!pick_gpu_init(E.dev, ms_init)) return false;
  if (g_pgGeo.maxTotalThreadsPerThreadgroup < 256) { fprintf(stderr, "geo_res48: maxTotalThreadsPerThreadgroup %lu < 256\n", (unsigned long)g_pgGeo.maxTotalThreadsPerThreadgroup); return false; }
  id<MTLDevice> dev = E.dev;
  std::vector<uint32_t> cm; cm.reserve((size_t)GP.npacked);
  for (size_t r = 0; r < GP.st.size(); r++) if (GP.st[r] == 1) cm.push_back(GP.moff[r]);
  const size_t nc = cm.size(); nll.assign(nc, NAN); bad.assign(nc, 1); cond.assign(nc, NAN);
  if (!nc) return true;
  const auto t0 = std::chrono::steady_clock::now();
  chunk = std::max(1, std::min<int>(chunk, (int)nc));
  id<MTLBuffer> bP = [dev newBufferWithBytes:GP.proto.data() length:GP.proto.size() * 4 options:MTLResourceStorageModeShared];
  id<MTLBuffer> bM = [dev newBufferWithBytes:GP.mrec.data() length:GP.mrec.size() * 4 options:MTLResourceStorageModeShared];
  id<MTLBuffer> bO = [dev newBufferWithBytes:cm.data() length:nc * 4 options:MTLResourceStorageModeShared];
  id<MTLBuffer> bY = [dev newBufferWithBytes:Yv length:96 * 4 options:MTLResourceStorageModeShared];
  id<MTLBuffer> bR = post_sicis_new(dev, nc * 48 * 4), bB = post_sicis_new(dev, nc * 4), bS = post_sicis_new(dev, (size_t)chunk * 48 * 1024 * 2), bK = post_sicis_new(dev, nc * 4);
  id<MTLBuffer> bC = [dev newBufferWithBytes:cf.data() length:cf.size() * 4 options:MTLResourceStorageModeShared];
  id<MTLBuffer> bN = post_sicis_new(dev, nc * 4), bF = post_sicis_new(dev, nc * 4);
  id<MTLCommandBuffer> cb = [E.cq commandBuffer];
  id<MTLComputeCommandEncoder> en = [cb computeCommandEncoder];
  [en setComputePipelineState:g_pgGeo];
  [en setBuffer:bP offset:0 atIndex:0]; [en setBuffer:bM offset:0 atIndex:1]; [en setBuffer:bO offset:0 atIndex:2]; [en setBuffer:bY offset:0 atIndex:3];
  [en setBuffer:bR offset:0 atIndex:4]; [en setBuffer:bB offset:0 atIndex:5]; [en setBuffer:bS offset:0 atIndex:6]; [en setBuffer:bK offset:0 atIndex:8];
  for (size_t base = 0; base < nc; base += (size_t)chunk) {
    struct { uint32_t base, nchunk, sgn, pad; float gsx, gsy, p2, p3; } ga = {(uint32_t)base, (uint32_t)std::min<size_t>((size_t)chunk, nc - base), (uint32_t)sgn_tan, 0u, gsx, gsy, 0.f, 0.f};
    [en setBytes:&ga length:sizeof ga atIndex:7];
    [en dispatchThreadgroups:MTLSizeMake(ga.nchunk, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
  }
  [en endEncoding];
  en = [cb computeCommandEncoder];
  [en setComputePipelineState:g_pgScore];
  uint32_t nn = (uint32_t)nc;
  [en setBuffer:bR offset:0 atIndex:0]; [en setBuffer:bC offset:0 atIndex:1]; [en setBuffer:bN offset:0 atIndex:2]; [en setBuffer:bF offset:0 atIndex:3]; [en setBytes:&nn length:4 atIndex:4];
  [en dispatchThreadgroups:MTLSizeMake((nc + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
  [en endEncoding];
  [cb commit]; [cb waitUntilCompleted];
  if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "geo_res48 command buffer: %s\n", cb.error.description.UTF8String); return false; }
  const float* pn = (const float*)bN.contents; const uint32_t* pb = (const uint32_t*)bB.contents; const uint32_t* pf = (const uint32_t*)bF.contents;
  const float* pk_ = (const float*)bK.contents;
  for (size_t c = 0; c < nc; c++) { cond[c] = pk_[c]; nll[c] = pn[c]; bad[c] = (pb[c] != 0u || pf[c] != 0u || !std::isfinite(pn[c]) || pn[c] >= 1.0e29f) ? 1 : 0; }
  ms_gpu += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return true;
}

static void post16_h117_fin(Engine& E, PipeStep& PS, StepBatch& SB, float* err, float* rx, float* ry, Times& T, bool cpu_here) {
  h301_mark("post fin: alloc/bind/grid enc");
  if (g_h675) h670_join();
  h301_mark("post fin: m1 h670 join done");
  auto tw = clk::now();
  PostWall14 pw{T, tw};
  H53Scr& HS = g_h53s;
  id<MTLDevice> dev = E.dev;
  const std::vector<int>& ct = SB.ct;
  struct OpenGuard { bool& o; ~OpenGuard() { o = false; } } open_guard{HS.open};
  Loc16Reset loc_reset{g_loc16, HS.ulist};
  if (!HS.open || HS.K0.empty()) { h731_wait(HS, nullptr); return; }
  if (cpu_here) h117_cpu(E, PS, SB, HS, T);
  h301_mark("post fin: m2 h117_cpu done");
  const std::vector<int>& K = HS.K; const int nK = (int)K.size();
  T.p_cand += nK;
  if (nK == 0) { h731_wait(HS, nullptr); return; }
  const size_t h731_g0 = g_h731 ? HS.h731_gb : 0, h731_k0 = g_h731 ? HS.h731_kb : 0;
  const size_t NL = HS.ulist.size();
  id<MTLComputePipelineState> gps = E.pipe(E.Lgr, "grid_h26");
  id<MTLComputePipelineState> pg = E.pipe(E.Lb, "post_gather_i14"), pb = E.pipe(E.Lb, "band_lane"), pr = E.pipe(E.Lb, "post_repair_i14");
  if (!gps || !pg || !pb || !pr) { fprintf(stderr, "H117: grid/band pipeline missing\n"); exit(1); }
  static band::Tables BT = band::make_tables();
  static float STR[2 * band::NROW];
  for (int i = 0; i < band::NROW; i++) { STR[2 * i] = E.tg[0].yf[i]; STR[2 * i + 1] = E.tg[1].yf[i]; }
  const size_t nTile = ((size_t)nK + 31) / 32;
  const uint32_t Pk = (uint32_t)nK;
  id<MTLBuffer> bGout = g_h731 ? h731_gout(dev, HS, NL * 128 * 4) : h900_gout(dev, HS, std::max<size_t>(NL * 128 * 4, 256));
  id<MTLBuffer> bRows = HS.rows.b;
  const auto ta0 = clk::now();
  id<MTLBuffer> bCurve, bOT, bOS, bOO, bTV;
  if (g_h731) {
    const H731BandBufs B7 = h731_bufs(E, HS, (size_t)PS.n, NL);
    bCurve = B7.curve; bOT = B7.ot; bOS = B7.os; bOO = B7.oo; bTV = B7.tv;
  } else if (g_h131 & 1) {
    bCurve = h131_get(dev, g_h131p.curve, nTile * 4096 * 4);
    bOT = h131_get(dev, g_h131p.ot, (size_t)nK * 48 * 4);
    bOS = h131_get(dev, g_h131p.os, (size_t)nK * 48);
    bOO = h131_get(dev, g_h131p.oo, (size_t)nK * 4);
    bTV = h131_get(dev, g_h131p.tv, (size_t)nK * 48 * 4 + 64);
  } else {
    bCurve = h900_nbl(dev, __LINE__, nTile * 4096 * 4, MSH);
    std::memset(bCurve.contents, 0, nTile * 4096 * 4);
    bOT = h900_nbl(dev, __LINE__, (size_t)nK * 48 * 4, MSH);
    bOS = h900_nbl(dev, __LINE__, (size_t)nK * 48, MSH);
    bOO = h900_nbl(dev, __LINE__, (size_t)nK * 4, MSH);
    bTV = h900_nbl(dev, __LINE__, (size_t)nK * 48 * 4 + 64, MSH);
  }
  h301_mark("post fin: m3 buffers");
  const double h131_alloc_ms = ms_since(ta0);
  const auto tb0 = clk::now();
  id<MTLBuffer> bStr, bTb, bTd, bTv;
  if (g_h131 & 1) {
    bStr = h131_put(dev, g_h131p.str, STR, sizeof STR); bTb = h131_put(dev, g_h131p.tb, BT.base, sizeof BT.base);
    bTd = h131_put(dev, g_h131p.td, BT.del, sizeof BT.del); bTv = h131_put(dev, g_h131p.tvl, BT.valid, sizeof BT.valid);
  } else {
    bStr = E.B(STR, sizeof STR); bTb = E.B(BT.base, sizeof BT.base); bTd = E.B(BT.del, sizeof BT.del); bTv = E.B(BT.valid, sizeof BT.valid);
  }
  P14Buf SBf[2];
  for (int s = 0; s < 2; s++) {
    SBf[s].hdr = HS.phb[s].b; SBf[s].segs = HS.psgb[s].b; SBf[s].th = HS.pthb[s].b; SBf[s].rs = HS.prsb[s].b;
    SBf[s].a = HS.pab[s].b; SBf[s].stio = HS.stiob[s].b;
    SBf[s].yv = (g_h131 & 1) ? h131_put(dev, g_h131p.yv[s], E.tg[s].yf, sizeof E.tg[s].yf) : E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  }
  const double h131_bind_ms = ms_since(tb0);
  T.h131_alloc += h131_alloc_ms; T.h131_bind += h131_bind_ms;
  auto enc_grid = [&](id<MTLCommandBuffer> cb) {
    if (g_h731) { h731_enc_grid(E, HS, cb, h731_g0, bGout, gps); return; }
    if (HS.gbins.empty()) return;
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    size_t gi = 0;
    for (const H117G& g : HS.gbins) {
      id<MTLBuffer> tl = (g_h131 & 1) ? h131_put(dev, g_h131p.tl.size() > gi ? g_h131p.tl[gi] : (g_h131p.tl.resize(gi + 1), g_h131p.tl[gi]), g.tl.data(), g.tl.size() * 4)
                                      : E.B(g.tl.data(), g.tl.size() * 4);
      gi++;
      [en setComputePipelineState:gps];
      id<MTLBuffer> bs[10] = {HS.tok.b, HS.bHdr, HS.sp.b, HS.bSg, HS.bTh[g.side], HS.bAa[g.side], HS.bQ, tl, HS.bOid[g.side], bGout};
      for (int k = 0; k < 10; k++) [en setBuffer:bs[k] offset:0 atIndex:k];
      [en setBytes:g.pd length:8 atIndex:10]; [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
    }
    [en endEncoding];
  };
  auto enc_band = [&](id<MTLCommandBuffer> cb) {
    if (g_h731) { const H731BandBufs B7{bGout, bRows, bCurve, bStr, bTb, bTd, bTv, bOT, bOS, bOO, bTV}; h731_enc_band(cb, B7, h731_k0, (size_t)nK, pg, pb, pr); return; }
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pg];
    [en setBuffer:bGout offset:0 atIndex:0]; [en setBuffer:bRows offset:0 atIndex:1]; [en setBuffer:bCurve offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK * 128 + 127) / 128, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeSerial);
    [en setComputePipelineState:pb];
    [en setBuffer:bCurve offset:0 atIndex:0]; [en setBuffer:bStr offset:0 atIndex:1]; [en setBuffer:bTb offset:0 atIndex:2]; [en setBuffer:bTd offset:0 atIndex:3];
    [en setBuffer:bTv offset:0 atIndex:4]; [en setBuffer:bOT offset:0 atIndex:5]; [en setBuffer:bOS offset:0 atIndex:6]; [en setBuffer:bOO offset:0 atIndex:7]; [en setBytes:&Pk length:4 atIndex:8];
    [en dispatchThreadgroups:MTLSizeMake(nTile, 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [en endEncoding];
    en = tl_enc(cb, MTLDispatchTypeConcurrent);
    [en setComputePipelineState:pr];
    [en setBuffer:bOT offset:0 atIndex:0]; [en setBuffer:bOS offset:0 atIndex:1]; [en setBuffer:bTV offset:0 atIndex:2]; [en setBytes:&Pk length:4 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake(((size_t)nK + 63) / 64, 1, 1) threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [en endEncoding];
  };
  auto enc_pol = [&](id<MTLCommandBuffer> cb) {
    g_h899_fin = &g_h117_sl;
    p14_encode_h91(E, cb, HS.pbins, SBf, HS.ptok.b, HS.psp.b, bTV, 3.0f, 0.0f, g_h117_sl);
    if (g_h894 & 1) g_h894_need = g_h894_last_eo;
    if (g_h891 && g_h891_gpu_us > 0) {
      h680_init(E.dev);
      const uint32_t nit = (uint32_t)std::max(1.0, g_h891_gpu_us * 1e-3 / g_h680_ms_per_iter);
      id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeSerial);
      [en setComputePipelineState:g_h680ps]; [en setBuffer:g_h680o offset:0 atIndex:0]; [en setBytes:&nit length:4 atIndex:1];
      [en setBuffer:bTV offset:0 atIndex:2];
      [en dispatchThreadgroups:MTLSizeMake(256, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)]; [en endEncoding];
      g_h891_gpu_n += 1; g_h891_gpu_req_ms += g_h891_gpu_us * 1e-3; g_h891_gpu_cal_ms += nit * g_h680_ms_per_iter;
    }
  };
  double gpu_t[3] = {0, 0, 0}, wall_t[3] = {0, 0, 0};

  const bool h898_def = (g_h898 & 5) && !HS.h898_rng.empty();
  const bool h898_ovl = h898_def && !((g_h117 & 2) && !g_h240b) && !(g_h670 && g_h670_st > 0 && g_h240b) && (g_h257 & 1) && g_h240b;
  if (h898_def && !h898_ovl) h898_flush(E, SB, HS, T);
  auto run_cb = [&](id<MTLCommandBuffer> cb, int which) {
    auto t0 = clk::now();
    [cb commit]; [cb waitUntilCompleted];
    if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H117 command buffer failed\n"); h182_dump(cb, "H117"); gpu_fail(); }
    gpu_t[which] += (cb.GPUEndTime - cb.GPUStartTime) * 1e3; wall_t[which] += ms_since(t0);
  };

  T.h131_pre += h131_alloc_ms + h131_bind_ms;
  auto h131_enc_t = [&](const std::function<void()>& f) { const auto q = clk::now(); f(); const double ms = ms_since(q); T.h131_enc += ms; T.h131_pre += ms; };
  bool h670_on = false;

  const float* rs0d = (const float*)SBf[0].rs.contents; const float* rs1d = (const float*)SBf[1].rs.contents;
  if (g_dump14.size() && !g_dumpf) g_dumpf = fopen(g_dump14.c_str(), "w");
  if (g_dumpf && HS.call == 1) fprintf(g_dumpf, "call,child,rx8,ry8,rxf,ryf,acc0,acc1\n");
  struct OutCnt { double viol = 0, up = 0, a0 = 0, a1 = 0, fb = 0, pchk = 0, pchkn = 0, pfb = 0; };


  static const int g_pfr = getenv("POSTFIN_REUSE") ? atoi(getenv("POSTFIN_REUSE")) : 1;
  static const size_t PFG = getenv("K_PFG") ? (size_t)std::max(16, atoi(getenv("K_PFG"))) : 512;

  static const size_t jblk = getenv("J_BLK") ? (size_t)std::max(256, atoi(getenv("J_BLK"))) : 4096;
  const size_t HB = g_post_sicis_gpu ? 4096 : jblk;
  std::vector<OutCnt> oc(std::max({(size_t)h84_nb((size_t)nK, 4096), ((size_t)nK + HB - 1) / HB + 1, ((size_t)nK + PFG - 1) / PFG}));
  const bool want_dump = g_dumpf != nullptr;
  std::vector<float> ORW(want_dump ? (size_t)nK * 4 : 0);
  std::vector<uint8_t> OAC(want_dump ? (size_t)nK * 2 : 0);

  std::vector<PostSicisOut> psout(std::max<size_t>(1, ((size_t)nK + 4095) / 4096));
  std::vector<float> psCf;
  const bool psOn = g_post_sicis_gpu && (!g_yfx || g_post_sicis_yfx) && g_sicis_en && g_sicis_en->present && !HS.pbins.empty() && post_sicis_cfc(*g_sicis_en, psCf) && post_sicis_gpu_init(dev);
  static const int jptr = getenv("J_PTR") ? atoi(getenv("J_PTR")) : 1;
  const float* jthp[2] = {(const float*)SBf[0].th.contents, (const float*)SBf[1].th.contents}; const float* jtvp = (const float*)bTV.contents;
  auto outputs = [&](size_t lo, size_t hi, int bi) {
    const bool acc_off = g_h660 || h670_on;
    OutCnt& cn = oc[(size_t)bi];
#if defined(__aarch64__)

    const bool yfxb = g_yfxfast && g_yfx && jptr && sicis::fastEnabled() && NPT == 48 && !getenv("FJN_S1") && g_sicis_en && g_sicis_en->present;
    thread_local std::vector<float> ypre_res, ypre_ce; thread_local std::vector<uint8_t> ypre_have, ypre_need;
    if (yfxb) { const size_t nn_ = hi - lo; ypre_res.resize(nn_ * NPT); ypre_ce.resize(nn_); ypre_have.assign(nn_, 0); ypre_need.assign(nn_, 0);
      for (size_t k = lo; k < hi; k++) { const PostSicisOut& po = psout[k / 4096]; const size_t q = k - po.lo; const bool gok = po.err && ((const uint32_t*)po.bad.contents)[q] == 0;
        ypre_need[k - lo] = (!gok || (g_post_sicis_check > 0 && k % (size_t)g_post_sicis_check == 0)) ? 1 : 0; }
      yfx_chunk(E, SB.UQ, ct, K, lo, hi, ypre_need.data(), SBf, bTV, jptr ? jthp : nullptr, jptr ? jtvp : nullptr, ypre_res.data(), ypre_ce.data(), ypre_have.data()); }
#else
    const bool yfxb = false; std::vector<float> ypre_res, ypre_ce; std::vector<uint8_t> ypre_have;
#endif
    for (size_t k = lo; k < hi; k++) {
      const int i = K[k];
      float fin[2], r8[2]; int ac[2];
      const float rf0 = rs0d[2 * k + 1], rf1 = rs1d[2 * k + 1];
      for (int s = 0; s < 2; s++) { const int j = ct[(size_t)i * 2 + s]; r8[s] = SB.F[s].rfin[j]; }
      const double rbx = rs0d[2 * k], rby = rs1d[2 * k];
      const double r8sum = (double)SB.F[0].rfin[ct[(size_t)i * 2]] + (double)SB.F[1].rfin[ct[(size_t)i * 2 + 1]];
      ac[0] = ac[1] = (!acc_off && std::isfinite(rbx) && std::isfinite(rby) && rbx + rby <= r8sum) ? 1 : 0;
      if (ac[0] && (!std::isfinite(rf0) || !std::isfinite(rf1))) { cn.fb += 1; ac[0] = ac[1] = 0; }
      fin[0] = ac[0] ? rf0 : r8[0]; fin[1] = ac[1] ? rf1 : r8[1];
      if ((double)fin[0] + fin[1] > (double)r8[0] + r8[1]) cn.viol += 1;
      if (fin[0] > r8[0]) cn.up += 1;
      if (fin[1] > r8[1]) cn.up += 1;
      if (ac[0]) cn.a0 += 1; if (ac[1]) cn.a1 += 1;
      rx[i] = fin[0]; ry[i] = fin[1];
      if (g_rec_pol && ac[0] && (size_t)i < g_rp_ok.size()) { uint8_t okb = 0; for (int s = 0; s < 2; s++) { const int nt = SB.UQ.U[(size_t)ct[(size_t)i * 2 + s]].nth; if (nt < 1 || nt > 32) continue;
          std::memcpy(&g_rp_th[s][(size_t)i * 32], jthp[s] + 32 * k, sizeof(float) * (size_t)nt); okb |= (uint8_t)(1 << s); } g_rp_ok[(size_t)i] = okb; }
      const double sum = (double)fin[0] + (double)fin[1];
      err[i] = (float)std::min(sum, (double)DEAD);
      if (g_sicis_en && g_sicis_en->present) {
        const PostSicisOut& po = psout[k / 4096];
        const size_t q = k - po.lo;
        const bool gpu_ok = po.err && ((const uint32_t*)po.bad.contents)[q] == 0;
        if (gpu_ok) err[i] = ((const float*)po.err.contents)[q];
        else cn.pfb += 1;
        if (!gpu_ok || (g_post_sicis_check > 0 && k % (size_t)g_post_sicis_check == 0)) {
          float res48[NPT];
          const uint64_t jt0 = g_jprof ? mach_absolute_time() : 0;
          const bool pre = yfxb && ypre_have[k - lo];
          const bool jok = pre ? true : sicis_res48(E, SB.UQ, ct, K, k, SBf, bTV, res48, jptr ? jthp : nullptr, jptr ? jtvp : nullptr);
          const uint64_t jt1 = g_jprof ? mach_absolute_time() : 0;
          if (jok) {
            const float ce = pre ? ypre_ce[k - lo] : (float)std::min(sicis::score48(res48, *g_sicis_en), 1.0e30);
            if (pre && g_yfxfast == 2) { float r2[NPT]; if (sicis_res48(E, SB.UQ, ct, K, k, SBf, bTV, r2, jptr ? jthp : nullptr, jptr ? jtvp : nullptr)) { const float c2 = (float)std::min(sicis::score48(r2, *g_sicis_en), 1.0e30);
                g_yf_chk.fetch_add(1, std::memory_order_relaxed); if (std::memcmp(r2, ypre_res.data() + (size_t)(k - lo) * NPT, sizeof r2)) g_yf_bad.fetch_add(1, std::memory_order_relaxed);
                g_yf_schk.fetch_add(1, std::memory_order_relaxed); if (std::memcmp(&c2, &ce, 4)) g_yf_sbad.fetch_add(1, std::memory_order_relaxed); } }
            if (g_jprof && !pre) { const uint64_t jt2 = mach_absolute_time(); g_jp_res.fetch_add(jt1 - jt0, std::memory_order_relaxed); g_jp_sc.fetch_add(jt2 - jt1, std::memory_order_relaxed); g_jp_n.fetch_add(1, std::memory_order_relaxed); }
            if (gpu_ok) { cn.pchkn += 1; cn.pchk = std::max(cn.pchk, (double)std::fabs(ce - err[i])); }
            else err[i] = ce;
          }
        }
        if (!std::isfinite(err[i]) || err[i] <= -1.0e30f) err[i] = DEAD;
      }
      if (want_dump) { ORW[4 * k] = r8[0]; ORW[4 * k + 1] = r8[1]; ORW[4 * k + 2] = fin[0]; ORW[4 * k + 3] = fin[1]; OAC[2 * k] = (uint8_t)ac[0]; OAC[2 * k + 1] = (uint8_t)ac[1]; }
    }
  };
  bool h117b_on = false;
  if (g_h117b && (g_h117 & 2) && !g_h240b && h84_nb((size_t)nK, HB) > 1 && !HS.pbins.empty() && !g_h731 && !g_h660 && !g_h670 && !(g_h891 && g_h891_gpu_us > 0)) {


    h117b_on = true;
    const size_t nb = (size_t)h84_nb((size_t)nK, HB);
    auto wait_blk = [&](id<MTLCommandBuffer> cb, int which) {
      [cb waitUntilCompleted];
      if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H117 command buffer failed\n"); h182_dump(cb, "H117"); gpu_fail(); }
      gpu_t[which] += (cb.GPUEndTime - cb.GPUStartTime) * 1e3;
    };
    h301_mark("post fin: m4 enc grid+band");
    const bool kasc = std::is_sorted(K.begin(), K.end());
    if (g_kverb) printf("K_ESEL fin: nK %d nb %zu kasc %d hook %d freet %d\n", nK, nb, (int)kasc, (int)(bool)g_ovl_hook, g_freet);
    id<MTLCommandBuffer> cbG = h182_cb(E.cq);
    h131_enc_t([&]{ enc_grid(cbG); });
    [cbG commit];
    id<MTLCommandBuffer> cbB = h182_cb(E.cq);
    h131_enc_t([&]{ enc_band(cbB); });
    [cbB commit];
    h301_mark("post fin: m5 split+enc polish");
    std::vector<std::vector<P14Bin>> blk = h117b_split(HS.pbins, HB, nb);
    if (g_h117b_sl.size() < nb) g_h117b_sl.resize(nb);
    std::vector<id<MTLCommandBuffer>> cbP(nb, nil);
    for (size_t b = 0; b < nb; b++) {
      if (blk[b].empty()) continue;
      cbP[b] = h182_cb(E.cq);
      h131_enc_t([&]{ p14_encode_h91(E, cbP[b], blk[b], SBf, HS.ptok.b, HS.psp.b, bTV, 3.0f, 0.0f, g_h117b_sl[b]); });
      if (psOn) { const size_t lo = b * HB, hi = std::min((size_t)nK, (b + 1) * HB); h131_enc_t([&]{ post_sicis_enc(E, cbP[b], psCf, SB, ct, K, lo, hi, SBf, bTV, psout[b]); }); }
      [cbP[b] commit];
    }
    h301_mark("post fin: m6 wait grid+band");
    wait_blk(cbG, 0); wait_blk(cbB, 1);
    for (size_t b = 0; b < nb; b++) {
      h301_mark("post fin: m7 wait polish blk");
      if (cbP[b]) wait_blk(cbP[b], 2);
      h301_mark("post fin: m8 outputs blk");
      const size_t blo = b * HB, bhi = std::min((size_t)nK, (b + 1) * HB);
      if (g_pfr) { const size_t ns = (bhi - blo + PFG - 1) / PFG; PL(E).run((int)ns, [&](int c) { const size_t lo = blo + (size_t)c * PFG; outputs(lo, std::min(bhi, lo + PFG), (int)(lo / PFG)); }); }
      else outputs(blo, bhi, (int)b);

      if (g_ovl_hook && !g_freet && b + 1 < nb && (b + 1) * HB < (size_t)nK && kasc) { h301_mark("K: early select"); g_ovl_hook(K[(b + 1) * HB]); }
    }
  } else if ((g_h117 & 2) && !g_h240b) {
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    h131_enc_t([&]{ enc_grid(cb); enc_band(cb); enc_pol(cb); });
    if (psOn) { for (size_t b = 0; b < psout.size(); b++) { const size_t lo = b * 4096, hi = std::min((size_t)nK, (b + 1) * 4096); h131_enc_t([&]{ post_sicis_enc(E, cb, psCf, SB, ct, K, lo, hi, SBf, bTV, psout[b]); }); } }
    run_cb(cb, 2);
    if (g_h731) h731_wait(HS, gpu_t);
  } else if (g_h670 && g_h670_st > 0 && g_h240b && !HS.pbins.empty()) {

    h670_join();
    h670_on = true;
    H670Job& J = g_h670j;
    id<MTLCommandBuffer> cbG = h182_cb(E.cq); h131_enc_t([&]{ enc_grid(cbG); }); [cbG commit];
    id<MTLCommandBuffer> cbB = h182_cb(E.cq); h131_enc_t([&]{ enc_band(cbB); }); [cbB commit];
    J.sp = h240b_mark(E, HS);
    if (!g_h670q) g_h670q = [dev newCommandQueue];
    id<MTLCommandBuffer> cbP = h182_cb(g_h670q);
    h131_enc_t([&]{ enc_pol(cbP); });
    h240b_restore(HS, J.sp);
    J.SBf[0] = SBf[0]; J.SBf[1] = SBf[1]; J.tok = HS.ptok.b; J.psp = HS.psp.b; J.tv = bTV;
    J.cbG = cbG; J.cbB = cbB; J.cbP = cbP; J.cpu_ms = 0; J.cs_ok = false; J.step = g_h670_st;
    J.K.assign(K.begin(), K.end()); J.r8.resize((size_t)nK * 2);
    for (int k = 0; k < nK; k++) { J.r8[2 * (size_t)k] = SB.F[0].rfin[ct[(size_t)K[k] * 2]]; J.r8[2 * (size_t)k + 1] = SB.F[1].rfin[ct[(size_t)K[k] * 2 + 1]]; }
    J.pending = true;
    J.th = std::thread([] { H401S h401s_("H670 polish job"); tl_h675_job = g_h675 != 0;
      if (g_h675) pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0);
      H670Job& J = g_h670j;
      [J.cbB waitUntilCompleted];
      [J.cbP commit];
      h240b_cpu(J.sp, J.SBf, J.tok, J.psp, J.tv, J.cpu_ms);
      [J.cbP waitUntilCompleted]; });
    if (g_h675) { g_h53cur = (g_h53cur == &g_h53sb[0]) ? &g_h53sb[1] : &g_h53sb[0]; g_h53p = g_h53cur; }
  } else if ((g_h257 & 1) && g_h240b && (h898_ovl || !HS.pbins.empty())) {


    auto wait_cb = [&](id<MTLCommandBuffer> cb, int which) {
      [cb waitUntilCompleted];
      if (cb.status != MTLCommandBufferStatusCompleted) { fprintf(stderr, "H117 command buffer failed\n"); h182_dump(cb, "H117"); gpu_fail(); }
      gpu_t[which] += (cb.GPUEndTime - cb.GPUStartTime) * 1e3;
    };
    id<MTLCommandBuffer> cbG = h182_cb(E.cq);
    h131_enc_t([&]{ enc_grid(cbG); });
    [cbG commit];
    id<MTLCommandBuffer> cbB = h182_cb(E.cq);
    h131_enc_t([&]{ enc_band(cbB); });
    [cbB commit];
    if (h898_ovl) { h301_mark("post fin: h898 prep"); h898_flush(E, SB, HS, T); h301_mark("post fin: polish mark+enc"); }
    if (HS.pbins.empty()) { wait_cb(cbG, 0); wait_cb(cbB, 1); } else {
    const auto tq = clk::now();
    H240bSplit sp = h240b_mark(E, HS);
    if (g_h894 & 16) h301_mark("pol: mark done");
    id<MTLCommandBuffer> cbP = h182_cb(E.cq);
    h131_enc_t([&]{ enc_pol(cbP); });
    h240b_restore(HS, sp);
    const double g0 = gpu_t[2];
    const auto tc = clk::now();
    if (g_h894 & 16) { h301_mark("pol: enc end"); [cbP addScheduledHandler:^(id<MTLCommandBuffer>) { if (g_h892 && g_h892_ring) { const size_t i = g_h892_n.fetch_add(1, std::memory_order_relaxed); const double t = h301_now(); if (i < g_h892_cap) g_h892_ring[i] = {6, 0, "cbP scheduled", t, t}; } }]; }
    [cbP commit];
    if (g_h894 & 16) h301_mark("pol: committed");
    wait_cb(cbG, 0); wait_cb(cbB, 1);
    if (g_h731) h731_wait(HS, gpu_t);
    double cpu_ms = 0;
    FitHandle cpuC = fw_submit_nospan265([&]{ h240b_cpu(sp, SBf, HS.ptok.b, HS.psp.b, bTV, cpu_ms); });
    wait_cb(cbP, 2);
    wall_t[2] += ms_since(tc);
    cpuC.join();
    const double gpu_ms = gpu_t[2] - g0, tail = ms_since(tq);
    const size_t ngpu = sp.ntot - sp.items.size();
    h240b_observe(sp.ntot, ngpu, gpu_ms, sp.items.size(), cpu_ms);
    printf("H240B: f_cpu %.3f cpu ms %.3f gpu ms %.3f tail wall %.3f (trees cpu %zu gpu %zu)\n", sp.f, cpu_ms, gpu_ms, tail, sp.items.size(), ngpu);
    }
  } else {
    { id<MTLCommandBuffer> cb = h182_cb(E.cq); h131_enc_t([&]{ enc_grid(cb); }); h301_mark("post fin: grid wait"); run_cb(cb, 0); }
    { h301_mark("post fin: band enc"); id<MTLCommandBuffer> cb = h182_cb(E.cq); h131_enc_t([&]{ enc_band(cb); }); h301_mark("post fin: band wait"); run_cb(cb, 1); }
    if (g_h731) h731_wait(HS, gpu_t);
    if (g_h660 && !HS.pbins.empty()) { if (g_h660 == 2) h660_burn(); }
    else if (!HS.pbins.empty() && g_h240b) {
      h301_mark("post fin: polish mark+enc");
      const auto tq = clk::now();
      H240bSplit sp = h240b_mark(E, HS);
      id<MTLCommandBuffer> cb = h182_cb(E.cq);
      h131_enc_t([&]{ enc_pol(cb); });
      h240b_restore(HS, sp);
      const double g0 = gpu_t[2];
      double cpu_ms = 0;
      FitHandle cpuC = fw_submit_nospan265([&]{ h240b_cpu(sp, SBf, HS.ptok.b, HS.psp.b, bTV, cpu_ms); });
      h301_mark("post fin: polish GPU wait");
      run_cb(cb, 2);
      h301_mark("post fin: polish CPU join");
      cpuC.join();
      h301_mark("post fin: outputs");
      const double gpu_ms = gpu_t[2] - g0, tail = ms_since(tq);
      const size_t ngpu = sp.ntot - sp.items.size();
      h240b_observe(sp.ntot, ngpu, gpu_ms, sp.items.size(), cpu_ms);
      printf("H240B: f_cpu %.3f cpu ms %.3f gpu ms %.3f tail wall %.3f (trees cpu %zu gpu %zu)\n", sp.f, cpu_ms, gpu_ms, tail, sp.items.size(), ngpu);
    }
    if (!HS.pbins.empty() && !g_h240b) { id<MTLCommandBuffer> cb = h182_cb(E.cq); h131_enc_t([&]{ enc_pol(cb); }); run_cb(cb, 2); }

  }
  h301_mark("post fin: m9 tail");
  T.p_grid += gpu_t[0]; T.p_band += gpu_t[1]; T.p_pb += gpu_t[2]; T.p_pb_wall += wall_t[2];
  auto to = clk::now();
  if (!h117b_on) h84_run(E, (size_t)nK, 4096, outputs);
  for (const OutCnt& cn : oc) { T.p_viol += cn.viol; T.p_side_up += cn.up; T.p_acc0 += cn.a0; T.p_acc1 += cn.a1; T.p_fb += cn.fb; }
  if (psOn && g_post_sicis_check > 0) { double pc = 0, pn = 0, pf = 0; for (const OutCnt& cn : oc) { pc = std::max(pc, cn.pchk); pn += cn.pchkn; pf += cn.pfb; }
    printf("POST_SICIS_CHECK call %d: n=%zu checked=%.0f max|gpu-cpu|=%.3g cpu-fallback=%.0f\n", HS.call, (size_t)nK, pn, pc, pf); }
  if (g_dumpf) for (int k = 0; k < nK; k++)
    fprintf(g_dumpf, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%d\n", HS.call, K[k], ORW[4 * k], ORW[4 * k + 1], ORW[4 * k + 2], ORW[4 * k + 3], OAC[2 * k], OAC[2 * k + 1]);
  T.p_cpu += ms_since(to);
  T.p_ebmb += g_ebmb_call; g_ebmb_call = 0;
  h894_pregrow(E.dev, g_h117_sl.pe);
}


static void post16_h53(Engine& E, PipeStep& PS, StepBatch& SB, const std::vector<int>& cand, float* err, float* rx, float* ry, Times& T) {
  post16_h53_prep(E, PS, SB, cand, T);
  post16_h53_fin(E, PS, SB, err, rx, ry, T);
}


static void scan_begin_p_h215(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;
  const int M = (int)UQ.U.size();
  B.M = M;
  Persist5& P = *B.P;
  id<MTLDevice> dev = E.dev;
  auto t0 = clk::now();
  Tick tk;
  if (!P.bct) { const float ctv[16] = {1,0,0,0, 0,1,1,0, 0,1,-1,0, 0,-1,-1,0}; P.bct = E.B(ctv, sizeof ctv); P.bsob = E.B(SOB_U, sizeof(SOB_U)); }
  for (int s = 0; s < 2; s++) if (!P.S[s].yb) P.S[s].yb = E.B(E.tg[s].yf, sizeof E.tg[s].yf);
  const int NCH = 4 * (int)WK(E).size();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s];
    if ((int)X.listed.size() < M) X.listed.resize(M, 0);
    for (int j : X.ls) X.listed[j] = 0;
    X.ls.assign(B.need[s].begin(), B.need[s].end()); std::sort(X.ls.begin(), X.ls.end());
    X.ls.erase(std::unique(X.ls.begin(), X.ls.end()), X.ls.end());
    for (int j : X.ls) X.listed[j] = 1;
    const int L = (int)X.ls.size();
    auto lspan = [L, NCH](int c, int& a, int& z) { a = (int)((long)L * c / NCH); z = (int)((long)L * (c + 1) / NCH); };
    X.sok.resize(M); X.ninst.resize(M); X.toff.resize(M); X.spoff.resize(M); X.thoff.resize(M); X.dkey.resize(M);
    T.h_s_alloc += tk.lap();

    std::vector<uint32_t> bti(NCH, 0), bsi(NCH, 0), bhi(NCH, 0), bns(NCH, 0);
    std::vector<int> bq(NCH, 0); std::vector<double> bfl(NCH, 0.0);
    PL(E).run(NCH, [&](int c) {
      int a, z; lspan(c, a, z);
      uint32_t ti = 0, si = 0, hi = 0, ns = 0; int q = 0; double fl = 0;
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (te.s_ok) { ti += (uint32_t)(te.nsp() / 4); ns += 1; q = std::max(q, UQ.U[j].Q); } else fl += 1;
        si += (uint32_t)tr.sp.size(); hi += (uint32_t)tr.th.size();
      }
      bti[c] = ti; bsi[c] = si; bhi[c] = hi; bns[c] = ns; bq[c] = q; bfl[c] = fl;
    });
    std::vector<uint32_t> oti(NCH), osi(NCH), ohi(NCH), osq(NCH);
    uint32_t TI = 0, SI = 0, HI = 0, NS = 0; int QM = 0;
    for (int c = 0; c < NCH; c++) {
      oti[c] = TI; osi[c] = SI; ohi[c] = HI; osq[c] = NS;
      TI += bti[c]; SI += bsi[c]; HI += bhi[c]; NS += bns[c]; QM = std::max(QM, bq[c]); T.scan_fail += bfl[c];
    }
    X.ti = TI; X.si = SI; X.hi = HI;
    X.seq.resize(NS);

    PL(E).run(NCH, [&](int c) {
      int a, z; lspan(c, a, z);
      uint32_t ti = oti[c], si = osi[c], hi = ohi[c], sq = osq[c];
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        X.toff[j] = ti; X.spoff[j] = si; X.thoff[j] = hi;
        X.sok[j] = 0; X.ninst[j] = 0; X.dkey[j] = 0;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (te.s_ok) {
          X.sok[j] = 1; X.ninst[j] = (uint32_t)(te.nsp() / 4); ti += X.ninst[j];
          const uint32_t cls = (tr.nS <= 4 ? 0u : tr.nS <= 8 ? 1u : tr.nS <= 16 ? 2u : 3u) + (tr.nth == 0 ? 4u : 0u);
          X.dkey[j] = (uint32_t)(((size_t)tr.Q * 41 + (size_t)tr.nslot) * 8 + cls);
          X.seq[sq++] = (uint32_t)j;
        }
        si += (uint32_t)tr.sp.size(); hi += (uint32_t)tr.th.size();
      }
    });
    const size_t DN = (size_t)(QM + 1) * 41 * 8;
    T.h_s_pre += tk.lap();
    sort_by_key5_p(E, X.seq, X.dkey, DN, X.start, X.sorted);
    T.h_s_sort += tk.lap();
    X.sbins.clear(); size_t cur = 0;
    for (size_t dk = 0; dk < DN; dk++) {
      const size_t c = X.start[dk + 1] - X.start[dk]; if (!c) continue;
      SBin5 b; b.Q = (int)(dk / (41 * 8)); b.ns = (int)((dk / 8) % 41); const int cls = (int)(dk % 8);
      b.s2 = cls >= 4; b.nc = 4 << (cls & 3); b.off = cur; b.cnt = c; b.sst = X.start[dk];
      X.sbins.push_back(b); cur = al5(cur + c, BIN_ALIGN);
    }
    X.bl_words = cur;
    T.h_s_bins += tk.lap();
    h900_ens(X.tok, "X.tok", dev, ((size_t)X.ti * 4 + 8) * 4, MSH);
    h900_ens(X.hdr, "X.hdr", dev, (size_t)M * 16, MSH);
    h900_ens(X.sp, "X.sp", dev, (size_t)X.si * 2, MSH);
    h900_ens(X.th, "X.th", dev, (size_t)X.hi * 4, MSH);
    h900_ens(X.seg, "X.seg", dev, (size_t)M * 144, MSH);
    h900_ens(X.out, "X.out", dev, (size_t)M * 64, MSH);
    h900_ens(X.bl, "X.bl", dev, X.bl_words * 4, MSH);
    T.h_s_alloc += tk.lap();
    uint32_t* tok = (uint32_t*)X.tok.b.contents; uint32_t* hdr = (uint32_t*)X.hdr.b.contents; uint16_t* sp = (uint16_t*)X.sp.b.contents;
    float* th = (float*)X.th.b.contents; uint32_t* seg = (uint32_t*)X.seg.b.contents; float* out = (float*)X.out.b.contents;
    uint32_t* bl = (uint32_t*)X.bl.b.contents;
    PL(E).run((int)X.sbins.size(), [&](int bi) {
      const SBin5& b = X.sbins[bi]; std::memcpy(bl + b.off, X.sorted.data() + b.sst, b.cnt * 4);
    });
    T.h_s_bins += tk.lap();
    PL(E).run(NCH, [&](int c) {
      int a, z; lspan(c, a, z);
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        uint32_t* h = hdr + 4 * (size_t)j; float* ob = out + 16 * (size_t)j; uint32_t* sg = seg + 36 * (size_t)j;
        h[0] = X.toff[j]; h[2] = X.spoff[j];
        for (int q = 0; q < 16; q++) ob[q] = INFINITY;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.ninst[j] | ((uint32_t)tr.nS << 16);
        h[3] = (uint32_t)tr.nth | (X.thoff[j] << 8);
        if (te.s_ok) { std::memcpy(sg, te.sg, 144); std::memcpy(tok + 4 * (size_t)X.toff[j], te.psp(), te.nsp() * 4); }
        else std::memset(sg, 0, 144);
        if (!tr.sp.empty()) std::memcpy(sp + X.spoff[j], tr.sp.data(), tr.sp.size() * 2);
        if (!tr.th.empty()) std::memcpy(th + X.thoff[j], tr.th.data(), tr.th.size() * 4);
      }
    });
    for (int z = 0; z < 8; z++) tok[4 * (size_t)X.ti + z] = 0;
    T.h_s_fill += tk.lap();
  }
  T.scan_enc += ms_since(t0);


  Tick tc;
  size_t nsb = 0; for (int s = 0; s < 2; s++) nsb += P.S[s].sbins.size();
  if (nsb) {
    h900_ens(P.scr, "P.scr", dev, (size_t)M * 16 * 16 * 4, MPV);
    size_t h139_tot = 0;
    for (int s = 0; s < 2; s++) { Side5& X = P.S[s]; X.h139n = 0; X.h139t = 0;
      if (g_h139) for (const SBin5& b : X.sbins) { const size_t t = h139_tgb(b);
        if (t > 32768) { h139_tot += ((t * b.cnt) + 255) & ~(size_t)255; X.h139n++; X.h139t += b.cnt; } } }
    if (h139_tot) h900_ens(P.h139scr, "P.h139scr", dev, h139_tot, MPV);
    if (g_h139 == 1) for (int s = 0; s < 2; s++) if (P.S[s].h139n) h900_ens(P.S[s].h139o, "P.S[s].h139o", dev, (size_t)M * 64, MSH);
    size_t h139_off = 0;
    id<MTLCommandBuffer> cb = h182_cb(E.cq); id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (int s = 0; s < 2; s++) {
      Side5& X = P.S[s];
      for (const SBin5& b : X.sbins) {
        Tick tp;
        const size_t tgb = h139_tgb(b); const bool over = tgb > 32768; const bool dpath = over && g_h139 == 2;
        id<MTLComputePipelineState> ps = dpath ? h139_ps(E, b.s2 ? 1 : 0, b.Q) : h93_scan_ps(E, b.s2 ? 1 : 0, b.Q);
        T.h_s_pipe += tp.lap();
        if (!ps) { fprintf(stderr, "missing scan kernel Q%d\n", b.Q); exit(1); }
        [en setComputePipelineState:ps];
        id<MTLBuffer> bs[9] = {X.tok.b, X.hdr.b, X.sp.b, X.th.b, P.bsob, X.yb, X.out.b, X.bl.b, X.seg.b};
        for (int k = 0; k < 9; k++) [en setBuffer:bs[k] offset:(k == 7 ? b.off * 4 : 0) atIndex:k];
        uint32_t ns = (uint32_t)b.ns, nc = (uint32_t)b.nc; float tau = b.Q * 9.0f * ldexpf(1.0f, -24);
        [en setBuffer:P.bct offset:0 atIndex:11]; [en setBuffer:P.scr.b offset:0 atIndex:12];
        [en setBytes:&tau length:4 atIndex:13]; [en setBytes:&ns length:4 atIndex:9]; [en setBytes:&nc length:4 atIndex:10];
        uint32_t sstr = (uint32_t)(((size_t)ns * 96 + 2 * nc + 59) & ~(size_t)3); [en setBytes:&sstr length:4 atIndex:14];
        const int thr = b.s2 ? 96 / 3 : (b.Q <= 5 ? 768 / 6 : 768 / 3);
        if (dpath) {
          const uint32_t stride = (uint32_t)(tgb / 4); const size_t reg = h139_off; h139_off += ((tgb * b.cnt) + 255) & ~(size_t)255;
          [en setBuffer:P.h139scr.b offset:reg atIndex:15]; [en setBytes:&stride length:4 atIndex:16];
          [en setThreadgroupMemoryLength:0 atIndex:0];
          [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)];
        } else {
          if (tgb > 32768) g_h139_viol += b.cnt;
          [en setThreadgroupMemoryLength:tgb atIndex:0];
          [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)];
          if (over && g_h139 == 1) {
            const uint32_t stride = (uint32_t)(tgb / 4); const size_t reg = h139_off; h139_off += ((tgb * b.cnt) + 255) & ~(size_t)255;
            [en setComputePipelineState:h139_ps(E, b.s2 ? 1 : 0, b.Q)];
            [en setBuffer:X.h139o.b offset:0 atIndex:6];
            [en setBuffer:P.h139scr.b offset:reg atIndex:15]; [en setBytes:&stride length:4 atIndex:16];
            [en setThreadgroupMemoryLength:0 atIndex:0];
            [en dispatchThreadgroups:MTLSizeMake(b.cnt, 1, 1) threadsPerThreadgroup:MTLSizeMake(thr, 1, 1)];
          }
        }
      }
    }
    [en endEncoding]; [cb commit];
    B.cbS = cb;
    if (g_h139 && (P.S[0].h139n || P.S[1].h139n))
      printf("H139 commit mode %d: over-limit bins %zu (TGs %zu), TG-memory violations so far %llu\n", g_h139,
             P.S[0].h139n + P.S[1].h139n, P.S[0].h139t + P.S[1].h139t, (unsigned long long)g_h139_viol);
  }
  T.h_s_cmd += tc.lap();
}

static void gn_begin_p_h215(Engine& E, PipeStep& PS, Batch& B, Times& T) {
  const Unique& UQ = PS.SB->UQ;
  const int M = B.M;
  if ((!g_h217 && (int)PS.SB->UQ.U.size() < M) || (int)B.P->S[0].listed.size() < M || (int)B.P->S[1].listed.size() < M) { fprintf(stderr, "pipeline: slot arrays do not match chunk %d\n", B.b0); exit(1); }
  Persist5& P = *B.P; id<MTLDevice> dev = E.dev;
  const int NCH = 4 * (int)WK(E).size();
  auto t0 = clk::now();
  Tick tk;
  P.gbins.clear();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    X.gok.resize(M); X.gninst.resize(M); X.gper.resize(M); X.gtoff.resize(M); X.gspoff.resize(M); X.keyl.resize(M); X.dkey.resize(M);
    const int L = (int)X.ls.size();
    auto lspan = [L, NCH](int c, int& a, int& z) { a = (int)((long)L * c / NCH); z = (int)((long)L * (c + 1) / NCH); };
    T.h_g_alloc += tk.lap();

    std::vector<uint32_t> bgi(NCH, 0), bgs(NCH, 0), bns(NCH, 0);
    std::vector<int> bqg(NCH, 0), bcg(NCH, 0), blm(NCH, 0);
    std::vector<double> bal(NCH, 0.0), bkm(NCH, 0.0);
    PL(E).run(NCH, [&](int c) {
      int a, z; lspan(c, a, z);
      uint32_t gi = 0, gs = 0, ns = 0; int qg = 0, cg = 0, lm = 0; double al = 0, km = 0;
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        X.gok[j] = 0; X.gninst[j] = 0; X.gper[j] = 0;
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        if (tr.nth >= 1 && tr.nth > GN_TH_MAX) km += 1;
        if (te.g_ok && isfinite(o.smin[j])) {
          X.gok[j] = 1; o.is_g[j] = 1; al += 1;
          X.gninst[j] = (uint32_t)(te.ngp() / 4); X.gper[j] = te.gper; gi += X.gninst[j];
          ns += 1; qg = std::max(qg, tr.Q); cg = std::max(cg, (int)((X.gper[j] + 63) / 64)); lm = std::max(lm, (int)h313a::nprog(tr));
        }
        gs += (uint32_t)tr.sp.size();
      }
      bgi[c] = gi; bgs[c] = gs; bns[c] = ns; bqg[c] = qg; bcg[c] = cg; blm[c] = lm; bal[c] = al; bkm[c] = km;
    });
    std::vector<uint32_t> ogt(NCH), ogs(NCH), osq(NCH);
    uint32_t GTI = 0, GSI = 0, NS = 0; int Qg = 0, Cg = 0, LM = 0;
    for (int c = 0; c < NCH; c++) {
      ogt[c] = GTI; ogs[c] = GSI; osq[c] = NS;
      GTI += bgi[c]; GSI += bgs[c]; NS += bns[c];
      Qg = std::max(Qg, bqg[c]); Cg = std::max(Cg, bcg[c]); LM = std::max(LM, blm[c]);
      T.g_alive += bal[c]; T.gn_kmiss += bkm[c];
    }
    X.gti = GTI; X.gsi = GSI;
    X.seq.resize(NS);
    if (Qg >= 100) { fprintf(stderr, "GN: Q %d >= 100 breaks the bin key\n", Qg); exit(1); }

    PL(E).run(NCH, [&](int c) {
      int a, z; lspan(c, a, z);
      uint32_t gt = ogt[c], gs = ogs[c], sq = osq[c];
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        X.gtoff[j] = gt; X.gspoff[j] = gs; X.keyl[j] = 0; X.dkey[j] = 0;
        const Tree& tr = UQ.U[j];
        if (X.gok[j]) { X.seq[sq++] = (uint32_t)j; X.keyl[j] = (uint32_t)(LM - (int)h313a::nprog(tr)); gt += X.gninst[j]; }
        gs += (uint32_t)tr.sp.size();
      }
    });
    T.h_g_pre += tk.lap();
    sort_by_key5_p(E, X.seq, X.keyl, (size_t)LM + 1, X.start, X.seq2);
    T.h_g_sort += tk.lap();
    const size_t DNg = (size_t)(Cg + 1) * (size_t)(Qg + 1) * 25;
    PL(E).run(NCH, [&](int c) {
      const size_t n2 = X.seq2.size(), a = n2 * c / NCH, z = n2 * (c + 1) / NCH;
      for (size_t i = a; i < z; i++) {
        const uint32_t j = X.seq2[i]; const Tree& tr = UQ.U[j]; const size_t cls = (X.gper[j] + 63) / 64;
        X.dkey[j] = (uint32_t)((cls * (size_t)(Qg + 1) + (size_t)tr.Q) * 25 + (size_t)tr.nth);
      }
    });
    T.h_g_plan += tk.lap();
    sort_by_key5_p(E, X.seq2, X.dkey, DNg, X.start, X.sorted);
    T.h_g_sort += tk.lap();

    std::vector<size_t> ne; for (size_t dk = 0; dk < DNg; dk++) if (X.start[dk + 1] > X.start[dk]) ne.push_back(dk);
    std::vector<size_t> perb(ne.size(), 0); std::vector<size_t> perb2(ne.size(), 0); std::vector<size_t> perb3(ne.size(), 0);
    PL(E).run(NCH, [&](int c) {
      const size_t a = ne.size() * c / NCH, z = ne.size() * (c + 1) / NCH;
      for (size_t bi = a; bi < z; bi++) {
        size_t per = 0;
        for (size_t k = X.start[ne[bi]]; k < X.start[ne[bi] + 1]; k++) per = std::max<size_t>(per, X.gper[X.sorted[k]]);
        size_t pen = 0, pen32 = 0; for (size_t k = X.start[ne[bi]]; k < X.start[ne[bi] + 1]; k++) { const size_t p32_ = h141_pe(UQ.enc[X.sorted[k]], X.gper[X.sorted[k]]); pen32 = std::max<size_t>(pen32, p32_); pen = std::max<size_t>(pen, h332_pe(UQ.enc[X.sorted[k]], (uint32_t)p32_)); }
        perb[bi] = per; perb2[bi] = pen; perb3[bi] = pen32;
      }
    });
    for (size_t bi = 0; bi < ne.size(); bi++) {
      const size_t dk = ne[bi];
      const size_t c = X.start[dk + 1] - X.start[dk];
      const size_t sst = X.start[dk];
      const int Tk = (int)(dk % 25); const size_t rem = dk / 25; const int Qk = (int)(rem % (size_t)(Qg + 1));
      Tick tp;
      id<MTLComputePipelineState> ps = h93_gn_ps(E, Qk, Tk);
      T.h_g_pipe += tp.lap();
      if (!ps) { T.gn_kmiss += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      const size_t per = perb[bi]; const size_t pen = perb2[bi]; const size_t pen32 = perb3[bi];
      const int NPk = std::max(4, (Tk + 3) / 4 * 4), QPk = std::max(4, (Qk + 3) / 4 * 4); const size_t WSk = (4 * NPk + 3 * QPk + 6) / 4 * 4;
      const size_t pers = WSk + 2 * per; const size_t shb = (size_t)(Tk <= 4 ? 1 + (Tk + 1) / 2 : 3) * 48;
      if (pers * 4 > 32768) { T.gn_drop += (double)c; for (size_t k = 0; k < c; k++) o.is_g[X.sorted[sst + k]] = 0; continue; }
      const size_t persN = WSk + 2 * pen;
      const int cap141 = (int)ps.maxTotalThreadsPerThreadgroup / 32;
      int tpg = h141_pick(shb, persN, c, cap141);
      h141_note(shb, pers, persN, c, cap141, tpg);
      const size_t xf713 = h713_gn(E, Qk, Tk, shb, persN, ps, tpg);
      GBin5 g{}; g.side = s; g.Q = Qk; g.nth = Tk; g.tpg = tpg; g.cnt = c;
      g.nt = (c + (size_t)tpg - 1) / (size_t)tpg * (size_t)tpg;
      g.ntg = xf713 ? g.nt : g.nt / (size_t)tpg;
      g.tgb = ((shb + persN * tpg + xf713) * 4 + 15) & ~(size_t)15; g.per = pen; g.per32 = pen32; g.sst = sst; g.ps = ps;
      g.eb_bytes = g.ntg * (size_t)tpg * 2 * (size_t)(Qk + 1) * 64 * 4;
      P.gbins.push_back(g);
    }
    T.h_g_plan += tk.lap();
  }

  if (h94_mode() >= 1) {
    std::vector<size_t> nat, tgb, cnt, cap; size_t before = 0, after = 0, ntr = 0;
    for (const GBin5& g : P.gbins) { nat.push_back(g.ntg); tgb.push_back(g.tgb); cnt.push_back(g.cnt); before += g.ntg; ntr += g.cnt; }
    cap.assign(nat.size(), 0); h94_cap(cap, nat, tgb, cnt);
    for (size_t b = 0; b < P.gbins.size(); b++) { GBin5& g = P.gbins[b]; g.ntg = cap[b]; g.eb_bytes = g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4; after += g.ntg; }
    fprintf(stderr, "H94 GN dispatch: bins %zu trees %zu TGs before %zu after %zu\n", P.gbins.size(), ntr, before, after);
  }
  size_t tlc[2] = {0, 0}, ebc = 0; uint32_t NB = 0;
  for (GBin5& g : P.gbins) {
    g.wix = NB++;
    g.tl_off = tlc[g.side]; tlc[g.side] = al5(tlc[g.side] + 1 + g.nt, BIN_ALIGN);
    g.eb_off = ebc; ebc = al5(ebc + g.eb_bytes + ebt_bytes(g.Q, g.nth, g.ntg, g.tpg), 256);
  }
  T.h_g_plan += tk.lap();
  for (int s = 0; s < 2; s++) {
    Side5& X = P.S[s]; FitOut& o = PS.SB->F[s];
    h900_ens(X.gtok, "X.gtok", dev, ((size_t)X.gti * 4 + 8) * 4, MSH);
    h900_ens(X.ghdr, "X.ghdr", dev, (size_t)M * 16, MSH);
    h900_ens(X.gsp, "X.gsp", dev, (size_t)X.gsi * 2, MSH);
    h900_ens(X.gseg, "X.gseg", dev, (size_t)M * 144, MSH);
    h900_ens(X.gth, "X.gth", dev, (size_t)M * 128, MSH);
    h900_ens(X.grs, "X.grs", dev, (size_t)M * 8, MSH);
    h900_ens(X.gst, "X.gst", dev, (size_t)M * 4, MSH);
    h900_ens(X.gai, "X.gai", dev, (size_t)M * 64, MSH);
    h900_ens(X.gbl, "X.gbl", dev, tlc[s] * 4, MSH);
    T.h_g_alloc += tk.lap();
    uint32_t* gtok = (uint32_t*)X.gtok.b.contents; uint32_t* ghdr = (uint32_t*)X.ghdr.b.contents; uint16_t* gsp = (uint16_t*)X.gsp.b.contents;
    uint32_t* gseg = (uint32_t*)X.gseg.b.contents; float* gth = (float*)X.gth.b.contents;
    const int L2 = (int)X.ls.size(); auto lspan2 = [L2, NCH](int c, int& a, int& z) { a = (int)((long)L2 * c / NCH); z = (int)((long)L2 * (c + 1) / NCH); };
    PL(E).run(NCH, [&](int c) {
      int a, z; lspan2(c, a, z);
      for (int k = a; k < z; k++) { const int j = X.ls[k];
        uint32_t* h = ghdr + 4 * (size_t)j; float* th = gth + 32 * (size_t)j; uint32_t* sg = gseg + 36 * (size_t)j;
        for (int q = 0; q < 32; q++) th[q] = 0.0f;
        h[0] = X.gtoff[j]; h[2] = X.gspoff[j];
        const Tree& tr = UQ.U[j]; const TreeEnc& te = UQ.enc[j];
        h[1] = X.gninst[j] | ((uint32_t)tr.nS << 16);
        if (X.gok[j]) {
          h[3] = te.ghw; std::memcpy(sg, te.gg, 144);
          std::memcpy(gtok + 4 * (size_t)X.gtoff[j], te.pgp(), te.ngp() * 4);
          double th2[32]; start_theta(tr.th, tr.nth, PS.SB->F[s].sarg[j], th2);
          for (int q = 0; q < tr.nth; q++) th[q] = (float)th2[q];
        } else { h[3] = 0; std::memset(sg, 0, 144); }
        if (!tr.sp.empty()) std::memcpy(gsp + X.gspoff[j], tr.sp.data(), tr.sp.size() * 2);
      }
    });
    for (int z = 0; z < 8; z++) gtok[4 * (size_t)X.gti + z] = 0;
    T.h_g_fill += tk.lap();
    uint32_t* gbl = (uint32_t*)X.gbl.b.contents;
    std::vector<const GBin5*> sb; for (const GBin5& g : P.gbins) if (g.side == s) sb.push_back(&g);
    if (!sb.empty()) PL(E).run((int)sb.size(), [&](int i) {
      const GBin5& g = *sb[i]; uint32_t* tl = gbl + g.tl_off; tl[0] = (uint32_t)g.per;
      for (size_t k = 0; k < g.cnt; k++) tl[1 + k] = X.sorted[g.sst + k];
      for (size_t k = g.cnt; k < g.nt; k++) tl[1 + k] = 0xFFFFFFFFu;
    });
    T.h_g_fill += tk.lap();
  }
  h900_ens(P.eb, "P.eb", dev, ebc, MPV);
  h900_ens(P.ctr, "P.ctr", dev, 4 * (size_t)NB + 16, MSH);
  T.h_g_alloc += tk.lap();
  T.gn_enc += ms_since(t0);


  Tick tc;
  if (NB) {
    std::memset(P.ctr.b.contents, 0, P.ctr.cap);
    id<MTLCommandBuffer> cb = h182_cb(E.cq);
    id<MTLComputeCommandEncoder> en = tl_enc(cb, MTLDispatchTypeConcurrent);
    for (size_t gj779 = 0; gj779 < P.gbins.size(); gj779++) { const GBin5& g = P.gbins[g_h779 ? P.gbins.size() - 1 - gj779 : gj779];
      Side5& X = P.S[g.side];
      if (h182_on()) [en pushDebugGroup:h182_gname(g)];
      [en setComputePipelineState:g.ps];
      id<MTLBuffer> bs[9] = {X.gtok.b, X.ghdr.b, X.gsp.b, X.gseg.b, X.yb, X.gth.b, X.grs.b, X.gst.b, X.gbl.b};
      for (int k = 0; k < 9; k++) [en setBuffer:bs[k] offset:(k == 8 ? g.tl_off * 4 : 0) atIndex:k];
      float tf[4] = {g.Q * 9.0f * ldexpf(1.0f, -24), E.tg[g.side].floor_, (g_freet ? 0.0f : 3.0f), 0.0f};
      [en setBytes:tf length:16 atIndex:9];
      [en setBuffer:X.gai.b offset:0 atIndex:10];
      [en setBuffer:P.eb.b offset:g.eb_off atIndex:11];
      if (g_h711) [en setBuffer:P.eb.b offset:g.eb_off + g.ntg * (size_t)g.tpg * 2 * (size_t)(g.Q + 1) * 64 * 4 atIndex:15];
      uint32_t wq[2] = {g.wix, (uint32_t)g.nt}; [en setBytes:wq length:8 atIndex:12];
      [en setBuffer:P.ctr.b offset:0 atIndex:13];
      [en setThreadgroupMemoryLength:g.tgb atIndex:0];
      if (h150c_on()) { g_h150.gnDisp++; g_h150.gnLaunch += (long)(g.ntg * (size_t)g.tpg); g_h150.gnOcc += (long)g.cnt; }
      [en dispatchThreadgroups:MTLSizeMake(g.ntg, 1, 1) threadsPerThreadgroup:MTLSizeMake(32 * (size_t)g.tpg, 1, 1)];
      if (h182_on()) [en popDebugGroup];
    }
    [en endEncoding]; [cb commit];
    B.cbG = cb;
  }
  T.h_g_cmd += tc.lap();
}


static void h265_bal_update(const PipeStep& PS, const CpuSpan265& R) {
  if (!g_h242d_bal) return;
  double gs = 1e300, ge = -1e300;
  for (const auto& q : PS.iv) { gs = std::min(gs, q.first); ge = std::max(ge, q.second); }
  const double eg = PS.iv.empty() ? 0.0 : (ge - gs) * 1e3;
  const double ec = R.any ? std::chrono::duration<double, std::milli>(R.e - R.s).count() : 0.0;
  if (ec > 0 && eg > 0) { const double L = std::log(ec / eg); g_bal_s = std::min(20.0, std::max(0.05, g_bal_s * std::exp(-0.5 * L))); }
  printf("H242D step %2d: s %.3f | CPU-fit end %.1f ms | GPU-busy end %.1f ms | next CPU shares = clamp(f*s, 0, 0.5)\n", ++g_h242d_steps, g_bal_s, ec, eg);
}


static const int g_tiny = getenv("TINY_FIX") ? atoi(getenv("TINY_FIX")) : 0;
static const int g_tiny_chk = getenv("TINY_CHK") ? atoi(getenv("TINY_CHK")) : 0;
static std::vector<uint8_t> g_tdone;
static const int g_tiny_early = getenv("TINY_EARLY") ? atoi(getenv("TINY_EARLY")) : 0;
static void tiny_score(Engine& E, PipeStep& PS, StepBatch& SB, float* err, int n, const std::vector<uint8_t>& admit, int lo_ = 0, int hi_ = -1) {
  if (hi_ < 0) hi_ = n; const bool fin = lo_ == 0 && hi_ == n;
  if (!g_tiny || g_yfx || !g_sf_axis || g_freet || !g_sicis_en || !g_sicis_en->present || NPT != 48) return;
  const auto tt0 = clk::now();
  static float tgf[48]; for (int p = 0; p < 48; p++) tgf[p] = (float)((double)p / 47.0);
  const std::vector<int>& ct = SB.ct; if ((int)ct.size() < 2 * n) return;
  std::vector<int> todo, chk; todo.reserve(65536);
  for (int i = lo_; i < hi_; i++) { if (!admit[(size_t)i] || ((size_t)i < g_tdone.size() && g_tdone[(size_t)i])) continue; const int a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1]; if (a < 0 || b < 0) continue;
    if (!std::isfinite(SB.F[0].rfin[(size_t)a]) || !std::isfinite(SB.F[1].rfin[(size_t)b])) continue;
    if (postlin_side(SB, 0, a) && postlin_side(SB, 1, b)) { if (fin && g_tiny_chk && std::isfinite(err[i]) && err[i] < 1e20f && (i % 97) == 0) chk.push_back(i); continue; }
    todo.push_back(i); }
  std::atomic<long> nsc{0}, ndead{0}, nboth0{0}, nmix{0}; long nskip = 0;


  const size_t NU = SB.UQ.U.size(); std::vector<int> slot[2] = {std::vector<int>(NU, -1), std::vector<int>(NU, -1)}; std::vector<int> rep[2]; std::vector<int8_t> cls(todo.size(), 0);
  for (size_t q = 0; q < todo.size(); q++) { const int i = todo[q]; const int u[2] = {ct[(size_t)i * 2], ct[(size_t)i * 2 + 1]}; int nlin = 0; bool bad = false;
    for (int sd = 0; sd < 2; sd++) { const Tree& tr = SB.UQ.U[(size_t)u[sd]]; if (tr.Q > 16 || tr.nth > 32) { bad = true; break; }
      if (tr.nth >= 1) { if (!SB.F[sd].is_g[(size_t)u[sd]]) { bad = true; break; } nlin++; } }
    if (bad) { cls[q] = -1; continue; }
    if (g_tiny == 2 && nlin) { nskip++; continue; }
    cls[q] = nlin ? 2 : 3;
    for (int sd = 0; sd < 2; sd++) if (slot[sd][(size_t)u[sd]] < 0) { slot[sd][(size_t)u[sd]] = (int)rep[sd].size(); rep[sd].push_back(i); } }
  const double tt1 = ms_since(tt0);
  std::vector<double> zc[2]; std::vector<uint8_t> zok[2];
  for (int sd = 0; sd < 2; sd++) { zc[sd].assign(rep[sd].size() * NPT, 0.0); zok[sd].assign(rep[sd].size(), 0); }
  const size_t nz0 = rep[0].size(), nzt = nz0 + rep[1].size(); const size_t BZ = 64, NBZ = (nzt + BZ - 1) / BZ;
  PL(E).run((int)NBZ, [&](int c) { thread_local std::vector<int> K1(1);
    for (size_t q = (size_t)c * BZ; q < std::min(nzt, (size_t)c * BZ + BZ); q++) { const int sd = q < nz0 ? 0 : 1; const size_t l = q < nz0 ? q : q - nz0; const int i = rep[sd][l]; const int uid = ct[(size_t)i * 2 + sd];
      const float* th[2] = {nullptr, nullptr}; if (SB.UQ.U[(size_t)uid].nth >= 1) th[sd] = F14.th[sd].data() + 32 * (size_t)uid;
      K1[0] = i; double zf[2][NPT]; if (!sicis_res48_fast(E, SB.UQ, ct, K1, 0, nullptr, nil, zf, th, tgf, sd)) continue;
      std::memcpy(&zc[sd][l * NPT], zf[sd], sizeof(double) * NPT); zok[sd][l] = 1; } });
  const double tt2 = ms_since(tt0);
  const size_t NB = (todo.size() + 255) / 256;
  PL(E).run((int)NB, [&](int c) { const size_t lo = (size_t)c * 256, hi = std::min(todo.size(), lo + 256);
    for (size_t q = lo; q < hi; q++) { const int i = todo[q]; if (cls[q] == 0) continue;
      double sc = 0; bool ok = cls[q] > 0;
      if (ok) { const int l0 = slot[0][(size_t)ct[(size_t)i * 2]], l1 = slot[1][(size_t)ct[(size_t)i * 2 + 1]]; ok = zok[0][(size_t)l0] && zok[1][(size_t)l1];
        if (ok) { const double* z0 = &zc[0][(size_t)l0 * NPT]; const double* z1 = &zc[1][(size_t)l1 * NPT]; float res48[NPT];
          if (g_nres) { for (int p = 0; p < NPT; p++) res48[p] = (float)(z0[p] * (double)g_nres_n[p] + z1[p] * (double)g_nres_n[NPT + p]); }
          else for (int p = 0; p < NPT; p++) { const double h = std::hypot(z0[p], z1[p]); res48[p] = (float)((z0[p] + z1[p] < 0) ? -h : h); }
          sc = std::min(sicis::score48(res48, *g_sicis_en), 1.0e30); } }
      if (!ok || !std::isfinite(sc) || sc <= -1.0e30) { err[i] = DEAD; ndead.fetch_add(1, std::memory_order_relaxed); continue; }
      err[i] = (float)sc; nsc.fetch_add(1, std::memory_order_relaxed); (cls[q] == 2 ? nmix : nboth0).fetch_add(1, std::memory_order_relaxed); } });
  if (!fin) for (int i = lo_; i < hi_; i++) if ((size_t)i < g_tdone.size()) g_tdone[(size_t)i] = 1;
  if (getenv("TINY_T")) printf("TINY_T: classify %.2f ms, side evals %.2f (%zu+%zu uids), total %.2f ms (%zu children)\n", tt1, tt2 - tt1, rep[0].size(), rep[1].size(), ms_since(tt0), todo.size());
  if (fin && g_tiny_chk) { double bl = 1e30, bm = 1e30; int il = -1, im = -1; for (int i : todo) { if (!(err[i] < 1e20f)) continue; const int a = ct[(size_t)i * 2], b = ct[(size_t)i * 2 + 1]; const bool lin = SB.UQ.U[(size_t)a].nth == 0 && SB.UQ.U[(size_t)b].nth == 0;
      if (lin && err[i] < bl) { bl = err[i]; il = i; } if (!lin && err[i] < bm) { bm = err[i]; im = i; } }
    printf("TINY_FIX: best newly scored: both-linear nll %.6g (child %d, len %d/%d) | nonlinear-side nll %.6g (child %d)\n", bl, il, il >= 0 ? PS.cands[il]->len[0] : 0, il >= 0 ? PS.cands[il]->len[1] : 0, bm, im); }
  if (fin && g_tiny_chk) { std::vector<double> gap; for (int i : chk) { const int u[2] = {ct[(size_t)i * 2], ct[(size_t)i * 2 + 1]}; const float* th[2] = {nullptr, nullptr}; for (int sd = 0; sd < 2; sd++) if (SB.UQ.U[(size_t)u[sd]].nth >= 1) th[sd] = F14.th[sd].data() + 32 * (size_t)u[sd];
      std::vector<int> K1(1, i); double zf[2][NPT]; float r48[NPT]; if (!sicis_res48_fast(E, SB.UQ, ct, K1, 0, nullptr, nil, zf, th, tgf)) continue;
      if (g_nres) { for (int p = 0; p < NPT; p++) r48[p] = (float)(zf[0][p] * (double)g_nres_n[p] + zf[1][p] * (double)g_nres_n[NPT + p]); } else for (int p = 0; p < NPT; p++) { const double h = std::hypot(zf[0][p], zf[1][p]); r48[p] = (float)((zf[0][p] + zf[1][p] < 0) ? -h : h); }
      const double sc = std::min(sicis::score48(r48, *g_sicis_en), 1.0e30); if (std::isfinite(sc)) gap.push_back(sc - (double)err[i]); }
    std::sort(gap.begin(), gap.end()); const size_t m = gap.size();
    printf("TINY_FIX: scored %ld (nonlinear side %ld, both linear %ld), DEAD %ld, skipped %ld of %zu not-in-K children | K sample %zu: grid-t stage-8 score minus post16 err: min %.4g median %.4g p90 %.4g max %.4g\n", nsc.load(), nmix.load(), nboth0.load(), ndead.load(), nskip, todo.size(), m,
           m ? gap[0] : 0.0, m ? gap[m / 2] : 0.0, m ? gap[m * 9 / 10] : 0.0, m ? gap[m - 1] : 0.0); }
  else { static int ncall = 0; if (ncall++ < 6 || getenv("TINY_V")) printf("TINY_FIX: scored %ld (nonlinear side %ld, both linear %ld), DEAD %ld, skipped %ld of %zu not-in-K children\n", nsc.load(), nmix.load(), nboth0.load(), ndead.load(), nskip, todo.size()); }
}


static const int g_recfix_g = getenv("REC_FIX") ? atoi(getenv("REC_FIX")) : 1;
static inline bool rec_pinpost() { return g_recfix_g && !g_freet && g_sf_axis && g_sf_pin && !g_pinfix; }
static inline void rec_set_pin() { fitx::g_lazy_retain = g_recfix_g && !g_freet && g_sf_axis && (getenv("REC_LAZY") ? atoi(getenv("REC_LAZY")) != 0 : true); fitx::g_pin_t = rec_pinpost() ? (g_sf_axis == 1 ? 0 : 1) : -1; fitx::g_pin_op = evo::TVAR; }
static inline bool rec_isdrv(const evo::Cand* c, int s, int drv) { return s == drv && (rec_pinpost() || (c->len[s] == 1 && c->op[s][0] == evo::TVAR)); }
static inline int s_drv_q(int drv, int s, const evo::Cand* c, const Tree& t) { return (s == drv && c->len[s] == 1 && c->op[s][0] == evo::TVAR) ? 1 : t.Q + t.nth; }
struct FitRefT { int Q, nth, nS; const uint32_t* prog; size_t nprog; const uint16_t* sp; size_t nsp; const uint16_t* segs; };


static float k_batch_gate(Engine& E, size_t n, const float* err, const std::vector<uint8_t>& admit) {
  const float gg = fitx::g_gate; const int tn = fitx::topn(0); if (tn <= 0) return gg; const size_t N = (size_t)tn; const size_t B = 4096, nb = (n + B - 1) / B;
  std::vector<size_t> cnt(nb, 0); std::vector<std::vector<float>> keep(nb);
  PL(E).run((int)nb, [&](int b) { const size_t lo = (size_t)b * B, hi = std::min(n, lo + B); size_t c = 0; std::vector<float>& kv = keep[(size_t)b];
    for (size_t k = lo; k < hi; k++) { if (!admit[k] || fitx::was_captured((int)k)) continue; const float v = err[k]; if (!std::isfinite(v)) continue; c++; if (v <= gg) kv.push_back(v); }
    cnt[(size_t)b] = c; });
  size_t tot = 0, nk = 0; for (size_t b = 0; b < nb; b++) { tot += cnt[b]; nk += keep[b].size(); }
  float g = gg;
  if (tot > N && nk >= N) { std::vector<float> e; e.reserve(nk); for (size_t b = 0; b < nb; b++) e.insert(e.end(), keep[b].begin(), keep[b].end());
    std::nth_element(e.begin(), e.begin() + (long)(N - 1), e.end()); g = std::min(g, e[N - 1]); }
  return g;
}
static void fit_capture_fixed(Engine& E, const PipeStep& PS, StepBatch& SB, const float* err, int n, const std::vector<uint8_t>& admit) {

  if ((int)SB.ct.size() < 2 * n || (!g_freet && !g_sf_axis)) return;
  const int drv = g_sf_axis == 1 ? 0 : g_sf_axis == 2 ? 1 : -1; double tg[48]; for (int i = 0; i < 48; i++) tg[i] = (double)i / 47.0;
  std::vector<uint64_t> hc[2] = {std::vector<uint64_t>(SB.UQ.U.size(), 0), std::vector<uint64_t>(SB.UQ.U.size(), 0)};
  std::vector<std::shared_ptr<const fitx::Side>> pcache[2] = {std::vector<std::shared_ptr<const fitx::Side>>(SB.UQ.U.size()), std::vector<std::shared_ptr<const fitx::Side>>(SB.UQ.U.size())};
  fitx::Rec r; static long s_unc = 0; long unc = 0, u_neg = 0, u_nog = 0, u_nog_fin = 0, u_cap = 0;
  auto rpth = [&](int sd, int i, int uid) -> const float* { if (g_rec_pol && (size_t)i < g_rp_ok.size() && ((g_rp_ok[(size_t)i] >> sd) & 1)) { g_rp_n.fetch_add(1, std::memory_order_relaxed); return &g_rp_th[sd][(size_t)i * 32]; } return &F14.th[sd][32 * (size_t)uid]; };
  for (int s = 0; s < 2; s++) std::memcpy(fitx::g_target[s], E.tg[s].yf, 48 * sizeof(float));
  const auto jbg0 = clk::now();

  const bool rfx = g_recfix_g && !g_freet;
  rec_set_pin();
  const bool g_tiny_cap = g_tiny && g_sf_axis && !rfx;
  float bg = rfx ? INFINITY : (g_kcap && !g_freet && !g_tiny_cap) ? k_batch_gate(E, (size_t)n, err, admit) : fitx::batch_gate((size_t)n, [&](size_t k) { return (admit[k] && !fitx::was_captured((int)k)) ? err[k] : INFINITY; });
  if (g_jprof) printf("J_PROF batch_gate %.2f ms (bg %g gate %g)\n", ms_since(jbg0), (double)bg, (double)fitx::g_gate);


  std::vector<float> tkey; float tbg = INFINITY; const float gpop = fitx::g_gate;
  if (g_tiny_cap) { static float PT[40][26][26]; static bool pti = false; if (!pti) { for (auto& a : PT) for (auto& b : a) for (float& c : b) c = NAN; pti = true; }
    tkey.assign((size_t)n, INFINITY); std::vector<int> kk((size_t)n, -1);
    for (int i = 0; i < n; i++) { if (!admit[(size_t)i] || !isfinite(err[i]) || fitx::was_captured(i)) continue; const int u0 = SB.ct[(size_t)i * 2], u1 = SB.ct[(size_t)i * 2 + 1]; if (u0 < 0 || u1 < 0) continue;
      const int k = s_drv_q(drv, 0, PS.cands[i], SB.UQ.U[(size_t)u0]) + s_drv_q(drv, 1, PS.cands[i], SB.UQ.U[(size_t)u1]); const int nx = PS.cands[i]->len[0], ny = PS.cands[i]->len[1]; kk[(size_t)i] = k;
      if (k >= 0 && k < 40 && nx <= 25 && ny <= 25) { float& pt = PT[k][nx][ny]; if (std::isnan(pt)) pt = (float)fp::fallback_terms(k, 48, nx, ny).total; } }
    PL(E).run((int)((n + 4095) / 4096), [&](int b) { for (int i = b * 4096; i < std::min(n, b * 4096 + 4096); i++) { if (kk[(size_t)i] < 0 && !(admit[(size_t)i] && isfinite(err[i]))) continue; if (!admit[(size_t)i] || !isfinite(err[i]) || fitx::was_captured(i)) continue;
      const int k = kk[(size_t)i]; if (k < 0) continue; const int nx = PS.cands[i]->len[0], ny = PS.cands[i]->len[1];
      tkey[(size_t)i] = (k < 40 && nx <= 25 && ny <= 25) ? err[i] + PT[k][nx][ny] : err[i]; } });
    std::vector<float> e; e.reserve(tkey.size()); for (float v : tkey) if (std::isfinite(v)) e.push_back(v);
    const size_t N = (size_t)std::max(1, fitx::topn(0)); if (e.size() > N) { std::nth_element(e.begin(), e.begin() + (long)(N - 1), e.end()); tbg = e[N - 1]; }
    bg = INFINITY; }
  static const int g_rdiag = getenv("REC_DIAG") ? atoi(getenv("REC_DIAG")) : 0;
  if (g_rdiag) { long na = 0, nf = 0, ng = 0, nwc = 0, nbad = 0; for (int i = 0; i < n; i++) { if (!admit[(size_t)i]) continue; na++; if (!isfinite(err[i])) { nbad++; continue; } nf++; if (fitx::was_captured(i)) nwc++; if (err[i] <= bg) ng++; }
    printf("REC_DIAG capture: n %d admitted %ld finite %ld nonfinite %ld wascap %ld err<=bg %ld (bg %g gate %g)\n", n, na, nf, nbad, nwc, ng, (double)bg, (double)fitx::g_gate); }
  if (rfx && !g_fit_chk) {
    std::vector<int> el; el.reserve((size_t)n); std::vector<uint8_t> need[2] = {std::vector<uint8_t>(SB.UQ.U.size(), 0), std::vector<uint8_t>(SB.UQ.U.size(), 0)}; std::vector<int> nl[2];
    {
      const int nb = (n + 4095) / 4096; struct CB { std::vector<int> el; long neg = 0, nog = 0; }; std::vector<CB> cb((size_t)nb);
      PL(E).run(nb, [&](int b) { CB& B = cb[(size_t)b]; for (int i = b * 4096; i < std::min(n, b * 4096 + 4096); i++) {
        if (!admit[(size_t)i] || !isfinite(err[i]) || fitx::was_captured(i)) continue;
        const int uid[2] = {SB.ct[(size_t)i * 2], SB.ct[(size_t)i * 2 + 1]}; if (uid[0] < 0 || uid[1] < 0) { B.neg++; continue; }
        bool drvs[2] = {rec_isdrv(PS.cands[i], 0, drv), rec_isdrv(PS.cands[i], 1, drv)}; bool ok = true;
        for (int sd = 0; sd < 2; sd++) { if (drvs[sd]) continue; if (SB.UQ.U[(size_t)uid[sd]].nth >= 1 && !SB.F[sd].is_g[(size_t)uid[sd]]) { ok = false; break; } }
        if (!ok) { B.nog++; continue; }
        for (int sd = 0; sd < 2; sd++) if (!drvs[sd]) need[sd][(size_t)uid[sd]] = 1;
        B.el.push_back(i); } });
      for (const CB& B : cb) { unc += B.neg + B.nog; u_neg += B.neg; u_nog += B.nog; el.insert(el.end(), B.el.begin(), B.el.end()); }
      for (int sd = 0; sd < 2; sd++) for (size_t u = 0; u < need[sd].size(); u++) if (need[sd][u]) nl[sd].push_back((int)u);
    }
    {
      auto fref = [&](int sd, int u) { const Tree& t = SB.UQ.U[(size_t)u]; return FitRefT{t.Q, t.nth, t.nS, h313a::prog(t), h313a::nprog(t), t.sp.data(), t.sp.size(), h313a::segs(t)}; };
      std::vector<uint64_t> hh[2]; for (int sd = 0; sd < 2; sd++) hh[sd].resize(nl[sd].size());
      const size_t n0 = nl[0].size(), nt = n0 + nl[1].size();
      PL(E).run((int)((nt + 1023) / 1024), [&](int b) { for (size_t q = (size_t)b * 1024; q < std::min(nt, (size_t)b * 1024 + 1024); q++) { const int sd = q < n0 ? 0 : 1; const size_t l = q < n0 ? q : q - n0; const FitRefT f = fref(sd, nl[sd][l]); hh[sd][l] = fitx::shash(f) | 1; } });
      fitx::Store& S = fitx::st(); std::lock_guard<std::mutex> lk(S.mx);
      struct Miss { uint64_t h; int sd, u; std::shared_ptr<const fitx::Side> r; }; std::vector<Miss> ms; std::unordered_map<uint64_t, size_t> mi;
      for (int sd = 0; sd < 2; sd++) for (size_t l = 0; l < nl[sd].size(); l++) { const uint64_t h = hh[sd][l]; const int u = nl[sd][l];
        auto ii = S.intern.find(h); std::shared_ptr<const fitx::Side> r; if (ii != S.intern.end()) r = ii->second.lock();
        if (r) { pcache[sd][(size_t)u] = r; continue; }
        auto m = mi.find(h); if (m == mi.end()) { mi.emplace(h, ms.size()); ms.push_back({h, sd, u, nullptr}); } }
      PL(E).run((int)((ms.size() + 255) / 256), [&](int b) { for (size_t q = (size_t)b * 256; q < std::min(ms.size(), (size_t)b * 256 + 256); q++) { float z[32] = {0}; auto sdp = std::make_shared<fitx::Side>(); const FitRefT f = fref(ms[q].sd, ms[q].u); if (fitx::build_side(f, z, z, *sdp)) ms[q].r = sdp; } });
      for (Miss& m : ms) { if (m.r) S.intern[m.h] = m.r; else S.nbuildfail++; }
      for (int sd = 0; sd < 2; sd++) for (size_t l = 0; l < nl[sd].size(); l++) { const int u = nl[sd][l]; if (pcache[sd][(size_t)u]) continue; auto m = mi.find(hh[sd][l]); if (m != mi.end()) pcache[sd][(size_t)u] = ms[m->second].r; }
    }
    { std::vector<int> el2; el2.reserve(el.size()); for (int i : el) { const int uid[2] = {SB.ct[(size_t)i * 2], SB.ct[(size_t)i * 2 + 1]}; bool good = true;
        for (int sd = 0; sd < 2; sd++) { if (rec_isdrv(PS.cands[i], sd, drv)) continue; if (!pcache[sd][(size_t)uid[sd]]) good = false; }
        if (good) el2.push_back(i); else { unc++; u_cap++; } } el.swap(el2); }
    const double rt1 = ms_since(jbg0);
    struct Blk { std::vector<fitx::Cap> c; std::vector<float> ar; std::vector<int> ch; }; const size_t NE = el.size(); const int nbk = (int)((NE + 4095) / 4096); std::vector<Blk> bl((size_t)nbk);
    static const std::shared_ptr<const fitx::Side> dp = fitx::drv_proto();
    PL(E).run(nbk, [&](int b) { Blk& B = bl[(size_t)b]; const size_t lo = (size_t)b * 4096, hi = std::min(NE, lo + 4096); B.c.reserve(hi - lo); B.ar.reserve((hi - lo) * 24); B.ch.reserve(hi - lo);
      for (size_t q = lo; q < hi; q++) { const int i = el[q]; const int uid[2] = {SB.ct[(size_t)i * 2], SB.ct[(size_t)i * 2 + 1]}; fitx::Cap c; c.key = fitx::gkey(*PS.cands[i]); c.err = err[i]; c.kind = 2;
        const float* th[2] = {nullptr, nullptr}; const float* a[2] = {nullptr, nullptr};
        for (int sd = 0; sd < 2; sd++) { const bool isd = rec_isdrv(PS.cands[i], sd, drv); c.p[sd] = isd ? dp : pcache[sd][(size_t)uid[sd]]; if (!isd) { th[sd] = &F14.th[sd][32 * (size_t)uid[sd]]; a[sd] = &F14.a[sd][16 * (size_t)uid[sd]]; if (g_rec_pol && (size_t)i < g_rp_ok.size() && ((g_rp_ok[(size_t)i] >> sd) & 1)) { th[sd] = &g_rp_th[sd][(size_t)i * 32]; g_rp_n.fetch_add(1, std::memory_order_relaxed); } } }
        c.off = B.ar.size();
        for (int sd = 0; sd < 2; sd++) if (!c.p[sd]->drv) { B.ar.insert(B.ar.end(), th[sd], th[sd] + c.p[sd]->nth); B.ar.insert(B.ar.end(), a[sd], a[sd] + c.p[sd]->Q); }
        B.c.push_back(std::move(c)); B.ch.push_back(i); } });
    { fitx::Store& S = fitx::st(); std::lock_guard<std::mutex> lk(S.mx); size_t tc = 0, ta = 0; for (Blk& B : bl) { tc += B.c.size(); ta += B.ar.size(); } S.cur.reserve(S.cur.size() + tc); S.ar.reserve(S.ar.size() + ta);
      for (Blk& B : bl) { const size_t base = S.ar.size(); S.ar.insert(S.ar.end(), B.ar.begin(), B.ar.end());
        for (size_t q = 0; q < B.c.size(); q++) { B.c[q].off += base; S.cur.push_back(std::move(B.c[q])); S.ncap++; const int ch = B.ch[q]; if ((size_t)ch < S.got.size()) S.got[(size_t)ch] = 1; } } }
    if (g_rdiag) printf("REC_DIAG captured %zu of %d (unexportable %ld) | prepass %.1f ms, pool+merge %.1f ms\n", NE, n, unc, rt1, ms_since(jbg0) - rt1);
    if (unc) s_unc += unc;
    return;
  }
  for (int i = 0; i < n; i++) {
    if (!admit[(size_t)i] || !isfinite(err[i]) || fitx::was_captured(i)) continue;
    if (g_tiny_cap ? (err[i] > gpop || tkey[(size_t)i] > tbg) : err[i] > bg) continue;
    const int uid[2] = {SB.ct[(size_t)i * 2], SB.ct[(size_t)i * 2 + 1]}; if (uid[0] < 0 || uid[1] < 0) { unc++; u_neg++; continue; }
    const Tree* T[2] = {&SB.UQ.U[(size_t)uid[0]], &SB.UQ.U[(size_t)uid[1]]}; FitRefT f[2]; bool isd[2] = {false, false}; bool ok = true;
    for (int s = 0; s < 2 && ok; s++) {
      if (rec_isdrv(PS.cands[i], s, drv)) { isd[s] = true; continue; }
      if (T[s]->nth >= 1 && !SB.F[s].is_g[(size_t)uid[s]]) ok = false;
      f[s] = FitRefT{T[s]->Q, T[s]->nth, T[s]->nS, h313a::prog(*T[s]), h313a::nprog(*T[s]), T[s]->sp.data(), T[s]->sp.size(), h313a::segs(*T[s])};
    }
    if (!ok) { unc++; u_nog++; if (g_fit_chk && isfinite(err[i]) && err[i] < 1e20) u_nog_fin++; continue; }
    auto hsh = [&](int sd) -> uint64_t { if (isd[sd]) return 0; uint64_t& h = hc[sd][(size_t)uid[sd]]; if (!h) h = fitx::shash(f[sd]) | 1; return h; };
    if (rfx) { const std::shared_ptr<const fitx::Side>* pp[2]; for (int sd = 0; sd < 2; sd++) { if (isd[sd]) { static const std::shared_ptr<const fitx::Side> dp = fitx::drv_proto(); pp[sd] = &dp; continue; } auto& slot = pcache[sd][(size_t)uid[sd]]; if (!slot) slot = fitx::proto_of(f[sd], hsh(sd)); pp[sd] = &slot; }
      if (!fitx::capture_p(i, fitx::gkey(*PS.cands[i]), err[i], 2, *pp[0], *pp[1], isd[0] ? nullptr : rpth(0, i, uid[0]), isd[0] ? nullptr : &F14.a[0][16 * (size_t)uid[0]], isd[1] ? nullptr : rpth(1, i, uid[1]), isd[1] ? nullptr : &F14.a[1][16 * (size_t)uid[1]], nullptr)) { unc++; u_cap++; continue; } }
    else if (!fitx::capture(i, fitx::gkey(*PS.cands[i]), err[i], 2, isd[0] ? (const FitRefT*)nullptr : &f[0], hsh(0), isd[0] ? nullptr : rpth(0, i, uid[0]), isd[0] ? nullptr : &F14.a[0][16 * (size_t)uid[0]],
                       isd[1] ? (const FitRefT*)nullptr : &f[1], hsh(1), isd[1] ? nullptr : rpth(1, i, uid[1]), isd[1] ? nullptr : &F14.a[1][16 * (size_t)uid[1]], nullptr)) { unc++; u_cap++; continue; }
    if (g_fit_chk) { fitx::materialize(fitx::st().cur.back(), r); for (int s = 0; s < 2; s++) { if (r.s[s].drv) continue; double yy = 0;
      { fitx::Side L = r.s[s]; fitx::ls_fit(L, tg, E.tg[s].yf); static long ln = 0, lbad = 0; static double lmax = 0; const double R2 = fit_recon_rss(L, tg, E.tg[s].yf), F = SB.F[s].rfin[(size_t)uid[s]]; double y2 = 0; for (int k = 0; k < 48; k++) y2 += (double)E.tg[s].yf[k] * E.tg[s].yf[k];
        const double rel = std::fabs(F - R2) / std::max(std::max(F, R2), 1e-5 * y2); ln++; if (rel > 1e-3) lbad++; lmax = std::max(lmax, rel); if (ln % 20000 == 0) printf("FIT_CHECK fixed-t LS-refit a (F14.th, a recomputed): %ld sides, rel>1e-3 %ld, max rel %.3g\n", ln, lbad, lmax); } for (int k = 0; k < 48; k++) yy += (double)E.tg[s].yf[k] * E.tg[s].yf[k];
      fit_chk_one((double)SB.F[s].rfin[(size_t)uid[s]], fit_recon_rss(r.s[s], tg, E.tg[s].yf), 1e-5 * yy); } }
  }
  if (g_rdiag) printf("REC_DIAG proto total %.1f ms of capture %.1f ms\n", fitx::g_ms_proto, ms_since(jbg0));
  if (unc) { s_unc += unc; if (g_fit_chk) printf("FIT_EXPORT: %ld admitted children without an exportable fit this call (total %ld): no tree %ld, nth>=1 without GN theta %ld, build failed %ld\n", unc, s_unc, u_neg, u_nog, u_cap); }
}


static void fit_set_theta(fitx::Side& S, const double* th) { for (int i = 0; i < S.nth; i++) S.th[i] = (float)th[i]; for (auto& tr : S.thref) S.n[tr.first].c = S.th[tr.second]; }
static double fit_rss_ls(fitx::Side& S, const double* th, const double* tg, const float* Y) {
  fit_set_theta(S, th); fitx::ls_fit(S, tg, Y); std::vector<double> v(S.n.size() + 1); double r2 = 0;
  for (int i = 0; i < 48; i++) { const double r = (double)Y[i] - fitx::eval_side(S, tg[i], v.data()); r2 += r * r; }
  return std::isfinite(r2) ? r2 : 1e300;
}
static bool fit_refit_side(const evo::Cand& g, int s, const float* Y, fitx::Side& out, double& rss) {
  static skn::Ctx* C = new skn::Ctx(); static skn::Prog* P = new skn::Prog(); skn::NFOut nf; double tg[48]; for (int i = 0; i < 48; i++) tg[i] = i / 47.0;
  if (!skn::nf_tree(*C, g.op[s], g.c[s], g.len[s], &nf) || !skn::compile(*C, nf, g.c[s], *P)) return false;
  skn::View V(P->w); ::Tree T; if (!compile_tree(V, T)) return false;
  FitRefT f{T.Q, T.nth, T.nS, h313a::prog(T), h313a::nprog(T), T.sp.data(), T.sp.size(), h313a::segs(T)};
  float z[32] = {0}; if (!fitx::build_side(f, z, z, out)) return false;
  const int n = out.nth; std::vector<double> th((size_t)n); for (int i = 0; i < n; i++) th[(size_t)i] = P->thd[i];
  double F = fit_rss_ls(out, th.data(), tg, Y), lam = 1e-3;


  { static const int sc = getenv("FP_REFIT_SCAN") ? atoi(getenv("FP_REFIT_SCAN")) : 1;
    for (int pass = 0; sc && pass < 2 && n > 0; pass++) for (int j = 0; j < n; j++) { const double keep = th[(size_t)j]; double bv = keep;
      for (int k = -2; k <= 128; k++) for (int sg = -1; sg <= 1; sg += 2) { if (k == 0) continue; th[(size_t)j] = sg * (k < 0 ? std::ldexp(1.0, k - 1) : 0.5 * k); const double Fc = fit_rss_ls(out, th.data(), tg, Y); if (Fc < F) { F = Fc; bv = th[(size_t)j]; } }
      th[(size_t)j] = bv; } }
  for (int it = 0; it < 80 && n > 0; it++) {
    std::vector<std::vector<double>> J((size_t)n, std::vector<double>(48)); std::vector<double> r0(48), r1(48), g((size_t)n, 0.0); std::vector<std::vector<double>> H((size_t)n, std::vector<double>((size_t)n, 0.0));
    auto resid = [&](const std::vector<double>& t_, std::vector<double>& r) { fitx::Side S2 = out; fit_rss_ls(S2, t_.data(), tg, Y); std::vector<double> v(S2.n.size() + 1); for (int i = 0; i < 48; i++) r[(size_t)i] = (double)Y[i] - fitx::eval_side(S2, tg[i], v.data()); };
    resid(th, r0);
    for (int j = 0; j < n; j++) { const double h = 1e-3 * (1.0 + std::fabs(th[(size_t)j])); std::vector<double> tp = th, tm = th; tp[(size_t)j] += h; tm[(size_t)j] -= h; resid(tp, r0.size() ? r1 : r1); std::vector<double> rm(48); resid(tm, rm); for (int i = 0; i < 48; i++) J[(size_t)j][(size_t)i] = (r1[(size_t)i] - rm[(size_t)i]) / (2 * h); }
    for (int a = 0; a < n; a++) { for (int i = 0; i < 48; i++) g[(size_t)a] += J[(size_t)a][(size_t)i] * r0[(size_t)i]; for (int b = 0; b < n; b++) for (int i = 0; i < 48; i++) H[(size_t)a][(size_t)b] += J[(size_t)a][(size_t)i] * J[(size_t)b][(size_t)i]; }
    bool acc = false;
    for (int tr = 0; tr < 12 && !acc; tr++) {
      std::vector<std::vector<double>> A = H; std::vector<double> d((size_t)n); for (int a = 0; a < n; a++) { A[(size_t)a][(size_t)a] += lam * (H[(size_t)a][(size_t)a] + 1e-12); d[(size_t)a] = -g[(size_t)a]; }
      bool okk = true; for (int c = 0; c < n && okk; c++) { int pv = c; for (int r = c + 1; r < n; r++) if (std::fabs(A[(size_t)r][(size_t)c]) > std::fabs(A[(size_t)pv][(size_t)c])) pv = r; if (std::fabs(A[(size_t)pv][(size_t)c]) < 1e-300) { okk = false; break; } std::swap(A[(size_t)c], A[(size_t)pv]); std::swap(d[(size_t)c], d[(size_t)pv]);
        for (int r = c + 1; r < n; r++) { const double m = A[(size_t)r][(size_t)c] / A[(size_t)c][(size_t)c]; for (int k = c; k < n; k++) A[(size_t)r][(size_t)k] -= m * A[(size_t)c][(size_t)k]; d[(size_t)r] -= m * d[(size_t)c]; } }
      if (okk) { for (int c = n - 1; c >= 0; c--) { for (int k = c + 1; k < n; k++) d[(size_t)c] -= A[(size_t)c][(size_t)k] * d[(size_t)k]; d[(size_t)c] /= A[(size_t)c][(size_t)c]; }
        std::vector<double> tn = th; for (int a = 0; a < n; a++) tn[(size_t)a] += d[(size_t)a]; fitx::Side S2 = out; const double Fn = fit_rss_ls(S2, tn.data(), tg, Y);
        if (Fn < F) { th = tn; F = Fn; lam = std::max(lam / 3, 1e-9); acc = true; } else lam *= 4; } else lam *= 10;
    }
    if (!acc) break;
  }
  F = fit_rss_ls(out, th.data(), tg, Y); rss = F; return true;
}

static bool fit_get(Engine& E, const evo::Cand& g, fitx::Rec& r, double* refit_rss = nullptr) {
  static const bool force = getenv("FIT_FORCE_REFIT") != nullptr;
  if (!force && fitx::find(g, r)) return true;
  if (!fitx::on()) return false;
  r = fitx::Rec(); r.key = fitx::gkey(g); r.kind = 3; double rs[2] = {0, 0};
  for (int s = 0; s < 2; s++) {
    if (g_sf_axis && s == (g_sf_axis == 1 ? 0 : 1) && g.len[s] == 1 && g.op[s][0] == evo::TVAR) { r.s[s].drv = true; r.s[s].Q = 1; continue; }
    if (!fit_refit_side(g, s, E.tg[s].yf, r.s[s], rs[s])) return false;
  }
  if (refit_rss) *refit_rss = rs[0] + rs[1];
  { static long nrf = 0; if (g_fit_chk && nrf++ < 20) printf("FIT_EXPORT: shown member without stored fit -> CPU refit rss_x+rss_y %.6g\n", rs[0] + rs[1]); }
  return true;
}


static void sp_fill_k(Engine& E, const PipeStep& PS, const StepBatch& SB, int n, const std::vector<uint8_t>& admit, std::vector<Meta>& M) {
  const int drv = g_sf_axis == 1 ? 0 : g_sf_axis == 2 ? 1 : -1;
  PL(E).run((n + 4095) / 4096, [&](int b) {
    long unk = 0;
    for (int i = b * 4096; i < std::min(n, b * 4096 + 4096); i++) {
      if (!admit[(size_t)i]) continue;
      const evo::Cand* c = PS.cands[i]; int k = 0; int nn[2] = {c->len[0], c->len[1]};
      for (int s = 0; s < 2; s++) { const int u = SB.ct[(size_t)i * 2 + s];
        if (s == drv && g_sf_pin) { k += 1; nn[s] = 1; continue; }
        if (u < 0) { k += fp::nml::const_leaves(c->op[s], c->len[s]); unk++; } else k += s_drv_q(drv, s, c, SB.UQ.U[(size_t)u]); }
      M[(size_t)i].kf = (uint8_t)std::min(k, 255); M[(size_t)i].nx = (uint8_t)nn[0]; M[(size_t)i].ny = (uint8_t)nn[1];
    }
    if (unk) g_sp_kunk.fetch_add(unk, std::memory_order_relaxed); });
}


static const int g_h418 = (getenv("H418") ? atoi(getenv("H418")) : 0) | (g_h874 & 1);
static const int g_h_sb_nofill = getenv("H_SB_NOFILL") ? atoi(getenv("H_SB_NOFILL")) : 0;
struct H418Store { FlatKeySet step, uidx; std::vector<FlatKeySet> stepSh, uidxSh; std::vector<std::vector<uint32_t>> shMap; std::vector<uint8_t> fitted[2]; };
static H418Store g_h418s;
static void h418_swap(PipeStep& P) {
  std::swap(P.step, g_h418s.step); std::swap(P.uidx, g_h418s.uidx); std::swap(P.stepSh, g_h418s.stepSh); std::swap(P.uidxSh, g_h418s.uidxSh);
  std::swap(P.shMap, g_h418s.shMap); std::swap(P.fitted[0], g_h418s.fitted[0]); std::swap(P.fitted[1], g_h418s.fitted[1]);
}
static void h418_take(PipeStep& P) { h418_swap(P); }
static void h418_put(PipeStep& P) { h418_swap(P); }
static void evaluate_v7(Engine& E, const evo::Cand* const* cands, const uint8_t* valid, int n, const FlatKeySet& popset,
                        float* err, float* rx, float* ry, std::vector<uint8_t>& admit, std::vector<Meta>& M, Times& T, StepBatch& SB) {
  auto tw = clk::now();
  if (fitx::on()) fitx::begin(n);
  rp_begin(n);
  CpuSpan265 bal_span; g_fit_span = &bal_span;
  if (g_h771) h771_begin(n);
  const int NCH = (n + 63) / 64;
  PipeStep PS; PS.n = n; PS.NCH = NCH; PS.cands = cands; PS.valid = valid; PS.popset = &popset;
  PS.admit = &admit; PS.M = &M; PS.SB = &SB;
  if (g_h264) { if ((int)M.size() != n) M.resize(n); } else M.assign(n, Meta());
  admit.assign(n, 0);
  if (g_h264) { if ((int)NFD.ch.size() != n) NFD.ch.resize(n); } else NFD.ch.assign(n, nf16::Child());
  NFD.sc.resize(n);
  if ((int)NFD.blk.size() < NCH) NFD.blk.resize(NCH);
  SB.ct.assign((size_t)n * 2, -1);
  SB.UQ.U.clear(); SB.UQ.hash.clear(); SB.UQ.enc.clear();
  if (g_h864 & 4) { SB.UQ.U.reserve((size_t)2 * n + 2); SB.UQ.hash.reserve((size_t)2 * n + 2); SB.UQ.enc.reserve((size_t)2 * n + 2); }
  for (int s = 0; s < 2; s++) { SB.F[s].init(2 * n + 2); SB.ids[s].clear(); if (!(g_h418 && g_h_sb_nofill)) PS.fitted[s].assign((size_t)2 * n + 2, 0); }
  F14.init((size_t)2 * n + 2);
  if (g_h870 & 2) { PS.rep.reserve((size_t)2 * n + 2); PS.fid.reserve((size_t)2 * n + 2); }
  struct H418Guard { PipeStep& P; bool on; ~H418Guard() { if (on) h418_put(P); } } h418g{PS, g_h418 != 0};
  if (g_h418) {
    h418_take(PS);
    if (!(g_h92 & 2)) { PS.step.clear(); PS.uidx.clear(); PS.step.reserve((size_t)n); PS.uidx.reserve((size_t)n * 2); }
    for (int s = 0; s < 2; s++) PS.fitted[s].assign((size_t)2 * n + 2, 0);
  } else
  PS.step.reserve((size_t)n), PS.uidx.reserve((size_t)n * 2);
  if (g_h92 & 2) {
    const int NSh = g_h410ns > 0 ? std::min(254, g_h410ns) : std::max(1, (int)WK(E).size());
    if (g_h418 && (int)PS.stepSh.size() == NSh) {
      for (int s = 0; s < NSh; s++) { PS.stepSh[(size_t)s].clear(); PS.uidxSh[(size_t)s].clear(); PS.shMap[(size_t)s].clear(); }
    } else {
    PS.stepSh.assign((size_t)NSh, FlatKeySet()); PS.uidxSh.assign((size_t)NSh, FlatKeySet()); PS.shMap.assign((size_t)NSh, std::vector<uint32_t>());
    }
    for (int s = 0; s < NSh; s++) { PS.stepSh[(size_t)s].reserve((size_t)n / (size_t)NSh + 1); PS.uidxSh[(size_t)s].reserve((size_t)n * 2 / (size_t)NSh + 1); }
  }

  if (g_h890 & 1) {
    PS.step.byref = PS.uidx.byref = true;
    for (auto& x : PS.stepSh) x.byref = true;
    for (auto& x : PS.uidxSh) x.byref = true;
  }

  const int K = std::max(1, std::min(g_chunks, NCH));
  const int CS = ((NCH + K - 1) / K) * 64;


  std::vector<int> cstart;
  if (g_h402 && K > 2 && NCH >= 2 * K) {
    const double rho = g_h402r > 0 ? g_h402r : g_h402_rho;
    std::vector<double> w((size_t)K); double W = 0;
    for (int q = 0; q < K; q++) { w[(size_t)q] = std::pow(rho, (double)std::min(q, K - 1 - q)); W += w[(size_t)q]; }
    double pre = 0; int last = 0; cstart.push_back(0);
    for (int q = 0; q + 1 < K; q++) { pre += w[(size_t)q]; int b = (int)std::llround((double)NCH * pre / W); if (b <= last) b = last + 1; if (b >= NCH) break; cstart.push_back(b * 64); last = b; }
    if (g_h402v) { printf("H402: rho %.3f K %d blocks:", rho, (int)cstart.size()); for (size_t q = 0; q < cstart.size(); q++) printf(" %d", ((q + 1 < cstart.size() ? cstart[q + 1] : NCH * 64) - cstart[q]) / 64); printf("\n"); }
  } else if (g_h301t && (g_h257 & 2) && K > 2 && NCH >= 2 * K) {
    const int m0 = std::max(1, NCH / 32), mid = NCH - 2 * m0, per = (mid + (K - 2) - 1) / (K - 2);
    cstart.push_back(0);
    for (int b = m0; b < NCH - m0; b += per) cstart.push_back(b * 64);
    cstart.push_back((NCH - m0) * 64);
  } else if ((g_h257 & 2) && K > 1 && NCH >= 2 * K) {
    const int m0 = std::max(1, NCH / 32), per = (NCH - m0 + (K - 1) - 1) / (K - 1);
    cstart.push_back(0);
    for (int b = m0; b < NCH; b += per) cstart.push_back(b * 64);
  } else {
    for (int c = 0; c * CS < n; c++) cstart.push_back(c * CS);
  }
  const int Kn = (int)cstart.size();
  if (g_h774) g_h774_cstart = cstart;
  std::vector<Batch> BT(Kn);
  const bool tiny_early = g_tiny && g_tiny_early && g_sf_axis && !g_freet && !g_h217 && g_sicis_en && g_sicis_en->present && Kn >= 2;
  std::deque<std::pair<int, int>> tq; if (g_tiny) g_tdone.assign((size_t)n, 0);
  if (tiny_early) g_tfill = [&]() -> bool { if (tq.empty()) return false; static const int sl = getenv("TINY_SLICE") ? atoi(getenv("TINY_SLICE")) : 2048; std::pair<int, int> r = tq.front(); if (r.second - r.first > sl) { tq.front().first = r.first + sl; r.second = r.first + sl; } else tq.pop_front(); tiny_score(E, PS, SB, err, n, admit, r.first, r.second); return !tq.empty(); };
  for (int c = 0; c < Kn; c++) { BT[c].b0 = cstart[c]; BT[c].b1 = c + 1 < Kn ? cstart[c + 1] : n; BT[c].ci = c; }
  h301_mark("K: ev chunks setup");
  if (g_h826_head) { h301_mark("K: breed head hook"); auto hk = std::move(g_h826_head); g_h826_head = nullptr; hk(BT[0].b1); }


  const bool h897on = (g_h897 & 1) != 0 && !g_h217 && !g_h250b && !g_h650d;
  int h897_n = 0; double h897_ms = 0;
  auto h897_early = [&](int c, Times& T) {
    if (c < 0) return;
    Batch& Bc = BT[c];
    if (Bc.egn || Bc.sdone || Bc.gbegun) return;
    if (c >= 2 && !BT[c - 2].gdone) return;
    if (!g_h882q->is_done(Bc.t882S)) return;
    if (Bc.cbS && Bc.cbS.status < MTLCommandBufferStatusCompleted) return;
    const auto e0 = clk::now();
    h301_mark("early scan_finish");
    scan_finish(E, PS, Bc, T);
    h301_mark("early gn_begin encode");
    if (g_h93) { gn_begin_p(E, PS, Bc, T);  } else gn_begin(E, PS, Bc, T);
    Bc.egn = true; h897_n++;
    h897_ms += ms_since(e0);
    h301_mark("front(c+1)");
  };

  auto front = [&](int c, Times& T) {
    Batch& B = BT[c]; B.P = &PP[g_h217 ? c % 3 : (c & 1)];
    if (c >= 2 && !BT[c - 2].sdone) { fprintf(stderr, "pipeline: scan slot %d reused before S(%d) finished\n", c & 1, c - 2); exit(1); }
    h401_sub("fe: NF"); if (g_h264) nf_range_h264(E, PS, B.b0, B.b1, T); else nf_range(E, PS, B.b0, B.b1, T);
    if (h897on) h897_early(c - 1, T);
    h401_sub("fe: dedup+unique"); if (g_h92 & 2) { dedup_h92(E, PS, B.b0, B.b1, T); if (g_h257 & 4) build_unique_h257(E, PS, B.b0, B.b1, T); else build_unique_h92(E, PS, B.b0, B.b1, T); }
    else { dedup_range(PS, B.b0, B.b1, T); build_unique_range(E, PS, B.b0, B.b1, T); }
    if (h897on) h897_early(c - 1, T);
    h401_sub("fe: schedule"); schedule_range(E, PS, B, T);
    h401_sub("fe: scan_begin"); if (g_h93) { scan_begin_p(E, PS, B, T);  } else scan_begin(E, PS, B, T);
    h401_sub("fe: end");
  };


  auto post_chunk = [&](const Batch& B) {
    if (g_freet) { freet_chunk(E, PS, B, err, rx, ry, T); return; }
    if (g_stage8_only || g_post16 != 1) return;
    std::vector<int> cc; cc.reserve((size_t)(B.b1 - B.b0));
    for (int i = B.b0; i < B.b1; i++) cc.push_back(i);
    if (g_h53) post16_h53(E, PS, SB, cc, err, rx, ry, T); else post16(E, PS, SB, cc, err, rx, ry, T);
  };


  const bool h84c = g_h84 >= 2 && g_h53 && g_post16 == 2;
  const bool h91c = h84c && g_h91 && !g_chk_band && !g_chk_fin && !g_time14;
  const bool h117c = h84c && !h91c && !g_h91 && (g_h117 & 3) && !g_chk_band && !g_chk_fin && !g_time14;
  const int n91 = g_h414r > 0 ? std::min(g_h414r, Kn) : Kn;
  if (h91c && (int)g_h91s.size() < n91) g_h91s.resize((size_t)n91);
  auto pre_chunk = [&](const Batch& B) {
    if (g_freet) return;
    if (!h84c || g_stage8_only) return;
    std::vector<int> cc; cc.reserve((size_t)(B.b1 - B.b0));
    for (int i = B.b0; i < B.b1; i++) cc.push_back(i);
    if (h91c) { const int ci = (int)(&B - BT.data()); H91Slot& SL = g_h91s[(size_t)(ci % n91)];
      if (SL.pend) { h301_mark("H414 slot reuse: collect"); post16_h91_col(E, SL, SB, err, rx, ry, T); }
      g_h53p = &SL.hs; post16_h53_prep(E, PS, SB, cc, T); post16_h91_enc(E, PS, SB, SL, T); g_h53p = &g_h53s; return; }
    post16_h53_prep(E, PS, SB, cc, T);
    if (h117c && (g_h117 & 1)) h117_cpu(E, PS, SB, g_h53s, T);
    if (g_h731 && h117c && (g_h117 & 1) && (int)(&B - BT.data()) < Kn - 1) h731_chunk(E, PS, g_h53s, T);
    if (g_h680) h680_inject(E.dev, Kn, (int)(&B - BT.data()) == Kn - 1);
  };
  if (g_h891 && g_h891_head_us > 0) g_h891_head_ms += h891_main_spin_ms(g_h891_head_us);
  h301_mark("head front(0)");
  auto tf = clk::now();
  double h402_front_ms = 0;


  struct ProdSync { std::mutex m; std::condition_variable cv; int ready = 0; } PSY;
  Times TP;
  std::thread prod;
  auto wait_ready = [&](int c) { std::unique_lock<std::mutex> lk(PSY.m); PSY.cv.wait(lk, [&] { return PSY.ready > c; }); };
  if (g_h217) {
    (void)h182_on();
    SB.UQ.U.reserve((size_t)2 * n + 2); SB.UQ.hash.reserve((size_t)2 * n + 2); SB.UQ.enc.reserve((size_t)2 * n + 2);
    PS.rep.reserve((size_t)2 * n + 2); PS.fid.reserve((size_t)2 * n + 2);
    prod = std::thread([&] {
      tl_pool = g_pool2; tl_W = &g_W2;
      for (int c = 0; c < Kn; c++) {


        if (c >= 2) { std::unique_lock<std::mutex> lk(PSY.m); PSY.cv.wait(lk, [&] { return BT[c - 2].sdone && BT[c - 2].gbegun && (c < 3 || BT[c - 3].pdone); }); }
        front(c, TP);
        { std::lock_guard<std::mutex> lk(PSY.m); PSY.ready = c + 1; }
        PSY.cv.notify_all();
      }
      tl_pool = nullptr; tl_W = nullptr;
    });
    wait_ready(0);
  } else {
    if (h650_is("front")) h650_spin_ns(g_h650ms * 1e6 / Kn);
    { const auto f0 = clk::now(); front(0, T); h402_front_ms += ms_since(f0); }
    if (!g_h313_defer && Kn == 1 && g_h301_after_front) { g_h301_after_front(); g_h301_after_front = nullptr; }
  }
  T.head += ms_since(tf);
  for (int c = 0; c < Kn; c++) {
    if (c == 0 && g_h826_rest) { auto rk = std::move(g_h826_rest); g_h826_rest = nullptr; rk(); }
    h301_mark("front(c+1)");
    if (h650_is("front") && !g_h217 && c + 1 < Kn) h650_spin_ns(g_h650ms * 1e6 / Kn);
    if (!g_h217 && c + 1 < Kn && !BT[c + 1].fdone) { const auto f0 = clk::now(); const double e0 = h897_ms; front(c + 1, T); const double fw = ms_since(f0) - (h897_ms - e0); h402_front_ms += fw; BT[c + 1].fdone = true; BT[c + 1].fms = fw; }
    if (!g_h313_defer && !g_h217 && c + 1 == Kn - 1 && g_h301_after_front) { g_h301_after_front(); g_h301_after_front = nullptr; }
    if (g_h217 && g_h302 && c > 0) { h301_mark("wait front(c)"); wait_ready(c); }
    if (!BT[c].egn) { h301_mark("scan_finish wait");
    scan_finish(E, PS, BT[c], T); }
    if (g_h217) {
      { std::lock_guard<std::mutex> lk(PSY.m); }
      PSY.cv.notify_all();
      if (!g_h302 && c + 1 < Kn) wait_ready(c + 1);
    }
    if (c >= 2 && !BT[c - 2].gdone) { fprintf(stderr, "pipeline: GN slot reuse at chunk %d\n", c); exit(1); }
    if (!BT[c].egn) h301_mark("gn_begin encode");
    if (h650_is("enc")) h650_spin_ns(g_h650ms * 1e6 / Kn);
    if (!BT[c].egn) { if (g_h93) { gn_begin_p(E, PS, BT[c], T);  } else gn_begin(E, PS, BT[c], T); }
    if (h650_is("gpumid")) { h680_init(E.dev); h680_spin_cb(E.cq, g_h650ms / Kn); }
    if (g_h217) { { std::lock_guard<std::mutex> lk(PSY.m); BT[c].gbegun = true; } PSY.cv.notify_all(); }
    if (g_kahead && !g_freet && !g_h217 && c >= 2 && c + 2 < Kn && !BT[c + 2].fdone && BT[c].sdone && BT[c - 2].gdone && BT[c - 2].gn_gpu > 0 && BT[c + 1].fms > 0) {


      const Batch& A = BT[c - 2];
      const double gest = A.gn_gpu / (double)(A.b1 - A.b0) * (double)(BT[c].b1 - BT[c].b0);
      const double fest = BT[c + 1].fms / (double)(BT[c + 1].b1 - BT[c + 1].b0) * (double)(BT[c + 2].b1 - BT[c + 2].b0);
      const bool gpend = !g_h882q->is_done(BT[c - 1].t882G) || (BT[c - 1].cbG && BT[c - 1].cbG.status < MTLCommandBufferStatusCompleted);
      if (g_kverb) printf("K_AHEAD c=%d: gpend %d fest %.2f gest %.2f\n", c, (int)gpend, fest, gest);
      if (gpend && fest <= gest) { h301_mark("front(c+2) ahead"); const auto f0 = clk::now(); front(c + 2, T); BT[c + 2].fdone = true; BT[c + 2].fms = ms_since(f0); }
    }
    if (c > 0) { h301_mark("gn_finish wait"); gn_finish(E, PS, BT[c - 1], T); h301_mark("glue+post prep"); glue_range(PS, BT[c - 1].b0, BT[c - 1].b1, err, rx, ry, T); if (tiny_early) tq.push_back({BT[c - 1].b0, BT[c - 1].b1}); post_chunk(BT[c - 1]); pre_chunk(BT[c - 1]);
      if (g_h217) { { std::lock_guard<std::mutex> lk(PSY.m); BT[c - 1].pdone = true; } PSY.cv.notify_all(); } }
  }
  g_tfill = nullptr;
  if (g_h217) { prod.join(); addT(T, TP); }
  if (g_h897 & 2) printf("H897: early GN chunks %d of %d, %.3f ms\n", h897_n, Kn, h897_ms);
  if (!g_h313_defer && g_h301_after_front) { g_h301_after_front(); g_h301_after_front = nullptr; }
  h301_mark("last gn_finish wait");
  gn_finish(E, PS, BT[Kn - 1], T);
  if (g_h242 && g_h413) { h301_mark("H413 put drain"); h242_drain(); }
  h301_mark("last glue+post prep");
  glue_range(PS, BT[Kn - 1].b0, BT[Kn - 1].b1, err, rx, ry, T);
  post_chunk(BT[Kn - 1]);
  pre_chunk(BT[Kn - 1]);
  if (g_ovl_hook && g_freet && g_freet_async && Kn >= 2) { auto hk = std::move(g_ovl_hook); g_ovl_hook = nullptr; h301_mark("STEP_OVL early select"); freet_drain_prior(); hk(BT[Kn - 1].b0); }
  if (g_freet) g_ovl_hook = nullptr;
  if (g_freet_async) freet_flush(T);
  if (h650_is("post")) h650_spin_ns(g_h650ms * 1e6);
  if (g_h313_defer && h84c && g_h301_after_front) { g_h301_after_front(); g_h301_after_front = nullptr; }
  h882_drain(); h882_report();
  h301_mark("union/post fin setup");

  for (int s = 0; s < 2; s++) SB.ids[s] = PS.allids[s];
  T.nuniq += (double)SB.UQ.U.size();
  if (g_h313dbg) { size_t np = 0, a2 = 0, a1s = 0, a1g = 0; for (const Tree& t : SB.UQ.U) { np++; a2 += t.progp != nullptr; } for (const TreeEnc& e : SB.UQ.enc) { a1s += e.sp4p != nullptr; a1g += e.gp4p != nullptr; }
    printf("H313 step trees %zu | progp (A2) %zu | sp4p (A1) %zu | gp4p (A1) %zu\n", np, a2, a1s, a1g); }
  {
    std::sort(PS.iv.begin(), PS.iv.end()); double busy = 0, a = -1, b = -1;
    for (auto& q : PS.iv) { if (q.first > b) { if (b >= a) busy += b - a; a = q.first; b = q.second; } else b = std::max(b, q.second); }
    if (b >= a && a >= 0) busy += b - a;
    T.gpu_union += busy * 1e3; g_h242d_union = busy * 1e3;
    if (g_h648) { const double mid = ms_since(tw), idle = mid - busy * 1e3; g_h648_idle = idle; g_h648_mid = mid;
      if (mid > 1) g_h648_s = std::min(1.5, std::max(0.0, g_h648_s * (1.0 - g_h648k * (idle - g_h648t) / mid))); } }
  if (g_freet) { static int frst = 0; freet_step_print(++frst); }
  if (!g_stage8_only && !g_freet) {
    if (g_post16 == 0) post14(E, PS, SB, n, err, rx, ry, T);
    else if (g_post16 == 2) {
      if (h91c) { for (int c = std::max(0, Kn - n91); c < Kn; c++) post16_h91_col(E, g_h91s[(size_t)(c % n91)], SB, err, rx, ry, T); }
      else if (h84c) { if (h117c) post16_h117_fin(E, PS, SB, err, rx, ry, T, (g_h117 & 1) == 0); else post16_h53_fin(E, PS, SB, err, rx, ry, T); }
      else {
        std::vector<int> all((size_t)n); for (int i = 0; i < n; i++) all[(size_t)i] = i;
        if (g_h53) post16_h53(E, PS, SB, all, err, rx, ry, T); else post16(E, PS, SB, all, err, rx, ry, T);
      }
    }
  }
  g_ovl_hook = nullptr;
  if (g_h313_defer && g_h301_after_front) { g_h301_after_front(); g_h301_after_front = nullptr; }
  if (g_h891 && g_h891_tail_us > 0) g_h891_tail_ms += h891_main_spin_ms(g_h891_tail_us);
  T.pipe_wall += ms_since(tw);
  if (g_h402 && n > 0 && h402_front_ms > 0) {
    const double f = h402_front_ms / n, gs = T.scan_gpu / n, gg = T.gn_gpu / n;
    const double r = (gs + std::sqrt(gs * gs + 4.0 * f * gg)) / (2.0 * f);
    if (std::isfinite(r)) g_h402_rho = std::min(4.0, std::max(1.0, r));
    if (g_h402v) printf("H402: rates f %.3f gS %.3f gG %.3f us/child -> rho %.3f\n", 1e3 * f, 1e3 * gs, 1e3 * gg, g_h402_rho); }
  k_flush();
  tiny_score(E, PS, SB, err, n, admit);
  if (g_sp_env || g_sp_diag) sp_fill_k(E, PS, SB, n, admit, M);
  if (fitx::on()) { const auto fx0 = clk::now(); fit_capture_fixed(E, PS, SB, err, n, admit); const double jcap = ms_since(fx0); fitx::commit(); fitx::g_ms += ms_since(fx0); if (g_jprof) printf("J_PROF capture: fit_capture_fixed %.2f commit+rest %.2f ms (n %d)\n", jcap, ms_since(fx0) - jcap, n); fit_chk_print(g_freet ? 1 : 2); }
  g_fit_span = nullptr; h265_bal_update(PS, bal_span);
}


static void memo_counts(const Engine& E, long& h, long& m, long& d) {
  h = m = 0; for (const std::vector<Worker>* ws : {&E.W, (const std::vector<Worker>*)&g_W2}) for (const Worker& w : *ws) { h += w.H->hits; m += w.H->misses; }
  d = (long)E.SH->drops.load();
}


static double g_h129_cap_mb = 768;
static long g_h129_clears = 0;
static double mb_of(double bytes) { return bytes / 1048576.0; }
static double rss_mb_now() {
  mach_task_basic_info_data_t info; mach_msg_type_number_t cnt = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &cnt) != KERN_SUCCESS) return -1;
  return mb_of((double)info.resident_size);
}
static double phys_fp_mb() {
  task_vm_info_data_t vi; mach_msg_type_number_t cnt = TASK_VM_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vi, &cnt) != KERN_SUCCESS) return -1;
  return mb_of((double)vi.phys_footprint);
}
static double cnt_clamp(uint64_t v, uint64_t cap) { return (double)std::min<uint64_t>(v, cap); }


static double memo_var_bytes(const Engine& E) {
  const double ne = cnt_clamp(E.CS->nent.load(), cc14::Shared::MAXENT);
  double b = (double)E.SH->bump.load();
  b += ne * sizeof(cc14::Ent) + cnt_clamp(E.CS->nkey.load(), cc14::Shared::MAXKEY) * 4.0
     + cnt_clamp(E.CS->nprog.load(), cc14::Shared::MAXPROG) * 4.0 + cnt_clamp(E.CS->nsp.load(), cc14::Shared::MAXSP) * 2.0;
  b += ne * sizeof(cmi6::Slot) + (double)E.MM->bump.load() * 4.0;
  return b;
}


static void memo_clear(Engine& E) {
  const uint32_t nent_cc = std::min<uint32_t>(E.CS->nent.load(), cc14::Shared::MAXENT);
  for (uint32_t i = 0; i < nf16::TB; i++) E.SH->tab[i].store(0, std::memory_order_relaxed);
  E.SH->bump.store(0); E.SH->nent.store(0);
  for (uint32_t i = 0; i < cc14::Shared::TB; i++) E.CS->tab[i].store(0, std::memory_order_relaxed);
  E.CS->nent.store(0); E.CS->nkey.store(0); E.CS->nprog.store(0); E.CS->nsp.store(0);
  if (nent_cc) memset((void*)E.MM->slot, 0, sizeof(cmi6::Slot) * nent_cc);
  E.MM->bump.store(0); E.MM->nst.store(0);
  for (std::vector<Worker>* ws : {&E.W, &g_W2}) for (Worker& w : *ws) {
    w.H->r_next = 0; w.H->r_end = 0;
    cc14::Handle& C = w.CH->C;
    C.r_ent_next = C.r_ent_end = 0; C.r_key_next = C.r_key_end = 0; C.r_prog_next = C.r_prog_end = 0; C.r_sp_next = C.r_sp_end = 0;
  }
  g_h129_clears++;
}

static void memo_line(const Engine& E, int st, double step_ms) {
  printf("H129 step %2d | step ms %.1f | RSS %.0f MB | phys %.0f MB | Metal alloc %.0f MB | NF entries %u arena %.1f MB | cc14 entries %u key %.1f prog %.1f sp %.1f MB | cmi6 arena %.1f MB | memo var %.1f MB | clears %ld\n",
         st, step_ms, rss_mb_now(), phys_fp_mb(), mb_of((double)[E.dev currentAllocatedSize]),
         E.SH->nent.load(), mb_of((double)E.SH->bump.load()),
         E.CS->nent.load(), mb_of(cnt_clamp(E.CS->nkey.load(), cc14::Shared::MAXKEY) * 4.0),
         mb_of(cnt_clamp(E.CS->nprog.load(), cc14::Shared::MAXPROG) * 4.0), mb_of(cnt_clamp(E.CS->nsp.load(), cc14::Shared::MAXSP) * 2.0),
         mb_of((double)E.MM->bump.load() * 4.0), mb_of(memo_var_bytes(E)), g_h129_clears);
}


struct H156Snap { long nf_h = 0, nf_m = 0, cc_h = 0, cc_m = 0, enc_h = 0, enc_m = 0, enc_l = 0; };
static H156Snap g_h156_s0;
static h156::RetainStats g_h156_ev;
static int g_h156_fallback = 0;
static long g_h156_compacts = 0;

static const int g_h886 = getenv("H886") ? atoi(getenv("H886")) : 0;
static std::atomic<int> g_h886_boost{0};
static std::atomic<int> g_h886_qos_seen{0};
static thread_local bool tl_h886_boosted = false;
static inline uint64_t h886_nq() { static const uint64_t n = (g_h886 & 8) ? std::max<uint64_t>(4, std::min<uint64_t>(8, std::thread::hardware_concurrency())) : 4; return n; }
static double g_h886_sl[4][3];
static inline void h886_poll() {
  if (!(g_h886 & 2) || tl_h886_boosted || !g_h886_boost.load(std::memory_order_relaxed)) return;
  tl_h886_boosted = true;
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED, 0);
  qos_class_t q; int rp; pthread_get_qos_class_np(pthread_self(), &q, &rp); g_h886_qos_seen.store((int)q, std::memory_order_relaxed);
}
struct H886Pend { h156::RetainStats ev; int fb = 0; long cmp = 0; bool valid = false; };
static H886Pend g_h886_pend; static std::mutex g_h886_mx;
static H156Snap h156_snap(const Engine& E) {
  H156Snap s;
  for (const std::vector<Worker>* ws : {&E.W, (const std::vector<Worker>*)&g_W2}) for (const Worker& w : *ws) {
    s.nf_h += w.H->hits; s.nf_m += w.H->misses;
    s.cc_h += w.CH->C.hits; s.cc_m += w.CH->C.misses;
    s.enc_h += w.CH->enc_hit; s.enc_m += w.CH->enc_miss; s.enc_l += w.CH->enc_local;
  }
  return s;
}
static void h156_step_begin(const Engine& E, int st) {
  h156::g_touch = (uint32_t)st;
  g_h156_s0 = h156_snap(E);
}


static void h874_pretouch(size_t n) {
  const size_t cap = 2 * n + 2; const auto t0 = clk::now(); size_t tot = 0;


  static const int mem_pre = getenv("MEM_PRETOUCH") ? atoi(getenv("MEM_PRETOUCH")) : 1;
  auto touch = [&](auto& v) { using T = typename std::decay<decltype(v)>::type::value_type; v.reserve(cap); if (!mem_pre) memset((void*)v.data(), 0, cap * sizeof(T)); tot += mem_pre ? 0 : cap * sizeof(T); };
  StepBatch* sp = new StepBatch();
  touch(sp->UQ.U); touch(sp->UQ.enc); touch(sp->UQ.hash);
  sp->ct.assign(2 * n, -1); sp->F[0].init((int)cap); sp->F[1].init((int)cap);
  g_h415_spare = sp;
  touch(g_h301r_t); touch(g_h301r_e);
  printf("H874: pretouch %.1f MB in %.1f ms (cap %zu)\n", (double)tot / 1048576.0, std::chrono::duration<double, std::milli>(clk::now() - t0).count(), cap);
}
static long h874_mf() { struct rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_minflt; }
static long g_h874_mf0 = -1;


static const int g_h853 = getenv("H853") ? atoi(getenv("H853")) : 0;
static const int g_prefault = getenv("PREFAULT") ? atoi(getenv("PREFAULT")) : 0;
static void h853_prefault(Engine& E) {
  if (!g_h853 || !(g_h129_cap_mb > 0)) return;
  const double X = getenv("H853_X") ? atof(getenv("H853_X")) : 1.1;
  const double cap = g_h129_cap_mb * X * 1048576.0;
  const auto t0 = clk::now(); size_t tot = 0;
  auto touch = [&](void* base, size_t from, double want, size_t limit) {
    const size_t to = std::min(limit, (size_t)want);
    if (to > from) { memset((char*)base + from, 0, to - from); tot += to - from; }
  };
  if (g_h853 & 1) {
    touch(E.SH->ar, (size_t)16 << 20, cap * 0.37, nf16::ARENA_BYTES);
    touch(E.MM->arena, 0, cap * 0.25, cmi6::ARENA_WORDS * sizeof(uint32_t));
    touch(E.MM->slot, 0, cap * 0.18, (size_t)cc14::Shared::MAXENT * sizeof(cmi6::Slot));
  }
  if (g_h853 & 2) {
    touch(E.CS->ent, 0, cap * 0.077, (size_t)cc14::Shared::MAXENT * sizeof(cc14::Ent));
    touch(E.CS->key, 0, cap * 0.087, cc14::Shared::MAXKEY * sizeof(uint32_t));
    touch(E.CS->prog, 0, cap * 0.037, cc14::Shared::MAXPROG * sizeof(uint32_t));
    touch(E.CS->sp, 0, cap * 0.0055, cc14::Shared::MAXSP * sizeof(uint16_t));
  }
  printf("H853: %d | prefault %.1f MB in %.1f ms (cap %.0f MB x %.2f)\n", g_h853, (double)tot / 1048576.0,
         std::chrono::duration<double, std::milli>(clk::now() - t0).count(), g_h129_cap_mb, X);
}

static const int g_h854 = getenv("H854") ? atoi(getenv("H854")) : 0;
static const int g_h854_qos = getenv("H854_QOS") ? atoi(getenv("H854_QOS")) : 2;
struct H854Rec { double launch = 0, start = 0, end = 0, hist = 0, nfc = 0, ccc = 0; bool ran = false; };
static H854Rec g_h854_rec;
static double g_h854_join_wait = 0, g_h854_join_done = 0; static bool g_h854_joined = false;
static struct rusage g_h854_ru0; static bool g_h854_ru0_ok = false;
static void h854_qos_self() {
  if (!(g_h854 & 2)) return;
  pthread_set_qos_class_self_np(g_h854_qos == 1 ? QOS_CLASS_BACKGROUND : g_h854_qos == 3 ? QOS_CLASS_DEFAULT : QOS_CLASS_UTILITY, 0);
}

template <class F> static void h854_par_ranges(uint64_t n, F f) {
  const uint64_t NQ = h886_nq(); std::vector<std::thread> th(NQ - 1);
  auto fc = [&](uint64_t lo, uint64_t hi) {
    if (!(g_h886 & 2)) { f(lo, hi); return; }
    for (uint64_t a = lo; a < hi; ) { const uint64_t b = std::min<uint64_t>(hi, a + (1u << 16)); h886_poll(); f(a, b); a = b; }
  };
  for (uint64_t q = 0; q < NQ - 1; q++) th[q] = std::thread([&, q] { fc(n * q / NQ, n * (q + 1) / NQ); });
  fc(n * (NQ - 1) / NQ, n);
  for (uint64_t q = 0; q < NQ - 1; q++) th[q].join();
}

static void h854_cc_live(cc14::Shared& C, uint32_t ne, std::vector<uint8_t>& live) {
  live.assign(ne, 0);
  h854_par_ranges(cc14::Shared::TB, [&](uint64_t lo, uint64_t hi) {
    for (uint64_t p = lo; p < hi; p++) { const uint32_t v = C.tab[p].load(std::memory_order_relaxed); if (v && v - 1 < ne) live[v - 1] = 1; }
  });
}

typedef std::vector<std::pair<uint32_t, uint32_t>> H886Rec;
static void h854_compact_nf(nf16::Shared& S, uint32_t thr, h156::RetainStats& out, const H886Rec* rec = nullptr) {
  const uint64_t NQ = h886_nq(); std::vector<std::vector<uint64_t>> kp(NQ); std::vector<long> lv(NQ, 0);
  {
    std::vector<std::thread> th(NQ - 1);
    auto scan = [&](uint64_t q) {
      const uint64_t lo = (uint64_t)nf16::TB * q / NQ, hi = (uint64_t)nf16::TB * (q + 1) / NQ;
      for (uint64_t p = lo; p < hi; p++) {
        if (!(p & 0xFFFF)) h886_poll();
        const uint64_t v = S.tab[p].load(std::memory_order_relaxed);
        if (!v) continue;
        lv[q]++;
        const uint64_t off = (uint64_t)((uint32_t)v - 1u);
        if (((const nf16::BHdr*)(S.ar + off))->touch >= thr) kp[q].push_back(off);
      }
    };
    auto scan2 = [&](uint64_t q) {
      if (rec) { lv[q] = (long)rec[q].size(); for (const auto& pr : rec[q]) if (pr.second >= thr) kp[q].push_back(pr.first); }
      else scan(q);
      if (g_h886 & 32) std::sort(kp[q].begin(), kp[q].end());
    };
    for (uint64_t q = 0; q < NQ - 1; q++) th[q] = std::thread(scan2, q);
    scan2(NQ - 1);
    for (uint64_t q = 0; q < NQ - 1; q++) th[q].join();
  }
  std::vector<uint64_t> keep; long live = 0;
  std::vector<size_t> bnd;
  for (uint64_t q = 0; q < NQ; q++) { live += lv[q]; keep.insert(keep.end(), kp[q].begin(), kp[q].end()); bnd.push_back(keep.size()); }
  if (g_h886 & 32) { for (size_t w = 1; w < bnd.size(); w++) std::inplace_merge(keep.begin(), keep.begin() + bnd[w - 1], keep.begin() + bnd[w]); }
  else std::sort(keep.begin(), keep.end());
  h854_par_ranges(nf16::TB, [&](uint64_t lo, uint64_t hi) { for (uint64_t p = lo; p < hi; p++) S.tab[p].store(0, std::memory_order_relaxed); });
  uint64_t dst = 0;
  uint64_t h886c = 0;
  for (size_t ki = 0; ki < keep.size(); ki++) {
    const uint64_t off = keep[ki];
    if (!(++h886c & 0x3FFF)) h886_poll();
    if (g_h886 & 16) {
      if (ki + 16 < keep.size()) __builtin_prefetch(S.ar + keep[ki + 16]);
      if (ki + 8 < keep.size()) { const nf16::BHdr* h8 = (const nf16::BHdr*)(S.ar + keep[ki + 8]); __builtin_prefetch((const void*)&S.tab[(uint32_t)(skn::NFMemo::hkey(h8->key) & (nf16::TB - 1))], 1); }
    }
    const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + off);
    const uint64_t bytes = h156::nf_block_bytes(hd);
    const uint64_t h = skn::NFMemo::hkey(hd->key);
    if (dst != off) memmove(S.ar + dst, S.ar + off, (size_t)bytes);
    uint32_t p = (uint32_t)(h & (nf16::TB - 1));
    while (S.tab[p].load(std::memory_order_relaxed) != 0) p = (p + 1) & (nf16::TB - 1);
    S.tab[p].store(((uint64_t)(uint32_t)(h >> 32) << 32) | (dst + 1), std::memory_order_relaxed);
    dst += bytes;
  }
  S.bump.store(dst, std::memory_order_relaxed);
  S.nent.store((uint32_t)keep.size(), std::memory_order_relaxed);
  out.nf_ev = live - (long)keep.size();
  out.nf_kept = (long)keep.size();
}

static void h854_compact_cc(cc14::Shared& C, cmi6::Memo& M, uint32_t thr, const std::vector<uint8_t>& live, h156::RetainStats& out) {
  const double q0 = h301_now();
  const uint32_t ne = std::min<uint32_t>(C.nent.load(std::memory_order_relaxed), cc14::Shared::MAXENT);
  const double q1 = h301_now();
  long nlive = 0, nlive_enc = 0;
  std::vector<uint32_t> keep;
  for (uint32_t i = 0; i < ne; i++) {
    if (!(i & 0x3FFF)) h886_poll();
    if (!live[i]) continue;
    nlive++;
    if (M.slot[i].st == 2) nlive_enc++;
    if (C.ent[i].touch >= thr) keep.push_back(i);
  }
  const uint32_t K = (uint32_t)keep.size();
  std::vector<uint32_t> nk_off(ne, 0), np_off(ne, 0), ns_off(ne, 0), nt_off(ne, 0);
  auto slide = [](int id, const std::vector<uint32_t>& lst, auto offf, auto cntf, auto* base, std::vector<uint32_t>& nw) -> uint64_t {
    const double s0 = h301_now(); double mvb = 0;
    std::vector<uint32_t> o(lst);
    auto lt = [&](uint32_t a, uint32_t b) { return offf(a) < offf(b); };
    if (!((g_h886 & 64) && std::is_sorted(o.begin(), o.end(), lt))) std::sort(o.begin(), o.end(), lt);
    const double s1 = h301_now();
    uint64_t dst = 0;
    uint64_t c = 0;
    for (uint32_t i : o) {
      if (!(++c & 0x3FFF)) h886_poll();
      const uint64_t src = offf(i), n = cntf(i);
      nw[i] = (uint32_t)dst;
      if (dst != src && n) { memmove(base + dst, base + src, (size_t)(n * sizeof(*base))); mvb += (double)(n * sizeof(*base)); }
      dst += n;
    }
    g_h886_sl[id][0] = (s1 - s0) * 1e3; g_h886_sl[id][1] = (h301_now() - s1) * 1e3; g_h886_sl[id][2] = mvb;
    return dst;
  };
  std::vector<uint32_t> enc_keep;
  for (uint32_t i : keep) if (M.slot[i].st == 2) enc_keep.push_back(i);
  uint64_t tot_k = 0, tot_p = 0, tot_s = 0, tot_t = 0;
  std::thread t1([&] { tot_k = slide(0, keep, [&](uint32_t i) { return (uint64_t)C.ent[i].key_off; }, [&](uint32_t i) { return (uint64_t)C.ent[i].nkey; }, C.key, nk_off); });
  std::thread t2([&] { tot_p = slide(1, keep, [&](uint32_t i) { return (uint64_t)C.ent[i].prog_off; }, [&](uint32_t i) { return (uint64_t)C.ent[i].nprog; }, C.prog, np_off); });
  std::thread t3([&] { tot_s = slide(2, keep, [&](uint32_t i) { return (uint64_t)C.ent[i].sp_off; }, [&](uint32_t i) { return (uint64_t)C.ent[i].nsp; }, C.sp, ns_off); });
  tot_t = slide(3, enc_keep, [&](uint32_t i) { return (uint64_t)M.slot[i].sp_off; }, [&](uint32_t i) { return (uint64_t)M.slot[i].nsp4 + M.slot[i].ngp4; }, M.arena, nt_off);
  t1.join(); t2.join(); t3.join();
  const double q2 = h301_now();
  for (uint32_t r = 0; r < K; r++) {
    const uint32_t i = keep[r];
    if (r != i) { C.ent[r] = C.ent[i]; M.slot[r] = M.slot[i]; }
    cc14::Ent& e = C.ent[r];
    e.key_off = nk_off[i]; e.prog_off = np_off[i]; e.sp_off = ns_off[i];
    cmi6::Slot& s = M.slot[r];
    if (s.st == 2) { s.sp_off = nt_off[i]; s.gp_off = nt_off[i] + s.nsp4; }
  }
  if (ne > K) memset((void*)(M.slot + K), 0, sizeof(cmi6::Slot) * (ne - K));
  h854_par_ranges(cc14::Shared::TB, [&](uint64_t lo, uint64_t hi) { for (uint64_t p = lo; p < hi; p++) C.tab[p].store(0, std::memory_order_relaxed); });
  for (uint32_t r = 0; r < K; r++) {
    if (!(r & 0x3FFF)) h886_poll();
    if ((g_h886 & 16) && r + 8 < K) __builtin_prefetch((const void*)&C.tab[(uint32_t)(C.ent[r + 8].h & (cc14::Shared::TB - 1))], 1);
    uint32_t p = (uint32_t)(C.ent[r].h & (cc14::Shared::TB - 1));
    while (C.tab[p].load(std::memory_order_relaxed) != 0) p = (p + 1) & (cc14::Shared::TB - 1);
    C.tab[p].store(r + 1, std::memory_order_relaxed);
  }
  C.nent.store(K, std::memory_order_relaxed);
  C.nkey.store((uint32_t)tot_k, std::memory_order_relaxed);
  C.nprog.store((uint32_t)tot_p, std::memory_order_relaxed);
  C.nsp.store((uint32_t)tot_s, std::memory_order_relaxed);
  M.bump.store(tot_t, std::memory_order_relaxed);
  long enc_kept = 0;
  for (uint32_t r = 0; r < K; r++) if (M.slot[r].st == 2) enc_kept++;
  M.nst.store((uint32_t)enc_kept, std::memory_order_relaxed);
  out.cc_ev = nlive - (long)K;
  out.cc_kept = (long)K;
  out.enc_ev = nlive_enc - enc_kept;
  out.enc_kept = enc_kept;
  if (g_h854 & 1) printf("H854 cc compact: live scan %.1f | keep+slides %.1f ... entries+table+rest %.1f ms (ne %u K %u) | slides sort/move ms, MB: key %.1f/%.1f/%.0f prog %.1f/%.1f/%.0f sp %.1f/%.1f/%.0f enc %.1f/%.1f/%.0f\n", (q1 - q0) * 1e3, (q2 - q1) * 1e3, (h301_now() - q2) * 1e3, ne, K,
    g_h886_sl[0][0], g_h886_sl[0][1], g_h886_sl[0][2] / 1048576.0, g_h886_sl[1][0], g_h886_sl[1][1], g_h886_sl[1][2] / 1048576.0, g_h886_sl[2][0], g_h886_sl[2][1], g_h886_sl[2][2] / 1048576.0, g_h886_sl[3][0], g_h886_sl[3][1], g_h886_sl[3][2] / 1048576.0);
}

static void h854_retain(nf16::Shared& S, cc14::Shared& C, cmi6::Memo& M, uint32_t st, double cap_bytes, h156::RetainStats& out) {
  const double ta = h301_now();
  const uint32_t NT = st + 2;
  std::vector<double> hist(NT, 0.0);
  const uint32_t ne = std::min<uint32_t>(C.nent.load(std::memory_order_relaxed), cc14::Shared::MAXENT);
  std::vector<std::pair<uint32_t, double>> ccadd;
  auto nf_hist = [&]() {
    for (uint32_t p = 0; p < nf16::TB; p++) {
      const uint64_t v = S.tab[p].load(std::memory_order_relaxed);
      if (!v) continue;
      const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + (uint64_t)((uint32_t)v - 1u));
      hist[h156::gen_of(hd->touch, NT)] += (double)h156::nf_block_bytes(hd);
    }
  };
  std::vector<H886Rec> nfrec(h886_nq()); const bool h886rec = (g_h886 & 4) && (g_h854 & 12) == 12 && (g_h854 & 8);
  std::vector<uint8_t> live;
  auto cc_collect = [&]() {
    if (g_h854 & 8) h854_cc_live(C, ne, live);
    else {
      live.assign(ne, 0);
      for (uint32_t p = 0; p < cc14::Shared::TB; p++) {
        const uint32_t v = C.tab[p].load(std::memory_order_relaxed);
        if (v && v - 1 < ne) live[v - 1] = 1;
      }
    }
    for (uint32_t i = 0; i < ne; i++) {
      if (!(i & 0x3FFF)) h886_poll();
      if (live[i]) ccadd.emplace_back(h156::gen_of(C.ent[i].touch, NT), h156::cc_entry_bytes(C.ent[i], M.slot[i]));
    }
  };
  if ((g_h854 & 12) == 12) {
    const uint32_t NQ = (uint32_t)h886_nq(); std::vector<std::vector<std::pair<uint32_t, double>>> add(NQ); std::vector<std::thread> th(NQ);
    for (uint32_t q = 0; q < NQ; q++) th[q] = std::thread([&, q] {
      const uint32_t lo = (uint32_t)((uint64_t)nf16::TB * q / NQ), hi = (uint32_t)((uint64_t)nf16::TB * (q + 1) / NQ);
      for (uint32_t p = lo; p < hi; p++) {
        if (!(p & 0xFFFF)) h886_poll();
        if ((g_h886 & 16) && p + 16 < hi) { const uint64_t w = S.tab[p + 16].load(std::memory_order_relaxed); if (w) __builtin_prefetch(S.ar + (uint64_t)((uint32_t)w - 1u)); }
        const uint64_t v = S.tab[p].load(std::memory_order_relaxed);
        if (!v) continue;
        const nf16::BHdr* hd = (const nf16::BHdr*)(S.ar + (uint64_t)((uint32_t)v - 1u));
        add[q].emplace_back(h156::gen_of(hd->touch, NT), (double)h156::nf_block_bytes(hd));
        if (h886rec) nfrec[q].emplace_back((uint32_t)v - 1u, (uint32_t)hd->touch);
      }
    });
    cc_collect();
    for (uint32_t q = 0; q < NQ; q++) { th[q].join(); for (const auto& pr : add[q]) hist[pr.first] += pr.second; }
  }
  else if (g_h854 & 4) { std::thread th(cc_collect); nf_hist(); th.join(); }
  else { nf_hist(); cc_collect(); }
  for (const auto& pr : ccadd) hist[pr.first] += pr.second;
  const double target = 0.5 * cap_bytes;
  uint32_t thr = st > 0 ? st - 1 : 0;
  double cum = 0;
  for (int64_t g = (int64_t)st - 1; g >= 0; g--) {
    cum += hist[(size_t)g];
    if (cum > target) break;
    thr = (uint32_t)g;
  }
  out.thr = thr;
  const double tb = h301_now();
  double tc = tb, td = tb;
  if (g_h854 & 4) {
    double tn = 0;
    std::thread th([&] { h854_qos_self(); const double t0 = h301_now(); if (g_h854 & 8) h854_compact_nf(S, thr, out, h886rec ? nfrec.data() : nullptr); else h156::compact_nf(S, thr, out); tn = h301_now() - t0; });
    if (g_h854 & 8) h854_compact_cc(C, M, thr, live, out); else h156::compact_cc(C, M, thr, out);
    td = h301_now(); th.join(); tc = tb + tn;
    g_h854_rec.nfc = tn * 1e3; g_h854_rec.ccc = (td - tb) * 1e3;
  } else {
    h156::compact_nf(S, thr, out); tc = h301_now();
    h156::compact_cc(C, M, thr, out); td = h301_now();
    g_h854_rec.nfc = (tc - tb) * 1e3; g_h854_rec.ccc = (td - tc) * 1e3;
  }
  g_h854_rec.hist = (tb - ta) * 1e3;
}

struct H886Pub { bool on; H886Pend p; ~H886Pub() { if (on) { std::lock_guard<std::mutex> lk(g_h886_mx); g_h886_pend.ev = p.ev; g_h886_pend.fb = p.fb; g_h886_pend.cmp += p.cmp; g_h886_pend.valid = true; } } };
static void memo_boundary(Engine& E, int st) {
  H886Pub h886pub; h886pub.on = (g_h886 & 1) != 0;
  h156::RetainStats& g_h156_ev = h886pub.on ? h886pub.p.ev : ::g_h156_ev; int& g_h156_fallback = h886pub.on ? h886pub.p.fb : ::g_h156_fallback; long& g_h156_compacts = h886pub.on ? h886pub.p.cmp : ::g_h156_compacts;
  g_h156_ev = h156::RetainStats(); g_h156_fallback = 0;
  if (g_h854 & 1) { g_h854_rec = H854Rec(); }
  if (!(g_h129_cap_mb > 0 && mb_of(memo_var_bytes(E)) > g_h129_cap_mb)) return;
  if (!h156::g_on) {
    g_h156_ev.nf_ev = (long)E.SH->nent.load(); g_h156_ev.cc_ev = (long)E.CS->nent.load(); g_h156_ev.enc_ev = (long)E.MM->nst.load();
    memo_clear(E); g_h156_fallback = 1;
    return;
  }
  if (g_h854 & 13) { g_h854_rec.ran = true; g_h854_rec.start = h301_now(); h854_retain(*E.SH, *E.CS, *E.MM, (uint32_t)st, g_h129_cap_mb * 1048576.0, g_h156_ev); }
  else h156::retain(*E.SH, *E.CS, *E.MM, (uint32_t)st, g_h129_cap_mb * 1048576.0, g_h156_ev);
  g_h156_compacts++;


  for (std::vector<Worker>* ws : {&E.W, &g_W2}) for (Worker& w : *ws) {
    w.H->r_next = 0; w.H->r_end = 0;
    cc14::Handle& C = w.CH->C;
    C.r_ent_next = C.r_ent_end = 0; C.r_key_next = C.r_key_end = 0; C.r_prog_next = C.r_prog_end = 0; C.r_sp_next = C.r_sp_end = 0;
  }
  if (mb_of(memo_var_bytes(E)) > g_h129_cap_mb) {
    g_h156_ev.nf_ev += (long)E.SH->nent.load(); g_h156_ev.cc_ev += (long)E.CS->nent.load(); g_h156_ev.enc_ev += (long)E.MM->nst.load();
    memo_clear(E); g_h156_fallback = 1;
  }
  if (g_h854 & 13) g_h854_rec.end = h301_now();
}
static H854Rec g_h886_rec;
static void h886_adopt() {
  if (!(g_h886 & 1)) return;
  { std::lock_guard<std::mutex> lk(g_h886_mx); if (g_h886_pend.valid) { g_h156_ev = g_h886_pend.ev; g_h156_fallback = g_h886_pend.fb; g_h156_compacts += g_h886_pend.cmp; g_h886_pend = H886Pend(); } }
  g_h886_rec = g_h854_rec;
}

static void h156_step_print(const Engine& E, int st) {
  const H156Snap s1 = h156_snap(E); const H156Snap& s0 = g_h156_s0;
  const long nh = s1.nf_h - s0.nf_h, nm = s1.nf_m - s0.nf_m;
  const long ch = s1.cc_h - s0.cc_h, cm = s1.cc_m - s0.cc_m;
  const long eh = s1.enc_h - s0.enc_h, em = s1.enc_m - s0.enc_m, el = s1.enc_l - s0.enc_l;
  const long nl = nh + nm, cl = ch + cm, xl = eh + em + el;
  const double ne = cnt_clamp(E.CS->nent.load(), cc14::Shared::MAXENT);
  const double cc_b = ne * sizeof(cc14::Ent) + cnt_clamp(E.CS->nkey.load(), cc14::Shared::MAXKEY) * 4.0
                    + cnt_clamp(E.CS->nprog.load(), cc14::Shared::MAXPROG) * 4.0 + cnt_clamp(E.CS->nsp.load(), cc14::Shared::MAXSP) * 2.0;
  const double enc_b = (double)E.MM->bump.load() * 4.0 + ne * sizeof(cmi6::Slot);
  printf("H156 step %2d | mode %s | NF lookups %ld hits %ld rate %.4f entries %u MB %.1f evicted %ld"
         " | cc14 lookups %ld hits %ld rate %.4f entries %u MB %.1f evicted %ld"
         " | cmi6 lookups %ld hits %ld rate %.4f entries %u MB %.1f evicted %ld"
         " | compactions %ld fallback %d thr %u\n",
         st, h156::g_on ? "retain" : "h129",
         nl, nh, nl ? (double)nh / (double)nl : 0.0, E.SH->nent.load(), mb_of((double)E.SH->bump.load()), g_h156_ev.nf_ev,
         cl, ch, cl ? (double)ch / (double)cl : 0.0, E.CS->nent.load(), mb_of(cc_b), g_h156_ev.cc_ev,
         xl, eh, xl ? (double)eh / (double)xl : 0.0, E.MM->nst.load(), mb_of(enc_b), g_h156_ev.enc_ev,
         g_h156_compacts, g_h156_fallback, g_h156_ev.thr);
}


#include "evo_loop_h33.hpp"
#include "evo_loop_h53.hpp"
#include <thread>
#include <functional>
static const int g_h419 = getenv("H419") ? atoi(getenv("H419")) : 0;
#include "evo_loop_i20_h709.hpp"
#include "evo_loop_h696_h709.hpp"


static void h826_breed(evo::Pool& pool, evo::State& S, const uint32_t* pidx, const uint32_t* cs, evo::BreedBuf& BB, uint64_t seed, uint64_t step,
                       uint32_t* rebreeds, uint32_t* failed, const evo::Cand** chp, uint8_t* cv, int d_lo, int d_hi) {
  using namespace evo;
  if (d_hi <= d_lo) return;
  if (g_h777 && g_h777_gp.size() != (size_t)NPOP) g_h777_gp.assign((size_t)NPOP, nullptr);
  const evo::Cand** Gp = g_h777 ? g_h777_gp.data() : nullptr;
  pool.run(d_hi - d_lo, [&](int k) {
    const int d = d_lo + k;
    if (g_h777) for (int i = d * DECK; i < (d + 1) * DECK; i++) { const Ind& P = XC(S, pidx[i]); Gp[i] = &P.g; if (g_h740) h740_rc(S, pidx[i], BB.prk[i], BB.pcw[i]); else { BB.prk[i] = P.rank; BB.pcw[i] = P.crowd; } }
    else for (int i = d * DECK; i < (d + 1) * DECK; i++) { const Ind& P = XC(S, pidx[i]); BB.G[i] = P.g; if (g_h740) h740_rc(S, pidx[i], BB.prk[i], BB.pcw[i]); else { BB.prk[i] = P.rank; BB.pcw[i] = P.crowd; } }
    Rng r(mix(mix(seed, step), (uint64_t)d));
    uint32_t rb = 0, fl = 0;
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      Ind& out = XA(S, cs[i]);
      const bool ok = g_h777 ? make_child_i12_p(BB, Gp, d * DECK, r, out, rb, g_link320 ? g_link320 + i : nullptr)
                             : make_child_i12(BB, d * DECK, r, out, rb, g_link320 ? g_link320 + i : nullptr);
      if (!ok) fl++;
      if (g_h740) g_h740e[cs[i]].valid = ok ? 1 : 0;
      if (g_sf_axis && g_sf_pin && g_pinfix) sf_pin_driver(out.g);
      chp[i] = &out.g; cv[i] = out.valid;
    }
    rebreeds[d] = rb; failed[d] = fl;
  });
}


static const int g_h320 = getenv("H320") ? atoi(getenv("H320")) : 0;
static const char* g_h320_path = getenv("H320_DUMP");
static const int g_h320_from = getenv("H320_DUMP_FROM") ? atoi(getenv("H320_DUMP_FROM")) : 1;
static std::vector<evo::ChildLink320> g_h320_link;


struct Rec320 { int32_t step, child; uint8_t side, n, lenA, lenB; uint8_t swap, xo, s, szi; uint8_t js, szj, k, pad; uint8_t op[24], opA[24], opB[24]; };
static_assert(sizeof(Rec320) == 92, "Rec320 is 92 bytes");
static FILE* g_h320_f = nullptr;
static void h320_flush_step(int st, const evo::State& S) {
  std::vector<nf308::DumpRec> R; { std::lock_guard<std::mutex> lk(nf308::g_h320_mx); R.swap(nf308::g_h320_recs); }
  std::sort(R.begin(), R.end(), [](const nf308::DumpRec& a, const nf308::DumpRec& b) { return a.child != b.child ? a.child < b.child : a.side < b.side; });
  if (st < g_h320_from || !g_h320_f) return;
  std::vector<Rec320> out(R.size());
  for (size_t q = 0; q < R.size(); q++) {
    const nf308::DumpRec& d = R[q]; Rec320& o = out[q]; std::memset(&o, 0, sizeof o);
    const evo::ChildLink320& L = g_h320_link[(size_t)d.child]; const int t = d.side;
    o.step = st; o.child = d.child; o.side = d.side; o.n = d.n;
    o.swap = L.swap; o.xo = L.xo[t]; o.s = L.s[t]; o.szi = L.szi[t]; o.js = L.js[t]; o.szj = L.szj[t]; o.k = L.k[t];
    std::memcpy(o.op, d.op, 24);
    const evo::Cand& A = g_h777 ? *g_h777_gp[(size_t)L.ia] : S.BB.G[(size_t)L.ia]; const evo::Cand& B = g_h777 ? *g_h777_gp[(size_t)L.ib] : S.BB.G[(size_t)L.ib];
    o.lenA = A.len[t]; std::memcpy(o.opA, A.op[t], 24);
    o.lenB = B.len[t ^ L.swap]; std::memcpy(o.opB, B.op[t ^ L.swap], 24);
  }
  int32_t hd[2] = {st, (int32_t)out.size()};
  fwrite(hd, 4, 2, g_h320_f); if (!out.empty()) fwrite(out.data(), sizeof(Rec320), out.size(), g_h320_f);
  fflush(g_h320_f);
  printf("H320 dump step %d: %zu NF misses written\n", st, out.size());
}
#if 0
static void select_v3(evo::Pool& pool, evo::State& S, const std::vector<Meta>& pm, const std::vector<Meta>& cm,
                      const uint8_t* admit, std::vector<Meta>& nm, Times& T) {
  using namespace evo;
  const int NP = NPOP;
  pool.run(NDEME, [&](int d) { migrants(d, S.par.data(), S.migtab.data() + d * MIGR, S.sel[d].mtmp); });
  std::vector<int> mstart(NDEME, 0);
  pool.run(NDEME, [&](int d) {
    SelScratch& s = S.sel[d];
    s.pidx.clear();
    for (int i = d * DECK; i < (d + 1) * DECK; i++) s.pidx.push_back((uint32_t)i);
    for (int i = d * DECK; i < (d + 1) * DECK; i++) {
      if (!admit[i]) continue;
      Ind& c = S.chi[i];
      float sum = 0.0f; for (int q = 0; q < NRES; q++) sum += c.res[q];
      c.err = sum; c.cx = (uint16_t)cplx(c.g);
      s.pidx.push_back((uint32_t)(NP + i));
    }
    mstart[d] = (int)s.pidx.size();
    int pd = (d + NDEME - 1) % NDEME;
    for (int k = 0; k < MIGR; k++) s.pidx.push_back(S.migtab[pd * MIGR + k]);
    const int m = (int)s.pidx.size();
    s.pe.resize(m); s.pc.resize(m); s.prank.resize(m); s.pcrowd.resize(m); s.surv.resize(m);
    for (int q = 0; q < m; q++) {
      uint32_t p = s.pidx[q];
      const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
      s.pe[q] = src.err; s.pc[q] = (uint16_t)src.cx;
    }
    rank_crowd(m, s.pe.data(), s.pc.data(), s.prank.data(), s.pcrowd.data(), s.rs);
    choose_survivors(m, s.pcrowd.data(), s.rs, s.surv.data(), s.tmp);
  });


  FlatKeySet placed; placed.reserve((size_t)NP);
  auto place = [&](int dst, int slot, int j, uint32_t pos) {
    const SelScratch& sj = S.sel[j];
    const uint32_t p = sj.pidx[pos];
    const Meta& mt = p < (uint32_t)NP ? pm[p] : cm[p - NP];
    const Ind& src = p < (uint32_t)NP ? S.par[p] : S.chi[p - NP];
    Ind& o = S.nxt[dst * DECK + slot];
    o = src; o.rank = sj.prank[pos]; o.crowd = sj.pcrowd[pos]; o.err = sj.pe[pos]; o.valid = 1;
    nm[dst * DECK + slot] = mt;
    if (p >= (uint32_t)NP) T.surv_child += 1;
  };
  std::vector<int> cur(NDEME, 0), fillA(NDEME, 0);
  for (int d = 0; d < NDEME; d++) {
    const SelScratch& s = S.sel[d];
    const int m = (int)s.pidx.size();
    int k = 0, q = 0;
    for (; q < m && k < DECK; q++) {
      const uint32_t pos = s.surv[q];
      const Meta& mt = s.pidx[pos] < (uint32_t)NP ? pm[s.pidx[pos]] : cm[s.pidx[pos] - NP];
      if (!placed.insert(mt.kh, mt.key.data(), (uint32_t)mt.key.size())) { T.sel_skip += 1; if (pos >= (uint32_t)mstart[d]) T.sel_skip_mig += 1; continue; }
      place(d, k, d, pos); k++;
    }
    cur[d] = q; fillA[d] = k;
  }
  for (int d = 0; d < NDEME; d++) {
    int slot = fillA[d];
    for (int j = 0; j < NDEME && slot < DECK; j++) {
      const SelScratch& sj = S.sel[j];
      const int mj = (int)sj.pidx.size();
      while (cur[j] < mj && slot < DECK) {
        const uint32_t pos = sj.surv[cur[j]++];
        const Meta& mt = sj.pidx[pos] < (uint32_t)NP ? pm[sj.pidx[pos]] : cm[sj.pidx[pos] - NP];
        if (!placed.insert(mt.kh, mt.key.data(), (uint32_t)mt.key.size())) { T.sel_skip += 1; if (pos >= (uint32_t)mstart[j]) T.sel_skip_mig += 1; continue; }
        place(d, slot, j, pos); slot++; T.xfill += 1;
      }
    }
    if (slot < DECK) { fprintf(stderr, "select: deme %d short by %d after both passes (impossible: parents are distinct)\n", d, DECK - slot); exit(1); }
  }
}

#endif
struct TreeBuf { uint8_t op[evo::MAXLEN]; float c[evo::MAXLEN]; int len; bool over; };
static void rgrow(evo::Rng& r, int depth, int maxd, bool full, TreeBuf& B) {
  if (B.over) return;
  if (B.len >= evo::MAXLEN) { B.over = true; return; }
  if (depth >= maxd || (!full && r.uni() < 0.3)) {
    B.op[B.len] = r.uni() < 0.5 ? (uint8_t)evo::CONST : (uint8_t)evo::TVAR;
    B.c[B.len] = (float)r.gauss(); B.len++; return;
  }
  if (r.uni() < 0.5) {
    rgrow(r, depth + 1, maxd, full, B); rgrow(r, depth + 1, maxd, full, B);
    if (B.over) return;
    if (B.len >= evo::MAXLEN) { B.over = true; return; }
    B.op[B.len] = (uint8_t)(evo::ADD + r.below(4)); B.c[B.len] = 0.0f; B.len++;
  } else {
    rgrow(r, depth + 1, maxd, full, B);
    if (B.over) return;
    if (B.len >= evo::MAXLEN) { B.over = true; return; }
    B.op[B.len] = (uint8_t)(evo::SIN + r.below(8)); B.c[B.len] = 0.0f; B.len++;
  }
}
static void random_cand(evo::Rng& r, evo::Cand& c) {
  for (;;) {
    std::memset(&c, 0, sizeof c);
    for (int t = 0; t < 2; t++) {
      bool done = false;
      for (int attempt = 0; attempt < 1000 && !done; attempt++) {
        TreeBuf B; std::memset(&B, 0, sizeof B);
        const int maxd = 2 + (int)r.below(5); const bool full = r.uni() < 0.5;
        rgrow(r, 0, maxd, full, B);
        if (!B.over && B.len >= 1 && B.len <= evo::MAXLEN) {
          c.len[t] = (uint8_t)B.len;
          std::memcpy(c.op[t], B.op, B.len); std::memcpy(c.c[t], B.c, B.len * sizeof(float));
          done = true;
        }
      }
      if (!done) { c.len[t] = 1; c.op[t][0] = evo::TVAR; c.c[t][0] = (float)r.gauss(); }
    }
    if (evo::cand_valid(c)) return;
  }
}

struct SeedSrc { FILE* f = nullptr; evo::Rng rng{20261008ull}; long file_read = 0, file_valid = 0, random_made = 0; bool file_done = false; };
static void next_batch(SeedSrc& S, int want, std::vector<evo::Cand>& out) {
  out.clear();
  while (S.f && !S.file_done && (int)out.size() < want) {
    evo::Cand c;
    if (fread(&c, sizeof c, 1, S.f) != 1) { S.file_done = true; break; }
    S.file_read++;
    if (evo::cand_valid(c)) { out.push_back(c); S.file_valid++; }
  }
  while ((int)out.size() < want) { evo::Cand c; random_cand(S.rng, c); out.push_back(c); S.random_made++; }
}


static void check_fit(const StepBatch& SB, int s, const Target& tg, int ncheck, uint64_t seed, long& n, long& miss, long& better, long& startmiss, std::vector<double>& ratios,
                      long& rmiss, long& rbetter, long& rstart, std::vector<double>& sdiff) {
  if (ncheck <= 0) return;
  std::vector<int> cand; for (int j : SB.ids[s]) if (SB.UQ.U[j].nth >= 1 && SB.F[s].is_g[j]) cand.push_back(j);
  std::mt19937_64 rng(seed); std::shuffle(cand.begin(), cand.end(), rng);
  if ((int)cand.size() > ncheck) cand.resize(ncheck);
  const double tau_base = 9.0 * ldexp(1.0, -24), floor_ = 1e-5 * tg.yy;
  for (int j : cand) {
    const Tree& tr = SB.UQ.U[j]; View V(tr.rec.data());
    double best = INFINITY; std::vector<double> bt;
    for (int q = 0; q < 16; q++) { double t2[32]; start_theta(tr.th, V.nth, q, t2); double r = eval_rss(V, t2, tg.y); if (r < best) { best = r; bt.assign(t2, t2 + V.nth); } }
    if (!std::isfinite(best)) continue;
    n++;
    if (fabs((double)SB.F[s].smin[j] - best) > 1e-3 * best + 1e-5 * tg.yy) startmiss++;
    if (fabs((double)SB.F[s].smin[j] - best) > 1e-3 * best) rstart++;
    sdiff.push_back(fabs((double)SB.F[s].smin[j] - best) / tg.yy);
    RefOut o = ref_gn(V, bt, tg.y, V.Q * tau_base, floor_);
    if (!std::isfinite(o.rssf)) continue;
    double g = SB.F[s].rfin[j];
    if (g > o.rssf * (1 + 1e-3) + 1e-5 * tg.yy) miss++;
    if (o.rssf > g * (1 + 1e-3) + 1e-5 * tg.yy) better++;
    if (g > o.rssf * (1 + 1e-3)) rmiss++;
    if (o.rssf > g * (1 + 1e-3)) rbetter++;
    ratios.push_back(g / std::max(o.rssf, 1e-30));
  }
}

#include <execinfo.h>
#include <mach-o/dyld.h>
extern "C" { typedef void (malloc_logger_t)(uint32_t, uintptr_t, uintptr_t, uintptr_t, uintptr_t, uint32_t); extern malloc_logger_t* malloc_logger; }
struct H900Rec { int st; uint32_t type; size_t size; int mainthr; void* fr[10]; int nfr; };
static H900Rec* g_h900_ring = nullptr; static std::atomic<size_t> g_h900_n{0}; static size_t g_h900_minb = (size_t)1 << 20;
static void h900_logger(uint32_t type, uintptr_t a1, uintptr_t a2, uintptr_t a3, uintptr_t res, uint32_t) {
  if (!(type & 2)) return;
  const size_t sz = (type & 4) ? (size_t)a3 : (size_t)a2;
  if (sz < g_h900_minb) return;
  const size_t i = g_h900_n.fetch_add(1, std::memory_order_relaxed); if (i >= 8192) return;
  H900Rec& r = g_h900_ring[i]; r.st = g_h900_st; r.type = type; r.size = sz; r.mainthr = pthread_main_np(); r.nfr = backtrace(r.fr, 10); (void)res;
}
static void h900_dump() {
  malloc_logger = nullptr; const size_t n = std::min<size_t>(g_h900_n.load(), 8192); const intptr_t slide = _dyld_get_image_vmaddr_slide(0);
  for (size_t i = 0; i < n; i++) { const H900Rec& r = g_h900_ring[i]; fprintf(stderr, "H900 host st=%d %s %.3f MB thr=%s @", r.st, (r.type & 4) ? "realloc" : "alloc", (double)r.size / 1048576.0, r.mainthr ? "main" : "other");
    for (int k = 2; k < r.nfr; k++) fprintf(stderr, " 0x%lx", (unsigned long)((intptr_t)r.fr[k] - slide)); fprintf(stderr, "\n"); }
  fprintf(stderr, "H900 host: %zu allocations >= %zu B\n", n, g_h900_minb);
}
static void h900_install() {
  g_h900_ring = (H900Rec*)calloc(8192, sizeof(H900Rec)); if (getenv("H900_MINB")) g_h900_minb = (size_t)atol(getenv("H900_MINB"));
  void* w[4]; backtrace(w, 4); atexit(h900_dump); malloc_logger = h900_logger;
}
volatile long g_thermal_state = -1;
static void* g_h917_pool = nullptr;
static std::chrono::steady_clock::time_point g_ui_t_start;
static double ui_ms_since_start() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - g_ui_t_start).count(); }
int engine_main(int argc, char** argv) { @autoreleasepool {
  if (getenv("MEM_DIAG")) { memdiag::init(); static id<MTLDevice> mdev = MTLCreateSystemDefaultDevice(); memdiag::start_sampler([]() -> double { return mb_of((double)[mdev currentAllocatedSize]); }); }
  g_ui_t_start = std::chrono::steady_clock::now();
  g_h917_pool = nullptr;
  printf("STATIC_FLAGS freet_async=%d h854=%d h257=%d h298q=%d h242d_fix=%.3f\n", (int)g_freet_async, g_h854, g_h257, g_h298q, g_h242d_fix);
  if (g_h900 & 1) h900_install();
  g_h870 = getenv("H870") ? atoi(getenv("H870")) : 0;
  if (argc < 6) { fprintf(stderr, "usage: %s <seeds.bin|random> kernel_dir outdir prims.txt label [steps=25] [threads=4] [check=200] [rseed=20261008]\n", argv[0]); return 2; }
  const std::string seedarg = argv[1], kdir = argv[2], outdir = argv[3], label = argv[5];
  g_psKdir = kdir;
  const int steps = argc > 6 ? atoi(argv[6]) : 25, nthr = argc > 7 ? atoi(argv[7]) : 4, ncheck = argc > 8 ? atoi(argv[8]) : 200;
  const uint64_t rseed = argc > 9 ? (uint64_t)strtoull(argv[9], nullptr, 10) : 20261008ull;
  if (const char* ev = getenv("S12_CHUNKS")) g_chunks = std::max(1, atoi(ev));
  printf("STATIC_FLAGS chunks=%d\n", g_chunks);
  if (const char* ev = getenv("H403")) evo::g_pool_h403 = atoi(ev);
  if (const char* ev = getenv("H744")) evo::g_pool_h744 = atoi(ev);
  if (const char* ev = getenv("H744S")) evo::g_pool_h744s = atof(ev);
  printf("H744: %d (spin-then-park pool, %.0f us window; needs H403)\n", evo::g_pool_h744, evo::g_pool_h744s);
  printf("H873: %d (persistent worker set; bit 1 h250a/h240b inner threads, bit 2 fw_submit)\n", g_h873);
  if (const char* ev = getenv("S14_STAGE8")) g_stage8_only = atoi(ev);
  if (const char* ev = getenv("S14_SRC")) g_src14 = atoi(ev);
  if (const char* ev = getenv("S14_CHK_BAND")) g_chk_band = atoi(ev);
  if (const char* ev = getenv("S14_DUMP")) g_dump14 = ev;
  if (const char* ev = getenv("S14_TIME")) g_time14 = atoi(ev);
  if (const char* ev = getenv("S14_CHK_FIN")) g_chk_fin = atoi(ev);
  fitx::reset(); rec_set_pin(); g_fit_cands = nullptr;
  g_sf_axis = 0; uib::g_xf_on = 0; g_fr_score96 = nullptr; g_sicis_en = nullptr; g_cf_on = false; g_cf_buf = nil; g_pp_seen.clear(); g_pp_hit = g_pp_all = 0;
  for (auto& pp : PP) for (int s2 = 0; s2 < 2; s2++) pp.S[s2].yb = nil;
  g_wb_par = nullptr; g_kgp_par = nullptr; g_h301_after_front = nullptr; g_ovl_hook = nullptr; g_kfill = nullptr; g_kfin = nullptr; g_kskey.clear(); g_kvalid.clear(); g_kcs_prev.clear(); g_h826_head = nullptr; g_h826_rest = nullptr;
  delete g_h415_spare; g_h415_spare = nullptr; g_pool2 = nullptr; g_h250_pending = nullptr;

  g_h301r_t.clear(); g_h301r_e.clear(); g_h301r_n.store(0); g_stp.clear(); g_st_chunk = 0; g_W2.clear();
  g_h670_res.clear(); g_h670_have = false; g_h670_cs.clear(); g_loc16.clear(); g_h91s.clear(); g_h774_cstart.clear(); g_h320_link.clear();
  g_fdq[0] = FreetDeferred(); g_fdq[1] = FreetDeferred(); g_fdi = 0; g_fr = FrStats(); g_frs = FrStats(); g_cfn = 0; g_cfdead = 0; g_cfmaxabs = 0; g_cfmaxrel = 0;
  g_h896p = nullptr; g_fit_span = nullptr; g_h53p = &g_h53sb[0]; g_h53cur = &g_h53sb[0]; g_h899_fin = nullptr;
  if (const char* ev = getenv("S12_FREET")) g_freet = atoi(ev);
  if (g_freet) {
    if (const char* ev = getenv("S12_FREET_MAXIT")) g_fr_maxit = std::max(1, atoi(ev));
    if (const char* ev = getenv("S12_FREET_DUMP")) g_fr_dump = ev;
    if (const char* ev = getenv("S12_FREET_DUMPAT")) g_fr_dumpat = atoi(ev);
    if (const char* ev = getenv("S12_FREET_DUMPN")) g_fr_dumpn = atoi(ev);
    if (const char* ev = getenv("S12_FREET_DUMPK")) g_fr_dumpk = std::max(1, atoi(ev));
    if (const char* ev = getenv("S12_FREET_VERBOSE")) g_fr_verbose = atoi(ev);
    if (g_h250a || g_h240b || g_h770 || g_h99 || g_h242 || getenv("H771") || getenv("H91") || getenv("H731") || getenv("H898") || getenv("H730") || getenv("H860")) { fprintf(stderr, "S12_FREET refuses H250A/H240B/H770/H771/H99/H242/H91/H731/H898/H730/H860 (fit paths replaced)\n"); return 1; }
    fjn::make_warp(g_fr_warp);
    printf("H950: S12_FREET=1 joint free-t pair fit (V1 LM, W=8 warp, maxit %d); GN itmax 0 (scan theta = fixed-t start); post16 replaced\n", g_fr_maxit);
  }
  if (const char* ev = getenv("S16_POST")) g_post16 = atoi(ev);
  if (const char* ev = getenv("H53")) g_h53 = atoi(ev);
  if (const char* ev = getenv("H92")) g_h92 = atoi(ev);
  if ((g_h890 & 1) && (g_h670 || g_h673 || g_h675 || g_h774 == 2)) { fprintf(stderr, "H890 bit 1 refused: H670/H673/H675/H774 carry child Metas across the next NF\n"); exit(1); }
  printf("H92: %d (bit 1 %s select, bit 2 %s dedup + unique ids)\n", g_h92, (g_h92 & 1) ? "move-mode (no record copies)" : "copy-mode", (g_h92 & 2) ? "sharded" : "serial");
  printf("H53: %s\n", g_h53 ? "on (post16_h53 + select_h53; I18 default)" : "off (I16 post16 + select_h33; opt-out)");
  if (const char* ev = getenv("H98")) g_h98 = atoi(ev);
  printf("H98: %d (%s)\n", g_h98, (g_h98 & 1) ? "record hash computed once, passed to prefetch and compile" : "off");
  if (const char* ev = getenv("H84")) g_h84 = atoi(ev);
  if (const char* ev = getenv("H91")) g_h91 = atoi(ev);
  if (!(g_h84 >= 2 && g_h53 && g_post16 == 2)) g_h91 = 0;
  printf("H91: %d (%s)\n", g_h91, g_h91 ? "per-chunk grid+band+polish encoded per chunk, waited at the end" : "off");
  if (const char* ev = getenv("H129")) g_h129_cap_mb = atof(ev);
  if (const char* ev = getenv("H310_CAP_MB")) { g_h129_cap_mb = atof(ev); printf("H310: memo cap override %.0f MB (config only; policy unchanged)\n", g_h129_cap_mb); }
  if (const char* ev = getenv("H912")) if (atoi(ev) != 0 && !getenv("H129") && !getenv("H310_CAP_MB")) {
    uint64_t ms = 0; size_t msz = sizeof ms; sysctlbyname("hw.memsize", &ms, &msz, nullptr, 0);
    double mem_mb = (double)ms / 1048576.0; if (const char* mv = getenv("H912_MEM_MB")) mem_mb = atof(mv);
    const double frac = 1.0 / 32.0, lo = 64.0, hi = 768.0;
    if (mem_mb > 0) { g_h129_cap_mb = std::min(hi, std::max(lo, frac * mem_mb));
      printf("H912: memo cap %.0f MB = clamp(%.5f x %.0f MB memsize, %.0f, %.0f)\n", g_h129_cap_mb, frac, mem_mb, lo, hi); }
    else printf("H912: hw.memsize unreadable, cap stays %.0f MB\n", g_h129_cap_mb);
  }
  {

    const int mc = getenv("MEM_CAP") ? atoi(getenv("MEM_CAP")) : 1;
    if (mc && g_h129_cap_mb > 0) { uint64_t ms = 0; size_t l = sizeof ms; double mem_mb = (sysctlbyname("hw.memsize", &ms, &l, nullptr, 0) == 0) ? (double)ms / 1048576.0 : 0.0;
      if (const char* mv = getenv("H912_MEM_MB")) mem_mb = atof(mv);
      if (mem_mb > 0) { const double lim = std::max(64.0, mem_mb / 64.0); if (g_h129_cap_mb > lim) { printf("MEM_CAP: memo cap %.0f MB -> %.0f MB (physical memory %.0f MB / 64, floor 64)\n", g_h129_cap_mb, lim, mem_mb); g_h129_cap_mb = lim; } } }
  }
  h156::g_on = true;
  if (const char* ev = getenv("H156")) h156::g_on = atoi(ev) != 0;
  if (const char* ev = getenv("H117")) g_h117 = atoi(ev);
  if (const char* ev = getenv("H131")) g_h131 = atoi(ev);
  if (const char* ev = getenv("H335")) g_h335 = atoi(ev);
  if (const char* ev = getenv("H333")) g_h333 = atoi(ev);
  if (g_h333 && !g_h335) { fprintf(stderr, "H333 needs H335=1 (the top-4 kernel writes the handoff)\n"); return 1; }
  if (const char* ev = getenv("H139")) g_h139 = atoi(ev);
  printf("H131: %d (bit 1 %s, bit 2 %s, bit 4 %s)\n", g_h131, (g_h131 & 1) ? "persistent post buffers" : "per-call post buffers",
         (g_h131 & 2) ? "parallel glue install" : "serial glue install", (g_h131 & 4) ? "dump + sets on a helper thread" : "serial dump + sets");
  g_h117_bit4 = (g_h117 & 4) ? 1 : 0;
  if (const char* ev = getenv("H128")) g_h128 = atoi(ev);
  printf("H117: %d (bit 1 %s, bit 2 %s, bit 4 %s, bit 8 %s)\n", g_h117, (g_h117 & 1) ? "tail CPU per chunk" : "tail CPU at tail",
         (g_h117 & 2) ? "one command buffer" : "command buffer per stage", (g_h117 & 4) ? "dense key ids" : "table", (g_h117 & 8) ? "popset on helper thread" : "serial");
  if (const char* ev = getenv("H93")) { g_h93 = atoi(ev) != 0; g_h93_print = true; }
  if (const char* ev = getenv("H215")) g_h215 = atoi(ev);
  if (const char* ev = getenv("H916")) g_h916 = atoi(ev);
  if (const char* ev = getenv("H217")) g_h217 = atoi(ev);
  if (const char* ev = getenv("SCANGN_OVL")) { if (atoi(ev) != 0) g_h897 |= 1; }
  if (g_h217 & ~1) { fprintf(stderr, "H217 accepts 0 or 1 (got %d)\n", g_h217); return 2; }
  if (g_h217 && (!g_h93 || g_post16 != 2 || g_h84 < 2 || g_h53 != 1)) { fprintf(stderr, "H217=1 needs H93 on, S16_POST=2, H84>=2, H53=1 (the default)\n"); return 2; }
  if (g_h826) {
    if (g_h826 != 1) { fprintf(stderr, "H826 accepts 0 or 1 (stage 2 is not built; got %d)\n", g_h826); return 2; }
    if (g_h217) { fprintf(stderr, "H826 refused with H217=1: front(c) runs on a producer thread, the rest-breed hook would race with it\n"); return 2; }
    if (!(g_h262 & 1)) { fprintf(stderr, "H826 needs H262 bit 1 (children read in place through chp = &XC(S,cs[i]).g, no ch[] copy)\n"); return 2; }
    if (g_h257 & 8) { fprintf(stderr, "H826 refused with H257 bit 8 (timing mode: the ch[] copy is inside T.breed)\n"); return 2; }
    printf("H826=1: breed demes of chunk 0 before front(0), the rest after S(0) is committed and before front(1) (%s breed)\n", g_h777 ? "H777" : "H92");
  }
  if (g_h217) printf("H217=1: producer-ahead front(c+1) on a dedicated thread (own pool, own workers, bounded to two parity slots)\n");
  if (g_h215 & 2) { fprintf(stderr, "H215 bit 2 (producer-ahead pipeline) is not built in loop_v11_i22_h215.mm; refusing H215=%d\n", g_h215); return 2; }
  if (g_h215) printf("H215=%d: bit 1 listed-only scan/GN encode (H93 %s)\n", g_h215, g_h93 ? "on" : "off");
  if (g_h916) {
    const char* why = nullptr;
    if (g_h916 != 1) why = "accepts 0 or 1";
    else if (!g_h93) why = "needs H93 on (scan_begin_p / gn_begin_p)";
    else if (!(g_h215 & 1)) why = "needs H215 bit 1 (the sorted listed set X.ls defines the local ids)";
    else if (g_h217) why = "refused with H217 (producer-ahead slots)";
    else if (g_h139 == 1) why = "refused with H139=1 (h139o compare is indexed by global id)";
    else if (g_h250b) why = "refused with H250B (CPU scan share reads the tables by global id)";
    else if (g_h333) why = "refused with H333 (hoff[tree])";
    else if (g_h335 == 2) why = "refused with H335=2 (check path reads out/h139o by global id)";
    else if (g_h770 || g_h771) why = "refused with H770/H771 (hit rows and ST[j] by global id)";
    if (why) { fprintf(stderr, "H916: %s\n", why); return 2; }
    printf("H916=1: chunk-local tree ids (per-tree GPU tables sized by the chunk's listed set)\n");
  }
  if (g_h263) printf("H263: mask %d (1 persistent scratch, 2 scan sort scattered into bl, 4 encoder binding cache)\n", g_h263);
  if (g_h257) printf("H257: bits %d (1 tail encode-ahead, 2 head taper, 4 build_unique on the pool, 8 timing honesty) H257F %.3f\n", g_h257, g_h257f);
  printf("H93: %s\n", g_h93 ? "on (scan_begin_p / gn_begin_p, parallel CPU encode)" : "off (serial scan_begin / gn_begin; fine timers only)");
  printf("H644: %d (pre-decoded token records)\n", (int)g644_on());
  printf("H660: %d D %s NS %.0f CPU %.0f\n", g_h660, g_h660d ? g_h660d : "-", g_h660ns, g_h660cpu);
  printf("H680: %d (D9 kill probe; CPU %.0f core-ms, GPU %.0f ms per step)\n", g_h680, g_h680cpu, g_h680g);
  printf("H731: %d (grid + band per chunk)\n", g_h731);
  printf("H770: %d (raw %d) VERIFY %d (needs H242)\n", g_h770, g_h770_raw, g_h770v);
  printf("H779: %d (bins heaviest first)\n", g_h779);
  printf("H777: %d (breed reads parent genomes by pointer)\n", g_h777);
  printf("H774: %d (phased select; 2 = verify vs select_h92)\n", g_h774);
  printf("H802: %d (scan pass 1 Q>=6 at 64 threads, G6/LN4)  H803: %d (scan pass 2 Q>=6 at 32 threads, G6/LN8)\n", g_h802, g_h803);
  printf("H790: %d (NaN-safe unary tangent, CPU + kbuild_h790)\n", g_h790);
  if (g_h680 == 2) h860_load();
  if (g_h731 && (g_h670 || !(g_h131 & 1))) { fprintf(stderr, "H731 needs H131 bit 1 and no H670\n"); return 1; }
  printf("H670: %d (one-step-lag polish)\n", g_h670);
  printf("H673: %d (revival of dropped polished children; needs H670)\n", g_h673);
  printf("H675: %d (lagged polish over the whole next step, %d utility threads; needs H670)\n", g_h675, g_h675t);
  if (g_h670 && !(g_h92 & 1)) { fprintf(stderr, "H670 needs move-mode select (H92 bit 1)\n"); exit(1); }
  printf("H657: (MGS in place) | H642: %d (dual exec p-outer) | H643: %d (MGS/dual setup in registers)\n", (int)g641_on(), (int)g642_on(), g643_bits());
  printf("H651: %d (lazy memo join) | H652: %d (best/check on a helper)\n", g_h651, g_h652);
  printf("H891: %d (delay injection; CPU %.0f us/batch-of-4 site %d, GPU %.0f us, head %.0f us, tail %.0f us)\n", g_h891, g_h891_cpu_us, g_h891_site, g_h891_gpu_us, g_h891_head_us, g_h891_tail_us);
  printf("H650: %s (M4 step-sensitivity probe; ms %.1f, fit ns/item %.0f)\n", g_h650d ? g_h650d : "off", g_h650ms, g_h650ns);
  printf("H421: %d (value-only unary ops in the CPU fit exec when C == 1)\n", g421_on() ? 1 : 0);
  printf("H419: %d (res[0]-only install, select sum shortcut)\n", g_h419);
  printf("H418: %d (per-step key sets kept across steps)\n", g_h418);
  if (g_h874) printf("H874: %d\n", g_h874);
  printf("H417: %d (CPU fit exec without operand/result staging)\n", g417_on() ? 1 : 0);
  printf("H313A: mask %d (A2 %d, A3 %d, A1 %d) defer %d X %d\n", g_h313a_mask, (g_h313a_mask & 1) != 0, (g_h313a_mask & 2) != 0, (g_h313a_mask & 4) != 0, (int)g_h313_defer, g_h313x);
  printf("H84: %d (%s)\n", g_h84, g_h84 >= 2 ? "parallel post CPU loops + per-chunk prep (S16_POST=2, H53)" : g_h84 == 1 ? "parallel post CPU loops" : "off");
  printf("H94: %d (%s)\n", h94_mode(), h94_mode() >= 1 ? "GN and polish bins at GPU-sized TG counts (claim loop, several trees per SIMD)" : "off");
  if (const char* ev = getenv("H99")) g_h99 = atoi(ev);
  if (const char* ev = getenv("H413")) g_h413 = atoi(ev);
  if (const char* ev = getenv("H242")) g_h242 = atoi(ev);
  if (g_h242) {
    g_h99 = 0;
    int l2 = 20;
    if (const char* ev = getenv("H242_LOG2")) l2 = atoi(ev);
    const char* vv = getenv("H242_VERIFY");
    G242.init((size_t)1 << l2, vv != nullptr && atoi(vv) != 0);
    const char* sy = getenv("H242_SYNC");
    G242.sync = sy != nullptr && atoi(sy) != 0;
    if (g_h770) {
      if (!g_h93) { fprintf(stderr, "H770 needs H93=1 (gn_begin_p)\n"); return 1; }
      if (g_h250b) { fprintf(stderr, "H770 does not support H250B (CPU scan share)\n"); return 1; }
      printf("H770: raw %d -> bits %d (1 side-keyed memo, 2 complete rows: CPU-fit a[0..Q), 4 virtual assignment; 1 and 4 coupled) | VERIFY %d | generation rotation (bit 5) not built: the table clears at 3/4 full as H242, clears are printed per step\n", g_h770_raw, g_h770, g_h770v);
    }
    printf("H242 memo: ON | mode %d | puts %s | slots %zu (%.0f MB, 128 B line each) | row cap %u floats | key = h99_key(record hash, t), 64-bit | verify %d | clear at 3/4 full\n",
           g_h242, G242.sync ? "sync (main thread)" : "async (one worker, barrier before each lookup)", G242.cap, (double)G242.cap * 128.0 / 1048576.0, (unsigned)H242_ROW, G242.verify ? 1 : 0);
  }
  if (g_h771) {
    if (g_h242) { fprintf(stderr, "H771 does not use the H242 memo: unset H242 (H413 is then a no-op)\n"); return 1; }
    if (g_h99) { fprintf(stderr, "H771 is refused with H99\n"); return 1; }
    if (!g_h93) { fprintf(stderr, "H771 needs H93=1 (gn_begin_p)\n"); return 1; }
    if (g_h250b) { fprintf(stderr, "H771 does not support H250B (CPU scan share)\n"); return 1; }
    int l2 = 4; while (((size_t)1 << l2) < (size_t)evo::NPOP * 3 / 2) l2++;
    if (const char* ev = getenv("H771_LOG2")) l2 = std::max(8, atoi(ev));
    const int ng = getenv("H771G") ? atoi(getenv("H771G")) : 2;
    h771_init((size_t)evo::NPOP, ng, l2);
    printf("H771: ON (implies H770 bits %d) | VERIFY %d | tables: %d generations x 2 sides x %zu slots x 128 B = %.1f MB, %d sub-tables each, no clear (stamped), flat (no lock, no atomics) | key = h99_key(NF record hash, t) | reads the %d previous generation(s), carries full hits forward | H242/H413 not used\n",
           g_h770, g_h770v, G771.ng, G771.subcap * H771_NSUB, (double)G771.ng * 2.0 * (double)(G771.subcap * H771_NSUB) * 128.0 / 1048576.0, H771_NSUB, G771.ng - 1);
  }
  { int le = 21, la = 26;
    if (const char* ev = getenv("H99_LOG2")) le = atoi(ev);
    if (const char* ev = getenv("H99_ARENA_LOG2")) la = atoi(ev);
    if (g_h99) G99.init((size_t)1 << le, (size_t)1 << la);
    printf("H99 memo: %s | entries %zu | key arena %zu words (%.0f MB) | results %.0f MB | key = full NF record (hash + exact compare), FIFO\n",
           g_h99 ? "ON" : "OFF", (size_t)(g_h99 ? G99.cap : 0), (size_t)(g_h99 ? G99.arena_w : 0), g_h99 ? G99.arena_w * 4.0 / 1048576.0 : 0.0,
           g_h99 ? G99.bytes() / 1048576.0 - G99.arena_w * 4.0 / 1048576.0 : 0.0); }
  if (![[NSFileManager defaultManager] createDirectoryAtPath:[NSString stringWithUTF8String:outdir.c_str()] withIntermediateDirectories:YES attributes:nil error:nil]) return 1;


  std::vector<Stroke> strokes; if (!load_strokes(argv[4], strokes)) { fprintf(stderr, "cannot read %s\n", argv[4]); return 1; }
  const Stroke* sk = nullptr; for (auto& s : strokes) if (s.label == label) { sk = &s; break; }
  if (!sk) { fprintf(stderr, "no stroke with label %s\n", label.c_str()); return 1; }
  std::vector<double> X, Y; to_npt(sk->xs, sk->ys, X, Y);
  size_t n_raw = sk->xs.size();
  if (const char* sin = getenv("STROKE_IN")) {
    std::vector<float> rx, ry;
    if (!stroke_in::load_raw(sin, rx, ry) || !stroke_in::resample_normalize(rx.data(), ry.data(), (int)rx.size(), NPT, X, Y)) {
      fprintf(stderr, "STROKE_IN: cannot load or normalize %s\n", sin); return 1; }
    n_raw = rx.size();
    printf("STROKE_IN: %s: %zu raw pts -> %d pts (arc-length resample + bbox-centre/span normalize, selection_core.hpp:95)\n", sin, rx.size(), NPT);
  }
  std::vector<double> XG, YG;

  int sf_axis = 0;
  if (getenv("SINGLE_FAST") && atoi(getenv("SINGLE_FAST")) != 0) {
    SingleFastDecision sfd = single_fast_detect(X, Y); sf_axis = sfd.axis;
    printf("SINGLE_FAST: axis %d (%s) margin %.3f eps %.4f backtrack x %.3f y %.3f closed %d\n", sf_axis,
           sf_axis == 1 ? "y=f(x), t:=x" : sf_axis == 2 ? "x=f(y), t:=y" : "not a function graph -> free t", sfd.margin, sfd.eps, sfd.backtrack_x, sfd.backtrack_y, (int)sfd.closed);
    if (sf_axis && g_yfx_flag && NPT == 48) {
      const std::vector<double>& A0 = sf_axis == 1 ? X : Y; std::vector<double> ax(A0), bx(sf_axis == 1 ? Y : X), tx(X), ty(Y);
      if (A0.front() > A0.back()) { std::reverse(ax.begin(), ax.end()); std::reverse(bx.begin(), bx.end()); std::reverse(tx.begin(), tx.end()); std::reverse(ty.begin(), ty.end()); }
      const double a0 = *std::min_element(A0.begin(), A0.end()), d = *std::max_element(A0.begin(), A0.end()) - a0;
      if (d > 0) { for (int p = 0; p < NPT; p++) { const double sd = (ax[p] - a0) / d; g_yfx_t[p] = sd; g_yfx_y[sf_axis == 1 ? 0 : 1][p] = sd; g_yfx_y[sf_axis == 1 ? 1 : 0][p] = bx[p]; }
        for (int s = 0; s < 2; s++) for (int p = 0; p < NPT; p++) g_yfx_f[s][p] = (float)g_yfx_y[s][p];
        for (int p = 0; p < NPT; p++) { const int a = p > 0 ? p - 1 : p, c = p < NPT - 1 ? p + 1 : p; const double ux = tx[c] - tx[a], uy = ty[c] - ty[a], h = std::hypot(ux, uy);
          const double nx = h > 0 ? -uy / h : 0.0, ny = h > 0 ? ux / h : 1.0;
          g_yfx_nt[p] = (float)nx; g_yfx_nt[NPT + p] = (float)ny;
          g_yfx_n[p] = sf_axis == 1 ? d * nx : nx; g_yfx_n[NPT + p] = sf_axis == 1 ? ny : d * ny; }
        g_yfx_txy[0] = tx; g_yfx_txy[1] = ty; g_yfx = true; g_yfx_epoch++; geop::g_sc[0] = sf_axis == 1 ? d : 1.0; geop::g_sc[1] = sf_axis == 1 ? 1.0 : d; } }
    if (g_yfx) { double mo = 0, so = 0; for (int p = 0; p < NPT; p++) { const double o = std::fabs(47.0 * g_yfx_t[p] - p); mo = std::max(mo, o); so += o; }
      printf("YFX_CELLS: on, driver span d %.4f, residuals at the 48 arc-length points (driver-grid offset max %.2f mean %.2f cells), energies from the same points\n", geop::g_sc[sf_axis == 1 ? 0 : 1], mo, so / NPT); }
    if (sf_axis) { std::vector<double> ox, oy; single_fast_resample(X, Y, sf_axis, NPT, ox, oy); X = ox; Y = oy; g_freet = 0; g_sf_axis = sf_axis; if (g_kesel && !getenv("H774") && g_ovl != 2) g_h774 = 1;
      XG = X; YG = Y; { std::vector<double>& D = sf_axis == 1 ? XG : YG; for (auto& v : D) v = g_sf_a0 + (g_sf_a1 - g_sf_a0) * v; } }
  }
  if (XG.empty()) { XG = X; YG = Y; }
  Sketch::CisPerCellEnergyInput sicis_energies;
  if (sicis::enabled()) {


    static const int sgrid = getenv("SICIS_GRID") ? atoi(getenv("SICIS_GRID")) : 1;
    std::vector<float> fx(sgrid ? XG.begin() : X.begin(), sgrid ? XG.end() : X.end()), fy(sgrid ? YG.begin() : Y.begin(), sgrid ? YG.end() : Y.end());
    if (g_yfx) { fx.assign(g_yfx_txy[0].begin(), g_yfx_txy[0].end()); fy.assign(g_yfx_txy[1].begin(), g_yfx_txy[1].end()); }
    bool const ok = sicis::strokeEnergies(fx.data(), fy.data(), (int)fx.size(), sicis_energies);
    printf("SICIS: ZI5 energies %s cells=%zu members=%zu\n", ok ? "ok" : "UNAVAILABLE (scoring off)",
           sicis_energies.perCellEnergy.size(), sicis_energies.memberEnergy.size());
    if (ok) g_sicis_en = &sicis_energies;
    if (ok && (getenv("FJN_RES48_HOST") || getenv("FJN_CISFUSE_HOST"))) {
      g_fr_score96 = [](const float* r) -> double {
        float c[NPT];
        for (int p = 0; p < NPT; p++) { const double zx = r[p], zy = r[NPT + p], h = std::hypot(zx, zy); c[p] = g_nres ? (float)(zx * (double)g_nres_n[p] + zy * (double)g_nres_n[NPT + p]) : (float)((zx + zy < 0) ? -h : h); }
        return sicis::score48(c, *g_sicis_en);
      };
      printf("SICIS: free-t objective = SI-CIS NLL on kernel res48 (FJN_RES48_HOST)\n");
    }
    if (ok && getenv("FJN_CISFUSE_HOST")) {
      const auto& me = sicis_energies.memberEnergy; bool good = !me.empty() && me.size() <= 8 && Sketch::cisEnsembleModeFromEnv() != 0;
      std::vector<float> cf(400, 0.0f); double K[8] = {0}, Kmin = 1e300;
      for (size_t m = 0; good && m < me.size(); m++) {
        if (me[m].size() != 48) { good = false; break; }
        for (size_t c = 0; c < 48; c++) { const double p = std::max(me[m][c], 0.0); if (!(p > 0.0) || !std::isfinite(p)) { good = false; break; }
          const double h = 0.5 / p; if (!std::isfinite((float)h)) { good = false; break; } cf[8 + m * 48 + c] = (float)h; K[m] += 0.5 * (std::log(2.0 * 3.14159265358979323846) + std::log(p)); }
        Kmin = std::min(Kmin, K[m]);
      }
      if (!good) { fprintf(stderr, "FJN_CISFUSE_HOST: ensemble energies unusable for the in-kernel score (need 1..8 members x 48 cells, p>0, SKETCH_CIS_ENSEMBLE!=0)\n"); return 1; }
      const bool joint = Sketch::cisJointCoarseFromEnv() != 0;
      if (joint) {
        const int blkN[7] = {24, 12, 6, 3, 1, 1, 1};
        Kmin = 1e300;
        for (size_t m = 0; m < me.size(); m++) {
          double Km = 0.0;
          for (int b = 0; b < 7; b++) { const double D = blkN[b]; Km += 0.5 * D * std::log(2.0 * 3.14159265358979323846) - std::lgamma(0.5 * D); }
          for (size_t c = 0; c < 48; c++) Km += 0.5 * std::log(std::max(me[m][c], 0.0));
          K[m] = Km; Kmin = std::min(Kmin, Km);
        }
      }
      const float Khi = (float)Kmin; cf[0] = (float)me.size(); cf[1] = g_cfchk > 0 ? 1.0f : 0.0f; cf[2] = Khi; cf[3] = (float)(Kmin - (double)Khi);
      cf[4] = joint ? 1.0f : 0.0f;
      cf[5] = getenv("FJN_S1") ? 1.0f : 0.0f;
      for (size_t m = 0; m < me.size(); m++) cf[392 + m] = (float)(K[m] - Kmin);
      g_cf_data = cf; g_cf_on = true; g_cf_buf = nil;
      printf("SICIS: CISFUSE ON: score in the fit kernel (members %zu, check every %d pairs)\n", me.size(), g_cfchk);
    }
  }
  printf("stroke %s (%s): %zu pts -> %d pts (arc-length resample%s)\n", sk->id.c_str(), sk->label.c_str(), sk->xs.size(), NPT,
         sk->xs.size() == (size_t)NPT ? ", skipped: already NPT" : "");
  if (!g_freet) g_fr_score96 = nullptr;
  if (getenv("STROKE_IN")) {
    const double s = stroke_in::g_span, d = g_sf_a1 - g_sf_a0;
    uib::g_map[0] = s; uib::g_map[1] = stroke_in::g_cx; uib::g_map[2] = s; uib::g_map[3] = stroke_in::g_cy;
    if (sf_axis == 1) { uib::g_map[0] = s * d; uib::g_map[1] = s * g_sf_a0 + stroke_in::g_cx; }
    if (sf_axis == 2) { uib::g_map[2] = s * d; uib::g_map[3] = s * g_sf_a0 + stroke_in::g_cy; }
  }


  uib::g_xf_on = 0; for (int i = 0; i < 4; i++) uib::g_smap[i] = uib::g_map[i];
  { static const int yfx = getenv("YFX_FRAME") ? atoi(getenv("YFX_FRAME")) : 1;
    const double d = g_sf_a1 - g_sf_a0;
    if (yfx && sf_axis && d > 0) {
      uib::g_xf_on = 1; uib::g_xf_drv = sf_axis == 1 ? 0 : 1; uib::g_xf_a = 1.0 / d; uib::g_xf_b = -g_sf_a0 / d;
      if (!getenv("STROKE_IN")) { uib::g_smap[0] = uib::g_smap[2] = 1; uib::g_smap[1] = uib::g_smap[3] = 0; if (sf_axis == 1) { uib::g_smap[0] = d; uib::g_smap[1] = g_sf_a0; } else { uib::g_smap[2] = d; uib::g_smap[3] = g_sf_a0; } }
      else { uib::g_map[0] = uib::g_map[2] = stroke_in::g_span; uib::g_map[1] = stroke_in::g_cx; uib::g_map[3] = stroke_in::g_cy; }
      printf("YFX_FRAME: %s equations in the drawing frame: t = %s over [%.6g, %.6g] (engine s = %.9g*t + %.9g), x_raw = %.6g*x + %.6g\n", sf_axis == 1 ? "y=f(x)" : "x=f(y)", sf_axis == 1 ? "x" : "y", g_sf_a0, g_sf_a1, uib::g_xf_a, uib::g_xf_b, uib::g_map[0], uib::g_map[1]);
    } }

  const int NP = evo::NPOP;
  SeedSrc src; src.rng = evo::Rng(rseed);
  const bool random_seeds = seedarg == "random";
  if (!random_seeds) { src.f = fopen(seedarg.c_str(), "rb"); if (!src.f) { perror("seeds"); return 1; } }
  printf("s12 loop v11 (I18: stages 8 + 10 + 9, K=%d chunks): P=%d children/step (%d demes x %d), steps %d, threads %d, check %d, seeds %s, rseed %llu\n",
         g_chunks, NP, evo::NDEME, evo::DECK, steps, nthr, ncheck, random_seeds ? "random (ramped half-and-half, len<=24)" : seedarg.c_str(), (unsigned long long)rseed);
  printf("I18: miss = f42::compile_fused (H42 SU order, fused encode), post mode S16_POST=%d (%s), select H53 (H33 if H53=0), breed H25 (breed_i12), NF caps x8\n",
         g_post16,
         g_post16 == 1 ? "post16 per chunk" : g_post16 == 2 ? "post16 once after all chunks" : "post14 serial (v9 reference)");

  Engine E;
  E.dev = MTLCreateSystemDefaultDevice(); E.cq = [E.dev newCommandQueue];
  if (g_h891 && g_h891_gpu_us > 0) {
    h680_init(E.dev);
    const uint32_t nit = (uint32_t)std::max(1.0, g_h891_gpu_us * 1e-3 / g_h680_ms_per_iter); double best = 1e30;
    for (int r = 0; r < 5; r++) { id<MTLCommandBuffer> cb = [g_h680q commandBuffer]; id<MTLComputeCommandEncoder> en = tl_enc(cb);
      [en setComputePipelineState:g_h680ps]; [en setBuffer:g_h680o offset:0 atIndex:0]; [en setBytes:&nit length:4 atIndex:1];
      [en dispatchThreadgroups:MTLSizeMake(256, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)]; [en endEncoding]; [cb commit]; [cb waitUntilCompleted];
      best = std::min(best, (cb.GPUEndTime - cb.GPUStartTime) * 1e6); }
    printf("H891: GPU spin %u iters, requested %.0f us, standalone CB GPU time %.1f us\n", nit, g_h891_gpu_us, best);
  }
  if (h178_on()) h178_open(E.dev);
  printf("I21 defaults: H141A %d H141B %d H150P %d H139 %d (env overrides: H141A=0 H141B=0 H150P=0 H139=0 restore i20)\n", h141a_on(), h141b_on(), h150p_on(), g_h139);
  printf("I22 defaults: H156 %d H182 %d H178 %d H178P %d H139 %d H141A %d H141B %d H150P %d (env off: H156=0 H182=0 H178=0 H178P=0)\n", h156::g_on ? 1 : 0, h182_on(), h178_on(), h178_pre(), g_h139, h141a_on(), h141b_on(), h150p_on());
  { NSError* e = nil; NSString* ps = [NSString stringWithUTF8String:(kdir + "/scan.metallib").c_str()];
    E.Ls = [E.dev newLibraryWithURL:[NSURL fileURLWithPath:ps] error:&e]; if (!E.Ls) { fprintf(stderr, "scan.metallib: %s\n", e.description.UTF8String); return 1; } }
  if (g_h335) { NSError* e = nil; g_h335L = [E.dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(kdir + (g_h333 ? "/kbuild_h335/scan_h333.metallib" : ((g_h802 || g_h803) ? (g_h802 && g_h803 ? "/kbuild_h802_h803/scan_h335.metallib" : g_h802 ? "/kbuild_h802/scan_h335.metallib" : "/kbuild_h803/scan_h335.metallib") : (getenv("H808F") && atoi(getenv("H808F")) >= 2) ? (atoi(getenv("H808F")) == 2 ? "/kbuild_h810_b1/scan_h335.metallib" : "/kbuild_h810_b1_bu/scan_h335.metallib") : (getenv("H808F") && atoi(getenv("H808F"))) ? (std::string("/kbuild_h808_f") + ((getenv("H805") && atoi(getenv("H805"))) ? "805" : "") + ((getenv("H806") && atoi(getenv("H806"))) ? ((getenv("H805") && atoi(getenv("H805"))) ? "_806" : "806") : "") + "/scan_h335.metallib").c_str() : getenv("H805") && atoi(getenv("H805")) && getenv("H806") && atoi(getenv("H806")) ? "/kbuild_h805_h806/scan_h335.metallib" : getenv("H805") && atoi(getenv("H805")) ? "/kbuild_h805/scan_h335.metallib" : getenv("H806") && atoi(getenv("H806")) ? "/kbuild_h806/scan_h335.metallib" : getenv("H622") && atoi(getenv("H622")) ? "/kbuild_h622/scan_h335.metallib" : getenv("H621") && atoi(getenv("H621")) ? "/kbuild_h621/scan_h335.metallib" : "/kbuild_h335/scan_h335.metallib"))).c_str()]] error:&e];
    if (!g_h335L) { fprintf(stderr, "H335: kbuild_h335/scan_h335.metallib: %s\n", e.description.UTF8String); return 1; }
    if ((g_h802 || g_h803) && ((getenv("H805") && atoi(getenv("H805"))) || (getenv("H806") && atoi(getenv("H806"))) || g_h333)) { fprintf(stderr, "H802/H803 need H335=1 and exclude H805/H806/H333\n"); return 1; }
    printf("H335=%d: 24-point scan ranking + exact top-4 (kbuild_h335/scan_h335.metallib)\n", g_h335);  if (getenv("H621") && atoi(getenv("H621"))) printf("H621: scan from kbuild_h621 (lane-parallel top-4)\n");  if (getenv("H622") && atoi(getenv("H622"))) printf("H622: scan from kbuild_h622 (H621 + all-bad MGS skip)\n");  if (getenv("H805") && atoi(getenv("H805"))) printf("H805: scan from kbuild_h805 (skip rank-dropped columns)\n"); if (getenv("H806") && atoi(getenv("H806"))) printf("H806: scan from kbuild_h806 (bad-tree early exit)\n"); }
  if (g_h139) { NSError* e = nil; g_h139L = [E.dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(kdir + "/scan_h139.metallib").c_str()]] error:&e];
    if (!g_h139L) { fprintf(stderr, "H139: scan_h139.metallib: %s\n", e.description.UTF8String); return 1; }
    printf("H139=%d: device-scratch scan kernels rcscand_* from scan_h139.metallib (over-limit = tgb > 32768 B)\n", g_h139); }

  auto load14 = [&](const char* base) -> id<MTLLibrary> {

    const bool h94k = h94_mode() >= 1 && (std::string(base) == "gn_i14" || std::string(base) == "pol_i14");
    const char* h331t = getenv("H331_TAUD");
    const bool h331k = !h94k && h331t && *h331t && (std::string(base) == "gn_i14" || std::string(base) == "pol_i14" || std::string(base) == "gnp_i14" || std::string(base) == "polp_i14");
    const bool h333k = !h94k && getenv("H333") && atoi(getenv("H333")) && std::string(base) == "gn_i14";
    const bool h332k = !h94k && g_h332 && std::string(base) == "gn_i14";
    const std::string dir = kdir + (h94k ? "/kbuild_h94/" : h332k ? ("/kbuild_h332" + std::string(h333k ? "_h333" : "") + (h331t && *h331t ? "_" + std::string(h331t) : std::string()) + "/") : h333k ? ("/kbuild_h333" + (h331t && *h331t ? "_" + std::string(h331t) : std::string()) + "/") : h331k ? (((g_h712 || g_h713) ? "/kbuild_h712b_" : g_h711 ? "/kbuild_h711b_" : g_h790 ? "/kbuild_h790_" : getenv("H710B") && atoi(getenv("H710B")) ? "/kbuild_h331b_" : "/kbuild_h331_") + std::string(h331t) + ((g_h712 && !g_h711) ? ((getenv("H710B") && atoi(getenv("H710B"))) ? "_n711" : "_n711_a") : "") + ((g_h713 && !g_h712) ? "_f713" : "") + "/") : "/kbuild_i14/") + base;
    id<MTLLibrary> L = nil;
    if (g_src14) {
      std::ifstream f(dir + ".metal"); if (!f) { fprintf(stderr, "missing %s.metal\n", dir.c_str()); return nil; }
      std::stringstream ss; ss << f.rdbuf(); NSError* e = nil; MTLCompileOptions* o = [MTLCompileOptions new];
      auto t0 = clk::now(); L = [E.dev newLibraryWithSource:[NSString stringWithUTF8String:ss.str().c_str()] options:o error:&e];
      printf("runtime source compile %s: %.0f ms\n", base, ms_since(t0));
      if (!L) { fprintf(stderr, "%s: %s\n", base, e.description.UTF8String); return nil; }
      return L; }
    NSError* e = nil; L = [E.dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(dir + ".metallib").c_str()]] error:&e];
    if (!L) { fprintf(stderr, "%s.metallib: %s\n", base, e.description.UTF8String); return nil; }
    return L; };
  if (g_h711) {
    const char* t711 = getenv("H331_TAUD");
    auto on711 = [](const char* n) { const char* v = getenv(n); return v && atoi(v); };
    if (!on711("H710B") || !t711 || !*t711) { fprintf(stderr, "H711 needs H710B=1 and H331_TAUD=T (kbuild_h711b_<T>/)\n"); return 1; }
    if (on711("H333") || on711("H332") || on711("H790") || on711("H94") || g_src14) { fprintf(stderr, "H711 is refused with H333/H332/H790/H94/S14_SRC\n"); return 1; }
    printf("H711=%d: fused mode-0 + first dual, gn/pol from kbuild_h711b_%s/, EBt at buffer 15\n", g_h711, t711); }
  if (g_h713 && !g_h711) { fprintf(stderr, "H713 needs H711=1 (kbuild_h712b_<T>/)\n"); return 1; }
  if (g_h712 && !g_h711 && (g_h790 || g_h713)) { fprintf(stderr, "H712 without H711 is refused with H790/H713\n"); return 1; }
  if (g_h712 || g_h713) printf("H712=%d H713=%d: gn/pol/gnp/polp from kbuild_h712b_%s%s%s/\n", g_h712, g_h713, getenv("H331_TAUD"), (g_h712 && !g_h711) ? ((getenv("H710B") && atoi(getenv("H710B"))) ? "_n711" : "_n711_a") : "", (g_h713 && !g_h712) ? "_f713" : "");
  E.Lg = load14("gn_i14"); E.Lp = load14("pol_i14"); E.Lgr = load14("grid_i14"); E.Lb = load14("band_i14");
  if (!E.Lg || !E.Lp || !E.Lgr || !E.Lb) return 1;
  if (g_freet) { if (const char* fl = getenv("FJN_LAY")) fjn::g_lay = atoi(fl); if (const char* fz = getenv("FJN_FUSE_HOST")) fjn::g_fuse = atoi(fz); NSError* fe = nil; g_fjnL = [E.dev newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:(kdir + "/fjn_h950.metallib").c_str()]] error:&fe];
    if (!g_fjnL) { fprintf(stderr, "fjn_h950.metallib: %s\n", fe.description.UTF8String); return 1; } }
  if (g_h713) { E.Lgp = load14("gnp_i14"); E.Lpp = load14("polp_i14"); if (!E.Lgp || !E.Lpp) return 1; }
  E.SH = new nf16::Shared(); E.SH->prefault((size_t)16 << 20); E.CS = new cc14::Shared(); if (!(g_h853 & 2)) E.CS->prefault((size_t)64 << 20);
  E.MM = new cmi6::Memo();
  struct PfJoin { std::thread t; ~PfJoin() { if (t.joinable()) t.join(); } } pf_thr;
  if (g_prefault) pf_thr.t = std::thread([&E] { h853_prefault(E); });
  else h853_prefault(E);
  if ((g_h874 & 2) && g_h415 && g_h301r && (g_h864 & 4)) h874_pretouch((size_t)NP);
  if (g_h874 & 16) g_h874_mf0 = h874_mf();
  WK(E).resize(std::max(1, nthr));
  for (Worker& w : WK(E)) { w.C = new Ctx(); w.H = new nf16::Handle(*E.SH, w.C); w.CH = new cmi6::Handle(*E.CS, *E.MM); }
  E.set_target(0, X); E.set_target(1, Y);
  if (g_nres) { nres_init(E.tg[0].y, E.tg[1].y); if (g_cf_on && g_cf_data.size() == 400) { g_cf_data.resize(496); memcpy(g_cf_data.data() + 400, g_nres_n, sizeof g_nres_n); g_cf_buf = nil; }
    printf("NORMAL_RES: SI-CIS residual = signed normal residual (stroke normal at the fit points)\n"); }
  if (getenv("PSO_PAR") && atoi(getenv("PSO_PAR")) != 0 && (!g_h178_arch || !g_h178_exists)) {
    const auto tpp = clk::now();
    std::vector<std::pair<id<MTLLibrary>, std::string>> L;
    if (g_freet) for (NSString* fnm in g_fjnL.functionNames) L.push_back({g_fjnL, std::string(fnm.UTF8String)});
    for (int q = 1; q <= 16; q++) { L.push_back({E.Ls, "rcscan_q" + std::to_string(q)}); L.push_back({E.Ls, "rcscan_s2_q" + std::to_string(q)}); }
    if (g_h139 && g_h139L) for (int q = 1; q <= 16; q++) { L.push_back({g_h139L, "rcscand_q" + std::to_string(q)}); L.push_back({g_h139L, "rcscand_s2_q" + std::to_string(q)}); }
    for (int q = 1; q <= 16; q++) for (int t = 1; t <= GN_TH_MAX; t++) {
      L.push_back({E.Lg, "gnw_q" + std::to_string(q) + "_t" + std::to_string(t)}); L.push_back({E.Lp, "pol_q" + std::to_string(q) + "_t" + std::to_string(t)}); }
    L.push_back({E.Lgr, "grid_h26"}); L.push_back({E.Lb, "post_gather_i14"}); L.push_back({E.Lb, "band_lane"}); L.push_back({E.Lb, "post_repair_i14"});
    Engine* EP = &E; auto* LP = &L;
    const long width = std::max(1, getenv("PSO_PAR_N") ? atoi(getenv("PSO_PAR_N")) : 3);
    dispatch_semaphore_t sem = dispatch_semaphore_create(width);
    std::atomic<long> done{0}; std::atomic<long>* DP = &done; const size_t NL = L.size();
    printf("PSO_PAR: start %zu pipelines, width %ld\n", NL, width); fflush(stdout);
    dispatch_apply(L.size(), dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t i) {
      dispatch_semaphore_wait(sem, DISPATCH_TIME_FOREVER);
      @autoreleasepool { EP->pipe((*LP)[i].first, (*LP)[i].second); }
      dispatch_semaphore_signal(sem);
      const long d = ++(*DP); if (d % 100 == 0 || d == (long)NL) { printf("PSO_PAR: %ld / %zu\n", d, NL); fflush(stdout); } });
    printf("PSO_PAR: %zu pipelines compiled concurrently in %.0f ms (created %ld absent %ld)\n", L.size(), ms_since(tpp), g_h178_created, g_h178_absent); fflush(stdout);
  }
  if (g_freet && !(getenv("S12_FREET_PRECREATE") && atoi(getenv("S12_FREET_PRECREATE")) == 0)) {
    const auto tpc = clk::now(); int npc = 0; for (NSString* fnm in g_fjnL.functionNames) { E.pipe(g_fjnL, std::string(fnm.UTF8String)); npc++; }
    printf("H950: pre-created %d fjn PSOs in %.1f ms (archive hits %ld misses %ld)\n", npc, ms_since(tpc), g_h178_hits, g_h178_miss); }
  if (g_freet) {
    double Lr = 0; for (int i = 1; i < NPT; i++) Lr += std::hypot(X[i] - X[i - 1], Y[i] - Y[i - 1]);
    const double gap = std::hypot(X[NPT - 1] - X[0], Y[NPT - 1] - Y[0]);


    static const bool cnr = getenv("CLOSED_NRAW") && atoi(getenv("CLOSED_NRAW")) != 0;
    g_fr_closed = gap <= Lr / std::sqrt((double)(cnr ? n_raw : sk->xs.size()));


    { static const bool cr = getenv("CLOSED_RULE") && atoi(getenv("CLOSED_RULE")) != 0; if (cr) g_fr_closed = single_fast_detect(X, Y).closed; }
    printf("H950: stroke arc length %.4f end gap %.4f raw points %zu -> %s\n", Lr, gap, cnr ? n_raw : sk->xs.size(), g_fr_closed ? "CLOSED (t0 free)" : "open (pinned ends)"); }
  printf("targets: x in [%.4f, %.4f], y in [%.4f, %.4f], yy_x %.4f yy_y %.4f\n",
         *std::min_element(X.begin(), X.end()), *std::max_element(X.begin(), X.end()),
         *std::min_element(Y.begin(), Y.end()), *std::max_element(Y.begin(), Y.end()), E.tg[0].yy, E.tg[1].yy);
  { auto w0 = clk::now();
    for (int q = 1; q <= 16; q++) { E.pipe(E.Ls, "rcscan_q" + std::to_string(q)); E.pipe(E.Ls, "rcscan_s2_q" + std::to_string(q)); }
    if (g_h139) for (int q = 1; q <= 16; q++) { h139_ps(E, 0, q); h139_ps(E, 1, q); }
    int npso = 0; for (int q = 1; q <= 16; q++) for (int t = 1; t <= GN_TH_MAX; t++) { if (E.pipe(E.Lg, "gnw_q" + std::to_string(q) + "_t" + std::to_string(t))) npso++; }
    if (g_h713) { int np713 = 0; for (int q = 1; q <= 16; q++) for (int t = 9; t <= GN_TH_MAX; t++) { if (E.pipe(E.Lgp, "gnp_q" + std::to_string(q) + "_t" + std::to_string(t))) np713++; if (E.pipe(E.Lpp, "polp_q" + std::to_string(q) + "_t" + std::to_string(t))) np713++; }
      printf("H713 pair pipelines: %d\n", np713); }
    printf("pipeline warm-up: %d GN + 32 scan pipelines in %.0f ms (one-time, excluded from step timings)\n", npso, ms_since(w0)); }
  if (h178_pre()) {
    auto w1 = clk::now(); long inuse = 0; const long a0 = g_h178_absent;
    auto add = [&](id<MTLLibrary> L, const std::string& nm) { if (E.pipe(L, nm)) inuse++; };
    for (int q = 1; q <= 16; q++) { add(E.Ls, "rcscan_q" + std::to_string(q)); add(E.Ls, "rcscan_s2_q" + std::to_string(q)); }
    if (g_h139 && g_h139L) for (int q = 1; q <= 16; q++) { add(g_h139L, "rcscand_q" + std::to_string(q)); add(g_h139L, "rcscand_s2_q" + std::to_string(q)); }
    for (int q = 1; q <= 16; q++) for (int t = 1; t <= GN_TH_MAX; t++) {
      add(E.Lg, "gnw_q" + std::to_string(q) + "_t" + std::to_string(t));
      add(E.Lp, "pol_q" + std::to_string(q) + "_t" + std::to_string(t)); }
    add(E.Lgr, "grid_h26"); add(E.Lb, "post_gather_i14"); add(E.Lb, "band_lane"); add(E.Lb, "post_repair_i14");
    printf("H178P=1 eager: %ld PSOs in use (absent names %ld) in %.1f ms | pipes created so far %ld | archive hits %ld misses %ld\n",
           inuse, g_h178_absent - a0, ms_since(w1), g_h178_created, g_h178_hits, g_h178_miss); }
  if (h178_on() || h178_pre()) printf("H178 startup: PSO table %zu | created %ld | absent names %ld | archive hits %ld misses %ld | PSO creation wall %.1f ms (total, all creations so far)\n",
                                      E.pso.size(), g_h178_created, g_h178_absent, g_h178_hits, g_h178_miss, g_h178_ms);
  if (h876_mask() & 4) printf("H876 mem init (before release): rss %.1f MB phys %.1f MB\n", rss_mb_now(), phys_fp_mb());
  if ((h876_mask() & 1) && g_h178_arch) {
    h178_flush(); g_h178_dirty = false;
    g_h178_arch = nil; g_h178_released = true;
    printf("H876: MTLBinaryArchive reference released\n"); }
  if (h876_mask() & 4) printf("H876 mem init end: rss %.1f MB phys %.1f MB\n", rss_mb_now(), phys_fp_mb());

  if (g_h298q & 2) pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
  printf("H299: %d (%s)\n", g_h299, g_h299 ? "per-chunk polish on the CPU after the chunk's grid + band (needs H91=1)" : "off");
  printf("H298Q: %d (bit1 fit threads UTILITY, bit2 main+pool USER_INTERACTIVE, bit4 tail polish UTILITY)\n", g_h298q);
  evo::Pool pool(nthr); E.pool = &pool; evo::State S;
  struct KGuard { ~KGuard() { g_kfill = nullptr; g_kfin = nullptr; } } kguard_;
  std::unique_ptr<evo::Pool> pool2;
  const int np2 = getenv("H412P2") ? std::max(1, atoi(getenv("H412P2"))) : std::max(1, nthr);
  if (g_h217) { pool2.reset(new evo::Pool(np2)); g_pool2 = pool2.get(); g_W2.resize(np2);
    for (Worker& w : g_W2) { w.C = new Ctx(); w.H = new nf16::Handle(*E.SH, w.C); w.CH = new cmi6::Handle(*E.CS, *E.MM); } }
  std::vector<Meta> pm, nm(NP), cmeta;
  std::vector<Meta> h673_cm;
  auto h673_next_children = [&](const uint32_t* pidx_, uint32_t* cs_) {
    const size_t NS3 = (size_t)3 * NP; std::vector<uint8_t>& u = g_sel_h92.used; u.assign(NS3, 0);
    for (int k = 0; k < NP; k++) u[pidx_[k]] = 1;
    for (int i = 0; i < NP; i++) if (!u[cs_[i]]) u[cs_[i]] = 2;
    size_t c = 0; for (size_t s = 0; s < NS3 && c < (size_t)NP; s++) if (!u[s]) cs_[c++] = (uint32_t)s;
    if (c != (size_t)NP) { fprintf(stderr, "H673 next_children: only %zu free slots\n", c); exit(1); }
  };
  if (getenv("H315A") && atoi(getenv("H315A")) != 0) { if (g_h308 < 1) { printf("H315A: needs H308>=1\n"); return 1; } nf308::g_h315a = true; }
  printf("H315A: %d\n", nf308::g_h315a ? 1 : 0);
  if (getenv("H315T") && atoi(getenv("H315T")) != 0) nf308::g_h315t = true;
  if (g_h320_path && !g_h320) { printf("H320_DUMP: refused (needs H320=1)\n"); return 1; }
  if (g_h320) {
    g_h320_link.assign((size_t)NP, evo::ChildLink320{}); evo::g_link320 = g_h320_link.data();
    if (g_h320_path) { g_h320_f = fopen(g_h320_path, "wb"); if (!g_h320_f) { printf("H320_DUMP: cannot open %s\n", g_h320_path); return 1; }
      char mg[16] = {'H', '3', '2', '0', 'D', '1', 0}; fwrite(mg, 1, 16, g_h320_f); nf308::g_h320_on = true; }
    printf("H320: 1 (breed links %zu B/child x %d, dump %s from step %d)\n", sizeof(evo::ChildLink320), NP, g_h320_path ? g_h320_path : "off", g_h320_from);
  } else printf("H320: 0 (off)\n");
  std::vector<uint32_t> pidx(NP), pidx2(NP), cs(NP);
  if (g_h740) {
    if (g_h774 == 2) { fprintf(stderr, "H740 is refused with H774=2 (the verify select writes the shared flat arrays)\n"); return 2; }
    g_h740e.assign((size_t)3 * NP, H740E{0.0f, 0, 0, 0}); g_h740r.assign((size_t)3 * NP, H740R{0.0, 0, {0, 0, 0}});
    printf("H740: %d (flat hot fields: %zu KiB err/cx/valid + %zu KiB rank/crowd)\n", g_h740, g_h740e.size() * sizeof(H740E) / 1024, g_h740r.size() * sizeof(H740R) / 1024);
  }
  g_sp_on = g_sp_env != 0;
  if (g_sp_on && (g_h740 != 1 || g_ovl == 2 || g_kesel)) { fprintf(stderr, "SEARCH_PEN: needs H740=1, STEP_OVL!=2, K_ESEL=0 (the child install writes the level); off\n"); g_sp_on = false; }
  { const bool pk_nml = getenv("FINAL_PICK") && g_sicis_en != nullptr && !(getenv("FINAL_PICK_NML") && atoi(getenv("FINAL_PICK_NML")) == 0), pk_kfit = getenv("FP_KFIT") && atoi(getenv("FP_KFIT")) != 0;
    if (g_sp_on && !(pk_nml && pk_kfit)) { fprintf(stderr, "SEARCH_PEN: off (the pick does not charge this penalty: needs FINAL_PICK, the SI-CIS NML pick, FP_KFIT=1; have FINAL_PICK %d NML %d FP_KFIT %d)\n", getenv("FINAL_PICK") != nullptr, (int)pk_nml, (int)pk_kfit); g_sp_on = false; } }
  evo::g_pen_tab = g_sp_on ? sp::tab().val.data() : nullptr; g_sp_kunk.store(0);
  printf("SEARCH_PEN: %d (2nd NSGA objective %s; diag %d)\n", (int)g_sp_on, g_sp_on ? "= the pick's penalty 0.5 k log(48/2pi) + struct(nx) + struct(ny), 16-bit levels, Fenwick sweep" : "= node count (6-bit levels)", g_sp_diag);
  FlatKeySet popset;
  std::vector<float> err(NP), rx(NP), ry(NP);
  FILE* traj = fopen((outdir + "/traj.csv").c_str(), "w");
  if (traj) fprintf(traj, "step,children,admitted,rej_within,rej_pop,shortfall,uniq_trees,fits_x,fits_y,sel_skip,sel_skip_mig,surv_child,xfill,"
                          "pop_keys,pop_nf_pairs,pop_trees,best_child_err,best_child_cx,best_pop_err,best_pop_cx,"
                          "t_breed,t_nf,t_dedup,t_comp,t_scan,t_gn,t_glue,t_select,t_timed,children_per_s\n");
  auto dump_gen = [&](int g) { char p[1024]; snprintf(p, sizeof p, "%s/pop_%03d.bin", outdir.c_str(), g); FILE* f = fopen(p, "wb"); if (!f) return;
    for (int i = 0; i < NP; i++) { uint64_t kh = KeyHash()(pm[i].key); uint32_t cx; float e; if (g_h740 == 1) { cx = g_h740e[pidx[i]].cx; e = g_h740e[pidx[i]].err; } else { const evo::Ind& PI = XC(S, pidx[i]); cx = PI.cx; e = PI.err; }
      fwrite(&kh, 8, 1, f); fwrite(&pm[i].r0, 8, 1, f); fwrite(&pm[i].r1, 8, 1, f); fwrite(&e, 4, 1, f); fwrite(&cx, 4, 1, f); }
    fclose(f); };


  {
    const long h178_g0 = g_h178_created;
    Times T0; std::vector<evo::Cand> popG; std::vector<float> popE, popRx, popRy; int nb_batches = 0; auto g0 = clk::now();
    while ((int)popG.size() < NP) {
      if (++nb_batches > 200) { fprintf(stderr, "gen 0: could not fill the population\n"); return 1; }
      std::vector<evo::Cand> bc; next_batch(src, NP, bc);
      for (auto& c : bc) sf_pin_driver(c);
      const int nb = (int)bc.size();
      std::vector<uint8_t> bv(nb, 1), adm; std::vector<float> be(nb), brx(nb), bry(nb); std::vector<Meta> bm;
      StepBatch SBx;
      std::vector<const evo::Cand*> bcp(nb); for (int i = 0; i < nb; i++) bcp[i] = &bc[i];
      evaluate_v7(E, bcp.data(), bv.data(), nb, popset, be.data(), brx.data(), bry.data(), adm, bm, T0, SBx);
      for (int i = 0; i < nb && (int)popG.size() < NP; i++) if (adm[i]) {
        popG.push_back(bc[i]); popE.push_back(be[i]); popRx.push_back(brx[i]); popRy.push_back(bry[i]); popset.insert(bm[i].kh, bm[i].key.data(), (uint32_t)bm[i].key.size()); pm.push_back(bm[i]); }
    }
    const double tot = ms_since(g0);
    int bi = (int)(std::min_element(popE.begin(), popE.end()) - popE.begin());
    printf("gen 0 (seeds %s): batches %d, file records read %ld (valid %ld), random made %ld | admitted %.0f, rej within %.0f, rej vs pop %.0f | population %d, distinct keys %zu | eval %.1f ms | best err %.6g (cx %d)\n",
           random_seeds ? "random" : "file", nb_batches, src.file_read, src.file_valid, src.random_made, T0.admitted, T0.rej_within, T0.rej_pop,
           NP, popset.size(), tot, popE[bi], evo::cplx(popG[bi]));
    {
      const char* gd = getenv("G0DUMP");
      if (gd && atoi(gd) != 0) {
        std::vector<int> ord(NP); for (int i = 0; i < NP; i++) ord[i] = i;
        std::vector<uint64_t> gk(NP); for (int i = 0; i < NP; i++) gk[i] = fitx::gkey(popG[i]);
        auto prn = [&](const char* tag, int rk, int i) { printf("G0DUMP %s %d idx %d key %016llx err %.9g rfin_x %.9g rfin_y %.9g cx %d\n", tag, rk, i, (unsigned long long)gk[i], popE[i], popRx[i], popRy[i], (int)evo::cplx(popG[i])); };
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return gk[a] != gk[b] ? gk[a] < gk[b] : a < b; });
        for (int r = 0; r < NP && r < 4096; r++) prn("key", r, ord[r]);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return popE[a] != popE[b] ? popE[a] < popE[b] : gk[a] < gk[b]; });
        for (int r = 0; r < NP && r < 50; r++) prn("best", r, ord[r]);
        printf("G0DUMP end NP %d\n", NP); fflush(stdout);
      }
    }
    if (h178_on() || h178_pre()) printf("H178 gen0: PSOs created %ld\n", g_h178_created - h178_g0);
    for (int i = 0; i < NP; i++) {
      evo::Ind& I = S.par[i]; I.g = popG[i]; I.valid = 1; std::memset(I.res, 0, sizeof I.res); I.res[0] = popE[i]; I.err = popE[i]; I.cx = sp_cx(I.g, pm[i]); }
    evo::init_rank(pool, S.par.data(), S.sel.data());
    if (g_h740) for (int i = 0; i < NP; i++) { const evo::Ind& I = S.par[i]; g_h740e[i] = H740E{I.err, I.cx, I.valid, 0}; g_h740r[i] = H740R{I.crowd, I.rank, {0, 0, 0}}; }
    for (int i = 0; i < NP; i++) pidx[i] = (uint32_t)i;
    next_children_h92(pidx.data(), cs.data());
    dump_gen(0);
  }

  const int steady_from = steps >= 15 ? 11 : 1;
  Times TT, TS; long nsteady = 0; double best_all = INFINITY;
  long h99_sl = 0, h99_sh = 0; double h99_su = 0;
  std::vector<evo::Cand> ch(NP); std::vector<uint8_t> cv(NP);
  std::vector<const evo::Cand*> chp(NP);
  H896Pop h896pop; g_h896p = &h896pop;
  std::thread pop_thr; std::function<void()> pend_print;
  std::thread h4_thr; std::shared_ptr<H4Job> Jh4;

  if ((g_h131 & 1) && !(g_freet && (evo::mem_lean() & 4))) h131_prealloc(E.dev, (size_t)NP);
  g_h301_main = std::this_thread::get_id(); (void)h401_tid();
  double h301_t0 = 0, h301_busy = 0, h301_idle = 0; long h301_n = 0; std::map<std::string, double> h301_acc;
  std::thread sb_thr, mb_thr;
  std::thread h652_thr; double h652_best = INFINITY, h652_join_ms = 0, h651_join_ms = 0;
  if (pf_thr.t.joinable()) pf_thr.t.join();
  uib::g_t0 = h301_now();
  if (getenv("UI_STEP_LOG") && atoi(getenv("UI_STEP_LOG"))) printf("UI_STEP gen0_done t_ms %.1f\n", ui_ms_since_start());
  for (int st = 1; st <= steps && !uib::g_stop; st++) {
    if (g_h313x) { const bool on = ((st + (g_h313x == 2 ? 1 : 0)) & 1) != 0; c308::g_a2 = c720::g_a2 = on && (g_h313a_mask & 1); cmi6::g_a1 = on && (g_h313a_mask & 4); }
    if (g_h308x) g_t4 = (g_h308 >= 3) && (((st + (g_h308x == 2 ? 1 : 0)) & 1) != 0);
    static double h401_t0 = 0;
    if (g_h892) { static const char* const L[] = {"", "H892 STEP 1 begin", "H892 STEP 2 begin", "H892 STEP 3 begin", "H892 STEP 4 begin"};
      if (st == 1) { g_h892_ring = (H401Rec*)calloc(g_h892_cap, sizeof(H401Rec)); evo::g_pool_trace_on = h401_pool_on; evo::g_pool_trace_now = h301_now; evo::g_pool_trace_rec = h401_pool_rec; evo::g_pool_trace_wait = h401_pool_wait; g_h892_t0 = h301_now(); g_h401_on.store(1); }
      if (st >= 1 && st <= 4) { const double t = h301_now(); h401_push(0, L[st], t, t); }
      else if (st == 5) { g_h401_on.store(0); const double t = h301_now(); h401_push(0, "H892 STEP 4 end", t, t); } }
    if (g_h401) { if (st == g_h401) { evo::g_pool_trace_on = h401_pool_on; evo::g_pool_trace_now = h301_now; evo::g_pool_trace_rec = h401_pool_rec; evo::g_pool_trace_wait = h401_pool_wait; h401_t0 = h301_now(); g_h401_on.store(1); }
      else if (st == g_h401 + 1) { g_h401_on.store(0); const double t = h301_now(); h401_push(0, "(H401 window end)", t, t); }
      else if (st == g_h401 + 2) h401_dump(h401_t0); }
    if (g_h900 & 1) h900_resident(g_h900_st);
    g_h670_st = st; g_h900_st = st;
    if (g_h917 & 1) { if (g_h917_pool) objc_autoreleasePoolPop(g_h917_pool); g_h917_pool = objc_autoreleasePoolPush(); }
    if ((g_h917 & 2) && st % g_h917_every == 0) { const double t0 = now_ms_h917(); malloc_zone_pressure_relief(nullptr, 0); g_h917_relief_ms += now_ms_h917() - t0; g_h917_relief_n++; }
    h301_mark("step boundary (memo)");
    h141_step_reset(); h150_step_reset();
    const long h178_c0 = g_h178_created;
    Times T;

    const bool h651_lazy = g_h651 && mb_thr.joinable();
    if (!h651_lazy) { if (mb_thr.joinable()) { const double tj0 = h301_now(); if (g_h886 & 2) g_h886_boost.store(1, std::memory_order_relaxed); mb_thr.join(); g_h854_join_wait = (h301_now() - tj0) * 1e3; g_h854_join_done = h301_now(); g_h854_joined = true; } else { const double tj0 = h301_now(); memo_boundary(E, st); if (g_h854 & 1) { g_h854_rec.launch = tj0; g_h854_join_wait = (h301_now() - tj0) * 1e3; g_h854_join_done = h301_now(); g_h854_joined = (st > 1); } }
    h886_adopt(); h156_step_begin(E, st); }
    const auto t_iter0 = clk::now();
    if (g_h301) { const double tn = h301_now();
      if (h301_t0 > 0) { std::map<std::string, double> acc; double idle = 0; const double busy = h301_window(h301_t0, tn, acc, idle);
        if (st - 1 >= steady_from) { h301_busy += busy; h301_idle += idle; h301_n++; for (auto& kv : acc) h301_acc[kv.first] += kv.second; }
        printf("H301 step %2d: window %.1f ms | GPU busy %.1f | idle %.1f\n", st - 1, (tn - h301_t0) * 1e3, busy, idle); }
      h301_t0 = tn; }
    h301_mark("breed + ch[]");
    long nf_h0 = 0, nf_m0 = 0, nf_d0 = 0; if (!h651_lazy) memo_counts(E, nf_h0, nf_m0, nf_d0);
    const long h99_l0 = G99.lookups, h99_h0 = G99.hits, h99_i0 = G99.inserts, h99_e0 = G99.evicts, h99_r0 = G99.refused, h99_d0 = G99.dups;
    auto t = clk::now();
    if (g_h826) {
      g_h826_head = [&, st](int b1_0) {
        const int d0 = std::min(evo::NDEME - 1, (b1_0 - 1) / evo::DECK);
        const auto th = clk::now();
        h826_breed(pool, S, pidx.data(), cs.data(), S.BB, 7, (uint64_t)st, S.rb.data(), S.fl.data(), chp.data(), cv.data(), 0, d0 + 1);
        T.breed += ms_since(th);
        if (d0 + 1 < evo::NDEME) g_h826_rest = [&, st, d0]() {
          const auto tr = clk::now();
          h826_breed(pool, S, pidx.data(), cs.data(), S.BB, 7, (uint64_t)st, S.rb.data(), S.fl.data(), chp.data(), cv.data(), d0 + 1, evo::NDEME);
          T.breed += ms_since(tr);
        };
      };
    } else {
    if (g_h777) breed_h777(pool, S, pidx.data(), cs.data(), S.BB, 7, (uint64_t)st, S.rb.data(), S.fl.data());
    else breed_h92(pool, S, pidx.data(), cs.data(), S.BB, 7, (uint64_t)st, S.rb.data(), S.fl.data());
    if (!(g_h257 & 8)) T.breed = ms_since(t);
    if (g_h262 & 1) {
      for (int i = 0; i < NP; i++) { chp[i] = &XC(S, cs[i]).g; cv[i] = g_h740 ? h740_valid(S, cs[i]) : XC(S, cs[i]).valid; }
    } else {
      for (int i = 0; i < NP; i++) { ch[i] = XC(S, cs[i]).g; chp[i] = &ch[i]; cv[i] = g_h740 ? h740_valid(S, cs[i]) : XC(S, cs[i]).valid; }
    }
    if (g_h257 & 8) T.breed = ms_since(t);
    }
    if (pop_thr.joinable()) { auto tj = clk::now(); pop_thr.join(); T.rebuild += ms_since(tj); }
    if (pend_print && !(g_h131 & 4)) { pend_print(); pend_print = nullptr; }
    std::vector<uint8_t> adm; std::unique_ptr<StepBatch> SBp;
    if (g_h415 && g_h415_spare) { if (sb_thr.joinable()) { auto tj = clk::now(); sb_thr.join(); g_h415_join_ms += ms_since(tj); } SBp.reset(g_h415_spare); g_h415_spare = nullptr; }
    else SBp.reset(new StepBatch());
    StepBatch& SB = *SBp;
    if (h652_thr.joinable()) { auto tj = clk::now(); h652_thr.join(); h652_join_ms = ms_since(tj); }
    if (h651_lazy) { h301_mark("H651 memo join"); auto tj = clk::now(); if (g_h886 & 2) g_h886_boost.store(1, std::memory_order_relaxed); mb_thr.join(); h651_join_ms = ms_since(tj); h886_adopt(); h156_step_begin(E, st); memo_counts(E, nf_h0, nf_m0, nf_d0); }
    h301_mark("evaluate entry");
    if (g_h301mb == 2 && st < steps) g_h301_after_front = [&E, &mb_thr, st] { { const double tl = h301_now(); g_h886_boost.store(0, std::memory_order_relaxed); mb_thr = std::thread([&E, st, tl] { h854_qos_self(); H401S h401s_("memo_boundary thread"); if (g_h854 & 1) { memo_boundary(E, st + 1); g_h854_rec.launch = tl; } else memo_boundary(E, st + 1); }); } };
    if (g_h320_f) { std::lock_guard<std::mutex> lk(nf308::g_h320_mx); nf308::g_h320_recs.clear(); }
    const long mf874a = (g_h874 & 16) ? h874_mf() : 0;
    int ovl_nd = 0;
    if ((g_ovl == 2 && g_h774 == 1 && g_freet && g_freet_async) || (g_kesel && g_h774 == 1 && !g_freet && g_sf_axis && g_h740)) {
      g_ovl_hook = [&](int b0) {
        const int nd = std::min(evo::NDEME, b0 / evo::DECK); if (nd <= ovl_nd) return; const int nd0 = ovl_nd;
        pool.run(nd - nd0, [&](int q) { const int d = nd0 + q; for (int i = d * evo::DECK; i < (d + 1) * evo::DECK; i++) {
          evo::Ind& c = XA(S, cs[i]); c.valid = adm[i]; if (!g_h419) std::memset(c.res, 0, sizeof c.res); c.res[0] = err[i]; c.cx = (uint16_t)evo::cplx(c.g);
          if (g_h740) { if (adm[i]) c.err = 0.0f + err[i]; H740E& fe = g_h740e[cs[i]]; fe.err = 0.0f + err[i]; fe.cx = c.cx; fe.valid = adm[i]; } } });
        if (nd0 == 0) sel774_begin(pool, S, pidx.data());
        sel774_demes(pool, S, pidx.data(), cs.data(), pm, cmeta, adm.data(), nd0, nd);
        ovl_nd = nd;
      };
    }
    h301_mark("K: pin drivers");
    if (g_sf_axis && g_sf_pin && !(g_h826 && g_pinfix)) {
      static const int kpin = getenv("K_PIN") ? atoi(getenv("K_PIN")) : 1;
      if (kpin && NP >= 8192) { const int nb = (NP + 4095) / 4096; pool.run(nb, [&](int b) { for (int i = b * 4096; i < std::min(NP, (b + 1) * 4096); i++) if (cv[i]) sf_pin_driver(const_cast<evo::Cand&>(*chp[i])); }); }
      else for (int i = 0; i < NP; i++) if (cv[i]) sf_pin_driver(const_cast<evo::Cand&>(*chp[i]));
    }
    h301_mark("K: evaluate_v7 setup");
    evaluate_v7(E, chp.data(), cv.data(), NP, popset, err.data(), rx.data(), ry.data(), adm, cmeta, T, SB);
    const long mf874b = (g_h874 & 16) ? h874_mf() : 0;
    if (g_h826 && (g_h826_head || g_h826_rest)) { fprintf(stderr, "H826: a breed hook was not consumed by evaluate_v7\n"); return 2; }
    if (g_h320_f) h320_flush_step(st, S);
    bool h870_early = false;
    if ((g_h870 & 4) && g_h415 && g_h301r && ncheck < 2) {
      std::swap(g_h301r_t, SB.UQ.U); std::swap(g_h301r_e, SB.UQ.enc); g_h301r_n.store(0);
      if (sb_thr.joinable()) sb_thr.join();
      StepBatch* p = SBp.release(); g_h415_spare = p;
      sb_thr = std::thread([p] { H401S h401s_("StepBatch clear thread"); p->UQ.U.clear(); p->UQ.enc.clear(); p->UQ.hash.clear(); });
      h870_early = true;
    }
    if (g_pincensus && g_sf_axis) { long nc = 0, nnd = 0; for (int i = 0; i < NP; i++) if (cv[i]) { nc++; if (!sf_is_driver(*chp[i])) nnd++; } printf("PIN_CENSUS step %d: evaluated children %ld, genotype changes after evaluation %ld (PIN_FIX=%d H826=%d)\n", st, nc, nnd, g_pinfix, g_h826); }
    long nf_h1 = 0, nf_m1 = 0, nf_d1 = 0; memo_counts(E, nf_h1, nf_m1, nf_d1);
    T.nf_hit = (double)(nf_h1 - nf_h0); T.nf_miss = (double)(nf_m1 - nf_m0); T.nf_drop = (double)(nf_d1 - nf_d0);
    if (g_h301mb == 1 && st < steps) mb_thr = std::thread([&E, st] { memo_boundary(E, st + 1); });
    T.shortfall = NP - T.admitted;
    h301_mark("install");
    t = clk::now();
    if (g_h131 & 2) {
      const int blk = (g_h405 || (g_h921 & 4)) ? 2048 : 8192, nbk = (NP + blk - 1) / blk;
      pool.run(nbk, [&](int b) {
        for (int i = std::max(b * blk, ovl_nd * evo::DECK); i < std::min(NP, (b + 1) * blk); i++) {
          if ((g_h921 & 4) && i + 8 < NP) { evo::Ind& nx = XA(S, cs[i + 8]); __builtin_prefetch(&nx.g.len[0], 1, 3); __builtin_prefetch(&nx.res[0], 1, 3); __builtin_prefetch(&nx.cx, 1, 3); __builtin_prefetch(&nx.valid, 1, 3); }
          evo::Ind& c = XA(S, cs[i]); c.valid = adm[i]; if (!g_h419) std::memset(c.res, 0, sizeof c.res); c.res[0] = err[i]; c.cx = sp_cx(c.g, cmeta[i]);
          if (g_h740) { if (adm[i]) c.err = 0.0f + err[i]; H740E& fe = g_h740e[cs[i]]; fe.err = 0.0f + err[i]; fe.cx = c.cx; fe.valid = adm[i]; } }
      });
      T.h131_inst += ms_since(t);
    } else {
      for (int i = ovl_nd * evo::DECK; i < NP; i++) { evo::Ind& c = XA(S, cs[i]); c.valid = adm[i]; if (!g_h419) std::memset(c.res, 0, sizeof c.res); c.res[0] = err[i]; c.cx = sp_cx(c.g, cmeta[i]);
          if (g_h740) { if (adm[i]) c.err = 0.0f + err[i]; H740E& fe = g_h740e[cs[i]]; fe.err = 0.0f + err[i]; fe.cx = c.cx; fe.valid = adm[i]; } }
      T.h131_inst += ms_since(t);
    }
    T.glue += ms_since(t);
    if (g_h131 & 4) {
      auto tj = clk::now(); if (h4_thr.joinable()) h4_thr.join(); T.h131_h4wait += ms_since(tj);
      if (Jh4) { T.h131_h4 += Jh4->dump_ms + Jh4->sets_ms; Jh4.reset(); }
      if (pend_print) { pend_print(); pend_print = nullptr; }
    }
    if (h652_thr.joinable()) { h652_thr.join(); }
    if (g_h670) {
      if (g_h670_have && g_h670_jst != g_h670_st - 1) { g_h670_have = false; printf("H670 apply: stale results (step %d at step %d) dropped\n", g_h670_jst, g_h670_st); }
      if (g_h670_have) { g_h670_have = false; const auto ta = clk::now();
        std::vector<uint8_t> f((size_t)3 * NP, 0); for (int i = 0; i < NP; i++) f[pidx[i]] = 1;
        long nap = 0, nlost = 0; const bool rev = g_h673 && g_h670_st <= steps;
        std::vector<std::vector<std::pair<float, int>>> rv(rev ? (size_t)evo::NDEME : 0);
        for (const auto& pr : g_h670_res) { best_all = std::min(best_all, (double)pr.second);
          const uint32_t slot = g_h670_cs[(size_t)pr.first]; if (!f[slot]) { nlost++; if (rev) rv[(size_t)(pr.first / evo::DECK)].push_back({pr.second, pr.first}); continue; }
          evo::Ind& c = XA(S, slot); c.res[0] = pr.second; c.err = 0.0f + pr.second; if (g_h740) g_h740e[slot].err = c.err; nap++; }
        if (rev) { long nrev = 0, ncand = 0;
          for (int d = 0; d < evo::NDEME; d++) { auto& v = rv[(size_t)d]; if (v.empty()) continue; ncand += (long)v.size();
            std::sort(v.begin(), v.end()); size_t t = 0;
            for (int i2 = d * evo::DECK; i2 < (d + 1) * evo::DECK && t < v.size(); i2++) {
              if (adm[(size_t)i2]) continue;
              const int i = v[t].second; const uint32_t slot = g_h670_cs[(size_t)i];
              evo::Ind& c = XA(S, slot); c.res[0] = v[t].first; c.err = 0.0f + v[t].first; c.valid = 1; if (g_h740) { g_h740e[slot].err = c.err; g_h740e[slot].valid = 1; }
              cs[(size_t)i2] = slot; adm[(size_t)i2] = 1; cmeta[(size_t)i2] = std::move(h673_cm[(size_t)i]);
              err[(size_t)i2] = v[t].first; chp[(size_t)i2] = &XA(S, slot).g; t++; nrev++; } }
          printf("H673 revive: dropped accepted %ld | revived %ld\n", ncand, nrev); }
        printf("H670 apply: accepted %zu | survivors updated %ld | not surviving %ld | %.2f ms\n", g_h670_res.size(), nap, nlost, ms_since(ta)); } }
    const bool kpop = fitx::on() && g_kpopkey && !g_freet && NP >= 4096 && (g_h92 & 1) && g_sf_axis;
    std::vector<uint32_t> kcsN; if (kpop) kcsN = cs;
    h301_mark("select");
    t = clk::now();
    g_h128_prof_on = (st >= steady_from);
    if (g_h774) h774_run(pool, S, pidx.data(), cs.data(), pidx2.data(), pm, cmeta, adm.data(), nm, T, (g_h92 & 1) != 0, st, g_h774_cstart, ovl_nd > 0 ? ovl_nd : -1);
    else select_h92(pool, S, pidx.data(), cs.data(), pidx2.data(), pm, cmeta, adm.data(), nm, T, (g_h92 & 1) != 0);
    g_h128_prof_on = false;
    const bool h673r = g_h673 && g_h670j.pending && !g_h670j.cs_ok;
    if (g_h670 && g_h670j.pending && !g_h670j.cs_ok) { g_h670j.cs = cs; g_h670j.cs_ok = true; }
    std::swap(pidx, pidx2);
    if (h673r) { std::swap(cmeta, h673_cm); h673_next_children(pidx.data(), cs.data()); } else next_children_h92(pidx.data(), cs.data());
    std::swap(pm, nm);
    if (g_h890 & 16) for (auto& blk : NFD.blk) for (auto& w : blk.w) w = 0xDEADBEEFu;
    T.select = ms_since(t);
    if (g_pincensus && g_sf_axis) { long nnd = 0; for (int i = 0; i < NP; i++) if (!sf_is_driver(XC(S, pidx[i]).g)) nnd++; printf("PIN_CENSUS step %d: population %d non-driver driver side %ld\n", st, NP, nnd); }
    h301_mark("after select: dump/popset/prints");
    if (h650_is("tail")) h650_spin_ns(g_h650ms * 1e6);

    auto J = std::make_shared<H4Job>();
    auto h4_body = [&, J, st]() {
      auto td = clk::now(); dump_gen(st); J->dump_ms = ms_since(td);
      auto ts = clk::now();
      std::set<std::pair<uint64_t, uint64_t>> pairs; std::unordered_set<uint64_t> trees;
      for (int i = 0; i < NP; i++) { pairs.insert({pm[i].r0, pm[i].r1}); trees.insert(pm[i].r0); trees.insert(pm[i].r1); }
      J->npair = pairs.size(); J->ntree = trees.size(); J->sets_ms = ms_since(ts);
    };
    if (getenv("H231") && atoi(getenv("H231")) == 1) {  }
    else if (g_h131 & 4) { h4_thr = std::thread([h4_body] { H401S h401s_("h4 dump+sets thread"); k_qos_set(2); h4_body(); }); Jh4 = J; }
    else { h4_body(); T.dump = J->dump_ms; T.h131_h4 += J->dump_ms + J->sets_ms; }

    t = clk::now();
    bool h921_hdl = false;
    if (g_h743 && (g_h921 & 3) && !(g_h896 & 1)) {
      auto& BK = g_h921_bk; if (BK.size() != 512) BK.resize(512);
      auto& BA = g_h921_ba; if (BA.size() != 256) BA.resize(256);
      auto& LA = g_h921_la; if (LA.size() != 16) LA.resize(16);
      auto& code = g_h921_code; auto& code2 = g_h921_code2; if (code.size() != (size_t)NP) { code.assign((size_t)NP, 0); code2.assign((size_t)NP, 0); }
      const bool inc1 = (g_h921 & 1) != 0; const bool one2 = (g_h921 & 3) != 0;
      auto& BW = g_h921_bw; if (BW.size() != 256) BW.assign(256, 0);
      const uint32_t EM = 0x0FFFFFFFu;
      auto shard_ins = [&](int s) { FlatKeySet& X = popset.sh[(size_t)s]; X.clear(); X.reserve((size_t)NP / 8);
        if (inc1) { size_t rows = 0, words = 0; for (int b = 0; b < 16; b++) { rows += BK[((size_t)b * 16 + (size_t)s) * 2].size(); words += BW[(size_t)b * 16 + (size_t)s]; }
          const size_t ne = rows * 5 / 2 + 16; X.eh.reserve(ne); X.eoff.reserve(ne); X.elen.reserve(ne); X.arena.reserve(words * 5 / 2 + 64); }
        auto ins = [&](int i) { bool f = false; const uint32_t e = X.find_or_add(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size(), f); code[(size_t)i] = ((uint32_t)s << 28) | e; };
        if (one2) { for (int b = 0; b < 16; b++) for (uint32_t i : BK[((size_t)b * 16 + (size_t)s) * 2]) ins((int)i); }
        else { for (int i = 0; i < NP; i++) if (FlatKeySet::shard_of(pm[i].kh) == (unsigned)s) ins(i); }
        if (inc1) { X.inc = true; X.dead.assign(X.eh.size(), 0); X.ndead = 0; X.dead.reserve(X.eh.capacity()); } };
      auto h921_full = [&]() {
        if (popset.sh.empty()) { popset.clear(); popset.sh.resize(16); }
        if (one2) PL(E).run(16, [&](int b) { for (int s = 0; s < 16; s++) { BK[((size_t)b * 16 + (size_t)s) * 2].clear(); BW[(size_t)b * 16 + (size_t)s] = 0; }
          const int lo = (int)((long)NP * b / 16), hi = (int)((long)NP * (b + 1) / 16);
          for (int i = lo; i < hi; i++) { const unsigned s = FlatKeySet::shard_of(pm[i].kh); BK[((size_t)b * 16 + s) * 2].push_back((uint32_t)i); BW[(size_t)b * 16 + s] += pm[i].key.size(); } });
        PL(E).run(16, shard_ins);
        g_h921_have = inc1; };
      bool done = false; size_t nadd = 0;
      if (inc1 && g_h921_have && !g_h774 && popset.sh.size() == 16 && g_sel_h92.plist.size() == (size_t)NP) {
        auto& keep = g_h921_keep; const auto& plist = g_sel_h92.plist;
        keep.assign((size_t)NP, 0);
        PL(E).run(16, [&](int b) { LA[(size_t)b].clear();
          for (size_t r = plist.size() * (size_t)b / 16; r < plist.size() * (size_t)(b + 1) / 16; r++) { const auto& q = plist[r];
            if (q.p < (uint32_t)NP) { keep[q.p] = 1; code2[q.dst_slot] = code[q.p]; } else LA[(size_t)b].push_back(q.dst_slot); } });
        for (int b = 0; b < 16; b++) nadd += LA[(size_t)b].size();
        size_t ent = 0; for (const auto& x : popset.sh) ent += x.eh.size();
        if (ent + nadd <= (size_t)NP * 5 / 2) {

          PL(E).run(16, [&](int b) { for (int s = 0; s < 16; s++) { BK[((size_t)b * 16 + (size_t)s) * 2].clear(); BA[(size_t)b * 16 + (size_t)s].clear(); }
            for (size_t p = (size_t)NP * (size_t)b / 16; p < (size_t)NP * (size_t)(b + 1) / 16; p++) if (!keep[p]) { const uint32_t c = code[p]; BK[((size_t)b * 16 + (c >> 28)) * 2].push_back(c & EM); }
            for (uint32_t idx : LA[(size_t)b]) { const uint64_t kh = pm[idx].kh; BA[(size_t)b * 16 + FlatKeySet::shard_of(kh)].push_back({kh, idx}); } });

          std::atomic<int> bad{0};
          PL(E).run(16, [&](int s) { FlatKeySet& X = popset.sh[(size_t)s];
            for (int b = 0; b < 16; b++) for (uint32_t e : BK[((size_t)b * 16 + (size_t)s) * 2]) { if (e >= X.dead.size() || X.dead[e]) bad.fetch_add(1); else { X.dead[e] = 1; X.ndead++; } }
            for (int b = 0; b < 16; b++) { const auto& v = BA[(size_t)b * 16 + (size_t)s]; const size_t nv = v.size();
              for (size_t j = 0; j < nv; j++) {
                if (j + 8 < nv) X.pf_slot(v[j + 8].first);
                if (j + 6 < nv) __builtin_prefetch(&pm[v[j + 6].second], 0, 3);
                if (j + 3 < nv) { const uint32_t* kp = pm[v[j + 3].second].key.data(); __builtin_prefetch(kp, 0, 3); __builtin_prefetch(kp + 16, 0, 3); __builtin_prefetch(kp + 32, 0, 3); }
                const Meta& m = pm[v[j].second]; bool f = false;
                const uint32_t e = X.find_or_add(m.kh, m.key.data(), (uint32_t)m.key.size(), f);
                if (f) X.dead.push_back(0); else if (X.dead[e]) { X.dead[e] = 0; X.ndead--; } else bad.fetch_add(1);
                code2[v[j].second] = ((uint32_t)s << 28) | e; } } });
          if (bad.load() == 0) { done = true; std::swap(code, code2); } else fprintf(stderr, "H921: step %d popset anomalies %d: full rebuild\n", st, bad.load());
        }
      }
      if (!done) h921_full();
      h921_hdl = true;
      if (g_h921 & 16) printf("H921 step %d: popset %s add %zu | %.2f ms\n", st, done ? "incremental" : "full", nadd, ms_since(t));
      if (g_h921 & 8) {
        FlatKeySet R; R.reserve((size_t)NP * 2); for (int i = 0; i < NP; i++) R.insert(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size());
        size_t miss = 0, badc = 0; for (int i = 0; i < NP; i++) { if (!popset.contains(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size())) miss++;
          if (inc1) { const FlatKeySet& X = popset.sh[code[(size_t)i] >> 28]; const uint32_t e = code[(size_t)i] & EM; if (e >= X.eh.size() || X.eh[e] != pm[i].kh || X.dead[e] || !X.same(e, pm[i].key.data(), (uint32_t)pm[i].key.size())) badc++; } }
        printf("H921 verify step %d: ref %zu popset %zu missing %zu bad codes %zu %s\n", st, R.size(), popset.size(), miss, badc, (R.size() == popset.size() && miss == 0 && badc == 0) ? "OK" : "BAD");
        if (!(R.size() == popset.size() && miss == 0 && badc == 0)) exit(3);
      }
    }
    if (h921_hdl) T.rebuild += ms_since(t); else
    if (g_h743) {
      if (popset.sh.empty()) { popset.clear(); popset.sh.resize(16); }
      if (g_h896 & 1) {
        h896pop.join();
        const int nt = (int)std::max(1u, std::min(15u, std::thread::hardware_concurrency() - 1));
        auto nxt = std::make_shared<std::atomic<int>>(0);
        for (int w = 0; w < nt; w++) h896pop.th.emplace_back([&popset, &pm, NP, nxt] { H401S h401s_("H896 popset thread");
          for (int s; (s = nxt->fetch_add(1)) < 16;) { FlatKeySet& X = popset.sh[(size_t)s]; X.clear(); X.reserve((size_t)NP / 8);
            for (int i = 0; i < NP; i++) if (FlatKeySet::shard_of(pm[i].kh) == (unsigned)s) X.insert(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size()); } });
      } else
      PL(E).run(16, [&](int s) { FlatKeySet& X = popset.sh[(size_t)s]; X.clear(); X.reserve((size_t)NP / 8);
        for (int i = 0; i < NP; i++) if (FlatKeySet::shard_of(pm[i].kh) == (unsigned)s) X.insert(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size()); });
      T.rebuild += ms_since(t);
    } else if (g_h117 & 8) {
      pop_thr = std::thread([&popset, &pm, NP]() {
    k_qos_set(2);
        popset.clear(); popset.reserve((size_t)NP * 2);
        for (int i = 0; i < NP; i++) popset.insert(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size());
      });
    } else {
      popset.clear(); popset.reserve((size_t)NP * 2);
      for (int i = 0; i < NP; i++) popset.insert(pm[i].kh, pm[i].key.data(), (uint32_t)pm[i].key.size());
      T.rebuild += ms_since(t);
    }
    h301_mark("after select: best/check");
    int bc = -1, bp = 0; const auto h921_tbc = clk::now();
    if (g_h652) {
      for (int i = 0; i < NP; i++) if (adm[i] && std::isfinite(err[i]) && err[i] < DEAD && (bc < 0 || err[i] < err[bc])) bc = i;
      const double bce = bc >= 0 ? (double)err[bc] : (double)INFINITY; const int bcx = bc >= 0 ? (int)evo::cplx(*chp[bc]) : -1;
      h652_best = bce;
      if (h652_thr.joinable()) h652_thr.join();
      h652_thr = std::thread([&S, &pidx, NP, bce, bcx, st]() {
    k_qos_set(2);
        int p = 0;
        for (int i = 1; i < NP; i++) if (h740_perr(S, pidx[i]) < h740_perr(S, pidx[p])) p = i;
        printf("        best child err %.6g (cx %d) | best pop err %.6g (cx %d) | (H652 helper, step %d)\n", bce, bcx, h740_perr(S, pidx[p]), (int)XC(S, pidx[p]).cx, st);
      });
    } else if (g_h405 || (g_h921 & 4)) {
      const int nb = 64; int bcs[nb], bps[nb]; std::atomic<int> nan_seen{0};
      pool.run(nb, [&](int b) {
        const int lo = (int)((long)NP * b / nb), hi = (int)((long)NP * (b + 1) / nb);
        int c = -1; for (int i = lo; i < hi; i++) if (adm[i] && std::isfinite(err[i]) && err[i] < DEAD && (c < 0 || err[i] < err[c])) c = i;
        int p = -1; for (int i = lo; i < hi; i++) { const double e = h740_perr(S, pidx[i]); if (std::isnan(e)) nan_seen.store(1); if (p < 0 || e < h740_perr(S, pidx[p])) p = i; }
        bcs[b] = c; bps[b] = p; });
      for (int b = 0; b < nb; b++) if (bcs[b] >= 0 && (bc < 0 || err[bcs[b]] < err[bc])) bc = bcs[b];
      if (!nan_seen.load()) { for (int b = 0; b < nb; b++) if (bps[b] >= 0 && h740_perr(S, pidx[bps[b]]) < h740_perr(S, pidx[bp])) bp = bps[b]; }
      else { bp = 0; for (int i = 1; i < NP; i++) if (h740_perr(S, pidx[i]) < h740_perr(S, pidx[bp])) bp = i; }
    } else {
    for (int i = 0; i < NP; i++) if (adm[i] && std::isfinite(err[i]) && err[i] < DEAD && (bc < 0 || err[i] < err[bc])) bc = i;
    for (int i = 1; i < NP; i++) if (h740_perr(S, pidx[i]) < h740_perr(S, pidx[bp])) bp = i;
    }
    best_all = std::min(best_all, bc >= 0 ? (double)err[bc] : (double)INFINITY);
    if (g_h921 & 16) printf("H921 step %d: install %.2f best/check %.2f ms\n", st, T.h131_inst, ms_since(h921_tbc));

    t = clk::now();
    long n[2] = {0, 0}, miss[2] = {0, 0}, better[2] = {0, 0}, smiss[2] = {0, 0}, rmiss = 0, rbetter = 0, rstart = 0; std::vector<double> ratios, sdiff;
    for (int s = 0; s < 2; s++) check_fit(SB, s, E.tg[s], ncheck / 2, 1000 + st * 2 + s, n[s], miss[s], better[s], smiss[s], ratios, rmiss, rbetter, rstart, sdiff);
    double ratio_med = 0; if (!ratios.empty()) { std::sort(ratios.begin(), ratios.end()); ratio_med = ratios[ratios.size() / 2]; }
    T.check = ms_since(t);
    h301_mark("J: after check");
    if (g_jchk) printf("J_CHK step %d: Ev bit mismatches %llu over %llu sides | QR z mismatches %llu over %llu\n", st, (unsigned long long)g_jchk_bad.load(), (unsigned long long)g_jchk_n.load(), (unsigned long long)g_jchk_qbad.load(), (unsigned long long)g_jchk_qn.load());
    if (g_jprof) { static mach_timebase_info_data_t tb = [] { mach_timebase_info_data_t t; mach_timebase_info(&t); return t; }(); const double f = (double)tb.numer / tb.denom * 1e-6;
      printf("J_PROF step %d: outputs lane-ms res48 %.2f score48 %.2f over %llu cands | eval %.2f qr %.2f | mean nodes %.1f Q %.1f\n", st, g_jp_res.exchange(0) * f, g_jp_sc.exchange(0) * f, (unsigned long long)g_jp_n.exchange(0), g_jp_ev.exchange(0) * f, g_jp_qr.exchange(0) * f, (double)g_jp_nn.load() / std::max<uint64_t>(1, g_jp_sd.load()), (double)g_jp_qq.load() / std::max<uint64_t>(1, g_jp_sd.load())); g_jp_nn = 0; g_jp_qq = 0; g_jp_sd = 0;  if (g_yfxfast) printf("YFX_FAST step %d: cache hit %llu miss %llu | check %llu mismatch %llu | score check %llu mismatch %llu | z check %llu mismatch %llu\n", st, (unsigned long long)g_yf_hit.exchange(0), (unsigned long long)g_yf_miss.exchange(0), (unsigned long long)g_yf_chk.load(), (unsigned long long)g_yf_bad.load(), (unsigned long long)g_yf_schk.load(), (unsigned long long)g_yf_sbad.load(), (unsigned long long)g_yf_zchk.load(), (unsigned long long)g_yf_zbad.load()); }
    const double timed = T.timed();
    printf("step %2d | children %d admitted %6.0f shortfall %5.0f | rej within %5.0f vs pop %5.0f | uniq trees %6.0f (fits x %5.0f y %5.0f)\n",
           st, NP, T.admitted, T.shortfall, T.rej_within, T.rej_pop, T.nuniq, T.nfit0, T.nfit1);
    h301_mark("after select: memo_line/prints");
    if (g_h648) printf("H648 step %2d: MID %.1f GPU idle %.1f -> s %.3f\n", st, g_h648_mid, g_h648_idle, g_h648_s);
    memo_line(E, st, timed); { static const int mfp = getenv("MEM_FP") ? atoi(getenv("MEM_FP")) : 0; if (mfp) { task_vm_info_data_t vi; mach_msg_type_number_t cn = TASK_VM_INFO_COUNT; if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vi, &cn) == KERN_SUCCESS) printf("MEMFP step %2d | phys %.0f MB | phys peak %.0f MB | compressed %.0f MB\n", st, mb_of((double)vi.phys_footprint), mb_of((double)vi.ledger_phys_footprint_peak), mb_of((double)vi.compressed)); } } h141_step_print(st); h150_step_print(st); h156_step_print(E, st); if (memdiag::tab) { memdiag::peak_report(st); memdiag::report(st, "step end"); memdiag::metal_report(st, mb_of((double)[E.dev currentAllocatedSize])); fprintf(stderr, "MEMDIAG st=%d h899 pe heap %.1f MB\n", st, mb_of((double)g_h899_pe.hsize)); fprintf(stderr, "MEMDIAG st=%d mmap NF arena used %.1f MB resident %.1f MB of reserve 1024\n", st, mb_of((double)E.SH->bump.load()), memdiag::region_res_mb(E.SH->ar, (size_t)E.SH->bump.load())); } if ((g_h917 & 8) && h917_at(st)) h917_report(st);
    h301_mark("J: print2 block");
    {
      const Times Tp = T; const double tsum = timed;
      const double tr_bce = bc >= 0 ? (double)err[bc] : (double)INFINITY; const int tr_bcx = bc >= 0 ? (int)evo::cplx(*chp[bc]) : -1;
      const double tr_bpe = g_h652 ? 0.0 : (double)XC(S, pidx[bp]).err; const int tr_bpx = g_h652 ? -1 : (int)XC(S, pidx[bp]).cx;
      auto print2 = [&popset, Tp, J, tsum, tr_bce, tr_bcx, tr_bpe, tr_bpx, st, NP, traj]() {
        h896_join();
        printf("        select skips %4.0f (migrant %4.0f) | survivors that are new children %5.0f | pass-B fill %4.0f | pop distinct keys %6zu, NF pairs %6zu, NF trees %6zu\n",
               Tp.sel_skip, Tp.sel_skip_mig, Tp.surv_child, Tp.xfill, popset.size(), J->npair, J->ntree);
        if (traj) fprintf(traj, "%d,%d,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%zu,%zu,%zu,%.9g,%d,%.9g,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f\n",
                          st, NP, Tp.admitted, Tp.rej_within, Tp.rej_pop, Tp.shortfall, Tp.nuniq, Tp.nfit0, Tp.nfit1, Tp.sel_skip, Tp.sel_skip_mig, Tp.surv_child, Tp.xfill,
                          popset.size(), J->npair, J->ntree,
                          tr_bce, tr_bcx, tr_bpe, tr_bpx,
                          Tp.breed, Tp.nf, Tp.dedup, Tp.comp, Tp.scan_enc + Tp.scan_wall, Tp.gn_enc + Tp.gn_wall, Tp.glue, Tp.select, tsum, NP / (tsum * 1e-3));
      };
      if ((g_h117 & 8) || (g_h131 & 4)) pend_print = print2; else print2();
    }
    h301_mark("J: best child print");
    if (!g_h652) printf("        best child err %.6g (cx %d) | best pop err %.6g (cx %d) | NF memo hit %.3f | alive S %.0f G %.0f | dead nf %.0f cmp %.0f scan %.0f gn %.0f\n",
           bc >= 0 ? err[bc] : INFINITY, bc >= 0 ? (int)evo::cplx(*chp[bc]) : -1, XC(S, pidx[bp]).err, (int)XC(S, pidx[bp]).cx,
           T.nf_hit / std::max(1.0, T.nf_hit + T.nf_miss), T.s_alive, T.g_alive, T.nf_fail, T.cmp_fail, T.scan_fail, T.gn_drop + T.gn_kmiss);
    if (g_h652) printf("        NF memo hit %.3f | alive S %.0f G %.0f | dead nf %.0f cmp %.0f scan %.0f gn %.0f\n", T.nf_hit / std::max(1.0, T.nf_hit + T.nf_miss), T.s_alive, T.g_alive, T.nf_fail, T.cmp_fail, T.scan_fail, T.gn_drop + T.gn_kmiss);
    printf("        ms: breed %.1f | NF %.1f (compile_nf %.1f) | dedup %.1f | compile %.1f (compile_tree %.1f) | scan enc %.1f wall %.1f gpu %.1f | GN enc %.1f wall %.1f gpu %.1f | glue %.1f | select %.1f | timed %.1f | children/s %.0f\n",
           T.breed, T.nf, T.nf_only, T.dedup, T.comp, T.cmp_only, T.scan_enc, T.scan_wall, T.scan_gpu, T.gn_enc, T.gn_wall, T.gn_gpu, T.glue, T.select, timed, NP / (timed * 1e-3));
    if (g_h890) { const double kn = (double)g_h890_key_ns.exchange(0) * 1e-6, kb = (double)g_h890_key_bytes.exchange(0) / 1048576.0;
      printf("        H890=%d: key handoff loop lane-ms %.2f (%.1f MB of keys %s)\n", g_h890, kn, kb, (g_h890 & 1) ? "referenced" : "copied"); }
    printf("        pipe (K=%d): step %.1f | head %.1f | GPU busy (union) %.1f = %.0f%% of step [CB spans sum %.1f: scan %.1f + GN %.1f] | main blocked on GPU %.1f (scan %.1f, GN %.1f) | CPU path %.1f | rebuild %.1f\n",
           g_chunks, timed, T.head, T.gpu_union, 100.0 * T.gpu_union / std::max(1e-9, timed), T.scan_gpu + T.gn_gpu, T.scan_gpu, T.gn_gpu,
           T.w_s + T.w_g, T.w_s, T.w_g, T.pipe_wall - (T.w_s + T.w_g), T.rebuild);
    printf("        I14 post wall %.2f ms, EB scratch %.1f MB\n", T.p_wall, T.p_ebmb);
    if (g_h898 & 2) { printf("H898: polish prep in-loop %.3f ms, in FIN window %.3f ms | arrays %.3f bins %.3f (of which grouping %.3f) (wall of calls, all sites)\n", g_h898_loop_ms, g_h898_fin_ms, g_h898_arr_ms, g_h898_bin_ms, g_h898_grp_ms); g_h898_loop_ms = g_h898_fin_ms = g_h898_arr_ms = g_h898_bin_ms = g_h898_grp_ms = 0; }
    printf("H92 step %2d: select %.2f | dedup %.2f | uq %.2f ms (bit1 %s, bit2 %s)\n", st, T.select, T.dedup, T.uq, (g_h92 & 1) ? "move" : "copy", (g_h92 & 2) ? "sharded" : "serial");
    printf("        I14 post: cand %.0f | enc+prep %.2f | grid gpu %.3f | band gpu %.3f | polish dispatch (itmax 3, incl. re-solve) gpu %.3f | re-solve alone gpu %.3f (S14_TIME) | cpu %.2f | accepted x %.0f y %.0f | viol(joint sum final>stage8) %.0f side-up %.0f fallback %.0f | dropped %.0f\n",
           T.p_cand, T.p_enc, T.p_grid, T.p_band, T.p_pb, T.p_pa, T.p_cpu, T.p_acc0, T.p_acc1, T.p_viol, T.p_side_up, T.p_fb, T.p_ndrop);
    if (g_chk_fin) printf("        I14 final-rss gate (double recompute at final theta, band t): n %.0f, max rel diff %.3g, mismatches (rel > 1e-3 and abs > 1e-5 yy) %.0f; of which in the rank-margin bucket (min keep ratio < 10 tau): %.0f of %.0f\n", T.p_fin_n, T.p_fin_max, T.p_fin_bad, T.p_fin_lobad, T.p_fin_lo);
    if (T.p_chk_n > 0) printf("        I14 band check: n %.0f, max rel obj err %.3g, non-tie path mismatches %.0f, repaired-monotone violations %.0f, corridor violations %.0f\n",
           T.p_chk_n, T.p_chk_objmax, T.p_chk_mis, T.p_chk_mono, T.p_chk_corr);
    printf("        check (%d x, %d y): GN final worse than CPU ref (1e-3) %ld, GPU better %ld, GPU/ref median ratio %.6f, check %.0f ms\n",
           ncheck / 2, ncheck / 2, miss[0] + miss[1], better[0] + better[1], ratio_med, T.check);
    if (g_h864 & 1) printf("H864 step %2d: p1 %.3f find %.3f pref %.3f ctor %.3f comp %.3f ids %.3f ct %.3f | serial p1+pref+ctor+ids+ct %.3f | uq %.3f comp(total) %.3f | H864=%d\n",
                           st, T.u_p1, T.u_find, T.u_pref, T.u_ctor, T.u_comp, T.u_ids, T.u_ct, T.u_p1 + T.u_pref + T.u_ctor + T.u_ids + T.u_ct, T.uq, T.comp, g_h864);
    if (g_h879 & 1) {
      printf("H879 step %2d: runs %.0f recs %.0f | run wall %.3f ms (lane-ms %.3f) | busy(compile) %.3f cpu %.3f | lane use busy %.3f cpu %.3f | start latency sum %.3f | tail idle sum %.3f | order build %.3f | H879=%d\n",
             st, T.c_runs, T.c_recs, T.c_wall, T.c_wall * (double)g_h879_lanes, T.c_busy, T.c_cpu, T.c_busy / std::max(1e-9, T.c_wall * (double)g_h879_lanes), T.c_cpu / std::max(1e-9, T.c_wall * (double)g_h879_lanes), T.c_lat, T.c_tail, T.c_ord, g_h879);
      { std::lock_guard<std::mutex> lg(g_h879_m); printf("H879 step %2d per task busy ms/cpu ms/claims:", st);
        for (int w = 0; w < 16 && w < 64; w++) if (g_h879_wn[w] || g_h879_wc[w] > 0) printf(" [%d %.2f %.2f %ld]", w, g_h879_wb[w], g_h879_wc[w], g_h879_wn[w]);
        printf("\n"); for (int w = 0; w < 64; w++) { g_h879_wb[w] = g_h879_wc[w] = 0; g_h879_wn[w] = 0; } }
    }
    if (g_h93_print) printf("H93 step %2d: scan enc %.3f [alloc %.3f pre %.3f sort %.3f bins %.3f fill %.3f] scan cmd %.3f [pipe %.3f] | GN enc %.3f [alloc %.3f pre %.3f sort %.3f plan %.3f fill %.3f] GN cmd %.3f [pipe %.3f] | H93=%d\n",
           st, T.scan_enc, T.h_s_alloc, T.h_s_pre, T.h_s_sort, T.h_s_bins, T.h_s_fill, T.h_s_cmd, T.h_s_pipe,
           T.gn_enc, T.h_g_alloc, T.h_g_pre, T.h_g_sort, T.h_g_plan, T.h_g_fill, T.h_g_cmd, T.h_g_pipe, g_h93 ? 1 : 0);
    if (g_h99) {
      const long dl = G99.lookups - h99_l0, dh = G99.hits - h99_h0;
      printf("H99 step %2d | lookups %6ld hits %6ld | hit rate %.4f | hits / unique trees %.4f | inserted %6ld evicted %6ld refused %ld dups %ld | entries %u of %zu | memo %.1f MB\n",
             st, dl, dh, dl > 0 ? (double)dh / (double)dl : 0.0, T.nuniq > 0 ? (double)dh / T.nuniq : 0.0,
             G99.inserts - h99_i0, G99.evicts - h99_e0, G99.refused - h99_r0, G99.dups - h99_d0, (unsigned)G99.nlive, G99.cap, G99.bytes() / 1048576.0);
      if (st >= steady_from) { h99_sl += dl; h99_sh += dh; h99_su += T.nuniq; } }
    h301_mark("J: diag prints");
    if (st == steps && pop_thr.joinable()) { auto tj = clk::now(); pop_thr.join(); T.rebuild += ms_since(tj); }
    if (st == steps && h4_thr.joinable()) { auto tj = clk::now(); h4_thr.join(); T.h131_h4wait += ms_since(tj); }
    if (st == steps && Jh4) { T.h131_h4 += Jh4->dump_ms + Jh4->sets_ms; Jh4.reset(); }
    if (st == steps && pend_print) { pend_print(); pend_print = nullptr; }
    if (st == steps && h652_thr.joinable()) { h652_thr.join(); }
    if (g_h651 || g_h652) printf("H651 step %2d: memo join %.3f ms | H652 join %.3f ms\n", st, h651_join_ms, h652_join_ms);
    h651_join_ms = h652_join_ms = 0;
    T.h131_iter = ms_since(t_iter0);
    addT(TT, T); fflush(stdout);
    if (h876_mask() & 4) printf("H876 mem step %2d: rss %.1f MB phys %.1f MB\n", st, rss_mb_now(), phys_fp_mb());
    if (h178_on() || h178_pre()) printf("H178 step %2d: PSOs created %ld | archive hits %ld misses %ld | PSO creation wall so far %.1f ms\n", st, g_h178_created - h178_c0, g_h178_hits, g_h178_miss, g_h178_ms);
    if (g_h650d && (h650_is("gncpu") || h650_is("polcpu"))) { printf("H650 step %d: injected fit core-ms %.1f (%llu items x %.0f ns)\n", st, g_h650_items.load() * g_h650ns * 1e-6, (unsigned long long)g_h650_items.load(), g_h650ns); g_h650_items = 0; }
    if (g_h891) { printf("H891 step %2d: cpu lane-ms %.2f (%llu batches, %.1f us avg) | gpu CBs %.0f requested %.3f ms calibrated %.3f ms | head %.3f ms tail %.3f ms\n", st,
        g_h891_cpu_ns.load() * 1e-6, (unsigned long long)g_h891_cpu_n.load(), g_h891_cpu_n.load() ? g_h891_cpu_ns.load() * 1e-3 / g_h891_cpu_n.load() : 0.0, g_h891_gpu_n, g_h891_gpu_req_ms, g_h891_gpu_cal_ms, g_h891_head_ms, g_h891_tail_ms);
      g_h891_cpu_ns = 0; g_h891_cpu_n = 0; g_h891_gpu_n = g_h891_gpu_req_ms = g_h891_gpu_cal_ms = 0; g_h891_head_ms = g_h891_tail_ms = 0; }
    printf("ITER WALL step %2d: %.1f ms (H131 whole-iteration wall; timed %.1f) thermal %ld\n", st, T.h131_iter, T.timed(), (long)g_thermal_state);
    h301_mark("J: iter wall/addT");
    const bool kdef = kpop;
    if (kdef) {
      k_flush();
      struct KPK { std::vector<uint64_t> pk; std::vector<uint8_t> pres; std::vector<float> jmin; int nbk = 0, next = 0; };
      auto K = std::make_shared<KPK>(); const int JB = 2048; K->nbk = (NP + JB - 1) / JB; K->pk.resize((size_t)NP); K->pres.assign((size_t)NP, 1); K->jmin.assign((size_t)K->nbk, INFINITY);


      const size_t NS3 = (size_t)3 * NP; if (g_kskey.size() != NS3) { g_kskey.assign(NS3, 0); g_kvalid.assign(NS3, 0); g_kcs_prev.clear(); }
      { std::vector<uint8_t> dirty(NS3, 0); for (uint32_t x : kcsN) dirty[x] = 1; if (!(g_h826 && g_pinfix && g_sf_pin)) for (uint32_t x : g_kcs_prev) dirty[x] = 1;
        const int nb2 = (NP + 4095) / 4096; std::atomic<long> nmis{0}, nrec{0};
        pool.run(nb2, [&](int b) { long rec = 0, mis = 0; for (int i = b * 4096; i < std::min(NP, (b + 1) * 4096); i++) { const uint32_t sl = pidx[i]; uint64_t k;
            if (dirty[sl] || !g_kvalid[sl]) { k = fitx::gkey(XC(S, sl).g); g_kskey[sl] = k; g_kvalid[sl] = 1; rec++; }
            else { k = g_kskey[sl]; if (g_kpopkey == 2) { const uint64_t k2 = fitx::gkey(XC(S, sl).g); if (k2 != k) { mis++; k = k2; g_kskey[sl] = k2; } } }
            K->pk[(size_t)i] = k; }
          nrec.fetch_add(rec); nmis.fetch_add(mis); });
        if (g_kpopkey == 2 || g_kverb) printf("K_POPKEY step %d: keys recomputed %ld of %d | verify mismatches %ld\n", st, nrec.load(), NP, nmis.load()); }
      g_kcs_prev = std::move(kcsN);
      g_kfill = [K, &pool, &S, &pidx, NP, JB]() -> bool {
        if (K->next >= K->nbk) return false;
        const int b0 = K->next, b1 = std::min(K->nbk, b0 + 8); K->next = b1;
        pool.run(b1 - b0, [&](int q) { const int b = b0 + q, lo = b * JB, hi = std::min(NP, lo + JB); float m = INFINITY;
          for (int i = lo; i < hi; i++) { const float e_ = XC(S, pidx[i]).err; if (e_ < m) m = e_; }
          K->jmin[(size_t)b] = m;
          fitx::present_mask(K->pk, (size_t)lo, (size_t)hi, K->pres.data() + lo); });
        return K->next < K->nbk;
      };
      g_kfin = [K, NP, &S, &pidx]() {
        const auto fx0 = clk::now(); float bpe = INFINITY;
        for (int b = 0; b < K->nbk; b++) if (K->jmin[(size_t)b] < bpe) bpe = K->jmin[(size_t)b];
        { static const int rd = getenv("REC_DIAG") ? atoi(getenv("REC_DIAG")) : 0; if (rd) { long np_ = 0; for (int i = 0; i < NP; i++) np_ += K->pres[(size_t)i]; std::unordered_set<uint64_t> u(K->pk.begin(), K->pk.end());
            int bi = -1; float be = INFINITY; long nin = 0; for (int i = 0; i < NP; i++) { const float e_ = XC(S, pidx[i]).err; if (e_ < be) { be = e_; bi = i; } }
            for (int i = 0; i < NP; i++) if (XC(S, pidx[i]).err <= be + 80.f) nin++;
            printf("REC_DIAG pop: NP %d with-record %ld distinct-keys %zu within-80-of-best %ld best err %g best-has-record %d table %zu\n", NP, np_, u.size(), nin, (double)be, bi >= 0 ? (int)K->pres[(size_t)bi] : -1, fitx::st().all.size()); } }
        fitx::g_gate = std::isfinite(bpe) ? bpe + (float)fitx::gate_margin() : INFINITY;
        std::vector<uint64_t>& pk = K->pk; size_t w = 0; for (int i = 0; i < NP; i++) if (K->pres[(size_t)i]) pk[w++] = pk[(size_t)i]; pk.resize(w);
        fitx::retain(pk); fitx::g_ms += ms_since(fx0);
      };
    }
    if (fitx::on() && !kdef) { const auto fx0 = clk::now(); std::vector<uint64_t> pk((size_t)NP); float bpe = INFINITY;
    static const int jfp = getenv("J_FITPAR") ? atoi(getenv("J_FITPAR")) : 1;
    const int JB = 2048, nbk = (NP + JB - 1) / JB; std::vector<float> jmin((size_t)nbk, INFINITY);
    if (jfp && NP >= 4096) pool.run(nbk, [&](int b) { const int lo = b * JB, hi = std::min(NP, lo + JB); float m = INFINITY; for (int i = lo; i < hi; i++) { const float e_ = XC(S, pidx[i]).err; if (e_ < m) m = e_; } jmin[(size_t)b] = m; });
    if (jfp && NP >= 4096) { for (int b = 0; b < nbk; b++) if (jmin[(size_t)b] < bpe) bpe = jmin[(size_t)b]; }
    else for (int i = 0; i < NP; i++) { const float e_ = XC(S, pidx[i]).err; if (e_ < bpe) bpe = e_; } fitx::g_gate = std::isfinite(bpe) ? bpe + (float)fitx::gate_margin() : INFINITY;
    const double jb = ms_since(fx0);

    std::vector<uint8_t> jpres;
    if (jfp && NP >= 4096) { jpres.assign((size_t)NP, 1); pool.run(nbk, [&](int b) { const int lo = b * JB, hi = std::min(NP, lo + JB); for (int i = lo; i < hi; i++) pk[(size_t)i] = fitx::gkey(XC(S, pidx[i]).g); fitx::present_mask(pk, (size_t)lo, (size_t)hi, jpres.data() + lo); }); }
    else for (int i = 0; i < NP; i++) pk[(size_t)i] = fitx::gkey(XC(S, pidx[i]).g);
    const double jc = ms_since(fx0);
    if (!jpres.empty()) { size_t w = 0; for (int i = 0; i < NP; i++) if (jpres[(size_t)i]) pk[w++] = pk[(size_t)i]; pk.resize(w); }
    fitx::retain(pk); if (g_jprof) printf("J_PROF fitx block: NP %d bpe %.2f gkey %.2f retain %.2f ms\n", NP, jb, jc - jb, ms_since(fx0) - jc); fitx::g_ms += ms_since(fx0); printf("FIT_EXPORT step %d: host ms on the main thread: writeback capture %.1f | fixed-capture+commit %.1f (commit %.1f) | retain %.1f | table %zu\n", st, fitx::g_ms_wb, fitx::g_ms, fitx::g_ms_commit, fitx::g_ms_retain, fitx::st().all.size()); fitx::g_ms = fitx::g_ms_wb = fitx::g_ms_commit = fitx::g_ms_retain = 0; }
    if (kdef) { printf("FIT_EXPORT step %d: host ms on the main thread: writeback capture %.1f | fixed-capture+commit %.1f (commit %.1f) | retain %.1f | table %zu\n", st, fitx::g_ms_wb, fitx::g_ms, fitx::g_ms_commit, fitx::g_ms_retain, fitx::st().all.size()); fitx::g_ms = fitx::g_ms_wb = fitx::g_ms_commit = fitx::g_ms_retain = 0; }
    h301_mark("J: fitx");
    int ui_pick = -1;
    {

      static const bool usp = !getenv("UI_STEP_PICK") || atoi(getenv("UI_STEP_PICK")) != 0;
      static const bool usl = getenv("UI_STEP_LOG") && atoi(getenv("UI_STEP_LOG")) != 0;
      const bool upk_ok = g_sicis_en != nullptr && fitx::on() && !(getenv("FINAL_PICK_NML") && atoi(getenv("FINAL_PICK_NML")) == 0);
      if ((usp && uib::g_ring || usl) && upk_ok) {
        const double t_end = ui_ms_since_start(); const auto up0 = clk::now();
        auto fexprU = [&](int i) { std::string o; fitx::Rec fr; if (i < 0 || !fit_get(E, XC(S, pidx[i]).g, fr)) return std::string("?"); for (int sd = 0; sd < 2; sd++) { uib::Tree T; std::memset(&T, 0, sizeof T);
            o += sd ? " | y=" : "x="; if (fitx::emit_tree(fr, sd, T)) { int pos = T.len - 1; o += fp::expr_of(T.op, T.c, T.len, pos); } else o += "?"; } return o; };
        static const double w = getenv("SIMPLICITY") ? atof(getenv("SIMPLICITY")) : 1.0;
        const bool t13_all = fp::t13_bit(0);
        static const bool trig_rule = !getenv("TRIG_RULE") || atoi(getenv("TRIG_RULE")) != 0;
        static const bool geo_pick = !getenv("GEO_PICK") || atoi(getenv("GEO_PICK")) != 0;
        static const bool kfit = getenv("FP_KFIT") && atoi(getenv("FP_KFIT")) != 0;


        const double ws = (w > 0 ? w : 1.0);
        std::vector<double> pe(NP), lb(NP); std::vector<int> pc(NP), pnx(NP), pny(NP), pkc(NP, 0); std::vector<fp::nml::FullTerms> fulls(NP);
        const int UB = 4096, unb = (NP + UB - 1) / UB;
        pool.run(unb, [&](int b) { const int lo = b * UB, hi = std::min(NP, lo + UB);
          for (int i = lo; i < hi; i++) { const evo::Ind& PI = XC(S, pidx[i]); pe[i] = PI.err; pc[i] = sp_show_cx(PI); pnx[i] = PI.g.len[0]; pny[i] = PI.g.len[1];
            lb[i] = std::isfinite(PI.err) ? (double)PI.err + ws * (fp::nml::struct_code_canon(pnx[i]) + fp::nml::struct_code_canon(pny[i])) : INFINITY; } });
        std::vector<int> oa; oa.reserve(NP); for (int i = 0; i < NP; i++) if (std::isfinite(lb[i])) oa.push_back(i);
        auto prep = [&](int i) { const evo::Ind& PI = XC(S, pidx[i]); int q_[2], n_[2];
          if (kfit) { if (fitx::qnth(PI.g, q_, n_, g_freet ? 1 : 0)) pkc[i] = q_[0] + q_[1] + n_[0] + n_[1]; else pkc[i] = -1; }
          else pkc[i] = fp::nml::const_leaves(PI.g.op[0], PI.g.len[0]) + fp::nml::const_leaves(PI.g.op[1], PI.g.len[1]);
          fulls[i] = fp::fallback_terms(pkc[i] < 0 ? 100000 : pkc[i], 48, pnx[i], pny[i]); if (pkc[i] < 0) pe[i] = INFINITY;
          if (w > 0 && w != 1.0) { fulls[i].paramTerm *= w; fulls[i].structTerm *= w; fulls[i].fisherTerm *= w; fulls[i].rangeTerm *= w; fulls[i].total *= w; } };
        auto tot0 = [&](int i) { return (double)pe[i] + fulls[i].total; };
        const double t_prep = std::chrono::duration<double, std::milli>(clk::now() - up0).count();
        std::vector<char> seen(NP, 0); long nvis = 0;
        {   std::priority_queue<double> hp; size_t done = 0, K = 64; bool stop = false;
          static const double sw = getenv("SIMPLICITY") ? std::max(atof(getenv("SIMPLICITY")), 1e-300) : 1.0;
          static const bool fa = !getenv("FREET_AFFINE") || atoi(getenv("FREET_AFFINE")) != 0;
          static const bool al = !getenv("ALIAS_RULE") || atoi(getenv("ALIAS_RULE")) != 0;
          double thrB = INFINITY;
          auto evalm = [&](int i) { prep(i); if (!std::isfinite(pe[i])) return; if (tot0(i) >= thrB) return; seen[i] = 1;
            if (!trig_rule) return;
            fitx::Rec fr; if (!fitx::find(XC(S, pidx[i]).g, fr)) return;
            if (fa && geop::freet_nonaffine(fr)) { pe[i] = INFINITY; return; }
            if (al && geop::aliased(fr)) { pe[i] = INFINITY; return; }
            const int ex = geop::trig_extra_k(fr);
            if (ex > 0) { const double dp = 0.5 * ex * std::log(48.0 / (2.0 * fp::nml::kPi)) * sw; fulls[i].paramTerm += dp; fulls[i].total += dp; } };
          while (!stop && done < oa.size()) {
            const size_t kk = std::min(oa.size(), done + K);
            std::partial_sort(oa.begin() + done, oa.begin() + kk, oa.end(), [&](int a, int b) { return lb[a] < lb[b]; });
            const double thr = (t13_all && hp.size() >= 64) ? hp.top() : INFINITY; thrB = thr;
            size_t ke = done; while (ke < kk && lb[oa[ke]] < thr) ke++;
            const int MB = 8, nb2 = (int)((ke - done + MB - 1) / MB);
            pool.run(nb2, [&](int b) { const size_t lo = done + (size_t)b * MB, hi = std::min(ke, lo + MB); for (size_t r = lo; r < hi; r++) evalm(oa[r]); });
            nvis += (long)(ke - done);
            for (size_t r = done; r < ke; r++) { const int i = oa[r]; if (t13_all && hp.size() >= 64 && lb[i] >= hp.top()) { stop = true; break; }
              if (!seen[i] || !std::isfinite(pe[i])) continue; hp.push(tot0(i)); if (hp.size() > 64) hp.pop(); }
            if (ke < kk) stop = true;
            done = kk; K *= 2; }
        }
        for (int i = 0; i < NP; i++) if (!seen[i]) pe[i] = INFINITY;
        const double t_walk = std::chrono::duration<double, std::milli>(clk::now() - up0).count();
        const double t_seen = std::chrono::duration<double, std::milli>(clk::now() - up0).count();
        int pk = fp::pick_index_nml_full(pe, pc, fulls, nullptr, t13_all);
        const double t_pick = std::chrono::duration<double, std::milli>(clk::now() - up0).count();
        if (geo_pick && pk >= 0) {
          std::vector<int> ord; for (int i = 0; i < NP; i++) if (std::isfinite(pe[i])) ord.push_back(i);
          static const bool pick_rescore_u = !getenv("PICK_RESCORE") || atoi(getenv("PICK_RESCORE")) != 0;   
          const int T64 = pick_rescore_u ? (int)ord.size() : std::min<int>(64, (int)ord.size());
          std::partial_sort(ord.begin(), ord.begin() + T64, ord.end(), [&](int a, int b) { return tot0(a) < tot0(b); });
          std::vector<int> sub; std::vector<double> se; std::vector<int> sc2; std::vector<fp::nml::FullTerms> sf;
          std::vector<double> gg((size_t)T64, NAN);
          pool.run(T64, [&](int r) { const int i = ord[r]; fitx::Rec fr; if (!fitx::find(XC(S, pidx[i]).g, fr)) return;
            const float* Y0 = fr.kind == 1 ? g_fit_yv : g_yfx ? g_yfx_f[0] : fitx::g_target[0]; const float* Y1 = fr.kind == 1 ? g_fit_yv + 48 : g_yfx ? g_yfx_f[1] : fitx::g_target[1];
            static const bool geo_fixed = getenv("GEO_FIXED") && atoi(getenv("GEO_FIXED")) != 0;
            gg[(size_t)r] = (fr.kind != 1 && !geo_fixed) ? (double)pe[i] : geop::geo_nll(fr, Y0, Y1, *g_sicis_en, g_nres ? (g_yfx ? g_yfx_nt : g_nres_n) : nullptr); });
          for (int r = 0; r < T64; r++) { const int i = ord[r]; const double g = gg[(size_t)r]; if (!std::isfinite(g)) continue;
            sub.push_back(i); se.push_back(g); sc2.push_back(pc[i]); sf.push_back(fulls[i]); }
          if (!sub.empty()) { const int q = fp::pick_index_nml_full(se, sc2, sf, nullptr, t13_all); if (q >= 0) { pe[sub[q]] = se[q]; pk = sub[q]; } }
        }
        if (usp) ui_pick = pk;
        if (usl) { printf("UI_STEP_MS step %d prep %.2f walk %.2f (visited %ld of %zu, charged %ld) seen-loop %.2f pick %.2f geo+rest %.2f total %.2f\n", st, t_prep, t_walk - t_prep, nvis, oa.size(), (long)std::count(seen.begin(), seen.end(), 1), t_seen - t_walk, t_pick - t_seen, std::chrono::duration<double, std::milli>(clk::now() - up0).count() - t_pick, std::chrono::duration<double, std::milli>(clk::now() - up0).count()); int rb = 0; for (int i = 1; i < NP; i++) if (XC(S, pidx[i]).err < XC(S, pidx[rb]).err) rb = i;
          printf("UI_STEP step %d t_ms %.1f pick_ms %.2f | raw best-err: err %.6g cx %d [%s] | nml pick: err %.6g cx %d total %.6g k %d [%s]\n", st, t_end, std::chrono::duration<double, std::milli>(clk::now() - up0).count(),
                 (double)XC(S, pidx[rb]).err, (int)XC(S, pidx[rb]).cx, fexprU(rb).c_str(), pk >= 0 ? pe[pk] : NAN, pk >= 0 ? pc[pk] : -1, pk >= 0 ? tot0(pk) : NAN, pk >= 0 ? pkc[pk] : -1, fexprU(pk).c_str()); }
      }
    }
    if (uib::g_ring) {
      const int np = std::min(NP, uib::PAIRS);
      uint32_t ord[uib::PAIRS]; float sc[uib::PAIRS]; uint16_t cxs[uib::PAIRS]; int nk = 0;
      for (int i = 0; i < NP; i++) {
        const float e = XC(S, pidx[i]).err;
        if (nk == np && e >= XC(S, ord[np - 1]).err) continue;
        int j = nk < np ? nk++ : np - 1;
        while (j > 0 && XC(S, ord[j - 1]).err > e) { ord[j] = ord[j - 1]; j--; }
        ord[j] = pidx[i];
      }
      if (ui_pick >= 0) {
        const uint32_t want = pidx[ui_pick]; int at = -1; for (int i = 0; i < nk; i++) if (ord[i] == want) { at = i; break; }
        if (at < 0) { at = nk < np ? nk++ : np - 1; ord[at] = want; }
        for (int i = at; i > 0; i--) ord[i] = ord[i - 1]; ord[0] = want; }
      for (int i = 0; i < np; i++) { const evo::Ind& PI = XC(S, ord[i]); sc[i] = PI.err; cxs[i] = (uint16_t)sp_show_cx(PI); }
      auto fillp = [&](int p, uib::Pair& P) { const evo::Ind& PI = XC(S, ord[p]); fitx::Rec fr; if (fit_get(E, PI.g, fr)) fitx::fill_pair(fr, P, PI.g); else uib::fill_from_cand(PI.g, P); };
      uib::push_burst(st, T.h131_iter, (int64_t)NP, np, sc, cxs, fillp, 0, 0.f, 1.f);
    }
    h657_step(st);
    if (uib::g_ring && uib::g_budget_s > 0 && (h301_now() - uib::g_t0) >= uib::g_budget_s) { printf("UIB budget reached after step %d\n", st); uib::g_stop = 1; }
    if (g_h874 & 16) { printf("H874 step %2d: minflt step %ld | evaluate_v7 %ld\n", st, g_h874_mf0 < 0 ? 0L : h874_mf() - g_h874_mf0, mf874b - mf874a); g_h874_mf0 = h874_mf(); }
    if (g_h854 & 1) {
      struct rusage r; getrusage(RUSAGE_SELF, &r);
      const double o = h301_now() - T.h131_iter * 1e-3;
      if (g_h854_ru0_ok) printf("H854 step %2d | minflt %ld majflt %ld | iter %.1f ms\n", st, r.ru_minflt - g_h854_ru0.ru_minflt, r.ru_majflt - g_h854_ru0.ru_majflt, T.h131_iter);
      g_h854_ru0 = r; g_h854_ru0_ok = true;
      if (g_h886 & 2) printf("H886 step %2d | boost qos seen %d\n", st, g_h886_qos_seen.load());
      if (g_h854_joined) {
        const H854Rec& q = (g_h886 & 1) ? g_h886_rec : g_h854_rec;
        printf("H854 step %2d | boundary(for this step): %s launch %.1f start %.1f end %.1f ms | hist %.1f nfcompact %.1f cccompact %.1f | join wait %.1f done %.1f (rel. iter start)\n",
               st, q.ran ? "RETAIN" : "no-op", (q.launch - o) * 1e3, (q.start - o) * 1e3, (q.end - o) * 1e3, q.hist, q.nfc, q.ccc, g_h854_join_wait, (g_h854_join_done - o) * 1e3);
        g_h854_joined = false;
      }
    }
    if (g_h770) h770_step_print(st);
    if (g_h771) h771_step_print(st);
    if (g_h413dbg && g_h242) { h242_drain(); printf("H413DBG step %2d: hits %ld (step %ld) hitkeys %016llx inserted rows %016llx (%ld)\n", st, G242.hits, G242.hits - g_h413_h0, (unsigned long long)g_h413_hk, (unsigned long long)g_h413_ik, G242.inserts); g_h413_h0 = G242.hits;
      printf("H413DBG fields step %2d: key %016llx meta %016llx rfin %016llx smin/sarg %016llx row %016llx steps %016llx refused %ld tail %ld la %ld lt %ld\n", st, (unsigned long long)g_h413_f[0], (unsigned long long)g_h413_f[1], (unsigned long long)g_h413_f[2], (unsigned long long)g_h413_f[3], (unsigned long long)g_h413_f[4], (unsigned long long)g_h413_f[5], G242.refused_wide, g_h413_tail.load(), g_h413_la.load(), g_h413_lt.load()); }
    if (st >= steady_from) { addT(TS, T); nsteady++; }
    h301_mark("loop tail (StepBatch free)");
    if (g_h301r && !h870_early && getenv("H301RDBG")) { size_t cap = 0, ne = 0; for (const Tree& t : g_h301r_t) { if (t.prog.capacity() || t.rec.capacity()) ne++; cap += t.prog.capacity() + t.sp.capacity() + t.th.capacity() + t.rec.capacity(); }
      printf("H301R step %2d: pool %zu takes %zu | old-pool entries with capacity %zu, capacity elems %zu | new pool %zu\n", st, g_h301r_t.size(), g_h301r_n.load(), ne, cap, SB.UQ.U.size()); }
    if (g_h415) {
      if (g_h301r && !h870_early) { std::swap(g_h301r_t, SB.UQ.U); std::swap(g_h301r_e, SB.UQ.enc); g_h301r_n.store(0); }
      h301_mark("loop tail: deleter join");
      if (!h870_early) {
      if (sb_thr.joinable()) sb_thr.join();
      StepBatch* p = SBp.release(); g_h415_spare = p;
      sb_thr = std::thread([p] { H401S h401s_("StepBatch clear thread"); p->UQ.U.clear(); p->UQ.enc.clear(); p->UQ.hash.clear(); });
      }
    } else {
    if (g_h301r && g_h405 && g_h301d) { std::swap(g_h301r_t, SB.UQ.U); std::swap(g_h301r_e, SB.UQ.enc); g_h301r_n.store(0); }
    else if (g_h301r) { g_h301r_t = std::move(SB.UQ.U); g_h301r_e = std::move(SB.UQ.enc); g_h301r_n.store(0); }
    h301_mark("loop tail: deleter join");
    if (g_h301d) { if (sb_thr.joinable()) sb_thr.join(); h301_mark("loop tail: after deleter"); sb_thr = std::thread([p = SBp.release()] { H401S h401s_("StepBatch delete thread"); delete p; }); }
    }
  }
  if (g_h670) { h670_join(); g_h670_st++; }
  if (g_h670) {
      if (g_h670_have && g_h670_jst != g_h670_st - 1) { g_h670_have = false; printf("H670 apply: stale results (step %d at step %d) dropped\n", g_h670_jst, g_h670_st); }
      if (g_h670_have) { g_h670_have = false; const auto ta = clk::now();
        std::vector<uint8_t> f((size_t)3 * NP, 0); for (int i = 0; i < NP; i++) f[pidx[i]] = 1;
        long nap = 0, nlost = 0;
        for (const auto& pr : g_h670_res) { best_all = std::min(best_all, (double)pr.second);
          const uint32_t slot = g_h670_cs[(size_t)pr.first]; if (!f[slot]) { nlost++; continue; }
          evo::Ind& c = XA(S, slot); c.res[0] = pr.second; c.err = 0.0f + pr.second; if (g_h740) g_h740e[slot].err = c.err; nap++; }
        printf("H670 apply: accepted %zu | survivors updated %ld | not surviving %ld | %.2f ms\n", g_h670_res.size(), nap, nlost, ms_since(ta)); } }
  if (g_h670) { int p = 0; for (int i = 1; i < NP; i++) if (XC(S, pidx[i]).err < XC(S, pidx[p]).err) p = i;
    printf("H670 final: joins %ld wait %.2f ms/join | best pop err %.6g (cx %d)\n", g_h670_n, g_h670_n ? g_h670_wait / g_h670_n : 0.0, XC(S, pidx[p]).err, (int)XC(S, pidx[p]).cx); }
  if (g_h415) { if (sb_thr.joinable()) sb_thr.join(); delete g_h415_spare; g_h415_spare = nullptr; printf("H415: spare-join wait %.2f ms total\n", g_h415_join_ms); }
  if (sb_thr.joinable()) sb_thr.join();
  if (mb_thr.joinable()) mb_thr.join();
  if (g_h139) printf("H139 final: TG-memory violations (threadgroup request > 32768 B on the threadgroup path) %llu\n", (unsigned long long)g_h139_viol);
  if (traj) fclose(traj);
  if (g_h320_f) fclose(g_h320_f);
  if (g_h301 && h301_n) {
    std::vector<std::pair<double, std::string>> v; for (auto& kv : h301_acc) v.push_back({kv.second / h301_n, kv.first});
    std::sort(v.rbegin(), v.rend());
    printf("STEADY H301 (%ld steps): GPU busy %.1f | GPU idle %.1f ms/step, by main-thread phase:", h301_n, h301_busy / h301_n, h301_idle / h301_n);
    for (auto& p : v) printf(" | %s %.1f", p.second.c_str(), p.first);
    printf("\n");
    std::vector<std::pair<double, std::string>> va; for (auto& kv : g_h301_all) va.push_back({kv.second / h301_n, kv.first}); std::sort(va.rbegin(), va.rend());
    printf("STEADY H301ALL (all wall ms/step by phase label; GPU busy or idle):"); for (auto& p : va) if (p.first >= 0.3) printf(" | %s %.1f", p.second.c_str(), p.first); printf("\n");
  }
  if (nsteady) {
    scaleT(TS, 1.0 / nsteady);
    if (nf308::g_h315t && g_h890) printf("H890 publish block (body assembly + publish) whole-run core-ms %.1f\n", nf308::g_t_pubcopy_h890.load() * 1e-6);
    if (nf308::g_h315t) printf("H315T (whole run, summed core-ms): chunk/group total %.1f | nf %.1f | compile %.1f | dag+publish %.1f | classify %.1f | sort %.1f | eval-loop %.1f | stage-out %.1f | rollback (inside nf) %.1f ms, %llu objects\n", nf308::g_t_all.load() * 1e-6, nf308::g_t_nf.load() * 1e-6, nf308::g_t_cmp.load() * 1e-6, nf308::g_t_pub.load() * 1e-6, nf308::g_t_class.load() * 1e-6, nf308::g_t_sort.load() * 1e-6, nf308::g_t_eval.load() * 1e-6, nf308::g_t_out.load() * 1e-6, nf308::g_t_rb.load() * 1e-6, (unsigned long long)nf308::g_n_rb.load());
    if (nf308::g_h315t && nf308::g_w_obj.load()) printf("H315T work proxies (completed tapes, whole run): objects created %.0f, reused %.0f = %.4f | term pairs (nvp) %.0f, reused %.0f = %.4f\n", (double)nf308::g_w_obj.load(), (double)nf308::g_w_objre.load(), (double)nf308::g_w_objre.load() / (double)nf308::g_w_obj.load(), (double)nf308::g_w_vp.load(), (double)nf308::g_w_vpre.load(), (double)nf308::g_w_vpre.load() / (double)nf308::g_w_vp.load());
    if (nf308::g_h315a) { const double mi = (double)nf308::g_h315_miss.load(), ni = (double)nf308::g_h315_inc.load(), nn = (double)nf308::g_h315_nodes.load(), nr = (double)nf308::g_h315_reused.load();
      printf("H315A COUNTERS (whole run): NF misses %.0f | incremental (reused >= 1 node) %.0f = %.4f | nodes %.0f reused %.0f = %.4f | duplicate-became-hit %llu | groups %llu\n", mi, ni, mi ? ni / mi : 0.0, nn, nr, nn ? nr / nn : 0.0, (unsigned long long)nf308::g_h315_stage_hit.load(), (unsigned long long)nf308::g_h315_groups.load()); }
    printf("HEADLINE step ms (H131 iter wall = whole step body, steady steps %d..%d, mean of %ld): %.1f ms per step | timed %.1f | untimed %.1f | H92=%d H117=%d H128=%d H129=%.0f H131=%d H84=%d S12_CHUNKS=%d\n",
           steady_from, steps, nsteady, TS.h131_iter, TS.timed(), TS.h131_iter - TS.timed(), g_h92, g_h117, g_h128, g_h129_cap_mb, g_h131, g_h84, g_chunks);
    printf("STEADY MEAN over steps %d..%d (%ld steps, ms per step): breed %.1f | NF %.1f (compile_nf %.1f) | dedup %.1f | compile %.1f (compile_tree %.1f) | scan enc %.1f wall %.1f gpu %.1f | GN enc %.1f wall %.1f gpu %.1f | glue %.1f | select %.1f | timed %.1f\n",
           steady_from, steps, nsteady, TS.breed, TS.nf, TS.nf_only, TS.dedup, TS.comp, TS.cmp_only, TS.scan_enc, TS.scan_wall, TS.scan_gpu, TS.gn_enc, TS.gn_wall, TS.gn_gpu, TS.glue, TS.select, TS.timed());
    printf("H92 STEADY (steps %d..%d): select %.2f | dedup %.2f | uq %.2f ms per step | H92=%d\n", steady_from, steps, TS.select, TS.dedup, TS.uq, g_h92);
    {
      const double nc = g_h128_prof.n > 0 ? (double)g_h128_prof.n : 1.0;
      printf("H128 PROF (calls %ld, H128=%d): migr %.3f | pool %.3f | keyids %.3f | passA %.3f | passB %.3f | freeF %.3f | copy %.3f | sum %.3f\n",
             g_h128_prof.n, g_h128, g_h128_prof.ms[0] / nc, g_h128_prof.ms[1] / nc, g_h128_prof.ms[2] / nc, g_h128_prof.ms[3] / nc,
             g_h128_prof.ms[4] / nc, g_h128_prof.ms[5] / nc, g_h128_prof.ms[6] / nc,
             (g_h128_prof.ms[0] + g_h128_prof.ms[1] + g_h128_prof.ms[2] + g_h128_prof.ms[3] + g_h128_prof.ms[4] + g_h128_prof.ms[5] + g_h128_prof.ms[6]) / nc);
      h128_print_sub(nc);
    }
    printf("STEADY MEAN v5: enc_tree (zconvS/zgconv per unique tree, inside compile) %.1f ms per step\n", TS.enc_tree);
    printf("STEADY MEAN: children/s %.0f | admitted %.0f | shortfall %.0f | uniq trees %.0f | rej within %.0f | rej vs pop %.0f | select skips %.0f\n",
           NP / (TS.timed() * 1e-3), TS.admitted, TS.shortfall, TS.nuniq, TS.rej_within, TS.rej_pop, TS.sel_skip);
  }
  if (nsteady) printf("STEADY I14 post (ms/step, steps %d..%d; TS holds per-step means): post wall %.2f | grid %.3f | band %.3f | polish dispatch %.3f (wall %.3f) | re-solve alone %.3f (wall %.3f, S14_TIME) | cpu %.3f | accepted x %.1f y %.1f | viol %.0f | CPU prep %.2f\n",
         steady_from, steps, TS.p_wall, TS.p_grid, TS.p_band, TS.p_pb, TS.p_pb_wall, TS.p_pa, TS.p_pa_wall, TS.p_cpu,
         TS.p_acc0, TS.p_acc1, TS.p_viol, TS.p_enc);
  const double H257DIV = (g_h257 & 8) ? 1.0 : (double)nsteady;
  if (nsteady && g_h93_print) printf("STEADY H93 MEAN (ms/step, steps %d..%d): scan enc %.3f [alloc %.3f pre %.3f sort %.3f bins %.3f fill %.3f] scan cmd %.3f [pipe %.3f] | GN enc %.3f [alloc %.3f pre %.3f sort %.3f plan %.3f fill %.3f] GN cmd %.3f [pipe %.3f]\n",
         steady_from, steps, TS.scan_enc, TS.h_s_alloc / H257DIV, TS.h_s_pre / H257DIV, TS.h_s_sort / H257DIV, TS.h_s_bins / H257DIV, TS.h_s_fill / H257DIV, TS.h_s_cmd / H257DIV, TS.h_s_pipe / H257DIV,
         TS.gn_enc, TS.h_g_alloc / H257DIV, TS.h_g_pre / H257DIV, TS.h_g_sort / H257DIV, TS.h_g_plan / H257DIV, TS.h_g_fill / H257DIV, TS.h_g_cmd / H257DIV, TS.h_g_pipe / H257DIV);
  if (nsteady && (g_h257 & 8)) printf("STEADY H257 CMD (ms/step, steps %d..%d): GN cmd encode %.3f (outside the GN enc total) | scan cmd encode %.3f | GN enc %.3f | scan enc %.3f | H257=%d\n", steady_from, steps, TS.h_g_cmd, TS.h_s_cmd, TS.gn_enc, TS.scan_enc, g_h257);
  if (nsteady) printf("STEADY PIPE (K=%d, steps %d..%d): step %.1f | head %.1f | GPU busy (union) %.1f = %.0f%% of step [CB spans sum %.1f: scan %.1f + GN %.1f] | main blocked on GPU %.1f (scan %.1f, GN %.1f) | CPU path %.1f | rebuild %.1f\n",
         g_chunks, steady_from, steps, TS.timed(), TS.head, TS.gpu_union, 100.0 * TS.gpu_union / std::max(1e-9, TS.timed()), TS.scan_gpu + TS.gn_gpu, TS.scan_gpu, TS.gn_gpu,
         TS.w_s + TS.w_g, TS.w_s, TS.w_g, TS.pipe_wall - (TS.w_s + TS.w_g), TS.rebuild);
  if (nsteady) printf("STEADY H131 iter wall %.1f | timed %.1f | untimed %.1f | post pre-CB parts: alloc %.2f | bind %.2f | enc %.2f | pre-CB total %.2f | glue install %.2f | h4 wait %.2f | h4 sets %.2f | H131=%d\n",
         TS.h131_iter, TS.timed(), TS.h131_iter - TS.timed(), TS.h131_alloc, TS.h131_bind, TS.h131_enc, TS.h131_pre, TS.h131_inst, TS.h131_h4wait, TS.h131_h4, g_h131);
  { long eh = 0, em = 0, el = 0; for (const Worker& w : WK(E)) { eh += w.CH->enc_hit; em += w.CH->enc_miss; el += w.CH->enc_local; }
    long nh = 0, nm = 0, nd = 0; memo_counts(E, nh, nm, nd);
    printf("MEMO v7: NF memo entries %u (hits %ld, misses %ld, refused %ld) | cc14 entries %u (hits %ld, misses %ld, neg hits %ld) | encode memo hits %ld, misses %ld, local %ld\n",
           (unsigned)E.SH->nent.load(), nh, nm, nd, (unsigned)E.CS->nent.load(),
           [&]{ long x = 0; for (const Worker& w : WK(E)) x += w.CH->C.hits; return x; }(),
           [&]{ long x = 0; for (const Worker& w : WK(E)) x += w.CH->C.misses; return x; }(),
           [&]{ long x = 0; for (const Worker& w : WK(E)) x += w.CH->C.neg_hits; return x; }(), eh, em, el); }
  if (g_h99) {
    if (nsteady) printf("H99 STEADY steps %d..%d: lookups %ld hits %ld | hit rate %.4f | hits / unique trees %.4f\n",
                        steady_from, steps, h99_sl, h99_sh, h99_sl > 0 ? (double)h99_sh / (double)h99_sl : 0.0, h99_su > 0 ? (double)h99_sh / h99_su : 0.0);
    printf("H99 MEMO total: lookups %ld hits %ld inserted %ld evicted %ld refused %ld dups %ld | entries %u of %zu | key arena %zu of %zu words | avg key words %.1f | memo %.1f MB\n",
           G99.lookups, G99.hits, G99.inserts, G99.evicts, G99.refused, G99.dups, (unsigned)G99.nlive, G99.cap, (size_t)G99.aptr, G99.arena_w,
           G99.inserts > 0 ? (double)G99.keyw / (double)G99.inserts : 0.0, G99.bytes() / 1048576.0);
  }
  if (g_h242) {
    h242_shutdown();
    if (g_h770) h770_total_print();
    printf("H242 TOTAL: lookups %ld hits %ld hit rate %.4f | lookup ms %.1f (main-thread wall, pool wait included) | apply ms %.1f | puts %ld inserted %ld dups %ld refused_wide %ld clears %ld | verify mismatches %ld collisions %ld | table %zu slots, %zu live\n",
           G242.lookups, G242.hits, G242.lookups > 0 ? (double)G242.hits / (double)G242.lookups : 0.0, G242.lookup_ms, G242.apply_ms,
           G242.nputs, G242.inserts, G242.dups, G242.refused_wide, G242.clears, G242.vmis, G242.vcoll, G242.cap, G242.nlive);
    printf("H242 PUT: stage %.1f ms (main thread, gn_finish) | worker insert %.1f ms over %ld batches (%s)\n",
           G242.put_ms, G242.worker_ms, G242.njobs, G242.sync ? "sync" : "async");
  }
  if (g_h771) h770_total_print();
  if (g_h892) { if (g_h401_on.load()) { g_h401_on.store(0); const double t = h301_now(); h401_push(0, "H892 STEP 4 end", t, t); } h892_dump(g_h892_t0); }
  printf("best admitted child err over run %.6g (rss_x + rss_y)\n", best_all);
  int fp_pick = -1;
  k_flush();
  if (fitx::on()) { std::vector<uint64_t> pk((size_t)NP); for (int i = 0; i < NP; i++) pk[(size_t)i] = fitx::gkey(XC(S, pidx[i]).g); fitx::retain(pk, true);
    long miss = 0, missfin = 0; double mine = 1e300, mine_all = 1e300; for (int i = 0; i < NP; i++) { const double e = XC(S, pidx[i]).err; mine_all = std::min(mine_all, e); int q[2], nt[2]; if (!fitx::qnth(XC(S, pidx[i]).g, q, nt)) { miss++; if (e < 1e20) missfin++; mine = std::min(mine, e); } }
    printf("FIT_EXPORT: population %d, members without a fit record %ld (finite err < 1e20: %ld, lowest err among them %.6g; lowest err in population %.6g)\n", NP, miss, missfin, mine, mine_all); }
  if (getenv("FINAL_PICK")) {
    std::vector<double> pe(NP); std::vector<int> pc(NP);
    for (int i = 0; i < NP; i++) { const evo::Ind& PI = XC(S, pidx[i]); pe[i] = PI.err; pc[i] = sp_show_cx(PI); }
    const char* why[1] = {""};

    const bool fp_nml = g_sicis_en != nullptr && !(getenv("FINAL_PICK_NML") && atoi(getenv("FINAL_PICK_NML")) == 0);
    int pk;
    double sp_pp = NAN, sp_ps = NAN, sp_pt = NAN;
    std::vector<double> sp_totpre;
    if (fp_nml) {
      std::vector<int> pnx(NP), pny(NP), pk_const(NP); std::vector<uint8_t> kmiss((size_t)NP, 0);
      for (int i = 0; i < NP; i++) {
        const evo::Ind& PI = XC(S, pidx[i]);
        pnx[i] = PI.g.len[0]; pny[i] = PI.g.len[1];
        pk_const[i] = fp::nml::const_leaves(PI.g.op[0], PI.g.len[0]) + fp::nml::const_leaves(PI.g.op[1], PI.g.len[1]);


        static const bool kfit = getenv("FP_KFIT") && atoi(getenv("FP_KFIT")) != 0;
        if (kfit) { int q_[2], n_[2]; if (fitx::qnth(PI.g, q_, n_, g_freet ? 1 : 0)) { pk_const[i] = 0; for (int sd = 0; sd < 2; sd++) pk_const[i] += q_[sd] + n_[sd]; } else { pk_const[i] = 0; kmiss[i] = 1; } }
      }
      const bool t13_all = fp::t13_bit(0);


      const int n_fit = 48;
      std::vector<fp::nml::FullTerms> fulls(NP);
      for (int i = 0; i < NP; i++) fulls[i] = fp::fallback_terms(pk_const[i], n_fit, pnx[i], pny[i]);
      {
        static const double w = getenv("SIMPLICITY") ? atof(getenv("SIMPLICITY")) : 1.0;
        if (w > 0 && w != 1.0) for (int i = 0; i < NP; i++) { fulls[i].paramTerm *= w; fulls[i].structTerm *= w; fulls[i].fisherTerm *= w; fulls[i].rangeTerm *= w; fulls[i].total *= w; }
      }
      {

        static const double w = getenv("SIMPLICITY") ? atof(getenv("SIMPLICITY")) : 1.0; const double wk = (w > 0 ? w : 1.0);
        double bestL = INFINITY; for (int i = 0; i < NP; i++) if (!kmiss[(size_t)i] && std::isfinite(pe[i])) bestL = std::min(bestL, pe[i] + fulls[i].total);
        long ncomp = 0, nexcl = 0;
        for (int i = 0; i < NP; i++) { if (!kmiss[(size_t)i]) continue;
          static const bool kmiss_on = !getenv("FP_KMISS") || atoi(getenv("FP_KMISS")) != 0;
          if (!kmiss_on || g_freet || !(pe[i] + fulls[i].total < bestL)) { fulls[i] = fp::fallback_terms(100000, n_fit, pnx[i], pny[i]); nexcl++; continue; }
          const evo::Cand& G = XC(S, pidx[i]).g; int k = 0, nthsum = 0; bool ok = true;
          for (int sd = 0; sd < 2 && ok; sd++) {
            if (g_sf_axis && sd == (g_sf_axis == 1 ? 0 : 1) && G.len[sd] == 1 && G.op[sd][0] == evo::TVAR) { k += 1; continue; }
            static skn::Ctx* C2 = new skn::Ctx(); static skn::Prog* P2 = new skn::Prog(); skn::NFOut nf;
            if (!skn::nf_tree(*C2, G.op[sd], G.c[sd], G.len[sd], &nf) || !skn::compile(*C2, nf, G.c[sd], *P2)) { ok = false; break; }
            skn::View V(P2->w); ::Tree T; if (!compile_tree(V, T)) { ok = false; break; } k += T.Q + T.nth; nthsum += T.nth; }
          if (ok && nthsum > 1) { fulls[i] = fp::fallback_terms(100000, n_fit, pnx[i], pny[i]); nexcl++; continue; }
          fulls[i] = fp::fallback_terms(ok ? k : 100000, n_fit, pnx[i], pny[i]); ncomp++;
          if (w > 0 && w != 1.0) { fulls[i].paramTerm *= w; fulls[i].structTerm *= w; fulls[i].fisherTerm *= w; fulls[i].rangeTerm *= w; fulls[i].total *= w; }
          if (ok) pk_const[i] = k; }
        (void)wk; printf("FINAL_PICK kfit: members without stored fit: compiled %ld, excluded %ld\n", ncomp, nexcl);
      }
      static const bool trig_rule = !getenv("TRIG_RULE") || atoi(getenv("TRIG_RULE")) != 0;
      static const bool geo_pick = !getenv("GEO_PICK") || atoi(getenv("GEO_PICK")) != 0;
      auto tot0 = [&](int i) { return (double)pe[i] + fulls[i].total; };

      struct D2Memo { int8_t done = 0, hasrec = 0, excl = 0; int ex = 0; };
      std::vector<D2Memo> d2memo((size_t)NP);
      auto d2eval = [&](int i, bool& hasrec, bool& excl, int& ex) {
        D2Memo& mm = d2memo[(size_t)i];
        if (!mm.done) {
          fitx::Rec fr; if (fitx::find(XC(S, pidx[i]).g, fr)) { mm.hasrec = 1;
            static const bool fa = !getenv("FREET_AFFINE") || atoi(getenv("FREET_AFFINE")) != 0;
            static const bool al = !getenv("ALIAS_RULE") || atoi(getenv("ALIAS_RULE")) != 0;
            if (fa && geop::freet_nonaffine(fr)) mm.excl = 1;
            else if (al && geop::aliased(fr)) mm.excl = 1;
            else mm.ex = geop::trig_extra_k(fr); }
          mm.done = 1; }
        hasrec = mm.hasrec != 0; excl = mm.excl != 0; ex = mm.ex; };
      const std::vector<double> pe_orig = pe;
      if (g_sp_diag) { sp_totpre.resize((size_t)NP); for (int i = 0; i < NP; i++) sp_totpre[(size_t)i] = fulls[i].total; }
      const int pk_old = fp::pick_index_nml_full(pe, pc, fulls, nullptr, t13_all);
      const auto t_geo0 = std::chrono::steady_clock::now();


      static const bool pick_rescore = !getenv("PICK_RESCORE") || atoi(getenv("PICK_RESCORE")) != 0;
      bool rescored_all = false;
      int pg_bestidx = -2; double pg_best_ex = NAN;
      static const bool pick_gpu = getenv("PICK_RESCORE_GPU") && atoi(getenv("PICK_RESCORE_GPU")) != 0;
      static const bool pick_gpu_chk = getenv("PICK_RESCORE_GPU_CHECK") && atoi(getenv("PICK_RESCORE_GPU_CHECK")) != 0;
      if (pick_rescore && geo_pick && g_sicis_en) {
        const auto t_r0 = std::chrono::steady_clock::now();
        std::vector<int> rc; rc.reserve((size_t)NP); for (int i = 0; i < NP; i++) if (std::isfinite(pe[i]) && fulls[i].paramTerm < 1e3) rc.push_back(i);
        std::vector<double> gg((size_t)rc.size(), NAN); std::vector<uint8_t> isk1(rc.size(), 0);
        const size_t RB = 32, nrb = (rc.size() + RB - 1) / RB;
        bool gpu_done = false;
        
        static long pg_calls = 0, pg_fallbacks = 0;
        if (pick_gpu && t13_all && !g_nres && !g_yfx) {
          pg_calls++;
          typedef std::chrono::steady_clock::time_point TP;
          auto msd = [](TP a, TP b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
          const size_t nr = rc.size();
          const TP tg0 = std::chrono::steady_clock::now();
          std::vector<uint64_t> keys(nr);
          PL(E).run((int)nrb, [&](int b) { for (size_t r = (size_t)b * RB; r < std::min(nr, (size_t)b * RB + RB); r++) keys[r] = fitx::gkey(XC(S, pidx[rc[r]]).g); });
          fitx::GeoPack GP; fitx::pack_geo(keys, 128, GP);
          const TP tg1 = std::chrono::steady_clock::now();
          std::vector<float> cf, gnll, gcond; std::vector<uint8_t> gbad; double ms_init = 0, ms_gpu = 0;
          static const int pg_chunk = getenv("PICK_RESCORE_GPU_CHUNK") ? std::max(1, atoi(getenv("PICK_RESCORE_GPU_CHUNK"))) : 512;
          static const bool pg_sgn = !getenv("GEO_SIGN") || atoi(getenv("GEO_SIGN")) != 0;
          const bool gok = post_sicis_cfc(*g_sicis_en, cf) && pick_gpu_run(E, GP, g_fit_yv, cf, pg_chunk, (float)geop::g_sc[0], (float)geop::g_sc[1], pg_sgn ? 1 : 0, gnll, gbad, gcond, ms_init, ms_gpu);
          if (gok) {
            const TP tg2 = std::chrono::steady_clock::now();
            std::vector<uint8_t> state(nr, 0); std::vector<double> ga(nr, NAN);
            long n_unk_kernel = 0, n_gbad = 0;
            { size_t c = 0; for (size_t r = 0; r < nr; r++) { if (GP.st[r] == 1) { if (gbad[c]) { state[r] = 2; n_gbad++; } else { state[r] = 1; ga[r] = gnll[c]; } c++; } else if (GP.st[r] == 2) { state[r] = 2; n_unk_kernel++; } } }
            std::vector<double> valv(pe); std::vector<uint8_t> exv((size_t)NP, 1);
            std::vector<int> rof((size_t)NP, -1); for (size_t r = 0; r < nr; r++) rof[(size_t)rc[r]] = (int)r;
            for (size_t r = 0; r < nr; r++) { if (state[r] == 1) { valv[(size_t)rc[r]] = ga[r]; exv[(size_t)rc[r]] = 0; } else if (state[r] == 2) exv[(size_t)rc[r]] = 0; }
            double eobs = 0, eover = 0, eunder = 0; long nex = 0, nex1 = 0;
            std::vector<int> bandv(nr, 0); double berr[20] = {0};
            auto exactify = [&](const std::vector<size_t>& rs) {
              std::vector<double> o(rs.size(), NAN); const size_t B8 = 8, nb8 = (rs.size() + B8 - 1) / B8;
              PL(E).run((int)nb8, [&](int b) { for (size_t k = (size_t)b * B8; k < std::min(rs.size(), (size_t)b * B8 + B8); k++) { const int i = rc[rs[k]];
                fitx::Rec fr; if (!fitx::find(XC(S, pidx[i]).g, fr)) { o[k] = pe[i]; continue; }
                o[k] = geop::geo_nll(fr, g_fit_yv, g_fit_yv + 48, *g_sicis_en, nullptr); } });
              for (size_t k = 0; k < rs.size(); k++) { const size_t r = rs[k]; const int i = rc[r]; double g = o[k]; if (!std::isfinite(g)) g = INFINITY;
                if (state[r] == 1 && std::isfinite(g)) { const double d = ga[r] - g; berr[bandv[r]] = std::max(berr[bandv[r]], std::fabs(d)); eobs = std::max(eobs, std::fabs(d)); eover = std::max(eover, d); eunder = std::max(eunder, -d); nex1++; }
                valv[(size_t)i] = g; exv[(size_t)i] = 1; nex++; } };
            auto pgd2dp = [&](int ex) { static const double swl = getenv("SIMPLICITY") ? std::max(atof(getenv("SIMPLICITY")), 1e-300) : 1.0; return 0.5 * ex * std::log((double)n_fit / (2.0 * fp::nml::kPi)) * swl; };
            static const double pg_M0 = getenv("PICK_RESCORE_GPU_MARGIN") ? atof(getenv("PICK_RESCORE_GPU_MARGIN")) : 0.05;
            static const size_t pg_S = getenv("PICK_RESCORE_GPU_S") ? (size_t)atoi(getenv("PICK_RESCORE_GPU_S")) : 64;
            static const size_t pg_K0 = getenv("PICK_RESCORE_GPU_K") ? (size_t)atoi(getenv("PICK_RESCORE_GPU_K")) : 256;
            const size_t pg_cap = std::max<size_t>(4096, nr / 8);
            auto key_u = [&](int i) { return valv[(size_t)i] + fulls[(size_t)i].total; };
            std::vector<size_t> need;
            { std::vector<size_t> c1; for (size_t r = 0; r < nr; r++) { if (state[r] == 2) need.push_back(r); else if (state[r] == 1 && std::isfinite(key_u(rc[r]))) c1.push_back(r); }
              double gmu = INFINITY; for (size_t r : c1) gmu = std::min(gmu, key_u(rc[r]));
              std::vector<std::vector<size_t>> bl(20);
              for (size_t r : c1) { const double G = key_u(rc[r]) - gmu; int bd = G < 1.0 ? 0 : std::min(19, 1 + (int)std::floor(std::log2(G))); bandv[r] = bd; bl[(size_t)bd].push_back(r); }
              std::sort(c1.begin(), c1.end(), [&](size_t a, size_t b) { return key_u(rc[a]) < key_u(rc[b]); });
              const size_t k0 = std::min(pg_K0, c1.size());
              need.insert(need.end(), c1.begin(), c1.begin() + (long)k0);

              for (auto& v : bl) { std::sort(v.begin(), v.end(), [&](size_t a, size_t b) { return key_u(rc[a]) < key_u(rc[b]); }); const size_t m_ = v.size(), S_ = std::min<size_t>(pg_S, m_); for (size_t k = 0; k < S_; k++) need.push_back(v[k * m_ / S_]); }
              std::sort(need.begin(), need.end()); need.erase(std::unique(need.begin(), need.end()), need.end()); }
            const TP tg3 = std::chrono::steady_clock::now();
            int passes = 0; bool valid = false; int bestidx = -1; double best_ex = INFINITY, Bm = pg_M0; long nblock_total = 0;
            for (; passes < 64; ) {
              if (need.size() + (size_t)nex > pg_cap) break;
              exactify(need); need.clear(); passes++;
              Bm = std::max(pg_M0, eobs);
              std::vector<int> El, Al; for (int i = 0; i < NP; i++) { if (!std::isfinite(valv[(size_t)i]) || !std::isfinite(key_u(i))) continue; (exv[(size_t)i] ? El : Al).push_back(i); }
              std::sort(El.begin(), El.end(), [&](int a, int b) { return key_u(a) < key_u(b); });
              std::sort(Al.begin(), Al.end(), [&](int a, int b) { return key_u(a) < key_u(b); });
              best_ex = INFINITY; bestidx = -1;
              for (int i : El) { const double u = key_u(i); if (u >= best_ex) break; bool hr, ex_; int ex; d2eval(i, hr, ex_, ex); if (ex_) continue; const double ct = u + (hr ? pgd2dp(ex) : 0.0); if (ct < best_ex) { best_ex = ct; bestidx = i; } }
              std::vector<int> blockers;
              for (int i : Al) { const double lb = key_u(i) - std::max(pg_M0, berr[bandv[(size_t)rof[(size_t)i]]]); if (lb >= best_ex) break; bool hr, ex_; int ex; d2eval(i, hr, ex_, ex); if (ex_) continue; if (lb + (hr ? pgd2dp(ex) : 0.0) < best_ex) blockers.push_back(i); }
              if (blockers.empty()) { valid = bestidx >= 0; break; }
              nblock_total += (long)blockers.size();
              for (size_t k = 0; k < blockers.size() && k < 2048; k++) need.push_back((size_t)rof[(size_t)blockers[k]]);
            }
            const TP tg4 = std::chrono::steady_clock::now();
            if (valid) {
              gpu_done = true;
              pg_bestidx = bestidx; pg_best_ex = best_ex;
              for (size_t r = 0; r < nr; r++) { const int i = rc[r]; if (state[r] == 0) { gg[r] = pe[i]; } else { gg[r] = valv[(size_t)i]; isk1[r] = 1; } }
              std::string bstr; for (int b = 0; b < 20; b++) if (berr[b] > 0) { char tb[48]; snprintf(tb, sizeof tb, "%d:%.3g ", b, berr[b]); bstr += tb; }
              printf("PICK_RESCORE_GPU: %zu members, packed %ld (not scorable by the kernel %ld, gpu-flagged %ld) | exact on CPU %ld (%.2f%%) in %d pass(es), %ld blockers | gpu-vs-cpu geo on the exact set (%ld): max|e| %.4g, gpu-over %.4g, gpu-under %.4g | per-band max|e| (band b = GPU total within 2^(b-1) nats of the best, b=0 within 1) [%s] | validated pick %d total %.6g | ms: pack %.1f (keys+pack) init %.1f gpu %.1f exact %.1f | fallbacks %ld of %ld\\n",
                     nr, GP.npacked, n_unk_kernel, n_gbad, nex, 100.0 * (double)nex / (double)nr, passes, nblock_total, nex1, eobs, eover, eunder, bstr.c_str(), bestidx, best_ex, msd(tg0, tg1), ms_init, ms_gpu, msd(tg3, tg4), pg_fallbacks, pg_calls);
              if (pick_gpu_chk) {
                const std::vector<double> valv_mixed = valv; const std::vector<uint8_t> exv_mixed = exv;
                std::vector<size_t> rest; for (size_t r = 0; r < nr; r++) if (state[r] == 1 && !exv[(size_t)rc[r]]) rest.push_back(r);
                const double eobs_top = eobs, eover_top = eover; exactify(rest);
                double mx = 0, mx_near = 0, mx_top = 0; long nbig = 0, nnear = 0; double bestg = INFINITY; for (size_t r = 0; r < nr; r++) if (state[r] == 1) bestg = std::min(bestg, valv[(size_t)rc[r]] + fulls[(size_t)rc[r]].total);
                for (size_t r = 0; r < nr; r++) if (state[r] == 1) { const int i = rc[r]; const double e = std::fabs(ga[r] - valv[(size_t)i]); if (!std::isfinite(e)) continue; mx = std::max(mx, e); if (e > pg_M0) nbig++; if (valv[(size_t)i] + fulls[(size_t)i].total < bestg + 30.0) { mx_near = std::max(mx_near, e); nnear++; } }
                (void)mx_top;
                if (const char* dp = getenv("PICK_RESCORE_GPU_DUMP")) { if (FILE* f = fopen(dp, "w")) { size_t c = 0; for (size_t r = 0; r < nr; r++) if (GP.st[r] == 1) { if (state[r] == 1) fprintf(f, "%.9g %.9g %.9g %.9g\n", (double)gcond[c], ga[r], valv[(size_t)rc[r]], fulls[(size_t)rc[r]].total); c++; } fclose(f); } }
                int fbest = -1; double fbest_ex = INFINITY; { std::vector<int> El; for (int i = 0; i < NP; i++) if (std::isfinite(valv[(size_t)i]) && std::isfinite(key_u(i))) El.push_back(i);
                  std::sort(El.begin(), El.end(), [&](int a, int b) { return key_u(a) < key_u(b); });
                  for (int i : El) { const double u = key_u(i); if (u >= fbest_ex) break; bool hr, ex_; int ex; d2eval(i, hr, ex_, ex); if (ex_) continue; const double ct = u + (hr ? pgd2dp(ex) : 0.0); if (ct < fbest_ex) { fbest_ex = ct; fbest = i; } } }
                printf("PICK_RESCORE_GPU_CHECK: all %zu free-t members exact | max|gpu-cpu| all %.4g (members above %.3g nats: %ld) | within 30 nats of the best total (%ld members): %.4g | exact-set error before the check %.4g (over %.4g) | full-CPU pick %d total %.12g vs validated pick %d total %.12g (diff %.3g) -> %s\\n",
                       nr, mx, pg_M0, nbig, nnear, mx_near, eobs_top, eover_top, fbest, fbest_ex, bestidx, best_ex, best_ex - fbest_ex, (fbest == bestidx) ? "SAME PICK" : (best_ex == fbest_ex ? "EXACT TIE, different member" : "PICK DIFFERS"));
                { long rk = -1; if (fbest >= 0 && rof[(size_t)fbest] >= 0 && state[(size_t)rof[(size_t)fbest]] == 1) { const double gu = ga[(size_t)rof[(size_t)fbest]] + fulls[(size_t)fbest].total; rk = 0; for (size_t r = 0; r < nr; r++) if (state[r] == 1 && ga[r] + fulls[(size_t)rc[r]].total < gu) rk++; }
                  printf("PICK_RESCORE_GPU_CHECK_RANK: the full-CPU argmin (member %d) has rank %ld (0 = first) in the GPU order of uncharged totals among %zu GPU-scored members (state %d); |exact set| was %ld\n", fbest, rk, nr, fbest >= 0 && rof[(size_t)fbest] >= 0 ? (int)state[(size_t)rof[(size_t)fbest]] : -1, nex - (long)rest.size()); }
                valv = valv_mixed; exv = exv_mixed;
                for (size_t r = 0; r < nr; r++) { const int i = rc[r]; if (state[r] != 0) gg[r] = valv[(size_t)i]; } }
            } else { pg_fallbacks++; printf("PICK_RESCORE_GPU: validation did not close within the exact-set cap (%zu members exact, %d passes): full CPU rescore | fallbacks %ld of %ld\\n", (size_t)nex, passes, pg_fallbacks, pg_calls); }
          } else { pg_fallbacks++; printf("PICK_RESCORE_GPU: GPU path unavailable (cfc/pipeline/command buffer): full CPU rescore | fallbacks %ld of %ld\\n", pg_fallbacks, pg_calls); }
        }

        if (!gpu_done) PL(E).run((int)nrb, [&](int b) { for (size_t r = (size_t)b * RB; r < std::min(rc.size(), (size_t)b * RB + RB); r++) { const int i = rc[r]; { int q_[2], n_[2]; if (!fitx::qnth(XC(S, pidx[i]).g, q_, n_, 1)) { gg[r] = pe[i]; continue; } }
            fitx::Rec fr; if (!fitx::find(XC(S, pidx[i]).g, fr)) { gg[r] = pe[i]; continue; }
            const float* Y0 = fr.kind == 1 ? g_fit_yv : g_yfx ? g_yfx_f[0] : fitx::g_target[0]; const float* Y1 = fr.kind == 1 ? g_fit_yv + 48 : g_yfx ? g_yfx_f[1] : fitx::g_target[1];
            static const bool geo_fixed_r = getenv("GEO_FIXED") && atoi(getenv("GEO_FIXED")) != 0;
            if (fr.kind != 1 && !geo_fixed_r) { gg[r] = pe[i]; continue; }
            isk1[r] = 1; gg[r] = geop::geo_nll(fr, Y0, Y1, *g_sicis_en, g_nres ? (g_yfx ? g_yfx_nt : g_nres_n) : nullptr); } });
        long nk1 = 0, ndrop = 0, nbetter = 0; std::vector<double> st_tot(NP, INFINITY); for (int i = 0; i < NP; i++) if (std::isfinite(pe[i])) st_tot[i] = tot0(i);
        for (size_t r = 0; r < rc.size(); r++) { const int i = rc[r]; nk1 += isk1[r]; if (!isk1[r]) continue; if (!std::isfinite(gg[r])) { pe[i] = INFINITY; ndrop++; continue; } if (gg[r] < pe[i]) nbetter++; pe[i] = gg[r]; }
        rescored_all = true;
        const int pk_r = fp::pick_index_nml_full(pe, pc, fulls, nullptr, t13_all);
        long rank_st = 0; if (pk_r >= 0) for (int i = 0; i < NP; i++) if (st_tot[i] < st_tot[pk_r]) rank_st++;
        const double ms_r = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_r0).count();
        printf("PICK_RESCORE: %zu pickable members, %ld re-scored with geo_nll (free-t records), %ld non-finite dropped, %ld lower than stored | %.1f ms | stored pick %d -> geo pick (pre-D2) %d (stored-total rank of the geo pick %ld)\n", rc.size(), nk1, ndrop, nbetter, ms_r, pk_old, pk_r, rank_st);
      }
      if (trig_rule) {

        std::vector<int> oa; for (int i = 0; i < NP; i++) if (std::isfinite(pe[i])) oa.push_back(i);
        std::sort(oa.begin(), oa.end(), [&](int a, int b) { return tot0(a) < tot0(b); });
        std::priority_queue<double> hp; long nch = 0, nrec = 0, nev = 0, nna = 0;
        static const double sw = getenv("SIMPLICITY") ? std::max(atof(getenv("SIMPLICITY")), 1e-300) : 1.0;
        for (int i : oa) { if (t13_all && hp.size() >= 64 && tot0(i) >= hp.top()) break; nev++;
          bool hasrec_ = false, excl_ = false; int ex = 0; d2eval(i, hasrec_, excl_, ex);
          if (hasrec_) { nrec++;
            if (excl_) { pe[i] = INFINITY; nna++; continue; }
            if (ex > 0) { nch++; const double dp = 0.5 * ex * std::log((double)n_fit / (2.0 * fp::nml::kPi)) * sw; fulls[i].paramTerm += dp; fulls[i].total += dp; } }
          hp.push(tot0(i)); if (hp.size() > 64) hp.pop(); }
        printf("TRIG_RULE: evaluated %ld of %d members (%ld with a fit record, %ld charged, %ld free-t non-affine excluded)\n", nev, NP, nrec, nch, nna);
      }
      pk = fp::pick_index_nml_full(pe, pc, fulls, nullptr, t13_all);
      if (pg_bestidx != -2) { static long pg_bad = 0; const double tk = pk >= 0 ? (double)pe[pk] + fulls[pk].total : NAN; const bool same = pk == pg_bestidx || tk == pg_best_ex;
        if (!same) pg_bad++; printf("PICK_RESCORE_GPU_FINAL: pick after the D2 walk %d total %.12g | validated %d total %.12g | %s%s\n", pk, tk, pg_bestidx, pg_best_ex, same ? "AGREE" : "DISAGREE", (same && pk != pg_bestidx) ? " (equal-total tie)" : ""); }
      if (geo_pick && g_sicis_en && !rescored_all) {
        std::vector<int> ord; for (int i = 0; i < NP; i++) if (std::isfinite(pe[i])) ord.push_back(i);
        const int T = std::min<int>(64, (int)ord.size());
        std::partial_sort(ord.begin(), ord.begin() + T, ord.end(), [&](int a, int b) { return tot0(a) < tot0(b); });
        std::vector<int> sub, exv; std::vector<double> se, onll; std::vector<int> sc; std::vector<fp::nml::FullTerms> sf;
        for (int r = 0; r < T; r++) { const int i = ord[r]; fitx::Rec fr; if (!fit_get(E, XC(S, pidx[i]).g, fr)) continue;
          const float* Y0 = fr.kind == 1 ? g_fit_yv : g_yfx ? g_yfx_f[0] : fitx::g_target[0]; const float* Y1 = fr.kind == 1 ? g_fit_yv + 48 : g_yfx ? g_yfx_f[1] : fitx::g_target[1];
          static const bool geo_fixed = getenv("GEO_FIXED") && atoi(getenv("GEO_FIXED")) != 0;
          const double g = (fr.kind != 1 && !geo_fixed) ? (double)pe[i] : geop::geo_nll(fr, Y0, Y1, *g_sicis_en, g_nres ? (g_yfx ? g_yfx_nt : g_nres_n) : nullptr); if (!std::isfinite(g)) continue;
          sub.push_back(i); exv.push_back(geop::trig_extra_k(fr)); se.push_back(g); onll.push_back(pe[i]); sc.push_back(pc[i]); sf.push_back(fulls[i]); }
        if (!sub.empty()) {
          const int q = fp::pick_index_nml_full(se, sc, sf, nullptr, t13_all);
          if (q >= 0) { int qo = -1; for (size_t r = 0; r < sub.size(); r++) if (sub[r] == pk) qo = (int)r;
            printf("GEO_PICK: %zu members re-scored | old pick %d nll %.6g -> geo nll %.6g | new pick %d nll %.6g -> geo nll %.6g | %s\n", sub.size(), pk, (double)pe[pk], qo >= 0 ? se[qo] : NAN, sub[q], onll[q], se[q], sub[q] == pk ? "unchanged" : "CHANGED");
            { double bp = 1e300, bu = 1e300; int iu = -1; for (size_t r = 0; r < sub.size(); r++) { const double L = se[r] + sf[r].total; bp = std::min(bp, L); if (exv[r] == 0 && L < bu) { bu = L; iu = (int)r; } }
              printf("GEO_PICK_GAP: pick charged k+%d | best plain member %d total %.4g vs pick total %.4g | gap %.4g nats (= %.2f params at 0.5 log(n/2pi) each)\n", exv[q], iu >= 0 ? sub[iu] : -1, bu, bp, bu - bp, (bu - bp) / (0.5 * std::log(48.0 / (2.0 * fp::nml::kPi)))); }
            pe[sub[q]] = se[q]; pk = sub[q]; }
        } else printf("GEO_PICK: no member re-scored\n");
      }
      if (pk >= 0) { sp_pp = fulls[pk].paramTerm; sp_ps = fulls[pk].structTerm; sp_pt = fulls[pk].total; }
      auto fexpr = [&](int i) { std::string o; fitx::Rec fr; if (i < 0 || !fit_get(E, XC(S, pidx[i]).g, fr)) return std::string("?"); for (int sd = 0; sd < 2; sd++) { uib::Tree T; std::memset(&T, 0, sizeof T);
          o += sd ? " | y=" : "x="; if (fitx::emit_tree(fr, sd, T)) { int pos = T.len - 1; o += fp::expr_of(T.op, T.c, T.len, pos); } else o += "?"; } return o; };
      if (getenv("NOISE_TIE") && atoi(getenv("NOISE_TIE")) != 0) {
        std::unordered_map<uint64_t, int> rep;
        for (int i = 0; i < NP; i++) { if (!std::isfinite(pe[i]) || pe[i] > 1e20) continue;
          const evo::Cand& G = XC(S, pidx[i]).g; uint64_t h = 1469598103934665603ull; for (int sd = 0; sd < 2; sd++) { h = (h ^ (uint64_t)(G.len[sd] + 1000 * sd)) * 1099511628211ull; for (int k = 0; k < G.len[sd]; k++) h = (h ^ (uint64_t)G.op[sd][k]) * 1099511628211ull; }
          auto it = rep.find(h); if (it == rep.end()) rep[h] = i; else if (pe[i] < pe[it->second]) it->second = i; }
        std::vector<int> o; for (auto& kv : rep) o.push_back(kv.second); std::sort(o.begin(), o.end(), [&](int a, int b) { const int na = pnx[a] + pny[a], nb = pnx[b] + pny[b]; return na < nb || (na == nb && (pe[a] < pe[b] || (pe[a] == pe[b] && a < b))); });
        printf("NOISE_TIE: simplest structures (nx+ny, nll):");
        for (size_t r = 0; r < o.size() && r < 12; r++) { const int i = o[r]; int q_[2], n_[2]; const bool fr0 = fitx::qnth(XC(S, pidx[i]).g, q_, n_, g_freet ? 1 : 0); printf(" [m%d nx %d ny %d nll %.6g k %d %s]", i, pnx[i], pny[i], (double)pe[i], fr0 ? q_[0] + q_[1] + n_[0] + n_[1] : -1, fr0 ? "fit" : "NO FIT"); }
        printf("\n"); }
      if (getenv("FP_EXPLAIN") && atoi(getenv("FP_EXPLAIN")) != 0) {
        const int NEX = getenv("FP_EXPLAIN_N") ? atoi(getenv("FP_EXPLAIN_N")) : 3000;
        const double half = 0.5 * std::log(48.0 / (2.0 * fp::nml::kPi));
        auto kfs = [&](int i) { return (double)pk_const[i] * half + fp::nml::struct_code_canon(pnx[i]) + fp::nml::struct_code_canon(pny[i]); };
        auto shash = [&](int i) { const evo::Cand& G = XC(S, pidx[i]).g; uint64_t h = 1469598103934665603ull; for (int sd = 0; sd < 2; sd++) { h = (h ^ (uint64_t)(G.len[sd] + 1000 * sd)) * 1099511628211ull; for (int k = 0; k < G.len[sd]; k++) h = (h ^ (uint64_t)G.op[sd][k]) * 1099511628211ull; } return h; };
        std::unordered_map<uint64_t, int> rt, rn; long nrec = 0, nnor = 0; double bnor = 1e300;
        for (int i = 0; i < NP; i++) { if (!std::isfinite(pe_orig[i]) || pe_orig[i] > 1e20) continue; int q_[2], n_[2];
          if (!fitx::qnth(XC(S, pidx[i]).g, q_, n_, g_freet ? 1 : 0)) { nnor++; bnor = std::min(bnor, (double)pe_orig[i]); continue; } nrec++;
          const uint64_t h = shash(i); auto it = rt.find(h); if (it == rt.end()) rt[h] = i; else if (pe_orig[i] + kfs(i) < pe_orig[it->second] + kfs(it->second)) it->second = i;
          auto jt = rn.find(h); if (jt == rn.end()) rn[h] = i; else if (pe_orig[i] < pe_orig[jt->second]) jt->second = i; }
        std::vector<int> a1, a2; for (auto& kv : rt) a1.push_back(kv.second); for (auto& kv : rn) a2.push_back(kv.second);
        std::sort(a1.begin(), a1.end(), [&](int a, int b) { const double ta = pe_orig[a] + kfs(a), tb = pe_orig[b] + kfs(b); return ta < tb || (ta == tb && a < b); });
        std::sort(a2.begin(), a2.end(), [&](int a, int b) { return pe_orig[a] < pe_orig[b] || (pe_orig[a] == pe_orig[b] && a < b); });
        std::vector<int> sel(a1.begin(), a1.begin() + std::min<size_t>(a1.size(), (size_t)NEX)); std::unordered_set<int> in(sel.begin(), sel.end());
        for (size_t r = 0; r < a2.size() && r < 300; r++) if (in.insert(a2[r]).second) sel.push_back(a2[r]);
        if (pk >= 0 && in.insert(pk).second) sel.push_back(pk);
        auto compk = [&](int i, int* q_, int* n_) -> bool { const evo::Cand& G = XC(S, pidx[i]).g; static skn::Ctx* C3 = new skn::Ctx(); static skn::Prog* P3 = new skn::Prog();
          for (int sd = 0; sd < 2; sd++) { if (g_sf_axis && sd == (g_sf_axis == 1 ? 0 : 1) && G.len[sd] == 1 && G.op[sd][0] == evo::TVAR) { q_[sd] = 1; n_[sd] = 0; continue; }
            skn::NFOut nf; if (!skn::nf_tree(*C3, G.op[sd], G.c[sd], G.len[sd], &nf) || !skn::compile(*C3, nf, G.c[sd], *P3)) return false; skn::View V(P3->w); ::Tree T; if (!compile_tree(V, T)) return false; q_[sd] = T.Q; n_[sd] = T.nth; } return true; };
        std::unordered_map<uint64_t, int> rnr; for (int i = 0; i < NP; i++) { if (!std::isfinite(pe_orig[i]) || pe_orig[i] > 1e20) continue; int q_[2], n_[2]; if (fitx::qnth(XC(S, pidx[i]).g, q_, n_, g_freet ? 1 : 0)) continue; const uint64_t h = shash(i); auto it = rnr.find(h); if (it == rnr.end()) rnr[h] = i; else if (pe_orig[i] < pe_orig[it->second]) it->second = i; }
        std::vector<int> a3; for (auto& kv : rnr) a3.push_back(kv.second); std::sort(a3.begin(), a3.end(), [&](int a, int b) { return pe_orig[a] < pe_orig[b] || (pe_orig[a] == pe_orig[b] && a < b); });
        const size_t NNR = getenv("FP_EXPLAIN_NOREC_N") ? (size_t)atol(getenv("FP_EXPLAIN_NOREC_N")) : 400;
        for (size_t r = 0; r < a3.size() && r < NNR; r++) if (in.insert(a3[r]).second) sel.push_back(a3[r]);
        printf("FP_EXPLAIN_HEAD: pop %d | members with fit record %ld (distinct structures %zu) | finite without record %ld (lowest nll %.6g; not pickable, k=100000) | printed %zu | free_t %d | kfit %d | engine pick %d | warp: %s | n_fit 48, half_log %.6g\n",
               NP, nrec, a1.size(), nnor, bnor, sel.size(), g_freet ? 1 : 0, getenv("FP_KFIT") ? atoi(getenv("FP_KFIT")) : 0, pk, g_freet ? "10 warp params (u0..u8, t0) NOT charged" : "n/a (fixed-t)", half);
        static const bool fa_ = !getenv("FREET_AFFINE") || atoi(getenv("FREET_AFFINE")) != 0, al_ = !getenv("ALIAS_RULE") || atoi(getenv("ALIAS_RULE")) != 0, geo_fixed_ = getenv("GEO_FIXED") && atoi(getenv("GEO_FIXED")) != 0;
        static const char* on_[] = {"", "c", "t", "+", "-", "*", "/", "sin", "cos", "exp", "log", "sq", "cube", "sqrt", "cbrt"};
        std::vector<std::pair<double, int>> tl;
        for (int i : sel) { const evo::Cand& G = XC(S, pidx[i]).g; int q_[2] = {-1, -1}, n_[2] = {-1, -1}; fitx::Rec fr;
          const bool hasrec = fitx::qnth(G, q_, n_, g_freet ? 1 : 0); int kuse = pk_const[i];
          const double c0 = fp::nml::struct_code_canon(pnx[i]), c1 = fp::nml::struct_code_canon(pny[i]);
          if (!hasrec) { const bool ok = !g_freet && compk(i, q_, n_); const int nths = ok ? n_[0] + n_[1] : -1; kuse = ok ? q_[0] + q_[1] + n_[0] + n_[1] : -1;
            if (!ok || nths > 1) { printf("FP_EXPLAIN_NOREC member %d pick %d | nll_stored %.6g | k_comp %d nth_sum %d (-1: free-t or not compilable) param_if_k %.6g | nx %d ny %d struct %.6g | elig 0 | X %s | Y %s\n", i, i == pk ? 1 : 0, (double)pe_orig[i], kuse, nths, kuse >= 0 ? kuse * half : -1.0, pnx[i], pny[i], c0 + c1, fp::tree_str(G, 0).c_str(), fp::tree_str(G, 1).c_str()); continue; } }
          if (!fit_get(E, G, fr)) continue;
          const int ex = geop::trig_extra_k(fr); const bool nonaff = geop::freet_nonaffine(fr), alias = geop::aliased(fr);
          const float* Y0 = fr.kind == 1 ? g_fit_yv : g_yfx ? g_yfx_f[0] : fitx::g_target[0]; const float* Y1 = fr.kind == 1 ? g_fit_yv + 48 : g_yfx ? g_yfx_f[1] : fitx::g_target[1];
          const double geo = (fr.kind != 1 && !geo_fixed_) ? (double)pe_orig[i] : geop::geo_nll(fr, Y0, Y1, *g_sicis_en, g_nres ? (g_yfx ? g_yfx_nt : g_nres_n) : nullptr);
          const double par = (double)kuse * half, d2 = 0.5 * ex * std::log(48.0 / (2.0 * fp::nml::kPi));
          int fl0 = 0, fl1 = 0; for (unsigned v = (unsigned)pnx[i]; v > 1U; v >>= 1U) ++fl0; for (unsigned v = (unsigned)pny[i]; v > 1U; v >>= 1U) ++fl1;
          const double g0 = (2.0 * fl0 + 1.0) * fp::nml::kLn2, g1 = (2.0 * fl1 + 1.0) * fp::nml::kLn2;
          std::string oh[2]; for (int sd = 0; sd < 2; sd++) { int cnt[16] = {0}; for (int k = 0; k < G.len[sd]; k++) if (G.op[sd][k] < 16) cnt[G.op[sd][k]]++; for (int o = 1; o < 15; o++) if (cnt[o]) { char b[24]; snprintf(b, sizeof b, "%s%s%d", oh[sd].empty() ? "" : ",", on_[o], cnt[o]); oh[sd] += b; } }
          std::string ex_[2] = {"?", "?"}; { uib::Tree DT[2]; uib::ds::Stats dst; if (fitx::display_pair(fr, DT, &dst, nullptr)) for (int sd = 0; sd < 2; sd++) { int pos = DT[sd].len - 1; ex_[sd] = fp::expr_of(DT[sd].op, DT[sd].c, DT[sd].len, pos); } }
          printf("FP_EXPLAIN member %d pick %d kind %d | nll_stored %.6g nll_geo %.6g | k %d (q %d+%d nth %d+%d leaves %d) param %.6g | nx %d ny %d struct %.6g (canon %.6g + gamma %.6g | canon %.6g + gamma %.6g) ops_x {%s} ops_y {%s} | trig_extra %d d2 %.6g | nonaffine %d alias %d excl_by_rule %d | total_noD2 %.6g total %.6g | elig %d | X %s | Y %s\n",
                 i, i == pk ? 1 : 0, fr.kind, (double)pe_orig[i], geo, kuse, q_[0], q_[1], n_[0], n_[1], fp::nml::const_leaves(G.op[0], G.len[0]) + fp::nml::const_leaves(G.op[1], G.len[1]), par, pnx[i], pny[i], c0 + c1, c0 - g0, g0, c1 - g1, g1, oh[0].c_str(), oh[1].c_str(), ex, d2,
                 nonaff ? 1 : 0, alias ? 1 : 0, ((fa_ && nonaff) || (al_ && alias)) ? 1 : 0, geo + par + c0 + c1, geo + par + c0 + c1 + d2, fulls[i].paramTerm < 1e3 ? 1 : 0, ex_[0].c_str(), ex_[1].c_str()); tl.push_back({geo + par + c0 + c1 + d2, i}); }
        {
          const int NJ = getenv("FP_EXPLAIN_J") ? atoi(getenv("FP_EXPLAIN_J")) : 600; std::sort(tl.begin(), tl.end());
          for (int r = 0; r < NJ && r < (int)tl.size(); r++) { const int i = tl[(size_t)r].second; fitx::Rec base; if (!fit_get(E, XC(S, pidx[i]).g, base)) continue;
            const float* Y0 = base.kind == 1 ? g_fit_yv : g_yfx ? g_yfx_f[0] : fitx::g_target[0]; const float* Y1 = base.kind == 1 ? g_fit_yv + 48 : g_yfx ? g_yfx_f[1] : fitx::g_target[1];
            fp::nml::Member mem; std::vector<std::array<int, 3>> mp;
            for (int sd = 0; sd < 2; sd++) { if (base.s[sd].drv) continue; for (int q = 0; q < base.s[sd].Q; q++) { mem.theta.push_back((double)base.s[sd].a[q]); mp.push_back({sd, 0, q}); } for (int j = 0; j < base.s[sd].nth; j++) { mem.theta.push_back((double)base.s[sd].th[j]); mp.push_back({sd, 1, j}); } }
            mem.roles.assign(mem.theta.size(), 0); mem.nx = 0; mem.nodesX = pnx[i]; mem.nodesY = pny[i];
            {
              for (size_t u = 0; u < mp.size(); u++) { if (mp[u][1] != 1) continue; const fitx::Side& SS = base.s[mp[u][0]]; const size_t N = SS.n.size(); int nd = -1; for (auto& tr : SS.thref) if (tr.second == mp[u][2]) { nd = tr.first; break; } if (nd < 0) continue;
                std::vector<int> par(N, -1); std::vector<char> hasT(N, 0);
                for (size_t x = 0; x < N; x++) { const fitx::Node& n = SS.n[x]; if (n.op == 2) hasT[x] = 1; else if (n.op >= 3 && n.op <= 6) { hasT[x] = hasT[n.a] || hasT[n.b]; if (par[n.a] < 0) par[n.a] = (int)x; if (par[n.b] < 0) par[n.b] = (int)x; } else if (n.op >= 7) { hasT[x] = hasT[n.a]; if (par[n.a] < 0) par[n.a] = (int)x; } }
                int cur = nd, A = -1; bool mulT = false, anyMul = false;
                while (par[cur] >= 0) { const int p = par[cur]; const fitx::Node& P = SS.n[p]; if (P.op == 3 || P.op == 4) { cur = p; continue; } if (P.op == 5) { anyMul = true; if (hasT[P.a == cur ? P.b : P.a]) mulT = true; cur = p; continue; } A = p; break; }
                int role = 0;
                if (A >= 0) { const int ao = SS.n[A].op; bool up = false; for (int c = A; c >= 0; c = par[c]) if (c != A && SS.n[c].op >= 7 && SS.n[c].op <= 9) up = true;
                  if (ao == 7 || ao == 8) role = mulT ? 1 : (!anyMul ? 2 : 4); else if (ao == 9) role = mulT ? 3 : 4; else role = up ? 4 : 0; }
                mem.roles[u] = role; } }
            fp::nml::ResidualFn rf = [&](const std::vector<double>& th, std::vector<double>& r48) -> bool { fitx::Rec R = base;
              for (size_t u = 0; u < mp.size(); u++) { fitx::Side& SS = R.s[mp[u][0]]; if (mp[u][1] == 0) SS.a[mp[u][2]] = (float)th[u]; else { SS.th[mp[u][2]] = (float)th[u]; for (auto& tr : SS.thref) if (tr.second == mp[u][2]) SS.n[tr.first].c = (float)th[u]; } }
              float c[48]; const double g = geop::geo_nll(R, Y0, Y1, *g_sicis_en, g_nres ? (g_yfx ? g_yfx_nt : g_nres_n) : nullptr, 1024, c); if (!std::isfinite(g)) return false; r48.assign(c, c + 48); return true; };
            fp::nml::Ranges rg; rg.extent = 1.0; rg.freq = 3.14159265358979323846 * 47.0; rg.phase = 2.0 * 3.14159265358979323846; rg.expRateOk = false; fp::nml::FullTerms ft;
            int nrole[5] = {0, 0, 0, 0, 0}; for (int ro : mem.roles) nrole[ro]++; const bool ok = fp::nml::full_terms(mem, rf, *g_sicis_en, rg, ft);
            printf("FP_EXPLAIN_J member %d | ok %d npar %zu kid %d dead %d collinear %d fisher %.6g range %.6g n_cells %d admissible %d out_of_range_id %d | roles extent %d freq %d phase %d exprate %d uncl %d\n", i, ok ? 1 : 0, mem.theta.size(), ft.kId, ft.nDead, ft.nCollinear, ft.fisherTerm, ft.rangeTerm, ft.n, ft.admissible ? 1 : 0, ft.nOutOfRangeId, nrole[0], nrole[1], nrole[2], nrole[3], nrole[4]); } }
      }
      { static double pe_old_save = 0; (void)pe_old_save; printf("GEO_PICK_SUMMARY: before: [%s] | after: [%s]\n", fexpr(pk_old).c_str(), fexpr(pk).c_str()); }
      if (g_rec_pol) printf("REC_POLISH: %ld (child, side) records stored with the polished theta over the run\n", g_rp_n.load());
      printf("GEO_PICK_MS: %.2f (D1+D2) | pick before D1/D2 %d\n", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_geo0).count(), pk_old);
      {
        std::vector<int> ord(NP);
        for (int i = 0; i < NP; i++) ord[i] = i;
        auto tot = [&](int i) { return (double)pe[i] + fulls[i].total; };
        std::partial_sort(ord.begin(), ord.begin() + std::min(10, NP), ord.end(), [&](int a, int b) { return tot(a) < tot(b); });
        for (int r = 0; r < std::min(10, NP); r++) {
          const int i = ord[r];
          const fp::nml::FullTerms& F = fulls[i];
          printf("FINAL_PICK member %d: nll %.6g param %.6g fisher %.6g range %.6g struct %.6g total %.6g admissible %d ok %d nx %d ny %d k %d\n",
                 i, (double)pe[i], F.paramTerm, F.fisherTerm, F.rangeTerm, F.structTerm, tot(i),
                 F.admissible ? 1 : 0, F.ok ? 1 : 0, pnx[i], pny[i], pk_const[i]);
        }
      }
      why[0] = t13_all ? "NML (T13): min(nll + param + fisher + range + struct) over all members"
                       : "NML: min(nll + param + fisher + range + struct) over skyline";
    } else {
      pk = fp::pick_index(pe, pc, 2 * (int)X.size() - 1, 1.96, why, g_sicis_en != nullptr);
    }
    fp_pick = pk;
    if (g_sp_diag && pk >= 0) {
      std::vector<int> ord; for (int i = 0; i < NP; i++) if (XC(S, pidx[i]).err < 1e19f) ord.push_back(i);
      auto nd = [&](int i) { return g_sp_on ? (int)pm[(size_t)i].nx + (int)pm[(size_t)i].ny : (int)XC(S, pidx[i]).cx; };
      auto pn = [&](int i) { return fp::search_pen(pm[(size_t)i].kf, 48, pm[(size_t)i].nx, pm[(size_t)i].ny); };
      std::sort(ord.begin(), ord.end(), [&](int a, int b) { const float ea = XC(S, pidx[a]).err, eb = XC(S, pidx[b]).err; if (ea != eb) return ea < eb; const int na = nd(a), nb = nd(b); return na != nb ? na < nb : a < b; });
      std::vector<int> rk((size_t)NP, 0); int F[64] = {0}; int R = 0;
      for (size_t a = 0; a < ord.size();) { size_t b = a + 1; const float e0 = XC(S, pidx[ord[a]]).err; const int c0 = nd(ord[a]); while (b < ord.size() && XC(S, pidx[ord[b]]).err == e0 && nd(ord[b]) == c0) b++;
        const int r = F[c0] + 1; for (size_t q = a; q < b; q++) rk[(size_t)ord[q]] = r; R = std::max(R, r); for (int c = c0; c < 64; c++) F[c] = std::max(F[c], r); a = b; }
      int bt = -1; double btv = INFINITY; for (int i : ord) { const double t = (double)XC(S, pidx[i]).err + pn(i); if (t < btv) { btv = t; bt = i; } }
      auto domP = [&](int i) { const double e = XC(S, pidx[i]).err, p = pn(i); for (int j : ord) { const double ej = XC(S, pidx[j]).err, pj = pn(j); if (ej <= e && pj <= p && (ej < e || pj < p)) return 1; } return 0; };
      auto domC = [&](int i) { return rk[(size_t)i] > 1 ? 1 : 0; };
      const long kunk = g_sp_kunk.load();
      if (!sp_totpre.empty()) {
        long ncmp = 0, ndiff = 0, nexcl = 0, ndiff_nodes = 0; for (int i : ord) { if (sp_totpre[(size_t)i] > 1e4) { nexcl++; continue; } ncmp++; const double d = std::fabs(sp_totpre[(size_t)i] - pn(i)); if (d > 1e-9) { ndiff++; const evo::Cand& g = XC(S, pidx[i]).g; if ((int)g.len[0] != (int)pm[(size_t)i].nx || (int)g.len[1] != (int)pm[(size_t)i].ny) ndiff_nodes++; } }
        const evo::Cand& gp = XC(S, pidx[pk]).g;
        printf("SEARCH_PEN_DIAG pick vs search penalty: members compared %ld (pick-excluded %ld) | differ by >1e-9: %ld (of which genotype nodes != the canonical (driver 1) nodes: %ld) | pick m%d: pick P %.6g, search P %.6g, genotype nodes %d+%d, canonical nodes %d+%d, pinned driver now %d\n",
               ncmp, nexcl, ndiff, ndiff_nodes, pk, sp_totpre[(size_t)pk], pn(pk), (int)gp.len[0], (int)gp.len[1], (int)pm[(size_t)pk].nx, (int)pm[(size_t)pk].ny, (int)(g_sf_axis ? sf_is_driver(gp) : 0)); }
      printf("SEARCH_PEN_DIAG final pop: members %zu, fronts %d under (nll, nodes) | pick m%d nll %.6g nodes %d k %d P %.4g: front %d (nll,nodes) dominated %d | dominated in (nll,P) %d | best by nll+P m%d total %.6g nodes %d k %d: front %d (nll,nodes), dominated in (nll,P) %d | pick==best %d | k fallback sides %ld\n",
             ord.size(), R, pk, (double)XC(S, pidx[pk]).err, nd(pk), (int)pm[(size_t)pk].kf, pn(pk), rk[(size_t)pk], domC(pk), domP(pk), bt, btv, nd(bt), (int)pm[(size_t)bt].kf, rk[(size_t)bt], domP(bt), (int)(bt == pk), kunk);
      printf("SEARCH_PEN_DIAG select (all steps, per deme pool): pools %ld | best by nll+P is on front 1: %ld, front 2: %ld, front>=3: %ld | ranked position >= DECK (culled by truncation): %ld | migrant copy: %ld\n",
             g_spd[0].load(), g_spd[1].load(), g_spd[2].load(), g_spd[3].load(), g_spd[4].load(), g_spd[5].load());
    }
    if (pk >= 0) {
      const evo::Ind& PI = XC(S, pidx[pk]);
      static uib::Pair* GP = new uib::Pair();
      auto xf_str = [&](uib::Tree T, int sd) -> std::string { if (!uib::to_xframe(T, sd)) return std::string("(too long for the x frame)"); int pos = T.len - 1; return fp::expr_of(T.op, T.c, T.len, pos); };
      if (uib::g_xf_on) {
        uib::fill_from_cand(PI.g, *GP);
        printf("FINAL_PICK frame: t = %s, the driver coordinate in the normalised stroke frame (bbox centre 0, one span for both axes; raw = span*t + centre); engine s = %.9g*t + %.9g\n", g_sf_axis == 1 ? "x" : "y", uib::g_xf_a, uib::g_xf_b);
        printf("FINAL_PICK x(t) = %s\n", xf_str(GP->x, 0).c_str());
        printf("FINAL_PICK y(t) = %s\n", xf_str(GP->y, 1).c_str());
        printf("FINAL_PICK s-frame x(t) = %s\n", fp::tree_str(PI.g, 0).c_str());
        printf("FINAL_PICK s-frame y(t) = %s\n", fp::tree_str(PI.g, 1).c_str());
      } else {
      printf("FINAL_PICK x(t) = %s\n", fp::tree_str(PI.g, 0).c_str());
      printf("FINAL_PICK y(t) = %s\n", fp::tree_str(PI.g, 1).c_str());
      }
      { fitx::Rec fr; double rfs = -1; if (fit_get(E, PI.g, fr, &rfs)) {
          static const char* nm[2] = {"x", "y"};
          printf("FINAL_PICK fitted (%s): x,y = sum_q a_q root_q(t; theta); t = %s\n", fr.kind == 1 ? "free-t joint fit" : fr.kind == 2 ? "fixed-t stage-8 GN theta, a by least squares" : "CPU refit of the genotype (no stored theta for this member; NOT the engine's stored fit)", fr.kind == 1 ? "t0 + sum_k S[i,k] u_k^2 / sum u_k^2 (warp)" : uib::g_xf_on ? "x_i = a0 + d*i/47 (the driver grid, display frame; the fit ran on s = i/47)" : "i/47");
          if (fr.kind == 3) printf("FINAL_PICK fitted refit: rss_x+rss_y %.6g vs stored err %.6g (err is %s)\n", rfs, (double)pe[pk], g_sicis_en ? "an SI-CIS nll, not comparable" : "rss_x+rss_y");
          uib::Tree DT[2], RT[2]; uib::ds::Stats dst; const bool dok = fitx::display_pair(fr, DT, &dst, nullptr), rok = dok && fitx::display_pair(fr, RT, nullptr, nullptr, false);
          for (int sd = 0; sd < 2; sd++) { uib::Tree T; std::memset(&T, 0, sizeof T); int pos = 0;
            if (dok) { T = DT[sd]; pos = T.len - 1; printf("FINAL_PICK fitted %s(t) = %s\n", nm[sd], fp::expr_of(T.op, T.c, T.len, pos).c_str());
              if (dst.applied && rok) { uib::Tree& U = RT[sd]; int upos = U.len - 1; printf("FINAL_PICK fitted %s(t) [unsimplified] = %s\n", nm[sd], fp::expr_of(U.op, U.c, U.len, upos).c_str()); } }
            else if (fitx::emit_tree(fr, sd, T)) { if (uib::g_xf_on) uib::to_xframe(T, sd); pos = T.len - 1; printf("FINAL_PICK fitted %s(t) = %s\n", nm[sd], fp::expr_of(T.op, T.c, T.len, pos).c_str()); }
            else printf("FINAL_PICK fitted %s(t): expression longer than %d nodes, see coefficients\n", nm[sd], uib::MAXLEN);
            printf("FINAL_PICK fitted %s: Q %d nth %d theta [", nm[sd], fr.s[sd].Q, fr.s[sd].nth); for (int k = 0; k < fr.s[sd].nth; k++) printf("%s%.9g", k ? ", " : "", (double)fr.s[sd].th[k]);
            printf("] a ["); for (int q = 0; q < fr.s[sd].Q; q++) printf("%s%.9g", q ? ", " : "", (double)fr.s[sd].a[q]); printf("]%s\n", fr.s[sd].drv ? " (driver axis: x = t)" : ""); }
          if (dok && uib::ds::on()) { if (dst.applied) printf("FINAL_PICK display simplify: mask %d (bit0 phase snap, bit1 dead terms/driver identity) | nodes x %d->%d y %d->%d | rss_x %.6g->%.6g rss_y %.6g->%.6g | rss total %.6g -> %.6g (ratio %.6f) | max curve change over the ink x %.3g y %.3g of the span %.4g\n", dst.mask, dst.len0[0], dst.len1[0], dst.len0[1], dst.len1[1], dst.rss0[0], dst.rss1[0], dst.rss0[1], dst.rss1[1], dst.rss0[0] + dst.rss0[1], dst.rss1[0] + dst.rss1[1], (dst.rss1[0] + dst.rss1[1]) / std::max(dst.rss0[0] + dst.rss0[1], 1e-300), dst.dmax[0], dst.dmax[1], dst.span);
            else printf("FINAL_PICK display simplify: not applied (guards failed or no change possible)\n"); }
          { double tt[48]; if (fr.kind == 1) fitx::warp_t(fr.u, 48, tt); else for (int i = 0; i < 48; i++) tt[i] = i / 47.0;
            const float* Y0 = fr.kind == 1 ? g_fit_yv : fitx::g_target[0]; const float* Y1 = fr.kind == 1 ? g_fit_yv + 48 : fitx::g_target[1];
            printf("FINAL_PICK fitted check: reconstructed rss_x %.6g rss_y %.6g at the 48 data points (model %s)\n", fit_recon_rss(fr.s[0], tt, Y0), fit_recon_rss(fr.s[1], tt, Y1), fr.kind == 1 ? "from the kernel's th/a/u" : fr.kind == 2 ? "th from F14, a recomputed by least squares" : "CPU refit"); }
          if (uib::g_xf_on && fr.kind != 1) {
            double rr[2] = {0, 0}; for (int sd = 0; sd < 2; sd++) { uib::Tree T = DT[sd]; if (!dok) { rr[sd] = NAN; continue; }
              for (int i = 0; i < 48; i++) { const double x = (1.0 / uib::g_xf_a) * (i / 47.0) - uib::g_xf_b / uib::g_xf_a, r = (sd == uib::g_xf_drv ? x : (double)fitx::g_target[sd][i]) - uib::eval_tree(T, x); rr[sd] += r * r; } }
            printf("FINAL_PICK fitted x-frame check: printed equations vs data at x_i = a0 + d*i/47: rss_x %.6g rss_y %.6g\n", rr[0], rr[1]); }
          if (fr.kind == 1) { printf("FINAL_PICK fitted warp u ["); for (int k = 0; k < 10; k++) printf("%s%.9g", k ? ", " : "", (double)fr.u[k]); printf("] (u[9] = t0)\n"); }
        } else printf("FINAL_PICK fitted: no fit record and no refit for the pick (genotype shown above)\n"); }
      printf("FINAL_PICK total: nll %.6g + pen %.6g (param %.6g + struct %.6g, incl. TRIG charge) = %.6g nats | k %d\n", (double)pe[pk], sp_pt, sp_pp, sp_ps, (double)pe[pk] + sp_pt, (int)(sp_totpre.empty() ? -1 : pm[(size_t)pk].kf));
      printf("FINAL_PICK err(%s) %.6g | cx %d | pop %d | reason: %s\n", g_sicis_en ? "SI-CIS nll" : "rss_x+rss_y", (double)pe[pk], pc[pk], NP, why[0]);
    }
  }
  if (uib::g_ring) {
    const int np = std::min(NP, uib::PAIRS); float sc[uib::PAIRS]; uint16_t cxs[uib::PAIRS]; uint32_t ord[uib::PAIRS]; int nk = 0;
    for (int i = 0; i < NP; i++) {
      const float e = XC(S, pidx[i]).err;
      if (nk == np && e >= XC(S, ord[np - 1]).err) continue;
      int j = nk < np ? nk++ : np - 1;
      while (j > 0 && XC(S, ord[j - 1]).err > e) { ord[j] = ord[j - 1]; j--; }
      ord[j] = pidx[i];
    }
    if (fp_pick >= 0) {
      const uint32_t want = pidx[fp_pick]; int at = -1; for (int i = 0; i < nk; i++) if (ord[i] == want) { at = i; break; }
      if (at < 0) at = (nk < np ? nk++ : np - 1);
      for (int i = at; i > 0; i--) ord[i] = ord[i - 1]; ord[0] = want;
    }
    for (int i = 0; i < np; i++) { const evo::Ind& PI = XC(S, ord[i]); sc[i] = PI.err; cxs[i] = (uint16_t)sp_show_cx(PI); }
    auto fillp = [&](int p, uib::Pair& P) { const evo::Cand& G = XC(S, ord[p]).g; fitx::Rec fr; if (fit_get(E, G, fr)) fitx::fill_pair(fr, P, G); else uib::fill_from_cand(G, P); };
    uib::push_burst(-1, h301_now() * 1e3 - uib::g_t0 * 1e3, (int64_t)NP, np, sc, cxs, fillp, 1, 0.f, 1.f);
  }
  if (pf_thr.t.joinable()) pf_thr.t.join();
  delete E.SH; E.SH = nullptr; delete E.CS; E.CS = nullptr; delete E.MM; E.MM = nullptr;
  return 0;
}}
