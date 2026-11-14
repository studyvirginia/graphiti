#import <Foundation/Foundation.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

extern int engine_main(int argc, char** argv);
extern "C" void* sketch_ui_ring_create(void);
extern "C" void sketch_ui_attach(void* ring);
extern "C" void sketch_ui_set_budget(double seconds);

#include <mutex>
static std::mutex g_engine_run;
static int run_engine(const float* xy, int n, double budget_s, void* ring, const char* steps);

#include <csignal>
#include <ctime>
#include <exception>
#include <algorithm>
#include <cerrno>
#include <execinfo.h>
#include <dirent.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <unistd.h>
static char g_logdir[1024] = "", g_session[64] = "";
static int g_console_fd = -1, g_tee_pipe[2] = {-1, -1}, g_engine_fd = -1;
static long g_engine_bytes = 0;
static const long kRotateBytes = 20L << 20;
static const int kKeepLogs = 5;
static volatile sig_atomic_t g_crash_written = 0;

static void write_all(int fd, const char* p, size_t n) {
  while (n > 0) { ssize_t w = write(fd, p, n); if (w < 0) { if (errno == EINTR) continue; return; } p += w; n -= (size_t)w; }
}
static void engine_log_name(char* out, size_t n, int seq) {
  if (seq == 0) snprintf(out, n, "%s/engine-%s.log", g_logdir, g_session);
  else snprintf(out, n, "%s/engine-%s.%d.log", g_logdir, g_session, seq);
}
static void prune_engine_logs(int keep) {
  DIR* d = opendir(g_logdir); if (!d) return;
  struct Ent { char name[256]; time_t mt; };
  std::vector<Ent> v;
  while (struct dirent* e = readdir(d)) {
    if (strncmp(e->d_name, "engine-", 7) != 0) continue;
    char p[1400]; snprintf(p, sizeof p, "%s/%s", g_logdir, e->d_name);
    struct stat st; if (stat(p, &st) != 0) continue;
    Ent x; snprintf(x.name, sizeof x.name, "%s", e->d_name); x.mt = st.st_mtime; v.push_back(x);
  }
  closedir(d);
  std::sort(v.begin(), v.end(), [](const Ent& a, const Ent& b) { return a.mt > b.mt; });
  for (size_t i = (size_t)keep; i < v.size(); i++) { char p[1400]; snprintf(p, sizeof p, "%s/%s", g_logdir, v[i].name); unlink(p); }
}
static int g_rot_seq = 0;
static void open_engine_log(void) {
  char p[1400]; engine_log_name(p, sizeof p, 0);
  g_engine_fd = open(p, O_WRONLY | O_CREAT | O_APPEND, 0644);
  g_engine_bytes = 0;
}
static void rotate_engine_log(void) {
  if (g_engine_fd >= 0) close(g_engine_fd);
  char cur[1400], rot[1400]; engine_log_name(cur, sizeof cur, 0);
  engine_log_name(rot, sizeof rot, ++g_rot_seq);
  rename(cur, rot);
  open_engine_log();
  prune_engine_logs(kKeepLogs);
}

// ---- per-solve performance log (reads the engine's own stdout lines in the tee thread; no engine change) ----
#include <condition_variable>
#include <chrono>
#include <sys/sysctl.h>
struct PerfStep { int st = 0; double iter_ms = 0, admitted = 0, uniq = 0, gpu_busy = 0, gpu_idle = 0, rss_mb = 0; long children = 0, lane_pairs = 0; int permille = -1; };
struct PerfRun { std::vector<PerfStep> steps; int axis = -1, closed = -1, cx = -1; double pso_ms = -1, final_err = NAN; long P = 0; };
static std::mutex g_perf_mu; static std::condition_variable g_perf_cv; static PerfRun g_perf; static int g_perf_sync = 0;
static PerfStep& perf_step(int st) { for (PerfStep& x : g_perf.steps) if (x.st == st) return x; g_perf.steps.push_back(PerfStep()); g_perf.steps.back().st = st; return g_perf.steps.back(); }
static void perf_line(const char* l) {
  int st = 0, a = 0; long c = 0, lp = 0; double x = 0, y = 0, z = 0;
  std::lock_guard<std::mutex> lk(g_perf_mu);
  if (sscanf(l, "step %d | children %ld admitted %lf shortfall %*lf | rej within %*lf vs pop %*lf | uniq trees %lf", &st, &c, &x, &y) == 4) { PerfStep& s = perf_step(st); s.children = c; s.admitted = x; s.uniq = y; }
  else if (sscanf(l, "ITER WALL step %d: %lf ms", &st, &x) == 2) perf_step(st).iter_ms = x;
  else if (sscanf(l, "COFIT step %d: permille %d cpu pairs %ld", &st, &a, &lp) == 3) { PerfStep& s = perf_step(st); s.permille = a; s.lane_pairs = lp; }
  else if (sscanf(l, "H301 step %d: window %lf ms | GPU busy %lf | idle %lf", &st, &x, &y, &z) == 4) { PerfStep& s = perf_step(st); s.gpu_busy = y; s.gpu_idle = z; }
  else if (sscanf(l, "H129 step %d | step ms %lf | RSS %lf MB", &st, &x, &y) == 3) perf_step(st).rss_mb = y;
  else if (sscanf(l, "SINGLE_FAST: axis %d", &a) == 1) { g_perf.axis = a; const char* q = strstr(l, "closed "); if (q) g_perf.closed = atoi(q + 7); }
  else if (strncmp(l, "H178 startup:", 13) == 0) { const char* q = strstr(l, "PSO creation wall "); if (q) g_perf.pso_ms = atof(q + 18); }
  else if (strncmp(l, "FINAL_PICK err(", 15) == 0) { const char* q = strchr(l, ')'); if (q && sscanf(q + 1, " %lf | cx %d", &x, &a) == 2) { g_perf.final_err = x; g_perf.cx = a; } }   // "err(rss_x+rss_y)" or "err(SI-CIS nll)" (noise model on)
  else if (strncmp(l, "s12 loop", 8) == 0) { const char* q = strstr(l, "P="); if (q) g_perf.P = atol(q + 2); }
  else if (sscanf(l, "[graphiti] PERF_SYNC %d", &a) == 1) { g_perf_sync = a; g_perf_cv.notify_all(); }
}
static void perf_feed(const char* buf, size_t n) {   // tee thread only
  static std::string carry; carry.append(buf, n);
  size_t b = 0;
  for (;;) { const size_t e = carry.find('\n', b); if (e == std::string::npos) break; carry[e] = 0; perf_line(carry.c_str() + b); b = e + 1; }
  carry.erase(0, b); if (carry.size() > 65536) carry.clear();
}
static void perf_reset(void) { std::lock_guard<std::mutex> lk(g_perf_mu); g_perf = PerfRun(); }
static const char* perf_thermal(void) { switch (NSProcessInfo.processInfo.thermalState) { case NSProcessInfoThermalStateNominal: return "nominal"; case NSProcessInfoThermalStateFair: return "fair"; case NSProcessInfoThermalStateSerious: return "serious"; default: return "critical"; } }
static void perf_report(const char* kind, int rc, int npts, double budget_s, double wall_ms, const char* therm0) {
  if (g_tee_pipe[1] < 0) return;   // output is not captured (graphiti_diag_init not called): nothing was parsed
  static int seq = 0; const int my = ++seq;
  printf("[graphiti] PERF_SYNC %d\n", my); fflush(stdout);
  PerfRun R;
  { std::unique_lock<std::mutex> lk(g_perf_mu); g_perf_cv.wait_for(lk, std::chrono::milliseconds(1500), [&] { return g_perf_sync >= my; }); R = g_perf; }
  std::vector<PerfStep> v; for (const PerfStep& s : R.steps) if (s.iter_ms > 0) v.push_back(s);
  std::sort(v.begin(), v.end(), [](const PerfStep& a, const PerfStep& b) { return a.st < b.st; });
  const size_t N = v.size();
  double sum = 0, ch = 0, ad = 0, uq = 0, lane = 0, gb = 0, gi = 0, rss = 0, first4 = 0; std::vector<double> ms;
  for (size_t i = 0; i < N; i++) { const PerfStep& s = v[i]; sum += s.iter_ms; ch += (double)s.children; ad += s.admitted; uq += s.uniq; lane += (double)s.lane_pairs; gb += s.gpu_busy; gi += s.gpu_idle; rss = std::max(rss, s.rss_mb); ms.push_back(s.iter_ms); if (i < 4) first4 += s.iter_ms; }
  std::sort(ms.begin(), ms.end());
  auto pct = [&](double q) { return N ? ms[std::min(N - 1, (size_t)(q * (double)(N - 1) + 0.5))] : 0.0; };
  const double sec = sum / 1000.0;
  char model[64] = "?"; size_t ml = sizeof model; sysctlbyname("hw.machine", model, &ml, nullptr, 0);
  const char* mode = R.axis < 0 ? "?" : R.axis == 0 ? (R.closed == 1 ? "free-t closed" : "free-t open") : "single-value";
  printf("[graphiti] PERF %s | %s | rc %d | pts %d | budget %.1f s | device %s | cores %d | thermal %s -> %s\n", kind, mode, rc, npts, budget_s, model, (int)NSProcessInfo.processInfo.activeProcessorCount, therm0, perf_thermal());
  printf("[graphiti] PERF wall %.0f ms | steps %zu | in steps %.0f ms | outside steps (start-up, final pick) %.0f ms | pipeline creation %.0f ms\n", wall_ms, N, sum, wall_ms - sum, R.pso_ms);
  if (N) {
    printf("[graphiti] PERF ms/step mean %.1f | median %.1f | min %.1f | p90 %.1f | max %.1f | first 4 steps %.0f ms\n", sum / (double)N, pct(0.5), ms.front(), pct(0.9), ms.back(), first4);
    printf("[graphiti] PERF candidates: evaluated %.0f (%.0f /s) | admitted and fitted %.0f (%.0f /s) | unique trees %.0f (%.0f /s) | population %ld\n", ch, ch / sec, ad, ad / sec, uq, uq / sec, R.P);
    printf("[graphiti] PERF CPU lane pairs %.0f (%.1f %% of admitted) | GPU busy %.0f ms (%.0f %% of step window) | peak RSS %.0f MB | final score %.6g (nll with the noise model on, else rss) cx %d\n", lane, ad > 0 ? 100.0 * lane / ad : 0.0, gb, gb + gi > 0 ? 100.0 * gb / (gb + gi) : 0.0, rss, R.final_err, R.cx);
    for (const PerfStep& s : v) printf("[graphiti] PERF step %3d | %7.1f ms | children %ld | admitted %.0f (%.0f /s) | uniq %.0f | lane pairs %ld (permille %d) | GPU busy %.1f idle %.1f | RSS %.0f MB\n",
                                       s.st, s.iter_ms, s.children, s.admitted, s.admitted / (s.iter_ms / 1000.0), s.uniq, s.lane_pairs, s.permille, s.gpu_busy, s.gpu_idle, s.rss_mb);
  }
  fflush(stdout);
  // files: Documents/perf/<stamp>.csv (one row per step) and Documents/perf/summary.csv (one row per solve)
  NSFileManager* fm = NSFileManager.defaultManager;
  NSString* dir = [[NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES) firstObject] stringByAppendingPathComponent:@"perf"];
  [fm createDirectoryAtPath:dir withIntermediateDirectories:YES attributes:nil error:nil];
  NSDateFormatter* f = [[NSDateFormatter alloc] init]; f.dateFormat = @"yyyyMMdd-HHmmss"; f.locale = [NSLocale localeWithLocaleIdentifier:@"en_US_POSIX"];
  NSString* stamp = [f stringFromDate:NSDate.date];
  NSMutableString* body = [NSMutableString stringWithFormat:@"# %s,%s,rc=%d,pts=%d,budget_s=%.1f,device=%s,thermal=%s->%s,wall_ms=%.0f,pipeline_ms=%.0f\nstep,ms,children,admitted,admitted_per_s,uniq_trees,lane_pairs,lane_permille,gpu_busy_ms,gpu_idle_ms,rss_mb\n", kind, mode, rc, npts, budget_s, model, therm0, perf_thermal(), wall_ms, R.pso_ms];
  for (const PerfStep& s : v) [body appendFormat:@"%d,%.1f,%ld,%.0f,%.0f,%.0f,%ld,%d,%.1f,%.1f,%.0f\n", s.st, s.iter_ms, s.children, s.admitted, s.admitted / (s.iter_ms / 1000.0), s.uniq, s.lane_pairs, s.permille, s.gpu_busy, s.gpu_idle, s.rss_mb];
  NSString* path = [dir stringByAppendingPathComponent:[NSString stringWithFormat:@"%@-%s.csv", stamp, kind]];
  [body writeToFile:path atomically:YES encoding:NSUTF8StringEncoding error:nil];
  NSString* sumPath = [dir stringByAppendingPathComponent:@"summary.csv"];
  if (![fm fileExistsAtPath:sumPath]) [@"stamp,kind,mode,rc,pts,budget_s,device,thermal_start,thermal_end,wall_ms,steps,ms_in_steps,ms_per_step_mean,ms_per_step_median,ms_per_step_max,first4_ms,evaluated_per_s,admitted_per_s,uniq_per_s,lane_pct,gpu_busy_pct,peak_rss_mb,pipeline_ms,final_score,cx\n" writeToFile:sumPath atomically:YES encoding:NSUTF8StringEncoding error:nil];
  NSString* row = [NSString stringWithFormat:@"%@,%s,%s,%d,%d,%.1f,%s,%s,%s,%.0f,%zu,%.0f,%.1f,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f,%.1f,%.0f,%.0f,%.0f,%.3g,%d\n", stamp, kind, mode, rc, npts, budget_s, model, therm0, perf_thermal(), wall_ms, N, sum,
                   N ? sum / (double)N : 0.0, pct(0.5), N ? ms.back() : 0.0, first4, N ? ch / sec : 0.0, N ? ad / sec : 0.0, N ? uq / sec : 0.0, ad > 0 ? 100.0 * lane / ad : 0.0, gb + gi > 0 ? 100.0 * gb / (gb + gi) : 0.0, rss, R.pso_ms, R.final_err, R.cx];
  if (NSFileHandle* h = [NSFileHandle fileHandleForWritingAtPath:sumPath]) { [h seekToEndOfFile]; [h writeData:[row dataUsingEncoding:NSUTF8StringEncoding]]; [h closeFile]; }
  printf("[graphiti] PERF log: %s\n", path.UTF8String); fflush(stdout);
}
static void* engine_tee_thread(void*) {
  char buf[8192];
  for (;;) {
    ssize_t n = read(g_tee_pipe[0], buf, sizeof buf);
    if (n < 0 && errno == EINTR) continue;
    if (n <= 0) break;
    if (g_console_fd >= 0) write_all(g_console_fd, buf, (size_t)n);
    perf_feed(buf, (size_t)n);
    if (g_engine_fd >= 0) {
      write_all(g_engine_fd, buf, (size_t)n); g_engine_bytes += n;
      if (g_engine_bytes > kRotateBytes) rotate_engine_log();
    }
  }
  return nullptr;
}
static void write_crash_report(const char* kind, const char* detail) {
  if (g_crash_written) return;
  g_crash_written = 1;
  void* bt[64]; const int n = backtrace(bt, 64);
  char p[1400]; snprintf(p, sizeof p, "%s/crash-%ld-%s.log", g_logdir, (long)time(nullptr), kind);
  int fd = open(p, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd < 0) return;
  char head[1200]; int m = snprintf(head, sizeof head, "GraphArt crash: %s %s\nsession %s\nbacktrace:\n", kind, detail, g_session);
  write_all(fd, head, (size_t)m);
  backtrace_symbols_fd(bt, n, fd);
  close(fd);
}
static void graphiti_signal(int sig) {
  char kind[32]; snprintf(kind, sizeof kind, "signal%d", sig);
  write_crash_report(kind, sig == SIGABRT ? "SIGABRT" : sig == SIGSEGV ? "SIGSEGV" : sig == SIGBUS ? "SIGBUS" : "SIGILL");
  signal(sig, SIG_DFL);
  raise(sig);
}
static void graphiti_terminate(void) { write_crash_report("terminate", "std::terminate"); abort(); }
static void graphiti_objc_exception(NSException* e) {
  char msg[700]; snprintf(msg, sizeof msg, "%s: %s", e.name.UTF8String, e.reason.UTF8String ?: "");
  write_crash_report("nsexception", msg);
}
extern "C" void graphiti_diag_init(const char* logdir, const char* session, int verbose_engine) {
  static bool once = false; if (once) return; once = true;
  snprintf(g_logdir, sizeof g_logdir, "%s", logdir); snprintf(g_session, sizeof g_session, "%s", session);
  g_console_fd = dup(1);
  if (pipe(g_tee_pipe) == 0) {
    prune_engine_logs(kKeepLogs - 1);
    if (verbose_engine) open_engine_log();
    dup2(g_tee_pipe[1], 1); dup2(g_tee_pipe[1], 2);
    pthread_t th; if (pthread_create(&th, nullptr, engine_tee_thread, nullptr) == 0) pthread_detach(th);
  }
  for (int sig : {SIGABRT, SIGSEGV, SIGBUS, SIGILL}) signal(sig, graphiti_signal);
  std::set_terminate(graphiti_terminate);
  NSSetUncaughtExceptionHandler(graphiti_objc_exception);
}

#include <execinfo.h>
static void graphiti_atexit() {
  void* bt[32]; const int n = backtrace(bt, 32); fprintf(stdout, "[graphiti] ATEXIT: process exiting, backtrace:\n"); fflush(stdout);
  backtrace_symbols_fd(bt, n, 1); fflush(stdout);
}
extern "C" int graphiti_prepare(void) {
  static bool once = false; if (!once) { once = true; atexit(graphiti_atexit); setvbuf(stdout, nullptr, _IOLBF, 0); }
  printf("[graphiti] prepare start\n");
  std::vector<float> xy; for (int i = 0; i < 48; i++) { const float a = 6.2831853f * i / 47.0f; xy.push_back(200 + 100 * std::cos(a)); xy.push_back(-(400 + 100 * std::sin(a))); }
  void* ring = sketch_ui_ring_create();
  const int rc = run_engine(xy.data(), 48, 0.0, ring, "0");
  printf("[graphiti] prepare rc=%d\n", rc); fflush(stdout);
  return rc;
}

extern "C" int graphiti_solve(const float* xy, int n, double budget_s, void* ring) { return run_engine(xy, n, budget_s, ring, "1000"); }

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <thread>
namespace stb {
constexpr int PAIRS = 9, PTS = 120, MAXLEN = 256;
struct Tree { uint8_t len; uint8_t op[MAXLEN]; int16_t cid[MAXLEN]; float c[MAXLEN]; };
struct Fit { uint8_t kind, ovf; uint8_t q[2], nth[2]; uint8_t drv[2]; uint8_t pad[3]; float u[10]; float th[2][32]; float a[2][16]; };
struct Pair { Tree x, y; float score; uint16_t cx; uint8_t champion; uint8_t fitted; float xs[PTS], ys[PTS]; Fit fit; };
struct Burst { uint32_t seq; int32_t step; int32_t n_pairs; int32_t final_pick; double wall_ms; int64_t candidates; float t_lo, t_hi; Pair pairs[PAIRS]; };
}
extern "C" int sketch_ui_pop(void* ring, void* out);
extern "C" void* sketch_ui_ring_create(void);
extern "C" void sketch_ui_ring_destroy(void* ring);
extern "C" unsigned long sketch_ui_burst_bytes(void);

struct SelftestStroke { std::string name; std::vector<float> xy; };

static std::vector<SelftestStroke> selftest_load(const std::string& path) {
  std::vector<SelftestStroke> out; std::ifstream f(path); std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { const size_t k = line.find_first_not_of(" \t", 1); out.push_back({k == std::string::npos ? std::string() : line.substr(k), {}}); continue; }
    if (out.empty()) continue;
    std::istringstream is(line); float x, y;
    if (is >> x >> y) { out.back().xy.push_back(x); out.back().xy.push_back(y); }
  }
  return out;
}

extern "C" int graphiti_selftest(void) {
  NSString* path = [NSBundle.mainBundle.resourcePath stringByAppendingPathComponent:@"in/selftest_strokes.txt"];
  auto strokes = selftest_load(path.UTF8String);
  printf("[graphiti] selftest: %zu strokes\n", strokes.size()); fflush(stdout);
  int failures = 0;
  const char* st_only = getenv("SELFTEST_ONLY");
  for (const auto& s : strokes) {
    if (st_only && s.name.find(st_only) == std::string::npos) continue;
    const int n = (int)(s.xy.size() / 2);
    void* ring = sketch_ui_ring_create();
    if (!ring) { failures++; continue; }
    struct Acc { int steps = 0; double ms_sum = 0; int final_seen = 0; float xs[stb::PTS] = {}, ys[stb::PTS] = {}; uint16_t cx = 0; float score = 0; } acc;
    std::atomic<bool> stop{false};
    std::thread drain([&] {
      std::vector<unsigned char> buf(sketch_ui_burst_bytes());
      for (;;) {
        const bool done = stop.load();
        const bool got = sketch_ui_pop(ring, buf.data()) == 1;
        if (got) {
          const stb::Burst* b = reinterpret_cast<const stb::Burst*>(buf.data());
          if (b->final_pick) {
            acc.final_seen = 1; acc.cx = b->pairs[0].cx; acc.score = b->pairs[0].score;
            std::memcpy(acc.xs, b->pairs[0].xs, sizeof(acc.xs)); std::memcpy(acc.ys, b->pairs[0].ys, sizeof(acc.ys));
          } else if (b->step >= 0) { acc.steps++; acc.ms_sum += b->wall_ms; }
        } else if (done) break;
        else std::this_thread::sleep_for(std::chrono::milliseconds(2));
      }
    });
    const char* st_steps = getenv("SELFTEST_STEPS"); const char* st_budget = getenv("SELFTEST_BUDGET");
    const int rc = run_engine(s.xy.data(), n, st_budget ? atof(st_budget) : 8.0, ring, st_steps ? st_steps : "1000");
    stop.store(true); drain.join();
    sketch_ui_ring_destroy(ring);
    const bool ok = rc == 0 && acc.final_seen;
    if (!ok) failures++;
    std::string xs, ys;
    char buf[64];
    for (int i = 0; i < stb::PTS; i += 12) { snprintf(buf, sizeof buf, "%s%.3g", i ? "," : "", acc.xs[i]); xs += buf; }
    for (int i = 0; i < stb::PTS; i += 12) { snprintf(buf, sizeof buf, "%s%.3g", i ? "," : "", acc.ys[i]); ys += buf; }
    const double ms = acc.steps ? acc.ms_sum / acc.steps : 0.0;
    printf("SELFTEST label=%s | mode free-t | FINAL_PICK cx %d score %.6g | x(t)=[%s] | y(t)=[%s] | ms/step %.1f | steps %d | rc %d\n",
           s.name.c_str(), (int)acc.cx, acc.score, xs.c_str(), ys.c_str(), ms, acc.steps, rc);
    fflush(stdout);
  }
  printf("[graphiti] selftest done: %zu strokes, %d failures\n", strokes.size(), failures); fflush(stdout);
  return failures;
}

static NSString* g_strokeLog = nil;
extern "C" const char* graphiti_stroke_log_path(void) { return g_strokeLog ? g_strokeLog.UTF8String : nullptr; }

static NSString* stroke_log_create(NSString* body) {
  NSFileManager* fm = NSFileManager.defaultManager;
  NSString* docs = [NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES) firstObject];
  NSString* dir = [docs stringByAppendingPathComponent:@"strokes"];
  [fm createDirectoryAtPath:dir withIntermediateDirectories:YES attributes:nil error:nil];
  NSDateFormatter* f = [[NSDateFormatter alloc] init];
  f.dateFormat = @"yyyyMMdd-HHmmss"; f.locale = [NSLocale localeWithLocaleIdentifier:@"en_US_POSIX"];
  NSString* stamp = [f stringFromDate:NSDate.date];
  NSString* path = [dir stringByAppendingPathComponent:[stamp stringByAppendingString:@".txt"]];
  for (int k = 2; [fm fileExistsAtPath:path]; k++) path = [dir stringByAppendingPathComponent:[NSString stringWithFormat:@"%@-%d.txt", stamp, k]];
  if (![body writeToFile:path atomically:YES encoding:NSUTF8StringEncoding error:nil]) return nil;
  return path;
}

static int run_engine(const float* xy, int n, double budget_s, void* ring, const char* steps) {
  std::lock_guard<std::mutex> lk(g_engine_run);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  g_strokeLog = nil;
  if (!xy || n < 4 || !ring) return 2;
  perf_reset(); const auto perf_t0 = std::chrono::steady_clock::now(); const char* perf_th0 = perf_thermal();
  NSBundle* b = NSBundle.mainBundle;
  NSString* res = b.resourcePath;
  NSString* work = NSTemporaryDirectory();
  NSString* strokePath = [work stringByAppendingPathComponent:@"graphiti_stroke.txt"];
  NSMutableString* s = [NSMutableString string];
  for (int i = 0; i < n; i++) [s appendFormat:@"%.6f %.6f\n", xy[2 * i], xy[2 * i + 1]];
  if (![s writeToFile:strokePath atomically:YES encoding:NSUTF8StringEncoding error:nil]) return 3;
  if (strcmp(steps, "0") != 0) {
    NSString* head = [NSString stringWithFormat:@"# GraphArt stroke: %d pts, budget %.2f s\n", n, budget_s];
    g_strokeLog = stroke_log_create([head stringByAppendingString:s]);
    printf("[graphiti] stroke log: %s\n", g_strokeLog ? g_strokeLog.UTF8String : "(write failed)"); fflush(stdout);
  }

  std::ifstream ef([res stringByAppendingPathComponent:@"in/env.txt"].UTF8String);
  std::string all, tok; std::getline(ef, all);
  std::istringstream is(all);
  while (is >> tok) { size_t e = tok.find('='); if (e != std::string::npos) setenv(tok.substr(0, e).c_str(), tok.substr(e + 1).c_str(), 1); }
  setenv("SICIS", "1", 1); setenv("FINAL_PICK", "1", 1); setenv("SINGLE_FAST", "1", 1);
  setenv("S12_FREET", "1", 1);
  if (!getenv("H242D_FIX")) setenv("H242D_FIX", "0.35", 1);
  NSString* arch = [[NSSearchPathForDirectoriesInDomains(NSCachesDirectory, NSUserDomainMask, YES) firstObject] stringByAppendingPathComponent:[NSString stringWithFormat:@"psos_%llu.bin", (unsigned long long)[[[NSFileManager defaultManager] attributesOfItemAtPath:[res stringByAppendingPathComponent:@"kd/fjn_h950.metallib"] error:nil] fileSize] ^ (unsigned long long)[[[NSFileManager defaultManager] attributesOfItemAtPath:[res stringByAppendingPathComponent:@"kd/scan.metallib"] error:nil] fileSize] << 20]];
  NSString* archShip = [res stringByAppendingPathComponent:@"in/psos_ios.bin"];
  NSFileManager* fm = NSFileManager.defaultManager;
  const bool haveRun = [fm fileExistsAtPath:arch];
  setenv("H178F", (haveRun || ![fm fileExistsAtPath:archShip] ? arch : archShip).UTF8String, 1);
  setenv("H178W", arch.UTF8String, 1);
  setenv("STROKE_IN", strokePath.UTF8String, 1);
  setenv("H930_MEMFIX", "1", 1);
  setenv("FJN_CISFUSE_HOST", "1", 1);
  {
    const NSInteger nm = [[NSUserDefaults standardUserDefaults] integerForKey:@"noiseModelMode"];
    if (getenv("GRAPHITI_NOISE_OFF")) { unsetenv("SICIS"); unsetenv("FJN_CISFUSE_HOST"); }
    if (nm == 1) setenv("NOISE_SCALE", "0.1", 1); else unsetenv("NOISE_SCALE");
  }
  if (!getenv("FJN_CISFUSE_CHECK")) setenv("FJN_CISFUSE_CHECK", "0", 1);
  if (!getenv("FJN_LPT")) setenv("FJN_LPT", "1", 1);
  setenv("SKETCH_SIGMA_V2_BIN", [res stringByAppendingPathComponent:@"in/sigma_v2_zi5.bin"].UTF8String, 1);
  sketch_ui_attach(ring);
  sketch_ui_set_budget(budget_s);
  NSString* kd = [res stringByAppendingPathComponent:@"kd"];
  NSString* out = [work stringByAppendingPathComponent:@"graphiti_out"];
  static const bool seed_file = !getenv("SEED_FILE") || atoi(getenv("SEED_FILE")) != 0;
  const std::string seeds = seed_file ? std::string([res stringByAppendingPathComponent:@"in/cands_r68late.bin"].UTF8String) : std::string("random");
  std::vector<std::string> a = { "engine",
    seeds, kd.UTF8String, out.UTF8String,
    [res stringByAppendingPathComponent:@"in/primitives6_flat.txt"].UTF8String, "circle", steps, "8", "0" };
  std::vector<char*> av; for (auto& x : a) av.push_back((char*)x.c_str());
  int rc = 4;
  try { rc = engine_main((int)av.size(), av.data()); }
  catch (const std::exception& e) { printf("[graphiti] engine run aborted: %s\n", e.what()); fflush(stdout); }
  catch (...) { printf("[graphiti] engine run aborted: unknown C++ exception\n"); fflush(stdout); }
  printf("[graphiti] engine run done rc=%d steps=%s n=%d budget=%.2f\n", rc, steps, n, budget_s); fflush(stdout);
  sketch_ui_attach(nullptr);
  perf_report(strcmp(steps, "0") == 0 ? "prepare" : "solve", rc, n, budget_s, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - perf_t0).count(), perf_th0);
  return rc;
}
