#import <Foundation/Foundation.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>

extern int engine_main(int argc, char** argv);
extern "C" { void* sketch_ui_ring_create(void); void sketch_ui_attach(void*); void sketch_ui_set_budget(double); int sketch_ui_pop(void*, void*); unsigned long sketch_ui_burst_bytes(void); }
namespace stb { constexpr int PAIRS = 9, PTS = 120, MAXLEN = 256;
struct Tree { uint8_t len; uint8_t op[MAXLEN]; int16_t cid[MAXLEN]; float c[MAXLEN]; };
struct Fit { uint8_t kind, ovf; uint8_t q[2], nth[2]; uint8_t drv[2]; uint8_t pad[3]; float u[10]; float th[2][32]; float a[2][16]; };
struct Pair { Tree x, y; float score; uint16_t cx; uint8_t champion; uint8_t fitted; float xs[PTS], ys[PTS]; Fit fit; };
struct Burst { uint32_t seq; int32_t step; int32_t n_pairs; int32_t final_pick; double wall_ms; int64_t candidates; float t_lo, t_hi; Pair pairs[PAIRS]; }; }

static void setenv_if_unset(const char* k, const char* v) { if (!getenv(k)) setenv(k, v, 1); }

int main(int argc, char** argv) {
  @autoreleasepool {
    if (const char* ef = getenv("APP_ENV_FILE")) {
      std::ifstream f(ef); std::string all, tok; std::getline(f, all);
      std::istringstream is(all);
      while (is >> tok) { size_t e = tok.find('='); if (e != std::string::npos) setenv_if_unset(tok.substr(0, e).c_str(), tok.substr(e + 1).c_str()); }
    }
    setenv_if_unset("SICIS", "1"); setenv("FINAL_PICK", "1", 1); setenv("SINGLE_FAST", "1", 1); setenv("S12_FREET", "1", 1);
    setenv_if_unset("H242D_FIX", "0.35");
    setenv("H930_MEMFIX", "1", 1); setenv_if_unset("FJN_CISFUSE_HOST", "1");
    setenv("FREET_QSPLIT", "1", 1); setenv_if_unset("H117B", "1"); setenv_if_unset("SINGLE_FAST_PIN", "1");
    setenv_if_unset("FJN_CISFUSE_CHECK", "256"); setenv_if_unset("PSO_PAR", "1"); setenv_if_unset("H178", "0");
    setvbuf(stdout, nullptr, _IOLBF, 0);
    if (argc < 9) { fprintf(stderr, "usage: %s cands kdir outdir prims label steps threads check\n", argv[0]); return 2; }
    std::vector<char*> av = { argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8] };
    if (getenv("RING_TEST")) {
      void* ring = sketch_ui_ring_create(); sketch_ui_attach(ring); sketch_ui_set_budget(0);
      static std::atomic<bool> stop{false}; const auto t0 = std::chrono::steady_clock::now();
      std::thread drain([&] { std::vector<unsigned char> buf(sketch_ui_burst_bytes());
        for (;;) { const bool done = stop.load();
          if (sketch_ui_pop(ring, buf.data()) == 1) { const stb::Burst* b = (const stb::Burst*)buf.data();
            printf("RING_BURST t_ms %.1f step %d final %d n_pairs %d champion: cx %d score %.6g champ_flag %d fitted %d | others cx:", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(), b->step, b->final_pick, b->n_pairs, b->pairs[0].cx, b->pairs[0].score, b->pairs[0].champion, b->pairs[0].fitted);
            for (int p = 1; p < b->n_pairs; p++) printf(" %d(%.4g)", b->pairs[p].cx, b->pairs[p].score); printf("\n"); }
          else if (done) break; else std::this_thread::sleep_for(std::chrono::milliseconds(1)); } });
      const int rc = engine_main((int)av.size(), av.data()); stop.store(true); drain.join(); return rc; }
    return engine_main((int)av.size(), av.data());
  }
}
