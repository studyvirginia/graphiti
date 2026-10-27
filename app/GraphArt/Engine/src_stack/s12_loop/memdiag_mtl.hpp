#pragma once
#import <Metal/Metal.h>
#include <vector>
namespace memdiag {
struct MEnt { __weak id<MTLBuffer> w; size_t n; uintptr_t pc[3]; bool priv; bool nocopy; };
static std::vector<MEnt> g_ment;
static void metal_note(id<MTLBuffer> b, size_t n, bool priv, bool nocopy = false) {
  if (!tab || !b) return;
  if (busy_get()) return; busy_set(1);
  void* fr[10]; const int nf = backtrace(fr, 10); MEnt e; e.w = b; e.n = n; e.priv = priv; e.nocopy = nocopy; e.pc[0] = e.pc[1] = e.pc[2] = 0; int k = 0, skip = 2;
  for (int i = 0; i < nf && k < 3; i++) { const uintptr_t a = (uintptr_t)fr[i] & 0x0000007fffffffffULL; if (a >= img_lo && a < img_hi) { if (skip > 0) { skip--; continue; } e.pc[k++] = a - slide; } }
  g_ment.push_back(e); busy_set(0);
}
static void metal_report(int st, double dev_alloc_mb) {
  busy_set(1);
  struct G { uintptr_t pc[3]; double mb; int n; bool priv, nocopy; };
  std::vector<G> g; double tot = 0, totp = 0, totnc = 0;
  for (auto& e : g_ment) { id<MTLBuffer> b = e.w; if (!b) continue; const double mb = (double)e.n / 1048576.0; tot += mb; if (e.priv) totp += mb; if (e.nocopy) totnc += mb;
    bool f = false; for (auto& x : g) if (x.pc[0] == e.pc[0] && x.pc[1] == e.pc[1] && x.pc[2] == e.pc[2] && x.priv == e.priv && x.nocopy == e.nocopy) { x.mb += mb; x.n++; f = true; break; }
    if (!f) g.push_back({{e.pc[0], e.pc[1], e.pc[2]}, mb, 1, e.priv, e.nocopy}); }
  std::sort(g.begin(), g.end(), [](const G& a, const G& b) { return a.mb > b.mb; });
  fprintf(stderr, "MEMDIAG st=%d metal live buffers %.1f MB (private %.1f, nocopy-host %.1f) device currentAllocatedSize %.1f MB\n", st, tot, totp, totnc, dev_alloc_mb);
  for (size_t i = 0; i < g.size() && i < 40; i++) { if (g[i].mb < 1.0) break;
    fprintf(stderr, "MEMDIAG st=%d msite %s live %.1f MB n %d pc 0x%lx 0x%lx 0x%lx\n", st, g[i].nocopy ? "nocopy" : g[i].priv ? "private" : "shared", g[i].mb, g[i].n, (unsigned long)g[i].pc[0], (unsigned long)g[i].pc[1], (unsigned long)g[i].pc[2]); }
  busy_set(0);
}
}
