#pragma once
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstddef>

namespace uib {

constexpr int RING = 64;
constexpr int PAIRS = 9;
constexpr int PTS = 120;
constexpr int GENLEN = 24;
constexpr int MAXLEN = 256;

struct Tree {
  uint8_t len;
  uint8_t op[MAXLEN];
  int16_t cid[MAXLEN];
  float c[MAXLEN];
};

struct Fit {
  uint8_t kind, ovf;
  uint8_t q[2], nth[2];
  uint8_t drv[2];
  uint8_t pad[3];
  float u[10];
  float th[2][32];
  float a[2][16];
};

struct Pair {
  Tree x, y;
  float score;
  uint16_t cx;
  uint8_t champion;
  uint8_t fitted;
  float xs[PTS], ys[PTS];
  Fit fit;
};

struct Burst {
  uint32_t seq;
  int32_t step;
  int32_t n_pairs;
  int32_t final_pick;
  double wall_ms;
  int64_t candidates;
  float t_lo, t_hi;
  Pair pairs[PAIRS];
};

struct Ring {
  Burst slot[RING];
  std::atomic<uint64_t> head{0};
  std::atomic<uint64_t> tail{0};
  std::atomic<uint64_t> dropped{0};
  uint32_t seq = 0;
};

inline Ring* g_ring = nullptr;
inline double g_budget_s = 0.0;
inline double g_t0 = 0.0;
inline double g_map[4] = {1, 0, 1, 0};
inline int g_xf_on = 0, g_xf_drv = 0;
inline double g_xf_a = 1, g_xf_b = 0;
inline double g_smap[4] = {1, 0, 1, 0};
inline volatile int g_stop = 0;

inline double op_eval(uint8_t o, double a, double b, double t, double c) {
  switch (o) {
    case 1: return c;
    case 2: return t;
    case 3: return a + b;
    case 4: return a - b;
    case 5: return a * b;
    case 6: return a / b;
    case 7: return std::sin(a);
    case 8: return std::cos(a);
    case 9: return std::exp(a);
    case 10: return std::log(std::fabs(a));
    case 11: return a * a;
    case 12: return a * a * a;
    case 13: return std::sqrt(std::fabs(a));
    case 14: return std::cbrt(a);
    case 15: return std::fabs(a);
    default: return NAN;
  }
}

inline double eval_tree(const Tree& g, double t) {
  double st[MAXLEN]; int sp = 0;
  for (int k = 0; k < g.len; k++) {
    const uint8_t o = g.op[k];
    if (o == 1 || o == 2) { st[sp++] = op_eval(o, 0, 0, t, g.c[k]); continue; }
    if (o >= 3 && o <= 6) { if (sp < 2) return NAN; const double r = st[--sp]; const double l = st[--sp]; st[sp++] = op_eval(o, l, r, t, 0); continue; }
    if (o >= 7 && o <= 15) { if (sp < 1) return NAN; const double a = st[--sp]; st[sp++] = op_eval(o, a, 0, t, 0); continue; }
    return NAN;
  }
  return sp == 1 ? st[0] : NAN;
}

inline bool to_xframe(Tree& T, int side) {
  if (!g_xf_on || T.len == 0 || T.len > MAXLEN - 1) return true;
  struct N { uint8_t op; int l, r; double c, A, B; bool af; };
  static thread_local N n[MAXLEN + 8]; int stk[MAXLEN + 8], sp = 0; uint8_t Op[MAXLEN + 8]; float Cc[MAXLEN + 8]; int L = T.len;
  std::memcpy(Op, T.op, T.len); std::memcpy(Cc, T.c, sizeof(float) * T.len);
  if (side == g_xf_drv) {
    Op[L] = 1; Cc[L++] = (float)(1.0 / g_xf_a); Op[L] = 5; Cc[L++] = 0; Op[L] = 1; Cc[L++] = (float)(-g_xf_b / g_xf_a); Op[L] = 3; Cc[L++] = 0; }
  for (int k = 0; k < L; k++) {
    const uint8_t o = Op[k]; N& q = n[k]; q = N{o, -1, -1, (double)Cc[k], 0, 0, false};
    if (o == 1) { q.af = true; q.B = Cc[k]; }
    else if (o == 2) { q.af = true; q.A = g_xf_a; q.B = g_xf_b; }
    else if (o >= 3 && o <= 6) { if (sp < 2) return false; q.r = stk[--sp]; q.l = stk[--sp]; const N& l = n[q.l]; const N& r = n[q.r];
      if (l.af && r.af) {
        if (o == 3 || o == 4) { const double sg = o == 3 ? 1.0 : -1.0; q.af = true; q.A = l.A + sg * r.A; q.B = l.B + sg * r.B; }
        else if (o == 5 && (l.A == 0 || r.A == 0)) { q.af = true; q.A = l.A * r.B + r.A * l.B; q.B = l.B * r.B; }
        else if (o == 6 && r.A == 0 && r.B != 0) { q.af = true; q.A = l.A / r.B; q.B = l.B / r.B; }
        if (q.af && !(std::isfinite(q.A) && std::isfinite(q.B) && std::isfinite((float)q.A) && std::isfinite((float)q.B))) q.af = false; } }
    else if (o >= 7 && o <= 15) { if (sp < 1) return false; q.l = stk[--sp]; }
    else return false;
    stk[sp++] = k;
  }
  if (sp != 1) return false;
  Tree R; R.len = 0; int len = 0; bool ovf = false;
  auto put = [&](uint8_t op, float c) { if (len >= 255) { ovf = true; return; } R.op[len] = op; R.c[len] = c; R.cid[len] = op == 1 ? (int16_t)(side * MAXLEN + len) : (int16_t)-1; len++; };
  struct Em { decltype(put)& put; const N* n; bool& ovf;
    void go(int i) { if (ovf) return; const N& q = n[i];
      if (q.af) { if (q.A == 0) { put(1, (float)q.B); return; }
        if (q.A == 1.0) put(2, 0); else { put(1, (float)q.A); put(2, 0); put(5, 0); }
        if (q.B != 0) { put(1, (float)q.B); put(3, 0); } return; }
      go(q.l); if (q.op <= 6) go(q.r); put(q.op, 0); } } em{put, n, ovf};
  em.go(sp ? stk[0] : 0);
  if (ovf) return false;
  R.len = (uint8_t)len; T = R; return true;
}

template <class Fill>
inline void push_burst(int step, double wall_ms, int64_t candidates, int n_pairs, const float* score,
                       const uint16_t* cx, Fill fill, int final_pick, float t_lo, float t_hi) {
  Ring* R = g_ring;
  if (!R) return;
  const uint64_t h = R->head.load(std::memory_order_relaxed);
  if (h - R->tail.load(std::memory_order_acquire) >= RING) { R->dropped.fetch_add(1, std::memory_order_relaxed); return; }
  Burst& B = R->slot[h % RING];
  B.seq = ++R->seq; B.step = step; B.n_pairs = n_pairs; B.final_pick = final_pick;
  B.wall_ms = wall_ms; B.candidates = candidates; B.t_lo = t_lo; B.t_hi = t_hi;
  const float dt = (t_hi - t_lo) / (PTS - 1);
  for (int p = 0; p < n_pairs && p < PAIRS; p++) {
    Pair& P = B.pairs[p];
    P.fitted = 0; std::memset(&P.fit, 0, sizeof P.fit);
    fill(p, P);
    P.score = score[p]; P.cx = cx[p]; P.champion = (p == 0);
    if (!P.fitted) for (int i = 0; i < PTS; i++) {
      const double t = t_lo + dt * i;
      P.xs[i] = (float)(g_smap[0] * eval_tree(P.x, t) + g_smap[1]);
      P.ys[i] = (float)(g_smap[2] * eval_tree(P.y, t) + g_smap[3]);
    }
    if (g_xf_on && !P.fit.pad[0]) { if (!to_xframe(P.x, 0) || !to_xframe(P.y, 1)) { P.fit.ovf = 2; static int w = 0; if (w++ < 5) fprintf(stderr, "YFX_FRAME: display tree too long after the x-frame rewrite (pair %d), shown in the engine t frame\n", p); } }
  }
  R->head.store(h + 1, std::memory_order_release);
}

template <class Cand>
inline void fill_from_cand(const Cand& g, Pair& P) {
  for (int t = 0; t < 2; t++) {
    Tree& T = t == 0 ? P.x : P.y;
    T.len = g.len[t];
    for (int k = 0; k < T.len && k < GENLEN; k++) {
      T.op[k] = g.op[t][k]; T.c[k] = g.c[t][k];
      T.cid[k] = (g.op[t][k] == 1) ? (int16_t)(t * MAXLEN + k) : (int16_t)-1;
    }
  }
}

inline bool pop(Ring* R, Burst& out) {
  const uint64_t t = R->tail.load(std::memory_order_relaxed);
  if (t == R->head.load(std::memory_order_acquire)) return false;
  std::memcpy(&out, &R->slot[t % RING], sizeof(Burst));
  R->tail.store(t + 1, std::memory_order_release);
  return true;
}

}

extern "C" {
void* sketch_ui_ring_create(void);
void sketch_ui_ring_destroy(void* ring);
void sketch_ui_attach(void* ring);
void sketch_ui_detach(void);
void sketch_ui_set_budget(double seconds);
void sketch_ui_map(double* out4);
void sketch_ui_layout(int* out16);
void sketch_ui_stop(void);
int sketch_ui_pop(void* ring, void* out);
unsigned long sketch_ui_burst_bytes(void);
unsigned long long sketch_ui_dropped(void* ring);
}

extern "C" {
void* sketch_ui_ring_create(void) { return new uib::Ring(); }
void sketch_ui_ring_destroy(void* ring) { delete static_cast<uib::Ring*>(ring); }
void sketch_ui_attach(void* ring) { uib::g_ring = static_cast<uib::Ring*>(ring); }
void sketch_ui_detach(void) { uib::g_ring = nullptr; }
void sketch_ui_set_budget(double seconds) { uib::g_budget_s = seconds; uib::g_stop = 0; }
void sketch_ui_stop(void) { uib::g_stop = 1; }
void sketch_ui_map(double* out4) { for (int i = 0; i < 4; i++) out4[i] = uib::g_map[i]; }
void sketch_ui_layout(int* o) {
  using uib::Burst; using uib::Pair; using uib::Fit; using uib::Tree;
  o[0] = (int)sizeof(Burst); o[1] = (int)offsetof(Burst, pairs); o[2] = (int)sizeof(Pair); o[3] = (int)sizeof(Tree);
  o[4] = (int)offsetof(Tree, op); o[5] = (int)offsetof(Tree, cid); o[6] = (int)offsetof(Tree, c); o[7] = (int)offsetof(Pair, score);
  o[8] = (int)offsetof(Pair, cx); o[9] = (int)offsetof(Pair, champion); o[10] = (int)offsetof(Pair, fitted); o[11] = (int)offsetof(Pair, xs);
  o[12] = (int)offsetof(Pair, ys); o[13] = (int)offsetof(Pair, fit); o[14] = (int)offsetof(Pair, y); o[15] = (int)sizeof(Fit);
}
int sketch_ui_pop(void* ring, void* out) { return uib::pop(static_cast<uib::Ring*>(ring), *static_cast<uib::Burst*>(out)) ? 1 : 0; }
unsigned long sketch_ui_burst_bytes(void) { return sizeof(uib::Burst); }
unsigned long long sketch_ui_dropped(void* ring) { return static_cast<uib::Ring*>(ring)->dropped.load(); }
}
