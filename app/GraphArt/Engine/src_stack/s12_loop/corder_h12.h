#pragma once
#include <cstring>
#include "../s02_nf/sk_nf_h2b.hpp"

namespace h11 {
constexpr int MAXN = 512, MAXQ = 64, MAXM = MAXN + MAXQ + 4;

struct Tab { const int* opnd; const int* OP; const int* CA; const int* CB; };
struct Seq { int* sn; int* sa; int* sb; };

static inline int emit_root(const Tab& T, unsigned char* dn, const Seq& S, int* rd, int* lr, int m, int r) {
  const int m0 = m;
  int stx[MAXN], sts[MAXN]; int sp = 0;
  stx[0] = r; sts[0] = 0; sp = 1;
  while (sp > 0) {
    const int t = sp - 1, x = stx[t];
    if (sts[t] == 0) {
      if (T.opnd[x] >= 0 || dn[x]) { --sp; continue; }
      sts[t] = 1; stx[sp] = T.CA[x]; sts[sp] = 0; ++sp; continue;
    }
    if (sts[t] == 1 && T.OP[x] >= 3 && T.OP[x] <= 6) { sts[t] = 2; stx[sp] = T.CB[x]; sts[sp] = 0; ++sp; continue; }
    const int A = T.CA[x], B = (sts[t] == 2) ? T.CB[x] : -1;
    S.sn[m] = x; S.sa[m] = A; S.sb[m] = B; rd[A]++; lr[A] = m; if (B >= 0) { rd[B]++; lr[B] = m; }
    ++m; dn[x] = 1; --sp;
  }
  if (m == m0) { S.sn[m] = -1; S.sa[m] = r; S.sb[m] = -1; rd[r]++; lr[r] = m; ++m; }
  return m;
}

static inline bool corder(const skn::View& V, const int* rep, const int* opnd, int n, int* out) {
  const int Q = V.Q;
  if (Q < 2 || Q > MAXQ || n < 1 || n > MAXN) return false;
  int OP[MAXN], CA[MAXN], CB[MAXN], rr[MAXQ];
  for (int q = 0; q < Q; q++) rr[q] = rep[V.root(q)];
  for (int x = 0; x < n; x++) {
    if (rep[x] != x || opnd[x] >= 0) { OP[x] = -1; CA[x] = -1; CB[x] = -1; continue; }
    const int o = V.op(x); OP[x] = o; CA[x] = rep[V.a(x)]; CB[x] = (o >= 3 && o <= 6) ? rep[V.b(x)] : -1;
  }
  const Tab T{opnd, OP, CA, CB};
  int sn[MAXM], sa[MAXM], sb[MAXM]; const Seq S{sn, sa, sb};
  unsigned char dn[MAXN];
  int heldst[MAXN];
  int rdP[MAXN], lrP[MAXN], rd[MAXN], lr[MAXN];
  for (int x = 0; x < n; x++) { dn[x] = 0; heldst[x] = 0; rdP[x] = 0; lrP[x] = -1; rd[x] = 0; lr[x] = -1; }
  int stamp = 0;
  int rem[MAXQ]; int nr = Q; for (int q = 0; q < Q; q++) rem[q] = q;
  int mP = 0;
  int no = 0;
  while (nr > 0) {
    int bc = 1 << 30, bi = 0;
    for (int i = 0; i < nr; i++) {
      int m = mP;
      m = emit_root(T, dn, S, rd, lr, m, rr[rem[i]]);
      for (int j = 0; j < nr; j++) if (j != i) m = emit_root(T, dn, S, rd, lr, m, rr[rem[j]]);
      if (m < bc) {
        ++stamp;
        int live = 0, mx = 0; bool cut = false;
        for (int k = 0; k < m; k++) {
          const int zn = sn[k], za = sa[k], zb = sb[k];
          int ch = 0;
          if (zn >= 0 && rd[zn] > 0) {
            const bool chained = (k + 1 < m) && (((sa[k + 1] == zn) + (sb[k + 1] == zn)) == rd[zn]);
            if (!chained) ch = 1;
          }
          if (opnd[za] < 0 && heldst[za] == stamp && lr[za] == k) { heldst[za] = 0; --live; }
          if (zb >= 0 && zb != za && opnd[zb] < 0 && heldst[zb] == stamp && lr[zb] == k) { heldst[zb] = 0; --live; }
          if (ch) {
            heldst[zn] = stamp;
            if (++live > mx) { mx = live; if (mx * 1000 + m >= bc) { cut = true; break; } }
          }
        }
        if (!cut) { const int c = mx * 1000 + m; if (c < bc) { bc = c; bi = i; } }
      }
      for (int k = mP; k < m; k++) {
        if (sn[k] >= 0) dn[sn[k]] = 0;
        rd[sa[k]] = rdP[sa[k]]; lr[sa[k]] = lrP[sa[k]];
        if (sb[k] >= 0) { rd[sb[k]] = rdP[sb[k]]; lr[sb[k]] = lrP[sb[k]]; }
      }
    }
    out[no++] = rem[bi];
    const int m2 = emit_root(T, dn, S, rdP, lrP, mP, rr[rem[bi]]);
    for (int k = mP; k < m2; k++) { rd[sa[k]] = rdP[sa[k]]; lr[sa[k]] = lrP[sa[k]]; if (sb[k] >= 0) { rd[sb[k]] = rdP[sb[k]]; lr[sb[k]] = lrP[sb[k]]; } }
    mP = m2;
    for (int j = bi; j + 1 < nr; j++) rem[j] = rem[j + 1];
    --nr;
  }
  return true;
}
}
