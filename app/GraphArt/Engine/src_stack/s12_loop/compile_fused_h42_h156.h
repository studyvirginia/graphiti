#pragma once
#include "compile_fused_h38_h156.h"
#include "prepare_h42.h"

namespace f42 {
using h14::MAXN; using h14::MAXQ; using h14::MAXM; using h14::NWD;
using f38::Grid; using f38::rec_fields; using f38::compile_composed; using f38::zgconv14_ref;
using f38::Scratch;

static bool compile_fused_stage(const View& V, Tree& T, cmi6::TreeEnc& te, f38::Grid& g, int upto, int* path) {
  Scratch& scr = f38::scr;
  const int Q = V.Q, n = V.n();
  if (path) *path = 0;
  if (Q > 17 || n > MAXN || (Q > 1 && n < 1)) { if (path) *path = 1; return compile_composed(V, T, te, g); }
  if (Q < 0) return false;
  T.Q = V.Q; T.nth = V.nth; T.nS = V.nS; T.prog.clear(); T.sp.clear(); T.th.clear();
  if (V.nS > NSV || V.nth > 32) return false;
  T.sp.reserve((size_t)V.nS);
  for (int i = 0; i < V.nS; i++) { int o = V.op(i), a = V.a(i), b = V.b(i); if (a > 31 || b > 31) return false; T.sp.push_back((uint16_t)(o | (a << 6) | (b << 11))); }

  h42::PrepBuf P;
  h42::prepare(V, P);
  const int* opnd = P.opnd;
  int ord[MAXQ];
  for (int q = 0; q < Q; q++) ord[q] = q;
  if (Q > 1 && !h42::order(h14::tab, ord)) { if (path) *path = 1; return compile_composed(V, T, te, g); }
  if (upto <= 1) return true;

  h14::Tab& C = h14::tab;
  const h14::Seq S{scr.sn, scr.sa, scr.sb};
  for (int x = 0; x <= n; x++) { scr.reads[x] = 0; scr.lastread[x] = -1; }
  uint64_t D[NWD]; std::memset(D, 0, sizeof D);
  int m = 0;
  for (int qi = 0; qi < Q; qi++) { m = h14::emit_q(C, ord[qi], D, S, scr.reads, scr.lastread, m); scr.segend[qi] = m; }
  scr.sn[m] = n; scr.sa[m] = n; scr.sb[m] = n;
  for (int k = 0; k < m; k++) { if (scr.sn[k] == n) scr.sn[k] = -1; if (scr.sb[k] == n) scr.sb[k] = -1; }
  for (int x = 0; x < n; x++) { scr.pos[x] = -1; scr.chain[x] = 0; }
  for (int k = 0; k < m; k++) if (scr.sn[k] >= 0) scr.pos[scr.sn[k]] = k;
  for (int k = 0; k + 1 < m; k++) {
    const int x = scr.sn[k];
    if (x < 0 || scr.reads[x] == 0) continue;
    const int r1 = (scr.sa[k + 1] == x) + (scr.sb[k + 1] == x);
    if (r1 == scr.reads[x]) scr.chain[x] = 1;
  }
  for (int x = 0; x < n; x++) scr.slot[x] = -1;
  uint16_t capt[18]; for (int q = 0; q < 18; q++) capt[q] = 1022;
  { int q = 0; for (int k = 0; k < m; k++) { scr.qof[k] = q; if (k + 1 == scr.segend[q]) q++; } }
  auto code = [&](int x, int k) -> int {
    if (opnd[x] >= 0) return opnd[x];
    if (scr.chain[x] && scr.pos[x] == k - 1) return 1022;
    return scr.slot[x];
  };
  int nsl = 0, nf = 0, pc = 0;
  bool decFail = false;
  scr.newend[0] = 0;
  for (int k = 0; k < m; k++) {
    const int node = scr.sn[k], A0 = scr.sa[k], B0 = scr.sb[k];
    const int A = code(A0, k), B = B0 >= 0 ? code(B0, k) : 0;
    if (A < 0 || B < 0) return false;
    if (opnd[A0] < 0 && scr.slot[A0] >= 0 && scr.lastread[A0] == k) { scr.freeS[nf++] = scr.slot[A0]; scr.slot[A0] = -2; }
    if (B0 >= 0 && B0 != A0 && opnd[B0] < 0 && scr.slot[B0] >= 0 && scr.lastread[B0] == k) { scr.freeS[nf++] = scr.slot[B0]; scr.slot[B0] = -2; }
    int d = 63;
    if (node >= 0 && scr.reads[node] > 0 && !scr.chain[node]) {
      d = nf > 0 ? scr.freeS[--nf] : nsl++;
      scr.slot[node] = d;
    }
    if (node < 0) { capt[scr.qof[k]] = (uint16_t)A; scr.newend[k + 1] = pc; continue; }
    const uint32_t op = (uint32_t)C.OP[node];
    scr.pk[pc] = op | (uint32_t)d << 6 | (uint32_t)A << 12 | (uint32_t)B << 22;
    if (!rec_fields(op, (uint32_t)d, (uint32_t)A, (uint32_t)B, scr.rec[pc])) decFail = true;
    ++pc;
    scr.newend[k + 1] = pc;
  }
  if (nsl > MAXSLOT) return false;
  T.prog.assign(scr.pk, scr.pk + pc);
  for (int q = 0; q < 36; q++) T.segs[q] = 0;
  for (int q = 0; q < Q; q++) T.segs[q + 1] = (uint16_t)scr.newend[scr.segend[q]];
  for (int q = 0; q < 18; q++) T.segs[18 + q] = 1022;
  for (int q = 0; q < Q; q++) T.segs[18 + q] = capt[q];
  T.nslot = nsl;
  T.th.reserve((size_t)V.nth);
  for (int s = 0; s < V.nth; s++) T.th.push_back(V.th(s));
  if (upto <= 2) return true;

  if (decFail) {
    if (path) *path = 2;
    cmi6::encode_tree_i6(T, te);
    g.ok = zgconv14_ref(T, g.p4, g.sg, g.hw) ? 1 : 0;
    return true;
  }

  const int SS = T.nth == 0 ? 2 : 16;
  const enc34::ScanE es = enc34::make_scan(T, SS);
  const bool gin = T.nth >= 1 && T.nth <= GN_TH_MAX;
  const int CD = T.nth <= 4 ? 1 + (T.nth + 1) / 2 : 3;
  const uint32_t NSC = (uint32_t)T.nS + 2;
  const uint32_t SB = (uint32_t)T.nslot * (uint32_t)CD * 48;
  const bool gok = gin && SB < 65536;
  const enc34::GnE eg = enc34::make_gn(T, CD);
  const enc34::Dec* R = scr.rec;
  uint32_t* os = scr.ob; uint32_t* og = scr.ob2;
  size_t k = 0, out = 0; int qp = 0;
  const int Qt = T.Q;
  while (k < (size_t)pc) {
    while (qp <= Qt && (size_t)T.segs[qp] == k) scr.segmap[qp++] = (uint32_t)out;
    const enc34::Dec& a = R[k];
    bool fuse = false;
    if (k + 1 < (size_t)pc && a.zo == 0 && a.d63) {
      const enc34::Dec& b = R[k + 1];
      fuse = b.zo >= 7 && b.zo <= 15 && b.A == 1022 && !(qp <= Qt && (size_t)T.segs[qp] == k + 1);
    }
    const enc34::Dec& z = fuse ? R[k + 1] : a;
    uint32_t* s = os + 4 * out;
    s[0] = z.zo; s[1] = z.d63 ? 0xFFFFFFFFu : es.dscale(z.d); s[2] = es.enc(a.A); s[3] = es.enc(a.B);
    if (gok) {
      uint32_t* gg = og + 4 * out;
      gg[0] = z.zo; gg[1] = z.d63 ? 0xFFFFFFFFu : eg.dscale(z.d); gg[2] = eg.enc(a.A); gg[3] = eg.enc(a.B);
    }
    k += fuse ? 2 : 1;
    ++out;
  }
  while (qp <= Qt && (size_t)T.segs[qp] == (size_t)pc) scr.segmap[qp++] = (uint32_t)out;
  if (qp != Qt + 1) {
    if (path) *path = 1;
    return compile_composed(V, T, te, g);
  }
  std::memset(te.sg, 0, sizeof te.sg); std::memset(te.gg, 0, sizeof te.gg);
  te.sp4.assign(os, os + 4 * out);
  enc34::segs_emit(te.sg, T, scr.segmap, es);
  te.s_ok = 1;
  if (gok) {
    te.gp4.assign(og, og + 4 * out);
    enc34::segs_emit(te.gg, T, scr.segmap, eg);
    te.ghw = NSC | (SB << 8) | ((uint32_t)T.nth << 24);
    te.gper = (uint32_t)((size_t)SB + 2 * (size_t)CD * NSC);
    te.g_ok = 1;
  }
  if (upto <= 3) return true;

  const uint32_t SB1 = (uint32_t)T.nslot * 48u;
  const uint32_t NSC1 = (uint32_t)T.nS + 2;
  g.ok = SB1 < 65536 ? 1 : 0;
  for (int q = 0; q < 18; q++) g.sg[q] = T.segs[q];
  const enc34::GnE e1 = enc34::make_gn(T, 1);
  for (int q = 0; q < 18; q++) g.sg[18 + q] = e1.enc(T.segs[18 + q]);
  if (g.ok) {
    g.hw = NSC1 | (SB1 << 8) | ((uint32_t)T.nth << 24);
    for (size_t i = 0; i < (size_t)pc; i++) {
      uint32_t* q4 = scr.gb + 4 * i;
      q4[0] = R[i].zo; q4[1] = R[i].d63 ? 0xFFFFFFFFu : e1.dscale(R[i].d); q4[2] = e1.enc(R[i].A); q4[3] = e1.enc(R[i].B);
    }
    g.p4.assign(scr.gb, scr.gb + 4 * (size_t)pc);
  } else {
    g.p4.clear();
  }
  if (path) *path = 0;
  return true;
}

static bool compile_fused(const View& V, Tree& T, cmi6::TreeEnc& te, f38::Grid& g, int* path = nullptr) {
  return f42::compile_fused_stage(V, T, te, g, 4, path);
}
}
