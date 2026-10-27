#pragma once
#include "loop_common_v2.h"
#include "corder_h12.h"

static bool compile_tree_corder_h12(const View& V, Tree& T) {
  T.Q = V.Q; T.nth = V.nth; T.nS = V.nS; T.prog.clear(); T.sp.clear(); T.th.clear();
  if (V.nS > NSV || V.nth > 32) return false;
  for (int i = 0; i < V.nS; i++) { int o = V.op(i), a = V.a(i), b = V.b(i); if (a > 31 || b > 31) return false; T.sp.push_back((uint16_t)(o | (a << 6) | (b << 11))); }
  const int n = V.n(); std::vector<int> opnd(n, -1);
  for (int i = 0; i < n; i++) { int o = V.op(i); if (i < V.nS) opnd[i] = 64 + i; else if (o == 2) opnd[i] = 900; else if (o == 16) opnd[i] = 901; }
  std::vector<int> rep(n); for (int i = 0; i < n; i++) rep[i] = i;
  { std::map<std::tuple<int,int,int>, int> hc;
    for (int i = V.nS; i < n; i++) { int o = V.op(i); if (o == 2 || o == 16) continue; int a = rep[V.a(i)], b = (o >= 3 && o <= 6) ? rep[V.b(i)] : -1;
      if ((o == 3 || o == 5) && a > b) std::swap(a, b); auto key = std::make_tuple(o, a, b); auto it = hc.find(key); if (it == hc.end()) hc[key] = i; else rep[i] = it->second; } }
  auto CA = [&](int x) { return rep[V.a(x)]; }; auto CB = [&](int x) { return rep[V.b(x)]; };
  struct I { int node, op, A, B; };
  std::vector<int> corder(V.Q); for (int q = 0; q < V.Q; q++) corder[q] = q;
  auto peak = [&](const std::vector<int>& ord) {
    std::vector<I> in2; std::vector<char> dn(n, 0);
    std::function<void(int)> df = [&](int x) { x = rep[x]; if (opnd[x] >= 0 || dn[x]) return; int o = V.op(x);
      if (o >= 3 && o <= 6) { df(CA(x)); df(CB(x)); in2.push_back({x, o, CA(x), CB(x)}); } else { df(CA(x)); in2.push_back({x, o, CA(x), -1}); } dn[x] = 1; };
    for (int q : ord) { int r = rep[V.root(q)]; size_t b0 = in2.size(); df(r); if (in2.size() == b0) in2.push_back({-1, 0, r, -1}); }
    const int mm = (int)in2.size(); std::vector<int> rd(n, 0), lr(n, -1); for (int k = 0; k < mm; k++) { rd[in2[k].A]++; lr[in2[k].A] = k; if (in2[k].B >= 0) { rd[in2[k].B]++; lr[in2[k].B] = k; } }
    int live = 0, mx = 0; std::vector<char> held(n, 0);
    for (int k = 0; k < mm; k++) { auto& z = in2[k]; int ch = 0;
      if (z.node >= 0 && rd[z.node] > 0) { bool chained = (k + 1 < mm) && ((in2[k + 1].A == z.node) + (in2[k + 1].B == z.node) == rd[z.node]); if (!chained) ch = 1; }
      if (opnd[z.A] < 0 && held[z.A] && lr[z.A] == k) { held[z.A] = 0; live--; } if (z.B >= 0 && z.B != z.A && opnd[z.B] < 0 && held[z.B] && lr[z.B] == k) { held[z.B] = 0; live--; }
      if (ch) { held[z.node] = 1; live++; mx = std::max(mx, live); } }
    return mx * 1000 + mm; };
  if (V.Q > 1) {
    int ob[h11::MAXQ];
    if (h11::corder(V, rep.data(), opnd.data(), n, ob)) corder.assign(ob, ob + V.Q);
    else {
    if (V.Q > 1) { std::vector<int> ord, rem = corder;
      while (!rem.empty()) { int bi = 0, bc = 1 << 30; for (size_t i = 0; i < rem.size(); i++) { std::vector<int> tr = ord; tr.push_back(rem[i]); for (size_t j = 0; j < rem.size(); j++) if (j != i) tr.push_back(rem[j]); int c = peak(tr); if (c < bc) { bc = c; bi = (int)i; } }
        ord.push_back(rem[bi]); rem.erase(rem.begin() + bi); }
      corder = ord; }
    }
  }
  std::vector<I> ins; std::vector<char> done(n, 0); std::vector<int> segend;
  std::function<void(int)> dfs = [&](int x) { x = rep[x]; if (opnd[x] >= 0 || done[x]) return; int o = V.op(x);
    if (o >= 3 && o <= 6) { dfs(CA(x)); dfs(CB(x)); ins.push_back({x, o, CA(x), CB(x)}); } else { dfs(CA(x)); ins.push_back({x, o, CA(x), -1}); } done[x] = 1; };
  for (int qi = 0; qi < V.Q; qi++) { int r = rep[V.root(corder[qi])]; size_t before = ins.size(); dfs(r);
    if (ins.size() == before) ins.push_back({-1, 0, r, -1}); segend.push_back((int)ins.size()); }
  const int m = (int)ins.size(); std::vector<int> reads(n, 0), lastread(n, -1);
  for (int k = 0; k < m; k++) { auto& z = ins[k]; reads[z.A]++; lastread[z.A] = k; if (z.B >= 0) { reads[z.B]++; lastread[z.B] = k; } }
  std::vector<int> pos(n, -1); for (int k = 0; k < m; k++) if (ins[k].node >= 0) pos[ins[k].node] = k;
  std::vector<char> chain(n, 0);
  for (int k = 0; k + 1 < m; k++) { int x = ins[k].node; if (x < 0 || reads[x] == 0) continue; auto& z = ins[k + 1]; int r1 = (z.A == x) + (z.B == x); if (r1 == reads[x]) chain[x] = 1; }
  std::vector<int> slot(n, -1), freeS; int nsl = 0;
  std::vector<uint16_t> capt(17, 1022); std::vector<int> qof(m, -1), newend(m + 1, 0); { int q = 0; for (int k = 0; k < m; k++) { qof[k] = q; if (k + 1 == segend[q]) q++; } }
  auto code = [&](int x, int k) { if (opnd[x] >= 0) return opnd[x]; if (chain[x] && pos[x] == k - 1) return 1022; return slot[x]; };
  T.prog.clear();
  for (int k = 0; k < m; k++) { auto& z = ins[k]; int A = code(z.A, k), B = z.B >= 0 ? code(z.B, k) : 0; if (A < 0 || B < 0) return false;
    if (opnd[z.A] < 0 && slot[z.A] >= 0 && lastread[z.A] == k) { freeS.push_back(slot[z.A]); slot[z.A] = -2; }
    if (z.B >= 0 && z.B != z.A && opnd[z.B] < 0 && slot[z.B] >= 0 && lastread[z.B] == k) { freeS.push_back(slot[z.B]); slot[z.B] = -2; }
    int d = 63; if (z.node >= 0 && reads[z.node] > 0 && !chain[z.node]) { if (!freeS.empty()) { d = freeS.back(); freeS.pop_back(); } else d = nsl++; slot[z.node] = d; }
    if (z.node < 0) { capt[qof[k]] = (uint16_t)A; newend[k + 1] = (int)T.prog.size(); continue; }
    T.prog.push_back((uint32_t)z.op | (uint32_t)d << 6 | (uint32_t)A << 12 | (uint32_t)B << 22); newend[k + 1] = (int)T.prog.size(); }
  for (int q = 0; q < 36; q++) T.segs[q] = 0;
  T.segs[0] = 0; for (int q = 0; q < V.Q; q++) T.segs[q + 1] = (uint16_t)newend[segend[q]];
  for (int q = 0; q < 18; q++) T.segs[18 + q] = 1022;
  for (int q = 0; q < V.Q; q++) T.segs[18 + q] = capt[q];
  T.nslot = nsl; if (nsl > MAXSLOT) return false;
  for (int s = 0; s < V.nth; s++) T.th.push_back(V.th(s)); return true;
}

