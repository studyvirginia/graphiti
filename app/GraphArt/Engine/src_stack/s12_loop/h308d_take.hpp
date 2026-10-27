#pragma once
#include <cstddef>
#include <vector>
#include <atomic>
namespace h308d {
static const size_t D = 2;
static inline void swap_vectors(Tree& t, cmi6::TreeEnc& e, Tree& pt, cmi6::TreeEnc& pe) {
  t.prog.swap(pt.prog); t.sp.swap(pt.sp); t.th.swap(pt.th); t.rec.swap(pt.rec);
  e.sp4.swap(pe.sp4); e.gp4.swap(pe.gp4); e.tk.swap(pe.tk);
  t.prog.clear(); t.sp.clear(); t.th.clear(); t.rec.clear();
#ifndef H308D_NEG_NOCLEAR
  e.sp4.clear(); e.gp4.clear(); e.tk.clear();
#endif
}
static inline void prefetch_tree(const Tree* t) { __builtin_prefetch((const char*)&t->prog, 1, 3); __builtin_prefetch((const char*)(&t->rec + 1) - 1, 1, 3); }
static inline void prefetch_enc(const cmi6::TreeEnc* e) { __builtin_prefetch((const char*)&e->sp4, 1, 3); __builtin_prefetch((const char*)(&e->gp4 + 1) - 1, 1, 3); }
struct Walk { size_t base = 0, cur = 0, n = 0; int bp = 0; };
static inline void begin(Walk& s, size_t nb, int bp = 0) { s.n = nb; s.cur = 0; s.bp = bp; }
static inline void pf(const std::vector<Tree>& PT, const std::vector<cmi6::TreeEnc>& PE, size_t r, Tree* t, cmi6::TreeEnc* e) {
  if (r < PT.size()) { prefetch_tree(&PT[r]); prefetch_enc(&PE[r]); } prefetch_tree(t); prefetch_enc(e);
}
static inline void take(Walk& s, std::atomic<size_t>& ctr, std::vector<Tree>& PT, std::vector<cmi6::TreeEnc>& PE, Tree& t, cmi6::TreeEnc& e) {
  const size_t j = s.cur++;
  if (j == 0) { s.base = ctr.fetch_add(s.n, std::memory_order_relaxed); for (size_t q = 0; q < D && q < s.n; q++) pf(PT, PE, s.base + q, &t + q, &e + q); }
  const size_t r = s.base + j;
  if (j + D < s.n) pf(PT, PE, r + D, &t + D, &e + D);
  if (r < PT.size()) swap_vectors(t, e, PT[r], PE[r]);
  if (s.bp && j + 1 < s.n && r + 1 < PT.size()) {
    const Tree& q = PT[r + 1]; const cmi6::TreeEnc& qe = PE[r + 1];
    __builtin_prefetch(q.prog.data(), 1, 3); __builtin_prefetch(q.sp.data(), 1, 3); __builtin_prefetch(q.th.data(), 1, 3);
    __builtin_prefetch(qe.sp4.data(), 1, 3); __builtin_prefetch(qe.gp4.data(), 1, 3); }
}
}
