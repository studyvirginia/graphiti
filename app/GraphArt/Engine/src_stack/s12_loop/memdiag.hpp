#pragma once
#include <execinfo.h>
#include <mach-o/dyld.h>
#include <malloc/malloc.h>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <vector>
extern "C" { extern void (*malloc_logger)(uint32_t, uintptr_t, uintptr_t, uintptr_t, uintptr_t, uint32_t); }
namespace memdiag {
static const size_t MD_MIN = 1;
static const size_t NTAB = 1u << 23;
static const int NSITE = 8192;
struct Slot { uintptr_t p; uint32_t sz; uint32_t site; };
static Slot* tab = nullptr;
struct Site { uintptr_t pc[3]; int64_t live; int64_t peak; int64_t nlive; uint32_t kind; int64_t res; int64_t swp; int64_t freed; int64_t freed_last; };
static Site sites[NSITE]; static int nsites = 0;
static std::atomic_flag lk = ATOMIC_FLAG_INIT;
#include <pthread.h>
static pthread_key_t bkey; static bool bkey_ok = false;
struct Busy { int v; Busy() : v(bkey_ok ? (int)(intptr_t)pthread_getspecific(bkey) : 1) {} };
static inline int busy_get() { return bkey_ok ? (int)(intptr_t)pthread_getspecific(bkey) : 1; }
static inline void busy_set(int v) { if (bkey_ok) pthread_setspecific(bkey, (void*)(intptr_t)v); }
static uintptr_t img_lo = 0, img_hi = 0, slide = 0;
static inline void lock() { while (lk.test_and_set(std::memory_order_acquire)) {} }
static inline void unlock() { lk.clear(std::memory_order_release); }
static int site_of(uintptr_t* fr, int nf, uint32_t kind) {
  uintptr_t key[3] = {0, 0, 0}; int k = 0;
  int skip = 2;
  for (int i = 0; i < nf && k < 3; i++) { const uintptr_t a = fr[i] & 0x0000007fffffffffULL; if (a >= img_lo && a < img_hi) { if (skip > 0) { skip--; continue; } key[k++] = a - slide; } }
  if (k == 0) { key[0] = 1; }
  uint64_t h = 1469598103934665603ULL; for (int i = 0; i < 3; i++) { h ^= key[i]; h *= 1099511628211ULL; } h ^= kind; h *= 1099511628211ULL;
  uint32_t s = (uint32_t)(h % NSITE);
  for (int probe = 0; probe < NSITE; probe++, s = (s + 1) % NSITE) {
    Site& S = sites[s];
    if (S.pc[0] == 0 && S.pc[1] == 0 && S.pc[2] == 0 && S.kind == 0 && S.nlive == 0 && S.peak == 0) { S.pc[0] = key[0]; S.pc[1] = key[1]; S.pc[2] = key[2]; S.kind = kind; nsites++; return (int)s; }
    if (S.pc[0] == key[0] && S.pc[1] == key[1] && S.pc[2] == key[2] && S.kind == kind) return (int)s;
  }
  return 0;
}
static void put(uintptr_t p, size_t sz, uint32_t kind) {
  void* fr[10]; const int nf = backtrace(fr, 10);
  lock();
  const int s = site_of((uintptr_t*)fr, nf, kind);
  size_t i = (p >> 4) * 11400714819323198485ULL >> 41 & (NTAB - 1);
  for (size_t n = 0; n < NTAB; n++, i = (i + 1) & (NTAB - 1)) { if (tab[i].p == 0 || tab[i].p == 1) { tab[i].p = p; tab[i].sz = (uint32_t)std::min<size_t>(sz, 0xffffffffu); tab[i].site = (uint32_t)s; break; } }
  sites[s].live += (int64_t)sz; sites[s].nlive++; if (sites[s].live > sites[s].peak) sites[s].peak = sites[s].live;
  unlock();
}
static bool take(uintptr_t p, size_t* szout = nullptr) {
  lock();
  size_t i = (p >> 4) * 11400714819323198485ULL >> 41 & (NTAB - 1);
  for (size_t n = 0; n < NTAB; n++, i = (i + 1) & (NTAB - 1)) {
    if (tab[i].p == 0) break;
    if (tab[i].p == p) { sites[tab[i].site].live -= tab[i].sz; sites[tab[i].site].nlive--; if (tab[i].sz >= 131072) sites[tab[i].site].freed += tab[i].sz; if (szout) *szout = tab[i].sz; tab[i].p = 1; unlock(); return true; }
  }
  unlock(); return false;
}
static void hook(uint32_t type, uintptr_t a1, uintptr_t a2, uintptr_t a3, uintptr_t res, uint32_t) {
  if (busy_get()) return; busy_set(1);
  const bool al = type & 2, de = type & 4, vm = type & (16 | 32);
  if (vm) {
    if ((type & 32) || (de && !al)) { take(a2); }
    else if (al && res && a2 >= MD_MIN) put(res, a2, 1);
  } else if (al && de) {
    if (a2) take(a2);
    if (res && a3 >= MD_MIN) put(res, a3, 0);
  } else if (al) { if (res && a2 >= MD_MIN) put(res, a2, 0); }
  else if (de) { if (a2) take(a2); }
  busy_set(0);
}
static void init() {
  pthread_key_create(&bkey, nullptr); bkey_ok = true;
  tab = (Slot*)calloc(NTAB, sizeof(Slot));
  const mach_header* mh = _dyld_get_image_header(0); slide = (uintptr_t)_dyld_get_image_vmaddr_slide(0);
  img_lo = (uintptr_t)mh; img_hi = img_lo + 0x4000000;
  void* w[4]; backtrace(w, 4);
  malloc_logger = hook;
}
#include <sys/mman.h>
#include <mach/mach.h>
#include <TargetConditionals.h>
#if !TARGET_OS_IPHONE
#include <mach/mach_vm.h>
#endif
static void residency() {
  for (int i = 0; i < NSITE; i++) { sites[i].res = 0; sites[i].swp = 0; }
  const size_t pg = 16384; static unsigned char vec[1 << 20];
  for (size_t i = 0; i < NTAB; i++) {
    Slot e = tab[i]; if (e.p <= 1) continue; Site& S = sites[e.site];
    if (e.sz < 131072) { S.res += e.sz; continue; }
    const uintptr_t a = e.p & ~(uintptr_t)(pg - 1); const size_t len = ((e.p + e.sz + pg - 1) & ~(uintptr_t)(pg - 1)) - a; const size_t np = len / pg;
    if (np > sizeof(vec)) { S.res += e.sz; continue; }
    int64_t r = 0; if (mincore((void*)a, len, (char*)vec) == 0) for (size_t k = 0; k < np; k++) r += (vec[k] & 1); S.res += r * (int64_t)pg;
#if !TARGET_OS_IPHONE
    mach_vm_address_t ra = a; mach_vm_size_t rs = 0; vm_region_extended_info_data_t ei; mach_msg_type_number_t c = VM_REGION_EXTENDED_INFO_COUNT; mach_port_t on = 0;
    if (mach_vm_region(mach_task_self(), &ra, &rs, VM_REGION_EXTENDED_INFO, (vm_region_info_t)&ei, &c, &on) == KERN_SUCCESS && ra <= a && rs <= len + pg) S.swp += (int64_t)ei.pages_swapped_out * (int64_t)pg;
#endif
  }
}
static void report(int st, const char* tag) {
  busy_set(1);
  static Site snap[NSITE]; lock(); residency(); memcpy(snap, sites, sizeof(sites)); unlock();
  int idx[NSITE]; int n = 0; for (int i = 0; i < NSITE; i++) if (snap[i].live > 0 || snap[i].peak > 0) idx[n++] = i;
  for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (snap[idx[b]].res > snap[idx[a]].res) std::swap(idx[a], idx[b]);
  double tot = 0, totv = 0; for (int i = 0; i < n; i++) { (snap[idx[i]].kind == 2 ? totv : tot) += (double)snap[idx[i]].res; }
  { malloc_statistics_t ms; malloc_zone_statistics(nullptr, &ms); fprintf(stderr, "MEMDIAG st=%d malloc_zone in-use %.1f MB blocks %u\n", st, ms.size_in_use / 1048576.0, ms.blocks_in_use); }
  fprintf(stderr, "MEMDIAG st=%d %s resident heap>=%zuB %.1f MB, mmap %.1f MB, sites %d\n", st, tag, MD_MIN, tot / 1048576.0, totv / 1048576.0, n);
  { double r = 0, w = 0; for (int i = 0; i < n; i++) if (snap[idx[i]].res < 256 * 1024) { r += (double)snap[idx[i]].res; w += (double)snap[idx[i]].swp; }
    fprintf(stderr, "MEMDIAG st=%d other sites under 0.25 MB resident: res %.1f MB swp %.1f MB\n", st, r / 1048576.0, w / 1048576.0); }
  { int ci[NSITE]; int m = 0; for (int i = 0; i < NSITE; i++) if (snap[i].freed - snap[i].freed_last > 0) ci[m++] = i;
    for (int a = 0; a < m; a++) for (int b = a + 1; b < m; b++) if (snap[ci[b]].freed - snap[ci[b]].freed_last > snap[ci[a]].freed - snap[ci[a]].freed_last) std::swap(ci[a], ci[b]);
    double tf = 0; for (int i = 0; i < m; i++) tf += (double)(snap[ci[i]].freed - snap[ci[i]].freed_last);
    fprintf(stderr, "MEMDIAG st=%d churn freed >=128KB blocks this step: %.1f MB\n", st, tf / 1048576.0);
    for (int i = 0; i < m && i < 25; i++) { const Site& S = snap[ci[i]]; if (S.freed - S.freed_last < 1048576) break;
      fprintf(stderr, "MEMDIAG st=%d churn %.1f MB pc 0x%lx 0x%lx 0x%lx\n", st, (double)(S.freed - S.freed_last) / 1048576.0, (unsigned long)S.pc[0], (unsigned long)S.pc[1], (unsigned long)S.pc[2]); }
    lock(); for (int i = 0; i < NSITE; i++) sites[i].freed_last = sites[i].freed; unlock(); }
  for (int i = 0; i < n && i < 400; i++) { const Site& S = snap[idx[i]]; if (S.res < 256 * 1024) break;
    fprintf(stderr, "MEMDIAG st=%d site %s live %.1f MB res %.1f MB swp %.1f MB n %lld peak %.1f MB pc 0x%lx 0x%lx 0x%lx\n", st, S.kind == 2 ? "mmap" : S.kind ? "vm" : "heap", (double)S.live / 1048576.0, (double)S.res / 1048576.0, (double)S.swp / 1048576.0, (long long)S.nlive, (double)S.peak / 1048576.0,
            (unsigned long)(S.pc[0] ), (unsigned long)(S.pc[1] ), (unsigned long)(S.pc[2] )); }
  busy_set(0);
}
}

namespace memdiag { static double region_res_mb(const void* p, size_t len) {
  const size_t pg = 16384; const uintptr_t a = (uintptr_t)p & ~(uintptr_t)(pg - 1); const size_t l = (((uintptr_t)p + len + pg - 1) & ~(uintptr_t)(pg - 1)) - a; const size_t np = l / pg;
  std::vector<char> v(np); if (mincore((void*)a, l, v.data()) != 0) return -1; size_t c = 0; for (char x : v) c += (x & 1); return (double)(c * pg) / 1048576.0; } }

#include <sys/mman.h>
namespace memdiag {
static inline void* mmap_w(void* a, size_t len, int prot, int fl, int fd, off_t off) {
  void* r = ::mmap(a, len, prot, fl, fd, off);
  if (tab && r != MAP_FAILED && (fl & MAP_ANON) && len >= (1u << 20) && !busy_get()) { busy_set(1); put((uintptr_t)r, len, 2); busy_set(0); }
  return r; }
static inline int munmap_w(void* a, size_t len) { if (tab && !busy_get()) { busy_set(1); take((uintptr_t)a); busy_set(0); } return ::munmap(a, len); }
}

#include <thread>
#include <chrono>
#include <mach/task_info.h>
namespace memdiag {
static std::atomic<double> pk_phys{0}; static std::atomic<double> pk_dev{0}; static std::atomic<long> pk_ms{0};
static Site* pk_snap = nullptr; static std::atomic<bool> samp_run{false}; static std::atomic<int> samp_reset{0};
static double (*dev_alloc_fn)() = nullptr;
static void sampler_main() {
  const auto t0 = std::chrono::steady_clock::now();
  busy_set(1);
  while (samp_run.load()) {
    task_vm_info_data_t vi; mach_msg_type_number_t cn = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vi, &cn) == KERN_SUCCESS) {
      const double ph = (double)vi.phys_footprint / 1048576.0;
      if (samp_reset.exchange(0)) { pk_phys.store(0); }
      if (ph > pk_phys.load()) { lock(); memcpy(pk_snap, sites, sizeof(Site) * NSITE); unlock(); pk_phys.store(ph); pk_dev.store(dev_alloc_fn ? dev_alloc_fn() : 0.0);
        pk_ms.store((long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count()); }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(4));
  }
}
static void start_sampler(double (*devfn)()) { if (samp_run.exchange(true)) return; dev_alloc_fn = devfn; pk_snap = (Site*)calloc(NSITE, sizeof(Site)); std::thread(sampler_main).detach(); }
static void peak_report(int st) {
  busy_set(1);
  fprintf(stderr, "MEMDIAG st=%d peak-in-step phys %.0f MB (Metal allocated %.0f MB) at t=%ld ms\n", st, pk_phys.load(), pk_dev.load(), pk_ms.load());
  static Site endsnap[NSITE]; lock(); memcpy(endsnap, sites, sizeof(endsnap)); unlock();
  int idx[NSITE]; int m = 0; for (int i = 0; i < NSITE; i++) if (pk_snap[i].live - endsnap[i].live > 1048576) idx[m++] = i;
  for (int a = 0; a < m; a++) for (int b = a + 1; b < m; b++) if (pk_snap[idx[b]].live - endsnap[idx[b]].live > pk_snap[idx[a]].live - endsnap[idx[a]].live) std::swap(idx[a], idx[b]);
  double tot = 0; for (int i = 0; i < m; i++) tot += (double)(pk_snap[idx[i]].live - endsnap[idx[i]].live);
  fprintf(stderr, "MEMDIAG st=%d transient live bytes at the peak above the step end (sum of sites over 1 MB): %.1f MB\n", st, tot / 1048576.0);
  for (int i = 0; i < m && i < 25; i++) { const Site& S = pk_snap[idx[i]];
    fprintf(stderr, "MEMDIAG st=%d trans %.1f MB kind %u pc 0x%lx 0x%lx 0x%lx\n", st, (double)(S.live - endsnap[idx[i]].live) / 1048576.0, S.kind, (unsigned long)S.pc[0], (unsigned long)S.pc[1], (unsigned long)S.pc[2]); }
  samp_reset.store(1);
  busy_set(0);
}
}
