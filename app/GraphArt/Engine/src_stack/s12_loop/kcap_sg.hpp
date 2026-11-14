#pragma once
#import <Metal/Metal.h>
#include "kcap_fmt.hpp"
#include "kcap_h.hpp"
#include <atomic>
#include <cstdio>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>
inline std::string g_kcap_kdir;
static inline bool kcap_sg_on() { static const bool on = kcap_dir() && !(getenv("FJN_KCAP_SG") && atoi(getenv("FJN_KCAP_SG")) == 0); return on; }
static inline bool kcap_sg_active() { return kcap_sg_on() && g_kcap_step == kcap_want_step(); }
struct KcapSgPso { std::string name, lib; uint32_t maxThr, execW; };
struct KcapSg {
  std::mutex m;
  std::map<const void*, KcapSgPso> reg;
  std::map<const void*, uint32_t> psoIdx; std::vector<KcapSgPso> psoV;
  std::map<const void*, uint32_t> privIdx; std::vector<uint64_t> privLen;
  std::vector<std::string> libs;
  FILE* f = nullptr; bool finished = false, atexitSet = false, warnedUnk = false; uint32_t nrec = 0, ord = 0; uint64_t bytesW = 0; uint32_t cnt[2] = {0, 0}; uint64_t privMax = 0;
  std::string recPath;
};
static inline KcapSg& kcap_sg() { static KcapSg* s = new KcapSg(); return *s; }
static inline void kcap_sg_reg(id<MTLComputePipelineState> ps, const std::string& name, id<MTLLibrary> L) {
  if (!ps || !kcap_sg_on()) return;
  std::string lib = L && L.label ? std::string(L.label.UTF8String) : std::string();
  if (!g_kcap_kdir.empty() && lib.compare(0, g_kcap_kdir.size() + 1, g_kcap_kdir + "/") == 0) lib = lib.substr(g_kcap_kdir.size() + 1);
  KcapSg& S = kcap_sg(); std::lock_guard<std::mutex> lk(S.m);
  S.reg[(__bridge const void*)ps] = {name, lib, (uint32_t)ps.maxTotalThreadsPerThreadgroup, (uint32_t)ps.threadExecutionWidth};
}
static inline id<MTLLibrary> kcap_newlib(id<MTLDevice> d, NSURL* u, NSError** e) {
  id<MTLLibrary> L = [d newLibraryWithURL:u error:e];
  if (L) L.label = u.path;
  return L;
}
static void kcap_sg_finish();
struct KcapRecState {
  uint32_t kind, dispType; std::vector<uint8_t> cmds; uint32_t ncmd = 0;
  std::map<const void*, uint32_t> shId; std::vector<id<MTLBuffer>> shBufs;
  std::vector<std::pair<uint32_t, uint64_t>> privs;
  std::set<uint32_t> privSeen;
  void put(const SgCmd& c) { const uint8_t* p = (const uint8_t*)&c; cmds.insert(cmds.end(), p, p + sizeof c); ncmd++; }
};
@interface KcapEnc : NSObject
- (instancetype)initWith:(id)real kind:(uint32_t)kind dispType:(uint32_t)dt;
@end
@implementation KcapEnc {
  id _real; KcapRecState* _st;
}
- (instancetype)initWith:(id)real kind:(uint32_t)kind dispType:(uint32_t)dt { self = [super init]; _real = real; _st = new KcapRecState(); _st->kind = kind; _st->dispType = dt; return self; }
- (void)dealloc { delete _st; }
- (id)forwardingTargetForSelector:(SEL)s {
  static std::mutex mm; static std::set<std::string> seen; const std::string n = sel_getName(s);
  static const std::set<std::string> ok = {"pushDebugGroup:", "popDebugGroup", "setLabel:", "label", "insertDebugSignpost:", "device", "retain", "release", "autorelease", "isKindOfClass:", "conformsToProtocol:"};
  if (!ok.count(n)) { std::lock_guard<std::mutex> lk(mm); if (seen.insert(n).second) fprintf(stderr, "KCAP_SG WARNING: encoder selector %s is forwarded but NOT recorded\n", n.c_str()); }
  return _real;
}
- (void)setComputePipelineState:(id<MTLComputePipelineState>)ps {
  KcapSg& S = kcap_sg(); uint32_t idx;
  { std::lock_guard<std::mutex> lk(S.m); const void* k = (__bridge const void*)ps; auto it = S.psoIdx.find(k);
    if (it != S.psoIdx.end()) idx = it->second;
    else { auto r = S.reg.find(k); KcapSgPso p; if (r != S.reg.end()) p = r->second; else { p = {"?", "", (uint32_t)ps.maxTotalThreadsPerThreadgroup, (uint32_t)ps.threadExecutionWidth}; if (!S.warnedUnk) { S.warnedUnk = true; fprintf(stderr, "KCAP_SG WARNING: pipeline not registered by Engine::pipe\n"); } }
      idx = (uint32_t)S.psoV.size(); S.psoV.push_back(p); S.psoIdx[k] = idx; } }
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_PS; c.u0 = idx; _st->put(c);
  [_real setComputePipelineState:ps];
}
- (void)setBuffer:(id<MTLBuffer>)b offset:(NSUInteger)o atIndex:(NSUInteger)i {
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_BUF; c.idx = (uint32_t)i; c.u1 = o;
  if (!b) c.u0 = 0xFFFFFFFFull;
  else if (b.storageMode == MTLStorageModePrivate) {
    KcapSg& S = kcap_sg(); std::lock_guard<std::mutex> lk(S.m); const void* k = (__bridge const void*)b; auto it = S.privIdx.find(k); uint32_t gi;
    if (it == S.privIdx.end()) { gi = (uint32_t)S.privLen.size(); S.privIdx[k] = gi; S.privLen.push_back(b.length); } else gi = it->second;
    S.privLen[gi] = std::max<uint64_t>(S.privLen[gi], b.length);
    if (_st->privSeen.insert(gi).second) _st->privs.push_back({gi, b.length});
    c.u0 = 0x80000000ull | gi;
  } else {
    const void* k = (__bridge const void*)b; auto it = _st->shId.find(k); uint32_t id;
    if (it == _st->shId.end()) { id = (uint32_t)_st->shBufs.size(); _st->shId[k] = id; _st->shBufs.push_back(b); } else id = it->second;
    c.u0 = id;
  }
  _st->put(c);
  [_real setBuffer:b offset:o atIndex:i];
}
- (void)setBytes:(const void*)p length:(NSUInteger)n atIndex:(NSUInteger)i {
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_BYTES; c.idx = (uint32_t)i; c.u0 = n; _st->put(c);
  const size_t pad = (n + 7) & ~(size_t)7; const size_t at = _st->cmds.size(); _st->cmds.resize(at + pad, 0); memcpy(_st->cmds.data() + at, p, n);
  [_real setBytes:p length:n atIndex:i];
}
- (void)setThreadgroupMemoryLength:(NSUInteger)n atIndex:(NSUInteger)i {
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_TG; c.idx = (uint32_t)i; c.u0 = n; _st->put(c);
  [_real setThreadgroupMemoryLength:n atIndex:i];
}
- (void)dispatchThreadgroups:(MTLSize)g threadsPerThreadgroup:(MTLSize)t {
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_DISP; c.g[0] = (uint32_t)g.width; c.g[1] = (uint32_t)g.height; c.g[2] = (uint32_t)g.depth; c.t[0] = (uint32_t)t.width; c.t[1] = (uint32_t)t.height; c.t[2] = (uint32_t)t.depth; _st->put(c);
  [_real dispatchThreadgroups:g threadsPerThreadgroup:t];
}
- (void)memoryBarrierWithScope:(MTLBarrierScope)s {
  SgCmd c; memset(&c, 0, sizeof c); c.op = SG_BAR; c.u0 = (uint64_t)s; _st->put(c);
  [_real memoryBarrierWithScope:s];
}
- (void)endEncoding {
  [_real endEncoding];
  KcapSg& S = kcap_sg(); std::lock_guard<std::mutex> lk(S.m);
  if (S.finished) return;
  if (!S.f) { S.recPath = std::string(kcap_dir()) + "/kcap_sg.recs"; S.f = fopen(S.recPath.c_str(), "wb"); if (!S.f) { fprintf(stderr, "KCAP_SG: cannot open %s\n", S.recPath.c_str()); S.finished = true; return; }
    if (!S.atexitSet) { S.atexitSet = true; atexit([] { kcap_sg_finish(); }); } }
  SgRec r; memset(&r, 0, sizeof r); r.magic = SG_REC_MAGIC; r.kind = _st->kind; r.step = g_kcap_step; r.ord = S.ord++; r.nbuf = (uint32_t)(_st->shBufs.size() + _st->privs.size()); r.ncmd = _st->ncmd; r.dispType = _st->dispType; r.cmdbytes = _st->cmds.size();
  fwrite(&r, sizeof r, 1, S.f); S.bytesW += sizeof r;
  for (size_t i = 0; i < _st->shBufs.size(); i++) { id<MTLBuffer> b = _st->shBufs[i]; SgBufH h; memset(&h, 0, sizeof h); h.id = (uint32_t)i; h.priv = 0; h.len = b.length;
    fwrite(&h, sizeof h, 1, S.f); fwrite(b.contents, 1, b.length, S.f); S.bytesW += sizeof h + b.length; }
  for (auto& pv : _st->privs) { SgBufH h; memset(&h, 0, sizeof h); h.id = 0x80000000u | pv.first; h.priv = 1; h.len = pv.second; fwrite(&h, sizeof h, 1, S.f); S.bytesW += sizeof h; }
  fwrite(_st->cmds.data(), 1, _st->cmds.size(), S.f); S.bytesW += _st->cmds.size();
  fflush(S.f); S.nrec++; S.cnt[_st->kind < 2 ? _st->kind : 0]++;
}
@end
static inline id kcap_sg_wrap(id real, uint32_t kind, uint32_t dispType) { return [[KcapEnc alloc] initWith:real kind:kind dispType:dispType]; }
static void kcap_sg_finish() {
  KcapSg& S = kcap_sg(); std::lock_guard<std::mutex> lk(S.m);
  if (S.finished) return; S.finished = true;
  if (!S.f) { if (kcap_sg_on()) fprintf(stderr, "KCAP_SG: no scan/GN record was written (capture step not reached?)\n"); return; }
  fclose(S.f); S.f = nullptr;
  std::map<std::string, uint32_t> libIdx; std::vector<SgPso> pv;
  for (const KcapSgPso& p : S.psoV) { uint32_t li = 0xFFFFFFFFu;
    if (!p.lib.empty()) { auto it = libIdx.find(p.lib); if (it == libIdx.end()) { li = (uint32_t)S.libs.size(); libIdx[p.lib] = li; S.libs.push_back(p.lib); } else li = it->second; }
    SgPso q; memset(&q, 0, sizeof q); strncpy(q.name, p.name.c_str(), sizeof q.name - 1); q.lib = li; q.maxThr = p.maxThr; q.execW = p.execW; pv.push_back(q); }
  SgHdr h; memset(&h, 0, sizeof h); memcpy(h.magic, "KCSG0001", 8); h.nrec = S.nrec; h.nlib = (uint32_t)S.libs.size(); h.npso = (uint32_t)pv.size();
  h.flags = (getenv("H333") && atoi(getenv("H333"))) ? 1u : 0u;
  { id<MTLDevice> d = MTLCreateSystemDefaultDevice(); if (d) strncpy(h.dev, d.name.UTF8String, sizeof h.dev - 1); }
  const std::string outp = std::string(kcap_dir()) + "/kcap_sg.bin";
  FILE* o = fopen(outp.c_str(), "wb"); FILE* in = fopen(S.recPath.c_str(), "rb");
  if (!o || !in) { fprintf(stderr, "KCAP_SG: cannot finalize %s\n", outp.c_str()); if (o) fclose(o); if (in) fclose(in); return; }
  fwrite(&h, sizeof h, 1, o);
  for (const std::string& l : S.libs) { SgLib x; memset(&x, 0, sizeof x); strncpy(x.path, l.c_str(), sizeof x.path - 1); fwrite(&x, sizeof x, 1, o); }
  fwrite(pv.data(), sizeof(SgPso), pv.size(), o);
  std::vector<char> buf(1 << 20); size_t n; while ((n = fread(buf.data(), 1, buf.size(), in)) > 0) fwrite(buf.data(), 1, n, o);
  fclose(in); fclose(o); remove(S.recPath.c_str());
  uint64_t pm = 0; for (uint64_t v : S.privLen) pm += v;
  fprintf(stderr, "KCAP_SG written %s: records %u (scan %u gn %u) pipelines %zu libs %zu bytes %llu private_buffers %zu private_MB %.1f h333 %u dev %s\n", outp.c_str(), S.nrec, S.cnt[0], S.cnt[1], pv.size(), S.libs.size(),
          (unsigned long long)S.bytesW, S.privLen.size(), pm / 1e6, h.flags & 1, h.dev);
}
static inline void kcap_sg_step(int st) { if (kcap_sg_on() && st > kcap_want_step()) kcap_sg_finish(); }
