#pragma once
#include <cstdint>
struct KcHdr { char magic[8]; uint32_t ncalls; uint32_t flags; };
struct KcCall {
  uint32_t magic; int32_t step; uint32_t npr, ndisp; uint32_t closed, maxit; float flr; uint32_t r48s;
  uint64_t sz[11];
  float yv[96]; float Smv[48 * 9 + 9]; uint32_t rpv[4]; float bwr; uint32_t pad_;
};
struct KcDisp { char name[32]; uint32_t Q, Tn; uint64_t ebo, plo, tgb; uint32_t ntg, tpg, nmem, per; };
static const uint32_t KC_CALL_MAGIC = 0x4C4C4143u;
struct SgHdr { char magic[8]; uint32_t nrec, nlib, npso, flags; char dev[64]; };
struct SgLib { char path[192]; };
struct SgPso { char name[48]; uint32_t lib, maxThr, execW, pad_; };
struct SgRec { uint32_t magic, kind; int32_t step; uint32_t ord, nbuf, ncmd, dispType, pad_; uint64_t cmdbytes; };
struct SgBufH { uint32_t id, priv; uint64_t len; };
enum { SG_PS = 1, SG_BUF = 2, SG_BYTES = 3, SG_TG = 4, SG_DISP = 5, SG_BAR = 6 };
struct SgCmd { uint32_t op, idx; uint64_t u0, u1; uint32_t g[3], t[3]; };
static const uint32_t SG_REC_MAGIC = 0x47534B43u;
