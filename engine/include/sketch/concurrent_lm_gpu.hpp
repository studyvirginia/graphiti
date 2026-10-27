#pragma once

#include "sketch/bytecode_evaluator.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#if defined(__APPLE__) && defined(__OBJC__)
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <mach/mach_time.h>
#endif

namespace Sketch {

constexpr std::size_t kCclmMaxTape = 64;
constexpr std::size_t kCclmMaxCoeffs = 16;
constexpr std::size_t kCclmMaxPoints = 256;

struct CclmInst {
    uint8_t op{0};
    uint8_t ar{0};
    uint16_t s1{0};
    uint16_t s2{0};
    int16_t slot{-1};
    float imm{0.0f};
    float w{1.0f};
};
static_assert(sizeof(CclmInst) == 16, "CclmInst must match the Metal struct layout");

struct CclmCandDesc {
    uint32_t offX{0}, lenX{0}, offY{0}, lenY{0};
    uint32_t nx{0}, n{0};
    uint32_t scrJraw{0}, scrJ{0}, scrR{0};
    uint32_t pad0{0}, pad1{0}, pad2{0};
};
static_assert(sizeof(CclmCandDesc) == 48, "CclmCandDesc must match the Metal struct layout");

struct CclmParams {
    uint32_t N{0}, M{0}, hasNoise{0}, ptClamp{0};
    uint32_t maxIter{8}, maxTrials{6}, screening{1}, pad0{0};
    float ftol{1e-6f}, lambda0{1e-3f}, lambdaUp{10.0f}, lambdaDown{0.1f};
    float lambdaMin{1e-12f}, lambdaMax{1e12f}, floorThr{0.0f}, pad1{0.0f};
};
static_assert(sizeof(CclmParams) == 64, "CclmParams must match the Metal struct layout");

inline auto cclmBuildTape(BytecodeProgram const& prog, int slotBase, std::vector<CclmInst>& out) -> bool
{
    out.clear();
    if (!prog.valid || prog.insts.empty() || prog.insts.size() > kCclmMaxTape) { return false; }
    std::vector<int> slotOf(prog.insts.size(), -1);
    for (std::size_t b = 0; b < prog.bindings.size(); ++b) {
        auto const ii = prog.bindings[b].instIndex;
        if (ii >= prog.insts.size()) { return false; }
        slotOf[ii] = slotBase + static_cast<int>(b);
    }
    std::vector<uint16_t> st;
    st.reserve(prog.insts.size());
    for (std::size_t k = 0; k < prog.insts.size(); ++k) {
        auto const& in = prog.insts[k];
        CclmInst g;
        g.op = static_cast<uint8_t>(in.op);
        g.ar = in.arity;
        g.imm = in.imm;
        g.w = in.weight;
        g.slot = static_cast<int16_t>(slotOf[k]);
        auto const idx = static_cast<uint16_t>(k);
        switch (in.op) {
            case BytecodeOp::Constant:
            case BytecodeOp::Variable:
                st.push_back(idx);
                break;
            case BytecodeOp::Add:
            case BytecodeOp::Mul:
                if (in.arity != 2 || st.size() < 2) { return false; }
                g.s1 = st.back(); st.pop_back(); g.s2 = st.back(); st.back() = idx;
                break;
            case BytecodeOp::Sub:
            case BytecodeOp::Div:
                if (in.arity == 1) {
                    if (st.empty()) { return false; }
                    g.s1 = st.back(); st.back() = idx;
                } else if (in.arity == 2) {
                    if (st.size() < 2) { return false; }
                    g.s1 = st.back(); st.pop_back(); g.s2 = st.back(); st.back() = idx;
                } else {
                    return false;
                }
                break;
            default:
                if (in.arity != 1 || st.empty()) { return false; }
                g.s1 = st.back(); st.back() = idx;
                break;
        }
        if (g.slot >= 0 && in.op != BytecodeOp::Constant && in.op != BytecodeOp::Variable) { return false; }
        out.push_back(g);
    }
    return st.size() == 1 && st[0] == prog.insts.size() - 1;
}

constexpr std::size_t kCclmResStride = 24;

inline auto cclmJRegion(std::size_t N, std::size_t M, std::size_t n) -> std::size_t {
    std::size_t const mpad = (N + M + 7) & ~std::size_t{7};
    return std::max((N + M) * n, mpad * 24);
}

struct CclmJob {
    std::vector<CclmInst> tapeX, tapeY;
    std::vector<float> coeffs;
    uint32_t nx{0};
};

struct CclmFitData {
    std::vector<float> t;
    std::vector<float> dataXY;
    std::vector<float> noise4;
    CclmParams params;
};

struct CclmOut {
    int status{0};
    int iters{0}, lossCalls{0}, builds{0};
    float loss{0.0f};
    std::vector<float> coeffs;
};

struct CclmTiming {
    double encodeMs{0}, swQueueMs{0}, driverMs{0}, hwQueueMs{0}, execMs{0}, notifyMs{0}, readbackMs{0}, totalMs{0};
    double submitHostS{0}, commitHostS{0}, gpuStartS{0}, gpuEndS{0}, waitDoneHostS{0};
};

using CclmCaptureFn = void (*)(CclmFitData const&, std::vector<CclmJob> const&, std::size_t P, void* user);
struct CclmCapture { CclmCaptureFn fn{nullptr}; void* user{nullptr}; };
inline auto cclmCaptureHook() -> CclmCapture& { static CclmCapture c; return c; }

#if defined(__APPLE__) && defined(__OBJC__)

static const char* kCclmShader = R"MSL(
#include <metal_stdlib>
using namespace metal;
#define MAXL 64
#define NMAX 16
#define MAXP 256
#define MAXSG 32

struct GInst { uchar op; uchar ar; ushort s1; ushort s2; short slot; float imm; float w; };
struct GCand { uint offX; uint lenX; uint offY; uint lenY; uint nx; uint n; uint scrJraw; uint scrJ; uint scrR; uint p0; uint p1; uint p2; };
struct GParams { uint N; uint M; uint hasNoise; uint ptClamp; uint maxIter; uint maxTrials; uint screening; uint pad0;
                 float ftol; float lambda0; float lambdaUp; float lambdaDown; float lambdaMin; float lambdaMax; float floorThr; uint pad1u; };

inline float cl_log1p(float x) {
    float u = 1.0f + x;
    if (u == 1.0f) { return x; }
    return log(u) * x / (u - 1.0f);
}
inline float pointLossG(float r2, float s2, float nu, uint clampM) {
    if (!(s2 > 0.0f) || !isfinite(s2)) { s2 = 1.0f; }
    float v;
    if (nu <= 0.0f) { v = r2 / s2; }
    else { v = 0.5f * (nu + 2.0f) * cl_log1p(r2 / (nu * s2)); }
    if (clampM != 0u) { v = max(v, 1.0f); }
    return v;
}
inline float weightOfG(float u, float s2, float nu, uint hasNoise) {
    if (hasNoise == 0u) { return 1.0f; }
    if (!(s2 > 0.0f)) { return 1.0f; }
    if (nu <= 0.0f) { return 1.0f / s2; }
    return 0.5f * (nu + 2.0f) / (nu * s2 + u);
}
inline float segD2(float px, float py, float ax, float ay, float bx, float by) {
    float dx = bx - ax, dy = by - ay;
    float len2 = dx * dx + dy * dy;
    float t = 0.0f;
    if (len2 > 1e-12f) { t = clamp(((px - ax) * dx + (py - ay) * dy) / len2, 0.0f, 1.0f); }
    float qx = ax + t * dx, qy = ay + t * dy;
    float ex = px - qx, ey = py - qy;
    return ex * ex + ey * ey;
}

inline float evalTapeOv(device const GInst* ins, uint L, float t, threadgroup const float* xc, thread float* v,
                        int ovSlot, float ovVal) {
    for (uint k = 0; k < L; ++k) {
        GInst I = ins[k];
        float r = 0.0f;
        switch (I.op) {
            case 0: r = (I.slot >= 0) ? ((I.slot == ovSlot) ? ovVal : xc[I.slot]) : I.imm; break;
            case 1: { float w = (I.slot >= 0) ? ((I.slot == ovSlot) ? ovVal : xc[I.slot]) : I.w; r = t * w; break; }
            case 2: r = (v[I.s2] + v[I.s1]) * I.w; break;
            case 3: r = ((I.ar == 1) ? -v[I.s1] : (v[I.s1] - v[I.s2])) * I.w; break;
            case 4: r = (v[I.s2] * v[I.s1]) * I.w; break;
            case 5: r = ((I.ar == 1) ? (1.0f / v[I.s1]) : (v[I.s1] / v[I.s2])) * I.w; break;
            case 6: r = sin(v[I.s1]) * I.w; break;
            case 7: r = cos(v[I.s1]) * I.w; break;
            case 8: r = exp(v[I.s1]) * I.w; break;
            case 9: r = log(v[I.s1]) * I.w; break;
            case 10: r = sqrt(v[I.s1]) * I.w; break;
            default: { float a = v[I.s1]; r = (a * a) * I.w; break; }
        }
        v[k] = r;
    }
    return v[L - 1];
}
inline float evalTape(device const GInst* ins, uint L, float t, threadgroup const float* xc, thread float* v) {
    return evalTapeOv(ins, L, t, xc, v, -1, 0.0f);
}
template <typename PF>
inline void fdRepair(device const GInst* ins, uint L, float t, threadgroup const float* xc, thread float* v,
                     PF grow) {
    for (uint k = 0; k < L; ++k) {
        int sl = ins[k].slot;
        if (sl < 0 || isfinite(grow[sl])) { continue; }
        float c = xc[sl];
        float h = 1e-4f * max(fabs(c), 1.0f);
        float fp = evalTapeOv(ins, L, t, xc, v, sl, c + h);
        float fm = evalTapeOv(ins, L, t, xc, v, sl, c - h);
        float g = (fp - fm) / (2.0f * h);
        grow[sl] = isfinite(g) ? g : 0.0f;
    }
}

template <typename PF>
inline bool gradTape(device const GInst* ins, uint L, float t, thread const float* v, thread float* adj,
                     PF grow) {
    for (uint k = 0; k < L; ++k) { adj[k] = 0.0f; }
    adj[L - 1] = 1.0f;
    bool ok = true;
    for (int kk = int(L) - 1; kk >= 0; --kk) {
        GInst I = ins[kk];
        float a = adj[kk];
        float w = I.w;
        switch (I.op) {
            case 0: if (I.slot >= 0) { grow[I.slot] = a; ok = ok && isfinite(a); } break;
            case 1: if (I.slot >= 0) { float g = a * t; grow[I.slot] = g; ok = ok && isfinite(g); } break;
            case 2: adj[I.s1] += a * w; adj[I.s2] += a * w; break;
            case 3: if (I.ar == 1) { adj[I.s1] -= a * w; } else { adj[I.s1] += a * w; adj[I.s2] -= a * w; } break;
            case 4: adj[I.s1] += a * w * v[I.s2]; adj[I.s2] += a * w * v[I.s1]; break;
            case 5:
                if (I.ar == 1) { float x = v[I.s1]; adj[I.s1] += -a * w / (x * x); }
                else { float num = v[I.s1], den = v[I.s2]; adj[I.s1] += a * w / den; adj[I.s2] += -a * w * num / (den * den); }
                break;
            case 6: adj[I.s1] += a * w * cos(v[I.s1]); break;
            case 7: adj[I.s1] += -a * w * sin(v[I.s1]); break;
            case 8: adj[I.s1] += a * w * exp(v[I.s1]); break;
            case 9: adj[I.s1] += a * w / v[I.s1]; break;
            case 10: adj[I.s1] += a * w * 0.5f / sqrt(v[I.s1]); break;
            default: adj[I.s1] += a * w * 2.0f * v[I.s1]; break;
        }
    }
    return ok;
}

inline float tgSum(float v, threadgroup float* red, threadgroup float* bc, uint tid, uint lane, uint sg, uint nsg) {
    float s = simd_sum(v);
    if (lane == 0) { red[sg] = s; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (tid == 0) { float tot = 0.0f; for (uint i = 0; i < nsg; ++i) { tot += red[i]; } bc[0] = tot; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float r = bc[0];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return r;
}

inline void argMinComb(thread float& v, thread uint& i, float ov, uint oi) { if (ov < v || (ov == v && oi < i)) { v = ov; i = oi; } }
inline void argMaxComb(thread float& v, thread uint& i, float ov, uint oi) { if (ov > v || (ov == v && oi < i)) { v = ov; i = oi; } }

struct Ctx {
    device const GInst* ix; uint lx;
    device const GInst* iy; uint ly;
    device const float* tv;
    device const float2* data;
    device const float4* noise;
    uint N; uint M; uint hasNoise; uint ptClamp;
    uint tid; uint T; uint lane; uint sg; uint nsg;
};

inline float lossAt(thread const Ctx& c, threadgroup const float* xc, threadgroup float* cxA, threadgroup float* cyA,
                    threadgroup float* red, threadgroup float* bc) {
    float vals[MAXL];
    float xlo = INFINITY, xhi = -INFINITY, ylo = INFINITY, yhi = -INFINITY;
    float bad = 0.0f;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float t = c.tv[j];
        float gx = evalTape(c.ix, c.lx, t, xc, vals);
        float gy = evalTape(c.iy, c.ly, t, xc, vals);
        if (!isfinite(gx) || !isfinite(gy)) { bad = 1.0f; }
        cxA[j] = gx; cyA[j] = gy;
        xlo = min(xlo, gx); xhi = max(xhi, gx); ylo = min(ylo, gy); yhi = max(yhi, gy);
    }
    xlo = simd_min(xlo); xhi = simd_max(xhi); ylo = simd_min(ylo); yhi = simd_max(yhi); bad = simd_max(bad);
    if (c.lane == 0) { red[c.sg * 5 + 0] = xlo; red[c.sg * 5 + 1] = xhi; red[c.sg * 5 + 2] = ylo; red[c.sg * 5 + 3] = yhi; red[c.sg * 5 + 4] = bad; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (c.tid == 0) {
        float a = INFINITY, b = -INFINITY, d = INFINITY, e = -INFINITY, f = 0.0f;
        for (uint s = 0; s < c.nsg; ++s) {
            a = min(a, red[s * 5 + 0]); b = max(b, red[s * 5 + 1]); d = min(d, red[s * 5 + 2]); e = max(e, red[s * 5 + 3]); f = max(f, red[s * 5 + 4]);
        }
        bc[1] = a; bc[2] = b; bc[3] = d; bc[4] = e; bc[5] = f;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    xlo = bc[1]; xhi = bc[2]; ylo = bc[3]; yhi = bc[4]; bad = bc[5];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (bad > 0.0f) { return FLT_MAX; }
    float cxc = (xhi + xlo) / 2.0f, cyc = (yhi + ylo) / 2.0f;
    float sp = max(max(xhi - xlo, yhi - ylo), 1e-9f);
    for (uint j = c.tid; j < c.M; j += c.T) { cxA[j] = (cxA[j] - cxc) / sp; cyA[j] = (cyA[j] - cyc) / sp; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float accD = 0.0f;
    for (uint i = c.tid; i < c.N; i += c.T) {
        float2 pt = c.data[i];
        float best = FLT_MAX;
        for (uint j = 0; j + 1 < c.M; ++j) {
            float d2 = segD2(pt.x, pt.y, cxA[j], cyA[j], cxA[j + 1], cyA[j + 1]);
            if (d2 < best) { best = d2; }
        }
        accD += (c.hasNoise != 0u) ? pointLossG(best, c.noise[i].x, c.noise[i].y, c.ptClamp) : best;
    }
    float accC = 0.0f;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float qx = cxA[j], qy = cyA[j];
        float best = FLT_MAX; uint bi = 0;
        for (uint i = 0; i < c.N; ++i) {
            float2 pt = c.data[i];
            float ddx = qx - pt.x, ddy = qy - pt.y;
            float d2 = ddx * ddx + ddy * ddy;
            if (d2 < best) { best = d2; bi = i; }
        }
        accC += (c.hasNoise != 0u) ? pointLossG(best, c.noise[bi].x, c.noise[bi].y, c.ptClamp) : best;
    }
    accD = tgSum(accD, red, bc, c.tid, c.lane, c.sg, c.nsg) / float(c.N);
    accC = tgSum(accC, red, bc, c.tid, c.lane, c.sg, c.nsg) / float(c.M);
    return 0.5f * (accD + accC);
}

template <typename PF>
inline int buildAt(thread const Ctx& c, uint nx, uint n, threadgroup const float* xc,
                   threadgroup float* cxA, threadgroup float* cyA,
                   threadgroup float* red, threadgroup uint* redI, threadgroup float* bc,
                   threadgroup float* dCx, threadgroup float* dCy, threadgroup float* dS,
                   threadgroup float* H, threadgroup float* gvec,
                   PF Jraw, PF J, PF R) {
    float vals[MAXL];
    float adj[MAXL];
    float bad = 0.0f, badG = 0.0f;
    float umin = INFINITY, umax = -INFINITY, vmin = INFINITY, vmax = -INFINITY;
    uint iumin = 0xFFFFFFFFu, iumax = 0xFFFFFFFFu, ivmin = 0xFFFFFFFFu, ivmax = 0xFFFFFFFFu;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float t = c.tv[j];
        PF row = Jraw + j * n;
        float u = evalTape(c.ix, c.lx, t, xc, vals);
        if (!gradTape(c.ix, c.lx, t, vals, adj, row)) { fdRepair(c.ix, c.lx, t, xc, vals, row); }
        float v = evalTape(c.iy, c.ly, t, xc, vals);
        if (!gradTape(c.iy, c.ly, t, vals, adj, row)) { fdRepair(c.iy, c.ly, t, xc, vals, row); }
        if (!isfinite(u) || !isfinite(v)) { bad = 1.0f; }
        cxA[j] = u; cyA[j] = v;
        argMinComb(umin, iumin, u, j); argMaxComb(umax, iumax, u, j);
        argMinComb(vmin, ivmin, v, j); argMaxComb(vmax, ivmax, v, j);
    }
    for (ushort off = 16; off > 0; off >>= 1) {
        argMinComb(umin, iumin, simd_shuffle_down(umin, off), simd_shuffle_down(iumin, off));
        argMaxComb(umax, iumax, simd_shuffle_down(umax, off), simd_shuffle_down(iumax, off));
        argMinComb(vmin, ivmin, simd_shuffle_down(vmin, off), simd_shuffle_down(ivmin, off));
        argMaxComb(vmax, ivmax, simd_shuffle_down(vmax, off), simd_shuffle_down(ivmax, off));
    }
    bad = simd_max(bad); badG = simd_max(badG);
    if (c.lane == 0) {
        red[c.sg * 6 + 0] = umin; red[c.sg * 6 + 1] = umax; red[c.sg * 6 + 2] = vmin; red[c.sg * 6 + 3] = vmax;
        red[c.sg * 6 + 4] = bad; red[c.sg * 6 + 5] = badG;
        redI[c.sg * 4 + 0] = iumin; redI[c.sg * 4 + 1] = iumax; redI[c.sg * 4 + 2] = ivmin; redI[c.sg * 4 + 3] = ivmax;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (c.tid == 0) {
        float a0 = INFINITY, a1 = -INFINITY, a2 = INFINITY, a3 = -INFINITY, fb = 0.0f, fg = 0.0f;
        uint i0 = 0xFFFFFFFFu, i1 = 0xFFFFFFFFu, i2 = 0xFFFFFFFFu, i3 = 0xFFFFFFFFu;
        for (uint s = 0; s < c.nsg; ++s) {
            argMinComb(a0, i0, red[s * 6 + 0], redI[s * 4 + 0]);
            argMaxComb(a1, i1, red[s * 6 + 1], redI[s * 4 + 1]);
            argMinComb(a2, i2, red[s * 6 + 2], redI[s * 4 + 2]);
            argMaxComb(a3, i3, red[s * 6 + 3], redI[s * 4 + 3]);
            fb = max(fb, red[s * 6 + 4]); fg = max(fg, red[s * 6 + 5]);
        }
        float status = 1.0f;
        if (fb > 0.0f) { status = 0.0f; }
        else if (fg > 0.0f) { status = 2.0f; }
        float Cx = 0.0f, Cy = 0.0f, invS = 0.0f;
        if (status == 1.0f) {
            Cx = 0.5f * (a1 + a0);
            Cy = 0.5f * (a3 + a2);
            float spanU = a1 - a0, spanV = a3 - a2;
            float S = max(max(spanU, spanV), 1e-9f);
            if (!(S > 0.0f) || !isfinite(S)) { status = 0.0f; }
            else {
                invS = 1.0f / S;
                for (uint k = 0; k < n; ++k) {
                    dCx[k] = (k < nx) ? 0.5f * (Jraw[i1 * n + k] + Jraw[i0 * n + k]) : 0.0f;
                    dCy[k] = (k >= nx) ? 0.5f * (Jraw[i3 * n + k] + Jraw[i2 * n + k]) : 0.0f;
                    dS[k] = 0.0f;
                }
                if (spanU >= spanV && spanU > 1e-9f) {
                    for (uint k = 0; k < nx; ++k) { dS[k] = Jraw[i1 * n + k] - Jraw[i0 * n + k]; }
                } else if (spanV > 1e-9f) {
                    for (uint k = nx; k < n; ++k) { dS[k] = Jraw[i3 * n + k] - Jraw[i2 * n + k]; }
                }
            }
        }
        bc[6] = status; bc[7] = Cx; bc[8] = Cy; bc[9] = invS;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float status = bc[6];
    float Cx = bc[7], Cy = bc[8], invS = bc[9];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (status != 1.0f) { return int(status); }
    for (uint j = c.tid; j < c.M; j += c.T) { cxA[j] = (cxA[j] - Cx) * invS; cyA[j] = (cyA[j] - Cy) * invS; }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float badR = 0.0f;
    float const fN = float(c.N), fM = float(c.M);
    for (uint i = c.tid; i < c.N; i += c.T) {
        float2 pt = c.data[i];
        float px = pt.x, py = pt.y;
        uint jb = 0; float tau = 0.0f, best = FLT_MAX;
        for (uint j = 0; j + 1 < c.M; ++j) {
            float ax = cxA[j], ay = cyA[j], bx = cxA[j + 1], by = cyA[j + 1];
            float ex = bx - ax, ey = by - ay;
            float len2 = ex * ex + ey * ey;
            float tt = 0.0f;
            if (len2 > 1e-12f) { tt = clamp(((px - ax) * ex + (py - ay) * ey) / len2, 0.0f, 1.0f); }
            float zx = ax + tt * ex, zy = ay + tt * ey;
            float d2 = (px - zx) * (px - zx) + (py - zy) * (py - zy);
            if (d2 < best) { best = d2; jb = j; tau = tt; }
        }
        float w = weightOfG(best, c.noise[i].z, c.noise[i].w, c.hasNoise);
        float scale = sqrt(max(w, 0.0f) / (2.0f * fN));
        float d = sqrt(best);
        float ri = scale * d;
        R[i] = ri;
        if (!isfinite(ri)) { badR = 1.0f; }
        PF row = J + i * n;
        if (d <= 1e-12f) {
            for (uint k = 0; k < n; ++k) { row[k] = 0.0f; }
            continue;
        }
        float ax = cxA[jb], ay = cyA[jb], bx = cxA[jb + 1], by = cyA[jb + 1];
        float zx = ax + tau * (bx - ax), zy = ay + tau * (by - ay);
        float nhx = (px - zx) / d, nhy = (py - zy) / d;
        PF r0 = Jraw + jb * n;
        PF r1 = Jraw + (jb + 1) * n;
        for (uint k = 0; k < n; ++k) {
            float gu0 = (k < nx) ? r0[k] : 0.0f, gv0 = (k >= nx) ? r0[k] : 0.0f;
            float gu1 = (k < nx) ? r1[k] : 0.0f, gv1 = (k >= nx) ? r1[k] : 0.0f;
            float dqx0 = (gu0 - dCx[k]) * invS - ax * dS[k] * invS, dqy0 = (gv0 - dCy[k]) * invS - ay * dS[k] * invS;
            float dqx1 = (gu1 - dCx[k]) * invS - bx * dS[k] * invS, dqy1 = (gv1 - dCy[k]) * invS - by * dS[k] * invS;
            float g0 = nhx * dqx0 + nhy * dqy0;
            float g1 = nhx * dqx1 + nhy * dqy1;
            float val = -scale * ((1.0f - tau) * g0 + tau * g1);
            row[k] = val;
            if (!isfinite(val)) { badR = 1.0f; }
        }
    }
    for (uint j = c.tid; j < c.M; j += c.T) {
        float qx = cxA[j], qy = cyA[j];
        uint kb = 0; float best = FLT_MAX;
        for (uint i = 0; i < c.N; ++i) {
            float2 pt = c.data[i];
            float ddx = qx - pt.x, ddy = qy - pt.y;
            float d2 = ddx * ddx + ddy * ddy;
            if (d2 < best) { best = d2; kb = i; }
        }
        float w = weightOfG(best, c.noise[kb].z, c.noise[kb].w, c.hasNoise);
        float scale = sqrt(max(w, 0.0f) / (2.0f * fM));
        float e = sqrt(best);
        float rj = scale * e;
        R[c.N + j] = rj;
        if (!isfinite(rj)) { badR = 1.0f; }
        PF row = J + (c.N + j) * n;
        if (e <= 1e-12f) {
            for (uint k = 0; k < n; ++k) { row[k] = 0.0f; }
            continue;
        }
        float2 pk = c.data[kb];
        float mhx = (qx - pk.x) / e, mhy = (qy - pk.y) / e;
        PF rr = Jraw + j * n;
        for (uint k = 0; k < n; ++k) {
            float gu = (k < nx) ? rr[k] : 0.0f, gv = (k >= nx) ? rr[k] : 0.0f;
            float dqx = (gu - dCx[k]) * invS - qx * dS[k] * invS;
            float dqy = (gv - dCy[k]) * invS - qy * dS[k] * invS;
            float val = scale * (mhx * dqx + mhy * dqy);
            row[k] = val;
            if (!isfinite(val)) { badR = 1.0f; }
        }
    }
    badR = simd_max(badR);
    if (c.lane == 0) { red[c.sg] = badR; }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (c.tid == 0) { float fb = 0.0f; for (uint s = 0; s < c.nsg; ++s) { fb = max(fb, red[s]); } bc[10] = fb; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float anyBad = bc[10];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (anyBad > 0.0f) { return 0; }

    uint const m = c.N + c.M;
    uint const npairs = n * (n + 1) / 2;
    for (uint q = c.tid; q < npairs + n; q += c.T) {
        if (q < npairs) {
            uint a = 0, rem = q;
            while (rem >= n - a) { rem -= n - a; ++a; }
            uint b = a + rem;
            float s = 0.0f;
            for (uint i = 0; i < m; ++i) { s += J[i * n + a] * J[i * n + b]; }
            H[a * n + b] = s; H[b * n + a] = s;
        } else {
            uint a = q - npairs;
            float s = 0.0f;
            for (uint i = 0; i < m; ++i) { s += J[i * n + a] * R[i]; }
            gvec[a] = s;
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return 1;
}

inline bool cholSolve(threadgroup const float* H, threadgroup const float* g, threadgroup const float* diag, float lambda, uint n,
                      threadgroup float* L, threadgroup float* y, threadgroup float* dx) {
    for (uint i = 0; i < n * n; ++i) { L[i] = 0.0f; }
    for (uint i = 0; i < n; ++i) {
        for (uint j = 0; j <= i; ++j) {
            float s = H[i * n + j];
            if (i == j) { s += lambda * diag[i]; }
            for (uint k = 0; k < j; ++k) { s -= L[i * n + k] * L[j * n + k]; }
            if (i == j) {
                if (!(s > 0.0f) || !isfinite(s)) { return false; }
                L[i * n + i] = sqrt(s);
            } else {
                L[i * n + j] = s / L[j * n + j];
            }
        }
    }
    for (uint i = 0; i < n; ++i) {
        float s = -g[i];
        for (uint k = 0; k < i; ++k) { s -= L[i * n + k] * y[k]; }
        y[i] = s / L[i * n + i];
    }
    for (int ii = int(n) - 1; ii >= 0; --ii) {
        float s = y[ii];
        for (uint k = uint(ii) + 1; k < n; ++k) { s -= L[k * n + uint(ii)] * dx[k]; }
        dx[ii] = s / L[uint(ii) * n + uint(ii)];
        if (!isfinite(dx[ii])) { return false; }
    }
    return true;
}

kernel void kCclm(device const GInst* insts [[buffer(0)]],
                  device const GCand* cands [[buffer(1)]],
                  device float* coeffs [[buffer(2)]],
                  device const float* tvals [[buffer(3)]],
                  device const float2* data [[buffer(4)]],
                  device const float4* noise [[buffer(5)]],
                  constant GParams& P [[buffer(6)]],
                  device float* scratch [[buffer(7)]],
                  device int* status [[buffer(8)]],
                  device float* lossOut [[buffer(9)]],
                  uint p [[threadgroup_position_in_grid]],
                  uint tid [[thread_position_in_threadgroup]],
                  uint T [[threads_per_threadgroup]],
                  uint lane [[thread_index_in_simdgroup]],
                  uint sg [[simdgroup_index_in_threadgroup]],
                  uint nsg [[simdgroups_per_threadgroup]])
{
    threadgroup float xs[NMAX], xn[NMAX], dxv[NMAX], gvec[NMAX], diag[NMAX], dCx[NMAX], dCy[NMAX], dS[NMAX], yv[NMAX];
    threadgroup float H[NMAX * NMAX];
    threadgroup float Lc[NMAX * NMAX];
    threadgroup float cxA[MAXP], cyA[MAXP];
    threadgroup float red[MAXSG * 6];
    threadgroup uint redI[MAXSG * 4];
    threadgroup float bc[16];
    threadgroup int go[1];

    GCand C = cands[p];
    uint const n = C.n, nx = C.nx;
    Ctx c;
    c.ix = insts + C.offX; c.lx = C.lenX; c.iy = insts + C.offY; c.ly = C.lenY;
    c.tv = tvals; c.data = data; c.noise = noise;
    c.N = P.N; c.M = P.M; c.hasNoise = P.hasNoise; c.ptClamp = P.ptClamp;
    c.tid = tid; c.T = T; c.lane = lane; c.sg = sg; c.nsg = nsg;
    device float* Jraw = scratch + C.scrJraw;
    device float* J = scratch + C.scrJ;
    device float* R = scratch + C.scrR;
    if (tid < n) { xs[tid] = coeffs[p * NMAX + tid]; }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int iters = 0, lossCalls = 0, builds = 0, st = 0;
    float f = lossAt(c, xs, cxA, cyA, red, bc);
    ++lossCalls;
    bool done = !isfinite(f) || (P.floorThr > 0.0f && f <= P.floorThr);
    float lambda = P.lambda0;
    for (uint it = 0; it < P.maxIter && !done; ++it) {
        ++builds;
        int b = buildAt(c, nx, n, xs, cxA, cyA, red, redI, bc, dCx, dCy, dS, H, gvec, Jraw, J, R);
        if (b == 2) { st = 1; break; }
        if (b == 0) { break; }
        if (tid < n) { float d = H[tid * n + tid]; diag[tid] = (d > 0.0f) ? d : 1.0f; }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        bool accepted = false;
        float rel = 0.0f;
        for (uint trial = 0; trial < P.maxTrials; ++trial) {
            if (tid == 0) {
                int ok = cholSolve(H, gvec, diag, lambda, n, Lc, yv, dxv) ? 1 : 0;
                if (ok != 0) {
                    for (uint a = 0; a < n; ++a) {
                        float v = xs[a] + dxv[a];
                        if (!isfinite(v)) { ok = 0; break; }
                        xn[a] = v;
                    }
                }
                go[0] = ok;
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            int g = go[0];
            threadgroup_barrier(mem_flags::mem_threadgroup);
            if (g == 0) { lambda = min(lambda * P.lambdaUp, P.lambdaMax); continue; }
            float fn = lossAt(c, xn, cxA, cyA, red, bc);
            ++lossCalls;
            if (isfinite(fn) && fn < f) {
                rel = (f - fn) / max(fabs(f), 1e-30f);
                if (tid < n) { xs[tid] = xn[tid]; }
                threadgroup_barrier(mem_flags::mem_threadgroup);
                f = fn;
                lambda = max(lambda * P.lambdaDown, P.lambdaMin);
                accepted = true;
                break;
            }
            lambda = min(lambda * P.lambdaUp, P.lambdaMax);
            if (P.screening != 0u && it == 0 && trial == 0 && f > 0.05f) { break; }
        }
        ++iters;
        if (!accepted) { break; }
        if (rel < P.ftol) { break; }
    }
    if (tid < n) { coeffs[p * NMAX + tid] = xs[tid]; }
    if (tid == 0) {
        status[p * 4 + 0] = st; status[p * 4 + 1] = iters; status[p * 4 + 2] = lossCalls; status[p * 4 + 3] = builds;
        lossOut[p] = f;
    }
}

constant bool kBB [[function_constant(0)]];
constant int kMicro [[function_constant(1)]];
constant bool kMicroOn = is_function_constant_defined(kMicro) && kMicro != 0;
#define BS 16
#define NCP 24
#define MAXB (MAXP / BS)

struct SC {
    device const GInst* ix; uint lx;
    device const GInst* iy; uint ly;
    device const float* tv;
    device const float2* data;
    device const float4* noise;
    device const float4* dbox; uint nbD;
    uint N; uint M; uint hasNoise; uint ptClamp;
    uint lane;
    threadgroup float2* cv;
    threadgroup float4* cbox;
};

inline float boxD2(float px, float py, float4 b) {
    float dx = max(max(b.x - px, px - b.z), 0.0f);
    float dy = max(max(b.y - py, py - b.w), 0.0f);
    return dx * dx + dy * dy;
}
inline bool boxSkip(float lb2, float best, float px, float py) {
    float m = sqrt(lb2) - 2e-6f * (1.0f + fabs(px) + fabs(py));
    return m > 0.0f && m * m > best;
}
inline void segScan(float px, float py, threadgroup const float2* cv, uint lo, uint hi,
                    thread float& best, thread uint& jb, thread float& tau) {
    for (uint j = lo; j < hi; ++j) {
        float2 a = cv[j], b = cv[j + 1];
        float ex = b.x - a.x, ey = b.y - a.y;
        float len2 = ex * ex + ey * ey;
        float tt = 0.0f;
        if (len2 > 1e-12f) { tt = clamp(((px - a.x) * ex + (py - a.y) * ey) / len2, 0.0f, 1.0f); }
        float zx = a.x + tt * ex, zy = a.y + tt * ey;
        float d2 = (px - zx) * (px - zx) + (py - zy) * (py - zy);
        if (d2 < best || (d2 == best && j < jb)) { best = d2; jb = j; tau = tt; }
    }
}
inline float nnSeg(thread const SC& c, float px, float py, uint nb, thread int& hint, thread uint& jb, thread float& tau) {
    float best = FLT_MAX; jb = 0xFFFFFFFFu; tau = 0.0f;
    uint const nseg = c.M - 1;
    if (!kBB) { segScan(px, py, c.cv, 0, nseg, best, jb, tau); }
    else {
        int const h = hint;
        if (h >= 0 && uint(h) < nb) { segScan(px, py, c.cv, uint(h) * BS, min(uint(h) * BS + BS, nseg), best, jb, tau); }
        for (uint k = 0; k < nb; ++k) {
            if (int(k) == h) { continue; }
            if (boxSkip(boxD2(px, py, c.cbox[k]), best, px, py)) { continue; }
            segScan(px, py, c.cv, k * BS, min(k * BS + BS, nseg), best, jb, tau);
        }
        hint = int(jb / BS);
    }
    if (jb == 0xFFFFFFFFu) { jb = 0; }
    return best;
}
inline void ptScan(float qx, float qy, device const float2* data, uint lo, uint hi, thread float& best, thread uint& kb) {
    for (uint i = lo; i < hi; ++i) {
        float2 pt = data[i];
        float ddx = qx - pt.x, ddy = qy - pt.y;
        float d2 = ddx * ddx + ddy * ddy;
        if (d2 < best || (d2 == best && i < kb)) { best = d2; kb = i; }
    }
}
inline float nnPt(thread const SC& c, float qx, float qy, thread int& hint, thread uint& kb) {
    float best = FLT_MAX; kb = 0xFFFFFFFFu;
    if (!kBB) { ptScan(qx, qy, c.data, 0, c.N, best, kb); }
    else {
        int const h = hint;
        uint const nb = c.nbD;
        if (h >= 0 && uint(h) < nb) { ptScan(qx, qy, c.data, uint(h) * BS, min(uint(h) * BS + BS, c.N), best, kb); }
        for (uint k = 0; k < nb; ++k) {
            if (int(k) == h) { continue; }
            if (boxSkip(boxD2(qx, qy, c.dbox[k]), best, qx, qy)) { continue; }
            ptScan(qx, qy, c.data, k * BS, min(k * BS + BS, c.N), best, kb);
        }
        hint = int(kb / BS);
    }
    if (kb == 0xFFFFFFFFu) { kb = 0; }
    return best;
}
inline uint segBlocks(thread const SC& c) {
    uint const nseg = c.M - 1;
    uint const nb = (nseg + BS - 1) / BS;
    if (kBB) {
        for (uint k = c.lane; k < nb; k += 32) {
            uint lo = k * BS, hi = min(lo + BS, nseg);
            float2 v = c.cv[lo];
            float4 b = float4(v.x, v.y, v.x, v.y);
            for (uint j = lo + 1; j <= hi; ++j) { float2 w = c.cv[j]; b.x = min(b.x, w.x); b.y = min(b.y, w.y); b.z = max(b.z, w.x); b.w = max(b.w, w.y); }
            c.cbox[k] = b;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    return nb;
}

inline float lossS(thread const SC& c, threadgroup const float* xc) {
    float vals[MAXL];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float xlo = INFINITY, xhi = -INFINITY, ylo = INFINITY, yhi = -INFINITY;
    bool bad = false;
    for (uint j = c.lane; j < c.M; j += 32) {
        float t = c.tv[j];
        float gx = evalTape(c.ix, c.lx, t, xc, vals);
        float gy = evalTape(c.iy, c.ly, t, xc, vals);
        if (!isfinite(gx) || !isfinite(gy)) { bad = true; }
        c.cv[j] = float2(gx, gy);
        xlo = min(xlo, gx); xhi = max(xhi, gx); ylo = min(ylo, gy); yhi = max(yhi, gy);
    }
    if (simd_any(bad)) { return FLT_MAX; }
    xlo = simd_min(xlo); xhi = simd_max(xhi); ylo = simd_min(ylo); yhi = simd_max(yhi);
    float cxc = (xhi + xlo) / 2.0f, cyc = (yhi + ylo) / 2.0f;
    float sp = max(max(xhi - xlo, yhi - ylo), 1e-9f);
    for (uint j = c.lane; j < c.M; j += 32) { float2 v = c.cv[j]; c.cv[j] = float2((v.x - cxc) / sp, (v.y - cyc) / sp); }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (kMicroOn && kMicro == 2) { return sp; }
    uint const nb = segBlocks(c);
    float accD = 0.0f;
    uint const chN = (c.N + 31) / 32;
    uint const i0 = min(c.lane * chN, c.N), i1 = min(i0 + chN, c.N);
    int hint = (i0 < c.N) ? int(((i0 * (c.M - 1)) / c.N) / BS) : -1;
    for (uint i = i0; i < i1; ++i) {
        float2 pt = c.data[i];
        uint jb; float tau;
        float best = nnSeg(c, pt.x, pt.y, nb, hint, jb, tau);
        accD += (c.hasNoise != 0u) ? pointLossG(best, c.noise[i].x, c.noise[i].y, c.ptClamp) : best;
    }
    float accC = 0.0f;
    uint const chM = (c.M + 31) / 32;
    uint const j0 = min(c.lane * chM, c.M), j1 = min(j0 + chM, c.M);
    int hintP = (j0 < c.M) ? int(((j0 * c.N) / c.M) / BS) : -1;
    for (uint j = j0; j < j1; ++j) {
        float2 q = c.cv[j];
        uint kb;
        float best = nnPt(c, q.x, q.y, hintP, kb);
        accC += (c.hasNoise != 0u) ? pointLossG(best, c.noise[kb].x, c.noise[kb].y, c.ptClamp) : best;
    }
    accD = simd_sum(accD) / float(c.N);
    accC = simd_sum(accC) / float(c.M);
    return 0.5f * (accD + accC);
}

#define CCLM_FOLD() do { \
        threadgroup_barrier(mem_flags::mem_threadgroup); \
        for (uint r0 = 0; r0 < 32; r0 += 8) { \
            for (uint a = 0; a < 3; ++a) { \
                if (a >= nt) { break; } \
                simdgroup_float8x8 A; \
                simdgroup_load(A, st, NCP, ulong2(a * 8, r0), true); \
                for (uint b = 0; b < 3; ++b) { \
                    if (b >= nt) { break; } \
                    simdgroup_float8x8 B; \
                    simdgroup_load(B, st, NCP, ulong2(b * 8, r0), false); \
                    simdgroup_multiply_accumulate(acc[a][b], A, B, acc[a][b]); \
                } \
            } \
        } \
        threadgroup_barrier(mem_flags::mem_threadgroup); \
    } while (0)

inline int buildS(thread const SC& c, uint nx, uint n, threadgroup const float* xc,
                  threadgroup float* dCx, threadgroup float* dCy, threadgroup float* dS,
                  threadgroup float* st, device float* Jraw) {
    float vals[MAXL];
    float adj[MAXL];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    bool bad = false;
    float umin = INFINITY, umax = -INFINITY, vmin = INFINITY, vmax = -INFINITY;
    uint iumin = 0xFFFFFFFFu, iumax = 0xFFFFFFFFu, ivmin = 0xFFFFFFFFu, ivmax = 0xFFFFFFFFu;
    for (uint j = c.lane; j < c.M; j += 32) {
        float t = c.tv[j];
        device float* row = Jraw + j * n;
        float u = evalTape(c.ix, c.lx, t, xc, vals);
        if (!gradTape(c.ix, c.lx, t, vals, adj, row)) { fdRepair(c.ix, c.lx, t, xc, vals, row); }
        float v = evalTape(c.iy, c.ly, t, xc, vals);
        if (!gradTape(c.iy, c.ly, t, vals, adj, row)) { fdRepair(c.iy, c.ly, t, xc, vals, row); }
        if (!isfinite(u) || !isfinite(v)) { bad = true; }
        c.cv[j] = float2(u, v);
        argMinComb(umin, iumin, u, j); argMaxComb(umax, iumax, u, j);
        argMinComb(vmin, ivmin, v, j); argMaxComb(vmax, ivmax, v, j);
    }
    if (simd_any(bad)) { return 0; }
    if (kMicroOn && kMicro == 4) { return 1; }
    for (ushort off = 16; off > 0; off >>= 1) {
        argMinComb(umin, iumin, simd_shuffle_xor(umin, off), simd_shuffle_xor(iumin, off));
        argMaxComb(umax, iumax, simd_shuffle_xor(umax, off), simd_shuffle_xor(iumax, off));
        argMinComb(vmin, ivmin, simd_shuffle_xor(vmin, off), simd_shuffle_xor(ivmin, off));
        argMaxComb(vmax, ivmax, simd_shuffle_xor(vmax, off), simd_shuffle_xor(ivmax, off));
    }
    float const Cx = 0.5f * (umax + umin);
    float const Cy = 0.5f * (vmax + vmin);
    float const spanU = umax - umin, spanV = vmax - vmin;
    float const S = max(max(spanU, spanV), 1e-9f);
    if (!(S > 0.0f) || !isfinite(S)) { return 0; }
    float const invS = 1.0f / S;
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (c.lane < n) {
        uint const k = c.lane;
        dCx[k] = (k < nx) ? 0.5f * (Jraw[iumax * n + k] + Jraw[iumin * n + k]) : 0.0f;
        dCy[k] = (k >= nx) ? 0.5f * (Jraw[ivmax * n + k] + Jraw[ivmin * n + k]) : 0.0f;
        float ds = 0.0f;
        if (spanU >= spanV && spanU > 1e-9f) { if (k < nx) { ds = Jraw[iumax * n + k] - Jraw[iumin * n + k]; } }
        else if (spanV > 1e-9f) { if (k >= nx) { ds = Jraw[ivmax * n + k] - Jraw[ivmin * n + k]; } }
        dS[k] = ds;
    }
    for (uint j = c.lane; j < c.M; j += 32) { float2 v = c.cv[j]; c.cv[j] = float2((v.x - Cx) * invS, (v.y - Cy) * invS); }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint const nb = segBlocks(c);

    uint const nt = (n + 1 + 7) / 8;     // 8-wide tiles covering [J | r]
    uint const ncol = nt * 8;
    simdgroup_float8x8 acc[3][3];
    for (uint a = 0; a < 3; ++a) { for (uint b = 0; b < 3; ++b) { acc[a][b] = make_filled_simdgroup_matrix<float, 8, 8>(0.0f); } }
    bool badR = false;
    float const fN = float(c.N), fM = float(c.M);
    threadgroup float* srow = st + c.lane * NCP;
    uint const chN = (c.N + 31) / 32;
    int hint = int(((min(c.lane * chN, c.N - 1) * (c.M - 1)) / c.N) / BS);
    for (uint ch = 0; ch < chN; ++ch) {
        uint const i = c.lane * chN + ch;
        for (uint k = 0; k < ncol; ++k) { srow[k] = 0.0f; }
        if (i < c.N) {
            float2 pt = c.data[i];
            float px = pt.x, py = pt.y;
            uint jb; float tau;
            float best;
            if (kMicroOn && kMicro == 6) { jb = 0; tau = 0.5f; float2 a0 = c.cv[0]; best = (px - a0.x) * (px - a0.x) + (py - a0.y) * (py - a0.y); }
            else { best = nnSeg(c, px, py, nb, hint, jb, tau); }
            float w = weightOfG(best, c.noise[i].z, c.noise[i].w, c.hasNoise);
            float scale = sqrt(max(w, 0.0f) / (2.0f * fN));
            float d = sqrt(best);
            float ri = scale * d;
            srow[n] = ri;
            if (!isfinite(ri)) { badR = true; }
            if (d > 1e-12f) {
                float2 A2 = c.cv[jb], B2 = c.cv[jb + 1];
                float ax = A2.x, ay = A2.y, bx = B2.x, by = B2.y;
                float zx = ax + tau * (bx - ax), zy = ay + tau * (by - ay);
                float nhx = (px - zx) / d, nhy = (py - zy) / d;
                device const float* r0 = Jraw + jb * n;
                device const float* r1 = Jraw + (jb + 1) * n;
                for (uint k = 0; k < n; ++k) {
                    float gu0 = (k < nx) ? r0[k] : 0.0f, gv0 = (k >= nx) ? r0[k] : 0.0f;
                    float gu1 = (k < nx) ? r1[k] : 0.0f, gv1 = (k >= nx) ? r1[k] : 0.0f;
                    float dqx0 = (gu0 - dCx[k]) * invS - ax * dS[k] * invS, dqy0 = (gv0 - dCy[k]) * invS - ay * dS[k] * invS;
                    float dqx1 = (gu1 - dCx[k]) * invS - bx * dS[k] * invS, dqy1 = (gv1 - dCy[k]) * invS - by * dS[k] * invS;
                    float g0 = nhx * dqx0 + nhy * dqy0;
                    float g1 = nhx * dqx1 + nhy * dqy1;
                    float val = -scale * ((1.0f - tau) * g0 + tau * g1);
                    srow[k] = val;
                    if (!isfinite(val)) { badR = true; }
                }
            }
        }
        CCLM_FOLD();
    }
    uint const chM = (c.M + 31) / 32;
    int hintP = int(((min(c.lane * chM, c.M - 1) * c.N) / c.M) / BS);
    for (uint ch = 0; ch < chM; ++ch) {
        uint const j = c.lane * chM + ch;
        for (uint k = 0; k < ncol; ++k) { srow[k] = 0.0f; }
        if (j < c.M) {
            float2 q = c.cv[j];
            float qx = q.x, qy = q.y;
            uint kb;
            float best;
            if (kMicroOn && kMicro == 6) { kb = 0; float2 d0 = c.data[0]; best = (qx - d0.x) * (qx - d0.x) + (qy - d0.y) * (qy - d0.y); }
            else { best = nnPt(c, qx, qy, hintP, kb); }
            float w = weightOfG(best, c.noise[kb].z, c.noise[kb].w, c.hasNoise);
            float scale = sqrt(max(w, 0.0f) / (2.0f * fM));
            float e = sqrt(best);
            float rj = scale * e;
            srow[n] = rj;
            if (!isfinite(rj)) { badR = true; }
            if (e > 1e-12f) {
                float2 pk = c.data[kb];
                float mhx = (qx - pk.x) / e, mhy = (qy - pk.y) / e;
                device const float* rr = Jraw + j * n;
                for (uint k = 0; k < n; ++k) {
                    float gu = (k < nx) ? rr[k] : 0.0f, gv = (k >= nx) ? rr[k] : 0.0f;
                    float dqx = (gu - dCx[k]) * invS - qx * dS[k] * invS;
                    float dqy = (gv - dCy[k]) * invS - qy * dS[k] * invS;
                    float val = scale * (mhx * dqx + mhy * dqy);
                    srow[k] = val;
                    if (!isfinite(val)) { badR = true; }
                }
            }
        }
        CCLM_FOLD();
    }
    if (simd_any(badR)) { return 0; }
    for (uint a = 0; a < 3; ++a) {
        if (a >= nt) { break; }
        for (uint b = 0; b < 3; ++b) {
            if (b >= nt) { break; }
            simdgroup_store(acc[a][b], st, NCP, ulong2(b * 8, a * 8));
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return 1;
}

inline bool cholS(thread const float* h, float gi, float dg, float lambda, uint n, uint lane,
                  threadgroup float* Ls, thread float& dxOut) {
    float a[NMAX];
    for (uint k = 0; k < NMAX; ++k) { a[k] = h[k]; if (k == lane) { a[k] += lambda * dg; } }
    for (uint j = 0; j < NMAX; ++j) {
        if (j >= n) { break; }
        float ajj = simd_shuffle(a[j], ushort(j));
        if (!(ajj > 0.0f) || !isfinite(ajj)) { return false; }
        float ljj = sqrt(ajj);
        float lij = (lane == j) ? ljj : a[j] / ljj;
        a[j] = lij;
        for (uint cc = j + 1; cc < NMAX; ++cc) {
            if (cc >= n) { break; }
            float lcj = simd_shuffle(lij, ushort(cc));
            if (lane >= cc) { a[cc] -= lij * lcj; }
        }
    }
    float val = -gi;
    for (uint k = 0; k < NMAX; ++k) {
        if (k >= n) { break; }
        float yk = simd_shuffle(val / a[k], ushort(k));
        if (lane == k) { val = yk; } else if (lane > k) { val -= a[k] * yk; }
    }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    for (uint k = 0; k < NMAX; ++k) { if (k < n && lane < n && k <= lane) { Ls[lane * NCP + k] = a[k]; } }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    float const lii = (lane < n) ? Ls[lane * NCP + lane] : 1.0f;
    for (int ii = int(n) - 1; ii >= 0; --ii) {
        float dxi = simd_shuffle(val / lii, ushort(ii));
        if (int(lane) == ii) { val = dxi; } else if (int(lane) < ii) { val -= Ls[uint(ii) * NCP + lane] * dxi; }
    }
    dxOut = val;
    bool ok = (lane >= n) || isfinite(val);
    return simd_all(ok);
}

kernel void kCclmS(device const GInst* insts [[buffer(0)]],
                   device const GCand* cands [[buffer(1)]],
                   device float* coeffs [[buffer(2)]],
                   device const float* tvals [[buffer(3)]],
                   device const float2* data [[buffer(4)]],
                   device const float4* noise [[buffer(5)]],
                   constant GParams& P [[buffer(6)]],
                   device float* scratch [[buffer(7)]],
                   device int* status [[buffer(8)]],
                   device float* lossOut [[buffer(9)]],
                   device const float4* dbox [[buffer(10)]],
                   uint p [[threadgroup_position_in_grid]],
                   uint lane [[thread_index_in_simdgroup]])
{
    threadgroup float xs[NMAX], xn[NMAX], dCx[NMAX], dCy[NMAX], dS[NMAX];
    threadgroup float2 cv[MAXP];
    threadgroup float4 cbox[MAXB];
    threadgroup float st[32 * NCP];

    GCand C = cands[p];
    uint const n = C.n, nx = C.nx;
    SC c;
    c.ix = insts + C.offX; c.lx = C.lenX; c.iy = insts + C.offY; c.ly = C.lenY;
    c.tv = tvals; c.data = data; c.noise = noise; c.dbox = dbox; c.nbD = P.pad0;
    c.N = P.N; c.M = P.M; c.hasNoise = P.hasNoise; c.ptClamp = P.ptClamp;
    c.lane = lane; c.cv = cv; c.cbox = cbox;
    device float* Jraw = scratch + C.scrJraw;
    if (lane < n) { xs[lane] = coeffs[p * NMAX + lane]; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (kMicroOn) {
        float accm = 0.0f;
        uint const R = P.maxIter;       // repetitions
        for (uint r = 0; r < R; ++r) {
            if (kMicro == 1 || kMicro == 2) { accm += lossS(c, xs); }
            else if (kMicro == 3 || kMicro == 4 || kMicro == 6) { accm += float(buildS(c, nx, n, xs, dCx, dCy, dS, st, Jraw)); }
            else if (kMicro == 5) {
                float hm[NMAX];
                for (uint k = 0; k < NMAX; ++k) { hm[k] = (lane < n && k < n) ? ((k == lane) ? 4.0f : 0.5f / float(1 + k + lane)) : 0.0f; }
                float dx = 0.0f;
                accm += cholS(hm, 1.0f, 1.0f, 1e-3f * float(r + 1), n, lane, st, dx) ? dx : 0.0f;
            }
        }
        if (lane == 0) { lossOut[p] = accm; status[p * 4 + 0] = 0; status[p * 4 + 1] = 0; status[p * 4 + 2] = 0; status[p * 4 + 3] = 0; }
        return;
    }

    int iters = 0, lossCalls = 0, builds = 0;
    float f = lossS(c, xs);
    ++lossCalls;
    bool done = !isfinite(f) || (P.floorThr > 0.0f && f <= P.floorThr);
    float lambda = P.lambda0;
    for (uint it = 0; it < P.maxIter && !done; ++it) {
        ++builds;
        int b = buildS(c, nx, n, xs, dCx, dCy, dS, st, Jraw);
        if (b == 0) { break; }
        float h[NMAX];
        for (uint k = 0; k < NMAX; ++k) { h[k] = (lane < n && k < n) ? st[lane * NCP + k] : 0.0f; }
        float const gi = (lane < n) ? st[lane * NCP + n] : 0.0f;
        float const hii = (lane < n) ? st[lane * NCP + lane] : 1.0f;
        float const dg = (hii > 0.0f) ? hii : 1.0f;
        bool accepted = false;
        float rel = 0.0f;
        for (uint trial = 0; trial < P.maxTrials; ++trial) {
            float dx = 0.0f;
            bool ok = cholS(h, gi, dg, lambda, n, lane, st, dx);
            if (ok) {
                bool fin = true;
                if (lane < n) { float v = xs[lane] + dx; fin = isfinite(v); xn[lane] = v; }
                ok = simd_all(fin);
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            if (!ok) { lambda = min(lambda * P.lambdaUp, P.lambdaMax); continue; }
            float fn = lossS(c, xn);
            ++lossCalls;
            if (isfinite(fn) && fn < f) {
                rel = (f - fn) / max(fabs(f), 1e-30f);
                if (lane < n) { xs[lane] = xn[lane]; }
                threadgroup_barrier(mem_flags::mem_threadgroup);
                f = fn;
                lambda = max(lambda * P.lambdaDown, P.lambdaMin);
                accepted = true;
                break;
            }
            lambda = min(lambda * P.lambdaUp, P.lambdaMax);
            if (P.screening != 0u && it == 0 && trial == 0 && f > 0.05f) { break; }
        }
        ++iters;
        if (!accepted) { break; }
        if (rel < P.ftol) { break; }
    }
    if (lane < n) { coeffs[p * NMAX + lane] = xs[lane]; }
    if (lane == 0) {
        status[p * 4 + 0] = 0; status[p * 4 + 1] = iters; status[p * 4 + 2] = lossCalls; status[p * 4 + 3] = builds;
        lossOut[p] = f;
    }
}
#define MAXG 8
#define RSTRIDE 24
struct GC {
    device const GInst* ix; uint lx;
    device const GInst* iy; uint ly;
    device const float* tv;
    device const float2* data;
    device const float4* noise;
    device const float4* dbox; uint nbD;
    uint N; uint M; uint hasNoise; uint ptClamp;
    uint tid; uint T; uint lane; uint sg; uint G;
    threadgroup float2* cv;
    threadgroup float4* cbox;
    threadgroup float* red;      // MAXG * 16
    threadgroup uint* redI;      // MAXG * 4
};

inline float nnSegG(thread const GC& c, float px, float py, uint nb, thread uint& jb, thread float& tau) {
    float best = FLT_MAX; jb = 0xFFFFFFFFu; tau = 0.0f;
    uint const nseg = c.M - 1;
    if (!kBB) { segScan(px, py, c.cv, 0, nseg, best, jb, tau); }
    else {
        uint k0 = 0; float lb0 = FLT_MAX;
        for (uint k = 0; k < nb; ++k) { float lb = boxD2(px, py, c.cbox[k]); if (lb < lb0) { lb0 = lb; k0 = k; } }
        segScan(px, py, c.cv, k0 * BS, min(k0 * BS + BS, nseg), best, jb, tau);
        for (uint k = 0; k < nb; ++k) {
            if (k == k0) { continue; }
            if (boxSkip(boxD2(px, py, c.cbox[k]), best, px, py)) { continue; }
            segScan(px, py, c.cv, k * BS, min(k * BS + BS, nseg), best, jb, tau);
        }
    }
    if (jb == 0xFFFFFFFFu) { jb = 0; }
    return best;
}
inline float nnPtG(thread const GC& c, float qx, float qy, thread uint& kb) {
    float best = FLT_MAX; kb = 0xFFFFFFFFu;
    if (!kBB) { ptScan(qx, qy, c.data, 0, c.N, best, kb); }
    else {
        uint const nb = c.nbD;
        uint k0 = 0; float lb0 = FLT_MAX;
        for (uint k = 0; k < nb; ++k) { float lb = boxD2(qx, qy, c.dbox[k]); if (lb < lb0) { lb0 = lb; k0 = k; } }
        ptScan(qx, qy, c.data, k0 * BS, min(k0 * BS + BS, c.N), best, kb);
        for (uint k = 0; k < nb; ++k) {
            if (k == k0) { continue; }
            if (boxSkip(boxD2(qx, qy, c.dbox[k]), best, qx, qy)) { continue; }
            ptScan(qx, qy, c.data, k * BS, min(k * BS + BS, c.N), best, kb);
        }
    }
    if (kb == 0xFFFFFFFFu) { kb = 0; }
    return best;
}
inline uint segBlocksG(thread const GC& c) {
    uint const nseg = c.M - 1;
    uint const nb = (nseg + BS - 1) / BS;
    if (kBB) {
        for (uint k = c.tid; k < nb; k += c.T) {
            uint lo = k * BS, hi = min(lo + BS, nseg);
            float2 v = c.cv[lo];
            float4 b = float4(v.x, v.y, v.x, v.y);
            for (uint j = lo + 1; j <= hi; ++j) { float2 w = c.cv[j]; b.x = min(b.x, w.x); b.y = min(b.y, w.y); b.z = max(b.z, w.x); b.w = max(b.w, w.y); }
            c.cbox[k] = b;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    return nb;
}

inline float lossG(thread const GC& c, threadgroup const float* xc) {
    float vals[MAXL];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float xlo = INFINITY, xhi = -INFINITY, ylo = INFINITY, yhi = -INFINITY, bad = 0.0f;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float t = c.tv[j];
        float gx = evalTape(c.ix, c.lx, t, xc, vals);
        float gy = evalTape(c.iy, c.ly, t, xc, vals);
        if (!isfinite(gx) || !isfinite(gy)) { bad = 1.0f; }
        c.cv[j] = float2(gx, gy);
        xlo = min(xlo, gx); xhi = max(xhi, gx); ylo = min(ylo, gy); yhi = max(yhi, gy);
    }
    xlo = simd_min(xlo); xhi = simd_max(xhi); ylo = simd_min(ylo); yhi = simd_max(yhi); bad = simd_max(bad);
    if (c.lane == 0) { c.red[0 * MAXG + c.sg] = xlo; c.red[1 * MAXG + c.sg] = xhi; c.red[2 * MAXG + c.sg] = ylo; c.red[3 * MAXG + c.sg] = yhi; c.red[4 * MAXG + c.sg] = bad; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    xlo = INFINITY; xhi = -INFINITY; ylo = INFINITY; yhi = -INFINITY; bad = 0.0f;
    for (uint g = 0; g < c.G; ++g) {
        xlo = min(xlo, c.red[0 * MAXG + g]); xhi = max(xhi, c.red[1 * MAXG + g]);
        ylo = min(ylo, c.red[2 * MAXG + g]); yhi = max(yhi, c.red[3 * MAXG + g]); bad = max(bad, c.red[4 * MAXG + g]);
    }
    if (bad > 0.0f) { return FLT_MAX; }
    float cxc = (xhi + xlo) / 2.0f, cyc = (yhi + ylo) / 2.0f;
    float sp = max(max(xhi - xlo, yhi - ylo), 1e-9f);
    for (uint j = c.tid; j < c.M; j += c.T) { float2 v = c.cv[j]; c.cv[j] = float2((v.x - cxc) / sp, (v.y - cyc) / sp); }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (kMicroOn && kMicro == 2) { return sp; }
    uint const nb = segBlocksG(c);
    float accD = 0.0f;
    for (uint i = c.tid; i < c.N; i += c.T) {
        float2 pt = c.data[i];
        uint jb; float tau;
        float best = nnSegG(c, pt.x, pt.y, nb, jb, tau);
        accD += (c.hasNoise != 0u) ? pointLossG(best, c.noise[i].x, c.noise[i].y, c.ptClamp) : best;
    }
    float accC = 0.0f;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float2 q = c.cv[j];
        uint kb;
        float best = nnPtG(c, q.x, q.y, kb);
        accC += (c.hasNoise != 0u) ? pointLossG(best, c.noise[kb].x, c.noise[kb].y, c.ptClamp) : best;
    }
    accD = simd_sum(accD); accC = simd_sum(accC);
    if (c.lane == 0) { c.red[5 * MAXG + c.sg] = accD; c.red[6 * MAXG + c.sg] = accC; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float sD = 0.0f, sC = 0.0f;
    for (uint g = 0; g < c.G; ++g) { sD += c.red[5 * MAXG + g]; sC += c.red[6 * MAXG + g]; }
    return 0.5f * (sD / float(c.N) + sC / float(c.M));
}

inline int buildG(thread const GC& c, uint nx, uint n, threadgroup const float* xc,
                  threadgroup float* dCx, threadgroup float* dCy, threadgroup float* dS,
                  threadgroup float* Hs, device float* Jraw, device float* Jd) {
    float vals[MAXL];
    float adj[MAXL];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float bad = 0.0f;
    float umin = INFINITY, umax = -INFINITY, vmin = INFINITY, vmax = -INFINITY;
    uint iumin = 0xFFFFFFFFu, iumax = 0xFFFFFFFFu, ivmin = 0xFFFFFFFFu, ivmax = 0xFFFFFFFFu;
    for (uint j = c.tid; j < c.M; j += c.T) {
        float t = c.tv[j];
        device float* row = Jraw + j * n;
        float u = evalTape(c.ix, c.lx, t, xc, vals);
        if (!gradTape(c.ix, c.lx, t, vals, adj, row)) { fdRepair(c.ix, c.lx, t, xc, vals, row); }
        float v = evalTape(c.iy, c.ly, t, xc, vals);
        if (!gradTape(c.iy, c.ly, t, vals, adj, row)) { fdRepair(c.iy, c.ly, t, xc, vals, row); }
        if (!isfinite(u) || !isfinite(v)) { bad = 1.0f; }
        c.cv[j] = float2(u, v);
        argMinComb(umin, iumin, u, j); argMaxComb(umax, iumax, u, j);
        argMinComb(vmin, ivmin, v, j); argMaxComb(vmax, ivmax, v, j);
    }
    for (ushort off = 16; off > 0; off >>= 1) {
        argMinComb(umin, iumin, simd_shuffle_xor(umin, off), simd_shuffle_xor(iumin, off));
        argMaxComb(umax, iumax, simd_shuffle_xor(umax, off), simd_shuffle_xor(iumax, off));
        argMinComb(vmin, ivmin, simd_shuffle_xor(vmin, off), simd_shuffle_xor(ivmin, off));
        argMaxComb(vmax, ivmax, simd_shuffle_xor(vmax, off), simd_shuffle_xor(ivmax, off));
    }
    bad = simd_max(bad);
    if (c.lane == 0) {
        c.red[8 * MAXG + c.sg] = umin; c.red[9 * MAXG + c.sg] = umax; c.red[10 * MAXG + c.sg] = vmin; c.red[11 * MAXG + c.sg] = vmax;
        c.red[12 * MAXG + c.sg] = bad;
        c.redI[0 * MAXG + c.sg] = iumin; c.redI[1 * MAXG + c.sg] = iumax; c.redI[2 * MAXG + c.sg] = ivmin; c.redI[3 * MAXG + c.sg] = ivmax;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    umin = INFINITY; umax = -INFINITY; vmin = INFINITY; vmax = -INFINITY; bad = 0.0f;
    iumin = 0xFFFFFFFFu; iumax = 0xFFFFFFFFu; ivmin = 0xFFFFFFFFu; ivmax = 0xFFFFFFFFu;
    for (uint g = 0; g < c.G; ++g) {
        argMinComb(umin, iumin, c.red[8 * MAXG + g], c.redI[0 * MAXG + g]);
        argMaxComb(umax, iumax, c.red[9 * MAXG + g], c.redI[1 * MAXG + g]);
        argMinComb(vmin, ivmin, c.red[10 * MAXG + g], c.redI[2 * MAXG + g]);
        argMaxComb(vmax, ivmax, c.red[11 * MAXG + g], c.redI[3 * MAXG + g]);
        bad = max(bad, c.red[12 * MAXG + g]);
    }
    if (bad > 0.0f) { return 0; }
    if (kMicroOn && kMicro == 4) { return 1; }
    float const Cx = 0.5f * (umax + umin);
    float const Cy = 0.5f * (vmax + vmin);
    float const spanU = umax - umin, spanV = vmax - vmin;
    float const S = max(max(spanU, spanV), 1e-9f);
    if (!(S > 0.0f) || !isfinite(S)) { return 0; }
    float const invS = 1.0f / S;
    if (c.tid < n) {
        uint const k = c.tid;
        dCx[k] = (k < nx) ? 0.5f * (Jraw[iumax * n + k] + Jraw[iumin * n + k]) : 0.0f;
        dCy[k] = (k >= nx) ? 0.5f * (Jraw[ivmax * n + k] + Jraw[ivmin * n + k]) : 0.0f;
        float ds = 0.0f;
        if (spanU >= spanV && spanU > 1e-9f) { if (k < nx) { ds = Jraw[iumax * n + k] - Jraw[iumin * n + k]; } }
        else if (spanV > 1e-9f) { if (k >= nx) { ds = Jraw[ivmax * n + k] - Jraw[ivmin * n + k]; } }
        dS[k] = ds;
    }
    for (uint j = c.tid; j < c.M; j += c.T) { float2 v = c.cv[j]; c.cv[j] = float2((v.x - Cx) * invS, (v.y - Cy) * invS); }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint const nb = segBlocksG(c);

    uint const nt = (n + 1 + 7) / 8;
    uint const ncol = nt * 8;
    uint const m = c.N + c.M;
    uint const mpad = (m + 7u) & ~7u;
    float badR = 0.0f;
    float const fN = float(c.N), fM = float(c.M);
    for (uint i = c.tid; i < c.N; i += c.T) {
        device float* row = Jd + i * NCP;
        for (uint k = 0; k < ncol; ++k) { row[k] = 0.0f; }
        float2 pt = c.data[i];
        float px = pt.x, py = pt.y;
        uint jb; float tau;
        float best;
        if (kMicroOn && kMicro == 6) { jb = 0; tau = 0.5f; float2 a0 = c.cv[0]; best = (px - a0.x) * (px - a0.x) + (py - a0.y) * (py - a0.y); }
        else { best = nnSegG(c, px, py, nb, jb, tau); }
        float w = weightOfG(best, c.noise[i].z, c.noise[i].w, c.hasNoise);
        float scale = sqrt(max(w, 0.0f) / (2.0f * fN));
        float d = sqrt(best);
        float ri = scale * d;
        row[n] = ri;
        if (!isfinite(ri)) { badR = 1.0f; }
        if (d > 1e-12f) {
            float2 A2 = c.cv[jb], B2 = c.cv[jb + 1];
            float ax = A2.x, ay = A2.y, bx = B2.x, by = B2.y;
            float zx = ax + tau * (bx - ax), zy = ay + tau * (by - ay);
            float nhx = (px - zx) / d, nhy = (py - zy) / d;
            device const float* r0 = Jraw + jb * n;
            device const float* r1 = Jraw + (jb + 1) * n;
            for (uint k = 0; k < n; ++k) {
                float gu0 = (k < nx) ? r0[k] : 0.0f, gv0 = (k >= nx) ? r0[k] : 0.0f;
                float gu1 = (k < nx) ? r1[k] : 0.0f, gv1 = (k >= nx) ? r1[k] : 0.0f;
                float dqx0 = (gu0 - dCx[k]) * invS - ax * dS[k] * invS, dqy0 = (gv0 - dCy[k]) * invS - ay * dS[k] * invS;
                float dqx1 = (gu1 - dCx[k]) * invS - bx * dS[k] * invS, dqy1 = (gv1 - dCy[k]) * invS - by * dS[k] * invS;
                float g0 = nhx * dqx0 + nhy * dqy0;
                float g1 = nhx * dqx1 + nhy * dqy1;
                float val = -scale * ((1.0f - tau) * g0 + tau * g1);
                row[k] = val;
                if (!isfinite(val)) { badR = 1.0f; }
            }
        }
    }
    for (uint j = c.tid; j < c.M; j += c.T) {
        device float* row = Jd + (c.N + j) * NCP;
        for (uint k = 0; k < ncol; ++k) { row[k] = 0.0f; }
        float2 q = c.cv[j];
        float qx = q.x, qy = q.y;
        uint kb;
        float best;
        if (kMicroOn && kMicro == 6) { kb = 0; float2 d0 = c.data[0]; best = (qx - d0.x) * (qx - d0.x) + (qy - d0.y) * (qy - d0.y); }
        else { best = nnPtG(c, qx, qy, kb); }
        float w = weightOfG(best, c.noise[kb].z, c.noise[kb].w, c.hasNoise);
        float scale = sqrt(max(w, 0.0f) / (2.0f * fM));
        float e = sqrt(best);
        float rj = scale * e;
        row[n] = rj;
        if (!isfinite(rj)) { badR = 1.0f; }
        if (e > 1e-12f) {
            float2 pk = c.data[kb];
            float mhx = (qx - pk.x) / e, mhy = (qy - pk.y) / e;
            device const float* rr = Jraw + j * n;
            for (uint k = 0; k < n; ++k) {
                float gu = (k < nx) ? rr[k] : 0.0f, gv = (k >= nx) ? rr[k] : 0.0f;
                float dqx = (gu - dCx[k]) * invS - qx * dS[k] * invS;
                float dqy = (gv - dCy[k]) * invS - qy * dS[k] * invS;
                float val = scale * (mhx * dqx + mhy * dqy);
                row[k] = val;
                if (!isfinite(val)) { badR = 1.0f; }
            }
        }
    }
    for (uint r = m + c.tid; r < mpad; r += c.T) { device float* row = Jd + r * NCP; for (uint k = 0; k < ncol; ++k) { row[k] = 0.0f; } }
    badR = simd_max(badR);
    if (c.lane == 0) { c.red[13 * MAXG + c.sg] = badR; }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    float anyBad = 0.0f;
    for (uint g = 0; g < c.G; ++g) { anyBad = max(anyBad, c.red[13 * MAXG + g]); }
    if (anyBad > 0.0f) { return 0; }
    uint const npairs = nt * (nt + 1) / 2;
    for (uint q = c.sg; q < npairs; q += c.G) {
        uint a = 0, rem = q;
        while (rem >= nt - a) { rem -= nt - a; ++a; }
        uint const b = a + rem;
        simdgroup_float8x8 acc = make_filled_simdgroup_matrix<float, 8, 8>(0.0f);
        for (uint r0 = 0; r0 < mpad; r0 += 8) {
            simdgroup_float8x8 A, B;
            simdgroup_load(A, Jd, NCP, ulong2(a * 8, r0), true);
            simdgroup_load(B, Jd, NCP, ulong2(b * 8, r0), false);
            simdgroup_multiply_accumulate(acc, A, B, acc);
        }
        simdgroup_store(acc, Hs, NCP, ulong2(b * 8, a * 8));
        if (a != b) { simdgroup_store(acc, Hs, NCP, ulong2(a * 8, b * 8), true); }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return 1;
}

kernel void kCclmG(device const GInst* insts [[buffer(0)]],
                   device const GCand* cands [[buffer(1)]],
                   device float* coeffs [[buffer(2)]],
                   device const float* tvals [[buffer(3)]],
                   device const float2* data [[buffer(4)]],
                   device const float4* noise [[buffer(5)]],
                   constant GParams& P [[buffer(6)]],
                   device float* scratch [[buffer(7)]],
                   device int* status [[buffer(8)]],
                   device float* lossOut [[buffer(9)]],
                   device const float4* dbox [[buffer(10)]],
                   device atomic_int* doneFlag [[buffer(11)]],
                   device atomic_uint* resOut [[buffer(12)]],
                   uint p [[threadgroup_position_in_grid]],
                   uint tid [[thread_position_in_threadgroup]],
                   uint T [[threads_per_threadgroup]],
                   uint lane [[thread_index_in_simdgroup]],
                   uint sg [[simdgroup_index_in_threadgroup]],
                   uint nsg [[simdgroups_per_threadgroup]])
{
    threadgroup float xs[NMAX], xn[NMAX], dCx[NMAX], dCy[NMAX], dS[NMAX];
    threadgroup float2 cv[MAXP];
    threadgroup float4 cbox[MAXB];
    threadgroup float red[MAXG * 16];
    threadgroup uint redI[MAXG * 4];
    threadgroup float Hs[NCP * NCP];
    threadgroup float Ls[NMAX * NCP];
    threadgroup int go[2];

    GCand C = cands[p];
    uint const n = C.n, nx = C.nx;
    GC c;
    c.ix = insts + C.offX; c.lx = C.lenX; c.iy = insts + C.offY; c.ly = C.lenY;
    c.tv = tvals; c.data = data; c.noise = noise; c.dbox = dbox; c.nbD = P.pad0;
    c.N = P.N; c.M = P.M; c.hasNoise = P.hasNoise; c.ptClamp = P.ptClamp;
    c.tid = tid; c.T = T; c.lane = lane; c.sg = sg; c.G = nsg;
    c.cv = cv; c.cbox = cbox; c.red = red; c.redI = redI;
    device float* Jraw = scratch + C.scrJraw;
    device float* Jd = scratch + C.scrJ;
    if (tid < n) { xs[tid] = coeffs[p * NMAX + tid]; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (kMicroOn) {
        float accm = 0.0f;
        uint const R = P.maxIter;
        for (uint r = 0; r < R; ++r) {
            if (kMicro == 1 || kMicro == 2) { accm += lossG(c, xs); }
            else if (kMicro == 3 || kMicro == 4 || kMicro == 6) { accm += float(buildG(c, nx, n, xs, dCx, dCy, dS, Hs, Jraw, Jd)); }
            else if (kMicro == 5 && sg == 0) {
                float hm[NMAX];
                for (uint k = 0; k < NMAX; ++k) { hm[k] = (lane < n && k < n) ? ((k == lane) ? 4.0f : 0.5f / float(1 + k + lane)) : 0.0f; }
                float dx = 0.0f;
                accm += cholS(hm, 1.0f, 1.0f, 1e-3f * float(r + 1), n, lane, Ls, dx) ? dx : 0.0f;
            }
        }
        if (tid == 0) { lossOut[p] = accm; status[p * 4 + 0] = 0; status[p * 4 + 1] = 0; status[p * 4 + 2] = 0; status[p * 4 + 3] = 0; }
        return;
    }

    int iters = 0, lossCalls = 0, builds = 0;
    float f = lossG(c, xs);
    ++lossCalls;
    bool done = !isfinite(f) || (P.floorThr > 0.0f && f <= P.floorThr);
    float lambda = P.lambda0;
    for (uint it = 0; it < P.maxIter && !done; ++it) {
        ++builds;
        int b = buildG(c, nx, n, xs, dCx, dCy, dS, Hs, Jraw, Jd);
        if (b == 0) { break; }
        float h[NMAX];
        float gi = 0.0f, dg = 1.0f;
        if (sg == 0) {
            for (uint k = 0; k < NMAX; ++k) { h[k] = (lane < n && k < n) ? Hs[lane * NCP + k] : 0.0f; }
            gi = (lane < n) ? Hs[lane * NCP + n] : 0.0f;
            float const hii = (lane < n) ? Hs[lane * NCP + lane] : 1.0f;
            dg = (hii > 0.0f) ? hii : 1.0f;
        }
        bool accepted = false;
        float rel = 0.0f;
        for (uint trial = 0; trial < P.maxTrials; ++trial) {
            if (sg == 0) {
                float dx = 0.0f;
                bool ok = cholS(h, gi, dg, lambda, n, lane, Ls, dx);
                if (ok) {
                    bool fin = true;
                    if (lane < n) { float v = xs[lane] + dx; fin = isfinite(v); xn[lane] = v; }
                    ok = simd_all(fin);
                }
                if (lane == 0) { go[trial & 1u] = ok ? 1 : 0; }
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            int const g = go[trial & 1u];
            if (g == 0) { lambda = min(lambda * P.lambdaUp, P.lambdaMax); continue; }
            float fn = lossG(c, xn);
            ++lossCalls;
            if (isfinite(fn) && fn < f) {
                rel = (f - fn) / max(fabs(f), 1e-30f);
                if (tid < n) { xs[tid] = xn[tid]; }
                threadgroup_barrier(mem_flags::mem_threadgroup);
                f = fn;
                lambda = max(lambda * P.lambdaDown, P.lambdaMin);
                accepted = true;
                break;
            }
            lambda = min(lambda * P.lambdaUp, P.lambdaMax);
            if (P.screening != 0u && it == 0 && trial == 0 && f > 0.05f) { break; }
        }
        ++iters;
        if (!accepted) { break; }
        if (rel < P.ftol) { break; }
    }
    if (tid == 0) {
        device atomic_uint* r = resOut + p * RSTRIDE;
        uint w[NMAX + 6];
        for (uint k = 0; k < NMAX; ++k) { w[k] = (k < n) ? as_type<uint>(xs[k]) : 0u; }
        w[NMAX + 0] = 0u; w[NMAX + 1] = uint(iters); w[NMAX + 2] = uint(lossCalls); w[NMAX + 3] = uint(builds);
        w[NMAX + 4] = as_type<uint>(f); w[NMAX + 5] = P.pad1u;
        uint h = 2166136261u ^ p;
        for (uint k = 0; k < NMAX + 6; ++k) { h = (h ^ w[k]) * 16777619u; atomic_store_explicit(&r[k], w[k], memory_order_relaxed); }
        atomic_store_explicit(&r[NMAX + 6], h, memory_order_relaxed);
        atomic_thread_fence(mem_flags::mem_device, memory_order_seq_cst, thread_scope_device);
        atomic_store_explicit(&doneFlag[p], 1, memory_order_relaxed);
        coeffs[p * NMAX] = coeffs[p * NMAX];   // (no-op; keeps the binding used)
        status[p * 4 + 0] = 0; status[p * 4 + 1] = iters; status[p * 4 + 2] = lossCalls; status[p * 4 + 3] = builds;
        lossOut[p] = f;
    }
}
kernel void kSpin(device float* out [[buffer(0)]], constant uint& iters [[buffer(1)]],
                  uint gid [[thread_position_in_grid]]) {
    float a = float(gid) * 1e-3f, b = 1.0001f;
    for (uint i = 0; i < iters; ++i) { a = fma(a, b, 1e-7f); b = fma(b, 0.9999f, 1e-4f); }
    if (a == 12345.678f) { out[0] = a + b; }
}
)MSL";

inline auto cclmMicroPipelines() -> bool& { static bool b = false; return b; }

inline auto cclmHostSeconds() -> double {
    static double const k = [] { mach_timebase_info_data_t tb; mach_timebase_info(&tb); return static_cast<double>(tb.numer) / static_cast<double>(tb.denom) * 1e-9; }();
    return static_cast<double>(mach_absolute_time()) * k;
}

struct CclmGpuRuntime {
    id<MTLDevice> dev{nil};
    id<MTLCommandQueue> queue{nil};
    id<MTLComputePipelineState> pso{nil};
    CclmTiming last{};
    double submitHostS{0}, commitHostS{0};
    id<MTLBuffer> bInsts{nil}, bCands{nil}, bCoeffs{nil}, bT{nil}, bData{nil}, bNoise{nil}, bParams{nil},
        bScratch{nil}, bStatus{nil}, bLoss{nil};
    id<MTLCommandBuffer> inflight{nil};
    std::size_t inflightP{0};
    std::vector<uint32_t> inflightN;
    std::shared_ptr<std::atomic<int64_t>> doneNs{std::make_shared<std::atomic<int64_t>>(0)};
    std::size_t allocEvents{0};

    auto ensure(__strong id<MTLBuffer>& b, std::size_t need) -> bool {
        std::size_t const cap = (b == nil) ? 0 : static_cast<std::size_t>([b length]);
        if (b != nil && need <= cap) { return true; }
        std::size_t nc = std::max<std::size_t>({need, cap * 2, std::size_t{256}});
        b = [dev newBufferWithLength:nc options:MTLResourceStorageModeShared];
        if (b == nil) { return false; }
        ++allocEvents;
        return true;
    }

    static constexpr int kNumProd = 5;
    static constexpr int kNumVariants = kNumProd + 24;
    static constexpr int kDefaultVariant = 4;
    static auto spinPsoRef() -> __strong id<MTLComputePipelineState>& { static id<MTLComputePipelineState> s = nil; return s; }
    static auto sharedPipelines(__strong id<MTLDevice>& devOut) -> std::vector<id<MTLComputePipelineState>> const& {
        static id<MTLDevice> sDev = nil;
        static std::vector<id<MTLComputePipelineState>> sPso;
        static std::once_flag once;
        std::call_once(once, [] {
            @autoreleasepool {
                sPso.assign(kNumVariants, nil);
                sDev = MTLCreateSystemDefaultDevice();
                if (sDev == nil) { return; }
                NSError* err = nil;
                MTLCompileOptions* opts = [MTLCompileOptions new];
                if (@available(macOS 15.0, iOS 18.0, *)) { opts.mathMode = MTLMathModeSafe; } else { opts.fastMathEnabled = NO; }
                id<MTLLibrary> lib = [sDev newLibraryWithSource:[NSString stringWithUTF8String:kCclmShader] options:opts error:&err];
                if (lib == nil) { NSLog(@"[cclm] shader compile failed: %@", err); return; }
                auto make = [&](NSString* name, int bb, int micro = 0) -> id<MTLComputePipelineState> {
                    id<MTLFunction> fn = nil;
                    if (bb < 0) { fn = [lib newFunctionWithName:name]; }
                    else {
                        MTLFunctionConstantValues* cv = [MTLFunctionConstantValues new];
                        bool v = bb != 0;
                        [cv setConstantValue:&v type:MTLDataTypeBool atIndex:0];
                        int mi = micro;
                        [cv setConstantValue:&mi type:MTLDataTypeInt atIndex:1];
                        NSError* e2 = nil;
                        fn = [lib newFunctionWithName:name constantValues:cv error:&e2];
                        if (fn == nil) { NSLog(@"[cclm] function %@ failed: %@", name, e2); }
                    }
                    if (fn == nil) { return nil; }
                    NSError* e3 = nil;
                    id<MTLComputePipelineState> ps = [sDev newComputePipelineStateWithFunction:fn error:&e3];
                    if (ps == nil) { NSLog(@"[cclm] pipeline %@ failed: %@", name, e3); }
                    return ps;
                };
                sPso[0] = make(@"kCclm", -1);
                sPso[1] = make(@"kCclmS", 0);
                sPso[2] = make(@"kCclmS", 1);
                sPso[3] = make(@"kCclmG", 0);
                sPso[4] = make(@"kCclmG", 1);
                spinPsoRef() = make(@"kSpin", -1);
                for (int mi = 1; mi <= 6 && cclmMicroPipelines(); ++mi) {
                    for (int bb = 0; bb < 2; ++bb) {
                        sPso[static_cast<std::size_t>(kNumProd + 2 * (mi - 1) + bb)] = make(@"kCclmS", bb, mi);
                        sPso[static_cast<std::size_t>(kNumProd + 12 + 2 * (mi - 1) + bb)] = make(@"kCclmG", bb, mi);
                    }
                }
            }
        });
        devOut = sDev;
        return sPso;
    }

    int variant{kDefaultVariant};
    int forceG{0};
    std::size_t lastThreadsPerGroup{0};
    std::vector<id<MTLComputePipelineState>> psos;
    id<MTLBuffer> bDBox{nil};
    id<MTLBuffer> bDone{nil};
    id<MTLBuffer> bResOut{nil};
    bool inflightStreams{false};

    auto init() -> bool {
        if (pso != nil) { return true; }
        psos = sharedPipelines(dev);
        if (dev == nil || psos.size() != kNumVariants) { return false; }
        for (int i = 0; i < kNumProd; ++i) { if (psos[static_cast<std::size_t>(i)] == nil) { return false; } }
        pso = psos[0];
        queue = [dev newCommandQueue];
        bParams = [dev newBufferWithLength:sizeof(CclmParams) options:MTLResourceStorageModeShared];
        return queue != nil && bParams != nil;
    }

    auto reserve(std::size_t P) -> bool {
        if (!init() || P == 0) { return false; }
        std::size_t const NM = kCclmMaxPoints, n = kCclmMaxCoeffs;
        return ensure(bInsts, P * 2 * kCclmMaxTape * sizeof(CclmInst))
            && ensure(bCands, P * sizeof(CclmCandDesc))
            && ensure(bCoeffs, P * kCclmMaxCoeffs * sizeof(float))
            && ensure(bStatus, P * 4 * sizeof(int32_t))
            && ensure(bLoss, P * sizeof(float))
            && ensure(bT, NM * sizeof(float))
            && ensure(bData, 2 * NM * sizeof(float))
            && ensure(bNoise, 4 * NM * sizeof(float))
            && ensure(bDBox, (NM / 16 + 1) * 4 * sizeof(float))
            && ensure(bDone, P * sizeof(int32_t))
            && ensure(bResOut, P * kCclmResStride * sizeof(uint32_t))
            && ensure(bScratch, P * (NM * n + cclmJRegion(NM, NM, n) + 2 * NM) * sizeof(float));
    }

    auto submit(CclmFitData const& fit, std::vector<CclmJob> const& jobs, std::size_t P) -> bool {
        if (!init() || P == 0 || P > jobs.size() || inflight != nil) { return false; }
        submitHostS = cclmHostSeconds();
        if (cclmCaptureHook().fn != nullptr) { cclmCaptureHook().fn(fit, jobs, P, cclmCaptureHook().user); }
        std::size_t const N = fit.params.N, M = fit.params.M;
        std::size_t totalInsts = 0, totalScratch = 0;
        for (std::size_t p = 0; p < P; ++p) {
            auto const& j = jobs[p];
            totalInsts += j.tapeX.size() + j.tapeY.size();
            std::size_t const n = j.coeffs.size();
            totalScratch += M * n + cclmJRegion(N, M, n) + (N + M);
        }
        if (!ensure(bInsts, totalInsts * sizeof(CclmInst))) { return false; }
        if (!ensure(bCands, P * sizeof(CclmCandDesc))) { return false; }
        if (!ensure(bCoeffs, P * kCclmMaxCoeffs * sizeof(float))) { return false; }
        if (!ensure(bStatus, P * 4 * sizeof(int32_t))) { return false; }
        if (!ensure(bLoss, P * sizeof(float))) { return false; }
        if (!ensure(bT, M * sizeof(float))) { return false; }
        if (!ensure(bData, 2 * N * sizeof(float))) { return false; }
        if (!ensure(bNoise, 4 * N * sizeof(float))) { return false; }
        if (!ensure(bScratch, totalScratch * sizeof(float))) { return false; }
        std::size_t const nbD = (N + 15) / 16;
        if (!ensure(bDBox, (nbD + 1) * 4 * sizeof(float))) { return false; }
        if (!ensure(bDone, P * sizeof(int32_t))) { return false; }
        if (!ensure(bResOut, P * kCclmResStride * sizeof(uint32_t))) { return false; }
        std::memset([bResOut contents], 0xFF, P * kCclmResStride * sizeof(uint32_t));
        std::memset([bDone contents], 0, P * sizeof(int32_t));

        auto* ins = static_cast<CclmInst*>([bInsts contents]);
        auto* cd = static_cast<CclmCandDesc*>([bCands contents]);
        auto* cf = static_cast<float*>([bCoeffs contents]);
        inflightN.resize(P);
        std::size_t io = 0, so = 0;
        for (std::size_t p = 0; p < P; ++p) {
            auto const& j = jobs[p];
            CclmCandDesc d;
            d.offX = static_cast<uint32_t>(io);
            d.lenX = static_cast<uint32_t>(j.tapeX.size());
            std::memcpy(ins + io, j.tapeX.data(), j.tapeX.size() * sizeof(CclmInst)); io += j.tapeX.size();
            d.offY = static_cast<uint32_t>(io);
            d.lenY = static_cast<uint32_t>(j.tapeY.size());
            std::memcpy(ins + io, j.tapeY.data(), j.tapeY.size() * sizeof(CclmInst)); io += j.tapeY.size();
            std::size_t const n = j.coeffs.size();
            d.nx = j.nx;
            d.n = static_cast<uint32_t>(n);
            d.scrJraw = static_cast<uint32_t>(so); so += M * n;
            d.scrJ = static_cast<uint32_t>(so); so += cclmJRegion(N, M, n);
            d.scrR = static_cast<uint32_t>(so); so += (N + M);
            cd[p] = d;
            std::memcpy(cf + p * kCclmMaxCoeffs, j.coeffs.data(), n * sizeof(float));
            inflightN[p] = static_cast<uint32_t>(n);
        }
        std::memcpy([bT contents], fit.t.data(), M * sizeof(float));
        std::memcpy([bData contents], fit.dataXY.data(), 2 * N * sizeof(float));
        std::memcpy([bNoise contents], fit.noise4.data(), 4 * N * sizeof(float));
        {
            CclmParams prm = fit.params;
            prm.pad0 = static_cast<uint32_t>(nbD);
            ++dispatchSeq;
            std::memcpy(&prm.pad1, &dispatchSeq, sizeof dispatchSeq);
            std::memcpy([bParams contents], &prm, sizeof(CclmParams));
            auto* db = static_cast<float*>([bDBox contents]);
            for (std::size_t k = 0; k < nbD; ++k) {
                std::size_t const lo = k * 16, hi = std::min(lo + 16, N);
                float xlo = fit.dataXY[2 * lo], xhi = xlo, ylo = fit.dataXY[2 * lo + 1], yhi = ylo;
                for (std::size_t i = lo + 1; i < hi; ++i) {
                    xlo = std::min(xlo, fit.dataXY[2 * i]); xhi = std::max(xhi, fit.dataXY[2 * i]);
                    ylo = std::min(ylo, fit.dataXY[2 * i + 1]); yhi = std::max(yhi, fit.dataXY[2 * i + 1]);
                }
                db[4 * k + 0] = xlo; db[4 * k + 1] = ylo; db[4 * k + 2] = xhi; db[4 * k + 3] = yhi;
            }
        }

        int vv = (variant >= 0 && variant < kNumVariants) ? variant : kDefaultVariant;
        if (psos[static_cast<std::size_t>(vv)] == nil) { vv = kDefaultVariant; }
        id<MTLComputePipelineState> ps = psos[static_cast<std::size_t>(vv)];
        NSUInteger want = 32;
        if (vv == 0) {
            NSUInteger const maxT = [ps maxTotalThreadsPerThreadgroup];
            NSUInteger const w = [ps threadExecutionWidth];
            want = static_cast<NSUInteger>(std::max(N, M));
            want = ((want + w - 1) / w) * w;
            want = std::min<NSUInteger>(std::min<NSUInteger>(want, 256), maxT);
            want = std::max<NSUInteger>(want, w);
        } else if (vv == 3 || vv == 4 || vv >= kNumProd + 12) {
            NSUInteger const maxG = std::max<NSUInteger>(1, std::min<NSUInteger>(8, [ps maxTotalThreadsPerThreadgroup] / 32));
            NSUInteger g = (static_cast<NSUInteger>(std::max(N, M)) + 31) / 32;
            if (forceG > 0) { g = static_cast<NSUInteger>(forceG); }
            want = 32 * std::max<NSUInteger>(1, std::min(g, maxG));
        }
        lastThreadsPerGroup = want;
        inflightStreams = (vv == 3 || vv == 4);

        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:ps];
        [enc setBuffer:bInsts offset:0 atIndex:0];
        [enc setBuffer:bCands offset:0 atIndex:1];
        [enc setBuffer:bCoeffs offset:0 atIndex:2];
        [enc setBuffer:bT offset:0 atIndex:3];
        [enc setBuffer:bData offset:0 atIndex:4];
        [enc setBuffer:bNoise offset:0 atIndex:5];
        [enc setBuffer:bParams offset:0 atIndex:6];
        [enc setBuffer:bScratch offset:0 atIndex:7];
        [enc setBuffer:bStatus offset:0 atIndex:8];
        [enc setBuffer:bLoss offset:0 atIndex:9];
        [enc setBuffer:bDBox offset:0 atIndex:10];
        [enc setBuffer:bDone offset:0 atIndex:11];
        [enc setBuffer:bResOut offset:0 atIndex:12];
        [enc dispatchThreadgroups:MTLSizeMake(P, 1, 1) threadsPerThreadgroup:MTLSizeMake(want, 1, 1)];
        [enc endEncoding];
        doneNs->store(0);
        auto done = doneNs;
        [cb addCompletedHandler:^(id<MTLCommandBuffer>) {
            done->store(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now().time_since_epoch()).count());
        }];
        [cb commit];
        commitHostS = cclmHostSeconds();
        inflight = cb;
        inflightP = P;
        return true;
    }

    id<MTLBuffer> bSpin{nil};
    auto pulse(uint32_t iters, uint32_t tgs) -> id<MTLCommandBuffer> {
        if (!init() || spinPsoRef() == nil) { return nil; }
        if (bSpin == nil) { bSpin = [dev newBufferWithLength:16 options:MTLResourceStorageModeShared]; if (bSpin == nil) { return nil; } }
        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:spinPsoRef()];
        [enc setBuffer:bSpin offset:0 atIndex:0];
        [enc setBytes:&iters length:sizeof(iters) atIndex:1];
        [enc dispatchThreadgroups:MTLSizeMake(std::max<uint32_t>(1, tgs), 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
        [enc endEncoding];
        [cb commit];
        return cb;
    }

    auto wait(std::vector<CclmOut>& out, int64_t& completionNs, double& gpuKernelMs) -> bool {
        if (inflight == nil) { return false; }
        [inflight waitUntilCompleted];
        double const waitDone = cclmHostSeconds();
        bool const ok = [inflight status] == MTLCommandBufferStatusCompleted;
        gpuKernelMs = ([inflight GPUEndTime] - [inflight GPUStartTime]) * 1e3;
        {
            CclmTiming& t = last;
            t = CclmTiming{};
            t.submitHostS = submitHostS; t.commitHostS = commitHostS;
            t.gpuStartS = [inflight GPUStartTime]; t.gpuEndS = [inflight GPUEndTime]; t.waitDoneHostS = waitDone;
            double const ks = [inflight kernelStartTime], ke = [inflight kernelEndTime];
            t.encodeMs = (commitHostS - submitHostS) * 1e3;
            t.swQueueMs = (ks - commitHostS) * 1e3;
            t.driverMs = (ke - ks) * 1e3;
            t.hwQueueMs = (t.gpuStartS - ke) * 1e3;
            t.execMs = (t.gpuEndS - t.gpuStartS) * 1e3;
            t.notifyMs = (waitDone - t.gpuEndS) * 1e3;
        }
        while (doneNs->load() == 0 && ok) { std::this_thread::yield(); }
        completionNs = doneNs->load();
        inflight = nil;
        if (!ok) { return false; }
        if (inflightStreams) {
            out.resize(inflightP);
            for (std::size_t p = 0; p < inflightP; ++p) { if (!readRecord(p, out[p])) { out[p] = CclmOut{}; out[p].status = 1; } }
            double const rbDone0 = cclmHostSeconds();
            last.readbackMs = (rbDone0 - waitDone) * 1e3;
            last.totalMs = (rbDone0 - submitHostS) * 1e3;
            return true;
        }
        auto const* cf = static_cast<float const*>([bCoeffs contents]);
        auto const* stv = static_cast<int32_t const*>([bStatus contents]);
        auto const* lo = static_cast<float const*>([bLoss contents]);
        out.resize(inflightP);
        for (std::size_t p = 0; p < inflightP; ++p) {
            auto& o = out[p];
            o.status = stv[p * 4 + 0];
            o.iters = stv[p * 4 + 1];
            o.lossCalls = stv[p * 4 + 2];
            o.builds = stv[p * 4 + 3];
            o.loss = lo[p];
            o.coeffs.assign(cf + p * kCclmMaxCoeffs, cf + p * kCclmMaxCoeffs + inflightN[p]);
        }
        double const rbDone = cclmHostSeconds();
        last.readbackMs = (rbDone - waitDone) * 1e3;
        last.totalMs = (rbDone - submitHostS) * 1e3;
        return true;
    }

    auto streams() const -> bool { return inflightStreams; }
    uint32_t dispatchSeq{0};
    auto readRecord(std::size_t q, CclmOut& o) const -> bool {
        auto* r = static_cast<uint32_t*>([bResOut contents]) + q * kCclmResStride;
        constexpr std::size_t NW = kCclmMaxCoeffs + 6;
        uint32_t w[NW + 1];
        for (std::size_t k = 0; k <= NW; ++k) { w[k] = __atomic_load_n(r + k, __ATOMIC_ACQUIRE); }
        if (w[kCclmMaxCoeffs + 5] != dispatchSeq) { return false; }
        uint32_t h = 2166136261u ^ static_cast<uint32_t>(q);
        for (std::size_t k = 0; k < NW; ++k) { h = (h ^ w[k]) * 16777619u; }
        if (h != w[NW]) { return false; }
        auto asF = [](uint32_t u) { float f; std::memcpy(&f, &u, sizeof f); return f; };
        o.coeffs.resize(inflightN[q]);
        for (std::size_t k = 0; k < inflightN[q]; ++k) { o.coeffs[k] = asF(w[k]); }
        o.status = static_cast<int>(w[kCclmMaxCoeffs + 0]);
        o.iters = static_cast<int>(w[kCclmMaxCoeffs + 1]);
        o.lossCalls = static_cast<int>(w[kCclmMaxCoeffs + 2]);
        o.builds = static_cast<int>(w[kCclmMaxCoeffs + 3]);
        o.loss = asF(w[kCclmMaxCoeffs + 4]);
        return true;
    }
    auto pollOne(std::size_t q, CclmOut& o) const -> bool {
        if (bDone == nil || q >= inflightP) { return false; }
        auto* fl = static_cast<int32_t*>([bDone contents]);
        if (__atomic_load_n(fl + q, __ATOMIC_ACQUIRE) == 0) { return false; }
        return readRecord(q, o);
    }
    std::mutex collectMu;
    std::atomic<bool> collected{false};
    int64_t collectedDoneNs{0};
    auto tryCollect(std::vector<CclmOut>& out) -> bool {
        if (collected.load(std::memory_order_acquire)) { return true; }
        if (doneNs->load() == 0) { return false; }
        std::lock_guard<std::mutex> lk(collectMu);
        if (collected.load(std::memory_order_relaxed)) { return true; }
        if (inflight == nil) { return false; }
        double km = 0.0;
        bool const ok = wait(out, collectedDoneNs, km);
        if (!ok) { out.clear(); }
        collected.store(true, std::memory_order_release);
        return true;
    }
    void beginRace() { collected.store(false); collectedDoneNs = 0; }
    auto discardPending() -> double {
        if (inflight == nil) { return 0.0; }
        double const t0 = cclmHostSeconds();
        std::vector<CclmOut> sink;
        int64_t ns = 0; double km = 0.0;
        (void)wait(sink, ns, km);
        return (cclmHostSeconds() - t0) * 1e3;
    }
    id<MTLCommandQueue> pulseQueue{nil};
    auto pulseAsync(uint32_t iters, uint32_t tgs) -> bool {
        if (!init() || spinPsoRef() == nil) { return false; }
        if (pulseQueue == nil) { pulseQueue = [dev newCommandQueue]; if (pulseQueue == nil) { return false; } }
        if (bSpin == nil) { bSpin = [dev newBufferWithLength:16 options:MTLResourceStorageModeShared]; if (bSpin == nil) { return false; } }
        id<MTLCommandBuffer> cb = [pulseQueue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:spinPsoRef()];
        [enc setBuffer:bSpin offset:0 atIndex:0];
        [enc setBytes:&iters length:sizeof(iters) atIndex:1];
        [enc dispatchThreadgroups:MTLSizeMake(std::max<uint32_t>(1, tgs), 1, 1) threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
        [enc endEncoding];
        [cb commit];
        return true;
    }
};

inline auto cclmTryCollect(CclmGpuRuntime& rt, std::vector<CclmOut>& out) -> bool { return rt.tryCollect(out); }
inline auto cclmStreams(CclmGpuRuntime const& rt) -> bool { return rt.streams(); }
inline auto cclmPollOne(CclmGpuRuntime const& rt, std::size_t q, CclmOut& o) -> bool { return rt.pollOne(q, o); }
inline void cclmBeginRace(CclmGpuRuntime& rt) { rt.beginRace(); }
inline auto cclmDiscardPending(CclmGpuRuntime& rt) -> double { return rt.discardPending(); }
inline auto cclmCollectedDoneNs(CclmGpuRuntime const& rt) -> int64_t { return rt.collectedDoneNs; }
inline auto cclmPulse(CclmGpuRuntime& rt, uint32_t iters, uint32_t tgs) -> bool { return rt.pulseAsync(iters, tgs); }

inline auto cclmMakeRuntime() -> std::shared_ptr<CclmGpuRuntime> {
    auto rt = std::make_shared<CclmGpuRuntime>();
    if (!rt->init()) { return nullptr; }
    return rt;
}
inline auto cclmSubmit(CclmGpuRuntime& rt, CclmFitData const& fit, std::vector<CclmJob> const& jobs, std::size_t P) -> bool {
    return rt.submit(fit, jobs, P);
}
inline auto cclmWait(CclmGpuRuntime& rt, std::vector<CclmOut>& out, int64_t& completionNs, double& gpuKernelMs) -> bool {
    return rt.wait(out, completionNs, gpuKernelMs);
}
inline auto cclmAllocEvents(CclmGpuRuntime const& rt) -> std::size_t { return rt.allocEvents; }
inline auto cclmReserve(CclmGpuRuntime& rt, std::size_t P) -> bool { return rt.reserve(P); }
inline void cclmSetVariant(CclmGpuRuntime& rt, int v) { rt.variant = v; }
inline auto cclmLastTiming(CclmGpuRuntime const& rt) -> CclmTiming { return rt.last; }

#else

struct CclmGpuRuntime {};
inline auto cclmMakeRuntime() -> std::shared_ptr<CclmGpuRuntime> { return nullptr; }
inline auto cclmSubmit(CclmGpuRuntime&, CclmFitData const&, std::vector<CclmJob> const&, std::size_t) -> bool { return false; }
inline auto cclmWait(CclmGpuRuntime&, std::vector<CclmOut>&, int64_t&, double&) -> bool { return false; }
inline auto cclmAllocEvents(CclmGpuRuntime const&) -> std::size_t { return 0; }
inline auto cclmReserve(CclmGpuRuntime&, std::size_t) -> bool { return false; }
inline void cclmSetVariant(CclmGpuRuntime&, int) {}
inline auto cclmLastTiming(CclmGpuRuntime const&) -> CclmTiming { return {}; }
inline auto cclmTryCollect(CclmGpuRuntime&, std::vector<CclmOut>&) -> bool { return false; }
inline auto cclmStreams(CclmGpuRuntime const&) -> bool { return false; }
inline auto cclmPollOne(CclmGpuRuntime const&, std::size_t, CclmOut&) -> bool { return false; }
inline void cclmBeginRace(CclmGpuRuntime&) {}
inline auto cclmDiscardPending(CclmGpuRuntime&) -> double { return 0.0; }
inline auto cclmCollectedDoneNs(CclmGpuRuntime const&) -> int64_t { return 0; }
inline auto cclmPulse(CclmGpuRuntime&, uint32_t, uint32_t) -> bool { return false; }

#endif

}
