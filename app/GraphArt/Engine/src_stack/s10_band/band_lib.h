#pragma once
#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>

namespace band {

constexpr int NROW = 48, NLEV = 64, NSEG = 63, WCOR = 12, SLOTS = 25;

struct Tables {
    int32_t base[NROW];
    int32_t del[NROW];
    uint32_t valid[NROW];
};

static inline long floordiv_l(long a, long b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
static inline long ceildiv_l(long a, long b) { return -floordiv_l(-a, b); }

static inline bool corridor_ok(int s, int i) {
    if (s < 0 || s > 62) return false;
    long d = 47L * s - 63L * i;
    return d <= 564 && d >= -564;
}

static inline Tables make_tables() {
    Tables T{};
    for (int i = 0; i < NROW; i++) T.base[i] = (int32_t)ceildiv_l(63L * i - 564L, 47L);
    for (int i = 0; i < NROW; i++) {
        T.del[i] = (i == 0) ? 1 : T.base[i] - T.base[i - 1];
        uint32_t v = 0;
        for (int k = 0; k < SLOTS; k++)
            if (corridor_ok(T.base[i] + k, i)) v |= 1u << k;
        T.valid[i] = v;
    }
    return T;
}

static inline void make_stroke(float S[NROW][2]) {
    for (int i = 0; i < NROW; i++) {
        double s = (double)i / 47.0;
        S[i][0] = (float)std::cos(2.0 * M_PI * s);
        S[i][1] = (float)std::sin(2.0 * M_PI * s);
    }
}

static inline void gen_candidates_rng(std::mt19937& rng, float* aos, int P);

static inline void gen_candidates(float* aos, int P, uint32_t seed) {
    std::mt19937 rng(seed);
    gen_candidates_rng(rng, aos, P);
}

static inline void gen_candidates_rng(std::mt19937& rng, float* aos, int P) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int c = 0; c < P; c++) {
        double a1 = 0.3 * U(rng), a2 = 0.3 * U(rng);
        double f1 = 0.5 + 3.5 * U(rng), f2 = 0.5 + 3.5 * U(rng);
        double p1 = 2.0 * M_PI * U(rng), p2 = 2.0 * M_PI * U(rng);
        for (int l = 0; l < NLEV; l++) {
            double th = 2.0 * M_PI * (double)l / 63.0;
            aos[c * 128 + 2 * l + 0] = (float)(std::cos(th) + a1 * std::cos(f1 * th + p1));
            aos[c * 128 + 2 * l + 1] = (float)(std::sin(th) + a2 * std::sin(f2 * th + p2));
        }
    }
}

static inline void aos_to_soa(const float* aos, float* soa, int P) {
    for (int c = 0; c < P; c++) {
        int tile = c >> 5, lane = c & 31;
        for (int j = 0; j < 128; j++) soa[(size_t)tile * 4096 + (size_t)j * 32 + lane] = aos[(size_t)c * 128 + j];
    }
}

static const char* msl_source() {
    return R"MSL(
#include <metal_stdlib>
using namespace metal;

constant float BIGV = 1.0e30f;

// Mapping A: one candidate per thread. curve layout SoA: [tile][j][lane], j = 2*level + coord.
kernel void band_lane(device const float* curve [[buffer(0)]],
                      constant float2* stroke [[buffer(1)]],
                      constant int* tbase [[buffer(2)]],
                      constant int* tdel [[buffer(3)]],
                      constant uint* tvalid [[buffer(4)]],
                      device float* outT [[buffer(5)]],
                      device uchar* outS [[buffer(6)]],
                      device float* outObj [[buffer(7)]],
                      constant uint& P [[buffer(8)]],
                      uint gid [[thread_position_in_grid]])
{
    if (gid >= P) return;
    const uint tile = gid >> 5;
    const uint lane = gid & 31u;
    device const float* cb = curve + tile * 4096u + lane;
    float PM[25];
    #pragma unroll
    for (int k = 0; k < 25; k++) PM[k] = 0.0f;
    uint masks[48];
    for (int i = 0; i < 48; i++) {
        const int base = tbase[i];
        const int dl = tdel[i];
        const uint vm = tvalid[i];
        const float2 p = stroke[i];
        float2 q[26];
        #pragma unroll
        for (int j = 0; j < 26; j++) {
            const int idx = clamp(base + j, 0, 63);
            q[j] = float2(cb[idx * 64], cb[idx * 64 + 32]);
        }
        float D[25];
        #pragma unroll
        for (int k = 0; k < 25; k++) {
            const float2 A = q[k];
            const float2 E = q[k + 1] - A;
            const float L2 = max(E.x * E.x + E.y * E.y, 1e-30f);
            const float inv = 1.0f / L2;
            const float2 ex = E * inv;
            const float h = A.x * ex.x + A.y * ex.y;
            const float dt = fma(p.x, ex.x, fma(p.y, ex.y, -h));
            const float u = clamp(dt, 0.0f, 1.0f);
            const float qx = fma(-u, E.x, p.x - A.x);
            const float qy = fma(-u, E.y, p.y - A.y);
            float cost = fma(qx, qx, qy * qy);
            cost = ((vm >> k) & 1u) ? cost : BIGV;
            const float Qk = (dl == 1) ? PM[min(k + 1, 24)] : PM[min(k + 2, 24)];
            D[k] = cost + Qk;
        }
        uint mask = 1u;
        float run = D[0];
        PM[0] = D[0];
        #pragma unroll
        for (int m = 1; m < 25; m++) {
            if (D[m] < run) mask |= (1u << m);
            run = min(run, D[m]);
            PM[m] = run;
        }
        masks[i] = mask;
    }
    const float obj = PM[24];
    int kk = 31 - int(clz(masks[47]));
    for (int i = 47; i >= 0; i--) {
        const int s = clamp(tbase[i] + kk, 0, 62);
        const float2 A = float2(cb[s * 64], cb[s * 64 + 32]);
        const float2 B = float2(cb[(s + 1) * 64], cb[(s + 1) * 64 + 32]);
        const float2 p = stroke[i];
        const float2 E = B - A;
        const float L2 = max(E.x * E.x + E.y * E.y, 1e-30f);
        const float inv = 1.0f / L2;
        const float u = clamp(((p.x - A.x) * E.x + (p.y - A.y) * E.y) * inv, 0.0f, 1.0f);
        outT[(uint)i * P + gid] = (float(s) + u) * (1.0f / 63.0f);
        outS[(uint)i * P + gid] = uchar(s);
        if (i > 0) {
            const uint m = uint(min(kk + tdel[i], 24));
            const uint mk = masks[i - 1] & ((2u << m) - 1u);
            kk = 31 - int(clz(mk));
        }
    }
    outObj[gid] = obj;
}

// Mapping B: one candidate per SIMD-group (32 lanes = 25 DP slots + idle). curve AoS [cand][128].
kernel void band_simd(device const float* curve [[buffer(0)]],
                      constant float2* stroke [[buffer(1)]],
                      constant int* tbase [[buffer(2)]],
                      constant int* tdel [[buffer(3)]],
                      constant uint* tvalid [[buffer(4)]],
                      device float* outT [[buffer(5)]],
                      device uchar* outS [[buffer(6)]],
                      device float* outObj [[buffer(7)]],
                      constant uint& P [[buffer(8)]],
                      uint tg [[threadgroup_position_in_grid]],
                      uint sg [[simdgroup_index_in_threadgroup]],
                      uint lane [[thread_index_in_simdgroup]],
                      uint nsg [[simdgroups_per_threadgroup]])
{
    threadgroup float2 pts[8][64];
    threadgroup float4 segA[8][63];
    threadgroup float4 segB[8][63];
    threadgroup uint msk[8][48];
    const uint cand = tg * nsg + sg;
    if (cand >= P) return;
    device const float4* src = (device const float4*)(curve + (size_t)cand * 128);
    const float4 v = src[lane];
    pts[sg][2 * lane] = v.xy;
    pts[sg][2 * lane + 1] = v.zw;
    simdgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = lane; s < 63u; s += 32u) {
        const float2 A = pts[sg][s];
        const float2 E = pts[sg][s + 1] - A;
        const float L2 = max(E.x * E.x + E.y * E.y, 1e-30f);
        const float inv = 1.0f / L2;
        const float2 ex = E * inv;
        const float h = A.x * ex.x + A.y * ex.y;
        segA[sg][s] = float4(A.x, A.y, E.x, E.y);
        segB[sg][s] = float4(ex.x, ex.y, h, 0.0f);
    }
    simdgroup_barrier(mem_flags::mem_threadgroup);

    float PMv = 0.0f;
    for (int i = 0; i < 48; i++) {
        const int base = tbase[i];
        const int dl = tdel[i];
        const uint vm = tvalid[i];
        const bool valid = (lane < 25u) && (((vm >> lane) & 1u) != 0u);
        const int s = clamp(base + int(lane), 0, 62);
        const float4 a = segA[sg][s];
        const float4 b = segB[sg][s];
        const float2 p = stroke[i];
        const float dt = fma(p.x, b.x, fma(p.y, b.y, -b.z));
        const float u = clamp(dt, 0.0f, 1.0f);
        const float qx = fma(-u, a.z, p.x - a.x);
        const float qy = fma(-u, a.w, p.y - a.y);
        float cost = fma(qx, qx, qy * qy);
        cost = valid ? cost : BIGV;
        float Q = 0.0f;
        if (i > 0) Q = simd_shuffle(PMv, ushort(min(lane + uint(dl), 24u)));
        const float D = cost + Q;
        float v2 = D;
        for (ushort d = 1; d < 32; d <<= 1) {
            const float t2 = simd_shuffle_up(v2, d);
            v2 = (lane >= uint(d)) ? min(v2, t2) : v2;
        }
        const float pe = simd_shuffle_up(v2, ushort(1));
        const bool rec = (lane == 0u) || (D < pe);
        const uint mask = uint(simd_vote::vote_t(simd_ballot(rec)));
        if (lane == 0u) msk[sg][i] = mask;
        PMv = v2;
    }
    simdgroup_barrier(mem_flags::mem_threadgroup);
    const float objv = simd_shuffle(PMv, ushort(24));
    int kk = 31 - int(clz(msk[sg][47] & 0x1FFFFFFu));
    for (int i = 47; i >= 0; i--) {
        const int s = clamp(tbase[i] + kk, 0, 62);
        const float4 a = segA[sg][s];
        const float4 b = segB[sg][s];
        const float2 p = stroke[i];
        const float dt = fma(p.x, b.x, fma(p.y, b.y, -b.z));
        const float u = clamp(dt, 0.0f, 1.0f);
        if (lane == 0u) {
            outT[(uint)i * P + cand] = (float(s) + u) * (1.0f / 63.0f);
            outS[(uint)i * P + cand] = uchar(s);
        }
        if (i > 0) {
            const uint m = uint(min(kk + tdel[i], 24));
            kk = 31 - int(clz(msk[sg][i - 1] & ((2u << m) - 1u)));
        }
    }
    if (lane == 0u) outObj[cand] = objv;
}
)MSL";
}

struct Engine {
    id<MTLDevice> dev;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> pLane;
    id<MTLComputePipelineState> pSimd;
    Engine() {
        dev = MTLCreateSystemDefaultDevice();
        queue = [dev newCommandQueue];
        NSError* err = nil;
        MTLCompileOptions* o = [MTLCompileOptions new];
        id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:msl_source()] options:o error:&err];
        if (!lib) { fprintf(stderr, "compile error: %s\n", [[err description] UTF8String]); exit(2); }
        pLane = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"band_lane"] error:&err];
        if (!pLane) { fprintf(stderr, "pso lane: %s\n", [[err description] UTF8String]); exit(2); }
        pSimd = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"band_simd"] error:&err];
        if (!pSimd) { fprintf(stderr, "pso simd: %s\n", [[err description] UTF8String]); exit(2); }
    }
};

struct Bufs {
    id<MTLBuffer> curve, stroke, tb, td, tv, oT, oS, oO;
};

static inline id<MTLBuffer> new_buf(Engine& E, size_t bytes) {
    return [E.dev newBufferWithLength:bytes options:MTLResourceStorageModeShared];
}

static inline Bufs alloc_bufs(Engine& E, int P) {
    Bufs b;
    b.curve = new_buf(E, (size_t)P * 128 * sizeof(float));
    b.stroke = new_buf(E, NROW * sizeof(float) * 2);
    b.tb = new_buf(E, NROW * sizeof(int32_t));
    b.td = new_buf(E, NROW * sizeof(int32_t));
    b.tv = new_buf(E, NROW * sizeof(uint32_t));
    b.oT = new_buf(E, (size_t)P * NROW * sizeof(float));
    b.oS = new_buf(E, (size_t)P * NROW * sizeof(uint8_t));
    b.oO = new_buf(E, (size_t)P * sizeof(float));
    return b;
}

static inline void fill_tables(Bufs& b, const Tables& T, const float S[NROW][2]) {
    memcpy([b.stroke contents], S, NROW * sizeof(float) * 2);
    memcpy([b.tb contents], T.base, sizeof(T.base));
    memcpy([b.td contents], T.del, sizeof(T.del));
    memcpy([b.tv contents], T.valid, sizeof(T.valid));
}

static inline double run_kernel(Engine& E, bool simd, const Bufs& b, uint32_t P) {
    id<MTLCommandBuffer> cb = [E.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:(simd ? E.pSimd : E.pLane)];
    [enc setBuffer:b.curve offset:0 atIndex:0];
    [enc setBuffer:b.stroke offset:0 atIndex:1];
    [enc setBuffer:b.tb offset:0 atIndex:2];
    [enc setBuffer:b.td offset:0 atIndex:3];
    [enc setBuffer:b.tv offset:0 atIndex:4];
    [enc setBuffer:b.oT offset:0 atIndex:5];
    [enc setBuffer:b.oS offset:0 atIndex:6];
    [enc setBuffer:b.oO offset:0 atIndex:7];
    [enc setBytes:&P length:sizeof(P) atIndex:8];
    if (simd)
        [enc dispatchThreadgroups:MTLSizeMake(P / 8, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
    else
        [enc dispatchThreads:MTLSizeMake(P, 1, 1) threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];
    if (cb.error) { fprintf(stderr, "gpu error: %s\n", [[cb.error description] UTF8String]); exit(3); }
    return cb.GPUEndTime - cb.GPUStartTime;
}

}
