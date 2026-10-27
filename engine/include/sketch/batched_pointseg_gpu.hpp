#pragma once


#include "sketch/point_segment_kernel.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#if defined(__APPLE__) && defined(__OBJC__)
#include <chrono>
#include <mach/mach_time.h>
#endif
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>

#if defined(__APPLE__) && defined(__OBJC__)
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#endif

namespace Sketch {

inline auto pointSegFp16ModeFromEnv() -> int
{
    static int const v = [] {
        char const* e = std::getenv("SKETCH_PS_FP16");
        int const m = (e != nullptr) ? std::atoi(e) : 0;
        return (m == 1 || m == 2 || m == 3) ? m : 0;
    }();
    return v;
}

inline auto pointSegFp16ModeForFloatBuffers() -> int
{
    int const m = pointSegFp16ModeFromEnv();
    return (m == 2 || m == 3) ? m : 0;
}

struct BatchedPointSegResult {
    bool ok{false};
    std::vector<Operon::Scalar> fitness;
    std::vector<Operon::Scalar> residualFlat;
    std::vector<Operon::Scalar> tauFlat;
    int fp16ModeUsed{0};
};

struct BatchedPointSegReverseResult {
    bool ok{false};
    std::vector<Operon::Scalar> fitness;
};

auto evaluatePopulationPointSegReverseGpu(
    std::vector<Vec2> const& strokePointsDrawn,
    std::vector<std::vector<Vec2>> const& candidatePolylines,
    std::vector<Operon::Scalar> const& s2,
    std::vector<Operon::Scalar> const& nu,
    bool ptClamp) -> BatchedPointSegReverseResult;

auto evaluatePopulationPointSegGpu(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::vector<Vec2>> const& candidatePolylines,
    bool unsignedMag = false) -> BatchedPointSegResult;

auto evaluatePopulationPointSegGpuResident(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::vector<Vec2>> const& candidatePolylines) -> BatchedPointSegResult;

#if defined(__APPLE__) && defined(__OBJC__)
inline constexpr bool kPointSegGpuCompiledIn = true;
#else
inline constexpr bool kPointSegGpuCompiledIn = false;
#endif

struct PointSegGpuCounters {
    std::atomic<long long> gpuBatches{0};
    std::atomic<long long> gpuCandidates{0};
    std::atomic<long long> cpuInGpuBatch{0};
    std::atomic<long long> fallbackBatches{0};
};
inline auto pointSegGpuCounters() -> PointSegGpuCounters&
{
    static PointSegGpuCounters c;
    return c;
}

struct PointSegSchedDiag {
    std::atomic<long long> firstPassNs{0};
};
inline auto pointSegSchedDiag() -> PointSegSchedDiag&
{
    static PointSegSchedDiag d;
    return d;
}
inline auto pointSegSchedDiagEnabled() -> bool
{
    static bool const v = std::getenv("SKETCH_GPU_UMA_TIMING") != nullptr;
    return v;
}

template <class Pt>
inline auto meanSegLenOfSpan(Pt const* poly, std::size_t n) -> float
{
    constexpr float kEps = std::numeric_limits<float>::epsilon();
    float const lenSqEps = kEps * 100.0f;
    double sum = 0.0; std::size_t cnt = 0;
    for (std::size_t j = 0; j + 1 < n; ++j) {
        float dx = static_cast<float>(poly[j + 1].x - poly[j].x);
        float dy = static_cast<float>(poly[j + 1].y - poly[j].y);
        float lenSq = dx * dx + dy * dy;
        if (lenSq > lenSqEps) { sum += std::sqrt(static_cast<double>(lenSq)); ++cnt; }
    }
    return cnt > 0 ? static_cast<float>(sum / static_cast<double>(cnt)) : 0.0f;
}

#if defined(__APPLE__) && defined(__OBJC__)

namespace detail_gpu_pointseg {

#include "sketch/pointseg_metal_source.hpp"
#define kShader kShaderGenerated

struct Params { uint32_t P, M, Nc; };
struct ReverseParams { uint32_t P, Nc, Ns, K, ptClamp; };
struct Vec2f { float x, y; };
struct Vec2h { _Float16 x, y; };

inline auto meanSegLenOf(std::vector<Vec2> const& poly) -> float
{
    constexpr float kEps = std::numeric_limits<float>::epsilon();
    float const lenSqEps = kEps * 100.0f;
    double sum = 0.0; std::size_t n = 0;
    for (std::size_t j = 0; j + 1 < poly.size(); ++j) {
        float dx = static_cast<float>(poly[j + 1].x - poly[j].x);
        float dy = static_cast<float>(poly[j + 1].y - poly[j].y);
        float lenSq = dx * dx + dy * dy;
        if (lenSq > lenSqEps) { sum += std::sqrt(static_cast<double>(lenSq)); ++n; }
    }
    return n > 0 ? static_cast<float>(sum / static_cast<double>(n)) : 0.0f;
}

struct Gpu {
    id<MTLDevice> dev{nil};
    id<MTLCommandQueue> queue{nil};
    id<MTLComputePipelineState> fusedTg{nil};
    id<MTLComputePipelineState> fusedTgUnsigned{nil};
    id<MTLComputePipelineState> fusedTgFp16Storage{nil};
    id<MTLComputePipelineState> fusedTgFp16StorageUnsigned{nil};
    id<MTLComputePipelineState> fusedTgFp16Arith{nil};
    id<MTLComputePipelineState> fusedTgFp16ArithUnsigned{nil};
    id<MTLComputePipelineState> fusedTgFp16Local{nil};
    id<MTLComputePipelineState> fusedTgFp16LocalUnsigned{nil};
    id<MTLComputePipelineState> reverseHalfTg{nil};
    bool built{false};
    bool buildAttempted{false};
};

inline auto gpuSingleton() -> Gpu&
{
    static Gpu g = [] {
        Gpu out;
        out.buildAttempted = true;
        out.dev = MTLCreateSystemDefaultDevice();
        if (!out.dev) { return out; }
        id<MTLLibrary> lib = nil;
#ifdef SKETCH_POINTSEG_METALLIB_PATH
        {
            char const* mlEnv = std::getenv("SKETCH_METALLIB");
            if (mlEnv != nullptr && std::atoi(mlEnv) != 0) {
                NSString* path = [NSString stringWithUTF8String:SKETCH_POINTSEG_METALLIB_PATH];
                if ([[NSFileManager defaultManager] fileExistsAtPath:path]) {
                    NSURL* url = [NSURL fileURLWithPath:path];
                    NSError* mlErr = nil;
                    lib = [out.dev newLibraryWithURL:url error:&mlErr];
                }
            }
        }
#endif
        if (!lib) {
            NSError* err = nil;
            MTLCompileOptions* opts = [MTLCompileOptions new];
            if (@available(macOS 15.0, iOS 18.0, *)) { opts.mathMode = MTLMathModeSafe; }
            else { opts.fastMathEnabled = NO; }
            lib = [out.dev newLibraryWithSource:[NSString stringWithUTF8String:kShader] options:opts error:&err];
            if (!lib && std::getenv("SKETCH_METAL_COMPILE_DEBUG") != nullptr) {
                std::fprintf(stderr, "[metal-compile-error] %s\n",
                              err ? [[err localizedDescription] UTF8String] : "(no NSError)");
            }
        }
        if (!lib) { return out; }
        id<MTLFunction> fn = [lib newFunctionWithName:@"fusedThreadgroup"];
        if (!fn) { return out; }
        NSError* perr = nil;
        out.fusedTg = [out.dev newComputePipelineStateWithFunction:fn error:&perr];
        if (!out.fusedTg) { return out; }
        if (id<MTLFunction> fnU = [lib newFunctionWithName:@"fusedThreadgroupUnsigned"]) {
            NSError* uerr = nil;
            out.fusedTgUnsigned = [out.dev newComputePipelineStateWithFunction:fnU error:&uerr];
        }
        if (id<MTLFunction> fn1 = [lib newFunctionWithName:@"fusedThreadgroupFp16Storage"]) {
            NSError* e1 = nil;
            out.fusedTgFp16Storage = [out.dev newComputePipelineStateWithFunction:fn1 error:&e1];
        }
        if (id<MTLFunction> fn1u = [lib newFunctionWithName:@"fusedThreadgroupFp16StorageUnsigned"]) {
            NSError* e1u = nil;
            out.fusedTgFp16StorageUnsigned = [out.dev newComputePipelineStateWithFunction:fn1u error:&e1u];
        }
        if (id<MTLFunction> fn2 = [lib newFunctionWithName:@"fusedThreadgroupFp16Arith"]) {
            NSError* e2 = nil;
            out.fusedTgFp16Arith = [out.dev newComputePipelineStateWithFunction:fn2 error:&e2];
        }
        if (id<MTLFunction> fn2u = [lib newFunctionWithName:@"fusedThreadgroupFp16ArithUnsigned"]) {
            NSError* e2u = nil;
            out.fusedTgFp16ArithUnsigned = [out.dev newComputePipelineStateWithFunction:fn2u error:&e2u];
        }
        if (id<MTLFunction> fn3 = [lib newFunctionWithName:@"fusedThreadgroupFp16Local"]) {
            NSError* e3 = nil;
            out.fusedTgFp16Local = [out.dev newComputePipelineStateWithFunction:fn3 error:&e3];
        }
        if (id<MTLFunction> fn3u = [lib newFunctionWithName:@"fusedThreadgroupFp16LocalUnsigned"]) {
            NSError* e3u = nil;
            out.fusedTgFp16LocalUnsigned = [out.dev newComputePipelineStateWithFunction:fn3u error:&e3u];
        }
        if (id<MTLFunction> fnRev = [lib newFunctionWithName:@"reverseHalfThreadgroup"]) {
            NSError* revErr = nil;
            out.reverseHalfTg = [out.dev newComputePipelineStateWithFunction:fnRev error:&revErr];
        }
        out.queue = [out.dev newCommandQueue];
        out.built = (out.queue != nil);
        return out;
    }();
    return g;
}

inline auto pipelineFor(Gpu const& g, bool unsignedMag) -> id<MTLComputePipelineState>
{
    return unsignedMag ? g.fusedTgUnsigned : g.fusedTg;
}

inline auto pipelineFor(Gpu const& g, bool unsignedMag, int fp16Mode) -> id<MTLComputePipelineState>
{
    if (fp16Mode == 1) {
        id<MTLComputePipelineState> const p = unsignedMag ? g.fusedTgFp16StorageUnsigned : g.fusedTgFp16Storage;
        if (p != nil) { return p; }
    } else if (fp16Mode == 2) {
        id<MTLComputePipelineState> const p = unsignedMag ? g.fusedTgFp16ArithUnsigned : g.fusedTgFp16Arith;
        if (p != nil) { return p; }
    } else if (fp16Mode == 3) {
        id<MTLComputePipelineState> const p = unsignedMag ? g.fusedTgFp16LocalUnsigned : g.fusedTgFp16Local;
        if (p != nil) { return p; }
    }
    return pipelineFor(g, unsignedMag);
}

inline auto fp16PipelineAvailable(Gpu const& g, bool unsignedMag, int fp16Mode) -> bool
{
    if (fp16Mode == 1) { return (unsignedMag ? g.fusedTgFp16StorageUnsigned : g.fusedTgFp16Storage) != nil; }
    if (fp16Mode == 2) { return (unsignedMag ? g.fusedTgFp16ArithUnsigned : g.fusedTgFp16Arith) != nil; }
    if (fp16Mode == 3) { return (unsignedMag ? g.fusedTgFp16LocalUnsigned : g.fusedTgFp16Local) != nil; }
    return true;
}

inline auto gpuTimingEnabled() -> bool
{
    static bool const v = std::getenv("SKETCH_GPU_UMA_TIMING") != nullptr;
    return v;
}

struct GpuTimingStats {
    double packMs{0.0};
    double memcpyMs{0.0};
    double allocMs{0.0};
    double encodeMs{0.0};
    double waitMs{0.0};
    double queueMs{0.0};
    double gpuExecMs{0.0};
    double wakeMs{0.0};
    double readbackMs{0.0};
    long calls{0};
    long streamChunks{0};
    std::atomic<long long> streamWorkerWaitNs{0};
    std::atomic<long> streamWorkerWaits{0};
    GpuTimingStats() = default;
    GpuTimingStats(GpuTimingStats const&) = delete;
    auto operator=(GpuTimingStats const&) -> GpuTimingStats& = delete;
};

inline auto gpuTimingStats() -> GpuTimingStats&
{
    static GpuTimingStats s;
    return s;
}

inline auto hostSecondsNow() -> double
{
    static double const scale = [] {
        mach_timebase_info_data_t tb{};
        mach_timebase_info(&tb);
        return static_cast<double>(tb.numer) / static_cast<double>(tb.denom) * 1e-9;
    }();
    return static_cast<double>(mach_absolute_time()) * scale;
}

inline void dumpGpuTimingStats(char const* label)
{
    if (!gpuTimingEnabled()) { return; }
    auto const& s = gpuTimingStats();
    if (s.calls == 0) { std::fprintf(stderr, "[gpu-timing] %s: no calls recorded\n", label); return; }
    std::fprintf(stderr,
        "[gpu-timing] %s: wait split per-call avg: queue(commit->GPUStart)=%.4f gpuExec(GPUStart->GPUEnd)=%.4f "
        "wake(GPUEnd->CPU)=%.4f readback=%.4f (ms) | stream chunks=%ld workerWaits=%ld workerWaitMs(tot)=%.3f\n",
        label, s.queueMs / s.calls, s.gpuExecMs / s.calls, s.wakeMs / s.calls, s.readbackMs / s.calls,
        s.streamChunks, s.streamWorkerWaits.load(), static_cast<double>(s.streamWorkerWaitNs.load()) * 1e-6);
    std::fprintf(stderr,
        "[gpu-timing] %s: calls=%ld pack=%.4fms(tot) memcpy=%.4fms(tot) alloc=%.4fms(tot) "
        "encode=%.4fms(tot) wait=%.4fms(tot) | per-call avg: pack=%.4f memcpy=%.4f alloc=%.4f "
        "encode=%.4f wait=%.4f (ms)\n",
        label, s.calls, s.packMs, s.memcpyMs, s.allocMs, s.encodeMs, s.waitMs,
        s.packMs / s.calls, s.memcpyMs / s.calls, s.allocMs / s.calls,
        s.encodeMs / s.calls, s.waitMs / s.calls);
}

struct ReverseHalfTimingStats {
    double totalMs{0.0};
    long calls{0};
};
inline auto reverseHalfTimingStats() -> ReverseHalfTimingStats&
{
    static ReverseHalfTimingStats s;
    return s;
}
inline void dumpReverseHalfTimingStats(char const* label)
{
    if (!gpuTimingEnabled()) { return; }
    auto const& s = reverseHalfTimingStats();
    if (s.calls == 0) { std::fprintf(stderr, "[gpu-timing] %s: no calls recorded\n", label); return; }
    std::fprintf(stderr, "[gpu-timing] %s: calls=%ld total=%.4fms avg=%.4fms/call\n",
                 label, s.calls, s.totalMs, s.totalMs / s.calls);
}

inline void resetGpuTimingStats()
{
    auto& s = gpuTimingStats();
    s.packMs = s.memcpyMs = s.allocMs = s.encodeMs = s.waitMs = 0.0;
    s.queueMs = s.gpuExecMs = s.wakeMs = s.readbackMs = 0.0;
    s.calls = 0;
    s.streamChunks = 0;
    s.streamWorkerWaitNs.store(0);
    s.streamWorkerWaits.store(0);
}

inline void accumulateCbSplit(id<MTLCommandBuffer> cb, double commitHostS, double wakeHostS)
{
    auto& s = gpuTimingStats();
    double const gs = cb.GPUStartTime, ge = cb.GPUEndTime;
    if (gs > 0.0 && ge >= gs) {
        s.queueMs += (gs - commitHostS) * 1e3;
        s.gpuExecMs += (ge - gs) * 1e3;
        if (wakeHostS > 0.0) { s.wakeMs += (wakeHostS - ge) * 1e3; }
    }
}

}

inline auto evaluatePopulationPointSegGpu(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::vector<Vec2>> const& candidatePolylines,
    bool unsignedMag) -> BatchedPointSegResult
{
    using namespace detail_gpu_pointseg;
    BatchedPointSegResult out;
  @autoreleasepool {
    auto const P = candidatePolylines.size();
    auto const N = strokePoints.size();
    if (P == 0 || N < 2) { return out; }
    auto const Nc = candidatePolylines[0].size();
    if (Nc < 2) { return out; }
    for (auto const& poly : candidatePolylines) {
        if (poly.size() != Nc) { return out; }
    }

    auto& g = gpuSingleton();
    if (!g.built) { return out; }
    bool const timing = gpuTimingEnabled();
    using Clock = std::chrono::steady_clock;

    int const fp16Mode = pointSegFp16ModeFromEnv();
    bool const fp16Storage = fp16Mode == 1 && fp16PipelineAvailable(g, unsignedMag, 1);

    Clock::time_point t0; if (timing) { t0 = Clock::now(); }
    std::vector<Vec2f> dataF(N);
    for (std::size_t i = 0; i < N; ++i) { dataF[i] = {static_cast<float>(strokePoints[i].x), static_cast<float>(strokePoints[i].y)}; }

    std::vector<Vec2f> candFlat(P * Nc);
    std::vector<float> meanSegLen(P);
    for (std::size_t p = 0; p < P; ++p) {
        auto const& poly = candidatePolylines[p];
        for (std::size_t j = 0; j < Nc; ++j) {
            candFlat[p * Nc + j] = {static_cast<float>(poly[j].x), static_cast<float>(poly[j].y)};
        }
        meanSegLen[p] = meanSegLenOf(poly);
    }
    std::vector<Vec2h> dataH, candH;
    if (fp16Storage) {
        dataH.resize(N);
        for (std::size_t i = 0; i < N; ++i) { dataH[i] = {static_cast<_Float16>(dataF[i].x), static_cast<_Float16>(dataF[i].y)}; }
        candH.resize(P * Nc);
        for (std::size_t k = 0; k < candFlat.size(); ++k) { candH[k] = {static_cast<_Float16>(candFlat[k].x), static_cast<_Float16>(candFlat[k].y)}; }
    }
    Clock::time_point t1; if (timing) { t1 = Clock::now(); }

    Params prm{static_cast<uint32_t>(P), static_cast<uint32_t>(N), static_cast<uint32_t>(Nc)};

    id<MTLBuffer> bufCand = fp16Storage
        ? [g.dev newBufferWithBytes:candH.data() length:candH.size() * sizeof(Vec2h) options:MTLResourceStorageModeShared]
        : [g.dev newBufferWithBytes:candFlat.data() length:candFlat.size() * sizeof(Vec2f) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufData = fp16Storage
        ? [g.dev newBufferWithBytes:dataH.data() length:dataH.size() * sizeof(Vec2h) options:MTLResourceStorageModeShared]
        : [g.dev newBufferWithBytes:dataF.data() length:dataF.size() * sizeof(Vec2f) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufMeanSegLen = [g.dev newBufferWithBytes:meanSegLen.data() length:meanSegLen.size() * sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufResidual = [g.dev newBufferWithLength:static_cast<std::size_t>(P) * N * sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufFitness = [g.dev newBufferWithLength:static_cast<std::size_t>(P) * sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufTau = [g.dev newBufferWithLength:static_cast<std::size_t>(P) * N * sizeof(float) options:MTLResourceStorageModeShared];
    if (!bufCand || !bufData || !bufMeanSegLen || !bufResidual || !bufFitness || !bufTau) { return out; }
    Clock::time_point t2; if (timing) { t2 = Clock::now(); }

    id<MTLComputePipelineState> const pso = fp16Storage ? pipelineFor(g, unsignedMag, 1)
        : pipelineFor(g, unsignedMag, (fp16Mode == 2 || fp16Mode == 3) ? fp16Mode : 0);
    if (pso == nil) { return out; }
    NSUInteger tgThreadsWanted = ((static_cast<NSUInteger>(N) + 31) / 32) * 32;
    NSUInteger tgThreads = std::min(tgThreadsWanted, pso.maxTotalThreadsPerThreadgroup);
    if (tgThreads < static_cast<NSUInteger>(N)) { return out; }

    id<MTLCommandBuffer> cb = [g.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:bufCand offset:0 atIndex:0];
    [enc setBuffer:bufData offset:0 atIndex:1];
    [enc setBytes:&prm length:sizeof prm atIndex:2];
    [enc setBuffer:bufMeanSegLen offset:0 atIndex:3];
    [enc setBuffer:bufResidual offset:0 atIndex:4];
    [enc setBuffer:bufFitness offset:0 atIndex:5];
    [enc setBuffer:bufTau offset:0 atIndex:6];
    [enc setThreadgroupMemoryLength:32 * sizeof(float) atIndex:0];
    [enc dispatchThreadgroups:MTLSizeMake(P, 1, 1) threadsPerThreadgroup:MTLSizeMake(tgThreads, 1, 1)];
    [enc endEncoding];
    [cb commit];
    Clock::time_point t3; if (timing) { t3 = Clock::now(); }
    [cb waitUntilCompleted];
    Clock::time_point t4; if (timing) { t4 = Clock::now(); }
    if (cb.status == MTLCommandBufferStatusError) { return out; }

    out.fitness.assign(static_cast<float*>(bufFitness.contents), static_cast<float*>(bufFitness.contents) + P);
    out.residualFlat.assign(static_cast<float*>(bufResidual.contents), static_cast<float*>(bufResidual.contents) + static_cast<std::size_t>(P) * N);
    out.tauFlat.assign(static_cast<float*>(bufTau.contents), static_cast<float*>(bufTau.contents) + static_cast<std::size_t>(P) * N);
    out.ok = true;
    out.fp16ModeUsed = fp16Storage ? 1
        : (((fp16Mode == 2 || fp16Mode == 3) && fp16PipelineAvailable(g, unsignedMag, fp16Mode)) ? fp16Mode : 0);
    if (timing) {
        auto& s = gpuTimingStats();
        s.packMs += std::chrono::duration<double, std::milli>(t1 - t0).count();
        s.allocMs += std::chrono::duration<double, std::milli>(t2 - t1).count();
        s.encodeMs += std::chrono::duration<double, std::milli>(t3 - t2).count();
        s.waitMs += std::chrono::duration<double, std::milli>(t4 - t3).count();
        ++s.calls;
    }
    return out;
  }
}

inline auto evaluatePopulationPointSegReverseGpu(
    std::vector<Vec2> const& strokePointsDrawn,
    std::vector<std::vector<Vec2>> const& candidatePolylines,
    std::vector<Operon::Scalar> const& s2,
    std::vector<Operon::Scalar> const& nu,
    bool ptClamp) -> BatchedPointSegReverseResult
{
    using namespace detail_gpu_pointseg;
    BatchedPointSegReverseResult out;
    bool const timing = gpuTimingEnabled();
    double const t0 = timing ? hostSecondsNow() : 0.0;
  @autoreleasepool {
    auto const P = candidatePolylines.size();
    auto const Ns = strokePointsDrawn.size();
    if (P == 0 || Ns < 2) { return out; }
    auto const Nc = candidatePolylines[0].size();
    if (Nc < 2) { return out; }
    for (auto const& poly : candidatePolylines) {
        if (poly.size() != Nc) { return out; }
    }
    auto const K = s2.size();
    if (K != nu.size()) { return out; }

    auto& g = gpuSingleton();
    if (!g.built || g.reverseHalfTg == nil) { return out; }

    std::vector<Vec2f> drawnF(Ns);
    for (std::size_t i = 0; i < Ns; ++i) {
        drawnF[i] = {static_cast<float>(strokePointsDrawn[i].x), static_cast<float>(strokePointsDrawn[i].y)};
    }
    std::vector<Vec2f> candFlat(P * Nc);
    for (std::size_t p = 0; p < P; ++p) {
        auto const& poly = candidatePolylines[p];
        for (std::size_t j = 0; j < Nc; ++j) {
            candFlat[p * Nc + j] = {static_cast<float>(poly[j].x), static_cast<float>(poly[j].y)};
        }
    }
    std::vector<float> s2F(K > 0 ? K : 1, 1.0f), nuF(K > 0 ? K : 1, 0.0f);
    for (std::size_t k = 0; k < K; ++k) { s2F[k] = static_cast<float>(s2[k]); nuF[k] = static_cast<float>(nu[k]); }

    ReverseParams prm{static_cast<uint32_t>(P), static_cast<uint32_t>(Nc), static_cast<uint32_t>(Ns),
                       static_cast<uint32_t>(K), ptClamp ? 1u : 0u};

    id<MTLBuffer> bufCandVerts = [g.dev newBufferWithBytes:candFlat.data() length:candFlat.size() * sizeof(Vec2f) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufDrawn = [g.dev newBufferWithBytes:drawnF.data() length:drawnF.size() * sizeof(Vec2f) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufS2 = [g.dev newBufferWithBytes:s2F.data() length:s2F.size() * sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufNu = [g.dev newBufferWithBytes:nuF.data() length:nuF.size() * sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufReverseFitness = [g.dev newBufferWithLength:static_cast<std::size_t>(P) * sizeof(float) options:MTLResourceStorageModeShared];
    if (!bufCandVerts || !bufDrawn || !bufS2 || !bufNu || !bufReverseFitness) { return out; }

    id<MTLComputePipelineState> const pso = g.reverseHalfTg;
    NSUInteger tgThreadsWanted = ((static_cast<NSUInteger>(Nc) + 31) / 32) * 32;
    NSUInteger tgThreads = std::min(tgThreadsWanted, pso.maxTotalThreadsPerThreadgroup);
    if (tgThreads < static_cast<NSUInteger>(Nc)) { return out; }

    id<MTLCommandBuffer> cb = [g.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:bufCandVerts offset:0 atIndex:0];
    [enc setBuffer:bufDrawn offset:0 atIndex:1];
    [enc setBytes:&prm length:sizeof prm atIndex:2];
    [enc setBuffer:bufS2 offset:0 atIndex:3];
    [enc setBuffer:bufNu offset:0 atIndex:4];
    [enc setBuffer:bufReverseFitness offset:0 atIndex:5];
    [enc setThreadgroupMemoryLength:32 * sizeof(float) atIndex:0];
    [enc dispatchThreadgroups:MTLSizeMake(P, 1, 1) threadsPerThreadgroup:MTLSizeMake(tgThreads, 1, 1)];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];
    if (cb.status == MTLCommandBufferStatusError) { return out; }

    out.fitness.assign(static_cast<float*>(bufReverseFitness.contents), static_cast<float*>(bufReverseFitness.contents) + P);
    out.ok = true;
    if (timing) {
        auto& s = reverseHalfTimingStats();
        s.totalMs += (hostSecondsNow() - t0) * 1e3;
        ++s.calls;
    }
    return out;
  }
}

namespace detail_gpu_pointseg {

struct ResidentBuffers {
    id<MTLBuffer> bufCand{nil};
    id<MTLBuffer> bufData{nil};
    id<MTLBuffer> bufMeanSegLen{nil};
    id<MTLBuffer> bufResidual{nil};
    id<MTLBuffer> bufFitness{nil};
    id<MTLBuffer> bufTau{nil};
    std::size_t capCand{0};
    std::size_t capData{0};
    std::size_t capMeanSegLen{0};
    std::size_t capFitness{0};
    std::size_t capResidual{0};
    std::size_t capTau{0};
};


inline void ensureCapacity(id<MTLDevice> dev, __strong id<MTLBuffer>& buf, std::size_t& cap,
                            std::size_t neededElems, std::size_t elemSize)
{
    if (neededElems <= cap && buf != nil) { return; }
    std::size_t const newCap = std::max(neededElems, cap == 0 ? neededElems : cap * 2);
    buf = [dev newBufferWithLength:newCap * elemSize options:MTLResourceStorageModeShared];
    cap = newCap;
}

}

class PointSegGpuResident {
public:
    detail_gpu_pointseg::ResidentBuffers rb;
    bool unsignedMag{false};
};

inline auto evaluatePopulationPointSegGpuResident(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::vector<Vec2>> const& candidatePolylines,
    PointSegGpuResident& owner) -> BatchedPointSegResult
{
    using namespace detail_gpu_pointseg;
    BatchedPointSegResult out;
    bool const unsignedMag = owner.unsignedMag;
  @autoreleasepool {
    auto const P = candidatePolylines.size();
    auto const N = strokePoints.size();
    if (P == 0 || N < 2) { return out; }
    auto const Nc = candidatePolylines[0].size();
    if (Nc < 2) { return out; }
    for (auto const& poly : candidatePolylines) {
        if (poly.size() != Nc) { return out; }
    }

    auto& g = gpuSingleton();
    if (!g.built) { return out; }
    auto& rb = owner.rb;
    bool const timing = gpuTimingEnabled();
    using Clock = std::chrono::steady_clock;

    Clock::time_point t0; if (timing) { t0 = Clock::now(); }
    std::vector<Vec2f> dataF(N);
    for (std::size_t i = 0; i < N; ++i) { dataF[i] = {static_cast<float>(strokePoints[i].x), static_cast<float>(strokePoints[i].y)}; }

    std::vector<Vec2f> candFlat(P * Nc);
    std::vector<float> meanSegLen(P);
    for (std::size_t p = 0; p < P; ++p) {
        auto const& poly = candidatePolylines[p];
        for (std::size_t j = 0; j < Nc; ++j) {
            candFlat[p * Nc + j] = {static_cast<float>(poly[j].x), static_cast<float>(poly[j].y)};
        }
        meanSegLen[p] = meanSegLenOf(poly);
    }
    Clock::time_point t1; if (timing) { t1 = Clock::now(); }

    Params prm{static_cast<uint32_t>(P), static_cast<uint32_t>(N), static_cast<uint32_t>(Nc)};

    ensureCapacity(g.dev, rb.bufCand, rb.capCand, P * Nc, sizeof(Vec2f));
    ensureCapacity(g.dev, rb.bufData, rb.capData, N, sizeof(Vec2f));
    ensureCapacity(g.dev, rb.bufMeanSegLen, rb.capMeanSegLen, P, sizeof(float));
    ensureCapacity(g.dev, rb.bufFitness, rb.capFitness, P, sizeof(float));
    ensureCapacity(g.dev, rb.bufResidual, rb.capResidual, static_cast<std::size_t>(P) * N, sizeof(float));
    ensureCapacity(g.dev, rb.bufTau, rb.capTau, static_cast<std::size_t>(P) * N, sizeof(float));
    if (!rb.bufCand || !rb.bufData || !rb.bufMeanSegLen || !rb.bufResidual || !rb.bufFitness || !rb.bufTau) { return out; }
    Clock::time_point t1b; if (timing) { t1b = Clock::now(); }

    std::memcpy(rb.bufCand.contents, candFlat.data(), candFlat.size() * sizeof(Vec2f));
    std::memcpy(rb.bufData.contents, dataF.data(), dataF.size() * sizeof(Vec2f));
    std::memcpy(rb.bufMeanSegLen.contents, meanSegLen.data(), meanSegLen.size() * sizeof(float));
    Clock::time_point t2; if (timing) { t2 = Clock::now(); }

    id<MTLComputePipelineState> const pso = pipelineFor(g, unsignedMag, pointSegFp16ModeForFloatBuffers());
    if (pso == nil) { return out; }
    NSUInteger tgThreadsWanted = ((static_cast<NSUInteger>(N) + 31) / 32) * 32;
    NSUInteger tgThreads = std::min(tgThreadsWanted, pso.maxTotalThreadsPerThreadgroup);
    if (tgThreads < static_cast<NSUInteger>(N)) { return out; }

    id<MTLCommandBuffer> cb = [g.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:rb.bufCand offset:0 atIndex:0];
    [enc setBuffer:rb.bufData offset:0 atIndex:1];
    [enc setBytes:&prm length:sizeof prm atIndex:2];
    [enc setBuffer:rb.bufMeanSegLen offset:0 atIndex:3];
    [enc setBuffer:rb.bufResidual offset:0 atIndex:4];
    [enc setBuffer:rb.bufFitness offset:0 atIndex:5];
    [enc setBuffer:rb.bufTau offset:0 atIndex:6];
    [enc setThreadgroupMemoryLength:32 * sizeof(float) atIndex:0];
    [enc dispatchThreadgroups:MTLSizeMake(P, 1, 1) threadsPerThreadgroup:MTLSizeMake(tgThreads, 1, 1)];
    [enc endEncoding];
    [cb commit];
    Clock::time_point t3; if (timing) { t3 = Clock::now(); }
    double const commitHostS = timing ? hostSecondsNow() : 0.0;
    [cb waitUntilCompleted];
    Clock::time_point t4; if (timing) { t4 = Clock::now(); }
    double const wakeHostS = timing ? hostSecondsNow() : 0.0;
    if (cb.status == MTLCommandBufferStatusError) { return out; }

    out.fitness.assign(static_cast<float*>(rb.bufFitness.contents), static_cast<float*>(rb.bufFitness.contents) + P);
    out.residualFlat.assign(static_cast<float*>(rb.bufResidual.contents), static_cast<float*>(rb.bufResidual.contents) + static_cast<std::size_t>(P) * N);
    out.tauFlat.assign(static_cast<float*>(rb.bufTau.contents), static_cast<float*>(rb.bufTau.contents) + static_cast<std::size_t>(P) * N);
    out.ok = true;
    {
        int const wantMode = pointSegFp16ModeForFloatBuffers();
        out.fp16ModeUsed = (wantMode != 0 && fp16PipelineAvailable(g, unsignedMag, wantMode)) ? wantMode : 0;
    }
    if (timing) {
        auto& s = gpuTimingStats();
        s.readbackMs += std::chrono::duration<double, std::milli>(Clock::now() - t4).count();
        accumulateCbSplit(cb, commitHostS, wakeHostS);
        s.packMs += std::chrono::duration<double, std::milli>(t1 - t0).count();
        s.allocMs += std::chrono::duration<double, std::milli>(t1b - t1).count();
        s.memcpyMs += std::chrono::duration<double, std::milli>(t2 - t1b).count();
        s.encodeMs += std::chrono::duration<double, std::milli>(t3 - t2).count();
        s.waitMs += std::chrono::duration<double, std::milli>(t4 - t3).count();
        ++s.calls;
    }
    return out;
  }
}

inline auto evaluatePopulationPointSegGpuResident(
    std::vector<Vec2> const& strokePoints,
    std::vector<std::vector<Vec2>> const& candidatePolylines) -> BatchedPointSegResult
{
    static PointSegGpuResident shared;
    return evaluatePopulationPointSegGpuResident(strokePoints, candidatePolylines, shared);
}

struct PointSegGpuStreamOptions {
    bool unsignedMag{false};
    bool unretainedRefs{false};
    bool spinWait{false};
    bool residencySet{false};
};

class PointSegGpuStream {
public:
    using Vec2f = detail_gpu_pointseg::Vec2f;
    PointSegGpuStream() = default;
    PointSegGpuStream(PointSegGpuStream const&) = delete;
    auto operator=(PointSegGpuStream const&) -> PointSegGpuStream& = delete;
    ~PointSegGpuStream() { drain(); }

    auto begin(Operon::Scalar const* sx, Operon::Scalar const* sy, std::size_t n, std::size_t pMax,
               std::size_t nc, PointSegGpuStreamOptions const& opts) -> bool
    {
        using namespace detail_gpu_pointseg;
        drain();
        ok_ = false;
        opts_ = opts;
        auto& g = gpuSingleton();
        if (!g.built || n < 2 || nc < 2 || pMax == 0) { return false; }
        id<MTLComputePipelineState> const pso = pipelineFor(g, opts.unsignedMag, pointSegFp16ModeForFloatBuffers());
        if (pso == nil) { return false; }
        NSUInteger const w = pso.threadExecutionWidth > 0 ? pso.threadExecutionWidth : 32;
        NSUInteger const tgWanted = ((static_cast<NSUInteger>(n) + w - 1) / w) * w;
        tgThreads_ = std::min(tgWanted, pso.maxTotalThreadsPerThreadgroup);
        if (tgThreads_ < static_cast<NSUInteger>(n)) { return false; }
        bool grew = false;
        MTLResourceOptions const bufOpts = opts_.unretainedRefs
            ? (MTLResourceStorageModeShared | MTLResourceHazardTrackingModeUntracked)
            : MTLResourceStorageModeShared;
        if (opts_.unretainedRefs != untracked_) { capCand_ = capData_ = capMean_ = capFit_ = capRes_ = capTau_ = 0; }
        untracked_ = opts_.unretainedRefs;
        auto grow = [&](__strong id<MTLBuffer>& buf, std::size_t& cap, std::size_t need, std::size_t elem) {
            if (need <= cap && buf != nil) { return; }
            std::size_t const nc2 = std::max(need, cap == 0 ? need : cap * 2);
            buf = [g.dev newBufferWithLength:nc2 * elem options:bufOpts];
            cap = nc2;
            grew = true;
        };
        grow(bufCand_, capCand_, pMax * nc, sizeof(Vec2f));
        grow(bufData_, capData_, n, sizeof(Vec2f));
        grow(bufMean_, capMean_, pMax, sizeof(float));
        grow(bufFit_, capFit_, pMax, sizeof(float));
        grow(bufRes_, capRes_, pMax * n, sizeof(float));
        grow(bufTau_, capTau_, pMax * n, sizeof(float));
        if (!bufCand_ || !bufData_ || !bufMean_ || !bufFit_ || !bufRes_ || !bufTau_) { return false; }
        if (opts_.residencySet && (grew || residency_ == nil)) { rebuildResidency(g.dev); }
        auto* d = static_cast<Vec2f*>(bufData_.contents);
        for (std::size_t i = 0; i < n; ++i) { d[i] = {static_cast<float>(sx[i]), static_cast<float>(sy[i])}; }
        n_ = n; nc_ = nc; pMax_ = pMax;
        chunkOf_.assign(pMax, 0u);
        committed_ = 0;
        ok_ = true;
        if (gpuTimingEnabled()) { ++gpuTimingStats().calls; }
        return true;
    }

    auto candSlot(std::size_t p) -> Vec2f* { return static_cast<Vec2f*>(bufCand_.contents) + p * nc_; }
    void setMeanSegLen(std::size_t p, float v) { static_cast<float*>(bufMean_.contents)[p] = v; }

    auto commitRange(std::size_t p0, std::size_t p1) -> bool
    {
        using namespace detail_gpu_pointseg;
        if (!ok_ || p1 <= p0 || p1 > pMax_) { return false; }
      @autoreleasepool {
        bool const timing = gpuTimingEnabled();
        auto const t0 = timing ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        auto& g = gpuSingleton();
        id<MTLCommandBuffer> cb = opts_.unretainedRefs ? [g.queue commandBufferWithUnretainedReferences]
                                                       : [g.queue commandBuffer];
        if (cb == nil) { ok_ = false; return false; }
        if (opts_.residencySet && residency_ != nil) {
            if (@available(macOS 15.0, iOS 18.0, *)) { [cb useResidencySet:(id<MTLResidencySet>)residency_]; }
        }
        Params prm{static_cast<uint32_t>(p1 - p0), static_cast<uint32_t>(n_), static_cast<uint32_t>(nc_)};
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:pipelineFor(g, opts_.unsignedMag, pointSegFp16ModeForFloatBuffers())];
        [enc setBuffer:bufCand_ offset:p0 * nc_ * sizeof(Vec2f) atIndex:0];
        [enc setBuffer:bufData_ offset:0 atIndex:1];
        [enc setBytes:&prm length:sizeof prm atIndex:2];
        [enc setBuffer:bufMean_ offset:p0 * sizeof(float) atIndex:3];
        [enc setBuffer:bufRes_ offset:p0 * n_ * sizeof(float) atIndex:4];
        [enc setBuffer:bufFit_ offset:p0 * sizeof(float) atIndex:5];
        [enc setBuffer:bufTau_ offset:p0 * n_ * sizeof(float) atIndex:6];
        [enc setThreadgroupMemoryLength:32 * sizeof(float) atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(p1 - p0, 1, 1) threadsPerThreadgroup:MTLSizeMake(tgThreads_, 1, 1)];
        [enc endEncoding];
        Chunk c;
        c.p0 = p0; c.p1 = p1; c.cb = cb;
        c.state = std::make_shared<std::atomic<int>>(0);
        auto st = c.state;
        [cb addCompletedHandler:^(id<MTLCommandBuffer> b) {
            st->store(b.status == MTLCommandBufferStatusCompleted ? 1 : 2, std::memory_order_release);
        }];
        auto const idx = static_cast<uint32_t>(chunks_.size());
        for (std::size_t p = p0; p < p1; ++p) { chunkOf_[p] = idx; }
        [cb commit];
        c.commitHostS = timing ? hostSecondsNow() : 0.0;
        chunks_.push_back(std::move(c));
        committed_ = std::max(committed_, p1);
        if (timing) {
            auto& s = gpuTimingStats();
            s.encodeMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            ++s.streamChunks;
        }
        return true;
      }
    }

    [[nodiscard]] auto committedUpTo() const -> std::size_t { return committed_; }
    [[nodiscard]] auto busy() const -> bool
    {
        return !chunks_.empty() && chunks_.back().state->load(std::memory_order_acquire) == 0;
    }

    auto waitAll() -> bool
    {
        bool const timing = detail_gpu_pointseg::gpuTimingEnabled();
        auto const t0 = timing ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        bool all = true;
        for (std::size_t k = 0; k < chunks_.size(); ++k) { all = waitChunk(k) && all; }
        if (timing) {
            auto& s = detail_gpu_pointseg::gpuTimingStats();
            s.waitMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            if (!chunks_.empty()) { s.wakeMs += (detail_gpu_pointseg::hostSecondsNow() - lastGpuEnd()) * 1e3; }
        }
        return all;
    }

    auto waitFor(std::size_t p) const -> bool
    {
        if (!ok_ || p >= committed_) { return false; }
        auto const k = chunkOf_[p];
        if (chunks_[k].state->load(std::memory_order_acquire) == 1) { return true; }
        bool const timing = detail_gpu_pointseg::gpuTimingEnabled();
        auto const t0 = timing ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        bool const r = waitChunk(k);
        if (timing) {
            auto& s = detail_gpu_pointseg::gpuTimingStats();
            s.streamWorkerWaitNs.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t0).count(), std::memory_order_relaxed);
            s.streamWorkerWaits.fetch_add(1, std::memory_order_relaxed);
        }
        return r;
    }

    struct Split { double gpuStart{0.0}, gpuEnd{0.0}; };
    [[nodiscard]] auto lastSplit() const -> Split
    {
        Split s;
        for (auto const& c : chunks_) {
            double const a = c.cb.GPUStartTime, b = c.cb.GPUEndTime;
            if (s.gpuStart == 0.0 || a < s.gpuStart) { s.gpuStart = a; }
            s.gpuEnd = std::max(s.gpuEnd, b);
        }
        return s;
    }

    [[nodiscard]] auto fitness() const -> float const* { return static_cast<float const*>(bufFit_.contents); }
    [[nodiscard]] auto residual() const -> float const* { return static_cast<float const*>(bufRes_.contents); }
    [[nodiscard]] auto tau() const -> float const* { return static_cast<float const*>(bufTau_.contents); }

private:
    using Params = detail_gpu_pointseg::Params;
    struct Chunk {
        std::size_t p0{0}, p1{0};
        id<MTLCommandBuffer> cb{nil};
        std::shared_ptr<std::atomic<int>> state;
        double commitHostS{0.0};
    };

    auto waitChunk(std::size_t k) const -> bool
    {
        auto const& c = chunks_[k];
        int st = c.state->load(std::memory_order_acquire);
        if (st != 0) { return st == 1; }
        if (opts_.spinWait) {
            while ((st = c.state->load(std::memory_order_acquire)) == 0) {
#if defined(__aarch64__)
                __asm__ __volatile__("yield");
#endif
            }
            return st == 1;
        }
        [c.cb waitUntilCompleted];
        return c.cb.status == MTLCommandBufferStatusCompleted;
    }

    [[nodiscard]] auto lastGpuEnd() const -> double
    {
        double e = 0.0;
        for (auto const& c : chunks_) { e = std::max(e, static_cast<double>(c.cb.GPUEndTime)); }
        return e;
    }

    void drain()
    {
        if (chunks_.empty()) { return; }
        for (std::size_t k = 0; k < chunks_.size(); ++k) { (void)waitChunk(k); }
        for (auto const& c : chunks_) { [c.cb waitUntilCompleted]; }
        if (detail_gpu_pointseg::gpuTimingEnabled()) {
            for (auto const& c : chunks_) { detail_gpu_pointseg::accumulateCbSplit(c.cb, c.commitHostS, 0.0); }
        }
        chunks_.clear();
    }

    void rebuildResidency(id<MTLDevice> dev)
    {
        if (@available(macOS 15.0, iOS 18.0, *)) {
            if (residency_ == nil) {
                MTLResidencySetDescriptor* d = [MTLResidencySetDescriptor new];
                d.initialCapacity = 6;
                NSError* err = nil;
                residency_ = [dev newResidencySetWithDescriptor:d error:&err];
                if (residency_ == nil) { return; }
            }
            id<MTLResidencySet> rs = (id<MTLResidencySet>)residency_;
            [rs removeAllAllocations];
            [rs addAllocation:bufCand_]; [rs addAllocation:bufData_]; [rs addAllocation:bufMean_];
            [rs addAllocation:bufFit_]; [rs addAllocation:bufRes_]; [rs addAllocation:bufTau_];
            [rs commit];
            [rs requestResidency];
        }
    }

    PointSegGpuStreamOptions opts_{};
    id<MTLBuffer> bufCand_{nil}, bufData_{nil}, bufMean_{nil}, bufFit_{nil}, bufRes_{nil}, bufTau_{nil};
    std::size_t capCand_{0}, capData_{0}, capMean_{0}, capFit_{0}, capRes_{0}, capTau_{0};
    id residency_{nil};
    NSUInteger tgThreads_{0};
    std::size_t n_{0}, nc_{0}, pMax_{0}, committed_{0};
    bool ok_{false};
    bool untracked_{false};
    std::vector<Chunk> chunks_;
    std::vector<uint32_t> chunkOf_;
};

inline constexpr bool kPointSegGpuStreamAvailable = true;

#else

struct PointSegGpuStreamOptions {
    bool unsignedMag{false};
    bool unretainedRefs{false};
    bool spinWait{false};
    bool residencySet{false};
};

class PointSegGpuStream {
public:
    struct Vec2f { float x, y; };
    auto begin(Operon::Scalar const*, Operon::Scalar const*, std::size_t, std::size_t, std::size_t,
               PointSegGpuStreamOptions const&) -> bool { return false; }
    auto candSlot(std::size_t) -> Vec2f* { return nullptr; }
    void setMeanSegLen(std::size_t, float) {}
    auto commitRange(std::size_t, std::size_t) -> bool { return false; }
    [[nodiscard]] auto committedUpTo() const -> std::size_t { return 0; }
    [[nodiscard]] auto busy() const -> bool { return false; }
    auto waitAll() -> bool { return false; }
    auto waitFor(std::size_t) const -> bool { return false; }
    [[nodiscard]] auto fitness() const -> float const* { return nullptr; }
    [[nodiscard]] auto residual() const -> float const* { return nullptr; }
    [[nodiscard]] auto tau() const -> float const* { return nullptr; }
};

inline constexpr bool kPointSegGpuStreamAvailable = false;

inline auto evaluatePopulationPointSegGpu(
    std::vector<Vec2> const&,
    std::vector<std::vector<Vec2>> const&,
    bool) -> BatchedPointSegResult
{
    return {};
}

inline auto evaluatePopulationPointSegReverseGpu(
    std::vector<Vec2> const&,
    std::vector<std::vector<Vec2>> const&,
    std::vector<Operon::Scalar> const&,
    std::vector<Operon::Scalar> const&,
    bool) -> BatchedPointSegReverseResult
{
    return {};
}

inline auto evaluatePopulationPointSegGpuResident(
    std::vector<Vec2> const&,
    std::vector<std::vector<Vec2>> const&) -> BatchedPointSegResult
{
    return {};
}

class PointSegGpuResident { public: bool unsignedMag{false}; };

inline auto evaluatePopulationPointSegGpuResident(
    std::vector<Vec2> const&,
    std::vector<std::vector<Vec2>> const&,
    PointSegGpuResident&) -> BatchedPointSegResult
{
    return {};
}

#endif

}
