#include "sketch/bracket_newton_runtime.hpp"
#include "sketch/fp16_hot_path.hpp"
#include "sketch/fp16_gram_runtime.hpp"

#include <algorithm>
#include <cstdint>

#if defined(__APPLE__)
#import <Metal/Metal.h>
#endif

namespace Sketch {

struct MetalRuntimeImpl {
#if defined(__APPLE__)
    id<MTLDevice> device{nil};
    id<MTLCommandQueue> queue{nil};
    MetalFp16GramRuntime* gramRuntime{nullptr};
    id<MTLLibrary> library{nil};
    id<MTLComputePipelineState> pipeDiagonal{nil};
    id<MTLComputePipelineState> pipeRidge2x2{nil};
    id<MTLComputePipelineState> pipeRidge4x4{nil};
    id<MTLComputePipelineState> pipeRidge8x8{nil};
    id<MTLComputePipelineState> pipeRidge8x8Half{nil};
    id<MTLComputePipelineState> pipeFootpoints{nil};
#endif
};

#if defined(__APPLE__)
static NSString* const kDiagonalUpdateShader = @R"MSL(
#include <metal_stdlib>
using namespace metal;
kernel void diagonal_ridge_update(device const float* diagonal [[buffer(0)]],
                                  device const float* rhs [[buffer(1)]],
                                  device float* update [[buffer(2)]],
                                  constant float& ridge [[buffer(3)]],
                                  uint i [[thread_position_in_grid]]) {
    float d = diagonal[i] + ridge;
    update[i] = (isfinite(d) && d > 0.0f && isfinite(rhs[i])) ? rhs[i] / d : 0.0f;
}
kernel void coupled_ridge_2x2(device const float* gram [[buffer(0)]],
                               device const float* rhs [[buffer(1)]],
                               device float* update [[buffer(2)]],
                               constant float& ridge [[buffer(3)]],
                               uint i [[thread_position_in_grid]]) {
    if (i != 0) return;
    float a = gram[0] + ridge;
    float b = gram[1];
    float d = gram[3] + ridge;
    float det = a * d - b * b;
    if (!isfinite(det) || fabs(det) <= 1.0e-12f) {
        update[0] = 0.0f; update[1] = 0.0f; return;
    }
    update[0] = (d * rhs[0] - b * rhs[1]) / det;
    update[1] = (a * rhs[1] - b * rhs[0]) / det;
}
kernel void coupled_ridge_4x4(device const float* gram [[buffer(0)]],
                               device const float* rhs [[buffer(1)]],
                               device float* update [[buffer(2)]],
                               constant float& ridge [[buffer(3)]],
                               uint i [[thread_position_in_grid]]) {
    if (i != 0) return;
    float a[16]; float b[4];
    for (uint r = 0; r < 4; ++r) {
        b[r] = rhs[r];
        for (uint c = 0; c < 4; ++c) { a[r * 4 + c] = gram[r * 4 + c]; }
        a[r * 4 + r] += ridge;
    }
    for (uint col = 0; col < 4; ++col) {
        uint pivot = col;
        for (uint row = col + 1; row < 4; ++row) {
            if (fabs(a[row * 4 + col]) > fabs(a[pivot * 4 + col])) { pivot = row; }
        }
        if (!isfinite(a[pivot * 4 + col]) || fabs(a[pivot * 4 + col]) <= 1.0e-12f) {
            for (uint k = 0; k < 4; ++k) { update[k] = 0.0f; }
            return;
        }
        if (pivot != col) {
            for (uint c = col; c < 4; ++c) {
                float t = a[col * 4 + c]; a[col * 4 + c] = a[pivot * 4 + c]; a[pivot * 4 + c] = t;
            }
            float t = b[col]; b[col] = b[pivot]; b[pivot] = t;
        }
        for (uint row = col + 1; row < 4; ++row) {
            float factor = a[row * 4 + col] / a[col * 4 + col];
            for (uint c = col + 1; c < 4; ++c) { a[row * 4 + c] -= factor * a[col * 4 + c]; }
            b[row] -= factor * b[col];
        }
    }
    for (int r = 3; r >= 0; --r) {
        float value = b[r];
        for (uint c = uint(r + 1); c < 4; ++c) { value -= a[uint(r) * 4 + c] * update[c]; }
        update[uint(r)] = value / a[uint(r) * 4 + uint(r)];
    }
}
kernel void coupled_ridge_8x8(device const float* gram [[buffer(0)]],
                               device const float* rhs [[buffer(1)]],
                               device float* update [[buffer(2)]],
                               constant float& ridge [[buffer(3)]],
                               uint i [[thread_position_in_grid]]) {
    if (i != 0) return;
    float a[64]; float b[8];
    for (uint r = 0; r < 8; ++r) {
        b[r] = rhs[r];
        for (uint c = 0; c < 8; ++c) { a[r * 8 + c] = gram[r * 8 + c]; }
        a[r * 8 + r] += ridge;
    }
    for (uint col = 0; col < 8; ++col) {
        uint pivot = col;
        for (uint row = col + 1; row < 8; ++row) {
            if (fabs(a[row * 8 + col]) > fabs(a[pivot * 8 + col])) { pivot = row; }
        }
        if (!isfinite(a[pivot * 8 + col]) || fabs(a[pivot * 8 + col]) <= 1.0e-12f) {
            for (uint k = 0; k < 8; ++k) { update[k] = 0.0f; }
            return;
        }
        if (pivot != col) {
            for (uint c = col; c < 8; ++c) {
                float t = a[col * 8 + c]; a[col * 8 + c] = a[pivot * 8 + c]; a[pivot * 8 + c] = t;
            }
            float t = b[col]; b[col] = b[pivot]; b[pivot] = t;
        }
        for (uint row = col + 1; row < 8; ++row) {
            float factor = a[row * 8 + col] / a[col * 8 + col];
            for (uint c = col + 1; c < 8; ++c) { a[row * 8 + c] -= factor * a[col * 8 + c]; }
            b[row] -= factor * b[col];
        }
    }
    for (int r = 7; r >= 0; --r) {
        float value = b[r];
        for (uint c = uint(r + 1); c < 8; ++c) { value -= a[uint(r) * 8 + c] * update[c]; }
        update[uint(r)] = value / a[uint(r) * 8 + uint(r)];
    }
}
kernel void coupled_ridge_8x8_half(device const half* gram [[buffer(0)]],
                                    device const half* rhs [[buffer(1)]],
                                    device float* update [[buffer(2)]],
                                    constant float& ridge [[buffer(3)]],
                                    uint i [[thread_position_in_grid]]) {
    if (i != 0) return;
    float a[64]; float b[8];
    for (uint r = 0; r < 8; ++r) {
        b[r] = float(rhs[r]);
        for (uint c = 0; c < 8; ++c) { a[r * 8 + c] = float(gram[r * 8 + c]); }
        a[r * 8 + r] += ridge;
    }
    for (uint col = 0; col < 8; ++col) {
        uint pivot = col;
        for (uint row = col + 1; row < 8; ++row) {
            if (fabs(a[row * 8 + col]) > fabs(a[pivot * 8 + col])) { pivot = row; }
        }
        if (!isfinite(a[pivot * 8 + col]) || fabs(a[pivot * 8 + col]) <= 1.0e-12f) {
            for (uint k = 0; k < 8; ++k) { update[k] = 0.0f; }
            return;
        }
        if (pivot != col) {
            for (uint c = col; c < 8; ++c) {
                float t = a[col * 8 + c]; a[col * 8 + c] = a[pivot * 8 + c]; a[pivot * 8 + c] = t;
            }
            float t = b[col]; b[col] = b[pivot]; b[pivot] = t;
        }
        for (uint row = col + 1; row < 8; ++row) {
            float factor = a[row * 8 + col] / a[col * 8 + col];
            for (uint c = col + 1; c < 8; ++c) { a[row * 8 + c] -= factor * a[col * 8 + c]; }
            b[row] -= factor * b[col];
        }
    }
    for (int r = 7; r >= 0; --r) {
        float value = b[r];
        for (uint c = uint(r + 1); c < 8; ++c) { value -= a[uint(r) * 8 + c] * update[c]; }
        update[uint(r)] = value / a[uint(r) * 8 + uint(r)];
    }
}
struct FootpointParams { uint samples; uint segments; uint warm; uint radius; };
kernel void recompute_footpoints(device const float2* target [[buffer(0)]],
                                 device const float2* curve [[buffer(1)]],
                                 device uint* segmentOut [[buffer(2)]],
                                 device float* tOut [[buffer(3)]],
                                 constant FootpointParams& prm [[buffer(4)]],
                                 device const uint* seedSegment [[buffer(5)]],
                                 device uint* evaluationOut [[buffer(6)]],
                                 uint i [[thread_position_in_grid]]) {
    if (i >= prm.samples || prm.segments == 0u) return;
    float best = INFINITY; uint bestSegment = 0u; float bestT = 0.0f;
    float2 p = target[i];
    uint begin = 0u, end = prm.segments;
    if (prm.warm != 0u) {
        uint center = min(seedSegment[i], prm.segments - 1u);
        begin = center > prm.radius ? center - prm.radius : 0u;
        end = min(prm.segments, center + prm.radius + 1u);
    }
    evaluationOut[i] = end - begin;
    for (uint j = begin; j < end; ++j) {
        float2 a = curve[j], v = curve[j + 1u] - a;
        float vv = dot(v, v);
        float t = vv > 0.0f ? dot(p - a, v) / vv : 0.0f;
        t = clamp(t, 0.0f, 1.0f);
        float2 d = p - (a + t * v);
        float d2 = dot(d, d);
        if (isfinite(d2) && d2 < best) { best = d2; bestSegment = j; bestT = t; }
    }
    segmentOut[i] = bestSegment; tOut[i] = bestT;
}
)MSL";
#endif

MetalBracketNewtonRuntime::MetalBracketNewtonRuntime() : impl_(new MetalRuntimeImpl{}) {
#if defined(__APPLE__)
    auto* impl = static_cast<MetalRuntimeImpl*>(impl_);
    impl->device = MTLCreateSystemDefaultDevice();
    if (impl->device != nil) {
        impl->queue = [impl->device newCommandQueue];
        NSError* error = nil;
        impl->library = [impl->device newLibraryWithSource:kDiagonalUpdateShader options:nil error:&error];
        if (impl->library != nil) {
            auto buildPipe = [&](NSString* name) -> id<MTLComputePipelineState> {
                id<MTLFunction> fn = [impl->library newFunctionWithName:name];
                return fn != nil ? [impl->device newComputePipelineStateWithFunction:fn error:nil] : nil;
            };
            impl->pipeDiagonal = buildPipe(@"diagonal_ridge_update");
            impl->pipeRidge2x2 = buildPipe(@"coupled_ridge_2x2");
            impl->pipeRidge4x4 = buildPipe(@"coupled_ridge_4x4");
            impl->pipeRidge8x8 = buildPipe(@"coupled_ridge_8x8");
            impl->pipeRidge8x8Half = buildPipe(@"coupled_ridge_8x8_half");
            impl->pipeFootpoints = buildPipe(@"recompute_footpoints");
        }
    }
    impl->gramRuntime = new MetalFp16GramRuntime{};
#endif
}

MetalBracketNewtonRuntime::~MetalBracketNewtonRuntime() {
#if defined(__APPLE__)
    auto* impl = static_cast<MetalRuntimeImpl*>(impl_);
    impl->pipeDiagonal = nil;
    impl->pipeRidge2x2 = nil;
    impl->pipeRidge4x4 = nil;
    impl->pipeRidge8x8 = nil;
    impl->pipeRidge8x8Half = nil;
    impl->pipeFootpoints = nil;
    impl->library = nil;
    impl->queue = nil;
    impl->device = nil;
    delete impl->gramRuntime;
    impl->gramRuntime = nullptr;
#endif
    delete static_cast<MetalRuntimeImpl*>(impl_);
    impl_ = nullptr;
}

auto MetalBracketNewtonRuntime::available() const -> bool {
#if defined(__APPLE__)
    auto const* impl = static_cast<MetalRuntimeImpl const*>(impl_);
    return impl != nullptr && impl->device != nil && impl->queue != nil;
#else
    return false;
#endif
}

auto MetalBracketNewtonRuntime::submitGeneration(BracketNewtonConfig const& cfg,
                                                  BracketNewtonGeneration const& generation)
    -> BracketNewtonResult {
    BracketNewtonResult result;
    BracketNewtonGeneration effective = generation;
    bool cpuGramFallback = false;
#if defined(__APPLE__)
    if (cfg.useFp16Storage && !generation.fp16Basis.empty()) {
        Fp16GramConfig gramCfg{cfg.samples, cfg.parameters, true};
        Fp16GramBatch batch{generation.fp16Basis, generation.fp16Residual};
        Fp16GramResult gram = static_cast<MetalRuntimeImpl*>(impl_)->gramRuntime != nullptr
            ? static_cast<MetalRuntimeImpl*>(impl_)->gramRuntime->submit(gramCfg, batch)
            : Fp16GramResult{};
        if (!gram.ok) {
            gram = accumulateFp16GramCpu(gramCfg, batch);
            cpuGramFallback = gram.ok;
        }
        if (!gram.ok) { result.usedFallback = true; return result; }
        effective.gram = std::move(gram.gram);
        effective.rhs = std::move(gram.rhs);
    }
#endif
    result.ridge = strictTraceRidge(effective.gram, cfg.parameters, cfg.ridgeFraction);
    if (!validGeneration(cfg, effective) || !std::isfinite(result.ridge)) {
        result.usedFallback = true;
        return result;
    }
    result.nextPacingState = generation.parentPacingState;
    result.warmStartStatus = cfg.warmStart
        ? WarmStartStatus::BlockedMissingCorrespondenceState
        : WarmStartStatus::NotRequested;
#if defined(__APPLE__)
    auto* impl = static_cast<MetalRuntimeImpl*>(impl_);
    if (impl == nullptr || impl->queue == nil) {
        result.usedFallback = true;
        return result;
    }
    id<MTLComputePipelineState> pipeline = nil;
    if (cfg.useFp16Storage && cfg.parameters == 8U) {
        pipeline = impl->pipeRidge8x8Half;
    } else if (cfg.parameters == 2U) {
        pipeline = impl->pipeRidge2x2;
    } else if (cfg.parameters == 4U) {
        pipeline = impl->pipeRidge4x4;
    } else if (cfg.parameters == 8U) {
        pipeline = impl->pipeRidge8x8;
    } else {
        pipeline = impl->pipeDiagonal;
    }
    if (pipeline == nil) {
        result.usedFallback = true;
        return result;
    }
    bool const coupled4 = cfg.parameters == 4U;
    bool const coupled8 = cfg.parameters == 8U;
    std::vector<float> diagonal((cfg.parameters == 2U || coupled4 || coupled8)
                                    ? cfg.parameters * cfg.parameters : cfg.parameters);
    std::vector<float> rhs(cfg.parameters);
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
        diagonal[i] = static_cast<float>(effective.gram[i]);
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) { rhs[i] = static_cast<float>(effective.rhs[i]); }
    std::vector<std::uint16_t> diagonalHalf, rhsHalf;
    if (cfg.useFp16Storage && cfg.parameters == 8U) {
        diagonalHalf.resize(diagonal.size()); rhsHalf.resize(rhs.size());
        for (std::size_t i = 0; i < diagonal.size(); ++i) {
            if (!fp16::encodeFinite(diagonal[i], diagonalHalf[i])) { result.usedFallback = true; return result; }
        }
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            if (!fp16::encodeFinite(rhs[i], rhsHalf[i])) { result.usedFallback = true; return result; }
        }
    }
    bool const halfIngress = cfg.useFp16Storage && cfg.parameters == 8U;
    void const* gramData = halfIngress ? static_cast<void const*>(diagonalHalf.data())
                                       : static_cast<void const*>(diagonal.data());
    void const* rhsData = halfIngress ? static_cast<void const*>(rhsHalf.data())
                                      : static_cast<void const*>(rhs.data());
    NSUInteger const gramBytes = halfIngress ? diagonalHalf.size() * sizeof(std::uint16_t)
                                             : diagonal.size() * sizeof(float);
    NSUInteger const rhsBytes = halfIngress ? rhsHalf.size() * sizeof(std::uint16_t)
                                            : rhs.size() * sizeof(float);
    id<MTLBuffer> gramBuffer = [impl->device newBufferWithBytes:gramData
                                                            length:gramBytes
                                                           options:MTLResourceStorageModeShared];
    id<MTLBuffer> rhsBuffer = [impl->device newBufferWithBytes:rhsData
                                                           length:rhsBytes
                                                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> updateBuffer = [impl->device newBufferWithLength:rhs.size() * sizeof(float)
                                                             options:MTLResourceStorageModeShared];
    if (gramBuffer == nil || rhsBuffer == nil || updateBuffer == nil) {
        result.usedFallback = true;
        return result;
    }
    float ridge = static_cast<float>(result.ridge);
    id<MTLCommandBuffer> command = [impl->queue commandBuffer];
    if (command == nil) {
        result.usedFallback = true;
        return result;
    }
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:gramBuffer offset:0 atIndex:0];
    [encoder setBuffer:rhsBuffer offset:0 atIndex:1];
    [encoder setBuffer:updateBuffer offset:0 atIndex:2];
    [encoder setBytes:&ridge length:sizeof(ridge) atIndex:3];
    MTLSize grid = MTLSizeMake((cfg.parameters == 2U || coupled4 || coupled8) ? 1U : cfg.parameters, 1, 1);
    NSUInteger width = std::max<NSUInteger>(1, pipeline.threadExecutionWidth);
    MTLSize threads = MTLSizeMake(std::min<NSUInteger>(width, cfg.parameters), 1, 1);
    [encoder dispatchThreads:grid threadsPerThreadgroup:threads];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status == MTLCommandBufferStatusError) { return result; }
    float const* values = static_cast<float const*>(updateBuffer.contents);
    result.update.resize(cfg.parameters);
    for (std::size_t i = 0; i < cfg.parameters; ++i) {
        if (!std::isfinite(values[i])) {
            result.update.clear();
            result.usedFallback = true;
            return result;
        }
        result.update[i] = values[i];
    }
    result.ok = true;
    result.usedFallback = cpuGramFallback;
    if (!generation.targetPoints.empty() && generation.candidateCurveSamples.size() >= 2U) {
        struct HostFloat2 { float x; float y; };
        std::vector<HostFloat2> target, curve;
        target.reserve(generation.targetPoints.size()); curve.reserve(generation.candidateCurveSamples.size());
        for (auto const& p : generation.targetPoints) target.push_back({static_cast<float>(p.x), static_cast<float>(p.y)});
        for (auto const& p : generation.candidateCurveSamples) curve.push_back({static_cast<float>(p.x), static_cast<float>(p.y)});
        NSUInteger const sampleBytes = target.size() * sizeof(HostFloat2);
        id<MTLComputePipelineState> footPipeline = impl->pipeFootpoints;
        id<MTLBuffer> targetBuffer = footPipeline != nil ? [impl->device newBufferWithBytes:target.data()
            length:target.size() * sizeof(HostFloat2) options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> curveBuffer = footPipeline != nil ? [impl->device newBufferWithBytes:curve.data()
            length:curve.size() * sizeof(HostFloat2) options:MTLResourceStorageModeShared] : nil;
#if 0
        id<MTLBuffer> segBuffer = footPipeline != nil ? [impl->device newBufferWithLength(target.size() * sizeof(std::uint32_t)
            options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> tBuffer = footPipeline != nil ? [impl->device newBufferWithLength(target.size() * sizeof(float)
            options:MTLResourceStorageModeShared] : nil;
        id<MTLCommandBuffer> footCommand = footPipeline != nil ? [impl->queue commandBuffer] : nil;
#endif
        id<MTLBuffer> segBuffer2 = footPipeline != nil
            ? [impl->device newBufferWithLength:sampleBytes options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> tBuffer2 = footPipeline != nil
            ? [impl->device newBufferWithLength:sampleBytes options:MTLResourceStorageModeShared] : nil;
        std::vector<std::uint32_t> seedSegments(target.size(), 0U);
        bool const useWarmFootpoints = cfg.warmStart &&
            generation.footpoints.segmentIndex.size() == target.size() &&
            generation.footpoints.segmentT.size() == target.size();
        if (useWarmFootpoints) {
            for (std::size_t i = 0; i < target.size(); ++i)
                seedSegments[i] = static_cast<std::uint32_t>(generation.footpoints.segmentIndex[i]);
        }
        id<MTLBuffer> seedBuffer = footPipeline != nil ? [impl->device newBufferWithBytes:seedSegments.data()
            length:seedSegments.size() * sizeof(std::uint32_t) options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> evaluationBuffer = footPipeline != nil ? [impl->device newBufferWithLength:target.size() * sizeof(std::uint32_t)
            options:MTLResourceStorageModeShared] : nil;
        id<MTLCommandBuffer> footCommand2 = footPipeline != nil ? [impl->queue commandBuffer] : nil;
        if (footCommand2 != nil && targetBuffer != nil && curveBuffer != nil && segBuffer2 != nil &&
            tBuffer2 != nil && seedBuffer != nil && evaluationBuffer != nil) {
            id<MTLComputeCommandEncoder> footEncoder = [footCommand2 computeCommandEncoder];
            [footEncoder setComputePipelineState:footPipeline];
            [footEncoder setBuffer:targetBuffer offset:0 atIndex:0];
            [footEncoder setBuffer:curveBuffer offset:0 atIndex:1];
            [footEncoder setBuffer:segBuffer2 offset:0 atIndex:2];
            [footEncoder setBuffer:tBuffer2 offset:0 atIndex:3];
            struct FootpointParamsHost { std::uint32_t samples; std::uint32_t segments; std::uint32_t warm; std::uint32_t radius; } params{
                static_cast<std::uint32_t>(target.size()), static_cast<std::uint32_t>(curve.size() - 1U),
                useWarmFootpoints ? 1U : 0U, static_cast<std::uint32_t>(cfg.warmSegmentRadius)};
            [footEncoder setBytes:&params length:sizeof(params) atIndex:4];
            [footEncoder setBuffer:seedBuffer offset:0 atIndex:5];
            [footEncoder setBuffer:evaluationBuffer offset:0 atIndex:6];
            NSUInteger width = std::max<NSUInteger>(1U, footPipeline.threadExecutionWidth);
            [footEncoder dispatchThreads:MTLSizeMake(target.size(), 1, 1)
                threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(width, target.size()), 1, 1)];
            [footEncoder endEncoding]; [footCommand2 commit]; [footCommand2 waitUntilCompleted];
            auto const* seg = static_cast<std::uint32_t const*>(segBuffer2.contents);
            auto const* ts = static_cast<float const*>(tBuffer2.contents);
            auto const* evaluations = static_cast<std::uint32_t const*>(evaluationBuffer.contents);
            if (footCommand2.status != MTLCommandBufferStatusError) {
                for (std::size_t i = 0; i < target.size(); ++i)
                    result.correspondenceEvaluations += evaluations[i];
                result.nextFootpoints.segmentIndex.resize(target.size());
                result.nextFootpoints.segmentT.resize(target.size());
                for (std::size_t i = 0; i < target.size(); ++i) {
                    if (!std::isfinite(ts[i]) || ts[i] < 0.0f || ts[i] > 1.0f || seg[i] >= params.segments) {
                        result.nextFootpoints.segmentIndex.clear(); result.nextFootpoints.segmentT.clear(); break;
                    }
                    result.nextFootpoints.segmentIndex[i] = seg[i]; result.nextFootpoints.segmentT[i] = ts[i];
                }
                if (!result.nextFootpoints.segmentIndex.empty()) {
                    result.nextFootpoints.completedSweeps = std::min(cfg.sweeps,
                        generation.footpoints.completedSweeps + 1U);
                    result.completedSweeps = result.nextFootpoints.completedSweeps;
                    if (useWarmFootpoints) result.warmStartStatus = WarmStartStatus::StateAccepted;
                }
            }
        }
    }
    if (result.nextFootpoints.segmentIndex.empty() && !generation.footpoints.segmentIndex.empty()) {
        result.nextFootpoints = generation.footpoints;
        result.nextFootpoints.completedSweeps = std::min(
            cfg.sweeps, generation.footpoints.completedSweeps + 1U);
        result.completedSweeps = result.nextFootpoints.completedSweeps;
    }
#else
    result.usedFallback = true;
#endif
    return result;
}

#if defined(__APPLE__)
auto MetalBracketNewtonRuntime::submitGenerationAsync(
    BracketNewtonConfig const& cfg, BracketNewtonGeneration const& generation,
    std::function<void(BracketNewtonResult)> completion) -> bool {
    auto fail = [&completion]() {
        BracketNewtonResult result;
        result.usedFallback = true;
        completion(std::move(result));
        return false;
    };
    auto* impl = static_cast<MetalRuntimeImpl*>(impl_);
    double const ridge = strictTraceRidge(generation.gram, cfg.parameters, cfg.ridgeFraction);
    if (impl == nullptr || impl->device == nil || impl->queue == nil ||
        !validGeneration(cfg, generation) || !std::isfinite(ridge) || cfg.useFp16Storage) return fail();
    id<MTLComputePipelineState> pipeline = nil;
    if (cfg.parameters == 2U) {
        pipeline = impl->pipeRidge2x2;
    } else if (cfg.parameters == 4U) {
        pipeline = impl->pipeRidge4x4;
    } else if (cfg.parameters == 8U) {
        pipeline = impl->pipeRidge8x8;
    } else {
        pipeline = impl->pipeDiagonal;
    }
    if (pipeline == nil) return fail();
    std::vector<float> gram(generation.gram.size()), rhs(generation.rhs.size());
    for (std::size_t i = 0; i < gram.size(); ++i) gram[i] = static_cast<float>(generation.gram[i]);
    for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] = static_cast<float>(generation.rhs[i]);
    id<MTLBuffer> gramBuffer = [impl->device newBufferWithBytes:gram.data()
                                                            length:gram.size() * sizeof(float)
                                                           options:MTLResourceStorageModeShared];
    id<MTLBuffer> rhsBuffer = [impl->device newBufferWithBytes:rhs.data()
                                                           length:rhs.size() * sizeof(float)
                                                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> updateBuffer = [impl->device newBufferWithLength:rhs.size() * sizeof(float)
                                                             options:MTLResourceStorageModeShared];
    id<MTLCommandBuffer> command = [impl->queue commandBuffer];
    if (gramBuffer == nil || rhsBuffer == nil || updateBuffer == nil || command == nil) return fail();
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:gramBuffer offset:0 atIndex:0];
    [encoder setBuffer:rhsBuffer offset:0 atIndex:1];
    [encoder setBuffer:updateBuffer offset:0 atIndex:2];
    float ridgeFloat = static_cast<float>(ridge);
    std::vector<double> parentPacing = generation.parentPacingState;
    BracketNewtonFootpointState inputFootpoints = generation.footpoints;
    bool const warmStart = cfg.warmStart;
    auto callback = std::move(completion);
    [encoder setBytes:&ridgeFloat length:sizeof(ridgeFloat) atIndex:3];
    MTLSize grid = MTLSizeMake((cfg.parameters == 2U || cfg.parameters == 4U || cfg.parameters == 8U)
                                   ? 1U : cfg.parameters, 1, 1);
    NSUInteger width = std::max<NSUInteger>(1U, pipeline.threadExecutionWidth);
    [encoder dispatchThreads:grid threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(width, cfg.parameters), 1, 1)];
    [encoder endEncoding];
    [command addCompletedHandler:^(id<MTLCommandBuffer> completed) {
        BracketNewtonResult result;
        result.ridge = ridge;
        result.nextPacingState = parentPacing;
        result.warmStartStatus = warmStart ? WarmStartStatus::BlockedMissingCorrespondenceState
                                                : WarmStartStatus::NotRequested;
        if (completed.status == MTLCommandBufferStatusError) {
            result.usedFallback = true;
            callback(std::move(result));
            return;
        }
        float const* values = static_cast<float const*>(updateBuffer.contents);
        result.update.resize(cfg.parameters);
        for (std::size_t i = 0; i < cfg.parameters; ++i) {
            if (!std::isfinite(values[i])) {
                result.update.clear(); result.usedFallback = true;
                callback(std::move(result)); return;
            }
            result.update[i] = values[i];
        }
        result.ok = true;
        if (!inputFootpoints.segmentIndex.empty()) {
            result.nextFootpoints = inputFootpoints;
            result.nextFootpoints.completedSweeps = std::min(
                cfg.sweeps, inputFootpoints.completedSweeps + 1U);
            result.completedSweeps = result.nextFootpoints.completedSweeps;
        }
        callback(std::move(result));
    }];
    [command commit];
    return true;
}
#endif

}
