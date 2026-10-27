#include "sketch/fp16_gram_runtime.hpp"

#include <Metal/Metal.h>

namespace Sketch {
namespace {
NSString* const kShader = @R"MSL(
#include <metal_stdlib>
using namespace metal;
kernel void fp16_gram(device const half* basis [[buffer(0)]], device const half* residual [[buffer(1)]],
                      device float* output [[buffer(2)]], constant uint& rows [[buffer(3)]],
                      constant uint& cols [[buffer(4)]], uint gid [[thread_position_in_grid]]) {
    uint total = cols * cols + cols;
    if (gid >= total) return;
    if (gid < cols * cols) {
        uint r = gid / cols, c = gid % cols; float sum = 0.0f;
        for (uint i = 0; i < rows; ++i) { sum += float(basis[i * cols + r]) * float(basis[i * cols + c]); }
        output[gid] = sum;
    } else {
        uint r = gid - cols * cols; float sum = 0.0f;
        for (uint i = 0; i < rows; ++i) { sum += float(basis[i * cols + r]) * float(residual[i]); }
        output[gid] = sum;
    }
}
kernel void fp32_gram(device const float* basis [[buffer(0)]], device const float* residual [[buffer(1)]],
                      device float* output [[buffer(2)]], constant uint& rows [[buffer(3)]],
                      constant uint& cols [[buffer(4)]], uint gid [[thread_position_in_grid]]) {
    uint total = cols * cols + cols;
    if (gid >= total) return;
    if (gid < cols * cols) {
        uint r = gid / cols, c = gid % cols; float sum = 0.0f;
        for (uint i = 0; i < rows; ++i) { sum += basis[i * cols + r] * basis[i * cols + c]; }
        output[gid] = sum;
    } else {
        uint r = gid - cols * cols; float sum = 0.0f;
        for (uint i = 0; i < rows; ++i) { sum += basis[i * cols + r] * residual[i]; }
        output[gid] = sum;
    }
}
)MSL";
struct Impl { id<MTLDevice> device{nil}; id<MTLCommandQueue> queue{nil}; };
}

MetalFp16GramRuntime::MetalFp16GramRuntime() : impl_(new Impl{}) {
    auto* impl = static_cast<Impl*>(impl_);
    impl->device = MTLCreateSystemDefaultDevice();
    if (impl->device != nil) { impl->queue = [impl->device newCommandQueue]; }
}

MetalFp16GramRuntime::~MetalFp16GramRuntime() {
    auto* impl = static_cast<Impl*>(impl_);
    impl->queue = nil; impl->device = nil; delete impl; impl_ = nullptr;
}

auto MetalFp16GramRuntime::available() const -> bool {
    auto const* impl = static_cast<Impl const*>(impl_);
    return impl != nullptr && impl->device != nil && impl->queue != nil;
}

auto MetalFp16GramRuntime::submit(Fp16GramConfig const& cfg, Fp16GramBatch const& batch)
    -> Fp16GramResult {
    Fp16GramResult out;
    if (!validFp16Gram(cfg, batch) || !available()) { out.usedFallback = true; return out; }
    auto* impl = static_cast<Impl*>(impl_);
    std::vector<std::uint16_t> basisHalf, residualHalf;
    std::vector<float> basisFloat, residualFloat;
    if (cfg.useFp16Storage) {
        basisHalf.resize(batch.basis.size()); residualHalf.resize(batch.residual.size());
        for (std::size_t i = 0; i < basisHalf.size(); ++i) {
            if (!fp16::encodeFinite(static_cast<float>(batch.basis[i]), basisHalf[i])) { out.usedFallback = true; return out; }
        }
        for (std::size_t i = 0; i < residualHalf.size(); ++i) {
            if (!fp16::encodeFinite(static_cast<float>(batch.residual[i]), residualHalf[i])) { out.usedFallback = true; return out; }
        }
    } else {
        basisFloat.assign(batch.basis.begin(), batch.basis.end());
        residualFloat.assign(batch.residual.begin(), batch.residual.end());
    }
    bool const half = cfg.useFp16Storage;
    NSError* error = nil;
    id<MTLLibrary> library = [impl->device newLibraryWithSource:kShader options:nil error:&error];
    id<MTLFunction> function = [library newFunctionWithName:(half ? @"fp16_gram" : @"fp32_gram")];
    id<MTLComputePipelineState> pipeline = function != nil
        ? [impl->device newComputePipelineStateWithFunction:function error:&error] : nil;
    if (pipeline == nil) { out.usedFallback = true; return out; }
    std::size_t const total = cfg.columns * cfg.columns + cfg.columns;
    void const* basisData = half ? static_cast<void const*>(basisHalf.data()) : static_cast<void const*>(basisFloat.data());
    void const* residualData = half ? static_cast<void const*>(residualHalf.data()) : static_cast<void const*>(residualFloat.data());
    NSUInteger const basisBytes = half ? basisHalf.size() * sizeof(std::uint16_t) : basisFloat.size() * sizeof(float);
    NSUInteger const residualBytes = half ? residualHalf.size() * sizeof(std::uint16_t) : residualFloat.size() * sizeof(float);
    id<MTLBuffer> basisBuffer = [impl->device newBufferWithBytes:basisData
        length:basisBytes options:MTLResourceStorageModeShared];
    id<MTLBuffer> residualBuffer = [impl->device newBufferWithBytes:residualData
        length:residualBytes options:MTLResourceStorageModeShared];
    id<MTLBuffer> outputBuffer = [impl->device newBufferWithLength:total * sizeof(float)
        options:MTLResourceStorageModeShared];
    if (basisBuffer == nil || residualBuffer == nil || outputBuffer == nil) { out.usedFallback = true; return out; }
    std::uint32_t rows = static_cast<std::uint32_t>(cfg.rows), cols = static_cast<std::uint32_t>(cfg.columns);
    id<MTLCommandBuffer> command = [impl->queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline]; [encoder setBuffer:basisBuffer offset:0 atIndex:0];
    [encoder setBuffer:residualBuffer offset:0 atIndex:1]; [encoder setBuffer:outputBuffer offset:0 atIndex:2];
    [encoder setBytes:&rows length:sizeof(rows) atIndex:3]; [encoder setBytes:&cols length:sizeof(cols) atIndex:4];
    NSUInteger requestedWidth = cfg.threadgroupWidth == 0U
        ? pipeline.threadExecutionWidth : static_cast<NSUInteger>(cfg.threadgroupWidth);
    NSUInteger threadWidth = std::max<NSUInteger>(1,
        std::min<NSUInteger>({requestedWidth, pipeline.maxTotalThreadsPerThreadgroup, total}));
    [encoder dispatchThreads:MTLSizeMake(total, 1, 1)
       threadsPerThreadgroup:MTLSizeMake(threadWidth, 1, 1)];
    [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted) { out.usedFallback = true; return out; }
    auto* values = static_cast<float*>(outputBuffer.contents);
    out.gram.resize(cfg.columns * cfg.columns); out.rhs.resize(cfg.columns);
    for (std::size_t i = 0; i < out.gram.size(); ++i) { out.gram[i] = values[i]; }
    for (std::size_t i = 0; i < out.rhs.size(); ++i) { out.rhs[i] = values[out.gram.size() + i]; }
    for (double value : out.gram) { if (!std::isfinite(value)) { out = {}; out.usedFallback = true; return out; } }
    for (double value : out.rhs) { if (!std::isfinite(value)) { out = {}; out.usedFallback = true; return out; } }
    out.ok = true;
    return out;
}
}
