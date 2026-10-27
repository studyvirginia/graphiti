#include <metal_stdlib>
using namespace metal;
struct Params { uint P; uint M; uint Nc; };

inline void segEval(float2 p, float2 a, float2 v, float lenSq, float lenSqEps,
                     thread float& d2Out, thread float2& diffOut, thread float& tOut) {
    float t = 0.0f;
    if (lenSq > lenSqEps) {
        float traw = ((p.x - a.x) * v.x + (p.y - a.y) * v.y) / lenSq;
        t = clamp(traw, 0.0f, 1.0f);
    }
    float2 pa = p - a;
    float paSq = pa.x * pa.x + pa.y * pa.y;
    float paDotV = pa.x * v.x + pa.y * v.y;
    float d2 = paSq - 2.0f * t * paDotV + t * t * lenSq;
    d2Out = max(d2, 0.0f);
    diffOut = pa - t * v;
    tOut = t;
}

inline void segEvalHalfArith(float2 pF, float2 aF, float2 bF, float lenSqEpsF,
                              thread float& d2Out, thread half2& diffOut, thread float& tOut,
                              thread half2& normalOut) {
    half2 p = half2(pF), a = half2(aF), b = half2(bF);
    half2 v = b - a;
    half lenSq = v.x * v.x + v.y * v.y;
    half lenSqEps = half(lenSqEpsF);
    half t = 0.0h;
    if (lenSq > lenSqEps) {
        half traw = ((p.x - a.x) * v.x + (p.y - a.y) * v.y) / lenSq;
        t = clamp(traw, half(0.0h), half(1.0h));
    }
    half2 pa = p - a;
    half paSq = pa.x * pa.x + pa.y * pa.y;
    half paDotV = pa.x * v.x + pa.y * v.y;
    half d2 = paSq - half(2.0h) * t * paDotV + t * t * lenSq;
    d2Out = float(max(d2, half(0.0h)));
    diffOut = pa - t * v;
    tOut = float(t);
    half2 normal = half2(0.0h, 1.0h);
    if (lenSq > lenSqEps) {
        half invLen = half(1.0h) / sqrt(lenSq);
        half2 tangent = v * invLen;
        normal = half2(-tangent.y, tangent.x);
    }
    normalOut = normal;
}

inline void segEvalLocalHalf(float2 vF, float2 dF, float lenSqEpsF,
                              thread float& d2Out, thread half2& diffOut, thread float& tOut,
                              thread half2& normalOut) {
    half2 v = half2(vF);
    half2 d = half2(dF);
    half lenSq = v.x * v.x + v.y * v.y;
    half lenSqEps = half(lenSqEpsF);
    half t = 0.0h;
    if (lenSq > lenSqEps) {
        half traw = (d.x * v.x + d.y * v.y) / lenSq;
        t = clamp(traw, half(0.0h), half(1.0h));
    }
    half paSq = d.x * d.x + d.y * d.y;
    half paDotV = d.x * v.x + d.y * v.y;
    half d2 = paSq - half(2.0h) * t * paDotV + t * t * lenSq;
    d2Out = float(max(d2, half(0.0h)));
    diffOut = d - t * v;
    tOut = float(t);
    half2 normal = half2(0.0h, 1.0h);
    if (lenSq > lenSqEps) {
        half invLen = half(1.0h) / sqrt(lenSq);
        half2 tangent = v * invLen;
        normal = half2(-tangent.y, tangent.x);
    }
    normalOut = normal;
}

kernel void fusedThreadgroup(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float2 v = b - a;
            float lenSq = v.x * v.x + v.y * v.y;
            float d2; float2 diff; float tc;
            segEval(pt, a, v, lenSq, lenSqEps, d2, diff, tc);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            float2 normal = float2(0.0f, 1.0f);
            if (lenSq > lenSqEps) {
                float invLen = 1.0f / sqrt(lenSq);
                float2 tangent = v * invLen;
                normal = float2(-tangent.y, tangent.x);
            }
            float signedVal = diff.x * normal.x + diff.y * normal.y;
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
        }
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupUnsigned(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float blendedAbs = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float2 v = b - a;
            float lenSq = v.x * v.x + v.y * v.y;
            float d2; float2 diff; float tc;
            segEval(pt, a, v, lenSq, lenSqEps, d2, diff, tc);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            float2 normal = float2(0.0f, 1.0f);
            if (lenSq > lenSqEps) {
                float invLen = 1.0f / sqrt(lenSq);
                float2 tangent = v * invLen;
                normal = float2(-tangent.y, tangent.x);
            }
            float signedVal = diff.x * normal.x + diff.y * normal.y;
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                blendedAbs *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
            blendedAbs += e * fabs(signedVal);
        }
        float u = kEps / 2.0f;
        float nRound = 2.0f * float(prm.Nc);
        float gammaN = (nRound * u) / (1.0f - nRound * u);
        float sgnI = (fabs(blended) > gammaN * blendedAbs && blended < 0.0f) ? -1.0f : 1.0f;
        blended = sgnI * blendedAbs;
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16Storage(device const half2* cand [[buffer(0)]], device const half2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const half2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = float2(data[i]);
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = float2(c[j]), b = float2(c[j + 1]);
            float2 v = b - a;
            float lenSq = v.x * v.x + v.y * v.y;
            float d2; float2 diff; float tc;
            segEval(pt, a, v, lenSq, lenSqEps, d2, diff, tc);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            float2 normal = float2(0.0f, 1.0f);
            if (lenSq > lenSqEps) {
                float invLen = 1.0f / sqrt(lenSq);
                float2 tangent = v * invLen;
                normal = float2(-tangent.y, tangent.x);
            }
            float signedVal = diff.x * normal.x + diff.y * normal.y;
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
        }
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16StorageUnsigned(device const half2* cand [[buffer(0)]], device const half2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const half2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = float2(data[i]);
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float blendedAbs = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = float2(c[j]), b = float2(c[j + 1]);
            float2 v = b - a;
            float lenSq = v.x * v.x + v.y * v.y;
            float d2; float2 diff; float tc;
            segEval(pt, a, v, lenSq, lenSqEps, d2, diff, tc);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            float2 normal = float2(0.0f, 1.0f);
            if (lenSq > lenSqEps) {
                float invLen = 1.0f / sqrt(lenSq);
                float2 tangent = v * invLen;
                normal = float2(-tangent.y, tangent.x);
            }
            float signedVal = diff.x * normal.x + diff.y * normal.y;
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                blendedAbs *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
            blendedAbs += e * fabs(signedVal);
        }
        float u = kEps / 2.0f;
        float nRound = 2.0f * float(prm.Nc);
        float gammaN = (nRound * u) / (1.0f - nRound * u);
        float sgnI = (fabs(blended) > gammaN * blendedAbs && blended < 0.0f) ? -1.0f : 1.0f;
        blended = sgnI * blendedAbs;
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16Arith(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float d2; half2 diffH; float tc; half2 normalH;
            segEvalHalfArith(pt, a, b, lenSqEps, d2, diffH, tc, normalH);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            half signedValH = diffH.x * normalH.x + diffH.y * normalH.y;
            float signedVal = float(signedValH);
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
        }
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16ArithUnsigned(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float blendedAbs = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float d2; half2 diffH; float tc; half2 normalH;
            segEvalHalfArith(pt, a, b, lenSqEps, d2, diffH, tc, normalH);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            half signedValH = diffH.x * normalH.x + diffH.y * normalH.y;
            float signedVal = float(signedValH);
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                blendedAbs *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
            blendedAbs += e * fabs(signedVal);
        }
        float u = kEps / 2.0f;
        float nRound = 2.0f * float(prm.Nc);
        float gammaN = (nRound * u) / (1.0f - nRound * u);
        float sgnI = (fabs(blended) > gammaN * blendedAbs && blended < 0.0f) ? -1.0f : 1.0f;
        blended = sgnI * blendedAbs;
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16Local(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float2 v = b - a;
            float2 dLocal = pt - a;
            float d2; half2 diffH; float tc; half2 normalH;
            segEvalLocalHalf(v, dLocal, lenSqEps, d2, diffH, tc, normalH);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            half signedValH = diffH.x * normalH.x + diffH.y * normalH.y;
            float signedVal = float(signedValH);
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
        }
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

kernel void fusedThreadgroupFp16LocalUnsigned(device const float2* cand [[buffer(0)]], device const float2* data [[buffer(1)]],
                              constant Params& prm [[buffer(2)]], device const float* meanSegLen [[buffer(3)]],
                              device float* residual [[buffer(4)]], device float* fitness [[buffer(5)]],
                              device float* tauOut [[buffer(6)]],
                              uint p [[threadgroup_position_in_grid]],
                              uint i [[thread_position_in_threadgroup]],
                              uint tgSize [[threads_per_threadgroup]],
                              uint simdLane [[thread_index_in_simdgroup]],
                              uint simdGroupId [[simdgroup_index_in_threadgroup]],
                              threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    device const float2* c = cand + p * prm.Nc;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;
    float softMinScale = max(meanSegLen[p], kEps);

    float contrib = 0.0f;
    if (i < prm.M) {
        float2 pt = data[i];
        float mRun = FLT_MAX;
        float sumExp = 0.0f, blended = 0.0f, blendedTau = 0.0f;
        float blendedAbs = 0.0f;
        float hardMin = FLT_MAX;
        for (uint j = 0; j + 1 < prm.Nc; ++j) {
            float2 a = c[j], b = c[j + 1];
            float2 v = b - a;
            float2 dLocal = pt - a;
            float d2; half2 diffH; float tc; half2 normalH;
            segEvalLocalHalf(v, dLocal, lenSqEps, d2, diffH, tc, normalH);
            hardMin = min(hardMin, d2);
            float dVal = sqrt(max(d2, 0.0f));

            half signedValH = diffH.x * normalH.x + diffH.y * normalH.y;
            float signedVal = float(signedValH);
            float tauVal = float(j) + tc;

            if (dVal < mRun) {
                float rescale = exp(-(mRun - dVal) / softMinScale);
                sumExp *= rescale;
                blended *= rescale;
                blendedTau *= rescale;
                blendedAbs *= rescale;
                mRun = dVal;
            }
            float e = exp(-(dVal - mRun) / softMinScale);
            sumExp += e;
            blended += e * signedVal;
            blendedTau += e * tauVal;
            blendedAbs += e * fabs(signedVal);
        }
        float u = kEps / 2.0f;
        float nRound = 2.0f * float(prm.Nc);
        float gammaN = (nRound * u) / (1.0f - nRound * u);
        float sgnI = (fabs(blended) > gammaN * blendedAbs && blended < 0.0f) ? -1.0f : 1.0f;
        blended = sgnI * blendedAbs;
        residual[p * prm.M + i] = (sumExp > 0.0f) ? (blended / sumExp) : 0.0f;
        tauOut[p * prm.M + i] = (sumExp > 0.0f) ? (blendedTau / sumExp) : 0.0f;
        contrib = hardMin;
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { fitness[p] = total; }
    }
}

struct ReverseParams { uint P; uint Nc; uint Ns; uint K; uint ptClamp; };

inline float pointLossMetal(float r2, float s2, float nu, bool ptClamp) {
    if (!(s2 > 0.0f) || !isfinite(s2)) { s2 = 1.0f; }
    float val;
    if (nu <= 0.0f) { val = r2 / s2; }
    else { val = 0.5f * (nu + 2.0f) * log(1.0f + r2 / (nu * s2)); }
    if (ptClamp) { val = max(val, 1.0f); }
    return val;
}

kernel void reverseHalfThreadgroup(device const float2* candVerts [[buffer(0)]],
                                    device const float2* drawn [[buffer(1)]],
                                    constant ReverseParams& prm [[buffer(2)]],
                                    device const float* s2Buf [[buffer(3)]],
                                    device const float* nuBuf [[buffer(4)]],
                                    device float* reverseFitness [[buffer(5)]],
                                    uint p [[threadgroup_position_in_grid]],
                                    uint i [[thread_position_in_threadgroup]],
                                    uint tgSize [[threads_per_threadgroup]],
                                    uint simdLane [[thread_index_in_simdgroup]],
                                    uint simdGroupId [[simdgroup_index_in_threadgroup]],
                                    threadgroup float* partial [[threadgroup(0)]]) {
    if (p >= prm.P) return;
    constexpr float kEps = 1.1920929e-7f;
    constexpr float lenSqEps = kEps * 100.0f;

    float contrib = 0.0f;
    if (i < prm.Nc) {
        float2 q = candVerts[p * prm.Nc + i];
        float hardMin = FLT_MAX;
        uint nearestVertex = 0;
        for (uint j = 0; j + 1 < prm.Ns; ++j) {
            float2 a = drawn[j], b = drawn[j + 1];
            float2 v = b - a;
            float lenSq = v.x * v.x + v.y * v.y;
            float d2; float2 diff; float tc;
            segEval(q, a, v, lenSq, lenSqEps, d2, diff, tc);
            if (d2 < hardMin) {
                hardMin = d2;
                nearestVertex = (tc < 0.5f) ? j : (j + 1);
            }
        }
        if (prm.K == 0) {
            contrib = hardMin;
        } else {
            uint k = min((nearestVertex * prm.K) / (prm.Ns - 1u), prm.K - 1u);
            contrib = pointLossMetal(hardMin, s2Buf[k], nuBuf[k], prm.ptClamp != 0);
        }
    }

    float simdSum = simd_sum(contrib);
    uint numSimdGroups = (tgSize + 31) / 32;
    if (simdLane == 0) { partial[simdGroupId] = simdSum; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (simdGroupId == 0) {
        float v = (simdLane < numSimdGroups) ? partial[simdLane] : 0.0f;
        float total = simd_sum(v);
        if (simdLane == 0) { reverseFitness[p] = total; }
    }
}
