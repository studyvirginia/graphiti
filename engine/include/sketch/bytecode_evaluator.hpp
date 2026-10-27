#ifndef SKETCH_BYTECODE_EVALUATOR_HPP
#define SKETCH_BYTECODE_EVALUATOR_HPP

#include <operon/operon/core/tree.hpp>
#include <operon/operon/core/node.hpp>
#include <operon/operon/core/types.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <span>
#include <cstdint>
#include <cstdlib>
#include <bit>

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#if defined(OPERON_MATH_EVE)
#include <operon/operon/interpreter/backend/eve/functions.hpp>
#endif

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif

namespace Sketch {

enum class BytecodeOp : uint8_t {
    Constant,
    Variable,
    Add,
    Sub,
    Mul,
    Div,
    Sin,
    Cos,
    Exp,
    Log,
    Sqrt,
    Square
};

struct BytecodeInst {
    BytecodeOp op;
    uint8_t arity{2};
    float imm{0.0f};
    float weight{1.0f};
};

class BytecodeProgram {
public:
    std::vector<BytecodeInst> insts;
    bool valid{false};

    struct CoeffBinding {
        uint16_t instIndex;
        bool isConstant;
    };
    std::vector<CoeffBinding> bindings;

    void UpdateCoefficients(std::span<Operon::Scalar const> coeffs)
    {
        for (std::size_t i = 0; i < bindings.size() && i < coeffs.size(); ++i) {
            auto const& b = bindings[i];
            if (b.isConstant) {
                insts[b.instIndex].imm = static_cast<float>(coeffs[i]);
            } else {
                insts[b.instIndex].weight = static_cast<float>(coeffs[i]);
            }
        }
    }

    static auto Compile(Operon::Tree const& tree, std::span<Operon::Scalar const> coeffs) -> BytecodeProgram
    {
        BytecodeProgram prog;
        auto const& nodes = tree.Nodes();
        if (nodes.empty()) { return prog; }

        prog.insts.reserve(nodes.size());
        std::size_t coeffIdx = 0;

        for (auto const& n : nodes) {
            BytecodeInst inst{};
            inst.arity = static_cast<uint8_t>(n.Arity);
            inst.weight = static_cast<float>(n.Value);

            if (n.IsConstant()) {
                inst.op = BytecodeOp::Constant;
                if (n.Optimize) {
                    prog.bindings.push_back({static_cast<uint16_t>(prog.insts.size()), true});
                    if (coeffIdx < coeffs.size()) {
                        inst.imm = static_cast<float>(coeffs[coeffIdx++]);
                    } else {
                        inst.imm = static_cast<float>(n.Value);
                    }
                } else {
                    inst.imm = static_cast<float>(n.Value);
                }
                prog.insts.push_back(inst);
            } else if (n.IsVariable()) {
                inst.op = BytecodeOp::Variable;
                if (n.Optimize) {
                    prog.bindings.push_back({static_cast<uint16_t>(prog.insts.size()), false});
                    if (coeffIdx < coeffs.size()) {
                        inst.weight = static_cast<float>(coeffs[coeffIdx++]);
                    } else {
                        inst.weight = static_cast<float>(n.Value);
                    }
                } else {
                    inst.weight = static_cast<float>(n.Value);
                }
                prog.insts.push_back(inst);
            } else if (n.IsFunction()) {
                using BO = Operon::BuiltinOp;
                auto h = n.HashValue;
                if (h == static_cast<Operon::Hash>(BO::Add)) {
                    inst.op = BytecodeOp::Add;
                } else if (h == static_cast<Operon::Hash>(BO::Sub)) {
                    inst.op = BytecodeOp::Sub;
                } else if (h == static_cast<Operon::Hash>(BO::Mul)) {
                    inst.op = BytecodeOp::Mul;
                } else if (h == static_cast<Operon::Hash>(BO::Div)) {
                    inst.op = BytecodeOp::Div;
                } else if (h == static_cast<Operon::Hash>(BO::Sin)) {
                    inst.op = BytecodeOp::Sin;
                } else if (h == static_cast<Operon::Hash>(BO::Cos)) {
                    inst.op = BytecodeOp::Cos;
                } else if (h == static_cast<Operon::Hash>(BO::Exp)) {
                    inst.op = BytecodeOp::Exp;
                } else if (h == static_cast<Operon::Hash>(BO::Log)) {
                    inst.op = BytecodeOp::Log;
                } else if (h == static_cast<Operon::Hash>(BO::Sqrt)) {
                    inst.op = BytecodeOp::Sqrt;
                } else if (h == static_cast<Operon::Hash>(BO::Square)) {
                    inst.op = BytecodeOp::Square;
                } else {
                    prog.valid = false;
                    return prog;
                }
                prog.insts.push_back(inst);
            } else {
                prog.valid = false;
                return prog;
            }
        }
        prog.valid = true;
        return prog;
    }

    static auto opMajorEnabled() -> bool
    {
        static thread_local int v = [] {
            char const* e = std::getenv("SKETCH_OP_MAJOR");
            return (e && e[0] == '0') ? 0 : 1;
        }();
        return v != 0;
    }

    void Evaluate(std::span<Operon::Scalar const> t_grid, std::span<Operon::Scalar> out) const
    {
        if (opMajorEnabled()) { EvaluateOpMajor(t_grid, out); return; }
        std::size_t const N = t_grid.size();
        std::size_t i = 0;

#if defined(__ARM_NEON) || defined(__aarch64__)
        for (; i + 4 <= N; i += 4) {
            float32x4_t stack[32];
            int sp = 0;
            float32x4_t vt = vld1q_f32(&t_grid[i]);

            for (auto const& inst : insts) {
                switch (inst.op) {
                    case BytecodeOp::Constant:
                        stack[sp++] = vdupq_n_f32(inst.imm);
                        break;
                    case BytecodeOp::Variable: {
                        auto v = vt;
                        if (inst.weight != 1.0f) {
                            v = vmulq_f32(v, vdupq_n_f32(inst.weight));
                        }
                        stack[sp++] = v;
                        break;
                    }
                    case BytecodeOp::Add: {
                        auto b = stack[--sp];
                        auto res = vaddq_f32(stack[sp - 1], b);
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Sub: {
                        if (inst.arity == 1) {
                            auto a = stack[sp - 1];
                            auto res = vnegq_f32(a);
                            if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                            stack[sp - 1] = res;
                        } else {
                            auto b = stack[--sp];
                            auto res = vsubq_f32(b, stack[sp - 1]);
                            if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                            stack[sp - 1] = res;
                        }
                        break;
                    }
                    case BytecodeOp::Mul: {
                        auto b = stack[--sp];
                        auto res = vmulq_f32(stack[sp - 1], b);
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Div: {
                        if (inst.arity == 1) {
                            auto a = stack[sp - 1];
                            auto res = vdivq_f32(vdupq_n_f32(1.0f), a);
                            if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                            stack[sp - 1] = res;
                        } else {
                            auto b = stack[--sp];
                            auto res = vdivq_f32(b, stack[sp - 1]);
                            if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                            stack[sp - 1] = res;
                        }
                        break;
                    }
                    case BytecodeOp::Square: {
                        auto a = stack[sp - 1];
                        auto res = vmulq_f32(a, a);
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Sqrt: {
                        auto a = stack[sp - 1];
                        auto res = vsqrtq_f32(a);
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Sin: {
#if defined(OPERON_MATH_EVE)
                        using W = eve::wide<float, eve::fixed<4>>;
                        auto w = std::bit_cast<W>(stack[sp - 1]);
                        auto s = Operon::Backend::detail::FastSinCos<true>(w);
                        auto res = std::bit_cast<float32x4_t>(s);
#elif defined(__APPLE__)
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        int n4 = 4;
                        vvsinf(tmp, tmp, &n4);
                        auto res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        tmp[0] = std::sin(tmp[0]); tmp[1] = std::sin(tmp[1]);
                        tmp[2] = std::sin(tmp[2]); tmp[3] = std::sin(tmp[3]);
                        auto res = vld1q_f32(tmp);
#endif
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Cos: {
#if defined(OPERON_MATH_EVE)
                        using W = eve::wide<float, eve::fixed<4>>;
                        auto w = std::bit_cast<W>(stack[sp - 1]);
                        auto c = Operon::Backend::detail::FastSinCos<false>(w);
                        auto res = std::bit_cast<float32x4_t>(c);
#elif defined(__APPLE__)
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        int n4 = 4;
                        vvcosf(tmp, tmp, &n4);
                        auto res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        tmp[0] = std::cos(tmp[0]); tmp[1] = std::cos(tmp[1]);
                        tmp[2] = std::cos(tmp[2]); tmp[3] = std::cos(tmp[3]);
                        auto res = vld1q_f32(tmp);
#endif
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Exp: {
#if defined(OPERON_MATH_EVE)
                        using W = eve::wide<float, eve::fixed<4>>;
                        auto w = std::bit_cast<W>(stack[sp - 1]);
                        auto e = Operon::Backend::detail::FastExp(w);
                        auto res = std::bit_cast<float32x4_t>(e);
#elif defined(__APPLE__)
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        int n4 = 4;
                        vvexpf(tmp, tmp, &n4);
                        auto res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        tmp[0] = std::exp(tmp[0]); tmp[1] = std::exp(tmp[1]);
                        tmp[2] = std::exp(tmp[2]); tmp[3] = std::exp(tmp[3]);
                        auto res = vld1q_f32(tmp);
#endif
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                    case BytecodeOp::Log: {
#if defined(OPERON_MATH_EVE)
                        using W = eve::wide<float, eve::fixed<4>>;
                        auto w = std::bit_cast<W>(stack[sp - 1]);
                        auto l = Operon::Backend::detail::FastLog(w);
                        auto res = std::bit_cast<float32x4_t>(l);
#elif defined(__APPLE__)
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        int n4 = 4;
                        vvlogf(tmp, tmp, &n4);
                        auto res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4];
                        vst1q_f32(tmp, stack[sp - 1]);
                        tmp[0] = std::log(tmp[0]); tmp[1] = std::log(tmp[1]);
                        tmp[2] = std::log(tmp[2]); tmp[3] = std::log(tmp[3]);
                        auto res = vld1q_f32(tmp);
#endif
                        if (inst.weight != 1.0f) { res = vmulq_f32(res, vdupq_n_f32(inst.weight)); }
                        stack[sp - 1] = res;
                        break;
                    }
                }
            }
            if (sp > 0) {
                vst1q_f32(&out[i], stack[0]);
            }
        }
#endif

        for (; i < N; ++i) {
            float stack[32];
            int sp = 0;
            float vt = t_grid[i];
            for (auto const& inst : insts) {
                switch (inst.op) {
                    case BytecodeOp::Constant:
                        stack[sp++] = inst.imm;
                        break;
                    case BytecodeOp::Variable: {
                        float v = vt * inst.weight;
                        stack[sp++] = v;
                        break;
                    }
                    case BytecodeOp::Add: {
                        auto b = stack[--sp];
                        stack[sp - 1] = (stack[sp - 1] + b) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Sub: {
                        if (inst.arity == 1) {
                            stack[sp - 1] = (-stack[sp - 1]) * inst.weight;
                        } else {
                            auto b = stack[--sp];
                            stack[sp - 1] = (b - stack[sp - 1]) * inst.weight;
                        }
                        break;
                    }
                    case BytecodeOp::Mul: {
                        auto b = stack[--sp];
                        stack[sp - 1] = (stack[sp - 1] * b) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Div: {
                        if (inst.arity == 1) {
                            stack[sp - 1] = (1.0f / stack[sp - 1]) * inst.weight;
                        } else {
                            auto b = stack[--sp];
                            stack[sp - 1] = (b / stack[sp - 1]) * inst.weight;
                        }
                        break;
                    }
                    case BytecodeOp::Square: {
                        auto a = stack[sp - 1];
                        stack[sp - 1] = (a * a) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Sqrt: {
                        auto a = stack[sp - 1];
                        stack[sp - 1] = std::sqrt(a) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Sin: {
                        stack[sp - 1] = std::sin(stack[sp - 1]) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Cos: {
                        stack[sp - 1] = std::cos(stack[sp - 1]) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Exp: {
                        stack[sp - 1] = std::exp(stack[sp - 1]) * inst.weight;
                        break;
                    }
                    case BytecodeOp::Log: {
                        stack[sp - 1] = std::log(stack[sp - 1]) * inst.weight;
                        break;
                    }
                }
            }
            if (sp > 0) {
                out[i] = stack[0];
            }
        }
    }

    void EvaluateOpMajor(std::span<Operon::Scalar const> t_grid, std::span<Operon::Scalar> out) const
    {
        int const N = static_cast<int>(t_grid.size());
        static thread_local std::vector<Operon::Scalar> buf;
        if (static_cast<int>(buf.size()) < 32 * N) { buf.resize(static_cast<std::size_t>(32) * N); }
        Operon::Scalar* const S = buf.data();
        auto slot = [&](int s) -> Operon::Scalar* { return S + static_cast<std::size_t>(s) * N; };
        int sp = 0;

        auto scale = [&](Operon::Scalar* p, float w) {
            if (w == 1.0f) { return; }
#if defined(__ARM_NEON) || defined(__aarch64__)
            int i = 0; float32x4_t vw = vdupq_n_f32(w);
            for (; i + 4 <= N; i += 4) { vst1q_f32(&p[i], vmulq_f32(vld1q_f32(&p[i]), vw)); }
            for (; i < N; ++i) { p[i] *= w; }
#else
            for (int i = 0; i < N; ++i) { p[i] *= w; }
#endif
        };
        auto transc = [&](Operon::Scalar* a, BytecodeOp op) {
            int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
            for (; i + 4 <= N; i += 4) {
                float32x4_t v = vld1q_f32(&a[i]); float32x4_t res;
                switch (op) {
                    case BytecodeOp::Sin: {
#if defined(OPERON_MATH_EVE)
                        using Wv = eve::wide<float, eve::fixed<4>>; auto w = std::bit_cast<Wv>(v);
                        res = std::bit_cast<float32x4_t>(Operon::Backend::detail::FastSinCos<true>(w));
#elif defined(__APPLE__)
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v); int n4 = 4; vvsinf(tmp, tmp, &n4); res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v);
                        tmp[0]=std::sin(tmp[0]);tmp[1]=std::sin(tmp[1]);tmp[2]=std::sin(tmp[2]);tmp[3]=std::sin(tmp[3]); res=vld1q_f32(tmp);
#endif
                        break; }
                    case BytecodeOp::Cos: {
#if defined(OPERON_MATH_EVE)
                        using Wv = eve::wide<float, eve::fixed<4>>; auto w = std::bit_cast<Wv>(v);
                        res = std::bit_cast<float32x4_t>(Operon::Backend::detail::FastSinCos<false>(w));
#elif defined(__APPLE__)
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v); int n4 = 4; vvcosf(tmp, tmp, &n4); res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v);
                        tmp[0]=std::cos(tmp[0]);tmp[1]=std::cos(tmp[1]);tmp[2]=std::cos(tmp[2]);tmp[3]=std::cos(tmp[3]); res=vld1q_f32(tmp);
#endif
                        break; }
                    case BytecodeOp::Exp: {
#if defined(OPERON_MATH_EVE)
                        using Wv = eve::wide<float, eve::fixed<4>>; auto w = std::bit_cast<Wv>(v);
                        res = std::bit_cast<float32x4_t>(Operon::Backend::detail::FastExp(w));
#elif defined(__APPLE__)
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v); int n4 = 4; vvexpf(tmp, tmp, &n4); res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v);
                        tmp[0]=std::exp(tmp[0]);tmp[1]=std::exp(tmp[1]);tmp[2]=std::exp(tmp[2]);tmp[3]=std::exp(tmp[3]); res=vld1q_f32(tmp);
#endif
                        break; }
                    default: {
#if defined(OPERON_MATH_EVE)
                        using Wv = eve::wide<float, eve::fixed<4>>; auto w = std::bit_cast<Wv>(v);
                        res = std::bit_cast<float32x4_t>(Operon::Backend::detail::FastLog(w));
#elif defined(__APPLE__)
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v); int n4 = 4; vvlogf(tmp, tmp, &n4); res = vld1q_f32(tmp);
#else
                        alignas(16) float tmp[4]; vst1q_f32(tmp, v);
                        tmp[0]=std::log(tmp[0]);tmp[1]=std::log(tmp[1]);tmp[2]=std::log(tmp[2]);tmp[3]=std::log(tmp[3]); res=vld1q_f32(tmp);
#endif
                        break; }
                }
                vst1q_f32(&a[i], res);
            }
#endif
            for (; i < N; ++i) {
                switch (op) {
                    case BytecodeOp::Sin: a[i] = std::sin(a[i]); break;
                    case BytecodeOp::Cos: a[i] = std::cos(a[i]); break;
                    case BytecodeOp::Exp: a[i] = std::exp(a[i]); break;
                    default: a[i] = std::log(a[i]); break;
                }
            }
        };

        for (auto const& inst : insts) {
            switch (inst.op) {
                case BytecodeOp::Constant: { Operon::Scalar* d = slot(sp++); for (int i = 0; i < N; ++i) { d[i] = inst.imm; } break; }
                case BytecodeOp::Variable: { Operon::Scalar* d = slot(sp++); float w = inst.weight; for (int i = 0; i < N; ++i) { d[i] = t_grid[i] * w; } break; }
                case BytecodeOp::Add: { Operon::Scalar* b = slot(--sp); Operon::Scalar* a = slot(sp - 1);
                    int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                    for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vaddq_f32(vld1q_f32(&a[i]), vld1q_f32(&b[i]))); }
#endif
                    for (; i < N; ++i) { a[i] += b[i]; } scale(a, inst.weight); break; }
                case BytecodeOp::Sub: { Operon::Scalar* a = slot(sp - 1);
                    if (inst.arity == 1) {
                        int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                        for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vnegq_f32(vld1q_f32(&a[i]))); }
#endif
                        for (; i < N; ++i) { a[i] = -a[i]; }
                    } else {
                        Operon::Scalar* b = slot(--sp); a = slot(sp - 1); int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                        for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vsubq_f32(vld1q_f32(&b[i]), vld1q_f32(&a[i]))); }
#endif
                        for (; i < N; ++i) { a[i] = b[i] - a[i]; }
                    }
                    scale(a, inst.weight); break; }
                case BytecodeOp::Mul: { Operon::Scalar* b = slot(--sp); Operon::Scalar* a = slot(sp - 1);
                    int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                    for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vmulq_f32(vld1q_f32(&a[i]), vld1q_f32(&b[i]))); }
#endif
                    for (; i < N; ++i) { a[i] *= b[i]; } scale(a, inst.weight); break; }
                case BytecodeOp::Div: { Operon::Scalar* a = slot(sp - 1);
                    if (inst.arity == 1) {
                        int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                        float32x4_t one = vdupq_n_f32(1.0f);
                        for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vdivq_f32(one, vld1q_f32(&a[i]))); }
#endif
                        for (; i < N; ++i) { a[i] = 1.0f / a[i]; }
                    } else {
                        Operon::Scalar* b = slot(--sp); a = slot(sp - 1); int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                        for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vdivq_f32(vld1q_f32(&b[i]), vld1q_f32(&a[i]))); }
#endif
                        for (; i < N; ++i) { a[i] = b[i] / a[i]; }
                    }
                    scale(a, inst.weight); break; }
                case BytecodeOp::Square: { Operon::Scalar* a = slot(sp - 1); int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                    for (; i + 4 <= N; i += 4) { float32x4_t v = vld1q_f32(&a[i]); vst1q_f32(&a[i], vmulq_f32(v, v)); }
#endif
                    for (; i < N; ++i) { a[i] = a[i] * a[i]; } scale(a, inst.weight); break; }
                case BytecodeOp::Sqrt: { Operon::Scalar* a = slot(sp - 1); int i = 0;
#if defined(__ARM_NEON) || defined(__aarch64__)
                    for (; i + 4 <= N; i += 4) { vst1q_f32(&a[i], vsqrtq_f32(vld1q_f32(&a[i]))); }
#endif
                    for (; i < N; ++i) { a[i] = std::sqrt(a[i]); } scale(a, inst.weight); break; }
                case BytecodeOp::Sin: case BytecodeOp::Cos: case BytecodeOp::Exp: case BytecodeOp::Log: {
                    Operon::Scalar* a = slot(sp - 1); transc(a, inst.op); scale(a, inst.weight); break; }
            }
        }
        if (sp > 0) { Operon::Scalar* r = slot(0); for (int i = 0; i < N; ++i) { out[i] = r[i]; } }
    }
};

inline auto computeAffineCoeffMask(BytecodeProgram const& prog, std::size_t nCoef) -> std::uint64_t
{
    if (!prog.valid || nCoef == 0 || prog.insts.empty()) { return 0; }
    std::vector<int> coefOfInst(prog.insts.size(), -1);
    for (std::size_t bi = 0; bi < prog.bindings.size(); ++bi) {
        auto const& b = prog.bindings[bi];
        if (b.instIndex < coefOfInst.size()) { coefOfInst[b.instIndex] = static_cast<int>(bi); }
    }
    auto bitOf = [](std::size_t k) -> std::uint64_t { return (k < 64) ? (std::uint64_t{1} << k) : 0ull; };
    struct Slot { std::uint64_t coeff{0}, aff{0}; };
    std::vector<Slot> stack;
    stack.reserve(prog.insts.size());
    for (std::size_t i = 0; i < prog.insts.size(); ++i) {
        auto const& inst = prog.insts[i];
        switch (inst.op) {
            case BytecodeOp::Constant:
            case BytecodeOp::Variable: {
                int const c = coefOfInst[i];
                Slot s;
                if (c >= 0) { s.coeff = bitOf(static_cast<std::size_t>(c)); s.aff = s.coeff; }
                stack.push_back(s);
                break;
            }
            case BytecodeOp::Add: {
                if (stack.size() < 2) { return 0; }
                Slot b = stack.back(); stack.pop_back();
                Slot a = stack.back(); stack.pop_back();
                stack.push_back({a.coeff | b.coeff, a.aff | b.aff});
                break;
            }
            case BytecodeOp::Sub: {
                if (inst.arity == 1) { break; }
                if (stack.size() < 2) { return 0; }
                Slot b = stack.back(); stack.pop_back();
                Slot a = stack.back(); stack.pop_back();
                stack.push_back({a.coeff | b.coeff, a.aff | b.aff});
                break;
            }
            case BytecodeOp::Mul: {
                if (stack.size() < 2) { return 0; }
                Slot b = stack.back(); stack.pop_back();
                Slot a = stack.back(); stack.pop_back();
                std::uint64_t const coeff = a.coeff | b.coeff;
                std::uint64_t const aff = (a.aff != 0 && b.aff != 0) ? 0ull : (a.aff | b.aff);
                stack.push_back({coeff, aff});
                break;
            }
            case BytecodeOp::Div: {
                if (inst.arity == 1) {
                    if (stack.empty()) { return 0; }
                    stack.back().aff = 0;
                    break;
                }
                if (stack.size() < 2) { return 0; }
                Slot b = stack.back(); stack.pop_back();
                Slot a = stack.back(); stack.pop_back();
                std::uint64_t const coeff = a.coeff | b.coeff;
                std::uint64_t const aff = (b.coeff != 0) ? 0ull : a.aff;
                stack.push_back({coeff, aff});
                break;
            }
            case BytecodeOp::Sin: case BytecodeOp::Cos: case BytecodeOp::Exp:
            case BytecodeOp::Log: case BytecodeOp::Sqrt: case BytecodeOp::Square: {
                if (stack.empty()) { return 0; }
                stack.back().aff = 0;
                break;
            }
        }
    }
    if (stack.empty()) { return 0; }
    return stack.back().aff;
}

}
#endif
