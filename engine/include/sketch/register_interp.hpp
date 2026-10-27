#ifndef SKETCH_REGISTER_INTERP_HPP
#define SKETCH_REGISTER_INTERP_HPP


#include "sketch/bytecode_evaluator.hpp"

#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace Sketch {

enum class RegOp : uint8_t {
    Identity,
    LoadConst,
    LoadVar,
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

struct RegInst {
    RegOp op{RegOp::Identity};
    uint8_t arity{1};
    uint16_t src1{0};
    uint16_t src2{0};
    uint16_t dst{0};
    Operon::Scalar imm{0};
    Operon::Scalar weight{1};
};

struct RegisterProgram {
    std::vector<RegInst> insts;
    int numRegs{0};
    bool valid{false};

    std::vector<BytecodeProgram::CoeffBinding> bindings;

    void UpdateCoefficients(std::span<Operon::Scalar const> coeffs)
    {
        for (std::size_t i = 0; i < bindings.size() && i < coeffs.size(); ++i) {
            auto const& b = bindings[i];
            if (b.isConstant) {
                insts[b.instIndex].imm = coeffs[i];
            } else {
                insts[b.instIndex].weight = coeffs[i];
            }
        }
    }

    static auto FromBytecode(BytecodeProgram const& bc) -> RegisterProgram
    {
        RegisterProgram rp;
        if (!bc.valid || bc.insts.empty()) { rp.valid = false; return rp; }

        std::vector<uint16_t> regStack;
        regStack.reserve(bc.insts.size());
        rp.insts.reserve(bc.insts.size());
        uint16_t nextReg = 0;

        auto mapOp = [](BytecodeOp op) -> RegOp {
            switch (op) {
                case BytecodeOp::Constant: return RegOp::LoadConst;
                case BytecodeOp::Variable: return RegOp::LoadVar;
                case BytecodeOp::Add: return RegOp::Add;
                case BytecodeOp::Sub: return RegOp::Sub;
                case BytecodeOp::Mul: return RegOp::Mul;
                case BytecodeOp::Div: return RegOp::Div;
                case BytecodeOp::Sin: return RegOp::Sin;
                case BytecodeOp::Cos: return RegOp::Cos;
                case BytecodeOp::Exp: return RegOp::Exp;
                case BytecodeOp::Log: return RegOp::Log;
                case BytecodeOp::Sqrt: return RegOp::Sqrt;
                default: return RegOp::Square;
            }
        };

        for (auto const& bi : bc.insts) {
            RegInst ri{};
            ri.op = mapOp(bi.op);
            ri.arity = bi.arity;
            ri.imm = bi.imm;
            ri.weight = bi.weight;

            bool isLeaf = (ri.op == RegOp::LoadConst || ri.op == RegOp::LoadVar);
            bool isUnaryMath = (ri.op == RegOp::Sin || ri.op == RegOp::Cos || ri.op == RegOp::Exp
                                 || ri.op == RegOp::Log || ri.op == RegOp::Sqrt || ri.op == RegOp::Square);

            if (isLeaf) {
            } else if (isUnaryMath || (bi.arity == 1)) {
                if (regStack.empty()) { rp.valid = false; return rp; }
                ri.src1 = regStack.back();
                regStack.pop_back();
                ri.arity = 1;
            } else {
                if (regStack.size() < 2) { rp.valid = false; return rp; }
                uint16_t b = regStack.back(); regStack.pop_back();
                uint16_t a = regStack.back(); regStack.pop_back();
                ri.src1 = a;
                ri.src2 = b;
                ri.arity = 2;
            }

            ri.dst = nextReg++;
            regStack.push_back(ri.dst);
            rp.insts.push_back(ri);
        }

        rp.numRegs = static_cast<int>(nextReg);
        rp.valid = (regStack.size() == 1);
        rp.bindings = bc.bindings;
        return rp;
    }

    static auto Compile(Operon::Tree const& tree, std::span<Operon::Scalar const> coeffs) -> RegisterProgram
    {
        auto bc = BytecodeProgram::Compile(tree, coeffs);
        return FromBytecode(bc);
    }
};

class RegisterInterpreter {
public:
    static void EvaluateSingle(RegisterProgram const& prog, std::span<Operon::Scalar const> t_grid,
                                std::span<Operon::Scalar> out)
    {
        std::size_t const N = t_grid.size();
        std::size_t const R = prog.insts.size();
        if (!prog.valid || R == 0) {
            for (auto& v : out) { v = Operon::Scalar{0}; }
            return;
        }
        std::vector<Operon::Scalar> regs(R);
        for (std::size_t i = 0; i < N; ++i) {
            Operon::Scalar t = t_grid[i];
            for (std::size_t k = 0; k < R; ++k) {
                regs[k] = EvalOne(prog.insts[k], regs, t);
            }
            out[i] = regs[R - 1];
        }
    }

private:
    static Operon::Scalar EvalOne(RegInst const& ri, std::vector<Operon::Scalar> const& regs, Operon::Scalar t)
    {
        switch (ri.op) {
            case RegOp::LoadConst: return ri.imm;
            case RegOp::LoadVar: return t * ri.weight;
            case RegOp::Identity: return regs[ri.src1] * ri.weight;
            case RegOp::Add: return (regs[ri.src1] + regs[ri.src2]) * ri.weight;
            case RegOp::Sub:
                return ri.arity == 1 ? (-regs[ri.src1]) * ri.weight
                                      : (regs[ri.src2] - regs[ri.src1]) * ri.weight;
            case RegOp::Mul: return (regs[ri.src1] * regs[ri.src2]) * ri.weight;
            case RegOp::Div:
                return ri.arity == 1 ? (Operon::Scalar{1} / regs[ri.src1]) * ri.weight
                                      : (regs[ri.src2] / regs[ri.src1]) * ri.weight;
            case RegOp::Sin: return std::sin(regs[ri.src1]) * ri.weight;
            case RegOp::Cos: return std::cos(regs[ri.src1]) * ri.weight;
            case RegOp::Exp: return std::exp(regs[ri.src1]) * ri.weight;
            case RegOp::Log: return std::log(regs[ri.src1]) * ri.weight;
            case RegOp::Sqrt: return std::sqrt(regs[ri.src1]) * ri.weight;
            case RegOp::Square: { auto a = regs[ri.src1]; return (a * a) * ri.weight; }
        }
        return Operon::Scalar{0};
    }
};

class BatchedRegisterInterp {
public:
    static void EvaluatePopulation(std::vector<RegisterProgram> const& progs,
                                    std::span<Operon::Scalar const> t_grid,
                                    std::vector<std::vector<Operon::Scalar>>& outPerProgram)
    {
        std::size_t const P = progs.size();
        std::size_t const M = t_grid.size();
        outPerProgram.assign(P, {});
        if (P == 0 || M == 0) { return; }

        std::size_t L = 0;
        for (auto const& pr : progs) { L = std::max(L, pr.insts.size()); }
        for (auto& o : outPerProgram) { o.assign(M, Operon::Scalar{0}); }
        if (L == 0) { return; }

        std::vector<RegInst> tape(L * P);
        std::vector<bool> progValid(P, true);

        for (std::size_t p = 0; p < P; ++p) {
            auto const& pr = progs[p];
            std::size_t const Lp = pr.insts.size();
            progValid[p] = pr.valid && Lp > 0;
            uint16_t finalDst = Lp > 0 ? pr.insts[Lp - 1].dst : 0;
            for (std::size_t k = 0; k < Lp; ++k) {
                tape[k * P + p] = pr.insts[k];
            }
            for (std::size_t k = Lp; k < L; ++k) {
                RegInst pad{};
                pad.op = RegOp::Identity;
                pad.arity = 1;
                pad.src1 = (k == Lp) ? finalDst : static_cast<uint16_t>(k - 1);
                pad.weight = Operon::Scalar{1};
                tape[k * P + p] = pad;
            }
        }

        std::vector<Operon::Scalar> regFile(L * P * M, Operon::Scalar{0});
        auto block = [&](std::size_t slot, std::size_t p) -> Operon::Scalar* {
            return regFile.data() + (slot * P + p) * M;
        };

        for (std::size_t k = 0; k < L; ++k) {
            for (std::size_t p = 0; p < P; ++p) {
                if (!progValid[p]) { continue; }
                RegInst const& ri = tape[k * P + p];
                Operon::Scalar* dst = block(k, p);
                Operon::Scalar const w = ri.weight;
                switch (ri.op) {
                    case RegOp::LoadConst: {
                        Operon::Scalar const c = ri.imm;
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = c; }
                        break;
                    }
                    case RegOp::LoadVar: {
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = t_grid[m] * w; }
                        break;
                    }
                    case RegOp::Identity: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = a[m] * w; }
                        break;
                    }
                    case RegOp::Add: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        Operon::Scalar const* b = block(ri.src2, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = (a[m] + b[m]) * w; }
                        break;
                    }
                    case RegOp::Sub: {
                        if (ri.arity == 1) {
                            Operon::Scalar const* a = block(ri.src1, p);
                            for (std::size_t m = 0; m < M; ++m) { dst[m] = (-a[m]) * w; }
                        } else {
                            Operon::Scalar const* a = block(ri.src1, p);
                            Operon::Scalar const* b = block(ri.src2, p);
                            for (std::size_t m = 0; m < M; ++m) { dst[m] = (b[m] - a[m]) * w; }
                        }
                        break;
                    }
                    case RegOp::Mul: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        Operon::Scalar const* b = block(ri.src2, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = (a[m] * b[m]) * w; }
                        break;
                    }
                    case RegOp::Div: {
                        if (ri.arity == 1) {
                            Operon::Scalar const* a = block(ri.src1, p);
                            for (std::size_t m = 0; m < M; ++m) { dst[m] = (Operon::Scalar{1} / a[m]) * w; }
                        } else {
                            Operon::Scalar const* a = block(ri.src1, p);
                            Operon::Scalar const* b = block(ri.src2, p);
                            for (std::size_t m = 0; m < M; ++m) { dst[m] = (b[m] / a[m]) * w; }
                        }
                        break;
                    }
                    case RegOp::Sin: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = std::sin(a[m]) * w; }
                        break;
                    }
                    case RegOp::Cos: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = std::cos(a[m]) * w; }
                        break;
                    }
                    case RegOp::Exp: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = std::exp(a[m]) * w; }
                        break;
                    }
                    case RegOp::Log: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = std::log(a[m]) * w; }
                        break;
                    }
                    case RegOp::Sqrt: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = std::sqrt(a[m]) * w; }
                        break;
                    }
                    case RegOp::Square: {
                        Operon::Scalar const* a = block(ri.src1, p);
                        for (std::size_t m = 0; m < M; ++m) { dst[m] = (a[m] * a[m]) * w; }
                        break;
                    }
                }
            }
        }

        for (std::size_t p = 0; p < P; ++p) {
            if (!progValid[p]) { continue; }
            Operon::Scalar const* r = block(L - 1, p);
            for (std::size_t m = 0; m < M; ++m) { outPerProgram[p][m] = r[m]; }
        }
    }
};

}
#endif
