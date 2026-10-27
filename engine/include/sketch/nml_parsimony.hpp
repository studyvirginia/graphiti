#pragma once


#include <Eigen/Dense>
#include <operon/operon/core/tree.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace Sketch::nml {

inline constexpr int kLeafSymbols = 2;
inline constexpr int kUnarySymbols = 6;
inline constexpr int kBinarySymbols = 4;
inline constexpr double kLn2 = 0.6931471805599453;
inline constexpr double kPi = 3.141592653589793;

inline auto logTreeCounts(int cmax, int leaves = kLeafSymbols, int unary = kUnarySymbols,
                          int binary = kBinarySymbols) -> std::vector<double>
{
    double const ninf = -std::numeric_limits<double>::infinity();
    std::vector<double> lt(static_cast<std::size_t>(std::max(cmax, 1)) + 1, ninf);
    if (cmax < 1 || leaves <= 0) { return lt; }
    lt[1] = std::log(static_cast<double>(leaves));
    double const lu = unary > 0 ? std::log(static_cast<double>(unary)) : ninf;
    double const lb = binary > 0 ? std::log(static_cast<double>(binary)) : ninf;
    auto lse = [](double a, double b) {
        if (a == -std::numeric_limits<double>::infinity()) { return b; }
        if (b == -std::numeric_limits<double>::infinity()) { return a; }
        double const m = std::max(a, b);
        return m + std::log1p(std::exp(std::min(a, b) - m));
    };
    for (int c = 2; c <= cmax; ++c) {
        double acc = lu + lt[static_cast<std::size_t>(c - 1)];
        double conv = ninf;
        for (int a = 1; a <= c - 2; ++a) {
            conv = lse(conv, lt[static_cast<std::size_t>(a)] + lt[static_cast<std::size_t>(c - 1 - a)]);
        }
        acc = lse(acc, lb + conv);
        lt[static_cast<std::size_t>(c)] = acc;
    }
    return lt;
}

inline auto treeCountsExactU64(int cmax, int leaves = kLeafSymbols, int unary = kUnarySymbols,
                               int binary = kBinarySymbols) -> std::vector<std::uint64_t>
{
    std::vector<std::uint64_t> t(static_cast<std::size_t>(std::max(cmax, 1)) + 1, 0);
    if (cmax < 1) { return t; }
    t[1] = static_cast<std::uint64_t>(leaves);
    for (int c = 2; c <= cmax; ++c) {
        std::uint64_t acc = static_cast<std::uint64_t>(unary) * t[static_cast<std::size_t>(c - 1)];
        for (int a = 1; a <= c - 2; ++a) {
            acc += static_cast<std::uint64_t>(binary) * t[static_cast<std::size_t>(a)] * t[static_cast<std::size_t>(c - 1 - a)];
        }
        t[static_cast<std::size_t>(c)] = acc;
    }
    return t;
}

inline auto logTreeCount(int c) -> double
{
    static std::vector<double> const table = logTreeCounts(512);
    if (c < 1) { return -std::numeric_limits<double>::infinity(); }
    if (static_cast<std::size_t>(c) < table.size()) { return table[static_cast<std::size_t>(c)]; }
    return logTreeCounts(c)[static_cast<std::size_t>(c)];
}

inline auto eliasGammaBits(int c) -> double
{
    if (c < 1) { return std::numeric_limits<double>::infinity(); }
    int fl = 0;
    for (unsigned v = static_cast<unsigned>(c); v > 1U; v >>= 1U) { ++fl; }
    return 2.0 * static_cast<double>(fl) + 1.0;
}

inline auto structCodeNats(int c) -> double
{
    return logTreeCount(c) + eliasGammaBits(c) * kLn2;
}

struct CanonOptions {
    bool commutativeAddMul{true};
    bool subViaNegation{true};
    bool subAlwaysRedundant{false};
    bool useFlags{true};
};

inline constexpr int kFlagN = 1;
inline constexpr int kFlagS = 2;

namespace detail {
struct ExactRing {
    using T = unsigned __int128;
    static auto zero() -> T { return 0; }
    static auto fromInt(long long k) -> T { return static_cast<T>(k); }
    static auto add(T a, T b) -> T { return a + b; }
    static auto mul(T a, T b) -> T { return a * b; }
    static auto pairSame(T a) -> T { return a * (a + 1) / 2; }
};
struct LogRing {
    using T = double;
    static auto zero() -> T { return -std::numeric_limits<double>::infinity(); }
    static auto fromInt(long long k) -> T { return k > 0 ? std::log(static_cast<double>(k)) : zero(); }
    static auto add(T a, T b) -> T
    {
        if (a == zero()) { return b; }
        if (b == zero()) { return a; }
        double const m = std::max(a, b);
        return m + std::log1p(std::exp(std::min(a, b) - m));
    }
    static auto mul(T a, T b) -> T { return (a == zero() || b == zero()) ? zero() : a + b; }
    static auto pairSame(T a) -> T { return a == zero() ? zero() : a + add(a, 0.0) - kLn2; }
};
}

template <class Ring>
auto canonCountsByFlag(int cmax, CanonOptions const& o) -> std::vector<std::array<typename Ring::T, 4>>
{
    using T = typename Ring::T;
    std::array<T, 4> zeros{};
    zeros.fill(Ring::zero());
    std::vector<std::array<T, 4>> n(static_cast<std::size_t>(std::max(cmax, 1)) + 1, zeros);
    if (cmax < 1) { return n; }
    auto fl = [&](int f) { return o.useFlags ? f : 0; };
    auto acc = [&](std::array<T, 4>& dst, int f, T v) { dst[static_cast<std::size_t>(fl(f))] = Ring::add(dst[static_cast<std::size_t>(fl(f))], v); };
    acc(n[1], kFlagN | kFlagS, Ring::fromInt(1));
    acc(n[1], kFlagN, Ring::fromInt(1));
    for (int c = 2; c <= cmax; ++c) {
        auto& out = n[static_cast<std::size_t>(c)];
        for (int g = 0; g < 4; ++g) {
            T const v = n[static_cast<std::size_t>(c - 1)][static_cast<std::size_t>(g)];
            if (v == Ring::zero()) { continue; }
            bool const gN = (g & kFlagN) != 0, gS = (g & kFlagS) != 0;
            acc(out, (gN || gS) ? kFlagN : 0, v);
            acc(out, gS ? kFlagN : 0, v);
            for (int u = 0; u < 4; ++u) { acc(out, 0, v); }
        }
        for (int a = 1; a <= c - 2; ++a) {
            int const b = c - 1 - a;
            auto const& A = n[static_cast<std::size_t>(a)];
            auto const& B = n[static_cast<std::size_t>(b)];
            for (int f = 0; f < 4; ++f) {
                if (A[static_cast<std::size_t>(f)] == Ring::zero()) { continue; }
                for (int g = 0; g < 4; ++g) {
                    if (B[static_cast<std::size_t>(g)] == Ring::zero()) { continue; }
                    bool const fN = (f & kFlagN) != 0, fS = (f & kFlagS) != 0;
                    bool const gN = (g & kFlagN) != 0, gS = (g & kFlagS) != 0;
                    T const ordered = Ring::mul(A[static_cast<std::size_t>(f)], B[static_cast<std::size_t>(g)]);
                    int const addF = ((fN && gN) ? kFlagN : 0) | ((fS || gS) ? kFlagS : 0);
                    int const mulF = (fN || gN) ? kFlagN : 0;
                    acc(out, mulF, ordered);
                    if (!o.subAlwaysRedundant && (!(o.subViaNegation && o.useFlags) || !gN)) {
                        acc(out, ((fN && gN) ? kFlagN : 0) | ((fS || gS) ? kFlagS : 0), ordered);
                    }
                    if (!o.commutativeAddMul) {
                        acc(out, addF, ordered);
                        acc(out, mulF, ordered);
                    } else if (a < b) {
                        acc(out, addF, ordered);
                        acc(out, mulF, ordered);
                    } else if (a == b) {
                        if (f < g) {
                            acc(out, addF, ordered);
                            acc(out, mulF, ordered);
                        } else if (f == g) {
                            T const pr = Ring::pairSame(A[static_cast<std::size_t>(f)]);
                            acc(out, addF, pr);
                            acc(out, mulF, pr);
                        }
                    }
                }
            }
        }
    }
    return n;
}

inline auto logCanonTreeCounts(int cmax, CanonOptions const& o = {}) -> std::vector<double>
{
    auto const byF = canonCountsByFlag<detail::LogRing>(cmax, o);
    std::vector<double> out(byF.size(), detail::LogRing::zero());
    for (std::size_t c = 0; c < byF.size(); ++c) {
        for (double v : byF[c]) { out[c] = detail::LogRing::add(out[c], v); }
    }
    return out;
}

inline auto logCanonTreeCount(int c) -> double
{
    static std::vector<double> const table = logCanonTreeCounts(128);
    if (c < 1) { return -std::numeric_limits<double>::infinity(); }
    if (static_cast<std::size_t>(c) < table.size()) { return table[static_cast<std::size_t>(c)]; }
    return logCanonTreeCounts(c)[static_cast<std::size_t>(c)];
}

inline auto structCodeCanonNats(int c) -> double
{
    return logCanonTreeCount(c) + eliasGammaBits(c) * kLn2;
}

inline auto foldToGrid(double w, double period) -> double
{
    if (!(period > 0.0) || !std::isfinite(w)) { return w; }
    return w - period * std::round(w / period);
}

struct LinGroup {
    bool trig{false};
    double lin{0.0};
    std::vector<int> coefs;
    int foldCoef{-1};
    double foldSign{1.0};
};
inline auto linearGroups(Operon::Tree const& tree) -> std::vector<LinGroup>
{
    auto const& nodes = tree.Nodes();
    std::vector<LinGroup> groups;
    if (nodes.empty()) { return groups; }
    std::vector<int> coefIdx(nodes.size(), -1);
    std::vector<char> hasVar(nodes.size(), 0);
    {
        int j = 0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].Optimize) { coefIdx[i] = j++; }
            char h = nodes[i].IsVariable() ? char{1} : char{0};
            if (nodes[i].Arity > 0) { for (auto c : tree.Indices(i)) { if (hasVar[c] != 0) { h = 1; } } }
            hasVar[i] = h;
        }
    }
    auto kids = [&](std::size_t i) {
        std::vector<std::size_t> ch;
        for (auto c : tree.Indices(i)) { ch.push_back(c); }
        return ch;
    };
    struct Aff { bool ok{false}; double a{0.0}; double b{0.0}; std::vector<int> coefs; };
    auto collectCoefs = [&](std::size_t i, std::vector<int>& out) {
        auto const len = static_cast<std::size_t>(nodes[i].Length);
        for (std::size_t q = i - len; q <= i; ++q) { if (coefIdx[q] >= 0) { out.push_back(coefIdx[q]); } }
    };
    std::function<Aff(std::size_t)> affine = [&](std::size_t i) -> Aff {
        Aff r;
        auto const& n = nodes[i];
        if (n.IsConstant()) { r.ok = true; r.b = static_cast<double>(n.Value); return r; }
        if (n.IsVariable()) { r.ok = true; r.a = static_cast<double>(n.Value); if (coefIdx[i] >= 0) { r.coefs.push_back(coefIdx[i]); } return r; }
        if (n.Arity != 2) { return r; }
        auto const ch = kids(i);
        if (n.IsAddition() || n.IsSubtraction()) {
            auto const l = affine(ch[0]);
            auto const q = affine(ch[1]);
            if (!l.ok || !q.ok) { return r; }
            double const s = n.IsSubtraction() ? -1.0 : 1.0;
            r.ok = true; r.a = l.a + s * q.a; r.b = l.b + s * q.b;
            r.coefs = l.coefs; r.coefs.insert(r.coefs.end(), q.coefs.begin(), q.coefs.end());
            return r;
        }
        if (n.IsMultiplication() || n.IsDivision()) {
            bool const lv = hasVar[ch[0]] != 0, qv = hasVar[ch[1]] != 0;
            if (lv && qv) { return r; }
            if (n.IsDivision() && qv) { return r; }
            auto const l = affine(ch[0]);
            auto const q = affine(ch[1]);
            if (!l.ok || !q.ok) { return r; }
            r.ok = true;
            if (n.IsMultiplication()) {
                r.a = lv ? l.a * q.b : l.b * q.a;
                r.b = l.b * q.b;
            } else {
                r.a = l.a / q.b;
                r.b = l.b / q.b;
            }
            if (lv || qv) {
                collectCoefs(ch[0], r.coefs);
                collectCoefs(ch[1], r.coefs);
            }
            return r;
        }
        return r;
    };
    for (std::size_t p = 0; p < nodes.size(); ++p) {
        auto const& n = nodes[p];
        bool const trig = n.IsSin() || n.IsCos();
        if (!trig && !n.IsExp()) { continue; }
        LinGroup g;
        g.trig = trig;
        double bestW = -1.0;
        std::function<void(std::size_t, double)> walk = [&](std::size_t i, double sign) {
            auto const& m = nodes[i];
            if (m.IsAddition() || m.IsSubtraction()) {
                auto const ch = kids(i);
                walk(ch[0], sign);
                walk(ch[1], m.IsSubtraction() ? -sign : sign);
                return;
            }
            auto const af = affine(i);
            if (!af.ok || hasVar[i] == 0) { return; }
            g.lin += sign * af.a;
            g.coefs.insert(g.coefs.end(), af.coefs.begin(), af.coefs.end());
            if (m.IsVariable() && coefIdx[i] >= 0 && std::abs(static_cast<double>(m.Value)) > bestW) {
                bestW = std::abs(static_cast<double>(m.Value));
                g.foldCoef = coefIdx[i];
                g.foldSign = sign;
            }
        };
        walk(p - 1, 1.0);
        std::sort(g.coefs.begin(), g.coefs.end());
        g.coefs.erase(std::unique(g.coefs.begin(), g.coefs.end()), g.coefs.end());
        if (!g.coefs.empty()) { groups.push_back(std::move(g)); }
    }
    return groups;
}

struct IdentLogDet {
    double logDet{0.0};
    int rank{0};
    int nDead{0};
    int nCollinear{0};
    std::vector<int> kept;
    bool ok{false};
};

inline auto identifiableLogDet(Eigen::MatrixXd const& J, double eps) -> IdentLogDet
{
    IdentLogDet out;
    auto const n = static_cast<int>(J.rows());
    auto const k = static_cast<int>(J.cols());
    if (k == 0) { out.ok = true; return out; }
    if (n == 0 || !J.allFinite()) { return out; }
    double const tol = static_cast<double>(std::max(n, k)) * eps;
    Eigen::VectorXd norms = J.colwise().norm().transpose();
    double const maxNorm = norms.maxCoeff();
    if (!(maxNorm > 0.0) || !std::isfinite(maxNorm)) { out.nDead = k; out.ok = true; return out; }
    std::vector<int> live;
    for (int j = 0; j < k; ++j) {
        if (norms(j) > tol * maxNorm) { live.push_back(j); } else { ++out.nDead; }
    }
    if (live.empty()) { out.ok = true; return out; }
    Eigen::MatrixXd Js(n, static_cast<Eigen::Index>(live.size()));
    for (std::size_t c = 0; c < live.size(); ++c) {
        Js.col(static_cast<Eigen::Index>(c)) = J.col(live[c]) / norms(live[c]);
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Js);
    qr.setThreshold(tol);
    int const r = static_cast<int>(qr.rank());
    auto const& R = qr.matrixR();
    auto const& perm = qr.colsPermutation().indices();
    double ld = 0.0;
    for (int i = 0; i < r; ++i) {
        ld += 2.0 * std::log(std::abs(R(i, i)));
        int const col = live[static_cast<std::size_t>(perm(i))];
        ld += 2.0 * std::log(norms(col));
        out.kept.push_back(col);
    }
    std::sort(out.kept.begin(), out.kept.end());
    out.rank = r;
    out.nCollinear = static_cast<int>(live.size()) - r;
    out.logDet = ld;
    out.ok = std::isfinite(ld);
    return out;
}

enum class ConstRole : int {
    Extent = 0,
    Frequency = 1,
    Phase = 2,
    ExpRate = 3,
    Unclassified = 4,
};

struct Terms {
    double nll{std::numeric_limits<double>::quiet_NaN()};
    double paramTerm{std::numeric_limits<double>::quiet_NaN()};
    double fisherTerm{std::numeric_limits<double>::quiet_NaN()};
    double rangeTerm{std::numeric_limits<double>::quiet_NaN()};
    double structTerm{std::numeric_limits<double>::quiet_NaN()};
    double total{std::numeric_limits<double>::quiet_NaN()};
    int nodesX{0}, nodesY{0};
    int k{0};
    int kIdentifiable{0};
    int nDead{0}, nCollinear{0};
    int n{0};
    int nFreq{0}, nPhase{0}, nExpRate{0}, nExtent{0}, nUnclassified{0};
    int nOutOfRange{0};
    double maxAbsConst{0.0};
    double maxAbsFreq{0.0};
    bool ok{false};

    double structOrdered{std::numeric_limits<double>::quiet_NaN()};
    double structCanon{std::numeric_limits<double>::quiet_NaN()};
    bool admissible{true};
    int nOutOfRangeId{0};
    int nFreqOut{0};
    int nExpOut{0};
    int nIndivOut{0};
    double maxAbsLinFreq{0.0};
    double maxAbsLinExp{0.0};
    int nExtentOverE{0};
    int nExtentOverEAll{0};
    int nPhaseOutAll{0};
    int nAliasInRange{0};
    int nAliasOutOfRange{0};
    int nAliasNotFoldable{0};
    double aliasNll{std::numeric_limits<double>::quiet_NaN()};
    double aliasGridMaxDiff{std::numeric_limits<double>::quiet_NaN()};
    bool polished{false};
    bool polishRolledBack{false};
    double polishMs{0.0};
    double shapeErrBefore{std::numeric_limits<double>::quiet_NaN()};
    double shapeErrAfter{std::numeric_limits<double>::quiet_NaN()};
    double nllPolishNoRollback{std::numeric_limits<double>::quiet_NaN()};

    double absNll{std::numeric_limits<double>::quiet_NaN()};
    double absSigmaBar{std::numeric_limits<double>::quiet_NaN()};
    double absLambda{std::numeric_limits<double>::quiet_NaN()};
    double absRangeTerm{std::numeric_limits<double>::quiet_NaN()};
    double absFisherTerm{std::numeric_limits<double>::quiet_NaN()};
    double absMaxWeight{std::numeric_limits<double>::quiet_NaN()};
    double absTotal{std::numeric_limits<double>::quiet_NaN()};
    bool absOk{false};
};

inline void assemble(Terms& t, double logDetGram, double sumLogRangeKept, bool useCanonStruct = false,
                     bool outOfRangeInadmissible = false)
{
    double const n = static_cast<double>(t.n);
    double const kId = static_cast<double>(t.kIdentifiable);
    t.paramTerm = 0.5 * kId * std::log(n / (2.0 * kPi));
    t.fisherTerm = 0.5 * (logDetGram - kId * std::log(n));
    t.rangeTerm = sumLogRangeKept;
    t.structOrdered = structCodeNats(t.nodesX) + structCodeNats(t.nodesY);
    t.structCanon = structCodeCanonNats(t.nodesX) + structCodeCanonNats(t.nodesY);
    t.structTerm = useCanonStruct ? t.structCanon : t.structOrdered;
    t.total = t.nll + t.paramTerm + t.fisherTerm + t.rangeTerm + t.structTerm;
    t.ok = std::isfinite(t.total);
    if (outOfRangeInadmissible && !t.admissible && t.ok) { t.total = std::numeric_limits<double>::infinity(); }
}

}
