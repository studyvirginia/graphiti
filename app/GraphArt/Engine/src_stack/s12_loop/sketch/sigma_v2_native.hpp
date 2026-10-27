#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif

namespace sketch {


struct SigmaV2Block {
    std::vector<float> c1W, c1B;
    std::vector<float> n1W, n1B;
    std::vector<float> f1gW, f1gB, f1bW, f1bB;
    std::vector<float> c2W, c2B;
    std::vector<float> n2W, n2B;
    std::vector<float> f2gW, f2gB, f2bW, f2bB;
};

struct SigmaV2Member {
    std::vector<float> h0;
    std::vector<float> invMean, invStd;
    std::vector<float> cond0W, cond0B;
    std::vector<float> cond1W, cond1B;
    std::vector<float> cond3W, cond3B;
    std::vector<float> cond4W, cond4B;
    std::vector<SigmaV2Block> blocks;
    std::vector<float> cell0W, cell0B;
    std::vector<float> cell2W, cell2B;
    std::vector<float> cell4W, cell4B;
};

struct SigmaV2Model {
    uint32_t version = 0, M = 0, H = 0, LAYERS = 0, N_INV = 0, L = 0, NC = 0, G = 0, GDIM = 0;
    std::vector<int32_t> cellLevel, cellIndex;
    std::vector<float> cellPos, cellNeg;
    std::vector<float> cellGeo;
    uint32_t nLevels = 0;
    std::vector<float> conformalQ;
    std::vector<SigmaV2Member> members;
};

namespace detail {

inline std::vector<float> readF32(std::ifstream& f, size_t n)
{
    std::vector<float> v(n);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(float)));
    if (!f) { throw std::runtime_error("sigma_v2_native: truncated file"); }
    for (float x : v) { if (!std::isfinite(x)) { throw std::runtime_error("sigma_v2_native: non-finite weight"); } }
    return v;
}

inline uint32_t readU32(std::ifstream& f)
{
    uint32_t v = 0;
    f.read(reinterpret_cast<char*>(&v), 4);
    if (!f) { throw std::runtime_error("sigma_v2_native: truncated header"); }
    return v;
}

}

inline SigmaV2Model loadSigmaV2Model(std::string const& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) { throw std::runtime_error("sigma_v2_native: cannot open " + path); }
    char magic[4];
    f.read(magic, 4);
    if (!f || std::memcmp(magic, "SGV2", 4) != 0) {
        throw std::runtime_error("sigma_v2_native: bad magic in " + path);
    }
    SigmaV2Model m;
    m.version = detail::readU32(f);
    m.M = detail::readU32(f);
    m.H = detail::readU32(f);
    m.LAYERS = detail::readU32(f);
    m.N_INV = detail::readU32(f);
    m.L = detail::readU32(f);
    m.NC = detail::readU32(f);
    m.G = detail::readU32(f);
    m.GDIM = detail::readU32(f);
    if (m.version != 1u) { throw std::runtime_error("sigma_v2_native: unsupported version " + std::to_string(m.version)); }
    if (m.M == 0u || m.H == 0u || (m.H % 2u) != 0u || m.LAYERS == 0u || m.N_INV != 28u || m.L < 2u || m.NC != m.L) {
        throw std::runtime_error("sigma_v2_native: header violates the model contract");
    }

    auto readI32 = [&](size_t n) {
        std::vector<int32_t> v(n);
        f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(int32_t)));
        if (!f) { throw std::runtime_error("sigma_v2_native: truncated geometry"); }
        return v;
    };
    m.cellLevel = readI32(m.NC);
    m.cellIndex = readI32(m.NC);
    m.cellPos = detail::readF32(f, static_cast<size_t>(m.NC) * m.L);
    m.cellNeg = detail::readF32(f, static_cast<size_t>(m.NC) * m.L);
    m.cellGeo = detail::readF32(f, static_cast<size_t>(m.NC) * m.G);
    m.nLevels = detail::readU32(f);
    m.conformalQ = detail::readF32(f, static_cast<size_t>(m.nLevels) * 3);

    size_t const H = m.H, N_INV = m.N_INV, cellIn = 2 * H + m.G + m.GDIM;
    for (uint32_t mi = 0; mi < m.M; ++mi) {
        SigmaV2Member mem;
        mem.h0 = detail::readF32(f, H);
        mem.invMean = detail::readF32(f, N_INV);
        mem.invStd = detail::readF32(f, N_INV);
        mem.cond0W = detail::readF32(f, H * N_INV); mem.cond0B = detail::readF32(f, H);
        mem.cond1W = detail::readF32(f, H); mem.cond1B = detail::readF32(f, H);
        mem.cond3W = detail::readF32(f, H * H); mem.cond3B = detail::readF32(f, H);
        mem.cond4W = detail::readF32(f, H); mem.cond4B = detail::readF32(f, H);
        mem.blocks.resize(m.LAYERS);
        for (uint32_t b = 0; b < m.LAYERS; ++b) {
            auto& bl = mem.blocks[b];
            bl.c1W = detail::readF32(f, H * H * 3); bl.c1B = detail::readF32(f, H);
            bl.n1W = detail::readF32(f, H); bl.n1B = detail::readF32(f, H);
            bl.f1gW = detail::readF32(f, H * H); bl.f1gB = detail::readF32(f, H);
            bl.f1bW = detail::readF32(f, H * H); bl.f1bB = detail::readF32(f, H);
            bl.c2W = detail::readF32(f, H * H * 3); bl.c2B = detail::readF32(f, H);
            bl.n2W = detail::readF32(f, H); bl.n2B = detail::readF32(f, H);
            bl.f2gW = detail::readF32(f, H * H); bl.f2gB = detail::readF32(f, H);
            bl.f2bW = detail::readF32(f, H * H); bl.f2bB = detail::readF32(f, H);
        }
        mem.cell0W = detail::readF32(f, H * cellIn); mem.cell0B = detail::readF32(f, H);
        mem.cell2W = detail::readF32(f, (H / 2) * H); mem.cell2B = detail::readF32(f, H / 2);
        mem.cell4W = detail::readF32(f, 1 * (H / 2)); mem.cell4B = detail::readF32(f, 1);
        m.members.push_back(std::move(mem));
    }
    if (f.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("sigma_v2_native: trailing bytes after the last member");
    }
    return m;
}

inline std::vector<std::array<double, 2>> sigmaV2ShapeSpaceRaw(
    float const* x, float const* y, int n, int L, std::array<double, 2>& outMean, double& outScale)
{
    std::vector<double> arc(static_cast<size_t>(n), 0.0);
    for (int i = 1; i < n; ++i) {
        double const dx = static_cast<double>(x[i]) - static_cast<double>(x[i - 1]);
        double const dy = static_cast<double>(y[i]) - static_cast<double>(y[i - 1]);
        arc[static_cast<size_t>(i)] = arc[static_cast<size_t>(i - 1)] + std::sqrt(dx * dx + dy * dy);
    }
    double const total = arc[static_cast<size_t>(n - 1)];
    if (n < 2 || !(total > 0.0) || !std::isfinite(total)) {
        throw std::invalid_argument("sigma_v2_native: degenerate stroke (zero or non-finite arc length)");
    }
    std::vector<double> u(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { u[static_cast<size_t>(i)] = arc[static_cast<size_t>(i)] / total; }

    std::vector<std::array<double, 2>> R(static_cast<size_t>(L));
    for (int j = 0; j < L; ++j) {
        double const s = static_cast<double>(j) / static_cast<double>(L - 1);
        int lo = 0;
        while (lo < n - 1 && u[static_cast<size_t>(lo + 1)] < s) { ++lo; }
        double px, py;
        if (s <= u[0]) { px = x[0]; py = y[0]; }
        else if (s >= u[static_cast<size_t>(n - 1)]) { px = x[n - 1]; py = y[n - 1]; }
        else {
            int hi = lo + 1;
            double const u0 = u[static_cast<size_t>(lo)], u1 = u[static_cast<size_t>(hi)];
            double const t = (u1 > u0) ? (s - u0) / (u1 - u0) : 0.0;
            px = static_cast<double>(x[lo]) + t * (static_cast<double>(x[hi]) - static_cast<double>(x[lo]));
            py = static_cast<double>(y[lo]) + t * (static_cast<double>(y[hi]) - static_cast<double>(y[lo]));
        }
        R[static_cast<size_t>(j)] = {px, py};
    }
    double mx = 0.0, my = 0.0;
    for (auto& p : R) { mx += p[0]; my += p[1]; }
    mx /= L; my /= L;
    outMean = {mx, my};
    std::vector<std::array<double, 2>> Rc(R);
    for (auto& p : Rc) { p[0] -= mx; p[1] -= my; }
    double ms = 0.0;
    for (auto& p : Rc) { ms += p[0] * p[0] + p[1] * p[1]; }
    ms /= L;
    outScale = std::sqrt(ms);
    if (!(outScale > 0.0) || !std::isfinite(outScale)) {
        throw std::invalid_argument("sigma_v2_native: degenerate stroke (zero RMS radius)");
    }
    return R;
}

inline std::vector<std::array<double, 2>> sigmaV2ShapeSpace(float const* x, float const* y, int n, int L)
{
    std::array<double, 2> m{};
    double sc = 0.0;
    std::vector<std::array<double, 2>> R = sigmaV2ShapeSpaceRaw(x, y, n, L, m, sc);
    for (auto& p : R) { p[0] = (p[0] - m[0]) / sc; p[1] = (p[1] - m[1]) / sc; }
    return R;
}

inline std::array<double, 28> sigmaV2IntrInvariants(std::vector<std::array<double, 2>> const& P)
{
    int const n = static_cast<int>(P.size());
    std::vector<float> x(static_cast<size_t>(n)), y(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { x[static_cast<size_t>(i)] = static_cast<float>(P[static_cast<size_t>(i)][0]); y[static_cast<size_t>(i)] = static_cast<float>(P[static_cast<size_t>(i)][1]); }

    std::vector<float> seg(static_cast<size_t>(n - 1)), ang(static_cast<size_t>(n - 1));
    for (int i = 0; i < n - 1; ++i) {
        float const dx = x[static_cast<size_t>(i + 1)] - x[static_cast<size_t>(i)];
        float const dy = y[static_cast<size_t>(i + 1)] - y[static_cast<size_t>(i)];
        seg[static_cast<size_t>(i)] = std::sqrt(dx * dx + dy * dy) + 1e-9f;
        ang[static_cast<size_t>(i)] = std::atan2(dy, dx);
    }
    int const nd = n - 2;
    std::vector<float> da(static_cast<size_t>(std::max(nd, 0)));
    constexpr float kPi = 3.14159265358979323846f;
    for (int i = 0; i < nd; ++i) {
        float v = ang[static_cast<size_t>(i + 1)] - ang[static_cast<size_t>(i)];
        v = std::fmod(v + kPi, 2 * kPi);
        if (v < 0) { v += 2 * kPi; }
        da[static_cast<size_t>(i)] = v - kPi;
    }
    std::vector<float> curv(static_cast<size_t>(n), 0.0f);
    for (int i = 0; i < nd; ++i) { curv[static_cast<size_t>(i + 1)] = std::max(-kPi, std::min(kPi, da[static_cast<size_t>(i)])); }
    std::vector<float> aca(static_cast<size_t>(std::max(nd, 0)));
    for (int i = 0; i < nd; ++i) { aca[static_cast<size_t>(i)] = std::fabs(da[static_cast<size_t>(i)]); }
    std::vector<float> dcurv(static_cast<size_t>(n - 1));
    for (int i = 0; i < n - 1; ++i) { dcurv[static_cast<size_t>(i)] = std::fabs(curv[static_cast<size_t>(i + 1)] - curv[static_cast<size_t>(i)]); }

    auto mean = [](std::vector<float> const& v) { float s = 0; for (float e : v) s += e; return v.empty() ? 0.0f : s / static_cast<float>(v.size()); };
    auto stddev = [&](std::vector<float> const& v) { float m = mean(v); float s = 0; for (float e : v) s += (e - m) * (e - m); return v.empty() ? 0.0f : std::sqrt(s / static_cast<float>(v.size())); };
    auto vmax = [](std::vector<float> const& v) { float b = v.empty() ? 0.0f : v[0]; for (float e : v) b = std::max(b, e); return b; };
    auto vsum = [](std::vector<float> const& v) { float s = 0; for (float e : v) s += e; return s; };
    auto percentile90 = [](std::vector<float> v) {
        if (v.empty()) { return 0.0f; }
        std::sort(v.begin(), v.end());
        float const pos = 0.90f * static_cast<float>(v.size() - 1);
        size_t const lo = static_cast<size_t>(std::floor(pos));
        size_t const hi = static_cast<size_t>(std::ceil(pos));
        float const t = pos - static_cast<float>(lo);
        return v[lo] + t * (v[hi] - v[lo]);
    };

    int inflections = 0;
    {
        std::vector<int> sgn(static_cast<size_t>(std::max(nd, 0)));
        for (int i = 0; i < nd; ++i) { sgn[static_cast<size_t>(i)] = (da[static_cast<size_t>(i)] > 0) - (da[static_cast<size_t>(i)] < 0); }
        for (int i = 0; i + 1 < nd; ++i) { if (std::abs(sgn[static_cast<size_t>(i + 1)] - sgn[static_cast<size_t>(i)]) > 0) { ++inflections; } }
    }

    float const cmean = mean(curv);
    int const nb = n / 2 + 1;
    std::vector<float> mag(static_cast<size_t>(nb));
    for (int k = 0; k < nb; ++k) {
        double re = 0.0, im = 0.0;
        for (int t = 0; t < n; ++t) {
            double const ang_ = -2.0 * 3.14159265358979323846 * static_cast<double>(k) * static_cast<double>(t) / static_cast<double>(n);
            double const v = static_cast<double>(curv[static_cast<size_t>(t)]) - static_cast<double>(cmean);
            re += v * std::cos(ang_);
            im += v * std::sin(ang_);
        }
        mag[static_cast<size_t>(k)] = static_cast<float>(std::sqrt(re * re + im * im));
    }
    float const magSum = vsum(mag) + 1e-9f;
    std::vector<float> cf(mag.size());
    for (size_t i = 0; i < mag.size(); ++i) { cf[i] = mag[i] / magSum; }
    float specEntropy = 0.0f;
    for (float c : cf) { specEntropy += -(c * std::log(c + 1e-12f)); }
    float const cf1 = cf.size() > 1 ? cf[1] : 0.0f, cf2 = cf.size() > 2 ? cf[2] : 0.0f,
                cf3 = cf.size() > 3 ? cf[3] : 0.0f, cf4v = cf.size() > 4 ? cf[4] : 0.0f;

    std::vector<float> rad(static_cast<size_t>(n));
    {
        float xm = mean(x), ym = mean(y);
        for (int i = 0; i < n; ++i) {
            float const dx = x[static_cast<size_t>(i)] - xm, dy = y[static_cast<size_t>(i)] - ym;
            rad[static_cast<size_t>(i)] = std::sqrt(dx * dx + dy * dy);
        }
    }
    float const arclen = vsum(seg);
    float const chordx = x[static_cast<size_t>(n - 1)] - x[0], chordy = y[static_cast<size_t>(n - 1)] - y[0];
    float const chord = std::sqrt(chordx * chordx + chordy * chordy);

    int const n2 = static_cast<int>(P.size());
    double xm2 = 0, ym2 = 0;
    for (auto const& p : P) { xm2 += p[0]; ym2 += p[1]; }
    xm2 /= n2; ym2 /= n2;
    double cxx = 0, cxy = 0, cyy = 0;
    for (auto const& p : P) {
        double const dx = p[0] - xm2, dy = p[1] - ym2;
        cxx += dx * dx; cxy += dx * dy; cyy += dy * dy;
    }
    cxx /= n2; cxy /= n2; cyy /= n2;
    double const tr = cxx + cyy, det = cxx * cyy - cxy * cxy;
    double const disc = std::sqrt(std::max(tr * tr / 4.0 - det, 0.0));
    double const evA = tr / 2.0 + disc, evB = tr / 2.0 - disc;
    double const l1 = std::max(std::max(evA, evB), 1e-12);
    double const l2 = std::max(std::min(evA, evB), 1e-12);
    double sa = 0.0;
    for (int i = 0; i < n2 - 1; ++i) {
        sa += P[static_cast<size_t>(i)][0] * P[static_cast<size_t>(i + 1)][1] - P[static_cast<size_t>(i + 1)][0] * P[static_cast<size_t>(i)][1];
    }
    sa = 0.5 * std::fabs(sa);

    std::array<double, 28> inv{};
    inv[0] = mean(aca);
    inv[1] = stddev(aca);
    inv[2] = vmax(aca);
    inv[3] = vsum(da);
    inv[4] = vsum(aca);
    inv[5] = static_cast<double>(inflections);
    inv[6] = std::sqrt(l2 / l1);
    inv[7] = chord;
    inv[8] = arclen;
    inv[9] = arclen / (chord + 1e-6f);
    inv[10] = std::sqrt(l1);
    inv[11] = std::sqrt(l2);
    inv[12] = mean(curv);
    inv[13] = stddev(curv);
    inv[14] = percentile90(aca);
    inv[15] = static_cast<double>(n);
    inv[16] = mean([&] { std::vector<float> c2(curv.size()); for (size_t i = 0; i < curv.size(); ++i) c2[i] = curv[i] * curv[i]; return c2; }());
    inv[17] = specEntropy;
    inv[18] = cf1; inv[19] = cf2; inv[20] = cf3; inv[21] = cf4v;
    inv[22] = mean(dcurv);
    inv[23] = vmax(dcurv);
    inv[24] = mean(rad);
    inv[25] = stddev(rad);
    inv[26] = sa / std::sqrt(l1 * l2);
    inv[27] = stddev(seg) / (mean(seg) + 1e-9f);
    for (auto& v : inv) {
        if (std::isnan(v)) { v = 0.0; }
        else if (v > 1e300) { v = 10.0; }
        else if (v < -1e300) { v = -10.0; }
    }
    return inv;
}

namespace detail {

inline float gelu(float x)
{
    return 0.5f * x * (1.0f + std::erf(x * 0.70710678118654752440f));
}

inline void linear(float const* W, float const* b, float const* x, float* y, int out, int in)
{
    for (int o = 0; o < out; ++o) {
        double acc = b[o];
        float const* row = W + static_cast<size_t>(o) * static_cast<size_t>(in);
        for (int i = 0; i < in; ++i) { acc += static_cast<double>(row[i]) * static_cast<double>(x[i]); }
        y[o] = static_cast<float>(acc);
    }
}

inline void layerNorm(float const* w, float const* b, float const* x, float* y, int n)
{
    double m = 0; for (int i = 0; i < n; ++i) { m += x[i]; } m /= n;
    double v = 0; for (int i = 0; i < n; ++i) { double d = x[i] - m; v += d * d; } v /= n;
    double const inv = 1.0 / std::sqrt(v + 1e-5);
    for (int i = 0; i < n; ++i) { y[i] = static_cast<float>((x[i] - m) * inv * w[i] + b[i]); }
}

inline void groupNorm(float const* w, float const* b, std::vector<float> const& x, std::vector<float>& y, int H, int L, int groups)
{
    y.resize(x.size());
    int const cpg = H / groups;
    for (int g = 0; g < groups; ++g) {
        double m = 0; size_t cnt = static_cast<size_t>(cpg) * static_cast<size_t>(L);
        for (int c = g * cpg; c < (g + 1) * cpg; ++c) {
            for (int l = 0; l < L; ++l) { m += x[static_cast<size_t>(c) * L + static_cast<size_t>(l)]; }
        }
        m /= static_cast<double>(cnt);
        double v = 0;
        for (int c = g * cpg; c < (g + 1) * cpg; ++c) {
            for (int l = 0; l < L; ++l) { double d = x[static_cast<size_t>(c) * L + static_cast<size_t>(l)] - m; v += d * d; }
        }
        v /= static_cast<double>(cnt);
        double const inv = 1.0 / std::sqrt(v + 1e-5);
        for (int c = g * cpg; c < (g + 1) * cpg; ++c) {
            for (int l = 0; l < L; ++l) {
                size_t const idx = static_cast<size_t>(c) * L + static_cast<size_t>(l);
                y[idx] = static_cast<float>((x[idx] - m) * inv * w[c] + b[c]);
            }
        }
    }
}

inline void conv1dSame(float const* W, float const* b, std::vector<float> const& x, std::vector<float>& y, int H, int L, int dilation)
{
    y.assign(x.size(), 0.0f);
    for (int oc = 0; oc < H; ++oc) {
        std::vector<double> acc(static_cast<size_t>(L), static_cast<double>(b[oc]));
        for (int ic = 0; ic < H; ++ic) {
            float const* wrow = W + (static_cast<size_t>(oc) * H + static_cast<size_t>(ic)) * 3;
            float const* xrow = x.data() + static_cast<size_t>(ic) * L;
            for (int l = 0; l < L; ++l) {
                double s = static_cast<double>(wrow[1]) * static_cast<double>(xrow[l]);
                int const lm = l - dilation, lp = l + dilation;
                if (lm >= 0) { s += static_cast<double>(wrow[0]) * static_cast<double>(xrow[lm]); }
                if (lp < L) { s += static_cast<double>(wrow[2]) * static_cast<double>(xrow[lp]); }
                acc[static_cast<size_t>(l)] += s;
            }
        }
        for (int l = 0; l < L; ++l) { y[static_cast<size_t>(oc) * L + static_cast<size_t>(l)] = static_cast<float>(acc[static_cast<size_t>(l)]); }
    }
}

inline auto b3Zi5AccelFromEnv() -> bool
{
    static bool const v = [] {
#if defined(__APPLE__)
        char const* e = std::getenv("SKETCH_B3_ZI5_ACCEL");
        return e == nullptr || std::atoi(e) != 0;
#else
        return false;
#endif
    }();
    return v;
}

#if defined(__APPLE__)

inline void linearAccel(float const* W, float const* b, float const* x, float* y, int out, int in)
{
    std::memcpy(y, b, static_cast<size_t>(out) * sizeof(float));
    cblas_sgemv(CblasRowMajor, CblasNoTrans, out, in, 1.0f, W, in, x, 1, 1.0f, y, 1);
}

inline void conv1dSameAccel(float const* W, float const* b, std::vector<float> const& x, std::vector<float>& y,
                             int H, int L, int dilation)
{
    y.assign(x.size(), 0.0f);
    size_t const HH = static_cast<size_t>(H) * static_cast<size_t>(H);
    size_t const HL = static_cast<size_t>(H) * static_cast<size_t>(L);
    for (int oc = 0; oc < H; ++oc) {
        for (int l = 0; l < L; ++l) { y[static_cast<size_t>(oc) * L + static_cast<size_t>(l)] = b[oc]; }
    }
    std::vector<float> Wk(HH);
    std::vector<float> Xs(HL);
    int const offsets[3] = {-dilation, 0, dilation};
    for (int k = 0; k < 3; ++k) {
        for (int oc = 0; oc < H; ++oc) {
            for (int ic = 0; ic < H; ++ic) {
                Wk[static_cast<size_t>(oc) * H + static_cast<size_t>(ic)] =
                    W[(static_cast<size_t>(oc) * H + static_cast<size_t>(ic)) * 3 + static_cast<size_t>(k)];
            }
        }
        int const off = offsets[k];
        for (int ic = 0; ic < H; ++ic) {
            float const* xrow = x.data() + static_cast<size_t>(ic) * L;
            float* xsrow = Xs.data() + static_cast<size_t>(ic) * L;
            for (int l = 0; l < L; ++l) {
                int const src = l + off;
                xsrow[l] = (src >= 0 && src < L) ? xrow[src] : 0.0f;
            }
        }
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, H, L, H, 1.0f, Wk.data(), H, Xs.data(), L,
                    1.0f, y.data(), L);
    }
}

#endif

inline void film(float const* gW, float const* gB, float const* bW, float const* bB, float const* z,
                 std::vector<float>& x, int H, int L)
{
    std::vector<float> g(static_cast<size_t>(H)), bb(static_cast<size_t>(H));
    linear(gW, gB, z, g.data(), H, H);
    linear(bW, bB, z, bb.data(), H, H);
    for (int c = 0; c < H; ++c) {
        for (int l = 0; l < L; ++l) {
            size_t const idx = static_cast<size_t>(c) * L + static_cast<size_t>(l);
            x[idx] = g[static_cast<size_t>(c)] * x[idx] + bb[static_cast<size_t>(c)];
        }
    }
}

inline std::vector<float> sigmaV2MemberForward(SigmaV2Model const& model, SigmaV2Member const& mem, float const* invStd)
{
    int const H = static_cast<int>(model.H), L = static_cast<int>(model.L), N_INV = static_cast<int>(model.N_INV);
#if defined(__APPLE__)
    bool const accel = b3Zi5AccelFromEnv();
    auto lin = accel ? &linearAccel : &linear;
    auto conv = accel ? &conv1dSameAccel : &conv1dSame;
#else
    auto lin = &linear;
    auto conv = &conv1dSame;
#endif
    std::vector<float> z0(static_cast<size_t>(H)), z0n(static_cast<size_t>(H)), z0g(static_cast<size_t>(H));
    lin(mem.cond0W.data(), mem.cond0B.data(), invStd, z0.data(), H, N_INV);
    layerNorm(mem.cond1W.data(), mem.cond1B.data(), z0.data(), z0n.data(), H);
    for (int i = 0; i < H; ++i) { z0g[static_cast<size_t>(i)] = gelu(z0n[static_cast<size_t>(i)]); }
    std::vector<float> z1(static_cast<size_t>(H)), z1n(static_cast<size_t>(H)), z(static_cast<size_t>(H));
    lin(mem.cond3W.data(), mem.cond3B.data(), z0g.data(), z1.data(), H, H);
    layerNorm(mem.cond4W.data(), mem.cond4B.data(), z1.data(), z1n.data(), H);
    for (int i = 0; i < H; ++i) { z[static_cast<size_t>(i)] = gelu(z1n[static_cast<size_t>(i)]); }

    std::vector<float> h(static_cast<size_t>(H) * static_cast<size_t>(L));
    for (int c = 0; c < H; ++c) { for (int l = 0; l < L; ++l) { h[static_cast<size_t>(c) * L + static_cast<size_t>(l)] = mem.h0[static_cast<size_t>(c)]; } }

    std::vector<float> t1, t1n, t2, t2n;
    for (uint32_t bi = 0; bi < model.LAYERS; ++bi) {
        auto const& bl = mem.blocks[bi];
        int const dilation = 1 << bi;
        conv(bl.c1W.data(), bl.c1B.data(), h, t1, H, L, dilation);
        groupNorm(bl.n1W.data(), bl.n1B.data(), t1, t1n, H, L, 8);
        film(bl.f1gW.data(), bl.f1gB.data(), bl.f1bW.data(), bl.f1bB.data(), z.data(), t1n, H, L);
        for (auto& v : t1n) { v = gelu(v); }
        conv(bl.c2W.data(), bl.c2B.data(), t1n, t2, H, L, dilation);
        groupNorm(bl.n2W.data(), bl.n2B.data(), t2, t2n, H, L, 8);
        film(bl.f2gW.data(), bl.f2gB.data(), bl.f2bW.data(), bl.f2bB.data(), z.data(), t2n, H, L);
        for (auto& v : t2n) { v = gelu(v); }
        for (size_t i = 0; i < h.size(); ++i) { h[i] = t2n[i] + h[i]; }
    }

    std::vector<float> gcat(model.GDIM);
    for (int c = 0; c < H; ++c) {
        double m = 0; for (int l = 0; l < L; ++l) { m += h[static_cast<size_t>(c) * L + static_cast<size_t>(l)]; }
        gcat[static_cast<size_t>(c)] = static_cast<float>(m / L);
    }
    for (int c = 0; c < H; ++c) { gcat[static_cast<size_t>(H + c)] = z[static_cast<size_t>(c)]; }

    int const NC = static_cast<int>(model.NC), G = static_cast<int>(model.G);
    size_t const cellIn = 2 * static_cast<size_t>(H) + static_cast<size_t>(G) + model.GDIM;
    int const H2 = H / 2;
    std::vector<float> feat(cellIn), t0(static_cast<size_t>(H)), tc(static_cast<size_t>(H2)), out1(1);
    std::vector<float> logE(static_cast<size_t>(NC));
    for (int c = 0; c < NC; ++c) {
        for (int ch = 0; ch < H; ++ch) {
            double hp = 0, hn = 0;
            for (int l = 0; l < L; ++l) {
                float const px = model.cellPos[static_cast<size_t>(c) * L + static_cast<size_t>(l)];
                float const nx = model.cellNeg[static_cast<size_t>(c) * L + static_cast<size_t>(l)];
                float const hv = h[static_cast<size_t>(ch) * L + static_cast<size_t>(l)];
                hp += static_cast<double>(px) * hv;
                hn += static_cast<double>(nx) * hv;
            }
            feat[static_cast<size_t>(ch)] = static_cast<float>(hp);
            feat[static_cast<size_t>(H + ch)] = static_cast<float>(hn);
        }
        for (int gi = 0; gi < G; ++gi) { feat[static_cast<size_t>(2 * H + gi)] = model.cellGeo[static_cast<size_t>(c) * G + static_cast<size_t>(gi)]; }
        for (uint32_t gi = 0; gi < model.GDIM; ++gi) { feat[static_cast<size_t>(2 * H) + G + gi] = gcat[gi]; }
        lin(mem.cell0W.data(), mem.cell0B.data(), feat.data(), t0.data(), H, static_cast<int>(cellIn));
        for (auto& v : t0) { v = gelu(v); }
        lin(mem.cell2W.data(), mem.cell2B.data(), t0.data(), tc.data(), H2, H);
        for (auto& v : tc) { v = gelu(v); }
        lin(mem.cell4W.data(), mem.cell4B.data(), tc.data(), out1.data(), 1, H2);
        logE[static_cast<size_t>(c)] = out1[0];
    }
    return logE;
}

}

struct SigmaV2CellOutput {
    std::vector<double> lnE;
    std::vector<double> sd;
    std::vector<int32_t> level, index;
    std::vector<std::vector<double>> memberLnE;
};

inline SigmaV2CellOutput sigmaV2Infer(SigmaV2Model const& model, float const* x, float const* y, int n)
{
    auto const P = sigmaV2ShapeSpace(x, y, n, static_cast<int>(model.L));
    auto const inv = sigmaV2IntrInvariants(P);

    int const NC = static_cast<int>(model.NC), N_INV = static_cast<int>(model.N_INV);
    std::vector<std::vector<float>> muMembers(model.M, std::vector<float>(static_cast<size_t>(NC)));
    for (uint32_t mi = 0; mi < model.M; ++mi) {
        auto const& mem = model.members[mi];
        std::vector<float> ivn(static_cast<size_t>(N_INV));
        for (int i = 0; i < N_INV; ++i) {
            double v = (inv[static_cast<size_t>(i)] - static_cast<double>(mem.invMean[static_cast<size_t>(i)])) / static_cast<double>(mem.invStd[static_cast<size_t>(i)]);
            v = std::max(-5.0, std::min(5.0, v));
            ivn[static_cast<size_t>(i)] = static_cast<float>(v);
        }
        auto const logE = detail::sigmaV2MemberForward(model, mem, ivn.data());
        for (int c = 0; c < NC; ++c) {
            double v = logE[static_cast<size_t>(c)];
            v = std::max(-30.0, std::min(10.0, v));
            muMembers[mi][static_cast<size_t>(c)] = static_cast<float>(v);
        }
    }

    SigmaV2CellOutput out;
    out.lnE.resize(static_cast<size_t>(NC));
    out.sd.resize(static_cast<size_t>(NC));
    out.level.assign(model.cellLevel.begin(), model.cellLevel.end());
    out.index.assign(model.cellIndex.begin(), model.cellIndex.end());
    out.memberLnE.assign(model.M, std::vector<double>(static_cast<size_t>(NC)));
    for (int c = 0; c < NC; ++c) {
        double sumExp = 0.0, sumMu = 0.0, sumMu2 = 0.0;
        for (uint32_t mi = 0; mi < model.M; ++mi) {
            double const mu = muMembers[mi][static_cast<size_t>(c)];
            out.memberLnE[mi][static_cast<size_t>(c)] = mu;
            sumExp += std::exp(mu);
            sumMu += mu;
            sumMu2 += mu * mu;
        }
        double const M = static_cast<double>(model.M);
        double const E = sumExp / M;
        double const meanMu = sumMu / M;
        double const varMu = std::max(0.0, sumMu2 / M - meanMu * meanMu);
        out.lnE[static_cast<size_t>(c)] = std::log(E);
        out.sd[static_cast<size_t>(c)] = std::sqrt(varMu);
    }
    return out;
}

namespace detail {
struct SigmaV2Cache {
    std::mutex m;
    bool attempted = false;
    bool ok = false;
    SigmaV2Model model;
};
inline auto sigmaV2CacheSingleton() -> SigmaV2Cache&
{
    static SigmaV2Cache cache;
    return cache;
}
}

inline auto resolveSigmaV2ModelPath() -> std::string
{
    char const* envPath = std::getenv("SKETCH_SIGMA_V2_BIN");
    if (envPath != nullptr && *envPath != '\0') { return std::string(envPath); }
    std::string const rel = "models/sigma_v2_zi5.bin";
    auto existsFile = [](std::string const& p) { std::ifstream f(p, std::ios::binary); return static_cast<bool>(f); };
    auto searchUp = [&](std::string dir) -> std::string {
        for (int depth = 0; depth < 64 && !dir.empty(); ++depth) {
            std::string const cand = dir + "/" + rel;
            if (existsFile(cand)) { return cand; }
            auto const slash = dir.find_last_of('/');
            if (slash == std::string::npos || slash == 0) { break; }
            dir = dir.substr(0, slash);
        }
        return {};
    };
    if (existsFile(rel)) { return rel; }
    if (char const* pwd = std::getenv("PWD"); pwd != nullptr) {
        if (auto r = searchUp(pwd); !r.empty()) { return r; }
    }
#if defined(__APPLE__)
    {
        char buf[4096];
        uint32_t sz = sizeof(buf);
        if (_NSGetExecutablePath(buf, &sz) == 0) {
            std::string exe(buf);
            auto const slash = exe.find_last_of('/');
            if (slash != std::string::npos) {
                if (auto r = searchUp(exe.substr(0, slash)); !r.empty()) { return r; }
            }
        }
    }
#endif
    return rel;
}

inline auto loadSigmaV2ModelCached() -> SigmaV2Model const*
{
    auto& cache = detail::sigmaV2CacheSingleton();
    std::lock_guard<std::mutex> lock(cache.m);
    if (!cache.attempted) {
        cache.attempted = true;
        std::string const path = resolveSigmaV2ModelPath();
        try {
            cache.model = loadSigmaV2Model(path);
            std::fprintf(stderr, "[sigma-net] ZI5 model loaded from '%s'\n", path.c_str());
            cache.ok = true;
        } catch (std::exception const& e) {
            cache.ok = false;
            std::fprintf(stderr, "[sigma-net] FAILED to load ZI5 model '%s': %s -- served sigma is DISABLED (set SKETCH_SIGMA_V2_BIN)\n", path.c_str(), e.what());
        }
    }
    return cache.ok ? &cache.model : nullptr;
}

}
