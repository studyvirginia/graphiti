#pragma once


#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

#include "frequency_discovery.hpp"

namespace Sketch {
namespace freqdisc {

constexpr double kPi = 3.141592653589793;

struct PencilPole {
    double omega{0.0};
    double damping{0.0};
    std::complex<double> amplitude{};
    double energy{0.0};
    double energyRe{0.0};
    double energyIm{0.0};
};

struct PencilResult {
    bool valid{false};
    std::size_t n{0};
    std::size_t L{0};
    std::size_t M{0};
    double dt{0.0};
    std::vector<PencilPole> poles;
    int seed{-1};
    std::vector<double> singularValues;

    [[nodiscard]] auto seedProposal(bool xShare) const -> FrequencyProposal
    {
        FrequencyProposal p;
        if (!valid || seed < 0) { return p; }
        auto const& pole = poles[static_cast<std::size_t>(seed)];
        p.valid = true;
        p.omega = std::fabs(pole.omega);
        p.magnitude = xShare ? pole.energyRe : pole.energyIm;
        p.damping = pole.damping;
        return p;
    }
};

inline auto pencilParameter(std::size_t n) -> std::size_t { return n / 3; }

inline auto mdlOrder(std::vector<double> const& sv, std::size_t Ns, std::size_t p, bool gram = false)
    -> std::size_t
{
    if (sv.empty() || Ns == 0 || p < 2) { return 0; }
    std::size_t const q = std::min(sv.size(), p);
    double const eps = std::numeric_limits<double>::epsilon();
    double const tol = gram ? sv[0] * std::sqrt(static_cast<double>(p) * eps)
                            : static_cast<double>(std::max(Ns, p)) * eps * sv[0];
    std::vector<double> lam(p, 0.0);
    for (std::size_t i = 0; i < p; ++i) {
        double const s = (i < q) ? std::max(sv[i], tol) : tol;
        lam[i] = s * s / static_cast<double>(Ns);
    }
    if (!(lam[0] > 0.0) || !std::isfinite(lam[0])) { return 0; }
    double const lnNs = std::log(static_cast<double>(Ns));
    std::size_t bestK = 0;
    double bestV = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < p; ++k) {
        std::size_t const r = p - k;
        double sumLog = 0.0, sum = 0.0;
        for (std::size_t i = k; i < p; ++i) { sumLog += std::log(lam[i]); sum += lam[i]; }
        double const lnGeo = sumLog / static_cast<double>(r);
        double const lnArith = std::log(sum / static_cast<double>(r));
        double const kd = static_cast<double>(k);
        double const v = -static_cast<double>(Ns) * static_cast<double>(r) * (lnGeo - lnArith)
                       + 0.5 * kd * (2.0 * static_cast<double>(p) - kd) * lnNs;
        if (v < bestV) { bestV = v; bestK = k; }
    }
    return bestK;
}

inline auto matrixPencil(std::vector<std::complex<double>> const& y, double span) -> PencilResult
{
    using CD = std::complex<double>;
    PencilResult out;
    std::size_t const n = y.size();
    out.n = n;
    if (n < 8 || !(span > kEps) || !std::isfinite(span)) { return out; }
    for (auto const& v : y) { if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) { return out; } }

    double const dt = span / static_cast<double>(n - 1);
    std::size_t const L = pencilParameter(n);
    std::size_t const rows = n - L, cols = L + 1;
    out.dt = dt; out.L = L;

    Eigen::MatrixXcd Y(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            Y(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = y[i + j];
        }
    }
    Eigen::MatrixXcd const G = Y.adjoint() * Y;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es0(G);
    if (es0.info() != Eigen::Success) { return out; }
    Eigen::Index const pc = static_cast<Eigen::Index>(cols);
    out.singularValues.resize(cols);
    for (Eigen::Index i = 0; i < pc; ++i) {
        out.singularValues[static_cast<std::size_t>(i)] = std::sqrt(std::max(0.0, es0.eigenvalues()(pc - 1 - i)));
    }
    if (out.singularValues.empty() || !(out.singularValues[0] > 0.0)) { return out; }

    std::size_t const M = mdlOrder(out.singularValues, rows, cols,true);
    out.M = M;
    if (M == 0) { return out; }

    Eigen::MatrixXcd W(pc, static_cast<Eigen::Index>(M));
    for (std::size_t m = 0; m < M; ++m) {
        W.col(static_cast<Eigen::Index>(m)) = es0.eigenvectors().col(pc - 1 - static_cast<Eigen::Index>(m)).conjugate();
    }
    Eigen::MatrixXcd const W1 = W.topRows(static_cast<Eigen::Index>(L));
    Eigen::MatrixXcd const W2 = W.bottomRows(static_cast<Eigen::Index>(L));
    Eigen::MatrixXcd const F = W1.colPivHouseholderQr().solve(W2);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(F,false);
    if (es.info() != Eigen::Success) { return out; }
    auto const& z = es.eigenvalues();

    Eigen::MatrixXcd V(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(M));
    for (std::size_t m = 0; m < M; ++m) {
        CD const zm = z(static_cast<Eigen::Index>(m));
        if (!std::isfinite(zm.real()) || !std::isfinite(zm.imag()) || std::abs(zm) == 0.0) { return out; }
        CD p{1.0, 0.0};
        for (std::size_t k = 0; k < n; ++k) {
            V(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(m)) = p;
            p *= zm;
        }
    }
    Eigen::VectorXcd yv(static_cast<Eigen::Index>(n));
    for (std::size_t k = 0; k < n; ++k) { yv(static_cast<Eigen::Index>(k)) = y[k]; }
    Eigen::VectorXcd const a = V.colPivHouseholderQr().solve(yv);

    out.poles.resize(M);
    for (std::size_t m = 0; m < M; ++m) {
        CD const zm = z(static_cast<Eigen::Index>(m));
        CD const s = std::log(zm) / dt;
        PencilPole pole;
        pole.omega = s.imag();
        pole.damping = s.real();
        pole.amplitude = a(static_cast<Eigen::Index>(m));
        double eRe = 0.0, eIm = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            CD const c = pole.amplitude * V(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(m));
            eRe += c.real() * c.real();
            eIm += c.imag() * c.imag();
        }
        pole.energyRe = eRe; pole.energyIm = eIm; pole.energy = eRe + eIm;
        if (!std::isfinite(pole.omega) || !std::isfinite(pole.damping) || !std::isfinite(pole.energy)) {
            return out;
        }
        out.poles[m] = pole;
    }
    std::stable_sort(out.poles.begin(), out.poles.end(),
                     [](PencilPole const& p1, PencilPole const& p2) { return p1.energy > p2.energy; });
    out.valid = true;
    double const minOmega = kPi / span;
    for (std::size_t m = 0; m < out.poles.size(); ++m) {
        if (std::fabs(out.poles[m].omega) >= minOmega && out.poles[m].energy > kEps) {
            out.seed = static_cast<int>(m);
            break;
        }
    }
    return out;
}

struct PencilSeedEvent {
    std::size_t n{0}, L{0}, M{0};
    bool pencilValid{false};
    double pencilOmega{0.0}, pencilDamping{0.0};
    bool pencilUseX{true};
    bool fftValidX{false}, fftValidY{false};
    double fftOmegaX{0.0}, fftOmegaY{0.0}, fftMagX{0.0}, fftMagY{0.0};
    double pencilNs{0.0}, fftNs{0.0};
    double sv0{0.0}, svLast{0.0};
    std::vector<double> const* resX{nullptr};
    std::vector<double> const* resY{nullptr};
};
using PencilSeedTraceFn = void (*)(PencilSeedEvent const&);
inline std::atomic<PencilSeedTraceFn> g_pencilSeedTrace{nullptr};

inline void proposeFrequencyPencilXY(std::vector<double> const& resX, std::vector<double> const& resY,
                                     double span, FrequencyProposal& propX, FrequencyProposal& propY)
{
    auto const trace = g_pencilSeedTrace.load(std::memory_order_relaxed);
    auto const t0 = std::chrono::steady_clock::now();
    std::vector<std::complex<double>> z(resX.size());
    bool const hasY = resY.size() == resX.size();
    for (std::size_t i = 0; i < resX.size(); ++i) { z[i] = {resX[i], hasY ? resY[i] : 0.0}; }
    auto const r = matrixPencil(z, span);
    propX = r.seedProposal(true);
    propY = r.seedProposal(hasY ? false : true);
    if (trace == nullptr) { return; }
    auto const t1 = std::chrono::steady_clock::now();
    PencilSeedEvent ev;
    ev.n = r.n; ev.L = r.L; ev.M = r.M;
    ev.resX = &resX; ev.resY = &resY;
    if (!r.singularValues.empty()) { ev.sv0 = r.singularValues.front(); ev.svLast = r.singularValues.back(); }
    ev.pencilValid = propX.valid;
    ev.pencilOmega = propX.omega; ev.pencilDamping = propX.damping;
    ev.pencilUseX = propX.magnitude >= propY.magnitude;
    ev.pencilNs = std::chrono::duration<double, std::nano>(t1 - t0).count();
    auto const fx = proposeFrequency(resX, span);
    auto const fy = hasY ? proposeFrequency(resY, span) : FrequencyProposal{};
    auto const t2 = std::chrono::steady_clock::now();
    ev.fftNs = std::chrono::duration<double, std::nano>(t2 - t1).count();
    ev.fftValidX = fx.valid; ev.fftOmegaX = fx.omega; ev.fftMagX = fx.magnitude;
    ev.fftValidY = fy.valid; ev.fftOmegaY = fy.omega; ev.fftMagY = fy.magnitude;
    trace(ev);
}

}
}
