#pragma once


#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "weighted_varpro.hpp"

namespace Sketch {
namespace freqdisc {

constexpr double kTwoPi = 6.283185307179586;
constexpr double kEps = 1e-12;

inline auto hannWindow(std::size_t n) -> std::vector<double>
{
    std::vector<double> w(n, 1.0);
    if (n < 2) { return w; }
    double const denom = static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        w[i] = 0.5 * (1.0 - std::cos(kTwoPi * static_cast<double>(i) / denom));
    }
    return w;
}

inline auto realDftMagnitude(std::vector<double> const& signal, std::vector<double>& outRe,
                              std::vector<double>& outIm) -> std::vector<double>
{
    std::size_t const n = signal.size();
    std::size_t const halfN = n / 2;
    outRe.assign(halfN + 1, 0.0);
    outIm.assign(halfN + 1, 0.0);
    std::vector<double> mag(halfN + 1, 0.0);
    if (n < 4) { return mag; }
    for (std::size_t k = 0; k <= halfN; ++k) {
        double re = 0.0, im = 0.0;
        double const w = -kTwoPi * static_cast<double>(k) / static_cast<double>(n);
        for (std::size_t t = 0; t < n; ++t) {
            double const angle = w * static_cast<double>(t);
            re += signal[t] * std::cos(angle);
            im += signal[t] * std::sin(angle);
        }
        outRe[k] = re; outIm[k] = im;
        double const m = std::sqrt(re * re + im * im);
        mag[k] = std::isfinite(m) ? m : 0.0;
    }
    return mag;
}

struct FrequencyProposal {
    bool valid{false};
    double omega{0.0};
    double magnitude{0.0};
    double damping{0.0};
};

inline auto proposeFrequency(std::vector<double> const& residual, double span) -> FrequencyProposal
{
    FrequencyProposal out;
    std::size_t const n = residual.size();
    if (n < 8 || !(span > kEps) || !std::isfinite(span)) { return out; }

    auto win = hannWindow(n);
    std::vector<double> windowed(n);
    double mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const v = residual[i];
        if (!std::isfinite(v)) { return out; }
        mean += v;
    }
    mean /= static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) { windowed[i] = (residual[i] - mean) * win[i]; }

    std::vector<double> re, im;
    auto mag = realDftMagnitude(windowed, re, im);
    if (mag.size() < 3) { return out; }

    std::size_t k0 = 0;
    double best = 0.0;
    for (std::size_t k = 1; k + 1 < mag.size(); ++k) {
        if (mag[k] > best) { best = mag[k]; k0 = k; }
    }
    if (k0 == 0 || !(best > kEps)) { return out; }

    double const ym1 = mag[k0 - 1], y0 = mag[k0], yp1 = mag[k0 + 1];
    double const denom = (ym1 - 2.0 * y0 + yp1);
    double delta = 0.0;
    if (std::fabs(denom) > kEps) {
        delta = 0.5 * (ym1 - yp1) / denom;
        if (!std::isfinite(delta) || std::fabs(delta) > 1.0) { delta = 0.0; }
    }
    double const kRefined = static_cast<double>(k0) + delta;
    double const freqCyclesPerSpan = kRefined;
    double const omega = kTwoPi * freqCyclesPerSpan / span;
    if (!std::isfinite(omega) || omega <= 0.0) { return out; }

    out.valid = true;
    out.omega = omega;
    out.magnitude = best;
    return out;
}

inline auto proposeFrequencyWeighted(std::vector<double> const& residual, double span,
                                      std::vector<double> const& sigma = {}) -> FrequencyProposal
{
    FrequencyProposal out;
    std::size_t const n = residual.size();
    if (n < 8 || !(span > kEps) || !std::isfinite(span)) { return out; }

    std::vector<double> weight;
    if (sigma.size() == n) {
        weight.resize(n);
        bool sigmaOk = true;
        double wmax = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double const s = sigma[i];
            double const s2 = s * s;
            if (!std::isfinite(s) || !(s2 > kEps)) { sigmaOk = false; break; }
            double const wi = 1.0 / s2;
            weight[i] = wi;
            if (wi > wmax) { wmax = wi; }
        }
        if (!sigmaOk || !(wmax > 0.0) || !std::isfinite(wmax)) {
            weight.clear();
        } else {
            for (std::size_t i = 0; i < n; ++i) { weight[i] /= wmax; }
        }
    }

    auto win = hannWindow(n);
    std::vector<double> windowed(n);
    double mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const v = residual[i];
        if (!std::isfinite(v)) { return out; }
        mean += v;
    }
    mean /= static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) { windowed[i] = (residual[i] - mean) * win[i]; }

    std::vector<double> weighted(n);
    for (std::size_t i = 0; i < n; ++i) {
        double const wi = weight.empty() ? 1.0 : weight[i];
        weighted[i] = wi * windowed[i];
    }

    std::vector<double> re, im;
    auto mag = realDftMagnitude(weighted, re, im);
    if (mag.size() < 3) { return out; }

    std::size_t k0 = 0;
    double best = 0.0;
    for (std::size_t k = 1; k + 1 < mag.size(); ++k) {
        if (mag[k] > best) { best = mag[k]; k0 = k; }
    }
    if (k0 == 0 || !(best > kEps)) { return out; }

    double const ym1 = mag[k0 - 1], y0 = mag[k0], yp1 = mag[k0 + 1];
    double const denom = (ym1 - 2.0 * y0 + yp1);
    double delta = 0.0;
    if (std::fabs(denom) > kEps) {
        delta = 0.5 * (ym1 - yp1) / denom;
        if (!std::isfinite(delta) || std::fabs(delta) > 1.0) { delta = 0.0; }
    }
    double const kRefined = static_cast<double>(k0) + delta;
    double const freqCyclesPerSpan = kRefined;
    double const omega = kTwoPi * freqCyclesPerSpan / span;
    if (!std::isfinite(omega) || omega <= 0.0) { return out; }

    out.valid = true;
    out.omega = omega;
    out.magnitude = best;
    return out;
}

struct QuadratureFit {
    bool valid{false};
    double a{0.0};
    double b{0.0};
    double amplitude{0.0};
    double phase{0.0};
    double residualNorm{0.0};
};

inline auto quadratureFit(std::vector<double> const& t, std::vector<double> const& target,
                           double omega) -> QuadratureFit
{
    QuadratureFit out;
    std::size_t const n = t.size();
    if (n < 2 || target.size() != n || !std::isfinite(omega)) { return out; }

    double Sss = 0.0, Scc = 0.0, Ssc = 0.0, Sys = 0.0, Syc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(t[i]) || !std::isfinite(target[i])) { return out; }
        double const s = std::sin(omega * t[i]);
        double const c = std::cos(omega * t[i]);
        Sss += s * s; Scc += c * c; Ssc += s * c;
        Sys += target[i] * s; Syc += target[i] * c;
    }
    double const det = Sss * Scc - Ssc * Ssc;
    if (!std::isfinite(det) || std::fabs(det) < kEps) { return out; }

    double const A = (Sys * Scc - Syc * Ssc) / det;
    double const B = (Sss * Syc - Ssc * Sys) / det;
    if (!std::isfinite(A) || !std::isfinite(B)) { return out; }

    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const s = std::sin(omega * t[i]);
        double const c = std::cos(omega * t[i]);
        double const r = target[i] - (A * s + B * c);
        sse += r * r;
    }

    out.valid = true;
    out.a = A;
    out.b = B;
    out.amplitude = std::sqrt(A * A + B * B);
    out.phase = std::atan2(B, A);
    out.residualNorm = std::sqrt(std::max(0.0, sse));
    return out;
}

struct FrequencyGnResult {
    double omega{0.0};
    QuadratureFit fit{};
    bool valid{false};
    bool accepted{false};
};

inline auto refineQuadratureFrequencyOneStep(std::vector<double> const& t,
                                              std::vector<double> const& target,
                                              double omega,
                                              std::vector<double> const& variance = {}) -> FrequencyGnResult
{
    FrequencyGnResult out;
    out.omega = omega;
    out.fit = quadratureFit(t, target, omega);
    out.valid = out.fit.valid;
    if (!out.valid) { return out; }

    std::size_t const n = t.size();
    std::vector<double> jac(n), residual(n);
    for (std::size_t i = 0; i < n; ++i) {
        double const arg = omega * t[i];
        double const s = std::sin(arg), c = std::cos(arg);
        double const yhat = out.fit.a * s + out.fit.b * c;
        residual[i] = target[i] - yhat;
        jac[i] = t[i] * (out.fit.a * c - out.fit.b * s);
    }
    WeightedVarProConfig cfg;
    cfg.ridgeFraction = 1e-3;
    cfg.useInverseVariance = !variance.empty();
    auto const step = solveWeightedVarPro(jac, n, 1, residual, variance, cfg);
    if (!step.ok || step.coefficients.size() != 1 || !std::isfinite(step.coefficients[0])) {
        return out;
    }

    double const trialOmega = omega + step.coefficients[0];
    if (!(trialOmega > 0.0) || !std::isfinite(trialOmega)) { return out; }
    auto trial = quadratureFit(t, target, trialOmega);
    if (trial.valid && trial.residualNorm < out.fit.residualNorm) {
        out.omega = trialOmega;
        out.fit = trial;
        out.accepted = true;
    }
    return out;
}

}
}
