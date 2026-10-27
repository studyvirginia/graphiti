#ifndef SKETCH_SHAPE_CONST_OPTIMIZER_HPP
#define SKETCH_SHAPE_CONST_OPTIMIZER_HPP

#include <operon/operon/core/types.hpp>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>
#include <Eigen/Dense>

namespace Sketch {

struct NelderMead {
    int maxIter = 200;
    Operon::Scalar xatol = 1e-3;
    Operon::Scalar fatol = 1e-6;

    auto Minimize(std::vector<Operon::Scalar> x0,
                  std::function<Operon::Scalar(std::vector<Operon::Scalar> const&)> const& loss,
                  Operon::Scalar& bestOut) const -> std::vector<Operon::Scalar>
    {
        auto const n = x0.size();
        if (n == 0) { bestOut = loss(x0); return x0; }

        std::vector<std::vector<Operon::Scalar>> simplex(n + 1, x0);
        std::vector<Operon::Scalar> f(n + 1);
        for (std::size_t i = 0; i < n; ++i) {
            auto& pt = simplex[i + 1];
            auto step = (std::abs(pt[i]) > 1e-6) ? Operon::Scalar{0.05} * pt[i]
                                                 : Operon::Scalar{0.05};
            pt[i] += step;
        }
        for (std::size_t i = 0; i <= n; ++i) { f[i] = loss(simplex[i]); }

        std::vector<std::size_t> order(n + 1);
        auto reorder = [&]() {
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](std::size_t a, std::size_t b) {
                          auto const fa = f[a], fb = f[b];
                          bool const finA = std::isfinite(fa), finB = std::isfinite(fb);
                          if (finA && finB) { return fa < fb; }
                          if (finA != finB) { return finA; }
                          return false;
                      });
        };

        constexpr Operon::Scalar ALPHA = 1.0, GAMMA = 2.0, RHO = 0.5, SIGMA = 0.5;

        auto shrink = [&]() {
            auto b = order.front();
            for (std::size_t i = 0; i <= n; ++i) {
                if (i == b) { continue; }
                for (std::size_t k = 0; k < n; ++k) {
                    simplex[i][k] = simplex[b][k] + SIGMA * (simplex[i][k] - simplex[b][k]);
                }
                f[i] = loss(simplex[i]);
            }
        };

        for (int iter = 0; iter < maxIter; ++iter) {
            reorder();
            auto best = order.front(), worst = order.back();
            auto second = (n >= 2) ? order[n - 1] : worst;

            if (std::abs(f[worst] - f[best]) <= fatol) {
                Operon::Scalar extent = 0;
                for (std::size_t i = 0; i <= n; ++i) {
                    if (i == best) { continue; }
                    for (std::size_t k = 0; k < n; ++k) {
                        auto d = std::abs(simplex[i][k] - simplex[best][k]);
                        if (d > extent) { extent = d; }
                    }
                }
                if (extent <= xatol) { break; }
            }

            static thread_local std::vector<Operon::Scalar> centroid, refl, expd, cont;
            centroid.assign(n, 0.0f);
            refl.resize(n);
            expd.resize(n);
            cont.resize(n);

            for (std::size_t i = 0; i <= n; ++i) {
                if (i == worst) { continue; }
                for (std::size_t k = 0; k < n; ++k) { centroid[k] += simplex[i][k]; }
            }
            for (auto& c : centroid) { c /= static_cast<Operon::Scalar>(n); }

            auto combineInto = [&](Operon::Scalar coeff, std::vector<Operon::Scalar>& out) {
                for (std::size_t k = 0; k < n; ++k) {
                    out[k] = centroid[k] + coeff * (centroid[k] - simplex[worst][k]);
                }
            };

            combineInto(ALPHA, refl);
            auto fr = loss(refl);
            if (fr < f[best]) {
                combineInto(GAMMA, expd);
                auto fe = loss(expd);
                if (fe < fr) { simplex[worst] = expd; f[worst] = fe; }
                else         { simplex[worst] = refl; f[worst] = fr; }
            } else if (fr < f[second]) {
                simplex[worst] = refl; f[worst] = fr;
            } else if (fr < f[worst]) {
                combineInto(RHO * ALPHA, cont);
                auto fc = loss(cont);
                if (fc <= fr) { simplex[worst] = cont; f[worst] = fc; }
                else          { shrink(); }
            } else {
                combineInto(-RHO, cont);
                auto fc = loss(cont);
                if (fc < f[worst]) { simplex[worst] = cont; f[worst] = fc; }
                else               { shrink(); }
            }
        }
        reorder();
        bestOut = f[order.front()];
        return simplex[order.front()];
    }
};

struct LevenbergMarquardt {
    int maxIter = 8;
    int maxTrials = 6;
    double ftol = 1e-6;
    double lambda0 = 1e-3;
    double lambdaUp = 10.0;
    double lambdaDown = 0.1;
    double lambdaMin = 1e-12;
    double lambdaMax = 1e12;

    double floorThreshold = 0.0;
    bool enableScreening = true;
    double screenLoss = 0.05;

    bool resolutionStop = false;
    double dataCount = 0.0;

    bool gaugeProject = false;
    double jacobianEps = DBL_EPSILON;
    mutable std::size_t residualRows = 0;
    mutable std::vector<int> gaugeNullDimHistory;
    static auto projectGaugeNull(std::vector<double> const& H, std::size_t n, double eps, std::size_t m,
                                  std::vector<double>& dx) -> int
    {
        if (n == 0) { return 0; }
        Eigen::MatrixXd Hm(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                Hm(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = H[i * n + j];
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hm);
        auto const& evals = es.eigenvalues();
        double const lambdaMaxH = evals(static_cast<Eigen::Index>(n) - 1);
        double const dim = static_cast<double>(std::max(m, n));
        double const rankTol = dim * eps;
        double const thresh = rankTol * rankTol * std::max(lambdaMaxH, 0.0);
        Eigen::VectorXd dxv(static_cast<Eigen::Index>(n));
        for (std::size_t i = 0; i < n; ++i) { dxv(static_cast<Eigen::Index>(i)) = dx[i]; }
        Eigen::VectorXd proj = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n));
        int nullDim = 0;
        for (Eigen::Index k = 0; k < static_cast<Eigen::Index>(n); ++k) {
            if (evals(k) <= thresh) {
                ++nullDim;
                Eigen::VectorXd const v = es.eigenvectors().col(k);
                proj += v * v.dot(dxv);
            }
        }
        dxv -= proj;
        for (std::size_t i = 0; i < n; ++i) { dx[i] = dxv(static_cast<Eigen::Index>(i)); }
        return nullDim;
    }

    mutable int itersRun = 0;
    mutable int lossCalls = 0;
    mutable int builds = 0;
    mutable int stopUnresolvableDamped = 0;
    mutable int stopRealizedSmall = 0;
    mutable int stopFixedPoint = 0;

    auto Minimize(std::vector<Operon::Scalar> x0,
                  std::function<bool(std::vector<Operon::Scalar> const&,
                                     std::vector<double>&, std::vector<double>&)> const& build,
                  std::function<Operon::Scalar(std::vector<Operon::Scalar> const&)> const& loss,
                  Operon::Scalar& bestOut) const -> std::vector<Operon::Scalar>
    {
        itersRun = 0; lossCalls = 0; builds = 0;
        stopUnresolvableDamped = 0; stopRealizedSmall = 0; stopFixedPoint = 0;
        if (gaugeProject) { gaugeNullDimHistory.clear(); }
        auto const n = x0.size();
        if (n == 0) { bestOut = loss(x0); ++lossCalls; return x0; }

        auto x = x0;
        auto f = loss(x); ++lossCalls;
        bestOut = f;
        if (!std::isfinite(static_cast<double>(f))) { return x; }

        if (floorThreshold > 0.0 && static_cast<double>(f) <= floorThreshold) {
            return x;
        }

        std::vector<double> r, J;
        std::vector<double> H(n * n), g(n), diag(n), dx(n);
        double lambda = lambda0;

        static thread_local std::vector<double> L_scratch, y_scratch;
        L_scratch.resize(n * n);
        y_scratch.resize(n);

        for (int it = 0; it < maxIter; ++it) {
            r.clear(); J.clear();
            ++builds;
            if (!build(x, r, J)) { break; }
            auto const m = r.size();
            if (m == 0 || J.size() != m * n) { break; }
            if (gaugeProject) { residualRows = m; }

            std::fill(H.begin(), H.end(), 0.0);
            std::fill(g.begin(), g.end(), 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                double const* __restrict row = J.data() + i * n;
                double const ri = r[i];
                for (std::size_t a = 0; a < n; ++a) {
                    double const ra = row[a];
                    g[a] += ra * ri;
                    double* __restrict Ha = H.data() + a * n;
                    for (std::size_t b = a; b < n; ++b) {
                        Ha[b] += ra * row[b];
                    }
                }
            }
            for (std::size_t a = 0; a < n; ++a) {
                for (std::size_t b = 0; b < a; ++b) { H[a * n + b] = H[b * n + a]; }
                diag[a] = H[a * n + a];
                if (!(diag[a] > 0.0)) { diag[a] = 1.0; }
            }

            if (resolutionStop) {
                double const fOld = static_cast<double>(f);
                static thread_local std::vector<Operon::Scalar> xBefore;
                xBefore = x;
                bool unresolvable = false;
                for (;;) {
                    if (!std::isfinite(lambda)) { unresolvable = true; break; }
                    if (!solveDamped(H, g, diag, lambda, n, dx, L_scratch, y_scratch)) {
                        lambda *= lambdaUp;
                        continue;
                    }
                    if (gaugeProject) { gaugeNullDimHistory.push_back(projectGaugeNull(H, n, jacobianEps, residualRows, dx)); }
                    double gdx = 0.0;
                    for (std::size_t a = 0; a < n; ++a) { gdx += g[a] * dx[a]; }
                    double dHd = 0.0;
                    for (std::size_t a = 0; a < n; ++a) {
                        double const* __restrict Ha = H.data() + a * n;
                        double hd = 0.0;
                        for (std::size_t b = 0; b < n; ++b) { hd += Ha[b] * dx[b]; }
                        dHd += dx[a] * hd;
                    }
                    double const pred = -2.0 * gdx - dHd;
                    bool resolvable;
                    if (pred <= 0.0) {
                        resolvable = false;
                    } else if (fOld - pred <= 0.0) {
                        resolvable = true;
                    } else {
                        double const gain = (dataCount / 2.0) * std::log(fOld / (fOld - pred));
                        resolvable = (gain >= 0.5);
                    }
                    if (!resolvable) { unresolvable = true; break; }

                    auto xn = x;
                    bool finite = true;
                    for (std::size_t a = 0; a < n; ++a) {
                        auto const v = static_cast<double>(x[a]) + dx[a];
                        if (!std::isfinite(v)) { finite = false; break; }
                        xn[a] = static_cast<Operon::Scalar>(v);
                    }
                    if (!finite) { lambda *= lambdaUp; continue; }
                    auto const fn = loss(xn); ++lossCalls;
                    if (std::isfinite(static_cast<double>(fn)) && fn < f) {
                        x = std::move(xn);
                        f = fn;
                        lambda = std::max(lambda * lambdaDown, lambdaMin);
                        break;
                    }
                    lambda *= lambdaUp;
                }
                ++itersRun;
                if (unresolvable) { ++stopUnresolvableDamped; break; }
                bool same = (x.size() == xBefore.size());
                for (std::size_t a = 0; same && a < x.size(); ++a) { same = same && (x[a] == xBefore[a]); }
                if (same) { ++stopFixedPoint; break; }
                double const gain = (dataCount / 2.0) * std::log(fOld / static_cast<double>(f));
                if (gain < 0.5) { ++stopRealizedSmall; break; }
                continue;
            }

            bool accepted = false;
            double rel = 0.0;
            for (int trial = 0; trial < maxTrials; ++trial) {
                if (!solveDamped(H, g, diag, lambda, n, dx, L_scratch, y_scratch)) { lambda = std::min(lambda * lambdaUp, lambdaMax); continue; }
                if (gaugeProject) { gaugeNullDimHistory.push_back(projectGaugeNull(H, n, jacobianEps, residualRows, dx)); }
                auto xn = x;
                bool finite = true;
                for (std::size_t a = 0; a < n; ++a) {
                    auto const v = static_cast<double>(x[a]) + dx[a];
                    if (!std::isfinite(v)) { finite = false; break; }
                    xn[a] = static_cast<Operon::Scalar>(v);
                }
                if (!finite) { lambda = std::min(lambda * lambdaUp, lambdaMax); continue; }
                auto const fn = loss(xn); ++lossCalls;
                if (std::isfinite(static_cast<double>(fn)) && fn < f) {
                    rel = (static_cast<double>(f) - static_cast<double>(fn))
                        / std::max(std::abs(static_cast<double>(f)), 1e-30);
                    x = std::move(xn);
                    f = fn;
                    lambda = std::max(lambda * lambdaDown, lambdaMin);
                    accepted = true;
                    break;
                }
                lambda = std::min(lambda * lambdaUp, lambdaMax);
                if (enableScreening && it == 0 && trial == 0 && static_cast<double>(f) > screenLoss) {
                    break;
                }
            }
            ++itersRun;
            if (!accepted) { break; }
            if (rel < ftol) { break; }
        }
        bestOut = f;
        return x;
    }

    auto MinimizeVarPro(std::vector<Operon::Scalar> x0,
                         std::function<bool(std::vector<Operon::Scalar> const&,
                                            std::vector<double>&, std::vector<double>&)> const& build,
                         std::function<Operon::Scalar(std::vector<Operon::Scalar> const&)> const& loss,
                         Operon::Scalar& bestOut,
                         std::uint64_t linMask) const -> std::vector<Operon::Scalar>
    {
        auto const n = x0.size();
        std::vector<std::size_t> linIdx, nlIdx;
        if (n > 0 && n <= 64 && resolutionStop) {
            for (std::size_t k = 0; k < n; ++k) {
                if ((linMask >> k) & 1ull) { linIdx.push_back(k); }
                else { nlIdx.push_back(k); }
            }
        }
        if (n == 0 || n > 64 || !resolutionStop || linIdx.empty() || nlIdx.empty()) {
            return Minimize(std::move(x0), build, loss, bestOut);
        }

        itersRun = 0; lossCalls = 0; builds = 0;
        stopUnresolvableDamped = 0; stopRealizedSmall = 0; stopFixedPoint = 0;
        if (gaugeProject) { gaugeNullDimHistory.clear(); }

        auto x = x0;
        auto f = loss(x); ++lossCalls;
        bestOut = f;
        if (!std::isfinite(static_cast<double>(f))) { return x; }
        if (floorThreshold > 0.0 && static_cast<double>(f) <= floorThreshold) { return x; }

        std::size_t const p = linIdx.size();
        std::size_t const q = nlIdx.size();
        std::vector<double> r, J;
        std::vector<double> H(n * n), g(n);
        std::vector<double> Hl(p * p), gl(p), diagl(p), dl(p);
        std::vector<double> Hnl(q * q), gnl(q), diagnl(q), dnl(q);
        double lambda = lambda0;

        static thread_local std::vector<double> Ll_scratch, yl_scratch, Lnl_scratch, ynl_scratch;

        auto formHg = [&](std::vector<double> const& rr, std::vector<double> const& JJ, std::size_t m) {
            std::fill(H.begin(), H.end(), 0.0);
            std::fill(g.begin(), g.end(), 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                double const* __restrict row = JJ.data() + i * n;
                double const ri = rr[i];
                for (std::size_t a = 0; a < n; ++a) {
                    double const ra = row[a];
                    g[a] += ra * ri;
                    double* __restrict Ha = H.data() + a * n;
                    for (std::size_t b = a; b < n; ++b) { Ha[b] += ra * row[b]; }
                }
            }
            for (std::size_t a = 0; a < n; ++a) {
                for (std::size_t b = 0; b < a; ++b) { H[a * n + b] = H[b * n + a]; }
            }
        };

        for (int it = 0; it < maxIter; ++it) {
            r.clear(); J.clear();
            ++builds;
            if (!build(x, r, J)) { break; }
            auto const m = r.size();
            if (m == 0 || J.size() != m * n) { break; }
            if (gaugeProject) { residualRows = m; }
            formHg(r, J, m);

            for (std::size_t a = 0; a < p; ++a) {
                for (std::size_t b = 0; b < p; ++b) { Hl[a * p + b] = H[linIdx[a] * n + linIdx[b]]; }
                gl[a] = g[linIdx[a]];
                diagl[a] = Hl[a * p + a];
                if (!(diagl[a] > 0.0)) { diagl[a] = 1.0; }
            }
            if (solveDamped(Hl, gl, diagl, 0.0, p, dl, Ll_scratch, yl_scratch)) {
                auto xTrial = x;
                bool finite = true;
                for (std::size_t a = 0; a < p; ++a) {
                    auto const v = static_cast<double>(x[linIdx[a]]) + dl[a];
                    if (!std::isfinite(v)) { finite = false; break; }
                    xTrial[linIdx[a]] = static_cast<Operon::Scalar>(v);
                }
                if (finite) {
                    auto const fTrial = loss(xTrial); ++lossCalls;
                    if (std::isfinite(static_cast<double>(fTrial)) && fTrial <= static_cast<double>(f)) {
                        x = std::move(xTrial);
                        f = fTrial;
                    }
                }
            }

            r.clear(); J.clear();
            ++builds;
            if (!build(x, r, J)) { break; }
            auto const m2 = r.size();
            if (m2 == 0 || J.size() != m2 * n) { break; }
            if (gaugeProject) { residualRows = m2; }
            formHg(r, J, m2);

            for (std::size_t a = 0; a < q; ++a) {
                for (std::size_t b = 0; b < q; ++b) { Hnl[a * q + b] = H[nlIdx[a] * n + nlIdx[b]]; }
                gnl[a] = g[nlIdx[a]];
                diagnl[a] = Hnl[a * q + a];
                if (!(diagnl[a] > 0.0)) { diagnl[a] = 1.0; }
            }

            double const fOld = static_cast<double>(f);
            static thread_local std::vector<Operon::Scalar> xBefore;
            xBefore = x;
            bool unresolvable = false;
            for (;;) {
                if (!std::isfinite(lambda)) { unresolvable = true; break; }
                if (!solveDamped(Hnl, gnl, diagnl, lambda, q, dnl, Lnl_scratch, ynl_scratch)) {
                    lambda *= lambdaUp;
                    continue;
                }
                if (gaugeProject) { gaugeNullDimHistory.push_back(projectGaugeNull(Hnl, q, jacobianEps, residualRows, dnl)); }
                double gdx = 0.0;
                for (std::size_t a = 0; a < q; ++a) { gdx += gnl[a] * dnl[a]; }
                double dHd = 0.0;
                for (std::size_t a = 0; a < q; ++a) {
                    double const* __restrict Ha = Hnl.data() + a * q;
                    double hd = 0.0;
                    for (std::size_t b = 0; b < q; ++b) { hd += Ha[b] * dnl[b]; }
                    dHd += dnl[a] * hd;
                }
                double const pred = -2.0 * gdx - dHd;
                bool resolvable;
                if (pred <= 0.0) {
                    resolvable = false;
                } else if (fOld - pred <= 0.0) {
                    resolvable = true;
                } else {
                    double const gain = (dataCount / 2.0) * std::log(fOld / (fOld - pred));
                    resolvable = (gain >= 0.5);
                }
                if (!resolvable) { unresolvable = true; break; }

                auto xn = x;
                bool finite = true;
                for (std::size_t a = 0; a < q; ++a) {
                    auto const v = static_cast<double>(x[nlIdx[a]]) + dnl[a];
                    if (!std::isfinite(v)) { finite = false; break; }
                    xn[nlIdx[a]] = static_cast<Operon::Scalar>(v);
                }
                if (!finite) { lambda *= lambdaUp; continue; }
                auto const fn = loss(xn); ++lossCalls;
                if (std::isfinite(static_cast<double>(fn)) && fn < f) {
                    x = std::move(xn);
                    f = fn;
                    lambda = std::max(lambda * lambdaDown, lambdaMin);
                    break;
                }
                lambda *= lambdaUp;
            }
            ++itersRun;
            if (unresolvable) { ++stopUnresolvableDamped; break; }
            bool same = (x.size() == xBefore.size());
            for (std::size_t a = 0; same && a < x.size(); ++a) { same = same && (x[a] == xBefore[a]); }
            if (same) { ++stopFixedPoint; break; }
            double const gain = (dataCount / 2.0) * std::log(fOld / static_cast<double>(f));
            if (gain < 0.5) { ++stopRealizedSmall; break; }
            continue;
        }
        bestOut = f;
        return x;
    }

    static auto solveDamped(std::vector<double> const& H, std::vector<double> const& g,
                            std::vector<double> const& diag, double lambda, std::size_t n,
                            std::vector<double>& dx,
                            std::vector<double>& L, std::vector<double>& y) -> bool
    {
        if (L.size() < n * n) { L.resize(n * n); }
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double s = H[i * n + j];
                if (i == j) { s += lambda * diag[i]; }
                for (std::size_t k = 0; k < j; ++k) { s -= L[i * n + k] * L[j * n + k]; }
                if (i == j) {
                    if (!(s > 0.0) || !std::isfinite(s)) { return false; }
                    L[i * n + i] = std::sqrt(s);
                } else {
                    L[i * n + j] = s / L[j * n + j];
                }
            }
        }
        dx.assign(n, 0.0);
        y.assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double s = -g[i];
            for (std::size_t k = 0; k < i; ++k) { s -= L[i * n + k] * y[k]; }
            y[i] = s / L[i * n + i];
        }
        for (std::size_t ii = n; ii-- > 0;) {
            double s = y[ii];
            for (std::size_t k = ii + 1; k < n; ++k) { s -= L[k * n + ii] * dx[k]; }
            dx[ii] = s / L[ii * n + ii];
            if (!std::isfinite(dx[ii])) { return false; }
        }
        return true;
    }
};

}
#endif
