#ifndef SKETCH_MULTITREE_GP_TYPES_HPP
#define SKETCH_MULTITREE_GP_TYPES_HPP

namespace Sketch {

#ifdef SKETCH_DUP_STATS
struct DupStats {
    static constexpr std::uint64_t kFnvBasis = 14695981039346656037ULL;
    static constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
    static void mix(std::uint64_t& h, void const* p, std::size_t n) {
        auto const* b = static_cast<unsigned char const*>(p);
        for (std::size_t i = 0; i < n; ++i) { h ^= b[i]; h *= kFnvPrime; }
    }
    static void mixTree(std::uint64_t& h, Operon::Tree const& t) {
        for (auto const& n : t.Nodes()) {
            auto ty = static_cast<std::uint8_t>(n.Type);
            mix(h, &ty, sizeof ty);
            mix(h, &n.HashValue, sizeof n.HashValue);
            mix(h, &n.Arity, sizeof n.Arity);
            mix(h, &n.Length, sizeof n.Length);
        }
        std::uint8_t sep = 0xFF; mix(h, &sep, sizeof sep);
    }
    static void mixCoefs(std::uint64_t& h, std::vector<Operon::Scalar> const& c) {
        std::size_t n = c.size(); mix(h, &n, sizeof n);
        if (n) { mix(h, c.data(), n * sizeof(Operon::Scalar)); }
    }
    mutable std::mutex mu;
    std::unordered_set<std::uint64_t> exact;
    std::unordered_set<std::uint64_t> structural;
    std::size_t total{0}, exactDup{0}, structOnlyDup{0}, empty{0};
    void reset() { std::lock_guard<std::mutex> g(mu); exact.clear(); structural.clear(); total = exactDup = structOnlyDup = empty = 0; }
    void record(std::uint64_t hs, std::uint64_t hx, bool hasCoefs) {
        std::lock_guard<std::mutex> g(mu);
        ++total;
        if (!hasCoefs) { ++empty; return; }
        bool const sSeen = !structural.insert(hs).second;
        bool const xSeen = !exact.insert(hx).second;
        if (xSeen) { ++exactDup; } else if (sSeen) { ++structOnlyDup; }
    }
    void print(char const* tag) const {
        std::lock_guard<std::mutex> g(mu);
        std::size_t const eff = total - empty;
        std::cerr << "[SKETCH_DUP_STATS " << tag << "] total_lm_calls=" << total
                  << " (no_coef=" << empty << ")"
                  << " exact_dup_calls=" << exactDup
                  << " exact_dup_rate=" << (eff ? double(exactDup) / double(eff) : 0.0)
                  << " structure_only_dup_calls=" << structOnlyDup
                  << " structure_only_dup_rate=" << (eff ? double(structOnlyDup) / double(eff) : 0.0)
                  << " unique_structures=" << structural.size()
                  << " unique_exact_keys=" << exact.size() << "\n";
    }
};
#endif

struct CisPointSegBatch {
    std::size_t N{0};
    std::vector<Operon::Scalar> fitness;
    std::vector<Operon::Scalar> residualFlat;
    std::vector<Operon::Scalar> tauFlat;
    std::vector<bool> valid;
    std::vector<Operon::Scalar> reverseFitness;
    int fp16ModeUsed{0};
    PointSegGpuStream const* stream{nullptr};
    std::size_t P{0};

    auto viewCandidate(std::size_t idx, std::size_t n, Operon::Scalar& fit, Operon::Scalar const*& res,
                       Operon::Scalar const*& tauOut) const -> bool
    {
        if (stream != nullptr) {
            if (idx >= P || !stream->waitFor(idx)) { return false; }
            fit = stream->fitness()[idx];
            res = stream->residual() + idx * n;
            tauOut = stream->tau() + idx * n;
            return true;
        }
        if (!(idx < fitness.size() && (idx + 1) * n <= residualFlat.size())) { return false; }
        fit = fitness[idx];
        res = residualFlat.data() + idx * n;
        tauOut = ((idx + 1) * n <= tauFlat.size()) ? tauFlat.data() + idx * n : nullptr;
        return true;
    }
};

struct GradMutErrModel { double mu{0}, b{0}, s0{1}, s1{0}, nu{3}, pi0{0}; };

namespace gmstat {
inline auto trigamma(double x) -> double
{
    double acc = 0.0;
    while (x < 6.0) { acc += 1.0 / (x * x); x += 1.0; }
    double const r = 1.0 / x, r2 = r * r;
    return acc + r + 0.5 * r2 + r * r2 * (1.0 / 6.0 - r2 * (1.0 / 30.0 - r2 * (1.0 / 42.0 - r2 * (1.0 / 30.0))));
}
inline auto logJeffreysNu(double nu) -> double
{
    if (!(nu > 0.0) || !std::isfinite(nu)) { return -std::numeric_limits<double>::infinity(); }
    if (nu >= 100.0) {
        return 0.5 * std::log(nu / (nu + 3.0)) + 0.5 * std::log(6.0 * (1.0 - 2.0 / nu)) - 2.0 * std::log(nu);
    }
    double const br = trigamma(0.5 * nu) - trigamma(0.5 * (nu + 1.0)) - 2.0 * (nu + 3.0) / (nu * (nu + 1.0) * (nu + 1.0));
    if (!(br > 0.0) || !std::isfinite(br)) { return -std::numeric_limits<double>::infinity(); }
    return 0.5 * std::log(nu / (nu + 3.0)) + 0.5 * std::log(br);
}
inline auto betacf(double a, double b, double x) -> double
{
    double const tiny = 1e-300;
    double const qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (std::abs(d) < tiny) { d = tiny; }
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= 100000; ++m) {
        double const dm = m, m2 = 2.0 * m;
        double aa = dm * (b - dm) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::abs(d) < tiny) { d = tiny; }
        c = 1.0 + aa / c; if (std::abs(c) < tiny) { c = tiny; }
        d = 1.0 / d; h *= d * c;
        aa = -(a + dm) * (qab + dm) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::abs(d) < tiny) { d = tiny; }
        c = 1.0 + aa / c; if (std::abs(c) < tiny) { c = tiny; }
        d = 1.0 / d;
        double const del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < 1e-15) { break; }
    }
    return h;
}
inline auto ibeta(double a, double b, double x) -> double
{
    if (!(x > 0.0)) { return 0.0; }
    if (!(x < 1.0)) { return 1.0; }
    double const lbt = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log1p(-x);
    if (x < (a + 1.0) / (a + b + 2.0)) { return std::exp(lbt) * betacf(a, b, x) / a; }
    return 1.0 - std::exp(lbt) * betacf(b, a, 1.0 - x) / b;
}
inline auto tCdf(double x, double nu) -> double
{
    if (std::isnan(x)) { return 0.5; }
    if (std::isinf(x)) { return x < 0.0 ? 0.0 : 1.0; }
    double const tail = 0.5 * ibeta(0.5 * nu, 0.5, nu / (nu + x * x));
    return x < 0.0 ? tail : 1.0 - tail;
}
template <class F>
inline auto nelderMead(F&& f, std::vector<double> x0, std::vector<double> const& step, int maxEval) -> std::vector<double>
{
    std::size_t const n = x0.size();
    std::vector<std::vector<double>> p(n + 1, x0);
    std::vector<double> fv(n + 1);
    for (std::size_t i = 0; i < n; ++i) { p[i + 1][i] += step[i]; }
    int ev = 0;
    for (std::size_t i = 0; i <= n; ++i) { fv[i] = f(p[i]); ++ev; }
    std::vector<double> c(n), xr(n), xe(n), xc(n);
    while (ev < maxEval) {
        std::vector<std::size_t> o(n + 1);
        std::iota(o.begin(), o.end(), 0);
        std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });
        auto const lo = o[0], hi = o[n], nh = o[n - 1];
        if (std::isfinite(fv[hi]) && std::abs(fv[hi] - fv[lo]) <= 1e-9 * (std::abs(fv[lo]) + 1e-300)) { break; }
        std::fill(c.begin(), c.end(), 0.0);
        for (std::size_t i = 0; i <= n; ++i) { if (i != hi) { for (std::size_t k = 0; k < n; ++k) { c[k] += p[i][k] / n; } } }
        for (std::size_t k = 0; k < n; ++k) { xr[k] = c[k] + (c[k] - p[hi][k]); }
        double const fr = f(xr); ++ev;
        if (fr < fv[lo]) {
            for (std::size_t k = 0; k < n; ++k) { xe[k] = c[k] + 2.0 * (c[k] - p[hi][k]); }
            double const fe = f(xe); ++ev;
            if (fe < fr) { p[hi] = xe; fv[hi] = fe; } else { p[hi] = xr; fv[hi] = fr; }
        } else if (fr < fv[nh]) {
            p[hi] = xr; fv[hi] = fr;
        } else {
            bool const outside = fr < fv[hi];
            for (std::size_t k = 0; k < n; ++k) { xc[k] = outside ? c[k] + 0.5 * (xr[k] - c[k]) : c[k] + 0.5 * (p[hi][k] - c[k]); }
            double const fc = f(xc); ++ev;
            if (fc < (outside ? fr : fv[hi])) { p[hi] = xc; fv[hi] = fc; }
            else {
                for (std::size_t i = 0; i <= n; ++i) {
                    if (i == lo) { continue; }
                    for (std::size_t k = 0; k < n; ++k) { p[i][k] = p[lo][k] + 0.5 * (p[i][k] - p[lo][k]); }
                    fv[i] = f(p[i]); ++ev;
                }
            }
        }
    }
    std::size_t best = 0;
    for (std::size_t i = 1; i <= n; ++i) { if (fv[i] < fv[best]) { best = i; } }
    return p[best];
}
inline auto fitErrMap(std::vector<double> const& p, std::vector<double> const& y, std::vector<double> const& dlt,
                      bool withPred, GradMutErrModel& m, bool warm) -> bool
{
    std::size_t const n = y.size();
    if (n == 0) { return false; }
    if (withPred) {
        bool vary = false;
        for (std::size_t i = 1; i < n && !vary; ++i) { vary = std::abs(p[i]) != std::abs(p[0]); }
        if (!vary) { withPred = false; }
    }
    auto nlp = [&](std::vector<double> const& th) -> double {
        double const mu = th[0];
        double const b = withPred ? th[1] : 0.0;
        double const ls0 = withPred ? th[2] : th[1];
        double const s1 = withPred ? std::abs(th[3]) : 0.0;
        double const lnu = withPred ? th[4] : th[2];
        double const s0 = std::exp(ls0), nu = std::exp(lnu);
        double const lp = logJeffreysNu(nu);
        if (!std::isfinite(lp) || !(s0 > 0.0) || !std::isfinite(s0)) { return std::numeric_limits<double>::infinity(); }
        double const C = std::lgamma(0.5 * (nu + 1.0)) - std::lgamma(0.5 * nu) - 0.5 * std::log(nu * 3.14159265358979323846);
        double const hp = 0.5 * (nu + 1.0);
        double ll = 0.0;
        double i00 = 0.0, i01 = 0.0, i11 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double const pi = withPred ? p[i] : 0.0;
            double const a = std::abs(pi);
            double const sc = s0 + s1 * a;
            double const z = (y[i] - (pi + mu + b * a)) / sc;
            ll += std::min(C - hp * std::log1p(z * z / nu) - std::log(sc) + std::log(dlt[i]), 0.0);
            double const w = 1.0 / (sc * sc);
            i00 += w; i01 += w * a; i11 += w * a * a;
        }
        double const lps = withPred ? 0.5 * std::log(std::max(i00 * i11 - i01 * i01, 0.0)) : -ls0;
        double const v = -(ll + lp + lps);
        return std::isfinite(v) ? v : std::numeric_limits<double>::infinity();
    };
    std::vector<double> th, st;
    double s0s = 1.0;
    if (warm) {
        std::vector<double> tw = withPred ? std::vector<double>{m.mu, m.b, std::log(m.s0), m.s1, std::log(m.nu)}
                                          : std::vector<double>{m.mu, std::log(m.s0), std::log(m.nu)};
        warm = std::isfinite(nlp(tw));
    }
    if (warm) {
        s0s = m.s0;
        th = withPred ? std::vector<double>{m.mu, m.b, std::log(m.s0), m.s1, std::log(m.nu)}
                      : std::vector<double>{m.mu, std::log(m.s0), std::log(m.nu)};
    } else {
        std::vector<double> d(n);
        for (std::size_t i = 0; i < n; ++i) { d[i] = y[i] - (withPred ? p[i] : 0.0); }
        auto med = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            std::size_t const k = v.size() / 2;
            return (v.size() % 2) ? v[k] : 0.5 * (v[k - 1] + v[k]);
        };
        double const md = med(d);
        std::vector<double> ad(n);
        for (std::size_t i = 0; i < n; ++i) { ad[i] = std::abs(d[i] - md); }
        s0s = std::max(1.4826 * med(ad), 1e-12 * (std::abs(md) + 1.0));
        th = withPred ? std::vector<double>{md, 0.0, std::log(s0s), 0.0, 0.0} : std::vector<double>{md, std::log(s0s), 0.0};
    }
    st = withPred ? std::vector<double>{s0s, 0.5, 1.0, 0.5, 1.0} : std::vector<double>{s0s, 1.0, 1.0};
    if (warm) {
        for (auto& v : st) { v /= std::sqrt(static_cast<double>(n)); }
    }
    auto best = nelderMead(nlp, th, st, 400 * static_cast<int>(th.size()));
    double fb = nlp(best);
    for (int rs = 0; rs < 20 && std::isfinite(fb); ++rs) {
        if (withPred) { st[0] = std::max(std::exp(best[2]), 1e-12); } else { st[0] = std::max(std::exp(best[1]), 1e-12); }
        auto nb = nelderMead(nlp, best, st, 400 * static_cast<int>(th.size()));
        double const fn = nlp(nb);
        bool const better = fn < fb - 1e-9 * (std::abs(fb) + 1e-300);
        if (fn < fb) { best = nb; fb = fn; }
        if (!better) { break; }
    }
    if (!std::isfinite(fb)) { return false; }
    {
        double const mu = best[0];
        double const b = withPred ? best[1] : 0.0;
        double const s0 = std::exp(withPred ? best[2] : best[1]);
        double const s1 = withPred ? std::abs(best[3]) : 0.0;
        double const nu = std::exp(withPred ? best[4] : best[2]);
        double const C = std::lgamma(0.5 * (nu + 1.0)) - std::lgamma(0.5 * nu) - 0.5 * std::log(nu * 3.14159265358979323846);
        for (std::size_t i = 0; i < n; ++i) {
            double const pi = withPred ? p[i] : 0.0;
            double const a = std::abs(pi);
            double const sc = s0 + s1 * a;
            double const z = (y[i] - (pi + mu + b * a)) / sc;
            if (C - 0.5 * (nu + 1.0) * std::log1p(z * z / nu) - std::log(sc) + std::log(dlt[i]) >= 0.0) { return false; }
        }
    }
    m.mu = best[0];
    if (withPred) { m.b = best[1]; m.s0 = std::exp(best[2]); m.s1 = std::abs(best[3]); m.nu = std::exp(best[4]); }
    else { m.b = 0.0; m.s0 = std::exp(best[1]); m.s1 = 0.0; m.nu = std::exp(best[2]); }
    return true;
}
inline auto medianOf(std::vector<double> v) -> double
{
    std::sort(v.begin(), v.end());
    std::size_t const k = v.size() / 2;
    return (v.size() % 2) ? v[k] : 0.5 * (v[k - 1] + v[k]);
}
inline void medmad(std::vector<double> const& p, std::vector<double> const& y, bool withPred,
                   double& med, double& mad)
{
    std::size_t const n = y.size();
    std::vector<double> d(n);
    for (std::size_t i = 0; i < n; ++i) { d[i] = y[i] - (withPred ? p[i] : 0.0); }
    med = medianOf(d);
    std::vector<double> ad(n);
    for (std::size_t i = 0; i < n; ++i) { ad[i] = std::abs(d[i] - med); }
    mad = medianOf(ad) * 1.4826;
}
inline auto betaBinomialTypeII(std::vector<double> const& nAll, std::vector<double> const& x0) -> std::pair<double, double>
{
    int kindsWithData = 0;
    for (double n : nAll) { if (n > 0.0) { ++kindsWithData; } }
    if (kindsWithData < 2) { return {0.5, 0.5}; }
    auto negMarginal = [&](std::vector<double> const& th) -> double {
        double const aa = std::exp(th[0]), bb = std::exp(th[1]);
        if (!(aa > 0.0 && bb > 0.0 && std::isfinite(aa) && std::isfinite(bb))) {
            return std::numeric_limits<double>::infinity();
        }
        double ll = 0.0;
        for (std::size_t m = 0; m < nAll.size(); ++m) {
            if (nAll[m] <= 0.0) { continue; }
            ll += std::lgamma(nAll[m] + 1.0) - std::lgamma(x0[m] + 1.0) - std::lgamma(nAll[m] - x0[m] + 1.0)
                + std::lgamma(aa + bb) - std::lgamma(aa) - std::lgamma(bb)
                + std::lgamma(x0[m] + aa) + std::lgamma(nAll[m] - x0[m] + bb) - std::lgamma(nAll[m] + aa + bb);
        }
        return std::isfinite(ll) ? -ll : std::numeric_limits<double>::infinity();
    };
    auto best = nelderMead(negMarginal, {std::log(0.5), std::log(0.5)}, {1.0, 1.0}, 4000);
    double const aa = std::exp(best[0]), bb = std::exp(best[1]);
    if (aa > 1e-6 && aa < 1e6 && bb > 1e-6 && bb < 1e6 && std::isfinite(negMarginal(best))) { return {aa, bb}; }
    return {0.5, 0.5};
}
inline auto dlPool(std::vector<double> const& est, std::vector<double> const& v, std::vector<char> const& has)
    -> std::pair<std::vector<double>, double>
{
    double sw = 0.0, swx = 0.0;
    int k = 0;
    for (std::size_t i = 0; i < has.size(); ++i) {
        if (has[i] != 0) { ++k; sw += 1.0 / v[i]; swx += est[i] / v[i]; }
    }
    double const fixed = sw > 0.0 ? swx / sw : 0.0;
    double tau2 = 0.0;
    if (k >= 2) {
        double q = 0.0, sw2 = 0.0;
        for (std::size_t i = 0; i < has.size(); ++i) {
            if (has[i] == 0) { continue; }
            double const wi = 1.0 / v[i];
            q += wi * (est[i] - fixed) * (est[i] - fixed);
            sw2 += wi * wi;
        }
        double const denom = sw - sw2 / sw;
        tau2 = (denom > 0.0) ? std::max(0.0, (q - (k - 1)) / denom) : 0.0;
    }
    double swp = 0.0, swxp = 0.0;
    for (std::size_t i = 0; i < has.size(); ++i) {
        if (has[i] == 0) { continue; }
        double const prec = 1.0 / (tau2 + v[i]);
        swp += prec; swxp += est[i] * prec;
    }
    double const poolOut = swp > 0.0 ? swxp / swp : 0.0;
    std::vector<double> out(has.size(), poolOut);
    for (std::size_t i = 0; i < has.size(); ++i) {
        if (has[i] == 0) { continue; }
        double const w = tau2 / (tau2 + v[i]);
        out[i] = w * est[i] + (1.0 - w) * poolOut;
    }
    return {out, poolOut};
}
}

struct MultitreeConfig {
    std::size_t populationSize = 1000;
    std::size_t generations = 50;
    std::size_t maxLength = 24;
    std::size_t maxDepth = 12;
    bool enforceStrictMaxLength = false;
    std::size_t gridSize = 200;
    double optFraction = 0.15;
    bool lmTargeted = true;
    bool useCascadeCull = false;
    bool useLiveCoarseSieve = false;
    double liveCoarseSieveThreshold = 3.379763530657003;
    double coeffSpread = 8.0;
    std::uint64_t seed = 42;
    bool verbose = false;
    bool useRotationGauge = false;
    bool useParametrizationGauge = false;
    bool useQuotientProjection = true;
    bool useQuotientFiniteMin = true;
    BracketNewtonRuntime* bracketNewtonRuntime = nullptr;
    BracketNewtonConfig bracketNewtonConfig = [] {
        BracketNewtonConfig c;
        c.useNoiseWeightedNewton = true;
        return c;
    }();
    bool useEngineAmplitudeRuntime = false;
    std::function<BracketNewtonGeneration(PairIndividual const&)> bracketNewtonGenerationBuilder{};
    std::function<void(PairIndividual&, BracketNewtonResult const&)> bracketNewtonResultConsumer{};
    std::function<BracketNewtonParentPacingMap(std::vector<PairIndividual> const&,
                                                std::vector<PairIndividual> const&)> bracketNewtonParentMapBuilder{};
    UmaBracketNewtonRunner* bracketNewtonUmaRunner = nullptr;
    std::function<bool(BracketNewtonConfig const&, BracketNewtonGeneration const&,
                       std::function<void(BracketNewtonResult)>)> bracketNewtonAsyncSubmit{};
    std::function<void()> bracketNewtonAsyncDrain{};
    std::function<void()> bracketNewtonAsyncSkip{};
    Operon::Scalar floorStop = 0;
    bool floorOnChamfer = false;
    Operon::Scalar arcLengthGuardRatio = 0;
    Operon::Scalar closedGap = 0.35;
    std::size_t alpsLayers = 0;
    std::size_t ageGap = 5;
    Operon::Scalar crossTreeProb = 0.3;
    Operon::Scalar orderWeight = 1.0;
    bool whitenOrderTerm = true;
    bool useArcOrderDistance = false;
    bool useBandedOrderDistance = true;
    int orderBandWidth = 8;
    bool useMonotonicityOrderGuard = true;
    bool useLm = false;
    int constOptimizer = 1;
    int lmBlockNn = 1;
    bool useLmBytecodeCache = true;
    std::size_t lmDownsampleK = 0;
    int lmDownsamplePolishIters = 2;
    double annealResolution = 0.75;
    std::size_t annealLowN = 64;
    double annealStage2Resolution = 0.0;
    std::size_t annealMidN = 0;
    std::size_t annealLowM = 64;
    std::size_t annealMidM = 0;
    bool semanticCrossover = false;
    bool useLexicase = true;
    bool useCoherentSelection = false;
    bool useNoiseWeightedLexicase = true;
    bool useCoherentLexicaseResidual = false;
    bool useResidualRetention = true;
    bool useResidualNormalRetention = false;
    bool useNoiseGuidedMutation = true;
    bool useSpectralStructuralMutation = false;
    bool useJacobianNodeSensitivity = false;
    int gradientMutationMode = 1;
    std::array<GradMutErrModel, 7> gradientMutationErr{};
    GradMutErrModel gradientMutationRandomArm{};
    bool gradientMutationUseErrModel = true;
    bool gradientMutationInRunErr = true;
    bool gradientMutationG6e = true;
    bool gradientMutationErrcalPImprove = false;
    bool gradientMutationBanditRefitReward = false;
    bool gradientMutationLog = false;
    std::size_t gradientMutationLogTopK = 5;
    bool gradientMutationEffLog = false;
    bool gradientMutationTrajectory = false;
    bool binaryInsertSubtree = false;
    bool mutationOpDiag = false;
    std::size_t tournamentSize = 7;
    bool useIslands = false;
    int numIslands = 1;
    int migrationInterval = 10;
    int migrationCount = 1;
    bool perTreeCredit = false;
    Operon::Scalar floorClamp = 0;
    std::size_t floorSimplifyGens = 0;
    double maxFitSeconds = 0;
    std::size_t constMultiStart = 1;
    bool constMultiStartDataSeed = true;
    bool polishFullObjective = false;
    bool guardConstOpt = false;
    bool rejectNonFiniteBeforeConstOpt = false;
    bool useArchive = false;
    Operon::Scalar curvatureWeight = 0;
    bool curvatureMultiplicative = false;
    Operon::Scalar curvatureFloorTurns = 0.75;
    Operon::Scalar curvatureFloorRel = 0.10;
    int curvatureHarmonics = 8;
    std::size_t nThreads = 1;

    bool useGridSdf = false;
    int gridSdfResolution = 128;

    bool useSymmetricFitness = false;

    bool useAccelerateTranscendentals = false;

    bool useBytecodeEvaluator = true;

    bool useRegisterInterp = (regEvalModeFromEnv() & 1) != 0;

    bool useBatchedRegisterEval = (regEvalModeFromEnv() & 2) != 0;

    bool useBatchedLmInnerLoop = false;

    bool useConcurrentCpuGpuLM = false;
    double concurrentLmForceGpuShare = -1.0;
    bool concurrentLmParityShadow = false;
    int concurrentLmKernelVariant = -1;
    int concurrentLmSchedule = 1;
    int concurrentLmKeepWarmIters = 0;
    int concurrentLmRaceGpuMax = 0;

    bool profileTiming = false;

    bool useFastNonDominatedSort = true;

    bool useExactNlogNNondominatedSort = true;

    bool pointClamp = false;

    int inLoopConstIter = 0;
    int polishConstIter = 0;

    bool useAmxChamfer = false;

    bool useCisFitness = true;

    bool useSiCisFitness = true;
    int cisGrid48Mode = cisGrid48ModeFromEnv();
    bool useB3Grid48Gpu = b3Grid48GpuFromEnv() != 0;
    bool useCoarseFourPointSieve = false;
    double coarseFourPointThreshold = std::numeric_limits<double>::infinity();
    double robustContinuousWeight = 0.0;
    double robustContinuousTangentWeight = 0.05;
    double robustContinuousCurvatureWeight = 0.01;

    bool useGlobalRobustNllFitness = false;
    double globalRobustNllNu = 2.4;

    bool useExternalPointSigma2 = false;

#if defined(__APPLE__) && defined(__OBJC__)
    bool useBatchedGpuPointSeg = true;
#else
    bool useBatchedGpuPointSeg = false;
#endif

    bool useParallelBreed = useParallelBreedFromEnv();

    bool useGpuReverseHalf = useGpuReverseHalfFromEnv();

    bool usePersistentGpuScoring = true;

    int oneResidualMode = oneResidualModeFromEnv();

    bool gpuPointSegZeroCopy = false;
    bool gpuPointSegOverlap = false;
    bool gpuPointSegUnretained = false;
    bool gpuPointSegSpinWait = false;
    bool gpuPointSegResidencySet = false;
    int gpuPointSegOverlapSchedule = 0;
    bool gpuPointSegOverlapHelper = false;

    bool legacyFrontRecommendation = false;

    bool useChowAbstention = false;
    bool reportNmlPick = false;
    bool nmlOutOfRangeInadmissible = true;
    bool nmlCanonicalStructCode = true;
    bool nmlPolishMembers = true;
    double nmlAbsSigma2 = std::numeric_limits<double>::quiet_NaN();
    bool computeStrictGateTelemetry = false;
    float chowAbstentionCost = 1.0F / 3.0F;

    bool useOperatorFamilyBias = true;

    bool useComplexityHallOfFame = false;

    bool useFrequencyDiscovery = false;
    bool useVarProFrequencyGn = false;
    bool useWeightedVarProFrequencyGn = false;
    bool useNoiseWeightedFrequencySeed = true;
    bool useMatrixPencilFrequencySeed = false;

    bool useTopologyGatedFrequency = false;

    bool useREscalationLadder = false;

    bool constMagnitudeComplexity = false;

    bool usePriorMinimalComplexity = false;

    double championAdequacyBand = 0.0;

    bool directFinalistEval = true;

    bool useVarProLinear = useVarProLinearFromEnv();
};

inline auto populationFloor(std::vector<PairIndividual> const& pop) -> Operon::Scalar
{
    Operon::Scalar worst = -std::numeric_limits<Operon::Scalar>::infinity();
    bool any = false;
    for (auto const& ind : pop) {
        if (std::isfinite(ind.shapeErr) && ind.shapeErr > worst) { worst = ind.shapeErr; any = true; }
    }
    return any ? worst : std::numeric_limits<Operon::Scalar>::infinity();
}

inline auto lmRefinementOrder(std::vector<PairIndividual> const& offspring, std::size_t k,
                              bool enableNoiseRouting) -> std::vector<std::size_t>
{
    std::vector<std::size_t> order(offspring.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    k = std::min(k, order.size());
    auto const supported = [&](PairIndividual const& ind) {
        return enableNoiseRouting && std::any_of(ind.regionErrCache.begin(), ind.regionErrCache.end(),
            [](Operon::Scalar e) { return std::isfinite(e) && e > Operon::Scalar{0}; });
    };
    std::partial_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(k), order.end(),
        [&](std::size_t lhs, std::size_t rhs) {
            auto const ls = supported(offspring[lhs]);
            auto const rs = supported(offspring[rhs]);
            if (ls != rs) { return ls; }
            return offspring[lhs].shapeErr < offspring[rhs].shapeErr;
        });
    return order;
}

struct FrontEntry {
    std::string exprX, exprY;
    Operon::Scalar shapeErr;
    int nodes;
    std::vector<Operon::Scalar> regionErr;
    std::vector<Operon::Scalar> curveX, curveY;
    Operon::Scalar orderErr{0};
    Operon::Scalar errX{0}, errY{0};
    double turning{0};
    double nmlAbsTotal{std::numeric_limits<double>::quiet_NaN()};
};

struct ConcurrentLmGenStat {
    int gen{0};
    std::size_t k{0};
    std::size_t nGpuTarget{0};
    std::size_t nGpu{0};
    std::size_t nCpu{0};
    std::size_t nFallback{0};
    double gpuShare{0.0};
    double prepMs{0.0};
    double cpuMs{0.0};
    double gpuMs{0.0};
    double gpuKernelMs{0.0};
    double idleCpuMs{0.0};
    double idleGpuMs{0.0};
    double postMs{0.0};
    double stageMs{0.0};
    double rateCpu{0.0};
    double rateGpu{0.0};
    std::size_t nGpuSubmitted{0};
    bool gpuUsed{false};
    double discardMs{0.0};
    double encodeMs{0.0}, queueMs{0.0}, execMs{0.0}, notifyMs{0.0};
    double gpuReadyMs{0.0};
};
struct ConcurrentLmParityGen {
    int gen{0};
    std::vector<float> mixed, cpuRef;
    std::vector<float> cpuUlp;
    std::vector<char> onGpu, fellBack;
};

struct MultitreeResult {
    std::vector<FrontEntry> front;
    std::vector<FrontEntry> hallOfFame;
    Operon::Tree x, y;
    std::string exprX, exprY;
    Operon::Scalar shapeError{std::numeric_limits<Operon::Scalar>::max()};
    int generations{0};
    bool truncated{false};
    std::size_t nEvalErrors{0};
    std::size_t nConstDescents{0};
    std::size_t nConstLossCalls{0};
    std::size_t nConstIters{0};
    std::size_t nLmStopUnresolvableDamped{0};
    std::size_t nLmStopRealizedSmall{0};
    std::size_t nLmStopFixedPoint{0};
    std::size_t nCascadeCullGens{0};
    std::size_t nCascadeCullSeen{0};
    std::size_t nCascadeCullSurvived{0};
    std::size_t nLiveCoarseSieveSeen{0};
    std::size_t nLiveCoarseSieveCulled{0};
    double liveCoarseSieveBoundMin{std::numeric_limits<double>::infinity()};
    double liveCoarseSieveBoundMax{-std::numeric_limits<double>::infinity()};
    double rotationGaugeAngle{0.0};
    ParametrizationGauge parametrizationGauge{};
    std::vector<int> parametrizationPerm{};
    double liveCoarseSieveBoundMean{std::numeric_limits<double>::quiet_NaN()};
    std::size_t nNonFiniteFitness{0};
    std::size_t nMutationOversizeRetries{0};
    std::size_t nMutationFallbackToParent{0};
    std::size_t nRejectedNonFiniteBeforeConstOpt{0};
    std::size_t nLmRecoveredFromNonFiniteStart{0};
    std::size_t nNonFiniteWithLog{0};
    std::size_t nNonFiniteWithSqrt{0};
    std::size_t nNonFiniteWithDiv{0};
    std::size_t nNonFiniteWithExp{0};
    std::size_t nNonFiniteWithNoneOfThese{0};
    double dataTurning{0};
    std::size_t simplifyGensRan{0};
    int annealSwitchGen{-1};
    double precomputeTimeMs{0.0};

    double profEvalMs{0.0};
    double profLmMs{0.0};
    double profScoreMs{0.0};
    double profBreedMs{0.0};
    double profSelectMs{0.0};
    double profGenMs{0.0};

    double profSelectCombineMs{0.0};
    double profSelectNdsMs{0.0};
    double profSelectCrowdMs{0.0};

    double profScoreDistMs{0.0};
    double profScoreCisMs{0.0};
    double profScoreOrderMs{0.0};
    double profScoreCurvMs{0.0};

    double profLmBuildMs{0.0};
    double profLmLossMs{0.0};

    double profCisPointSegMs{0.0};
    double profCisNllMs{0.0};
    double profCisReverseMs{0.0};
    double profCisGrid48Ms{0.0};

    double profGpuPrepassMs{0.0};

    double profParallelWallMs{0.0};

    double chowAbstainLeadBits{std::numeric_limits<double>::quiet_NaN()};
    double chowAbstainThresholdBits{std::numeric_limits<double>::quiet_NaN()};
    bool chowShouldAnswer{true};

    std::vector<Sketch::nml::Terms> nmlFront;
    int nmlPickIndex{-1};
    int nmlPickIndexAbs{-1};
    int nmlRecommendedIndex{-1};
    double nmlComputeMs{0.0};
    std::size_t nmlNMembers{0};
    std::size_t nmlNDistinct{0};
    double nmlFreqRange{std::numeric_limits<double>::quiet_NaN()};
    double nmlExpRateRange{std::numeric_limits<double>::quiet_NaN()};
    double nmlExtent{std::numeric_limits<double>::quiet_NaN()};
    std::vector<Sketch::nml::Terms> nmlFrontRaw;
    std::vector<std::string> nmlPolishedExprX;
    std::vector<std::string> nmlPolishedExprY;
    double nmlGridPeriod{std::numeric_limits<double>::quiet_NaN()};
    double nmlPolishMs{0.0};
    double nmlSigmaBar{std::numeric_limits<double>::quiet_NaN()};

    Operon::Scalar siCisRecommendedShapeNll{std::numeric_limits<Operon::Scalar>::quiet_NaN()};

    Operon::Scalar nullLineStrictShapeNll{std::numeric_limits<Operon::Scalar>::quiet_NaN()};
    std::vector<ConcurrentLmGenStat> cclmStats;
    std::vector<ConcurrentLmParityGen> cclmParity;
    std::size_t cclmGpuAllocEvents{0};
};

template <typename> struct DtwPathTestHook;
template <typename> struct QuotientFiniteMinTestHook;
template <typename> struct GradMutTestHook;

struct GradMutRecord {
    int gen{0};
    int event{0};
    int rankPred{0};
    int kind{0};
    int node{0};
    int treeLen{0};
    int nCands{0};
    bool isX{true};
    bool chosen{false};
    bool childIsParentClone{false};
    double predFixed{0}, predRefitTotal{0}, predRefitMut{0};
    double realFixed{0}, realRefitTotal{0}, realRefitMut{0};
    double predRefitMutNats{0}, realRefitMutNats{0};
    double L0{0}, Lref0{0};
    bool randomChosen{false};
    int dK{0};
    int randomLen{0};
    double predNml{0}, dStruct{0}, dParam{0};
    double realNmlMut{0};
    double randomRealNmlMut{0};
    double predFixedNats{0}, realFixedNats{0};
};

struct GradMutEff {
    int gen{0};
    int type{0};
    int kind{-1};
    int arm{-1};
    int lenDelta{0};
    double pred{0};
    double dVsPre{0};
    double dVsParent{0};
};

struct BinaryInsertSubtreeMutation final : public Operon::MutatorBase {
    BinaryInsertSubtreeMutation(Operon::CreatorBase const* creator, Operon::CoefficientInitializerBase const* coeffInit,
                                std::size_t maxDepth, std::size_t maxLength)
        : creator_(creator), coeffInit_(coeffInit), maxDepth_(maxDepth), maxLength_(maxLength) {}

    auto operator()(Operon::RandomGenerator& random, Operon::Tree tree) const -> Operon::Tree override
    {
        auto const& nodes = tree.Nodes();
        if (nodes.size() + 1 >= maxLength_) { return tree; }
        auto test = [](Operon::Node const& n) {
            return n.template IsOp<Operon::BuiltinOp::Add, Operon::BuiltinOp::Mul, Operon::BuiltinOp::Sub, Operon::BuiltinOp::Div>();
        };
        auto n = std::count_if(nodes.begin(), nodes.end(), test);
        if (n == 0) { return tree; }
        auto index = std::uniform_int_distribution<decltype(n)>(1, n)(random);
        std::size_t i = 0;
        for (; i < nodes.size(); ++i) {
            if (test(nodes[i]) && --index == 0) { break; }
        }
        std::size_t const availableLength = maxLength_ - nodes.size() - 1;
        std::size_t const availableDepth = std::max<std::size_t>(tree.Depth(), maxDepth_) - nodes[i].Level;
        auto const newLen = std::uniform_int_distribution<std::size_t>(1, availableLength)(random);
        auto subtree = (*creator_)(random, newLen, 1, std::max<std::size_t>(availableDepth, 1));
        (*coeffInit_)(random, subtree);

        Operon::Vector<Operon::Node> out;
        out.reserve(nodes.size() + subtree.Length() + 1);
        using Signed = std::make_signed_t<std::size_t>;
        auto const start = static_cast<Signed>(i - nodes[i].Length);
        std::copy(nodes.begin(), nodes.begin() + start, std::back_inserter(out));
        std::copy(subtree.Nodes().begin(), subtree.Nodes().end(), std::back_inserter(out));
        std::copy(nodes.begin() + start, nodes.begin() + static_cast<Signed>(i + 1), std::back_inserter(out));
        Operon::Node wrap = nodes[i];
        wrap.Arity = 2;
        out.push_back(wrap);
        std::copy(nodes.begin() + static_cast<Signed>(i + 1), nodes.end(), std::back_inserter(out));
        return Operon::Tree(out).UpdateNodes();
    }

private:
    Operon::CreatorBase const* creator_;
    Operon::CoefficientInitializerBase const* coeffInit_;
    std::size_t maxDepth_;
    std::size_t maxLength_;
};

enum OpDiagId : int {
    OpCrossoverX = 0,
    OpCrossoverY = 1,
    OpCrossTreeExchange = 2,
    OpMutOnePoint = 3,
    OpMutChangeFunc = 4,
    OpMutReplaceSub = 5,
    OpMutRemoveSub = 6,
    OpMutInsertSub = 7,
    OpMutNoiseGuided = 8,
    OpMutGradient = 9,
    OpDiagCount = 10
};
inline constexpr std::array<char const*, OpDiagCount> kOpDiagNames{
    "CrossoverX", "CrossoverY", "CrossTreeExchange",
    "Mut_OnePoint", "Mut_ChangeFunc", "Mut_ReplaceSub", "Mut_RemoveSub", "Mut_InsertSub",
    "Mut_NoiseGuided", "Mut_Gradient",
};
inline constexpr std::array<int, 5> kMutSlotToOpDiagId{
    OpMutOnePoint, OpMutChangeFunc, OpMutReplaceSub, OpMutRemoveSub, OpMutInsertSub,
};
struct OpDiagCounters {
    std::array<std::atomic<long long>, OpDiagCount> applications{};
    std::array<std::atomic<long long>, OpDiagCount> noop{};
    std::array<std::atomic<long long>, OpDiagCount> invalid{};
};
inline auto opDiagCounters() -> OpDiagCounters&
{
    static OpDiagCounters c;
    return c;
}
inline auto opDiagEnabled() -> bool
{
    static bool const v = std::getenv("SKETCH_OP_DIAG") != nullptr;
    return v;
}
inline void opDiagRecord(int id, bool noop, bool invalid)
{
    auto& c = opDiagCounters();
    c.applications[static_cast<std::size_t>(id)].fetch_add(1, std::memory_order_relaxed);
    if (noop) { c.noop[static_cast<std::size_t>(id)].fetch_add(1, std::memory_order_relaxed); }
    if (invalid) { c.invalid[static_cast<std::size_t>(id)].fetch_add(1, std::memory_order_relaxed); }
}
inline auto sameTreeStructure(Operon::Tree const& x, Operon::Tree const& y) -> bool
{
    auto const& a = x.Nodes();
    auto const& b = y.Nodes();
    if (a.size() != b.size()) { return false; }
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (a[k].HashValue != b[k].HashValue || a[k].Arity != b[k].Arity
            || std::memcmp(&a[k].Value, &b[k].Value, sizeof(a[k].Value)) != 0) { return false; }
    }
    return true;
}
struct NSafetyCounters {
    std::atomic<long long> nonFiniteMutatedConst{0};
    std::atomic<long long> nonFiniteLmConst{0};
};
inline auto nSafetyCounters() -> NSafetyCounters&
{
    static NSafetyCounters c;
    return c;
}

inline auto k1cMeasureMutex() -> std::mutex&
{
    static std::mutex m;
    return m;
}
inline auto k1cMeasureFile() -> std::ofstream*
{
    static std::ofstream* f = [] () -> std::ofstream* {
        char const* path = std::getenv("SKETCH_MEASURE_K1C");
        if (path == nullptr) { return nullptr; }
        auto* out = new std::ofstream(path, std::ios::app);
        return (out != nullptr && out->good()) ? out : nullptr;
    }();
    return f;
}

struct Op2GaugeCounters {
    static constexpr int kMaxNullDim = 16;
    std::array<std::atomic<long long>, kMaxNullDim + 1> nullDimHist{};
};
inline auto op2GaugeCounters() -> Op2GaugeCounters&
{
    static Op2GaugeCounters c;
    return c;
}
inline void op2GaugeRecordNullDim(int d)
{
    auto& c = op2GaugeCounters();
    std::size_t const idx = static_cast<std::size_t>(std::clamp(d, 0, Op2GaugeCounters::kMaxNullDim));
    c.nullDimHist[idx].fetch_add(1, std::memory_order_relaxed);
}
inline auto op2MaxAbsConstBitsRef() -> std::atomic<std::uint64_t>&
{
    static std::atomic<std::uint64_t> m{0};
    return m;
}
inline void op2GaugeRecordAbsConst(double v)
{
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof(bits));
    auto& ref = op2MaxAbsConstBitsRef();
    std::uint64_t cur = ref.load(std::memory_order_relaxed);
    while (bits > cur && !ref.compare_exchange_weak(cur, bits, std::memory_order_relaxed)) {}
}
inline auto op2MaxAbsConst() -> double
{
    std::uint64_t const bits = op2MaxAbsConstBitsRef().load(std::memory_order_relaxed);
    double v;
    std::memcpy(&v, &bits, sizeof(v));
    return v;
}
struct Op2ValidCounters {
    std::array<std::atomic<long long>, 7> events{};
    std::array<std::atomic<long long>, 7> draws{};
};
inline auto op2ValidCounters() -> Op2ValidCounters&
{
    static Op2ValidCounters c;
    return c;
}
inline constexpr std::array<char const*, 7> kOp2ValidNames{
    "OnePoint", "ChangeFunc", "ReplaceSub", "RemoveSub", "InsertSub", "CrossoverX", "CrossoverY",
};
inline void printOp2DiagTable(std::FILE* out)
{
    auto& vc = op2ValidCounters();
    std::fprintf(out, "[op2-valid] operator          events     draws  meanDraws/event\n");
    for (int i = 0; i < 7; ++i) {
        long long const e = vc.events[static_cast<std::size_t>(i)].load(std::memory_order_relaxed);
        long long const d = vc.draws[static_cast<std::size_t>(i)].load(std::memory_order_relaxed);
        double const mean = e > 0 ? static_cast<double>(d) / static_cast<double>(e) : 0.0;
        std::fprintf(out, "[op2-valid] %-15s %9lld %9lld %10.4f\n", kOp2ValidNames[static_cast<std::size_t>(i)], e, d, mean);
    }
    auto& gc = op2GaugeCounters();
    std::fprintf(out, "[op2-gauge] maxAbsConst=%.9g nullDimHist(0..%d+)=",
                 op2MaxAbsConst(), Op2GaugeCounters::kMaxNullDim);
    for (int i = 0; i <= Op2GaugeCounters::kMaxNullDim; ++i) {
        std::fprintf(out, "%lld ", gc.nullDimHist[static_cast<std::size_t>(i)].load(std::memory_order_relaxed));
    }
    std::fprintf(out, "\n");
}

inline void printOpDiagTable(std::FILE* out)
{
    auto& c = opDiagCounters();
    std::fprintf(out, "[op-diag] operator                apps      noop  noop%%   invalid  invalid%%\n");
    for (int i = 0; i < OpDiagCount; ++i) {
        long long const a = c.applications[static_cast<std::size_t>(i)].load(std::memory_order_relaxed);
        long long const n = c.noop[static_cast<std::size_t>(i)].load(std::memory_order_relaxed);
        long long const iv = c.invalid[static_cast<std::size_t>(i)].load(std::memory_order_relaxed);
        double const noopPct = a > 0 ? (100.0 * static_cast<double>(n) / static_cast<double>(a)) : 0.0;
        double const invPct = a > 0 ? (100.0 * static_cast<double>(iv) / static_cast<double>(a)) : 0.0;
        std::fprintf(out, "[op-diag] %-20s %9lld %9lld %6.2f %9lld %8.2f\n",
                     kOpDiagNames[static_cast<std::size_t>(i)], a, n, noopPct, iv, invPct);
    }
}

struct MutOpDiagCounters {
    std::array<std::atomic<std::size_t>, 5> calls{};
    std::array<std::atomic<std::size_t>, 5> noop{};
    std::atomic<std::size_t> mutateCalls{0};
    std::atomic<std::size_t> insertNoopAtCap{0};
    std::size_t maxLength{0};
};
inline thread_local int tlLastMutOp = -1;
struct CountingMutation final : public Operon::MutatorBase {
    CountingMutation(Operon::MutatorBase const* inner, int op, MutOpDiagCounters* ctr) : inner_(inner), op_(op), ctr_(ctr) {}
    auto operator()(Operon::RandomGenerator& random, Operon::Tree tree) const -> Operon::Tree override
    {
        Operon::Tree in = tree;
        auto out = (*inner_)(random, std::move(tree));
        ctr_->calls[static_cast<std::size_t>(op_)].fetch_add(1, std::memory_order_relaxed);
        auto const& a = in.Nodes();
        auto const& b = out.Nodes();
        bool same = a.size() == b.size();
        for (std::size_t k = 0; same && k < a.size(); ++k) {
            same = a[k].HashValue == b[k].HashValue && a[k].Arity == b[k].Arity
                && std::memcmp(&a[k].Value, &b[k].Value, sizeof(a[k].Value)) == 0;
        }
        if (same) {
            ctr_->noop[static_cast<std::size_t>(op_)].fetch_add(1, std::memory_order_relaxed);
            if (op_ == 4 && a.size() + 1 >= ctr_->maxLength) { ctr_->insertNoopAtCap.fetch_add(1, std::memory_order_relaxed); }
        }
        tlLastMutOp = op_;
        if (opDiagEnabled()) {
            opDiagRecord(kMutSlotToOpDiagId[static_cast<std::size_t>(op_)], same, false);
        }
        return out;
    }
private:
    Operon::MutatorBase const* inner_;
    int op_;
    MutOpDiagCounters* ctr_;
};

#if defined(__APPLE__)
template<typename T>
inline void RegisterAccelerateTranscendentals(Operon::DispatchTable<T>& dt) {
    constexpr auto S = Operon::Backend::BatchSize<T>;
    dt.template RegisterFunction<T>(
        static_cast<Operon::Hash>(Operon::BuiltinOp::Sin),
        [](Operon::Vector<Operon::Node> const& nodes, Operon::Backend::View<T, S> view, size_t i, Operon::Range ) {
            auto const w = nodes[i].Value;
            auto* res = Operon::Backend::Ptr<T, S>(view, i);
            auto const* arg = Operon::Backend::Ptr<T, S>(view, i - 1);
            int n = static_cast<int>(S);
            vvsinf(res, arg, &n);
            if (w != 1.0f) {
                for (int k = 0; k < n; ++k) res[k] *= w;
            }
        },
        Operon::Dispatch::MakeDiffCall<Operon::BuiltinOp::Sin, T, S>()
    );
    dt.template RegisterFunction<T>(
        static_cast<Operon::Hash>(Operon::BuiltinOp::Cos),
        [](Operon::Vector<Operon::Node> const& nodes, Operon::Backend::View<T, S> view, size_t i, Operon::Range ) {
            auto const w = nodes[i].Value;
            auto* res = Operon::Backend::Ptr<T, S>(view, i);
            auto const* arg = Operon::Backend::Ptr<T, S>(view, i - 1);
            int n = static_cast<int>(S);
            vvcosf(res, arg, &n);
            if (w != 1.0f) {
                for (int k = 0; k < n; ++k) res[k] *= w;
            }
        },
        Operon::Dispatch::MakeDiffCall<Operon::BuiltinOp::Cos, T, S>()
    );
}
#endif

#endif
