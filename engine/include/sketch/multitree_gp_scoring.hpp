
    auto symmetricFitnessError(std::vector<Operon::Scalar> const& gx,
                               std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto const M = gx.size();
        auto const N = dataX_.size();
        if (M < 2 || N < 2) { return BIG; }

        thread_local std::vector<Operon::Scalar> sfCx, sfCy;
        sfCx.resize(M);
        sfCy.resize(M);
        {
            Operon::Scalar gx0 = gx[0], gy0 = gy[0];
            if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return BIG; }
            Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
            sfCx[0] = gx0; sfCy[0] = gy0;
            for (std::size_t j = 1; j < M; ++j) {
                Operon::Scalar const vx = gx[j], vy = gy[j];
                if (!std::isfinite(vx) || !std::isfinite(vy)) { return BIG; }
                sfCx[j] = vx; sfCy[j] = vy;
                if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
                if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                sfCx[j] = (sfCx[j] - cxc) / sp;
                sfCy[j] = (sfCy[j] - cyc) / sp;
            }
        }

        thread_local std::vector<Vec2> sfStrokePts, sfCandPts;
        sfStrokePts.resize(N);
        sfCandPts.resize(M);
        for (std::size_t i = 0; i < N; ++i) { sfStrokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < M; ++j) { sfCandPts[j] = Vec2{sfCx[j], sfCy[j]}; }

        auto const sfRes = evaluateSymmetricFitness(sfStrokePts, sfCandPts, gridSdf_);
        if (!sfRes.ok || !std::isfinite(sfRes.symmetric)) { return BIG; }
        return static_cast<Operon::Scalar>(sfRes.symmetric);
    }

    auto oneResBit(int b) const -> bool { return (config_.oneResidualMode & b) != 0; }
    auto oneResOn() const -> bool { return oneResBit(1); }
    auto oneResUnsigned() const -> bool { return oneResidualUnsigned(config_.oneResidualMode); }

    auto cisGrid48Bit0() const -> bool { return (config_.cisGrid48Mode & 1) != 0; }
    auto cisGrid48Bit1() const -> bool { return (config_.cisGrid48Mode & 2) != 0; }
    auto cisGrid48Bit2() const -> bool { return (config_.cisGrid48Mode & 4) != 0; }

    struct Grid48PointSeg {
        PointSegmentResult ps;
        std::size_t strokeN{0};
        double invSc{1.0};
        bool ok{false};
    };
    auto cisGrid48Target() const -> std::vector<Vec2> const&
    {
        return cisGrid48Bit2() ? cisGrid48StrokeData_ : cisGrid48Stroke_;
    }
    auto cisGrid48InvSc() const -> double { return cisGrid48Bit2() ? (1.0 / cisGrid48Scale_) : 1.0; }
    auto cisGrid48PointSeg(std::vector<Vec2> const& candPts, bool computeTau, bool computeNormal) const
        -> Grid48PointSeg
    {
        Grid48PointSeg out;
        bool const useDataFrame = cisGrid48Bit2();
        auto const& stroke = cisGrid48Target();
        out.ps = evaluatePointSegmentOne(stroke, candPts, computeTau, computeNormal, oneResUnsigned());
        out.strokeN = stroke.size();
        out.invSc = cisGrid48InvSc();
        (void)useDataFrame;
        out.ok = out.ps.ok && out.ps.residual.size() == out.strokeN;
        return out;
    }

    auto pointSegOp(std::vector<Vec2> const& strokePts, std::vector<Vec2> const& candPts,
                    bool computeTau, bool computeNormal) const -> PointSegmentResult
    {
        if (oneResOn()) { return evaluatePointSegmentOne(strokePts, candPts, computeTau, computeNormal, oneResUnsigned()); }
        PointSegmentConfig psCfg;
        psCfg.computeTau = computeTau;
        psCfg.computeNormal = computeNormal;
        return evaluatePointSegment(strokePts, candPts, psCfg);
    }

    struct OneResCache {
        bool valid{false};
        void const* owner{nullptr};
        std::size_t N{0};
        std::vector<Operon::Scalar> gx, gy;
        std::vector<Operon::Scalar> residual, tau;
    };
    static inline thread_local OneResCache oneResCache_;
    mutable std::atomic<long long> nOneResCacheHit_{0};
    mutable std::atomic<long long> nOneResCacheMiss_{0};

    void oneResCacheStore(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                          std::size_t N, Operon::Scalar const* res, Operon::Scalar const* tau) const
    {
        auto& c = oneResCache_;
        c.valid = true;
        c.owner = this;
        c.N = N;
        c.gx.assign(gx.begin(), gx.end());
        c.gy.assign(gy.begin(), gy.end());
        c.residual.assign(res, res + N);
        if (tau != nullptr) { c.tau.assign(tau, tau + N); } else { c.tau.clear(); }
    }

    auto oneResCacheLookup(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy) const
        -> OneResCache const*
    {
        auto const& c = oneResCache_;
        if (c.valid && c.owner == this && c.N == dataX_.size()
            && c.gx.size() == gx.size() && c.gy.size() == gy.size()
            && std::equal(gx.begin(), gx.end(), c.gx.begin())
            && std::equal(gy.begin(), gy.end(), c.gy.begin())) {
            nOneResCacheHit_.fetch_add(1, std::memory_order_relaxed);
            return &c;
        }
        nOneResCacheMiss_.fetch_add(1, std::memory_order_relaxed);
        return nullptr;
    }

    auto oneResidualOf(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                       bool wantTau, std::vector<Operon::Scalar>& residual,
                       std::vector<Operon::Scalar>* tau) const -> bool
    {
        auto const N = dataX_.size();
        if (auto const* c = oneResCacheLookup(gx, gy)) {
            if (!wantTau || c->tau.size() == N) {
                residual = c->residual;
                if (wantTau && tau != nullptr) { *tau = c->tau; }
                return true;
            }
        }
        auto const M = gx.size();
        if (M < 2 || N < 2 || gy.size() != M) { return false; }
        thread_local std::vector<Operon::Scalar> oCx, oCy;
        oCx.resize(M);
        oCy.resize(M);
        {
            Operon::Scalar gx0 = gx[0], gy0 = gy[0];
            if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return false; }
            Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
            oCx[0] = gx0; oCy[0] = gy0;
            for (std::size_t j = 1; j < M; ++j) {
                Operon::Scalar const vx = gx[j], vy = gy[j];
                if (!std::isfinite(vx) || !std::isfinite(vy)) { return false; }
                oCx[j] = vx; oCy[j] = vy;
                if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
                if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                oCx[j] = (oCx[j] - cxc) / sp;
                oCy[j] = (oCy[j] - cyc) / sp;
            }
        }
        thread_local std::vector<Vec2> oStroke, oCand;
        oStroke.resize(N);
        oCand.resize(M);
        for (std::size_t i = 0; i < N; ++i) { oStroke[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < M; ++j) { oCand[j] = Vec2{oCx[j], oCy[j]}; }
        thread_local PointSegOneResult<Operon::Scalar> r;
        evaluatePointSegmentOneInto<Operon::Scalar>(oStroke.data(), N, oCand.data(), M, wantTau, false, r, oneResUnsigned());
        if (!r.ok || r.residual.size() != N) { return false; }
        residual = r.residual;
        if (wantTau && tau != nullptr) { *tau = r.tau; }
        return true;
    }

    auto oneResReverseHalfSum(std::vector<Operon::Scalar> const& candX, std::vector<Operon::Scalar> const& candY,
                           std::size_t M) const -> Operon::Scalar
    {
        auto const N = dataX_.size();
        thread_local std::vector<Vec2> drawn;
        drawn.resize(N);
        for (std::size_t i = 0; i < N; ++i) { drawn[originalArcIndex(i)] = Vec2{dataX_[i], dataY_[i]}; }
        auto const K = noise_.size();
        auto const& rlut = noiseRegionLUT(N, K);
        Operon::Scalar acc = 0;
        for (std::size_t j = 0; j < M; ++j) {
            auto const h = pointToPolylineSegOne(static_cast<double>(candX[j]), static_cast<double>(candY[j]),
                                                 drawn.data(), N);
            Operon::Scalar const d2 = static_cast<Operon::Scalar>(h.d2);
            if (K == 0) { acc += d2; }
            else {
                auto const k = rlut[h.nearestVertex];
                acc += pointLoss(d2, noise_.s2[k], noise_.nu[k]);
            }
        }
        return acc;
    }

    auto cisShapeError(std::vector<Operon::Scalar> const& gx,
                       std::vector<Operon::Scalar> const& gy,
                       bool directFin,
                       CisPointSegBatch const* psBatch = nullptr,
                       std::size_t psBatchIdx = 0,
                       CisPointSegBatch const* psBatch48 = nullptr) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        if (oneResOn()) { oneResCache_.valid = false; }
        auto const M = gx.size();
        auto const N = dataX_.size();
        if (M < 2 || N < 2) { return BIG; }

        thread_local std::vector<Operon::Scalar> cisCx, cisCy;
        cisCx.resize(M);
        cisCy.resize(M);
        {
            Operon::Scalar gx0 = gx[0], gy0 = gy[0];
            if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return BIG; }
            Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
            cisCx[0] = gx0; cisCy[0] = gy0;
            for (std::size_t j = 1; j < M; ++j) {
                Operon::Scalar const vx = gx[j], vy = gy[j];
                if (!std::isfinite(vx) || !std::isfinite(vy)) { return BIG; }
                cisCx[j] = vx; cisCy[j] = vy;
                if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
                if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                cisCx[j] = (cisCx[j] - cxc) / sp;
                cisCy[j] = (cisCy[j] - cyc) / sp;
            }
        }

        bool const profT12 = config_.profileTiming;

        thread_local std::vector<Vec2> strokePts, candPts;
        strokePts.resize(N);
        candPts.resize(M);
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < M; ++j) { candPts[j] = Vec2{cisCx[j], cisCy[j]}; }
        if (config_.useCoarseFourPointSieve
            && std::isfinite(config_.coarseFourPointThreshold)
            && coarseFourPointForwardLowerBound(strokePts, candPts) > config_.coarseFourPointThreshold) {
            return BIG;
        }
        std::chrono::steady_clock::time_point tPs0;
        if (profT12) { tPs0 = std::chrono::steady_clock::now(); }
        PointSegmentResult psRes;
        bool const psBatchValidHere = psBatch != nullptr
            && (psBatch->valid.empty() || (psBatchIdx < psBatch->valid.size() && psBatch->valid[psBatchIdx]));
        mgTauCacheValid_ = false;
        Operon::Scalar psViewFit{0};
        Operon::Scalar const* psViewRes = nullptr;
        Operon::Scalar const* psViewTau = nullptr;
        bool gpuFieldUsed = false;
        if (oneResOn()) {
            if (psBatchValidHere && psBatch->N == N
                && psBatch->viewCandidate(psBatchIdx, N, psViewFit, psViewRes, psViewTau)) {
                psRes.ok = true;
                psRes.fitness = psViewFit;
                psRes.residual.assign(psViewRes, psViewRes + N);
                oneResCacheStore(gx, gy, N, psRes.residual.data(), psViewTau);
                gpuFieldUsed = true;
            } else {
                psRes = evaluatePointSegmentOne(strokePts, candPts, true, false, oneResUnsigned());
                if (psRes.ok && psRes.residual.size() == N) {
                    oneResCacheStore(gx, gy, N, psRes.residual.data(),
                                     psRes.tau.size() == N ? psRes.tau.data() : nullptr);
                }
            }
        } else if (psBatchValidHere && psBatch->N == N
            && psBatch->viewCandidate(psBatchIdx, N, psViewFit, psViewRes, psViewTau)) {
            psRes.ok = true;
            psRes.fitness = psViewFit;
            psRes.residual.assign(psViewRes, psViewRes + N);
            gpuFieldUsed = true;
            if (config_.useMonotonicityOrderGuard && psViewTau != nullptr) {
                mgTauCache_.assign(psViewTau, psViewTau + N);
                mgTauCacheN_ = N;
                mgTauCacheM_ = M;
                mgTauCacheGx0_ = gx[0]; mgTauCacheGy0_ = gy[0];
                mgTauCacheGxL_ = gx[M - 1]; mgTauCacheGyL_ = gy[M - 1];
                mgTauCacheValid_ = true;
            }
        } else if (config_.useMonotonicityOrderGuard) {
            PointSegmentConfig psCfg;
            psCfg.computeTau = true;
            psRes = evaluatePointSegment(strokePts, candPts, psCfg);
            if (psRes.ok && psRes.tau.size() == N) {
                mgTauCache_ = psRes.tau;
                mgTauCacheN_ = N;
                mgTauCacheM_ = M;
                mgTauCacheGx0_ = gx[0]; mgTauCacheGy0_ = gy[0];
                mgTauCacheGxL_ = gx[M - 1]; mgTauCacheGyL_ = gy[M - 1];
                mgTauCacheValid_ = true;
            }
        } else {
            psRes = evaluatePointSegment(strokePts, candPts);
        }
        if (profT12) {
            auto const tPs1 = std::chrono::steady_clock::now();
            profCisPointSegNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tPs1 - tPs0).count(), std::memory_order_relaxed);
        }
        if (!psRes.ok) { return BIG; }

        std::vector<double> residualD(psRes.residual.size());
        for (std::size_t i = 0; i < psRes.residual.size(); ++i) {
            residualD[i] = static_cast<double>(psRes.residual[i]);
        }
        std::vector<double> sigma2D;
        if (!noise_.empty()) {
            sigma2D.assign(noise_.s2.begin(), noise_.s2.end());
        } else {
            sigma2D.assign(1, static_cast<double>(noiseMeanS2_));
        }
        CisScoringConfig cisCfg{};
        cisCfg.useScaleInvariantFitness = config_.useSiCisFitness;
        cisCfg.residualPrecision = gpuFieldUsed ? psBatch->fp16ModeUsed : 0;
        std::chrono::steady_clock::time_point tNll0;
        if (profT12) { tNll0 = std::chrono::steady_clock::now(); }

        CisScoringResult cisRes;
        bool grid48Used = false;
        if (cisGrid48Bit0() && cisGrid48Ok_) {
            PointSegmentResult ps48;
            double invSc48 = 1.0;
            bool ps48FromGpu = false;
            if (config_.useB3Grid48Gpu && psBatch48 != nullptr) {
                Operon::Scalar fit48{0};
                Operon::Scalar const* res48 = nullptr;
                Operon::Scalar const* tau48 = nullptr;
                if (psBatch48->viewCandidate(psBatchIdx, cisGrid48Stroke_.size(), fit48, res48, tau48)) {
                    ps48.ok = true;
                    ps48.residual.assign(res48, res48 + cisGrid48Stroke_.size());
                    invSc48 = cisGrid48InvSc();
                    ps48FromGpu = true;
                }
            }
            if (!ps48FromGpu) {
                static bool const b3g48DiagOn = std::getenv("SKETCH_B3_GRID48_DIAG") != nullptr;
                struct Grid48DiagTotals {
                    std::atomic<int64_t> ns{0};
                    std::atomic<int64_t> calls{0};
                    ~Grid48DiagTotals()
                    {
                        if (b3g48DiagOn) {
                            std::fprintf(stderr, "[B3G48DIAG] totalCpuGrid48Ns=%lld calls=%lld meanNs=%.1f\n",
                                         static_cast<long long>(ns.load()), static_cast<long long>(calls.load()),
                                         calls.load() > 0 ? static_cast<double>(ns.load()) / static_cast<double>(calls.load()) : 0.0);
                        }
                    }
                };
                static Grid48DiagTotals g48DiagTotals;
                std::chrono::steady_clock::time_point tG0d;
                if (b3g48DiagOn) { tG0d = std::chrono::steady_clock::now(); }
                std::chrono::steady_clock::time_point tG0;
                if (profT12) { tG0 = std::chrono::steady_clock::now(); }
                auto g48 = cisGrid48PointSeg(candPts, false, false);
                ps48 = std::move(g48.ps);
                ps48.ok = g48.ok;
                invSc48 = g48.invSc;
                if (profT12) {
                    auto const tG1 = std::chrono::steady_clock::now();
                    profCisGrid48Ns_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tG1 - tG0).count(), std::memory_order_relaxed);
                }
                if (b3g48DiagOn) {
                    auto const tG1d = std::chrono::steady_clock::now();
                    g48DiagTotals.ns.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tG1d - tG0d).count(), std::memory_order_relaxed);
                    g48DiagTotals.calls.fetch_add(1, std::memory_order_relaxed);
                }
            }
            if (ps48.ok && ps48.residual.size() == cisGrid48Stroke_.size()) {
                std::vector<double> residual48D(ps48.residual.size());
                for (std::size_t i = 0; i < ps48.residual.size(); ++i) {
                    residual48D[i] = static_cast<double>(ps48.residual[i]) * invSc48;
                }
                CisPerCellEnergyInput direct48;
                if (cisGrid48Bit1() && cisGrid48Zi5Ok_) {
                    direct48.perCellEnergy.resize(cisGrid48Zi5LnE_.size());
                    for (std::size_t c = 0; c < cisGrid48Zi5LnE_.size(); ++c) {
                        direct48.perCellEnergy[c] = std::exp(cisGrid48Zi5LnE_[c]);
                    }
                    direct48.memberEnergy.resize(cisGrid48Zi5MemberLnE_.size());
                    for (std::size_t m = 0; m < cisGrid48Zi5MemberLnE_.size(); ++m) {
                        auto const& memberLnE = cisGrid48Zi5MemberLnE_[m];
                        direct48.memberEnergy[m].resize(memberLnE.size());
                        for (std::size_t c = 0; c < memberLnE.size(); ++c) {
                            direct48.memberEnergy[m][c] = std::exp(memberLnE[c]);
                        }
                    }
                    direct48.present = true;
                }
                CisScoringConfig cisCfg48 = cisCfg;
                cisCfg48.residualPrecision = 0;
                bool mode4Applied = false;
                if (cisCfg48.exactLikelihood == 4 && cisGrid48Zi5Ok_ && direct48.present
                    && !direct48.memberEnergy.empty()) {
                    auto const dec48 = haarForward(residual48D);
                    auto cells48 = buildCisCells(dec48, std::vector<double>(1, 1.0), {}, true);
                    bool sizeOk = cells48.size() == cisGrid48Zi5LnE_.size();
                    for (auto const& me : direct48.memberEnergy) {
                        if (me.size() != cells48.size()) { sizeOk = false; break; }
                    }
                    if (sizeOk) {
                        std::vector<double> memberNll, memberW;
                        bool dataOk = true;
                        double const absVal = cisAbsDataNllEnsemble(cells48, direct48.memberEnergy, memberNll, memberW, dataOk);
                        if (dataOk && std::isfinite(absVal)) {
                            cisRes = CisScoringResult{};
                            cisRes.cells = std::move(cells48);
                            cisRes.shapeNll = absVal;
                            cisRes.scaleNll = 0.0;
                            cisRes.sHat = 1.0;
                            cisRes.fitness = absVal;
                            cisRes.ok = true;
                            mode4Applied = true;
                        }
                    }
                }
                if (!mode4Applied) {
                    cisRes = evaluateCisScoring(residual48D, sigma2D, {}, cisCfg48, direct48);
                }
                grid48Used = true;
            } else {
                cisGrid48Counters().residualFail.fetch_add(1, std::memory_order_relaxed);
            }
        }
        if (!grid48Used) {
            cisRes = evaluateCisScoring(residualD, sigma2D, {}, cisCfg);
        }
        if (profT12) {
            auto const tNll1 = std::chrono::steady_clock::now();
            profCisNllNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tNll1 - tNll0).count(), std::memory_order_relaxed);
        }
        if (!cisRes.ok || !std::isfinite(cisRes.fitness)) { return BIG; }

        double shapeFitness = static_cast<double>(cisRes.fitness);
        if (config_.useGlobalRobustNllFitness) {
            thread_local std::vector<Operon::Scalar> grSigma2;
            if (config_.useExternalPointSigma2 && !externalPointSigma2_.empty()) {
                auto const Ke = externalPointSigma2_.size();
                auto const& rlutE = noiseRegionLUT(N, Ke);
                grSigma2.resize(N);
                for (std::size_t i = 0; i < N; ++i) { grSigma2[i] = externalPointSigma2_[rlutE[originalArcIndex(i)]]; }
            } else {
                auto const Kg = noise_.size();
                auto const& rlutG = noiseRegionLUT(N, Kg);
                grSigma2.assign(1, noiseMeanS2_);
                if (Kg > 0) {
                    grSigma2.resize(N);
                    for (std::size_t i = 0; i < N; ++i) { grSigma2[i] = noise_.s2[rlutG[originalArcIndex(i)]]; }
                }
            }
            GlobalRobustNllConfig grCfg;
            grCfg.nu = config_.globalRobustNllNu;
            auto const grRes = evaluateGlobalRobustNll(psRes.residual, grSigma2, grCfg);
            if (grRes.ok) { shapeFitness = grRes.fitness; }
        }

        std::chrono::steady_clock::time_point tRev0;
        if (profT12) { tRev0 = std::chrono::steady_clock::now(); }
        Operon::Scalar accCurve = 0;
        auto const K = noise_.size();
        auto const& rlut = noiseRegionLUT(N, K);
        if (oneResBit(2)) {
            if (psBatchValidHere && psBatch->reverseFitness.size() > psBatchIdx) {
                accCurve = psBatch->reverseFitness[psBatchIdx];
            } else {
                accCurve = oneResReverseHalfSum(cisCx, cisCy, M);
            }
        } else if (config_.useGridSdf && !directFin && gridSdf_.ok) {
            for (std::size_t j = 0; j < M; ++j) {
                uint8_t reg = 0;
                Operon::Scalar const d2 = queryGridSdf(gridSdf_, cisCx[j], cisCy[j], reg);
                accCurve += (K == 0) ? d2 : pointLoss(d2, noise_.s2[reg], noise_.nu[reg]);
            }
        } else {
            int hint = -1;
            for (std::size_t j = 0; j < M; ++j) {
                auto const nn = nearestPointD(static_cast<double>(cisCx[j]), static_cast<double>(cisCy[j]),
                                              dataX_, dataY_, dataBlocks_, &hint);
                Operon::Scalar const d2 = static_cast<Operon::Scalar>(nn.d2);
                if (K == 0) { accCurve += d2; }
                else { auto const k = rlut[originalArcIndex(nn.i)]; accCurve += pointLoss(d2, noise_.s2[k], noise_.nu[k]); }
            }
        }
        if (profT12) {
            auto const tRev1 = std::chrono::steady_clock::now();
            profCisReverseNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tRev1 - tRev0).count(), std::memory_order_relaxed);
        }
        accCurve /= static_cast<Operon::Scalar>(M);
        if (!std::isfinite(accCurve)) { return BIG; }

        if (auto* k1cOut = k1cMeasureFile()) {
            std::lock_guard<std::mutex> lk(k1cMeasureMutex());
            (*k1cOut) << static_cast<double>(shapeFitness) << ' ' << static_cast<double>(accCurve) << '\n';
        }

        double robustPenalty = 0.0;
        if (config_.robustContinuousWeight > 0.0) {
            RobustContinuousConfig robustCfg;
            robustCfg.tangentWeight = config_.robustContinuousTangentWeight;
            robustCfg.curvatureWeight = config_.robustContinuousCurvatureWeight;
            auto const robust = evaluateRobustContinuous(
                strokePts, candPts, std::sqrt(std::max(static_cast<double>(noiseMeanS2_), 1e-12)), robustCfg);
            if (!robust.ok) { return BIG; }
            robustPenalty = config_.robustContinuousWeight * robust.fitness;
        }
        Operon::Scalar const result = static_cast<Operon::Scalar>(shapeFitness + robustPenalty) + accCurve;
        if (!std::isfinite(result)) { return BIG; }
        return result;
    }

    void buildQuotientFiniteMinBundles()
    {
        quotientFiniteMinBundles_.clear();
        auto const n = dataX_.size();
        if (n < 2) { return; }

        double mx = 0.0, my = 0.0;
        for (std::size_t i = 0; i < n; ++i) { mx += dataX_[i]; my += dataY_[i]; }
        mx /= static_cast<double>(n);
        my /= static_cast<double>(n);
        double cxx = 0.0, cyy = 0.0, cxy = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double const dx = static_cast<double>(dataX_[i]) - mx;
            double const dy = static_cast<double>(dataY_[i]) - my;
            cxx += dx * dx; cyy += dy * dy; cxy += dx * dy;
        }
        cxx /= static_cast<double>(n); cyy /= static_cast<double>(n); cxy /= static_cast<double>(n);
        double const tr = cxx + cyy;
        double const diff = cxx - cyy;
        double const disc = std::sqrt(diff * diff / 4.0 + cxy * cxy);
        double const lambda1 = tr / 2.0 + disc, lambda2 = tr / 2.0 - disc;
        double const denom = lambda1 + lambda2;
        double const qfmGapTol = sf2(10) ? 1.0 / static_cast<double>(n) : kQfmEigenGapRelTol;
        double const qfmThetaTol = sf2(10) ? 2.0 * M_PI / static_cast<double>(n) : kQfmThetaTol;
        bool const axisNearDegenerate = !(denom > 1e-12) || (lambda1 - lambda2) < qfmGapTol * denom;

        bool reversalNearTie = false;
        if (!dataClosed_ && n >= 3) {
            std::vector<float> pxF(n), pyF(n);
            for (std::size_t i = 0; i < n; ++i) {
                pxF[i] = static_cast<float>(dataX_[i]);
                pyF[i] = static_cast<float>(dataY_[i]);
            }
            double const theta = Sketch::detail::totalTurningAngle(pxF.data(), pyF.data(), static_cast<int>(n));
            reversalNearTie = std::abs(theta) < qfmThetaTol;
        }

        std::vector<int> tiedShifts;
        if (dataClosed_ && n >= 2) {
            std::vector<double> r2(n);
            double sumR2 = 0.0, maxR2 = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                double const dx = static_cast<double>(dataX_[i]) - mx;
                double const dy = static_cast<double>(dataY_[i]) - my;
                r2[i] = dx * dx + dy * dy;
                sumR2 += r2[i];
                if (r2[i] > maxR2) { maxR2 = r2[i]; }
            }
            double const rmsR = std::sqrt(sumR2 / static_cast<double>(n));
            double const scaleRef = rmsR > 0.0 ? rmsR : 1.0;
            double const epsR2 = kQfmRadiusRelTol * (scaleRef * scaleRef);
            for (std::size_t i = 1; i < n; ++i) {
                if (r2[i] >= maxR2 - epsR2) {
                    tiedShifts.push_back(static_cast<int>(i));
                    if (static_cast<int>(tiedShifts.size()) >= kQfmMaxCyclicAlts) { break; }
                }
            }
        }

        if (!axisNearDegenerate && !reversalNearTie && tiedShifts.empty()) {
            return;
        }

        struct Variant { bool signFlip; bool reversalFlip; int cyclicShift; };
        std::vector<Variant> variants;
        std::vector<bool> const signOpts = axisNearDegenerate ? std::vector<bool>{false, true} : std::vector<bool>{false};
        std::vector<bool> const revOpts = reversalNearTie ? std::vector<bool>{false, true} : std::vector<bool>{false};
        std::vector<int> shiftOpts{0};
        for (int s : tiedShifts) { shiftOpts.push_back(s); }
        for (bool sf : signOpts) {
            for (bool rf : revOpts) {
                for (int sh : shiftOpts) {
                    if (!sf && !rf && sh == 0) { continue; }
                    variants.push_back({sf, rf, sh});
                    if (variants.size() >= kQfmMaxBundles) { break; }
                }
                if (variants.size() >= kQfmMaxBundles) { break; }
            }
            if (variants.size() >= kQfmMaxBundles) { break; }
        }

        quotientFiniteMinBundles_.reserve(variants.size());
        for (auto const& v : variants) {
            QuotientFiniteMinBundle b;
            b.x = dataX_;
            b.y = dataY_;
            if (v.cyclicShift != 0 && static_cast<std::size_t>(v.cyclicShift) < n) {
                std::rotate(b.x.begin(), b.x.begin() + v.cyclicShift, b.x.end());
                std::rotate(b.y.begin(), b.y.begin() + v.cyclicShift, b.y.end());
            }
            if (v.reversalFlip) {
                std::reverse(b.x.begin(), b.x.end());
                std::reverse(b.y.begin(), b.y.end());
            }
            if (v.signFlip) {
                for (std::size_t i = 0; i < n; ++i) {
                    b.x[i] = static_cast<Operon::Scalar>(2.0 * mx - static_cast<double>(b.x[i]));
                    b.y[i] = static_cast<Operon::Scalar>(2.0 * my - static_cast<double>(b.y[i]));
                }
            }
            quotientFiniteMinBundles_.push_back(std::move(b));
        }
    }

    auto siCisScoreAgainstTarget(std::vector<Operon::Scalar> const& candX, std::vector<Operon::Scalar> const& candY,
                                 std::vector<Operon::Scalar> const& targetX, std::vector<Operon::Scalar> const& targetY) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto const M = candX.size();
        auto const N = targetX.size();
        if (M < 2 || N < 2) { return BIG; }
        std::vector<Vec2> strokePts(N), candPts(M);
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{targetX[i], targetY[i]}; }
        for (std::size_t j = 0; j < M; ++j) { candPts[j] = Vec2{candX[j], candY[j]}; }
        auto const psRes = pointSegOp(strokePts, candPts, false, false);
        if (!psRes.ok) { return BIG; }

        std::vector<double> residualD(psRes.residual.size());
        for (std::size_t i = 0; i < psRes.residual.size(); ++i) {
            residualD[i] = static_cast<double>(psRes.residual[i]);
        }
        std::vector<double> sigma2D;
        if (!noise_.empty()) {
            sigma2D.assign(noise_.s2.begin(), noise_.s2.end());
        } else {
            sigma2D.assign(1, static_cast<double>(noiseMeanS2_));
        }
        CisScoringConfig cfg{};
        cfg.useScaleInvariantFitness = true;
        auto const cisRes = evaluateCisScoring(residualD, sigma2D, {}, cfg);
        if (!cisRes.ok || !std::isfinite(cisRes.fitness)) { return BIG; }
        return static_cast<Operon::Scalar>(cisRes.fitness);
    }

    auto siCisRecommendationScore(PairIndividual const& ind, Operon::Dataset const& ds,
                                  Operon::Range range) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto& arena = getEvaluatorArena();
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        interpInto(ind.x, cx, ds, range, arena.gx);
        interpInto(ind.y, cy, ds, range, arena.gy);
        auto const& gx = arena.gx;
        auto const& gy = arena.gy;
        auto const M = gx.size();
        auto const N = dataX_.size();
        if (M < 2 || N < 2) { return BIG; }

        thread_local std::vector<Operon::Scalar> recCx, recCy;
        recCx.resize(M);
        recCy.resize(M);
        {
            Operon::Scalar gx0 = gx[0], gy0 = gy[0];
            if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return BIG; }
            Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
            recCx[0] = gx0; recCy[0] = gy0;
            for (std::size_t j = 1; j < M; ++j) {
                Operon::Scalar const vx = gx[j], vy = gy[j];
                if (!std::isfinite(vx) || !std::isfinite(vy)) { return BIG; }
                recCx[j] = vx; recCy[j] = vy;
                if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
                if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                recCx[j] = (recCx[j] - cxc) / sp;
                recCy[j] = (recCy[j] - cyc) / sp;
            }
        }

        Operon::Scalar best = siCisScoreAgainstTarget(recCx, recCy, dataX_, dataY_);
        if (config_.useQuotientFiniteMin) {
            for (auto const& b : quotientFiniteMinBundles_) {
                auto const s = siCisScoreAgainstTarget(recCx, recCy, b.x, b.y);
                if (s < best) { best = s; }
            }
        }
        return best;
    }

    auto buildLineNull(SearchState const& s) const -> PairIndividual
    {
        using BO = Operon::BuiltinOp;
        auto varNode = [&]() {
            Operon::Node n(Operon::NodeType::Variable);
            n.HashValue = s.inputs[0];
            n.Value = Operon::Scalar{1};
            return n;
        };
        auto buildAxis = [&](Operon::Scalar c0, Operon::Scalar c1) {
            std::vector<Operon::Node> nodes;
            nodes.reserve(5);
            nodes.push_back(Operon::Node::Constant(c1));
            nodes.push_back(varNode());
            nodes.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Mul), 2));
            nodes.push_back(Operon::Node::Constant(c0));
            nodes.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Add), 2));
            Operon::Tree tree(nodes);
            tree.UpdateNodes();
            return tree;
        };
        Operon::Scalar const x0 = dataX_.empty() ? Operon::Scalar{0} : dataX_.front();
        Operon::Scalar const x1 = dataX_.empty() ? Operon::Scalar{0} : (dataX_.back() - dataX_.front());
        Operon::Scalar const y0 = dataY_.empty() ? Operon::Scalar{0} : dataY_.front();
        Operon::Scalar const y1 = dataY_.empty() ? Operon::Scalar{0} : (dataY_.back() - dataY_.front());
        PairIndividual line;
        line.x = buildAxis(x0, x1);
        line.y = buildAxis(y0, y1);
        return line;
    }

    auto siCisRecommend(std::vector<PairIndividual*> const& front, Operon::Dataset const& ds,
                        Operon::Range range) const -> PairIndividual*
    {
        PairIndividual* best = front.back();
        Operon::Scalar bestScore = std::numeric_limits<Operon::Scalar>::max();
        for (auto* p : front) {
            auto const score = siCisRecommendationScore(*p, ds, range);
            if (std::isfinite(score) && score < bestScore) {
                bestScore = score;
                best = p;
            }
        }
        return best;
    }

    struct NmlRanges { double extent{1.0}; double freq{0.0}; double phase{0.0}; double expRate{0.0}; bool expRateOk{false}; double sigmaBar{0.0}; double gridPeriod{0.0}; };
    auto nmlRanges() const -> NmlRanges
    {
        NmlRanges r;
        if (!dataX_.empty()) {
            auto [xlo, xhi] = std::minmax_element(dataX_.begin(), dataX_.end());
            auto [ylo, yhi] = std::minmax_element(dataY_.begin(), dataY_.end());
            r.extent = std::max(static_cast<double>(*xhi - *xlo), static_cast<double>(*yhi - *ylo));
        }
        if (!(r.extent > 0.0)) { r.extent = 1.0; }
        double const N = static_cast<double>(dataX_.size());
        r.freq = Sketch::nml::kPi * std::max(N - 1.0, 1.0);
        r.phase = 2.0 * Sketch::nml::kPi;
        double const absS2 = config_.nmlAbsSigma2;
        bool const haveAbs = std::isfinite(absS2) && absS2 > 0.0;
        double const sigmaBar = haveAbs ? std::sqrt(absS2) : std::sqrt(std::max(static_cast<double>(noiseMeanS2_), 0.0));
        r.sigmaBar = sigmaBar;
        double const lam = (sigmaBar > 0.0) ? std::log(r.extent / sigmaBar) : 0.0;
        r.expRateOk = (haveAbs || !noise_.empty()) && std::isfinite(lam) && lam > 0.0;
        r.expRate = r.expRateOk ? 2.0 * lam : r.extent;
        return r;
    }

    static auto nmlConstRoles(Operon::Tree const& tree) -> std::vector<Sketch::nml::ConstRole>
    {
        using Sketch::nml::ConstRole;
        auto const& nodes = tree.Nodes();
        std::vector<ConstRole> roles;
        if (nodes.empty()) { return roles; }
        std::vector<char> hasVar(nodes.size(), 0);
        std::vector<long> parent(nodes.size(), -1);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            char h = nodes[i].IsVariable() ? char{1} : char{0};
            if (nodes[i].Arity > 0) {
                for (auto c : tree.Indices(i)) {
                    h = (h != 0 || hasVar[c] != 0) ? char{1} : char{0};
                    parent[c] = static_cast<long>(i);
                }
            }
            hasVar[i] = h;
        }
        auto isChainOp = [](Operon::Node const& n) { return n.IsAddition() || n.IsSubtraction() || n.IsMultiplication(); };
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (!nodes[i].Optimize) { continue; }
            bool const isVar = nodes[i].IsVariable();
            bool mulVar = false, anyMul = false;
            std::size_t cur = i;
            long p = parent[i];
            while (p >= 0 && isChainOp(nodes[static_cast<std::size_t>(p)])) {
                auto const pu = static_cast<std::size_t>(p);
                if (nodes[pu].IsMultiplication()) {
                    anyMul = true;
                    for (auto c : tree.Indices(pu)) {
                        if (c != cur && hasVar[c] != 0) { mulVar = true; }
                    }
                }
                cur = pu;
                p = parent[pu];
            }
            ConstRole role = ConstRole::Extent;
            if (p >= 0) {
                auto const& a = nodes[static_cast<std::size_t>(p)];
                if (a.IsSin() || a.IsCos()) {
                    if (isVar || mulVar) { role = ConstRole::Frequency; }
                    else if (!anyMul) { role = ConstRole::Phase; }
                    else { role = ConstRole::Unclassified; }
                } else if (a.IsExp()) {
                    role = (isVar || mulVar) ? ConstRole::ExpRate : ConstRole::Unclassified;
                } else {
                    long q = parent[static_cast<std::size_t>(p)];
                    while (q >= 0) {
                        auto const& b = nodes[static_cast<std::size_t>(q)];
                        if (b.IsSin() || b.IsCos() || b.IsExp()) { role = ConstRole::Unclassified; break; }
                        q = parent[static_cast<std::size_t>(q)];
                    }
                }
            }
            roles.push_back(role);
        }
        return roles;
    }

    auto nmlTermsFor(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range,
                     NmlRanges const& rg) const -> Sketch::nml::Terms
    {
        using Sketch::nml::ConstRole;
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        Sketch::nml::Terms t;
        t.nodesX = static_cast<int>(ind.x.Length());
        t.nodesY = static_cast<int>(ind.y.Length());
        auto const cx = ind.x.GetCoefficients();
        auto const cy = ind.y.GetCoefficients();
        std::size_t const nx = cx.size(), ny = cy.size(), k = nx + ny;
        t.k = static_cast<int>(k);

        auto rolesX = nmlConstRoles(ind.x);
        auto rolesY = nmlConstRoles(ind.y);
        if (rolesX.size() != nx || rolesY.size() != ny) { return t; }
        std::vector<ConstRole> roles = rolesX;
        roles.insert(roles.end(), rolesY.begin(), rolesY.end());
        std::vector<double> theta(k);
        for (std::size_t j = 0; j < nx; ++j) { theta[j] = static_cast<double>(cx[j]); }
        for (std::size_t j = 0; j < ny; ++j) { theta[nx + j] = static_cast<double>(cy[j]); }
        std::vector<double> rangeOf(k, rg.extent);
        for (std::size_t j = 0; j < k; ++j) {
            double const a = std::abs(theta[j]);
            t.maxAbsConst = std::max(t.maxAbsConst, a);
            switch (roles[j]) {
            case ConstRole::Frequency:
                ++t.nFreq; rangeOf[j] = rg.freq; t.maxAbsFreq = std::max(t.maxAbsFreq, a);
                if (a > rg.freq) { ++t.nOutOfRange; }
                break;
            case ConstRole::Phase: ++t.nPhase; rangeOf[j] = rg.phase; break;
            case ConstRole::ExpRate:
                if (rg.expRateOk) {
                    ++t.nExpRate; rangeOf[j] = rg.expRate;
                    if (a > 0.5 * rg.expRate) { ++t.nOutOfRange; }
                } else { ++t.nUnclassified; rangeOf[j] = rg.extent; }
                break;
            case ConstRole::Unclassified: ++t.nUnclassified; rangeOf[j] = rg.extent; break;
            case ConstRole::Extent: default: ++t.nExtent; rangeOf[j] = rg.extent; break;
            }
        }

        auto const nllS = siCisRecommendationScore(ind, ds, range);
        if (!std::isfinite(nllS) || nllS >= BIG) { return t; }
        t.nll = static_cast<double>(nllS);

        std::vector<Operon::Scalar> U, V;
        interpInto(ind.x, cx, ds, range, U);
        interpInto(ind.y, cy, ds, range, V);
        auto const M = U.size();
        auto const N = dataX_.size();
        if (M < 2 || N < 2 || V.size() != M) { return t; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(U[j]) || !std::isfinite(V[j])) { return t; }
        }
        std::vector<double> Ju, Jv;
        if (!treeJacobian(ind.x, cx, ds, range, M, Ju)) { return t; }
        if (!treeJacobian(ind.y, cy, ds, range, M, Jv)) { return t; }

        std::size_t iUmin = 0, iUmax = 0, iVmin = 0, iVmax = 0;
        for (std::size_t j = 1; j < M; ++j) {
            if (U[j] < U[iUmin]) { iUmin = j; }
            if (U[j] > U[iUmax]) { iUmax = j; }
            if (V[j] < V[iVmin]) { iVmin = j; }
            if (V[j] > V[iVmax]) { iVmax = j; }
        }
        double const Cx = 0.5 * (static_cast<double>(U[iUmax]) + static_cast<double>(U[iUmin]));
        double const Cy = 0.5 * (static_cast<double>(V[iVmax]) + static_cast<double>(V[iVmin]));
        double const spanU = static_cast<double>(U[iUmax]) - static_cast<double>(U[iUmin]);
        double const spanV = static_cast<double>(V[iVmax]) - static_cast<double>(V[iVmin]);
        double const S = std::max({spanU, spanV, 1e-9});
        if (!(S > 0.0) || !std::isfinite(S)) { return t; }
        std::vector<double> dCx(k, 0.0), dCy(k, 0.0), dS(k, 0.0);
        for (std::size_t j = 0; j < nx; ++j) { dCx[j] = 0.5 * (Ju[iUmax * nx + j] + Ju[iUmin * nx + j]); }
        for (std::size_t j = 0; j < ny; ++j) { dCy[nx + j] = 0.5 * (Jv[iVmax * ny + j] + Jv[iVmin * ny + j]); }
        if (spanU >= spanV && spanU > 1e-9) {
            for (std::size_t j = 0; j < nx; ++j) { dS[j] = Ju[iUmax * nx + j] - Ju[iUmin * nx + j]; }
        } else if (spanV > 1e-9) {
            for (std::size_t j = 0; j < ny; ++j) { dS[nx + j] = Jv[iVmax * ny + j] - Jv[iVmin * ny + j]; }
        }
        double const invS = 1.0 / S;
        std::vector<double> dqx(M * k), dqy(M * k);
        for (std::size_t j = 0; j < M; ++j) {
            double const ax = (static_cast<double>(U[j]) - Cx) * invS;
            double const ay = (static_cast<double>(V[j]) - Cy) * invS;
            for (std::size_t c = 0; c < k; ++c) {
                double const gu = (c < nx) ? Ju[j * nx + c] : 0.0;
                double const gv = (c >= nx) ? Jv[j * ny + (c - nx)] : 0.0;
                dqx[j * k + c] = (gu - dCx[c]) * invS - ax * dS[c] * invS;
                dqy[j * k + c] = (gv - dCy[c]) * invS - ay * dS[c] * invS;
            }
        }

        std::vector<Vec2> strokePts(N), candPts(M);
        {
            Operon::Scalar xlo = U[0], xhi = U[0], ylo = V[0], yhi = V[0];
            for (std::size_t j = 1; j < M; ++j) {
                if (U[j] < xlo) { xlo = U[j]; } else if (U[j] > xhi) { xhi = U[j]; }
                if (V[j] < ylo) { ylo = V[j]; } else if (V[j] > yhi) { yhi = V[j]; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) { candPts[j] = Vec2{(U[j] - cxc) / sp, (V[j] - cyc) / sp}; }
        }
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        auto const ps = pointSegOp(strokePts, candPts, true, true);
        if (!ps.ok || ps.residual.size() != N || ps.tau.size() != N || ps.normalAtClosest.size() != N) { return t; }

        Eigen::MatrixXd Jr(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(k));
        for (std::size_t i = 0; i < N; ++i) {
            double tau = static_cast<double>(ps.tau[i]);
            if (!std::isfinite(tau)) { return t; }
            tau = std::clamp(tau, 0.0, static_cast<double>(M - 1));
            auto s0 = static_cast<std::size_t>(std::floor(tau));
            if (s0 > M - 2) { s0 = M - 2; }
            double const f = std::clamp(tau - static_cast<double>(s0), 0.0, 1.0);
            double const nxv = static_cast<double>(ps.normalAtClosest[i].x);
            double const nyv = static_cast<double>(ps.normalAtClosest[i].y);
            for (std::size_t c = 0; c < k; ++c) {
                double const gx = (1.0 - f) * dqx[s0 * k + c] + f * dqx[(s0 + 1) * k + c];
                double const gy = (1.0 - f) * dqy[s0 * k + c] + f * dqy[(s0 + 1) * k + c];
                Jr(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(c)) = -(nxv * gx + nyv * gy);
            }
        }

        std::vector<double> residualD(N);
        for (std::size_t i = 0; i < N; ++i) { residualD[i] = static_cast<double>(ps.residual[i]); }
        std::vector<double> sigma2D;
        if (!noise_.empty()) { sigma2D.assign(noise_.s2.begin(), noise_.s2.end()); }
        else { sigma2D.assign(1, static_cast<double>(noiseMeanS2_)); }
        CisScoringConfig ccfg{};
        ccfg.useScaleInvariantFitness = true;
        auto const cis = evaluateCisScoring(residualD, sigma2D, {}, ccfg);
        if (!cis.ok || cis.cells.empty()) { return t; }
        std::size_t const nCells = cis.cells.size();
        t.n = static_cast<int>(nCells);
        bool const exactLaw = (ccfg.exactLikelihood >= 1 && ccfg.exactLikelihood <= 3);
        double const eps = exactLaw ? 0.0 : 1e-300;
        double sumRatio = 0.0;
        std::vector<double> pred(nCells);
        for (std::size_t c = 0; c < nCells; ++c) {
            pred[c] = std::max(cis.cells[c].predictedEnergy, eps);
            sumRatio += (std::max(cis.cells[c].energy, 0.0) + eps) / (pred[c] + eps);
        }
        double const sHat = exactLaw ? cis.sHat : sumRatio / static_cast<double>(nCells);
        if (!(sHat > 0.0) || !std::isfinite(sHat)) { return t; }

        Eigen::MatrixXd Jw(static_cast<Eigen::Index>(nCells), static_cast<Eigen::Index>(k));
        std::vector<double> col(N);
        for (std::size_t c = 0; c < k; ++c) {
            for (std::size_t i = 0; i < N; ++i) { col[i] = Jr(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(c)); }
            auto const flat = haarFlattenCoeffs(haarForward(col));
            if (flat.size() != nCells) { return t; }
            for (std::size_t r = 0; r < nCells; ++r) {
                Jw(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = flat[r] / std::sqrt(sHat * (pred[r] + eps));
            }
        }
        auto const id = Sketch::nml::identifiableLogDet(Jw, static_cast<double>(std::numeric_limits<Operon::Scalar>::epsilon()));
        if (!id.ok) { return t; }
        t.kIdentifiable = id.rank;
        t.nDead = id.nDead;
        t.nCollinear = id.nCollinear;
        double sumLogRange = 0.0;
        for (int c : id.kept) { sumLogRange += std::log(rangeOf[static_cast<std::size_t>(c)]); }

        std::vector<char> isKept(k, 0);
        for (int c : id.kept) { isKept[static_cast<std::size_t>(c)] = 1; }
        std::vector<char> inGroup(k, 0);
        double const lamBox = 0.5 * rg.expRate;
        std::vector<double> thetaAlias = theta;
        bool anyFold = false;
        for (int tr = 0; tr < 2; ++tr) {
            auto const groups = Sketch::nml::linearGroups(tr == 0 ? ind.x : ind.y);
            std::size_t const off = (tr == 0) ? 0 : nx;
            for (auto const& g : groups) {
                bool anyKept = false;
                for (int c : g.coefs) {
                    auto const gc = off + static_cast<std::size_t>(c);
                    if (gc < k) { inGroup[gc] = 1; if (isKept[gc] != 0) { anyKept = true; } }
                }
                if (!anyKept) { continue; }
                double const a = std::abs(g.lin);
                if (g.trig) {
                    t.maxAbsLinFreq = std::max(t.maxAbsLinFreq, a);
                    if (!(a <= rg.freq)) {
                        ++t.nFreqOut;
                        if (g.foldCoef >= 0 && rg.gridPeriod > 0.0) {
                            double const folded = Sketch::nml::foldToGrid(g.lin, rg.gridPeriod);
                            if (std::abs(folded) <= rg.freq) { ++t.nAliasInRange; } else { ++t.nAliasOutOfRange; }
                            auto const gc = off + static_cast<std::size_t>(g.foldCoef);
                            thetaAlias[gc] += g.foldSign * (folded - g.lin);
                            anyFold = true;
                        } else {
                            ++t.nAliasNotFoldable;
                        }
                    }
                } else if (rg.expRateOk) {
                    t.maxAbsLinExp = std::max(t.maxAbsLinExp, a);
                    if (!(a <= lamBox)) { ++t.nExpOut; ++t.nAliasNotFoldable; }
                }
            }
        }
        for (std::size_t j = 0; j < k; ++j) {
            if (roles[j] == ConstRole::Extent && std::abs(theta[j]) > rg.extent) { ++t.nExtentOverEAll; }
            if (roles[j] == ConstRole::Phase && std::abs(theta[j]) > rg.phase) { ++t.nPhaseOutAll; }
            if (isKept[j] == 0) { continue; }
            double const a = std::abs(theta[j]);
            if (inGroup[j] == 0) {
                if (roles[j] == ConstRole::Frequency && !(a <= rg.freq)) { ++t.nIndivOut; }
                if (roles[j] == ConstRole::ExpRate && rg.expRateOk && !(a <= lamBox)) { ++t.nIndivOut; }
            }
            if (roles[j] == ConstRole::Extent && a > rg.extent) { ++t.nExtentOverE; }
        }
        t.nOutOfRangeId = t.nFreqOut + t.nExpOut + t.nIndivOut;
        t.admissible = (t.nOutOfRangeId == 0);
        if (anyFold) {
            PairIndividual al = ind;
            std::vector<Operon::Scalar> ax(nx), ay(ny);
            for (std::size_t j = 0; j < nx; ++j) { ax[j] = static_cast<Operon::Scalar>(thetaAlias[j]); }
            for (std::size_t j = 0; j < ny; ++j) { ay[j] = static_cast<Operon::Scalar>(thetaAlias[nx + j]); }
            al.x.SetCoefficients(ax);
            al.y.SetCoefficients(ay);
            std::vector<Operon::Scalar> U2, V2;
            interpInto(al.x, ax, ds, range, U2);
            interpInto(al.y, ay, ds, range, V2);
            if (U2.size() == M && V2.size() == M) {
                double md = 0.0;
                for (std::size_t j = 0; j < M; ++j) {
                    md = std::max({md, std::abs(static_cast<double>(U2[j]) - static_cast<double>(U[j])),
                                   std::abs(static_cast<double>(V2[j]) - static_cast<double>(V[j]))});
                }
                t.aliasGridMaxDiff = md / S;
            }
            auto const an = siCisRecommendationScore(al, ds, range);
            if (std::isfinite(an) && an < BIG) { t.aliasNll = static_cast<double>(an); }
        }

        Sketch::nml::assemble(t, id.logDet, sumLogRange, config_.nmlCanonicalStructCode,
                              config_.nmlOutOfRangeInadmissible);

        if (Sketch::nmlAbsModeFromEnv() != 0 || (Sketch::sketchT13Mask() & 0x3) != 0) {
            bool absOk = cisGrid48Ok_ && cisGrid48Zi5Ok_ && !cisGrid48Zi5MemberLnE_.empty()
                         && !cisGrid48Zi5LnE_.empty();
            if (absOk) {
                auto const g48 = cisGrid48PointSeg(candPts, true, true);
                std::size_t const NC48 = g48.strokeN;
                absOk = g48.ok && g48.ps.tau.size() == NC48 && g48.ps.normalAtClosest.size() == NC48;
                if (absOk) {
                    std::vector<double> residual48(NC48);
                    for (std::size_t i = 0; i < NC48; ++i) {
                        residual48[i] = static_cast<double>(g48.ps.residual[i]) * g48.invSc;
                    }
                    auto const dec48 = haarForward(residual48);
                    auto cells48 = buildCisCells(dec48, std::vector<double>(1, 1.0), {}, true);
                    std::size_t const nCells48 = cells48.size();
                    absOk = nCells48 == cisGrid48Zi5LnE_.size();
                    std::size_t const Mm = cisGrid48Zi5MemberLnE_.size();
                    for (auto const& me : cisGrid48Zi5MemberLnE_) {
                        if (me.size() != nCells48) { absOk = false; break; }
                    }
                    if (absOk && Mm > 0) {
                        std::vector<std::vector<double>> memberEnergy(Mm, std::vector<double>(nCells48));
                        for (std::size_t m = 0; m < Mm; ++m) {
                            for (std::size_t c = 0; c < nCells48; ++c) {
                                memberEnergy[m][c] = std::exp(cisGrid48Zi5MemberLnE_[m][c]);
                            }
                        }
                        std::vector<double> memberNll, memberW;
                        bool dataOk = true;
                        double const absNllVal = cisAbsDataNllEnsemble(cells48, memberEnergy, memberNll, memberW, dataOk);
                        absOk = dataOk && std::isfinite(absNllVal);
                        if (absOk) {
                            t.absNll = absNllVal;
                            double maxW = 0.0;
                            for (double w : memberW) { maxW = std::max(maxW, w); }
                            t.absMaxWeight = maxW;

                            Operon::Scalar exlo = cisGrid48Stroke_[0].x, exhi = exlo;
                            Operon::Scalar eylo = cisGrid48Stroke_[0].y, eyhi = eylo;
                            for (auto const& v : cisGrid48Stroke_) {
                                exlo = std::min(exlo, v.x); exhi = std::max(exhi, v.x);
                                eylo = std::min(eylo, v.y); eyhi = std::max(eyhi, v.y);
                            }
                            double const E48 = std::max(static_cast<double>(exhi - exlo), static_cast<double>(eyhi - eylo));
                            int const dcLevel = static_cast<int>(dec48.details.size()) + 1;
                            double sumPbar = 0.0; std::size_t nPbar = 0;
                            for (std::size_t c = 0; c < nCells48; ++c) {
                                if (cells48[c].level == dcLevel) { continue; }
                                sumPbar += std::exp(cisGrid48Zi5LnE_[c]);
                                ++nPbar;
                            }
                            double const sigmaBar48 = (nPbar > 0) ? std::sqrt(std::max(sumPbar / static_cast<double>(nPbar), 0.0)) : 0.0;
                            t.absSigmaBar = sigmaBar48;
                            double const lam48 = (sigmaBar48 > 0.0 && E48 > 0.0) ? std::log(E48 / sigmaBar48) : 0.0;
                            bool const absExpRateOk = std::isfinite(lam48) && lam48 > 0.0;
                            t.absLambda = lam48;

                            double sumLogRangeAbs = 0.0;
                            for (int c : id.kept) {
                                double rr = rangeOf[static_cast<std::size_t>(c)];
                                if (roles[static_cast<std::size_t>(c)] == ConstRole::ExpRate && absExpRateOk) {
                                    rr = 2.0 * lam48;
                                }
                                sumLogRangeAbs += std::log(rr);
                            }
                            t.absRangeTerm = sumLogRangeAbs;

                            std::size_t const kId = id.kept.size();
                            if (kId > 0) {
                                std::vector<std::vector<double>> flatKept(kId);
                                bool haarOk = true;
                                for (std::size_t ci = 0; ci < kId && haarOk; ++ci) {
                                    int const c = id.kept[ci];
                                    std::vector<double> colv(NC48);
                                    for (std::size_t i = 0; i < NC48; ++i) {
                                        double tau = static_cast<double>(g48.ps.tau[i]);
                                        if (!std::isfinite(tau)) { haarOk = false; break; }
                                        tau = std::clamp(tau, 0.0, static_cast<double>(M - 1));
                                        auto s0 = static_cast<std::size_t>(std::floor(tau));
                                        if (s0 > M - 2) { s0 = M - 2; }
                                        double const f = std::clamp(tau - static_cast<double>(s0), 0.0, 1.0);
                                        double const nxv = static_cast<double>(g48.ps.normalAtClosest[i].x);
                                        double const nyv = static_cast<double>(g48.ps.normalAtClosest[i].y);
                                        double const gx = (1.0 - f) * dqx[s0 * k + static_cast<std::size_t>(c)]
                                                         + f * dqx[(s0 + 1) * k + static_cast<std::size_t>(c)];
                                        double const gy = (1.0 - f) * dqy[s0 * k + static_cast<std::size_t>(c)]
                                                         + f * dqy[(s0 + 1) * k + static_cast<std::size_t>(c)];
                                        colv[i] = -(nxv * gx + nyv * gy) * g48.invSc;
                                    }
                                    if (!haarOk) { break; }
                                    auto const flat = haarFlattenCoeffs(haarForward(colv));
                                    if (flat.size() != nCells48) { haarOk = false; break; }
                                    flatKept[ci] = std::move(flat);
                                }
                                if (haarOk) {
                                    Eigen::MatrixXd Jstack(static_cast<Eigen::Index>(Mm * nCells48), static_cast<Eigen::Index>(kId));
                                    bool wOk = true;
                                    for (std::size_t m = 0; m < Mm && wOk; ++m) {
                                        double const sw = std::sqrt(std::max(memberW[m], 0.0));
                                        for (std::size_t r = 0; r < nCells48; ++r) {
                                            double const pcm = memberEnergy[m][r];
                                            if (!(pcm > 0.0) || !std::isfinite(pcm)) { wOk = false; break; }
                                            double const invSqrtP = 1.0 / std::sqrt(pcm);
                                            for (std::size_t ci = 0; ci < kId; ++ci) {
                                                Jstack(static_cast<Eigen::Index>(m * nCells48 + r), static_cast<Eigen::Index>(ci))
                                                    = sw * flatKept[ci][r] * invSqrtP;
                                            }
                                        }
                                    }
                                    if (wOk) {
                                        auto const idAbs = Sketch::nml::identifiableLogDet(
                                            Jstack, static_cast<double>(std::numeric_limits<Operon::Scalar>::epsilon()));
                                        if (idAbs.ok && idAbs.rank == static_cast<int>(kId)) {
                                            t.absFisherTerm = 0.5 * (idAbs.logDet - static_cast<double>(kId) * std::log(static_cast<double>(nCells48)));
                                            t.absTotal = t.absNll + t.paramTerm + t.absFisherTerm + t.absRangeTerm + t.structTerm;
                                            t.absOk = std::isfinite(t.absTotal);
                                        }
                                    }
                                }
                            } else {
                                t.absFisherTerm = 0.0;
                                t.absRangeTerm = 0.0;
                                t.absTotal = t.absNll + t.structTerm;
                                t.absOk = std::isfinite(t.absTotal);
                            }
                        }
                    } else {
                        absOk = false;
                    }
                }
            }
            if (!t.absOk) {
                cisGrid48Counters().absFallback.fetch_add(1, std::memory_order_relaxed);
            }
        }
        return t;
    }

    auto shapeOnlyPartsFromGrid(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                                Operon::Scalar* chamferOut,
                                CisPointSegBatch const* psBatch = nullptr,
                                std::size_t psBatchIdx = 0,
                                CisPointSegBatch const* psBatch48 = nullptr) const -> Operon::Scalar
    {
        bool const directFin = directFinalistEvalMode();
        if (!config_.profileTiming) {
            bool const sdNeeded = (chamferOut != nullptr) || (!cisFitnessMode() && !config_.useSymmetricFitness);
            Operon::Scalar sd = 0;
            if (sdNeeded) {
                sd = shapeDistance(gx, gy, dataX_, dataY_, noise_, nullptr, 4, &dataBlocks_, (config_.useGridSdf && !directFin) ? &gridSdf_ : nullptr, (config_.useAmxChamfer && !directFin) ? &amxData_ : nullptr);
                if (chamferOut != nullptr) { *chamferOut = sd; }
            }
            if (cisFitnessMode()) {
                sd = cisShapeError(gx, gy, directFin, psBatch, psBatchIdx, psBatch48);
            } else if (config_.useSymmetricFitness) {
                sd = symmetricFitnessError(gx, gy);
            }
            if (config_.orderWeight > Operon::Scalar{0}) { sd += orderTerm(gx, gy); }
            if (config_.curvatureWeight > Operon::Scalar{0}) { sd = applyCurvature(sd, curvatureTerm(gx, gy)); }
            return sd;
        }
        auto tDist0 = std::chrono::steady_clock::now();
        auto sd = shapeDistance(gx, gy, dataX_, dataY_, noise_, nullptr, 4, &dataBlocks_, (config_.useGridSdf && !directFin) ? &gridSdf_ : nullptr, (config_.useAmxChamfer && !directFin) ? &amxData_ : nullptr);
        auto tDist1 = std::chrono::steady_clock::now();
        profScoreDistNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tDist1 - tDist0).count(), std::memory_order_relaxed);
        if (chamferOut != nullptr) { *chamferOut = sd; }
        if (cisFitnessMode()) {
            auto tCis0 = std::chrono::steady_clock::now();
            sd = cisShapeError(gx, gy, directFin, psBatch, psBatchIdx, psBatch48);
            auto tCis1 = std::chrono::steady_clock::now();
            profScoreCisNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tCis1 - tCis0).count(), std::memory_order_relaxed);
        } else if (config_.useSymmetricFitness) {
            sd = symmetricFitnessError(gx, gy);
        }
        if (config_.orderWeight > Operon::Scalar{0}) {
            auto tOrd0 = std::chrono::steady_clock::now();
            sd += orderTerm(gx, gy);
            auto tOrd1 = std::chrono::steady_clock::now();
            profScoreOrderNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tOrd1 - tOrd0).count(), std::memory_order_relaxed);
        }
        if (config_.curvatureWeight > Operon::Scalar{0}) {
            auto tCurv0 = std::chrono::steady_clock::now();
            sd = applyCurvature(sd, curvatureTerm(gx, gy));
            auto tCurv1 = std::chrono::steady_clock::now();
            profScoreCurvNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tCurv1 - tCurv0).count(), std::memory_order_relaxed);
        }
        return sd;
    }

    auto shapeOnlyParts(PairIndividual const& ind, Operon::Dataset const& ds,
                        Operon::Range range, Operon::Scalar* chamferOut) const -> Operon::Scalar
    {
        auto& arena = getEvaluatorArena();
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        if (config_.profileTiming) {
            auto t0 = std::chrono::steady_clock::now();
            interpInto(ind.x, cx, ds, range, arena.gx);
            interpInto(ind.y, cy, ds, range, arena.gy);
            auto t1 = std::chrono::steady_clock::now();
            auto result = shapeOnlyPartsFromGrid(arena.gx, arena.gy, chamferOut);
            auto t2 = std::chrono::steady_clock::now();
            profEvalNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count(), std::memory_order_relaxed);
            profScoreNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count(), std::memory_order_relaxed);
            return result;
        }
        interpInto(ind.x, cx, ds, range, arena.gx);
        interpInto(ind.y, cy, ds, range, arena.gy);
        return shapeOnlyPartsFromGrid(arena.gx, arena.gy, chamferOut);
    }

    auto shapeOnlyGuarded(PairIndividual const& ind, Operon::Dataset const& ds,
                          Operon::Range range) const -> Operon::Scalar
    {
        try {
            return shapeOnly(ind, ds, range);
        } catch (std::exception const&) {
            ++nEvalErrors_;
        } catch (...) {
            ++nEvalErrors_;
        }
        return std::numeric_limits<Operon::Scalar>::max();
    }

    static auto polishWorsened(Operon::Scalar after, Operon::Scalar before) -> bool
    {
        if (!std::isfinite(before)) { return false; }
        return !(after <= before);
    }

    static auto structurallyValid(PairIndividual const& ind) -> bool
    {
        return Sketch::structurallyValid(ind.x) && Sketch::structurallyValid(ind.y);
    }

    static auto allConstantsFinite(Operon::Tree const& t) -> bool
    {
        for (auto const& n : t.Nodes()) {
            if (n.IsLeaf() && !std::isfinite(static_cast<double>(n.Value))) { return false; }
        }
        return true;
    }
    static auto allConstantsFiniteVec(std::vector<Operon::Scalar> const& v) -> bool
    {
        for (auto const x : v) { if (!std::isfinite(static_cast<double>(x))) { return false; } }
        return true;
    }

    void attributeNonFiniteOps(PairIndividual const& ind) const
    {
        static bool const measureOps = std::getenv("SKETCH_MEASURE_NONFINITE_OPS") != nullptr;
        if (!measureOps) { return; }
        bool hasLog = false, hasSqrt = false, hasDiv = false, hasExp = false;
        auto scan = [&](Operon::Tree const& t) {
            for (auto const& n : t.Nodes()) {
                if (n.Type != Operon::NodeType::Function) { continue; }
                if (n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Log)
                    || n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Logabs)
                    || n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Log1p)) { hasLog = true; }
                if (n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Sqrt)
                    || n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Sqrtabs)) { hasSqrt = true; }
                if (n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Div)) { hasDiv = true; }
                if (n.HashValue == static_cast<Operon::Hash>(Operon::BuiltinOp::Exp)) { hasExp = true; }
            }
        };
        scan(ind.x);
        scan(ind.y);
        if (hasLog) { ++nNonFiniteWithLog_; }
        if (hasSqrt) { ++nNonFiniteWithSqrt_; }
        if (hasDiv) { ++nNonFiniteWithDiv_; }
        if (hasExp) { ++nNonFiniteWithExp_; }
        if (!hasLog && !hasSqrt && !hasDiv && !hasExp) { ++nNonFiniteWithNoneOfThese_; }
    }

    void evaluateGuarded(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        try {
            evaluate(ind, ds, range);
            if (!std::isfinite(ind.shapeErr)
                || (config_.perTreeCredit
                    && (!std::isfinite(ind.errX) || !std::isfinite(ind.errY)))) {
                ++nNonFiniteFitness_;
                attributeNonFiniteOps(ind);
            }
            return;
        } catch (std::exception const&) {
            ++nEvalErrors_;
        } catch (...) {
            ++nEvalErrors_;
        }
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        ind.shapeErr = BIG;
        ind.errX = BIG;
        ind.errY = BIG;
        ind.regionErrCache.clear();
        ind.complexity = static_cast<Operon::Scalar>(ind.totalLength());
        ind.constBits = constMagnitudeComplexityMode()
            ? (constBitsOf(ind.x) + constBitsOf(ind.y)) : Operon::Scalar{0};
    }

    void evaluateGuardedFromGrid(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range,
                                 std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                                 CisPointSegBatch const* psBatch = nullptr, std::size_t psBatchIdx = 0,
                                 CisPointSegBatch const* psBatch48 = nullptr) const
    {
        try {
            evaluateFromGrid(ind, ds, range, gx, gy, psBatch, psBatchIdx, psBatch48);
            if (!std::isfinite(ind.shapeErr)
                || (config_.perTreeCredit
                    && (!std::isfinite(ind.errX) || !std::isfinite(ind.errY)))) {
                ++nNonFiniteFitness_;
                attributeNonFiniteOps(ind);
            }
            return;
        } catch (std::exception const&) {
            ++nEvalErrors_;
        } catch (...) {
            ++nEvalErrors_;
        }
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        ind.shapeErr = BIG;
        ind.errX = BIG;
        ind.errY = BIG;
        ind.regionErrCache.clear();
        ind.complexity = static_cast<Operon::Scalar>(ind.totalLength());
        ind.constBits = constMagnitudeComplexityMode()
            ? (constBitsOf(ind.x) + constBitsOf(ind.y)) : Operon::Scalar{0};
    }

    public:
    void evaluateBatchedPrepass(std::vector<PairIndividual>& pop, Operon::Dataset const& ds,
                                Operon::Range range,
                                std::vector<std::vector<Operon::Scalar>>& gridsX,
                                std::vector<std::vector<Operon::Scalar>>& gridsY) const
    {
        std::size_t const P = pop.size();
        gridsX.assign(P, {});
        gridsY.assign(P, {});
        if (P == 0) { return; }

        auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());

        std::vector<Sketch::RegisterProgram> progsX;
        std::vector<Sketch::RegisterProgram> progsY;
        std::vector<std::size_t> mapX;
        std::vector<std::size_t> mapY;
        progsX.reserve(P); progsY.reserve(P);
        mapX.reserve(P); mapY.reserve(P);
        std::vector<bool> okX(P, false), okY(P, false);

        auto needsScalarPrecision = [](std::span<Operon::Scalar const> coeffs) {
            for (auto const v : coeffs) {
                float const narrowed = static_cast<float>(v);
                if (!std::isfinite(v) || !std::isfinite(narrowed)) { return true; }
                if (std::abs(v - static_cast<Operon::Scalar>(narrowed))
                    > Operon::Scalar{1e-12} * (Operon::Scalar{1} + std::abs(v))) { return true; }
            }
            return false;
        };
        auto requiresScalarMath = [](BytecodeProgram const& bc) {
            for (auto const& inst : bc.insts) {
                if (inst.op == BytecodeOp::Div || inst.op == BytecodeOp::Exp
                    || inst.op == BytecodeOp::Log || inst.op == BytecodeOp::Sqrt
                    || inst.op == BytecodeOp::Square || inst.op == BytecodeOp::Mul) { return true; }
            }
            return false;
        };

        for (std::size_t i = 0; i < P; ++i) {
            auto cx = pop[i].x.GetCoefficients();
            auto bcX = needsScalarPrecision({cx.data(), cx.size()})
                ? BytecodeProgram{} : BytecodeProgram::Compile(pop[i].x, {cx.data(), cx.size()});
            if (bcX.valid && !requiresScalarMath(bcX)) {
                auto rpX = Sketch::RegisterProgram::FromBytecode(bcX);
                if (rpX.valid) { okX[i] = true; mapX.push_back(i); progsX.push_back(std::move(rpX)); }
            }
            auto cy = pop[i].y.GetCoefficients();
            auto bcY = needsScalarPrecision({cy.data(), cy.size()})
                ? BytecodeProgram{} : BytecodeProgram::Compile(pop[i].y, {cy.data(), cy.size()});
            if (bcY.valid && !requiresScalarMath(bcY)) {
                auto rpY = Sketch::RegisterProgram::FromBytecode(bcY);
                if (rpY.valid) { okY[i] = true; mapY.push_back(i); progsY.push_back(std::move(rpY)); }
            }
        }

        std::vector<std::vector<Operon::Scalar>> outX, outY;
        if (!progsX.empty()) { Sketch::BatchedRegisterInterp::EvaluatePopulation(progsX, tValues, outX); }
        if (!progsY.empty()) { Sketch::BatchedRegisterInterp::EvaluatePopulation(progsY, tValues, outY); }

        for (std::size_t j = 0; j < mapX.size(); ++j) { gridsX[mapX[j]] = std::move(outX[j]); }
        for (std::size_t j = 0; j < mapY.size(); ++j) { gridsY[mapY[j]] = std::move(outY[j]); }

        for (std::size_t i = 0; i < P; ++i) {
            if (!okX[i]) {
                auto cx = pop[i].x.GetCoefficients();
                interpInto(pop[i].x, cx, ds, range, gridsX[i]);
            }
            if (!okY[i]) {
                auto cy = pop[i].y.GetCoefficients();
                interpInto(pop[i].y, cy, ds, range, gridsY[i]);
            }
        }
    }
