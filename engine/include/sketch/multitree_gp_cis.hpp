
    auto cisShapeErrorForTest(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                              bool directFin, CisPointSegBatch const* psBatch = nullptr,
                              std::size_t psBatchIdx = 0) const -> Operon::Scalar
    {
        return cisShapeError(gx, gy, directFin, psBatch, psBatchIdx);
    }

    auto reverseHalfCpuForTest(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy) const
        -> Operon::Scalar
    {
        auto const M = gx.size();
        if (M < 2) { return 0; }
        thread_local std::vector<Operon::Scalar> cisCx, cisCy;
        cisCx.resize(M); cisCy.resize(M);
        Operon::Scalar gx0 = gx[0], gy0 = gy[0];
        Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
        cisCx[0] = gx0; cisCy[0] = gy0;
        for (std::size_t j = 1; j < M; ++j) {
            Operon::Scalar const vx = gx[j], vy = gy[j];
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
        return oneResReverseHalfSum(cisCx, cisCy, M);
    }

    struct LastGenDumpForTest {
        bool valid{false};
        std::vector<std::vector<Operon::Scalar>> gridsX, gridsY;
        std::vector<Operon::Scalar> gpuReverseFitness;
    };
    [[nodiscard]] auto lastGenDumpForTest() const -> LastGenDumpForTest
    {
        return {dbgLastValid_, dbgLastGridsX_, dbgLastGridsY_, dbgLastReverseFitness_};
    }

    [[nodiscard]] auto cisGrid48StrokeSizeForTest() const -> std::size_t { return cisGrid48Stroke_.size(); }
    [[nodiscard]] auto cisGrid48OkForTest() const -> bool { return cisGrid48Ok_; }
    [[nodiscard]] auto cisGrid48Zi5OkForTest() const -> bool { return cisGrid48Zi5Ok_; }
    [[nodiscard]] auto cisGrid48Zi5CellCountForTest() const -> std::size_t { return cisGrid48Zi5LnE_.size(); }
    [[nodiscard]] auto cisGrid48StrokeForTest() const -> std::vector<Vec2> const& { return cisGrid48Stroke_; }
    [[nodiscard]] auto cisGrid48StrokeDataForTest() const -> std::vector<Vec2> const& { return cisGrid48StrokeData_; }
    [[nodiscard]] auto cisGrid48ScaleForTest() const -> double { return cisGrid48Scale_; }

    struct OneResidualSitesForTest {
        Operon::Scalar cis{0};
        bool cacheValid{false};
        std::vector<Operon::Scalar> fitRes, fitTau;
        bool consumerOk{false};
        std::vector<Operon::Scalar> consumerRes, consumerTau;
        long long hitsBefore{0}, hitsAfter{0};
        std::vector<Operon::Scalar> recRes;
        std::vector<Operon::Scalar> nmlRes;
        Operon::Scalar monotonicity{0};
    };
    auto oneResidualSitesForTest(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                                 CisPointSegBatch const* b = nullptr, std::size_t idx = 0) const
        -> OneResidualSitesForTest
    {
        OneResidualSitesForTest o;
        o.cis = cisShapeError(gx, gy, false, b, idx);
        o.cacheValid = oneResCache_.valid && oneResCache_.owner == this;
        if (o.cacheValid) { o.fitRes = oneResCache_.residual; o.fitTau = oneResCache_.tau; }
        o.hitsBefore = nOneResCacheHit_.load();
        o.consumerOk = oneResidualOf(gx, gy, true, o.consumerRes, &o.consumerTau);
        o.hitsAfter = nOneResCacheHit_.load();
        o.monotonicity = monotonicityOrderTerm(gx, gy);
        auto const M = gx.size();
        auto const N = dataX_.size();
        std::vector<Vec2> strokePts(N), candPts(M);
        {
            Operon::Scalar xlo = gx[0], xhi = gx[0], ylo = gy[0], yhi = gy[0];
            for (std::size_t j = 1; j < M; ++j) {
                if (gx[j] < xlo) { xlo = gx[j]; } else if (gx[j] > xhi) { xhi = gx[j]; }
                if (gy[j] < ylo) { ylo = gy[j]; } else if (gy[j] > yhi) { yhi = gy[j]; }
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) { candPts[j] = Vec2{(gx[j] - cxc) / sp, (gy[j] - cyc) / sp}; }
        }
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        o.recRes = pointSegOp(strokePts, candPts, false, false).residual;
        o.nmlRes = pointSegOp(strokePts, candPts, true, true).residual;
        return o;
    }

    struct OneResidualEvalForTest {
        Operon::Scalar shapeErr{0};
        std::vector<Operon::Scalar> fitRes;
        std::vector<Operon::Scalar> retained;
        std::vector<Operon::Scalar> lexCases;
        std::vector<std::size_t> arcIndex;
        Operon::Scalar lexS2{1};
    };
    auto oneResidualEvalForTest(PairIndividual ind) const -> OneResidualEvalForTest
    {
        OneResidualEvalForTest o;
        if (!state_ || !state_->fullDs) { return o; }
        Operon::Range const r{0, state_->M};
        CisFitnessModeGuard cisFitnessGuard(config_.useCisFitness);
        evaluateGuarded(ind, *state_->fullDs, r);
        o.shapeErr = ind.shapeErr;
        if (oneResCache_.valid && oneResCache_.owner == this) { o.fitRes = oneResCache_.residual; }
        o.retained = ind.perPointResid_;
        o.lexCases = ind.regionErrCache;
        o.lexS2 = config_.useNoiseWeightedLexicase ? noiseMeanS2_ : Operon::Scalar{1};
        o.arcIndex.resize(dataX_.size());
        for (std::size_t i = 0; i < dataX_.size(); ++i) { o.arcIndex[i] = originalArcIndex(i); }
        return o;
    }
    [[nodiscard]] auto oneResCacheCountersForTest() const -> std::pair<long long, long long>
    {
        return {nOneResCacheHit_.load(), nOneResCacheMiss_.load()};
    }

    auto orderTermForTest(std::vector<Operon::Scalar> const& gx,
                          std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        return orderTerm(gx, gy);
    }

    auto frontObjectiveFromGrid(std::vector<Operon::Scalar> const& gx,
                                std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        CreditModeGuard creditGuard(config_.perTreeCredit);
        FloorClampGuard clampGuard(config_.floorClamp);
        FloorModeGuard floorModeGuard(config_.floorOnChamfer);
        CisFitnessModeGuard cisFitnessGuard(config_.useCisFitness);
        PointClampGuard pointClampGuard(config_.pointClamp);
        ConstMagnitudeComplexityGuard constMagGuard(config_.constMagnitudeComplexity);
        PriorMinimalComplexityGuard priorMinimalGuard(config_.usePriorMinimalComplexity);
        FastNonDominatedSortGuard fastNdsGuard(config_.useFastNonDominatedSort);
        UseExactNlogNSortGuard exactNlogNGuard(config_.useExactNlogNNondominatedSort);
        ArcLengthGuardRatioGuard arcGuard(config_.arcLengthGuardRatio);
        DirectFinalistEvalGuard dfe(true);
        Operon::Scalar chamfer{0};
        return shapeOnlyPartsFromGrid(gx, gy, &chamfer);
    }

    struct RawFrontGridForTest {
        std::vector<Operon::Scalar> gx, gy;
        std::string exprX, exprY;
    };
    auto rawFrontGridForTest(Operon::Tree const& x, Operon::Tree const& y) const -> RawFrontGridForTest
    {
        RawFrontGridForTest out;
        if (!state_ || !state_->fullDs) { return out; }
        Operon::Range const r{0, state_->M};
        auto const cx = x.GetCoefficients();
        auto const cy = y.GetCoefficients();
        out.gx = interp(x, cx, *state_->fullDs, r);
        out.gy = interp(y, cy, *state_->fullDs, r);
        out.exprX = Operon::InfixFormatter::Format(x, *state_->fullDs, 4);
        out.exprY = Operon::InfixFormatter::Format(y, *state_->fullDs, 4);
        return out;
    }

    [[nodiscard]] auto dataXForTest() const -> std::vector<Operon::Scalar> const& { return dataX_; }
    [[nodiscard]] auto dataYForTest() const -> std::vector<Operon::Scalar> const& { return dataY_; }

    auto evaluateBatchedPointSegPrepass(std::vector<PairIndividual>& pop, Operon::Dataset const& ds,
                                        Operon::Range range,
                                        std::vector<std::vector<Operon::Scalar>>& gridsX,
                                        std::vector<std::vector<Operon::Scalar>>& gridsY,
                                        CisPointSegBatch& batch,
                                        CisPointSegBatch* batch48 = nullptr) const -> bool
    {
        return evaluateBatchedPointSegPrepassImpl(pop.size(),
            [&pop](std::size_t i) -> PairIndividual const& { return pop[i]; },
            ds, range, gridsX, gridsY, batch, batch48);
    }

    auto evaluateBatchedPointSegPrepassPtrs(std::vector<PairIndividual*> const& inds, Operon::Dataset const& ds,
                                            Operon::Range range,
                                            std::vector<std::vector<Operon::Scalar>>& gridsX,
                                            std::vector<std::vector<Operon::Scalar>>& gridsY,
                                            CisPointSegBatch& batch,
                                            CisPointSegBatch* batch48 = nullptr) const -> bool
    {
        return evaluateBatchedPointSegPrepassImpl(inds.size(),
            [&inds](std::size_t i) -> PairIndividual const& { return *inds[i]; },
            ds, range, gridsX, gridsY, batch, batch48);
    }

    template <class Get>
    auto evaluateBatchedPointSegPrepassImpl(std::size_t const P, Get const& get, Operon::Dataset const& ds,
                                            Operon::Range range,
                                            std::vector<std::vector<Operon::Scalar>>& gridsX,
                                            std::vector<std::vector<Operon::Scalar>>& gridsY,
                                            CisPointSegBatch& batch,
                                            CisPointSegBatch* batch48 = nullptr) const -> bool
    {
        if (batch48 != nullptr) { *batch48 = CisPointSegBatch{}; }
        bool const dbg = std::getenv("SKETCH_T15_DEBUG") != nullptr;
        auto fail = [&](char const* why) {
            if (dbg) { std::fprintf(stderr, "[T15] evaluateBatchedPointSegPrepass FAIL: %s\n", why); }
            if (config_.useBatchedGpuPointSeg) { pointSegGpuCounters().fallbackBatches.fetch_add(1, std::memory_order_relaxed); }
            return false;
        };

        batch = CisPointSegBatch{};
        gridsX.assign(P, {});
        gridsY.assign(P, {});
        if (!config_.useBatchedGpuPointSeg || P == 0) { return fail("flag off or empty pop"); }

        for (std::size_t i = 0; i < P; ++i) {
            auto cx = get(i).x.GetCoefficients();
            auto cy = get(i).y.GetCoefficients();
            interpInto(get(i).x, cx, ds, range, gridsX[i]);
            interpInto(get(i).y, cy, ds, range, gridsY[i]);
        }

        std::size_t const M = gridsX.empty() ? 0 : gridsX[0].size();
        if (M < 2) { return fail("M < 2"); }
        for (std::size_t i = 0; i < P; ++i) {
            if (gridsX[i].size() != M || gridsY[i].size() != M) { return fail("non-uniform M across population"); }
        }

        std::vector<std::vector<Vec2>> candPolys(P);
        std::vector<bool> validFlags(P, true);
        std::size_t numInvalid = 0;
        for (std::size_t i = 0; i < P; ++i) {
            auto const& gx = gridsX[i];
            auto const& gy = gridsY[i];
            bool ok = true;
            Operon::Scalar gx0 = gx[0], gy0 = gy[0];
            if (!std::isfinite(gx0) || !std::isfinite(gy0)) { ok = false; }
            Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
            std::vector<Operon::Scalar> cx(M), cy(M);
            cx[0] = ok ? gx0 : Operon::Scalar{0}; cy[0] = ok ? gy0 : Operon::Scalar{0};
            for (std::size_t j = 1; j < M && ok; ++j) {
                Operon::Scalar const vx = gx[j], vy = gy[j];
                if (!std::isfinite(vx) || !std::isfinite(vy)) { ok = false; break; }
                cx[j] = vx; cy[j] = vy;
                if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
                if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
            }
            candPolys[i].resize(M);
            if (!ok) {
                validFlags[i] = false;
                ++numInvalid;
                for (std::size_t j = 0; j < M; ++j) { candPolys[i][j] = Vec2{0, 0}; }
                continue;
            }
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                candPolys[i][j] = Vec2{(cx[j] - cxc) / sp, (cy[j] - cyc) / sp};
            }
        }
        if (numInvalid == P) { return fail("every candidate degenerate"); }

        std::size_t const N = dataX_.size();
        if (N < 2) { return fail("N < 2"); }
        std::vector<Vec2> strokePts(N);
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX_[i], dataY_[i]}; }

        if (config_.usePersistentGpuScoring && !psResident_) { psResident_ = std::make_unique<PointSegGpuResident>(); }
#ifdef SKETCH_P0_SHARED_RESIDENT_CONTROL
        auto gpuResult = config_.usePersistentGpuScoring
            ? evaluatePopulationPointSegGpuResident(strokePts, candPolys)
            : evaluatePopulationPointSegGpu(strokePts, candPolys);
#else
        if (psResident_) { psResident_->unsignedMag = oneResUnsigned(); }
        auto gpuResult = config_.usePersistentGpuScoring
            ? evaluatePopulationPointSegGpuResident(strokePts, candPolys, *psResident_)
            : evaluatePopulationPointSegGpu(strokePts, candPolys, oneResUnsigned());
#endif
        if (!gpuResult.ok) { return fail("evaluatePopulationPointSegGpu ok=false"); }

        if (batch48 != nullptr && config_.useB3Grid48Gpu && cisGrid48Bit0() && cisGrid48Ok_) {
            auto gpuResult48 = evaluatePopulationPointSegGpu(cisGrid48Target(), candPolys, oneResUnsigned());
            if (gpuResult48.ok) {
                batch48->N = cisGrid48Stroke_.size();
                batch48->fp16ModeUsed = gpuResult48.fp16ModeUsed;
                batch48->fitness = std::move(gpuResult48.fitness);
                batch48->residualFlat = std::move(gpuResult48.residualFlat);
                batch48->tauFlat = std::move(gpuResult48.tauFlat);
                batch48->valid = validFlags;
            }
        }

        batch.N = N;
        batch.fp16ModeUsed = gpuResult.fp16ModeUsed;
        batch.fitness = std::move(gpuResult.fitness);
        batch.residualFlat = std::move(gpuResult.residualFlat);
        batch.tauFlat = std::move(gpuResult.tauFlat);
        batch.valid = std::move(validFlags);

        if (config_.useGpuReverseHalf && kPointSegGpuCompiledIn && oneResBit(2)) {
            std::vector<Vec2> drawn(N);
            for (std::size_t i = 0; i < N; ++i) { drawn[originalArcIndex(i)] = Vec2{dataX_[i], dataY_[i]}; }
            auto const Kn = noise_.size();
            std::vector<Operon::Scalar> const s2Copy = Kn > 0 ? noise_.s2 : std::vector<Operon::Scalar>{};
            std::vector<Operon::Scalar> const nuCopy = Kn > 0 ? noise_.nu : std::vector<Operon::Scalar>{};
            auto revResult = evaluatePopulationPointSegReverseGpu(drawn, candPolys, s2Copy, nuCopy, pointClampMode());
            if (revResult.ok) { batch.reverseFitness = std::move(revResult.fitness); }
            static bool const dumpGrids = std::getenv("SKETCH_DUMP_LAST_GEN_GRIDS") != nullptr;
            if (dumpGrids) {
                dbgLastGridsX_ = gridsX;
                dbgLastGridsY_ = gridsY;
                dbgLastReverseFitness_ = batch.reverseFitness;
                dbgLastValid_ = true;
            }
        }
        {
            auto& ctr = pointSegGpuCounters();
            ctr.gpuBatches.fetch_add(1, std::memory_order_relaxed);
            ctr.gpuCandidates.fetch_add(static_cast<long long>(P - numInvalid), std::memory_order_relaxed);
            ctr.cpuInGpuBatch.fetch_add(static_cast<long long>(numInvalid), std::memory_order_relaxed);
        }
        if (dbg) {
            std::fprintf(stderr, "[T15] evaluateBatchedPointSegPrepass OK P=%zu M=%zu N=%zu invalid=%zu\n",
                         P, M, N, numInvalid);
        }
        return true;
    }

    struct PsStreamGen {
        bool active{false};
        bool overlap{false};
        std::size_t added{0};
        std::size_t numInvalid{0};
        std::vector<bool> valid;
    };

    struct PsOverlapHelper {
        std::function<void()> body;
        std::atomic<std::size_t> published{0};
        std::atomic<bool> done{false};
        std::atomic<bool> failed{false};
        bool started{false};
#ifndef _WIN32
        pthread_t tid{};
#else
        std::thread th;
#endif
        PsOverlapHelper() = default;
        PsOverlapHelper(PsOverlapHelper const&) = delete;
        auto operator=(PsOverlapHelper const&) -> PsOverlapHelper& = delete;
        ~PsOverlapHelper() { if (started) { (void)finish(); } }

        template <class F>
        void start(PairIndividual const* base, F fn)
        {
            body = [this, base, fn]() mutable {
#ifdef __APPLE__
                pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
                std::size_t next = 0;
                try {
                    for (;;) {
                        std::size_t const avail = published.load(std::memory_order_acquire);
                        if (next < avail) {
                            for (; next < avail; ++next) { fn(base[next], next); }
                            continue;
                        }
                        if (done.load(std::memory_order_acquire)) {
                            if (next >= published.load(std::memory_order_acquire)) { break; }
                            continue;
                        }
                        std::this_thread::yield();
                    }
                } catch (...) {
                    failed.store(true, std::memory_order_release);
                }
            };
#ifndef _WIN32
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setstacksize(&attr, std::size_t{16} << 20);
            started = pthread_create(&tid, &attr, &MultitreeGP::workerTrampoline, &body) == 0;
            pthread_attr_destroy(&attr);
#else
            th = std::thread(body);
            started = true;
#endif
        }
        [[nodiscard]] auto running() const -> bool { return started; }
        void publish(std::size_t n) { published.store(n, std::memory_order_release); }
        auto finish() -> bool
        {
            if (!started) { return !failed.load(); }
            done.store(true, std::memory_order_release);
#ifndef _WIN32
            pthread_join(tid, nullptr);
#else
            th.join();
#endif
            started = false;
            return !failed.load(std::memory_order_acquire);
        }
    };

    [[nodiscard]] auto psStreamWanted() const -> bool
    {
        return kPointSegGpuStreamAvailable && config_.useBatchedGpuPointSeg
            && (config_.gpuPointSegZeroCopy || config_.gpuPointSegOverlap) && cisFitnessMode();
    }

    auto psStreamBegin(std::size_t pExpected, Operon::Range range, PsStreamGen& st) const -> bool
    {
        st = PsStreamGen{};
        if (!psStreamWanted() || pExpected == 0) { return false; }
        std::size_t const N = dataX_.size();
        std::size_t const M = range.Size();
        if (N < 2 || M < 2) { return false; }
        if (!psStream_) { psStream_ = std::make_unique<PointSegGpuStream>(); }
        PointSegGpuStreamOptions o;
        o.unretainedRefs = config_.gpuPointSegUnretained;
        o.spinWait = config_.gpuPointSegSpinWait;
        o.residencySet = config_.gpuPointSegResidencySet;
        o.unsignedMag = oneResUnsigned();
        if (!psStream_->begin(dataX_.data(), dataY_.data(), N, pExpected, M, o)) { return false; }
        st.active = true;
        st.overlap = config_.gpuPointSegOverlap;
        st.valid.assign(pExpected, true);
        return true;
    }

    void psStreamAdd(PairIndividual const& ind, std::size_t i, Operon::Dataset const& ds, Operon::Range range,
                     std::vector<std::vector<Operon::Scalar>>& gridsX,
                     std::vector<std::vector<Operon::Scalar>>& gridsY, PsStreamGen& st) const
    {
        if (gridsX.size() <= i) { gridsX.resize(i + 1); gridsY.resize(i + 1); }
        auto cx0 = ind.x.GetCoefficients();
        auto cy0 = ind.y.GetCoefficients();
        interpInto(ind.x, cx0, ds, range, gridsX[i]);
        interpInto(ind.y, cy0, ds, range, gridsY[i]);
        auto const& gx = gridsX[i];
        auto const& gy = gridsY[i];
        std::size_t const M = gx.size();
        auto* slot = psStream_->candSlot(i);
        bool ok = true;
        Operon::Scalar gx0 = gx[0], gy0 = gy[0];
        if (!std::isfinite(gx0) || !std::isfinite(gy0)) { ok = false; }
        Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
        for (std::size_t j = 1; j < M && ok; ++j) {
            Operon::Scalar const vx = gx[j], vy = gy[j];
            if (!std::isfinite(vx) || !std::isfinite(vy)) { ok = false; break; }
            if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
            if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
        }
        if (!ok) {
            st.valid[i] = false;
            ++st.numInvalid;
            for (std::size_t j = 0; j < M; ++j) { slot[j] = {0.0f, 0.0f}; }
            psStream_->setMeanSegLen(i, 0.0f);
        } else {
            Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
            Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
            for (std::size_t j = 0; j < M; ++j) {
                slot[j] = {static_cast<float>((gx[j] - cxc) / sp), static_cast<float>((gy[j] - cyc) / sp)};
            }
            psStream_->setMeanSegLen(i, meanSegLenOfSpan(slot, M));
        }
        st.added = i + 1;
    }

    void psStreamPump(PsStreamGen& st) const
    {
        if (!st.active || !st.overlap) { return; }
        auto const c = psStream_->committedUpTo();
        if (st.added <= c) { return; }
        std::size_t const pending = st.added - c;
        std::size_t const remaining = st.valid.size() - st.added;
        bool commit = false;
        switch (config_.gpuPointSegOverlapSchedule) {
        case 1: commit = pending >= remaining; break;
        case 2: commit = c == 0 && st.added * 2 >= st.valid.size(); break;
        default: commit = !psStream_->busy(); break;
        }
        if (commit) { psStream_->commitRange(c, st.added); }
    }

    auto psStreamFinish(std::size_t P, PsStreamGen& st, CisPointSegBatch& batch) const -> bool
    {
        batch = CisPointSegBatch{};
        auto& ctr = pointSegGpuCounters();
        auto fail = [&] { ctr.fallbackBatches.fetch_add(1, std::memory_order_relaxed); return false; };
        if (!st.active || st.added != P || P == 0) { return fail(); }
        auto const c = psStream_->committedUpTo();
        if (P > c && !psStream_->commitRange(c, P)) { return fail(); }
        if (st.numInvalid == P) { return fail(); }
        bool const joinBeforeScoring = !st.overlap || config_.gpuPointSegOverlapSchedule == 3;
        if (joinBeforeScoring && !psStream_->waitAll()) { return fail(); }
        batch.N = dataX_.size();
        batch.P = P;
        batch.stream = psStream_.get();
        batch.fp16ModeUsed = pointSegFp16ModeForFloatBuffers();
        batch.valid = std::move(st.valid);
        ctr.gpuBatches.fetch_add(1, std::memory_order_relaxed);
        ctr.gpuCandidates.fetch_add(static_cast<long long>(P - st.numInvalid), std::memory_order_relaxed);
        ctr.cpuInGpuBatch.fetch_add(static_cast<long long>(st.numInvalid), std::memory_order_relaxed);
        return true;
    }

    private:
    mutable std::unique_ptr<PointSegGpuStream> psStream_;
    mutable std::unique_ptr<PointSegGpuResident> psResident_;

    mutable bool dbgLastValid_{false};
    mutable std::vector<std::vector<Operon::Scalar>> dbgLastGridsX_, dbgLastGridsY_;
    mutable std::vector<Operon::Scalar> dbgLastReverseFitness_;

    auto cclmStageOk(SearchState const& s) -> bool
    {
        if (!(config_.constOptimizer == 1 && !config_.useLm && !config_.guardConstOpt
              && config_.lmDownsampleK == 0 && config_.useBytecodeEvaluator && !config_.useGridSdf
              && !config_.useAmxChamfer && !config_.useLiveCoarseSieve
              && !(arcLengthGuardRatio() > Operon::Scalar{0}))) { return false; }
        std::size_t const N = dataX_.size(), M = s.range.Size();
        if (N < 2 || M < 2 || N > kCclmMaxPoints || M > kCclmMaxPoints || s.activeDs == nullptr) { return false; }
        if (!cclmGpuTried_) { cclmGpuTried_ = true; cclmGpu_ = cclmMakeRuntime(); }
        return cclmGpu_ != nullptr;
    }

    void runConcurrentLmStage(std::vector<PairIndividual>& offspring, std::vector<std::size_t> const& order,
                              std::size_t kEff, SearchState& s, std::size_t workers, int gen)
    {
        using clk = std::chrono::steady_clock;
        auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        auto const tStage = clk::now();
        ConcurrentLmGenStat st;
        st.gen = gen;
        st.k = kEff;
        auto const& ds = *s.activeDs;
        auto const range = s.range;

        double share = 0.5;
        if (config_.concurrentLmForceGpuShare >= 0.0) {
            share = std::min(1.0, config_.concurrentLmForceGpuShare);
        } else if (cclmRateCpu_ > 0.0 && cclmRateGpu_ > 0.0) {
            share = cclmRateGpu_ / (cclmRateCpu_ + cclmRateGpu_);
        }
        st.gpuShare = share;
        std::size_t const nGpuTarget = static_cast<std::size_t>(std::lround(share * static_cast<double>(kEff)));
        st.nGpuTarget = nGpuTarget;

        std::size_t const N = dataX_.size(), M = range.Size();
        auto& fit = cclmFit_;
        auto const tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
        fit.t.assign(tValues.begin(), tValues.end());
        fit.dataXY.resize(2 * N);
        for (std::size_t i = 0; i < N; ++i) {
            fit.dataXY[2 * i] = static_cast<float>(dataX_[i]);
            fit.dataXY[2 * i + 1] = static_cast<float>(dataY_[i]);
        }
        std::size_t const K = noise_.size();
        fit.noise4.assign(4 * N, 0.0f);
        if (K > 0) {
            auto const& rlut = noiseRegionLUT(N, K);
            for (std::size_t i = 0; i < N; ++i) {
                std::size_t const kl = rlut[i];
                std::size_t const kb = noiseIndexAtArc(originalArcIndex(i), N, K);
                fit.noise4[4 * i + 0] = static_cast<float>(noise_.s2[kl]);
                fit.noise4[4 * i + 1] = static_cast<float>(noise_.nu[kl]);
                fit.noise4[4 * i + 2] = static_cast<float>(noise_.s2[kb]);
                fit.noise4[4 * i + 3] = static_cast<float>(noise_.nu[kb]);
            }
        }
        LevenbergMarquardt const lmDefaults;
        auto& P = fit.params;
        P.N = static_cast<uint32_t>(N);
        P.M = static_cast<uint32_t>(M);
        P.hasNoise = K > 0 ? 1u : 0u;
        P.ptClamp = pointClampMode() ? 1u : 0u;
        P.maxIter = static_cast<uint32_t>(lmDefaults.maxIter);
        P.maxTrials = static_cast<uint32_t>(lmDefaults.maxTrials);
        P.screening = lmDefaults.enableScreening ? 1u : 0u;
        P.ftol = static_cast<float>(lmDefaults.ftol);
        P.lambda0 = static_cast<float>(lmDefaults.lambda0);
        P.lambdaUp = static_cast<float>(lmDefaults.lambdaUp);
        P.lambdaDown = static_cast<float>(lmDefaults.lambdaDown);
        P.lambdaMin = static_cast<float>(lmDefaults.lambdaMin);
        P.lambdaMax = static_cast<float>(lmDefaults.lambdaMax);
        P.floorThr = config_.floorClamp > 0.0f ? static_cast<float>(config_.floorClamp)
                   : (config_.floorStop > 0.0f ? static_cast<float>(config_.floorStop) : 0.0f);

        if (config_.concurrentLmSchedule == 1 && config_.concurrentLmForceGpuShare < 0.0) {
            runConcurrentLmRace(offspring, order, kEff, s, workers, st, tStage);
            return;
        }

        std::vector<char> onGpu(kEff, 0);
        cclmGpuPos_.clear();
        std::size_t nJobs = 0;
        for (std::size_t jj = kEff; jj-- > 0 && cclmGpuPos_.size() < nGpuTarget;) {
            auto const& ind = offspring[order[jj]];
            auto cx = ind.x.GetCoefficients();
            auto cy = ind.y.GetCoefficients();
            std::size_t const n = cx.size() + cy.size();
            if (n == 0 || n > kCclmMaxCoeffs) { continue; }
            auto progX = BytecodeProgram::Compile(ind.x, {cx.data(), cx.size()});
            auto progY = BytecodeProgram::Compile(ind.y, {cy.data(), cy.size()});
            if (!progX.valid || !progY.valid || progX.bindings.size() != cx.size() || progY.bindings.size() != cy.size()) { continue; }
            if (cclmJobs_.size() <= nJobs) { cclmJobs_.emplace_back(); }
            auto& job = cclmJobs_[nJobs];
            if (!cclmBuildTape(progX, 0, job.tapeX) || !cclmBuildTape(progY, static_cast<int>(cx.size()), job.tapeY)) { continue; }
            job.nx = static_cast<uint32_t>(cx.size());
            job.coeffs.clear();
            for (auto v : cx) { job.coeffs.push_back(static_cast<float>(v)); }
            for (auto v : cy) { job.coeffs.push_back(static_cast<float>(v)); }
            ++nJobs;
            cclmGpuPos_.push_back(jj);
            onGpu[jj] = 1;
        }
        std::vector<PairIndividual> shadow, shadowUlp;
        if (config_.concurrentLmParityShadow) {
            shadow.reserve(cclmGpuPos_.size());
            for (auto jj : cclmGpuPos_) { shadow.push_back(offspring[order[jj]]); }
            shadowUlp = shadow;
            for (auto& sh : shadowUlp) {
                auto cx = sh.x.GetCoefficients();
                auto cy = sh.y.GetCoefficients();
                for (auto& v : cx) { v = std::nextafter(v, std::numeric_limits<Operon::Scalar>::infinity()); }
                for (auto& v : cy) { v = std::nextafter(v, std::numeric_limits<Operon::Scalar>::infinity()); }
                sh.x.SetCoefficients(cx);
                sh.y.SetCoefficients(cy);
            }
        }
        bool submitted = false;
        if (config_.concurrentLmKernelVariant >= 0) { cclmSetVariant(*cclmGpu_, config_.concurrentLmKernelVariant); }
        if (nJobs > 0 && cclmAllocEvents(*cclmGpu_) == 0) { cclmReserve(*cclmGpu_, kEff); }
        if (nJobs > 0) {
            submitted = cclmSubmit(*cclmGpu_, fit, cclmJobs_, nJobs);
        }
        if (!submitted) {
            std::fill(onGpu.begin(), onGpu.end(), 0);
            cclmGpuPos_.clear();
            shadow.clear();
            shadowUlp.clear();
        }
        cclmCpuPos_.clear();
        for (std::size_t jj = 0; jj < kEff; ++jj) { if (onGpu[jj] == 0) { cclmCpuPos_.push_back(jj); } }
        st.nGpu = cclmGpuPos_.size();
        st.nCpu = cclmCpuPos_.size();
        auto const tSubmit = clk::now();
        st.prepMs = ms(tStage, tSubmit);

        parallelApply(cclmCpuPos_.size(), workers, [&](std::size_t q) {
            std::size_t const i = order[cclmCpuPos_[q]];
            optConstsInLoop(offspring[i], ds, range);
            evaluateGuarded(offspring[i], ds, range);
        });
        auto const tCpuDone = clk::now();
        st.cpuMs = ms(tSubmit, tCpuDone);

        std::vector<char> fellBack(kEff, 0);
        if (submitted) {
            int64_t doneNs = 0;
            double kernelMs = 0.0;
            bool const ok = cclmWait(*cclmGpu_, cclmOuts_, doneNs, kernelMs);
            auto const tGpuDone = ok ? clk::time_point(std::chrono::duration_cast<clk::duration>(std::chrono::nanoseconds(doneNs)))
                                     : clk::now();
            st.gpuMs = ms(tSubmit, tGpuDone);
            st.gpuKernelMs = kernelMs;
            if (tGpuDone <= tCpuDone) { st.idleGpuMs = ms(tGpuDone, tCpuDone); }
            else { st.idleCpuMs = ms(tCpuDone, tGpuDone); }
            for (std::size_t q = 0; q < cclmGpuPos_.size(); ++q) {
                std::size_t const jj = cclmGpuPos_[q];
                auto& ind = offspring[order[jj]];
                bool good = ok && q < cclmOuts_.size() && cclmOuts_[q].status == 0;
                if (good) {
                    auto const& o = cclmOuts_[q];
                    for (auto v : o.coeffs) { if (!std::isfinite(v)) { good = false; break; } }
                    if (good) {
                        std::size_t const nx = cclmJobs_[q].nx;
                        ind.x.SetCoefficients(std::vector<Operon::Scalar>(o.coeffs.begin(), o.coeffs.begin() + static_cast<long>(nx)));
                        ind.y.SetCoefficients(std::vector<Operon::Scalar>(o.coeffs.begin() + static_cast<long>(nx), o.coeffs.end()));
                        ++nConstDescents_;
                        nConstLossCalls_ += static_cast<std::size_t>(o.lossCalls);
                        nConstIters_ += static_cast<std::size_t>(o.iters);
                    }
                }
                if (!good) { fellBack[jj] = 1; ++st.nFallback; }
            }
            parallelApply(cclmGpuPos_.size(), workers, [&](std::size_t q) {
                std::size_t const jj = cclmGpuPos_[q];
                std::size_t const i = order[jj];
                if (fellBack[jj] != 0) { optConstsInLoop(offspring[i], ds, range); }
                evaluateGuarded(offspring[i], ds, range);
            });
            auto const tPost = clk::now();
            st.postMs = ms(tCpuDone > tGpuDone ? tCpuDone : tGpuDone, tPost);
            std::size_t const nGpuRan = st.nGpu - st.nFallback;
            if (st.gpuMs > 0.0 && nGpuRan > 0) { st.rateGpu = static_cast<double>(nGpuRan) / (st.gpuMs / 1e3); }
        }
        if (st.cpuMs > 0.0 && st.nCpu > 0) { st.rateCpu = static_cast<double>(st.nCpu) / (st.cpuMs / 1e3); }
        if (st.rateCpu > 0.0) { cclmRateCpu_ = st.rateCpu; }
        if (st.rateGpu > 0.0) { cclmRateGpu_ = st.rateGpu; }
        st.stageMs = ms(tStage, clk::now());
        cclmStats_.push_back(st);

        if (config_.concurrentLmParityShadow && !shadow.empty()) {
            parallelApply(shadow.size(), workers, [&](std::size_t q) {
                optConstsInLoop(shadow[q], ds, range);
                evaluateGuarded(shadow[q], ds, range);
                optConstsInLoop(shadowUlp[q], ds, range);
                evaluateGuarded(shadowUlp[q], ds, range);
            });
            ConcurrentLmParityGen pg;
            pg.gen = gen;
            pg.mixed.resize(kEff);
            pg.cpuRef.resize(kEff);
            pg.onGpu = onGpu;
            pg.fellBack = fellBack;
            for (std::size_t jj = 0; jj < kEff; ++jj) {
                pg.mixed[jj] = offspring[order[jj]].shapeErr;
                pg.cpuRef[jj] = offspring[order[jj]].shapeErr;
            }
            pg.cpuUlp = pg.cpuRef;
            for (std::size_t q = 0; q < cclmGpuPos_.size(); ++q) {
                pg.cpuRef[cclmGpuPos_[q]] = shadow[q].shapeErr;
                pg.cpuUlp[cclmGpuPos_[q]] = shadowUlp[q].shapeErr;
            }
            cclmParity_.push_back(std::move(pg));
        }
    }

    void runConcurrentLmRace(std::vector<PairIndividual>& offspring, std::vector<std::size_t> const& order,
                             std::size_t kEff, SearchState& s, std::size_t workers, ConcurrentLmGenStat& st,
                             std::chrono::steady_clock::time_point tStage)
    {
        using clk = std::chrono::steady_clock;
        auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        auto const& ds = *s.activeDs;
        auto const range = s.range;
        auto& rt = *cclmGpu_;
        st.discardMs = cclmDiscardPending(rt);
        if (cclmPrevUncollected_ && !cclmStats_.empty()) {
            auto const t = cclmLastTiming(rt);
            auto& pv = cclmStats_.back();
            pv.encodeMs = t.encodeMs; pv.queueMs = t.swQueueMs + t.driverMs + t.hwQueueMs; pv.execMs = t.execMs; pv.notifyMs = t.notifyMs;
            pv.gpuKernelMs = t.execMs;
        }
        cclmPrevUncollected_ = false;

        std::vector<int> jobOf(kEff, -1);
        cclmGpuPos_.clear();
        std::size_t nJobs = 0;
        std::size_t const gpuMax = config_.concurrentLmRaceGpuMax > 0 ? static_cast<std::size_t>(config_.concurrentLmRaceGpuMax) : kEff;
        for (std::size_t jj = (kEff > gpuMax ? kEff - gpuMax : 0); jj < kEff; ++jj) {
            auto const& ind = offspring[order[jj]];
            auto cx = ind.x.GetCoefficients();
            auto cy = ind.y.GetCoefficients();
            std::size_t const n = cx.size() + cy.size();
            if (n == 0 || n > kCclmMaxCoeffs) { continue; }
            auto progX = BytecodeProgram::Compile(ind.x, {cx.data(), cx.size()});
            auto progY = BytecodeProgram::Compile(ind.y, {cy.data(), cy.size()});
            if (!progX.valid || !progY.valid || progX.bindings.size() != cx.size() || progY.bindings.size() != cy.size()) { continue; }
            if (cclmJobs_.size() <= nJobs) { cclmJobs_.emplace_back(); }
            auto& job = cclmJobs_[nJobs];
            if (!cclmBuildTape(progX, 0, job.tapeX) || !cclmBuildTape(progY, static_cast<int>(cx.size()), job.tapeY)) { continue; }
            job.nx = static_cast<uint32_t>(cx.size());
            job.coeffs.clear();
            for (auto v : cx) { job.coeffs.push_back(static_cast<float>(v)); }
            for (auto v : cy) { job.coeffs.push_back(static_cast<float>(v)); }
            jobOf[jj] = static_cast<int>(nJobs);
            cclmGpuPos_.push_back(jj);
            ++nJobs;
        }
        std::vector<PairIndividual> shadow, shadowUlp;
        if (config_.concurrentLmParityShadow) {
            shadow.reserve(cclmGpuPos_.size());
            for (auto jj : cclmGpuPos_) { shadow.push_back(offspring[order[jj]]); }
            shadowUlp = shadow;
            for (auto& sh : shadowUlp) {
                auto cx = sh.x.GetCoefficients();
                auto cy = sh.y.GetCoefficients();
                for (auto& v : cx) { v = std::nextafter(v, std::numeric_limits<Operon::Scalar>::infinity()); }
                for (auto& v : cy) { v = std::nextafter(v, std::numeric_limits<Operon::Scalar>::infinity()); }
                sh.x.SetCoefficients(cx);
                sh.y.SetCoefficients(cy);
            }
        }
        bool submitted = false;
        if (nJobs > 0 && cclmAllocEvents(rt) == 0) { cclmReserve(rt, kEff); }
        cclmBeginRace(rt);
        if (nJobs > 0) { submitted = cclmSubmit(rt, cclmFit_, cclmJobs_, nJobs); }
        st.nGpuSubmitted = submitted ? nJobs : 0;
        auto const tSubmit = clk::now();
        st.prepMs = ms(tStage, tSubmit);

        std::vector<char> how(kEff, 0);
        std::vector<std::atomic<char>> claimed(kEff);
        for (auto& c : claimed) { c.store(0, std::memory_order_relaxed); }
        std::atomic<std::size_t> front{0};
        std::atomic<std::size_t> nGpuDesc{0}, nGpuLoss{0}, nGpuIter{0};
        std::atomic<int64_t> readyNs{0};
        bool const streaming = submitted && cclmStreams(rt);
        auto gpuResult = [&](std::size_t q, CclmOut& o) -> bool {
            if (!submitted) { return false; }
            if (streaming) { return cclmPollOne(rt, q, o); }
            if (!cclmTryCollect(rt, cclmOuts_) || q >= cclmOuts_.size()) { return false; }
            o = cclmOuts_[q];
            return true;
        };
        auto applyGpu = [&](std::size_t jj, CclmOut const& o) -> bool {
            int64_t expected = 0;
            readyNs.compare_exchange_strong(expected, clk::now().time_since_epoch().count());
            bool good = o.status == 0;
            for (auto v : o.coeffs) { if (!std::isfinite(v)) { good = false; break; } }
            if (!good) { how[jj] = 2; return false; }
            auto const q = static_cast<std::size_t>(jobOf[jj]);
            std::size_t const nx = cclmJobs_[q].nx;
            auto& ind = offspring[order[jj]];
            ind.x.SetCoefficients(std::vector<Operon::Scalar>(o.coeffs.begin(), o.coeffs.begin() + static_cast<long>(nx)));
            ind.y.SetCoefficients(std::vector<Operon::Scalar>(o.coeffs.begin() + static_cast<long>(nx), o.coeffs.end()));
            nGpuDesc.fetch_add(1, std::memory_order_relaxed);
            nGpuLoss.fetch_add(static_cast<std::size_t>(o.lossCalls), std::memory_order_relaxed);
            nGpuIter.fetch_add(static_cast<std::size_t>(o.iters), std::memory_order_relaxed);
            how[jj] = 1;
            return true;
        };
        auto claimLoop = [&](std::size_t) {
            CclmOut o;
            for (;;) {
                bool did = false;
                for (std::size_t q = cclmGpuPos_.size(); q-- > 0 && submitted;) {
                    std::size_t const jj = cclmGpuPos_[q];
                    if (claimed[jj].load(std::memory_order_relaxed) != 0) { continue; }
                    if (!gpuResult(q, o)) { if (streaming) { continue; } break; }
                    char e = 0;
                    if (!claimed[jj].compare_exchange_strong(e, 1)) { continue; }
                    if (!applyGpu(jj, o)) { optConstsInLoop(offspring[order[jj]], ds, range); }
                    evaluateGuarded(offspring[order[jj]], ds, range);
                    did = true;
                    break;
                }
                if (did) { continue; }
                std::size_t jj = kEff;
                for (;;) {
                    jj = front.fetch_add(1, std::memory_order_relaxed);
                    if (jj >= kEff) { break; }
                    char e = 0;
                    if (claimed[jj].compare_exchange_strong(e, 1)) { break; }
                }
                if (jj >= kEff) { break; }
                bool used = false;
                if (jobOf[jj] >= 0 && gpuResult(static_cast<std::size_t>(jobOf[jj]), o)) { used = applyGpu(jj, o); }
                if (!used) { optConstsInLoop(offspring[order[jj]], ds, range); }
                evaluateGuarded(offspring[order[jj]], ds, range);
            }
        };
        std::size_t const W = std::max<std::size_t>(1, workers);
        parallelApply(W * 4, W, claimLoop);
        auto const tDone = clk::now();
        st.cpuMs = ms(tSubmit, tDone);
        nConstDescents_ += nGpuDesc.load();
        nConstLossCalls_ += nGpuLoss.load();
        nConstIters_ += nGpuIter.load();
        for (std::size_t jj = 0; jj < kEff; ++jj) {
            if (how[jj] == 1) { ++st.nGpu; } else { ++st.nCpu; }
            if (how[jj] == 2) { ++st.nFallback; }
        }
        st.gpuUsed = readyNs.load() != 0;
        cclmPrevUncollected_ = submitted && (streaming || !st.gpuUsed);
        if (st.gpuUsed) { st.gpuReadyMs = ms(tStage, clk::time_point(clk::duration(readyNs.load()))); }
        if (st.gpuUsed && !streaming) {
            auto const t = cclmLastTiming(rt);
            st.encodeMs = t.encodeMs; st.queueMs = t.swQueueMs + t.driverMs + t.hwQueueMs; st.execMs = t.execMs; st.notifyMs = t.notifyMs;
            st.gpuKernelMs = t.execMs;
            int64_t const dn = cclmCollectedDoneNs(rt);
            if (dn > 0) { st.gpuMs = ms(tSubmit, clk::time_point(std::chrono::duration_cast<clk::duration>(std::chrono::nanoseconds(dn)))); }
        }
        if (config_.concurrentLmKeepWarmIters > 0) { cclmPulse(rt, static_cast<uint32_t>(config_.concurrentLmKeepWarmIters), 1); }
        st.stageMs = ms(tStage, clk::now());
        cclmStats_.push_back(st);

        if (config_.concurrentLmParityShadow && !shadow.empty()) {
            std::vector<std::size_t> used;
            for (std::size_t q = 0; q < cclmGpuPos_.size(); ++q) { if (how[cclmGpuPos_[q]] == 1) { used.push_back(q); } }
            parallelApply(used.size(), workers, [&](std::size_t u) {
                std::size_t const q = used[u];
                optConstsInLoop(shadow[q], ds, range);
                evaluateGuarded(shadow[q], ds, range);
                optConstsInLoop(shadowUlp[q], ds, range);
                evaluateGuarded(shadowUlp[q], ds, range);
            });
            ConcurrentLmParityGen pg;
            pg.gen = st.gen;
            pg.mixed.resize(kEff);
            pg.cpuRef.resize(kEff);
            pg.onGpu.assign(kEff, 0);
            pg.fellBack.assign(kEff, 0);
            for (std::size_t jj = 0; jj < kEff; ++jj) {
                pg.mixed[jj] = offspring[order[jj]].shapeErr;
                pg.cpuRef[jj] = offspring[order[jj]].shapeErr;
                pg.onGpu[jj] = how[jj] == 1 ? 1 : 0;
                pg.fellBack[jj] = how[jj] == 2 ? 1 : 0;
            }
            pg.cpuUlp = pg.cpuRef;
            for (auto q : used) {
                pg.cpuRef[cclmGpuPos_[q]] = shadow[q].shapeErr;
                pg.cpuUlp[cclmGpuPos_[q]] = shadowUlp[q].shapeErr;
            }
            cclmParity_.push_back(std::move(pg));
        }
    }

    void optConstsGuarded(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range,
                          bool fullObjective = false) const
    {
        try {
            optConsts(ind, ds, range, fullObjective);
        } catch (std::exception const&) {
            ++nEvalErrors_;
        } catch (...) {
            ++nEvalErrors_;
        }
    }

    void optConstsInLoop(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        if (config_.profileTiming) {
            auto t0 = std::chrono::steady_clock::now();
            optConstsInLoopImpl(ind, ds, range);
            auto t1 = std::chrono::steady_clock::now();
            profLmNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count(), std::memory_order_relaxed);
            return;
        }
        optConstsInLoopImpl(ind, ds, range);
    }

    void optConstsInLoopImpl(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        if (!config_.guardConstOpt) { optConstsGuarded(ind, ds, range); return; }
        auto const before = shapeOnlyGuarded(ind, ds, range);
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        optConstsGuarded(ind, ds, range);
        if (polishWorsened(shapeOnlyGuarded(ind, ds, range), before)) {
            ind.x.SetCoefficients(cx);
            ind.y.SetCoefficients(cy);
        }
    }

    static auto curveFinite(std::vector<Operon::Scalar> const& g) -> bool
    {
        for (auto v : g) { if (!std::isfinite(v)) { return false; } }
        return true;
    }

    auto preConstOptFinite(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range) const -> bool
    {
        thread_local std::vector<Operon::Scalar> gx, gy;
        try {
            interpInto(ind.x, ind.x.GetCoefficients(), ds, range, gx);
            interpInto(ind.y, ind.y.GetCoefficients(), ds, range, gy);
        } catch (std::exception const&) {
            return false;
        } catch (...) {
            return false;
        }
        return curveFinite(gx) && curveFinite(gy);
    }

    void optConstsInLoopGuardedFinite(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        if (config_.rejectNonFiniteBeforeConstOpt) {
            if (!preConstOptFinite(ind, ds, range)) {
                ++nRejectedNonFiniteBeforeConstOpt_;
                evaluateGuarded(ind, ds, range);
                return;
            }
            optConstsInLoop(ind, ds, range);
            evaluateGuarded(ind, ds, range);
            return;
        }
        static bool const measureRecovery = std::getenv("SKETCH_MEASURE_LM_NAN_RECOVERY") != nullptr;
        if (measureRecovery && !preConstOptFinite(ind, ds, range)) {
            ++nRejectedNonFiniteBeforeConstOpt_;
            optConstsInLoop(ind, ds, range);
            evaluateGuarded(ind, ds, range);
            if (std::isfinite(ind.shapeErr)) { ++nLmRecoveredFromNonFiniteStart_; }
            return;
        }
        optConstsInLoop(ind, ds, range);
        evaluateGuarded(ind, ds, range);
    }

    void optConstsInLoopDeferEval(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        if (config_.rejectNonFiniteBeforeConstOpt && !preConstOptFinite(ind, ds, range)) {
            ++nRejectedNonFiniteBeforeConstOpt_;
            return;
        }
        optConstsInLoop(ind, ds, range);
    }

    void evaluateTail(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        ind.complexity = static_cast<Operon::Scalar>(ind.totalLength());
        ind.constBits = constMagnitudeComplexityMode()
            ? (constBitsOf(ind.x) + constBitsOf(ind.y)) : Operon::Scalar{0};
        bool const splitCases = config_.useLexicase && config_.perTreeCredit;
        if (config_.useLexicase && !splitCases) {
            ind.regionErrCache = regionResid(ind, ds, range, sf2(6) ? static_cast<int>(dataX_.size()) : 16);
        }
        if (config_.perTreeCredit) {
            std::vector<Operon::Scalar> regs;
            auto [ex, ey] = perTreeResid(ind, ds, range, splitCases ? &regs : nullptr);
            ind.errX = ex;
            ind.errY = ey;
            if (splitCases) { ind.regionErrCache = std::move(regs); }
            if (config_.curvatureWeight > Operon::Scalar{0}
                && std::isfinite(ind.errX) && std::isfinite(ind.errY)) {
                auto gx = interp(ind.x, ind.x.GetCoefficients(), ds, range);
                auto gy = interp(ind.y, ind.y.GetCoefficients(), ds, range);
                auto const term = curvatureTerm(gx, gy);
                if (config_.curvatureMultiplicative) {
                    ind.errX = applyCurvature(ind.errX, term);
                    ind.errY = applyCurvature(ind.errY, term);
                } else {
                    ind.errX += Operon::Scalar{0.5} * term;
                    ind.errY += Operon::Scalar{0.5} * term;
                }
            }
        }
        if (config_.useResidualRetention) {
            ind.perPointResid_ = residualProfile(ind, ds, range);
        }
        if (config_.useResidualNormalRetention) {
            residualNormalProfile(ind, ds, range, ind.perPointResidSigned_, ind.perPointNormal_);
        }
    }

    auto residualProfile(PairIndividual const& ind, Operon::Dataset const& ds,
                         Operon::Range range) const -> std::vector<Operon::Scalar>
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        auto gx = interp(ind.x, cx, ds, range);
        auto gy = interp(ind.y, cy, ds, range);
        for (std::size_t j = 0; j < gx.size(); ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return {}; }
        }
        auto const N = dataX_.size();
        PointSegmentResult psRes;
        if (oneResOn()) {
            if (!oneResidualOf(gx, gy,false, psRes.residual, nullptr)) { return {}; }
        } else {
        std::vector<Operon::Scalar> nx(gx), ny(gy);
        normalizeInPlace(nx, ny);
        thread_local std::vector<Vec2> rpStrokePts, rpCandPts;
        rpStrokePts.resize(N);
        rpCandPts.resize(nx.size());
        for (std::size_t i = 0; i < N; ++i) { rpStrokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < nx.size(); ++j) { rpCandPts[j] = Vec2{nx[j], ny[j]}; }
        psRes = evaluatePointSegment(rpStrokePts, rpCandPts);
        }
        if (!psRes.ok || psRes.residual.size() != N) { return {}; }
        for (std::size_t i = 0; i < N; ++i) {
            if (!std::isfinite(psRes.residual[i])) { return {}; }
        }
        std::vector<Operon::Scalar> out(N);
        for (std::size_t i = 0; i < N; ++i) {
            out[originalArcIndex(i)] = std::abs(psRes.residual[i]);
        }
        return out;
    }

    void residualNormalProfile(PairIndividual const& ind, Operon::Dataset const& ds,
                                Operon::Range range, std::vector<Operon::Scalar>& outSigned,
                                std::vector<Vec2>& outNormal) const
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        auto gx = interp(ind.x, cx, ds, range);
        auto gy = interp(ind.y, cy, ds, range);
        for (std::size_t j = 0; j < gx.size(); ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return; }
        }
        std::vector<Operon::Scalar> nx(gx), ny(gy);
        normalizeInPlace(nx, ny);
        auto const N = dataX_.size();
        thread_local std::vector<Vec2> rnpStrokePts, rnpCandPts;
        rnpStrokePts.resize(N);
        rnpCandPts.resize(nx.size());
        for (std::size_t i = 0; i < N; ++i) { rnpStrokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < nx.size(); ++j) { rnpCandPts[j] = Vec2{nx[j], ny[j]}; }
        auto const psRes = pointSegOp(rnpStrokePts, rnpCandPts,false,true);
        if (!psRes.ok || psRes.residual.size() != N || psRes.normalAtClosest.size() != N) { return; }
        for (std::size_t i = 0; i < N; ++i) {
            if (!std::isfinite(psRes.residual[i])
                || !std::isfinite(psRes.normalAtClosest[i].x)
                || !std::isfinite(psRes.normalAtClosest[i].y)) {
                return;
            }
        }
        std::vector<Operon::Scalar> signedOut(N);
        std::vector<Vec2> normalOut(N);
        for (std::size_t i = 0; i < N; ++i) {
            auto const idx = originalArcIndex(i);
            signedOut[idx] = psRes.residual[i];
            normalOut[idx] = psRes.normalAtClosest[i];
        }
        outSigned = std::move(signedOut);
        outNormal = std::move(normalOut);
    }

    void evaluate(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range) const
    {
        if (config_.floorOnChamfer) {
            ind.shapeErr = shapeOnlyParts(ind, ds, range, &ind.chamferErr);
        } else {
            ind.shapeErr = shapeOnlyParts(ind, ds, range, nullptr);
        }
        evaluateTail(ind, ds, range);
    }

    auto liveCoarseBound(PairIndividual const& ind, Operon::Dataset const& ds,
                         Operon::Range range) const -> double
    {
        auto gx = interp(ind.x, ind.x.GetCoefficients(), ds, range);
        auto gy = interp(ind.y, ind.y.GetCoefficients(), ds, range);
        if (gx.size() < 2U || gx.size() != gy.size() || dataX_.size() < 2U) {
            return std::numeric_limits<double>::max();
        }
        for (std::size_t i = 0; i < gx.size(); ++i) {
            if (!std::isfinite(gx[i]) || !std::isfinite(gy[i])) {
                return std::numeric_limits<double>::max();
            }
        }
        normalizeInPlace(gx, gy);
        std::vector<Vec2> reference(dataX_.size());
        std::vector<Vec2> candidate(gx.size());
        for (std::size_t i = 0; i < reference.size(); ++i) {
            reference[i] = {dataX_[i], dataY_[i]};
        }
        for (std::size_t i = 0; i < candidate.size(); ++i) {
            candidate[i] = {gx[i], gy[i]};
        }
        return coarseFourPointForwardLowerBound(reference, candidate);
    }

    void evaluateFromGrid(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range,
                          std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                          CisPointSegBatch const* psBatch = nullptr, std::size_t psBatchIdx = 0,
                          CisPointSegBatch const* psBatch48 = nullptr) const
    {
        if (config_.floorOnChamfer) {
            ind.shapeErr = shapeOnlyPartsFromGrid(gx, gy, &ind.chamferErr, psBatch, psBatchIdx, psBatch48);
        } else {
            ind.shapeErr = shapeOnlyPartsFromGrid(gx, gy, nullptr, psBatch, psBatchIdx, psBatch48);
        }
        evaluateTail(ind, ds, range);
    }

    auto perTreeResid(PairIndividual const& ind, Operon::Dataset const& ds,
                      Operon::Range range,
                      std::vector<Operon::Scalar>* regionsOut = nullptr,
                      int R = 16) const
        -> std::pair<Operon::Scalar, Operon::Scalar>
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto xc = ind.x.GetCoefficients();
        auto yc = ind.y.GetCoefficients();
        auto gx = interp(ind.x, xc, ds, range);
        auto gy = interp(ind.y, yc, ds, range);
        auto const M = gx.size();
        if (M < 2 || dataX_.size() < 2) { return {BIG, BIG}; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return {BIG, BIG}; }
        }
        std::vector<Operon::Scalar> nx(gx), ny(gy);
        normalizeInPlace(nx, ny);

        constexpr int n = 64;
        std::vector<Operon::Scalar> px, py, bx, by;
        arcResample(dataX_, dataY_, n, px, py);
        arcResample(nx, ny, n, bx, by);

        auto const K = noise_.size();
        bool const whiten = K > 0;
        std::vector<Operon::Scalar> invS2;
        if (whiten) {
            invS2.resize(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                auto const k = std::min(static_cast<std::size_t>(i) * K
                                            / static_cast<std::size_t>(n), K - 1);
                auto const s2 = noise_.s2[k];
                invS2[static_cast<std::size_t>(i)] =
                    (std::isfinite(s2) && s2 > Operon::Scalar{0})
                        ? Operon::Scalar{1} / s2 : Operon::Scalar{1};
            }
        }

        std::vector<Operon::Scalar> sx(static_cast<std::size_t>(n));
        std::vector<Operon::Scalar> sy(static_cast<std::size_t>(n));
        std::vector<int> matchIdx;
        Operon::Scalar bestCost = BIG;
        Operon::Scalar outX = BIG, outY = BIG;
        std::vector<Operon::Scalar> regAcc;
        std::vector<int> regCnt;
        for (int dir = 0; dir < 2; ++dir) {
            std::vector<Operon::Scalar> rx(bx), ry(by);
            if (dir == 1) { std::reverse(rx.begin(), rx.end()); std::reverse(ry.begin(), ry.end()); }
            int const lo = dataClosed_ ? bestShift(px, py, rx, ry) - 1 : 0;
            int const hi = dataClosed_ ? lo + 2 : 0;
            for (int off = lo; off <= hi; ++off) {
                if (dataClosed_) {
                    int const r = ((off % n) + n) % n;
                    for (int k = 0; k < n; ++k) {
                        sx[static_cast<std::size_t>(k)] = rx[static_cast<std::size_t>((k + r) % n)];
                        sy[static_cast<std::size_t>(k)] = ry[static_cast<std::size_t>((k + r) % n)];
                    }
                } else {
                    sx = rx; sy = ry;
                }
                auto cost = dtwPath(px, py, sx, sy, matchIdx);
                if (cost < bestCost) {
                    bestCost = cost;
                    Operon::Scalar ax = 0, ay = 0;
                    if (regionsOut != nullptr) {
                        regAcc.assign(static_cast<std::size_t>(2 * R), Operon::Scalar{0});
                        regCnt.assign(static_cast<std::size_t>(R), 0);
                    }
                    for (int i = 0; i < n; ++i) {
                        auto const j = static_cast<std::size_t>(matchIdx[static_cast<std::size_t>(i)]);
                        auto ex = sx[j] - px[static_cast<std::size_t>(i)];
                        auto ey = sy[j] - py[static_cast<std::size_t>(i)];
                        if (whiten) {
                            auto const iv = invS2[static_cast<std::size_t>(i)];
                            auto const ex2 = ex * ex * iv;
                            auto const ey2 = ey * ey * iv;
                            ax += ex2;
                            ay += ey2;
                            if (regionsOut != nullptr) {
                                auto const k = static_cast<std::size_t>(std::min(i * R / n, R - 1));
                                regAcc[k] += ex2;
                                regAcc[static_cast<std::size_t>(R) + k] += ey2;
                                ++regCnt[k];
                            }
                        } else {
                            ax += ex * ex;
                            ay += ey * ey;
                            if (regionsOut != nullptr) {
                                auto const k = static_cast<std::size_t>(std::min(i * R / n, R - 1));
                                regAcc[k] += ex * ex;
                                regAcc[static_cast<std::size_t>(R) + k] += ey * ey;
                                ++regCnt[k];
                            }
                        }
                    }
                    outX = ax / static_cast<Operon::Scalar>(n);
                    outY = ay / static_cast<Operon::Scalar>(n);
                    if (regionsOut != nullptr) {
                        for (std::size_t k = 0; k < static_cast<std::size_t>(R); ++k) {
                            if (regCnt[k] == 0) { continue; }
                            auto const d = static_cast<Operon::Scalar>(regCnt[k]);
                            regAcc[k] /= d;
                            regAcc[static_cast<std::size_t>(R) + k] /= d;
                        }
                        *regionsOut = regAcc;
                    }
                }
            }
        }
        return {outX, outY};
    }

    auto sampledCurve(PairIndividual const& ind, Operon::Dataset const& ds,
                      Operon::Range range) const
        -> std::pair<std::vector<Operon::Scalar>, std::vector<Operon::Scalar>>
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        auto gx = interp(ind.x, cx, ds, range);
        auto gy = interp(ind.y, cy, ds, range);
        for (std::size_t j = 0; j < gx.size(); ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return {{}, {}}; }
        }
        std::vector<Operon::Scalar> nx(gx), ny(gy);
        normalizeInPlace(nx, ny);
        return {nx, ny};
    }

    auto regionResid(PairIndividual const& ind, Operon::Dataset const& ds,
                     Operon::Range range, int R = 16) const -> std::vector<Operon::Scalar>
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        auto gx = interp(ind.x, cx, ds, range);
        auto gy = interp(ind.y, cy, ds, range);
        for (std::size_t j = 0; j < gx.size(); ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return {}; }
        }
        std::vector<Operon::Scalar> nx(gx), ny(gy);
        normalizeInPlace(nx, ny);
        auto const N = dataX_.size();
        std::vector<Operon::Scalar> reg(static_cast<std::size_t>(R), Operon::Scalar{0});
        std::vector<int> cnt(static_cast<std::size_t>(R), 0);
        bool const nwLex = config_.useNoiseWeightedLexicase;
        std::vector<Operon::Scalar> const* extSig = nullptr;
        std::vector<std::size_t> const* sigLut = nullptr;
        Operon::Scalar constSig2 = Operon::Scalar{1};
        if (nwLex) {
            if (config_.useExternalPointSigma2 && !externalPointSigma2_.empty()) {
                extSig = &externalPointSigma2_;
                sigLut = &noiseRegionLUT(N, externalPointSigma2_.size());
            } else if (noise_.size() > 0) {
                extSig = &noise_.s2;
                sigLut = &noiseRegionLUT(N, noise_.size());
            } else {
                constSig2 = noiseMeanS2_;
            }
        }
        bool const cohLex = config_.useCoherentLexicaseResidual || (oneResOn() && oneResBit(8));
        std::vector<Operon::Scalar> cohRes;
        if (cohLex && oneResOn()) {
            std::vector<Operon::Scalar> oRes;
            if (oneResidualOf(gx, gy,false, oRes, nullptr) && oRes.size() == N) { cohRes = std::move(oRes); }
        } else if (cohLex) {
            thread_local std::vector<Vec2> g6cStrokePts, g6cCandPts;
            g6cStrokePts.resize(N);
            g6cCandPts.resize(nx.size());
            for (std::size_t i = 0; i < N; ++i) { g6cStrokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
            for (std::size_t j = 0; j < nx.size(); ++j) { g6cCandPts[j] = Vec2{nx[j], ny[j]}; }
            auto const g6cPsRes = evaluatePointSegment(g6cStrokePts, g6cCandPts);
            if (g6cPsRes.ok && g6cPsRes.residual.size() == N) {
                cohRes = g6cPsRes.residual;
            }
        }
        for (std::size_t i = 0; i < N; ++i) {
            auto r2 = (cohLex && cohRes.size() == N)
                ? (cohRes[i] * cohRes[i])
                : pointToPolyline2(dataX_[i], dataY_[i], nx, ny);
            if (nwLex) {
                Operon::Scalar const s2 = (extSig != nullptr && sigLut != nullptr)
                    ? (*extSig)[(*sigLut)[originalArcIndex(i)]] : constSig2;
                if (std::isfinite(s2) && s2 > Operon::Scalar{0}) { r2 /= s2; }
            }
            auto k = std::min<std::size_t>(i * static_cast<std::size_t>(R) / N,
                                           static_cast<std::size_t>(R) - 1);
            reg[k] += r2;
            ++cnt[k];
        }
        for (std::size_t k = 0; k < reg.size(); ++k) {
            if (cnt[k] > 0) { reg[k] /= static_cast<Operon::Scalar>(cnt[k]); }
        }
        return reg;
    }

    auto dtwPath(std::vector<Operon::Scalar> const& px, std::vector<Operon::Scalar> const& py,
                 std::vector<Operon::Scalar> const& qx, std::vector<Operon::Scalar> const& qy,
                 std::vector<int>& matchIdx) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max() / 4;
        int const n = static_cast<int>(px.size());
        int const m = static_cast<int>(qx.size());
        thread_local std::vector<Operon::Scalar> D;
        D.assign(static_cast<std::size_t>((n + 1)) * static_cast<std::size_t>(m + 1), BIG);
        auto const stride = static_cast<std::size_t>(m + 1);
        auto at = [stride](int i, int j) -> std::size_t {
            return static_cast<std::size_t>(i) * stride + static_cast<std::size_t>(j);
        };
        D[at(0, 0)] = 0;
        for (int i = 1; i <= n; ++i) {
            auto pix = px[static_cast<std::size_t>(i - 1)], piy = py[static_cast<std::size_t>(i - 1)];
            for (int j = 1; j <= m; ++j) {
                auto dx = qx[static_cast<std::size_t>(j - 1)] - pix, dy = qy[static_cast<std::size_t>(j - 1)] - piy;
                auto c = dx * dx + dy * dy;
                auto best = D[at(i - 1, j - 1)];
                best = std::min(best, D[at(i - 1, j)]);
                best = std::min(best, D[at(i, j - 1)]);
                D[at(i, j)] = c + best;
            }
        }
        matchIdx.assign(static_cast<std::size_t>(n), 0);
        int i = n, j = m;
        while (i > 0 && j > 0) {
            matchIdx[static_cast<std::size_t>(i - 1)] = j - 1;
            auto diag = D[at(i - 1, j - 1)];
            auto up   = D[at(i - 1, j)];
            auto left = D[at(i, j - 1)];
            if (diag <= up && diag <= left) { --i; --j; }
            else if (up <= left)            { --i; }
            else                            { --j; }
        }
        while (i > 0) { matchIdx[static_cast<std::size_t>(i - 1)] = 0; --i; }
        return D[at(n, m)];
    }

    auto footPointT(Operon::Tree const& xTree, Operon::Tree const& yTree,
                    Operon::Dataset const& gridDs, Operon::Range gridRange) const
        -> std::vector<Operon::Scalar>
    {
        auto const N = dataX_.size();
        auto cX = xTree.GetCoefficients();
        auto cY = yTree.GetCoefficients();
        auto xVals = interp(xTree, cX, gridDs, gridRange);
        auto yVals = interp(yTree, cY, gridDs, gridRange);
        auto const M = xVals.size();
        if (M < 2 || N < 2) { return {}; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(xVals[j]) || !std::isfinite(yVals[j])) { return {}; }
        }

        std::size_t const NS_CAP = 96;
        std::size_t const Ns = std::min(M, NS_CAP);
        std::vector<std::size_t> sampOrig(Ns);
        std::vector<Operon::Scalar> baseCx(Ns), baseCy(Ns);
        for (std::size_t s = 0; s < Ns; ++s) {
            std::size_t gi = (Ns == 1) ? 0 : (s * (M - 1)) / (Ns - 1);
            sampOrig[s] = gi;
            baseCx[s] = xVals[gi];
            baseCy[s] = yVals[gi];
        }

        std::vector<std::size_t> rots;
        if (dataClosed_) {
            std::size_t const R = 16;
            std::size_t const step = std::max<std::size_t>(1, Ns / R);
            for (std::size_t r = 0; r < Ns && rots.size() < R; r += step) { rots.push_back(r); }
        } else {
            rots.push_back(0);
        }

        std::vector<Operon::Scalar> qx(Ns), qy(Ns);
        std::vector<std::size_t> qOrig(Ns);
        std::vector<int> matchIdx;
        Operon::Scalar bestCost = std::numeric_limits<Operon::Scalar>::max();
        std::vector<Operon::Scalar> bestT(N, Operon::Scalar{0});
        bool found = false;
        for (int dir = 0; dir < 2; ++dir) {
            for (auto r : rots) {
                for (std::size_t k = 0; k < Ns; ++k) {
                    std::size_t idx = (k + r) % Ns;
                    std::size_t src = (dir == 0) ? idx : (Ns - 1 - idx);
                    qx[k] = baseCx[src]; qy[k] = baseCy[src]; qOrig[k] = sampOrig[src];
                }
                auto cost = dtwPath(dataX_, dataY_, qx, qy, matchIdx);
                if (cost < bestCost) {
                    bestCost = cost;
                    found = true;
                    for (std::size_t d = 0; d < N; ++d) {
                        std::size_t gi = qOrig[static_cast<std::size_t>(matchIdx[d])];
                        bestT[d] = static_cast<Operon::Scalar>(gi) / static_cast<Operon::Scalar>(M - 1);
                    }
                }
            }
        }
        if (!found) { return {}; }
        return bestT;
    }

    void optConstsLm(PairIndividual& ind, Operon::Dataset const& gridDs, Operon::Range gridRange) const
    {
        auto const N = dataX_.size();
        if (N < 2) { return; }
        if (ind.x.GetCoefficients().empty() && ind.y.GetCoefficients().empty()) { return; }

        std::vector<Operon::Scalar> flatX(N * 2, Operon::Scalar{0});
        std::vector<Operon::Scalar> flatY(N * 2, Operon::Scalar{0});
        for (std::size_t i = 0; i < N; ++i) {
            flatX[N + i] = dataX_[i];
            flatY[N + i] = dataY_[i];
        }
        Operon::Dataset dsX(gsl::not_null{flatX.data()}, static_cast<int>(N), 2);
        Operon::Dataset dsY(gsl::not_null{flatY.data()}, static_cast<int>(N), 2);

        Operon::Problem probX(&dsX);
        Operon::Problem probY(&dsY);
        auto varsX = dsX.GetVariables();
        auto varsY = dsY.GetVariables();
        probX.SetTrainingRange(0, static_cast<int>(N));
        probX.SetTarget(varsX[1].Hash);
        probX.SetInputs(std::vector{varsX[0].Hash});
        probY.SetTrainingRange(0, static_cast<int>(N));
        probY.SetTarget(varsY[1].Hash);
        probY.SetInputs(std::vector{varsY[0].Hash});

        using LMOpt = Operon::LevenbergMarquardtOptimizer<DTable>;
        LMOpt lmX(dtable_, &probX);
        LMOpt lmY(dtable_, &probY);
        lmX.SetIterations(10);
        lmY.SetIterations(10);
        Operon::RandomGenerator lmRng(config_.seed);

        auto const tVarX = varsX[0].Hash;
        auto const tVarY = varsY[0].Hash;

        auto const lmBefore = shapeOnlyGuarded(ind, gridDs, gridRange);
        auto const cx0 = ind.x.GetCoefficients();
        auto const cy0 = ind.y.GetCoefficients();
        constexpr int OUTER = 3;
        for (int it = 0; it < OUTER; ++it) {
            auto footT = footPointT(ind.x, ind.y, gridDs, gridRange);
            if (footT.empty()) { break; }
            dsX.SetValues(tVarX, Operon::Range{0, N}, {footT.data(), footT.size()});
            dsY.SetValues(tVarY, Operon::Range{0, N}, {footT.data(), footT.size()});
            auto ox = lmX.Optimize(lmRng, ind.x);
            if (ox.has_value() && allConstantsFiniteVec(ox->FinalParameters)) { ind.x.SetCoefficients(ox->FinalParameters); }
            else if (ox.has_value()) { ++nNonFiniteLmConst_; nSafetyCounters().nonFiniteLmConst.fetch_add(1, std::memory_order_relaxed); }
            auto oy = lmY.Optimize(lmRng, ind.y);
            if (oy.has_value() && allConstantsFiniteVec(oy->FinalParameters)) { ind.y.SetCoefficients(oy->FinalParameters); }
            else if (oy.has_value()) { ++nNonFiniteLmConst_; nSafetyCounters().nonFiniteLmConst.fetch_add(1, std::memory_order_relaxed); }
        }
        if (polishWorsened(shapeOnlyGuarded(ind, gridDs, gridRange), lmBefore)) {
            ind.x.SetCoefficients(cx0);
            ind.y.SetCoefficients(cy0);
        }
    }

    static auto frequencyParamMask(Operon::Tree const& tree) -> std::vector<char>
    {
        auto const& nodes = tree.Nodes();
        std::vector<char> isFreqNode(nodes.size(), 0);
        if (nodes.empty()) { return {}; }

        std::vector<char> hasVar(nodes.size(), 0);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            char h = nodes[i].IsVariable() ? char{1} : char{0};
            if (nodes[i].Arity > 0) {
                for (auto c : tree.Indices(i)) { h = h || hasVar[c] ? char{1} : char{0}; }
            }
            hasVar[i] = h;
        }

        std::vector<std::pair<std::size_t, char>> stack;
        stack.emplace_back(nodes.size() - 1, char{0});
        while (!stack.empty()) {
            auto [i, inTrig] = stack.back();
            stack.pop_back();
            auto const& n = nodes[i];
            if (inTrig != 0) {
                if (n.IsVariable()) {
                    isFreqNode[i] = 1;
                } else if (n.IsMultiplication() && n.Arity > 0) {
                    for (auto a : tree.Indices(i)) {
                        if (!nodes[a].IsConstant()) { continue; }
                        for (auto b : tree.Indices(i)) {
                            if (b != a && hasVar[b] != 0) { isFreqNode[a] = 1; break; }
                        }
                    }
                }
            }
            char const childInTrig = (inTrig != 0 || n.IsSin() || n.IsCos()) ? char{1} : char{0};
            if (n.Arity > 0) {
                for (auto c : tree.Indices(i)) { stack.emplace_back(c, childInTrig); }
            }
        }

        std::vector<char> mask;
        mask.reserve(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].Optimize) { mask.push_back(isFreqNode[i]); }
        }
        return mask;
    }

    static auto strokeTurningNumber(std::vector<Operon::Scalar> const& X,
                                    std::vector<Operon::Scalar> const& Y, bool closed) -> double
    {
        constexpr double TWO_PI = 6.283185307179586;
        constexpr double PI = 3.141592653589793;
        std::size_t const M = std::min(X.size(), Y.size());
        if (M < 5) { return 0.0; }

        std::size_t const h = std::min<std::size_t>(std::max<std::size_t>(M / 32, 1), (M - 1) / 2);
        std::vector<double> sx(M), sy(M);
        if (h == 0) {
            for (std::size_t i = 0; i < M; ++i) { sx[i] = static_cast<double>(X[i]); sy[i] = static_cast<double>(Y[i]); }
        } else {
            for (std::size_t i = 0; i < M; ++i) {
                double ax = 0.0, ay = 0.0;
                for (long j = static_cast<long>(i) - static_cast<long>(h);
                     j <= static_cast<long>(i) + static_cast<long>(h); ++j) {
                    std::size_t jj;
                    if (closed) {
                        long m = static_cast<long>(M);
                        jj = static_cast<std::size_t>(((j % m) + m) % m);
                    } else {
                        jj = (j < 0) ? 0 : (j >= static_cast<long>(M) ? M - 1 : static_cast<std::size_t>(j));
                    }
                    ax += static_cast<double>(X[jj]);
                    ay += static_cast<double>(Y[jj]);
                }
                double const c = static_cast<double>(2 * h + 1);
                sx[i] = ax / c; sy[i] = ay / c;
            }
        }

        double total = 0.0, prev = 0.0;
        bool have = false;
        for (std::size_t i = 1; i < M; ++i) {
            double dx = sx[i] - sx[i - 1];
            double dy = sy[i] - sy[i - 1];
            if (dx == 0.0 && dy == 0.0) { continue; }
            double ang = std::atan2(dy, dx);
            if (have) {
                double d = ang - prev;
                while (d > PI) { d -= TWO_PI; }
                while (d < -PI) { d += TWO_PI; }
                total += std::abs(d);
            }
            prev = ang;
            have = true;
        }
        return total / TWO_PI;
    }

