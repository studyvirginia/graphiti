
    template <typename> friend struct AmplitudeProviderTestHook;
    template <typename> friend struct DtwPathTestHook;
    template <typename> friend struct ConstOptTestHook;
    template <typename> friend struct QuotientFiniteMinTestHook;
    template <typename> friend struct GradMutTestHook;

    MultitreeGP(DTable const& dtable,
                std::vector<Operon::Scalar> dataX,
                std::vector<Operon::Scalar> dataY,
                MultitreeConfig config,
                NoiseModel noise = {})
        : dtable_{&dtable}
        , dataX_{std::move(dataX)}
        , dataY_{std::move(dataY)}
        , config_{config}
        , noise_{std::move(noise)}
    {
        if (dataX_.size() < 2 || dataX_.size() != dataY_.size()) {
            inputRejected_ = true;
            return;
        }
        for (std::size_t i = 0; i < dataX_.size(); ++i) {
            if (!std::isfinite(dataX_[i]) || !std::isfinite(dataY_[i])) {
                inputRejected_ = true;
                return;
            }
        }
        auto precomputeStart = std::chrono::steady_clock::now();
        ArcLengthGuardRatioGuard arcGuardBuildGuard(config_.arcLengthGuardRatio);
        normalizeInPlace(dataX_, dataY_);
        sf2ApplyBudgetRules();
        config_.gradientMutationInRunErr = config_.gradientMutationMode > 0;
        if (config_.useRotationGauge) {
            applyRotationGauge(dataX_, dataY_, rotationGaugeAngle_);
        }
        if (config_.useParametrizationGauge || config_.useQuotientProjection) {
            int const n = static_cast<int>(dataX_.size());
            if (n >= 2) {
                std::vector<float> px(static_cast<std::size_t>(n));
                std::vector<float> py(static_cast<std::size_t>(n));
                for (int i = 0; i < n; ++i) {
                    px[static_cast<std::size_t>(i)] = static_cast<float>(dataX_[static_cast<std::size_t>(i)]);
                    py[static_cast<std::size_t>(i)] = static_cast<float>(dataY_[static_cast<std::size_t>(i)]);
                }
                parametrizationPerm_.resize(static_cast<std::size_t>(n));
                parametrizationGaugeRecord_ = canonicalizeParametrization(
                    px.data(), py.data(), n, parametrizationPerm_.data());
                std::vector<Operon::Scalar> gx(static_cast<std::size_t>(n));
                std::vector<Operon::Scalar> gy(static_cast<std::size_t>(n));
                for (int k = 0; k < n; ++k) {
                    int const src = parametrizationPerm_[static_cast<std::size_t>(k)];
                    gx[static_cast<std::size_t>(k)] = dataX_[static_cast<std::size_t>(src)];
                    gy[static_cast<std::size_t>(k)] = dataY_[static_cast<std::size_t>(src)];
                }
                dataX_.swap(gx);
                dataY_.swap(gy);
            }
        }
        buildBlocks(dataX_, dataY_,false, dataBlocks_);
        arcResample(dataX_, dataY_, 64, dataArcX_, dataArcY_);
        if (cisGrid48Bit1() && !cisGrid48Bit0()) {
            cisGrid48Counters().bit1NoBit0.fetch_add(1, std::memory_order_relaxed);
        }
        if (cisGrid48Bit0()) {
            std::array<double, 2> shp48Mean{};
            double shp48Scale = 1.0;
            auto const shp48Raw = sketch::sigmaV2ShapeSpaceRaw(dataX_.data(), dataY_.data(),
                                                                static_cast<int>(dataX_.size()), 48,
                                                                shp48Mean, shp48Scale);
            std::vector<std::array<double, 2>> shp48(shp48Raw.size());
            for (std::size_t i = 0; i < shp48Raw.size(); ++i) {
                shp48[i] = {(shp48Raw[i][0] - shp48Mean[0]) / shp48Scale,
                            (shp48Raw[i][1] - shp48Mean[1]) / shp48Scale};
            }
            cisGrid48Stroke_.resize(shp48.size());
            for (std::size_t i = 0; i < shp48.size(); ++i) {
                cisGrid48Stroke_[i] = Vec2{static_cast<Operon::Scalar>(shp48[i][0]),
                                           static_cast<Operon::Scalar>(shp48[i][1])};
            }
            cisGrid48Ok_ = cisGrid48Stroke_.size() == 48;
            if (cisGrid48Bit2()) {
                cisGrid48StrokeData_.resize(shp48Raw.size());
                for (std::size_t i = 0; i < shp48Raw.size(); ++i) {
                    cisGrid48StrokeData_[i] = Vec2{static_cast<Operon::Scalar>(shp48Raw[i][0]),
                                                   static_cast<Operon::Scalar>(shp48Raw[i][1])};
                }
                cisGrid48Scale_ = shp48Scale;
            }
            if (!cisGrid48Ok_) {
                cisGrid48Counters().badStroke.fetch_add(1, std::memory_order_relaxed);
            }
            if (cisGrid48Ok_ && cisGrid48Bit1()) {
                sketch::SigmaV2Model const* model = sketch::loadSigmaV2ModelCached();
                if (model == nullptr) {
                    cisGrid48Counters().missingModel.fetch_add(1, std::memory_order_relaxed);
                } else {
                    auto const dec48 = haarForward(std::vector<double>(48, 0.0));
                    auto const expectCells = buildCisCells(dec48, std::vector<double>(1, 1.0), {}, true);
                    if (expectCells.size() != model->NC) {
                        cisGrid48Counters().cellMismatch.fetch_add(1, std::memory_order_relaxed);
                    } else {
                        static bool const zi5DiagOn = std::getenv("SKETCH_B3_ZI5_DIAG") != nullptr;
                        std::chrono::steady_clock::time_point zi5T0;
                        if (zi5DiagOn) { zi5T0 = std::chrono::steady_clock::now(); }
                        auto const zi5 = sketch::sigmaV2Infer(*model, dataX_.data(), dataY_.data(),
                                                              static_cast<int>(dataX_.size()));
                        if (zi5DiagOn) {
                            auto const ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - zi5T0).count();
                            std::fprintf(stderr, "[B3ZI5DIAG] sigmaV2Infer_ns=%lld\n", static_cast<long long>(ns));
                        }
                        cisGrid48Zi5LnE_ = zi5.lnE;
                        cisGrid48Zi5Sd_ = zi5.sd;
                        cisGrid48Zi5MemberLnE_ = zi5.memberLnE;
                        cisGrid48Zi5Ok_ = zi5.lnE.size() == expectCells.size();
                        if (!cisGrid48Zi5Ok_) {
                            cisGrid48Counters().cellMismatch.fetch_add(1, std::memory_order_relaxed);
                        }
                    }
                }
            }
        }
        if (config_.useGridSdf) {
            buildGridSdf(dataX_, dataY_, noise_, config_.gridSdfResolution, gridSdf_);
        }
#if defined(__APPLE__)
        if (config_.useAmxChamfer) {
            buildAmxData(dataX_, dataY_, amxData_);
        }
#endif
#if defined(__APPLE__)
        if (config_.useAccelerateTranscendentals) {
            accelDTable_ = *dtable_;
            RegisterAccelerateTranscendentals(*accelDTable_);
            dtable_ = &(*accelDTable_);
        }
#endif
        if (!noise_.empty()) {
            Operon::Scalar acc = 0;
            for (auto v : noise_.s2) { acc += v; }
            auto const m = acc / static_cast<Operon::Scalar>(noise_.size());
            if (std::isfinite(m) && m > Operon::Scalar{0}) { noiseMeanS2_ = m; }
        }
        if (dataX_.size() >= 2) {
            auto dx = dataX_.front() - dataX_.back();
            auto dy = dataY_.front() - dataY_.back();
            dataClosed_ = std::sqrt(dx*dx + dy*dy) < config_.closedGap;
            if (sf2(9)) {
                double const s0 = noise_.empty() ? sf2NoiseVar() : static_cast<double>(noise_.s2.front());
                double const sN = noise_.empty() ? sf2NoiseVar() : static_cast<double>(noise_.s2.back());
                double const g2 = static_cast<double>(dx) * dx + static_cast<double>(dy) * dy;
                dataClosed_ = g2 < (s0 + sN) * (-2.0 * std::log(0.05));
            }
        }
        if (config_.useQuotientFiniteMin && config_.useQuotientProjection) {
            buildQuotientFiniteMinBundles();
        }
        turningNumber_ = strokeTurningNumber(dataX_, dataY_, dataClosed_);
        dataTurningF_ = turningFourier(dataX_, dataY_, dataClosed_, config_.curvatureHarmonics);
        curvatureFloor_ = static_cast<double>(config_.curvatureFloorTurns)
                        + static_cast<double>(config_.curvatureFloorRel) * dataTurningF_;
        if (config_.lmDownsampleK > 0 && dataX_.size() > config_.lmDownsampleK) {
            auto const strokePts = [&] {
                std::vector<Vec2> pts(dataX_.size());
                for (std::size_t i = 0; i < dataX_.size(); ++i) { pts[i] = Vec2{dataX_[i], dataY_[i]}; }
                return pts;
            }();
            auto const sel = selectDownsampleIndices(strokePts, config_.lmDownsampleK);
            if (sel.ok && !sel.usedAllPoints && sel.indices.size() >= 4) {
                dsOrigIdx_ = sel.indices;
                dataXDs_.resize(dsOrigIdx_.size());
                dataYDs_.resize(dsOrigIdx_.size());
                for (std::size_t i = 0; i < dsOrigIdx_.size(); ++i) {
                    dataXDs_[i] = dataX_[dsOrigIdx_[i]];
                    dataYDs_[i] = dataY_[dsOrigIdx_[i]];
                }
                buildBlocks(dataXDs_, dataYDs_,false, dataBlocksDs_);
            }
        }
        if (sf2(3)) {
            auto const n = dataX_.size();
            std::size_t lowN = std::max<std::size_t>(32, sf2NyquistCount());
            if (lowN >= n) {
                config_.annealLowN = n;
                config_.annealLowM = config_.gridSize;
            } else {
                config_.annealLowN = lowN;
                config_.annealLowM = static_cast<std::size_t>(std::lround(
                    static_cast<double>(config_.gridSize) * static_cast<double>(lowN) / static_cast<double>(n)));
            }
            if (std::getenv("SKETCH_SF2_DIAG") != nullptr) {
                std::fprintf(stderr, "[SF2DIAG] n=%zu sigma2=%.6g kmax=%zu lowN=%zu lowM=%zu closed=%d\n", n, sf2NoiseVar(),
                             sf2KMax_, config_.annealLowN, config_.annealLowM, dataClosed_ ? 1 : 0);
            }
        }
        if (config_.annealResolution > 0.0 && config_.annealLowN >= 2
            && dataX_.size() > config_.annealLowN) {
            fullX_ = dataX_;
            fullY_ = dataY_;
            fullBlocks_ = dataBlocks_;
            if (config_.useGridSdf) {
                fullGridSdf_ = gridSdf_;
            }
#if defined(__APPLE__)
            if (config_.useAmxChamfer) {
                fullAmxData_ = amxData_;
            }
#endif
            if (config_.annealStage2Resolution > config_.annealResolution
                && config_.annealMidN >= 2
                && config_.annealMidN > config_.annealLowN
                && fullX_.size() > config_.annealMidN) {
                arcResample(fullX_, fullY_, static_cast<int>(config_.annealMidN), midX_, midY_);
                buildBlocks(midX_, midY_,false, midBlocks_);
#if defined(__APPLE__)
                if (config_.useAmxChamfer) {
                    buildAmxData(midX_, midY_, midAmxData_);
                }
#endif
                annealMidActive_ = true;
            }
            std::vector<Operon::Scalar> lx, ly;
            arcResample(dataX_, dataY_, static_cast<int>(config_.annealLowN), lx, ly);
            dataX_ = std::move(lx);
            dataY_ = std::move(ly);
            buildBlocks(dataX_, dataY_,false, dataBlocks_);
            if (config_.useGridSdf) {
                buildGridSdf(dataX_, dataY_, noise_, config_.gridSdfResolution, gridSdf_);
            }
#if defined(__APPLE__)
            if (config_.useAmxChamfer) {
                buildAmxData(dataX_, dataY_, amxData_);
            }
#endif
            annealActive_ = true;
        }
        auto precomputeEnd = std::chrono::steady_clock::now();
        precomputeTimeMs_ = std::chrono::duration<double, std::milli>(precomputeEnd - precomputeStart).count();
    }

    auto workerCount() const -> std::size_t
    {
        if (config_.nThreads != 0) { return std::max<std::size_t>(1, config_.nThreads); }
        auto const hc = std::thread::hardware_concurrency();
        return hc == 0 ? std::size_t{1} : static_cast<std::size_t>(hc);
    }

    struct CreditModeGuard {
        bool prev;
        explicit CreditModeGuard(bool v) : prev(perTreeCreditMode()) { perTreeCreditMode() = v; }
        CreditModeGuard(CreditModeGuard const&) = delete;
        CreditModeGuard(CreditModeGuard&&) = delete;
        auto operator=(CreditModeGuard const&) -> CreditModeGuard& = delete;
        auto operator=(CreditModeGuard&&) -> CreditModeGuard& = delete;
        ~CreditModeGuard() { perTreeCreditMode() = prev; }
    };

    struct FloorClampGuard {
        Operon::Scalar prev;
        explicit FloorClampGuard(Operon::Scalar v) : prev(errorClampValue()) { errorClampValue() = v; }
        FloorClampGuard(FloorClampGuard const&) = delete;
        FloorClampGuard(FloorClampGuard&&) = delete;
        auto operator=(FloorClampGuard const&) -> FloorClampGuard& = delete;
        auto operator=(FloorClampGuard&&) -> FloorClampGuard& = delete;
        ~FloorClampGuard() { errorClampValue() = prev; }
    };
    struct FloorModeGuard {
        bool prev;
        explicit FloorModeGuard(bool v) : prev(floorOnChamferMode()) { floorOnChamferMode() = v; }
        FloorModeGuard(FloorModeGuard const&) = delete;
        FloorModeGuard(FloorModeGuard&&) = delete;
        auto operator=(FloorModeGuard const&) -> FloorModeGuard& = delete;
        auto operator=(FloorModeGuard&&) -> FloorModeGuard& = delete;
        ~FloorModeGuard() { floorOnChamferMode() = prev; }
    };
    struct CisFitnessModeGuard {
        bool prev;
        explicit CisFitnessModeGuard(bool v) : prev(cisFitnessMode()) { cisFitnessMode() = v; }
        CisFitnessModeGuard(CisFitnessModeGuard const&) = delete;
        CisFitnessModeGuard(CisFitnessModeGuard&&) = delete;
        auto operator=(CisFitnessModeGuard const&) -> CisFitnessModeGuard& = delete;
        auto operator=(CisFitnessModeGuard&&) -> CisFitnessModeGuard& = delete;
        ~CisFitnessModeGuard() { cisFitnessMode() = prev; }
    };
    struct ArcLengthGuardRatioGuard {
        Operon::Scalar prev;
        explicit ArcLengthGuardRatioGuard(Operon::Scalar v) : prev(arcLengthGuardRatio()) { arcLengthGuardRatio() = v; }
        ArcLengthGuardRatioGuard(ArcLengthGuardRatioGuard const&) = delete;
        ArcLengthGuardRatioGuard(ArcLengthGuardRatioGuard&&) = delete;
        auto operator=(ArcLengthGuardRatioGuard const&) -> ArcLengthGuardRatioGuard& = delete;
        auto operator=(ArcLengthGuardRatioGuard&&) -> ArcLengthGuardRatioGuard& = delete;
        ~ArcLengthGuardRatioGuard() { arcLengthGuardRatio() = prev; }
    };
    struct DirectFinalistEvalGuard {
        bool prev;
        explicit DirectFinalistEvalGuard(bool v) : prev(directFinalistEvalMode()) { directFinalistEvalMode() = v; }
        DirectFinalistEvalGuard(DirectFinalistEvalGuard const&) = delete;
        DirectFinalistEvalGuard(DirectFinalistEvalGuard&&) = delete;
        auto operator=(DirectFinalistEvalGuard const&) -> DirectFinalistEvalGuard& = delete;
        auto operator=(DirectFinalistEvalGuard&&) -> DirectFinalistEvalGuard& = delete;
        ~DirectFinalistEvalGuard() { directFinalistEvalMode() = prev; }
    };
    struct PointClampGuard {
        bool prev;
        explicit PointClampGuard(bool v) : prev(pointClampMode()) { pointClampMode() = v; }
        PointClampGuard(PointClampGuard const&) = delete;
        PointClampGuard(PointClampGuard&&) = delete;
        auto operator=(PointClampGuard const&) -> PointClampGuard& = delete;
        auto operator=(PointClampGuard&&) -> PointClampGuard& = delete;
        ~PointClampGuard() { pointClampMode() = prev; }
    };
    struct ConstMagnitudeComplexityGuard {
        bool prev;
        explicit ConstMagnitudeComplexityGuard(bool v)
            : prev(constMagnitudeComplexityMode()) { constMagnitudeComplexityMode() = v; }
        ConstMagnitudeComplexityGuard(ConstMagnitudeComplexityGuard const&) = delete;
        ConstMagnitudeComplexityGuard(ConstMagnitudeComplexityGuard&&) = delete;
        auto operator=(ConstMagnitudeComplexityGuard const&) -> ConstMagnitudeComplexityGuard& = delete;
        auto operator=(ConstMagnitudeComplexityGuard&&) -> ConstMagnitudeComplexityGuard& = delete;
        ~ConstMagnitudeComplexityGuard() { constMagnitudeComplexityMode() = prev; }
    };
    struct PriorMinimalComplexityGuard {
        bool prev;
        explicit PriorMinimalComplexityGuard(bool v)
            : prev(priorMinimalComplexityMode()) { priorMinimalComplexityMode() = v; }
        PriorMinimalComplexityGuard(PriorMinimalComplexityGuard const&) = delete;
        PriorMinimalComplexityGuard(PriorMinimalComplexityGuard&&) = delete;
        auto operator=(PriorMinimalComplexityGuard const&) -> PriorMinimalComplexityGuard& = delete;
        auto operator=(PriorMinimalComplexityGuard&&) -> PriorMinimalComplexityGuard& = delete;
        ~PriorMinimalComplexityGuard() { priorMinimalComplexityMode() = prev; }
    };
    struct FastNonDominatedSortGuard {
        bool prev;
        explicit FastNonDominatedSortGuard(bool v)
            : prev(fastNonDominatedSortMode()) { fastNonDominatedSortMode() = v; }
        FastNonDominatedSortGuard(FastNonDominatedSortGuard const&) = delete;
        FastNonDominatedSortGuard(FastNonDominatedSortGuard&&) = delete;
        auto operator=(FastNonDominatedSortGuard const&) -> FastNonDominatedSortGuard& = delete;
        auto operator=(FastNonDominatedSortGuard&&) -> FastNonDominatedSortGuard& = delete;
        ~FastNonDominatedSortGuard() { fastNonDominatedSortMode() = prev; }
    };
    struct UseExactNlogNSortGuard {
        bool prev;
        explicit UseExactNlogNSortGuard(bool v)
            : prev(useExactNlogNSortMode()) { useExactNlogNSortMode() = v; }
        UseExactNlogNSortGuard(UseExactNlogNSortGuard const&) = delete;
        UseExactNlogNSortGuard(UseExactNlogNSortGuard&&) = delete;
        auto operator=(UseExactNlogNSortGuard const&) -> UseExactNlogNSortGuard& = delete;
        auto operator=(UseExactNlogNSortGuard&&) -> UseExactNlogNSortGuard& = delete;
        ~UseExactNlogNSortGuard() { useExactNlogNSortMode() = prev; }
    };

    struct FamilyProbeResult {
        Operon::Scalar polyCorr{0};
        Operon::Scalar trigCorr{0};
    };

    static auto familyProbe(std::vector<Operon::Scalar> const& x,
                             std::vector<Operon::Scalar> const& y) -> FamilyProbeResult
    {
        FamilyProbeResult out;
        auto const n = x.size();
        if (n < 4 || y.size() != n) { return out; }

        auto mean = [](std::vector<Operon::Scalar> const& v) -> double {
            double s = 0;
            for (auto e : v) { s += static_cast<double>(e); }
            return s / static_cast<double>(v.size());
        };
        double const mx = mean(x);
        double const my = mean(y);

        double sumPolyX = 0, sumPolyY = 0, sumSinX = 0, sumSinY = 0, sumCosX = 0, sumCosY = 0;
        double sumPoly2 = 0, sumSin2 = 0, sumCos2 = 0, sumRes2 = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double const t = static_cast<double>(i) / static_cast<double>(n - 1);
            double const rx = static_cast<double>(x[i]) - mx;
            double const ry = static_cast<double>(y[i]) - my;
            double const poly = (t * t) - (1.0 / 3.0);
            double const sn = std::sin(2.0 * M_PI * t);
            double const cs = std::cos(2.0 * M_PI * t);
            sumPolyX += poly * rx; sumPolyY += poly * ry;
            sumSinX  += sn * rx;   sumSinY  += sn * ry;
            sumCosX  += cs * rx;   sumCosY  += cs * ry;
            sumPoly2 += poly * poly;
            sumSin2  += sn * sn;
            sumCos2  += cs * cs;
            sumRes2  += rx * rx + ry * ry;
        }
        if (!(sumRes2 > 0.0)) { return out; }
        double const resNorm = std::sqrt(sumRes2);
        auto safeCorr = [&](double num2, double basisNorm2) -> double {
            if (!(basisNorm2 > 0.0)) { return 0.0; }
            double const num = std::sqrt(num2);
            double const denom = std::sqrt(basisNorm2) * resNorm;
            return (denom > 0.0) ? std::min(1.0, num / denom) : 0.0;
        };
        double const polyC = safeCorr(sumPolyX * sumPolyX + sumPolyY * sumPolyY, sumPoly2);
        double const sinC  = safeCorr(sumSinX * sumSinX + sumSinY * sumSinY, sumSin2);
        double const cosC  = safeCorr(sumCosX * sumCosX + sumCosY * sumCosY, sumCos2);

        out.polyCorr = std::isfinite(polyC) ? static_cast<Operon::Scalar>(polyC) : Operon::Scalar{0};
        double const trigC = std::max(sinC, cosC);
        out.trigCorr = std::isfinite(trigC) ? static_cast<Operon::Scalar>(trigC) : Operon::Scalar{0};
        return out;
    }

    struct FamilyWeights {
        Operon::Scalar poly{1};
        Operon::Scalar trig{1};
    };

    static constexpr Operon::Scalar kFamilyBiasCeiling = Operon::Scalar{0.25};

    static auto familyWeightsFromProbe(FamilyProbeResult const& probe, std::size_t sf2N = 0) -> FamilyWeights
    {
        Operon::Scalar strength = std::abs(probe.polyCorr - probe.trigCorr);
        if (!std::isfinite(strength)) { strength = 0; }
        if (sf2N > 0) {
            auto const se = static_cast<Operon::Scalar>(std::sqrt(3.0 / static_cast<double>(sf2N)));
            strength = std::max(Operon::Scalar{0}, strength - se);
        } else {
            strength = std::min(kFamilyBiasCeiling, std::max(Operon::Scalar{0}, strength));
        }

        FamilyWeights w;
        if (probe.polyCorr > probe.trigCorr) {
            w.poly = Operon::Scalar{1} + strength;
            w.trig = Operon::Scalar{1} - strength;
        } else if (probe.trigCorr > probe.polyCorr) {
            w.trig = Operon::Scalar{1} + strength;
            w.poly = Operon::Scalar{1} - strength;
        }
        return w;
    }

    static constexpr std::size_t kFamilyBiasBase = 1000;

    static void applyOperatorFamilyBias(Operon::PrimitiveSet& pset, FamilyWeights const& w)
    {
        auto const squareHash = static_cast<Operon::Hash>(Operon::BuiltinOp::Square);
        auto const sinHash = static_cast<Operon::Hash>(Operon::BuiltinOp::Sin);
        auto const cosHash = static_cast<Operon::Hash>(Operon::BuiltinOp::Cos);

        std::vector<Operon::Hash> hashes;
        hashes.reserve(pset.Primitives().size());
        for (auto const& kv : pset.Primitives()) { hashes.push_back(kv.first); }

        for (auto h : hashes) {
            if (!pset.IsEnabled(h) || pset.Frequency(h) == 0) { continue; }
            double weight = 1.0;
            if (h == squareHash) { weight = static_cast<double>(w.poly); }
            else if (h == sinHash || h == cosHash) { weight = static_cast<double>(w.trig); }
            auto freq = static_cast<std::size_t>(std::llround(static_cast<double>(kFamilyBiasBase) * weight));
            pset.SetFrequency(h, std::max<std::size_t>(1, freq));
        }
    }

    static auto noiseGuidedLeafPick(Operon::RandomGenerator& rng, Operon::Tree& tree,
                                     std::vector<Operon::Scalar> const& perPointResid)
        -> Operon::Node*
    {
        auto& nodes = tree.Nodes();
        std::vector<std::size_t> leafIdx;
        leafIdx.reserve(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].IsLeaf()) { leafIdx.push_back(i); }
        }
        if (leafIdx.empty()) { return nullptr; }
        if (leafIdx.size() == 1 || perPointResid.empty()) {
            std::uniform_int_distribution<std::size_t> ud(0, leafIdx.size() - 1);
            return &nodes[leafIdx[ud(rng)]];
        }
        double maxR = 0.0;
        double sumR = 0.0;
        bool allFinite = true;
        for (auto v : perPointResid) {
            double d = static_cast<double>(v);
            if (!std::isfinite(d)) { allFinite = false; break; }
            maxR = std::max(maxR, d);
            sumR += d;
        }
        if (!allFinite || !(sumR >= 0.0)) {
            std::uniform_int_distribution<std::size_t> ud(0, leafIdx.size() - 1);
            return &nodes[leafIdx[ud(rng)]];
        }
        double const meanR = sumR / static_cast<double>(perPointResid.size());
        double const severity = (meanR > 1e-12)
            ? std::clamp(maxR / meanR - 1.0, 0.0, 4.0)
            : 0.0;
        double const temp = 0.35 * severity;

        std::vector<double> w(leafIdx.size());
        double wsum = 0.0;
        for (std::size_t k = 0; k < leafIdx.size(); ++k) {
            double const logit = temp * static_cast<double>(nodes[leafIdx[k]].Level);
            w[k] = std::exp(logit);
            wsum += w[k];
        }
        if (!(std::isfinite(wsum) && wsum > 0.0)) {
            std::uniform_int_distribution<std::size_t> ud(0, leafIdx.size() - 1);
            return &nodes[leafIdx[ud(rng)]];
        }
        std::uniform_real_distribution<double> ud(0.0, wsum);
        double const r = ud(rng);
        double acc = 0.0;
        for (std::size_t k = 0; k < leafIdx.size(); ++k) {
            acc += w[k];
            if (r <= acc) { return &nodes[leafIdx[k]]; }
        }
        return &nodes[leafIdx.back()];
    }

    struct Sf2JeffreysCoeffInit : public Operon::CoefficientInitializerBase {
        double logLo{0.0};
        double logHi{0.0};
        auto operator()(Operon::RandomGenerator& random, Operon::Tree& tree) const -> void override
        {
            std::uniform_real_distribution<double> u(logLo, logHi);
            std::uniform_int_distribution<int> sgn(0, 1);
            for (auto& node : tree.Nodes()) {
                if (node.IsLeaf()) {
                    double const mag = std::exp(u(random));
                    auto const v = static_cast<Operon::Scalar>(sgn(random) == 0 ? mag : -mag);
                    if (std::isfinite(static_cast<double>(v))) { node.Value = v; }
                }
            }
        }
    };

    auto sf2LevelLeafPick(Operon::RandomGenerator& rng, Operon::Tree& tree) const -> Operon::Node*
    {
        auto& nodes = tree.Nodes();
        std::vector<std::size_t> leafIdx;
        leafIdx.reserve(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) { if (nodes[i].IsLeaf()) { leafIdx.push_back(i); } }
        if (leafIdx.empty()) { return nullptr; }
        std::vector<double> w(leafIdx.size());
        double wsum = 0.0;
        for (std::size_t k = 0; k < leafIdx.size(); ++k) {
            auto const L = static_cast<std::size_t>(nodes[leafIdx[k]].Level);
            double const sc = (L < sf2LvSucc_.size()) ? sf2LvSucc_[L] : 0.0;
            double const nt = (L < sf2LvTrials_.size()) ? sf2LvTrials_[L] : 0.0;
            w[k] = (sc + 1.0) / (nt + 2.0);
            wsum += w[k];
        }
        std::uniform_real_distribution<double> ud(0.0, wsum);
        double const r = ud(rng);
        double acc = 0.0;
        for (std::size_t k = 0; k < leafIdx.size(); ++k) {
            acc += w[k];
            if (r <= acc) { return &nodes[leafIdx[k]]; }
        }
        return &nodes[leafIdx.back()];
    }

    struct MutationSet {
        Operon::OnePointMutation<std::normal_distribution<Operon::Scalar>> onePoint;
        std::unique_ptr<Operon::ChangeFunctionMutation> changeFunc;
        std::unique_ptr<Operon::ReplaceSubtreeMutation> replaceSub;
        std::unique_ptr<Operon::RemoveSubtreeMutation> removeSub;
        std::unique_ptr<Operon::InsertSubtreeMutation> insertSub;
        std::unique_ptr<BinaryInsertSubtreeMutation> binInsertSub;
        Operon::MutatorBase const* insertOp{nullptr};
        std::array<std::unique_ptr<CountingMutation>, 5> diag;
        Operon::MultiMutation multi;
        bool insertIncluded{true};

        MutationSet(Operon::PrimitiveSet const& pset,
                    Operon::CreatorBase const& creator,
                    Operon::CoefficientInitializerBase const& coeffInit,
                    std::size_t maxDepth, std::size_t maxLength,
                    bool binaryInsert = false, MutOpDiagCounters* diagCtr = nullptr,
                    bool dropInsert = false)
        {
            bool const excludeInsert = dropInsert && !binaryInsert;
            insertIncluded = !excludeInsert;
            if (binaryInsert || diagCtr != nullptr || excludeInsert) {
                if (diagCtr != nullptr) { diagCtr->maxLength = maxLength; }
                onePoint.ParameterizeDistribution(Operon::Scalar{0}, Operon::Scalar{3});
                changeFunc = std::make_unique<Operon::ChangeFunctionMutation>(pset);
                replaceSub = std::make_unique<Operon::ReplaceSubtreeMutation>(
                    &creator, &coeffInit, maxDepth, maxLength);
                removeSub = std::make_unique<Operon::RemoveSubtreeMutation>(
                    &creator, &coeffInit, maxDepth);
                insertSub = std::make_unique<Operon::InsertSubtreeMutation>(
                    &creator, &coeffInit, maxDepth, maxLength);
                insertOp = insertSub.get();
                if (binaryInsert) {
                    binInsertSub = std::make_unique<BinaryInsertSubtreeMutation>(&creator, &coeffInit, maxDepth, maxLength);
                    insertOp = binInsertSub.get();
                }
                std::array<Operon::MutatorBase const*, 5> ops{&onePoint, changeFunc.get(), replaceSub.get(),
                                                              removeSub.get(), insertOp};
                for (int k = 0; k < 5; ++k) {
                    if (k == 4 && excludeInsert) { continue; }
                    Operon::MutatorBase const* op = ops[static_cast<std::size_t>(k)];
                    if (diagCtr != nullptr) {
                        diag[static_cast<std::size_t>(k)] = std::make_unique<CountingMutation>(op, k, diagCtr);
                        op = diag[static_cast<std::size_t>(k)].get();
                    }
                    multi.Add(op, 1.0);
                }
                return;
            }
            onePoint.ParameterizeDistribution(Operon::Scalar{0}, Operon::Scalar{3});
            changeFunc = std::make_unique<Operon::ChangeFunctionMutation>(pset);
            replaceSub = std::make_unique<Operon::ReplaceSubtreeMutation>(
                &creator, &coeffInit, maxDepth, maxLength);
            removeSub = std::make_unique<Operon::RemoveSubtreeMutation>(
                &creator, &coeffInit, maxDepth);
            insertSub = std::make_unique<Operon::InsertSubtreeMutation>(
                &creator, &coeffInit, maxDepth, maxLength);
            multi.Add(&onePoint, 1.0);
            multi.Add(changeFunc.get(), 1.0);
            multi.Add(replaceSub.get(), 1.0);
            multi.Add(removeSub.get(), 1.0);
            multi.Add(insertSub.get(), 1.0);
            insertOp = insertSub.get();
        }
    };

    struct SearchState {
        std::size_t M{0};
        std::size_t N{0};
        std::size_t lowM{0};
        std::size_t midM{0};

        std::vector<Operon::Scalar> fullFlat;
        std::unique_ptr<Operon::Dataset> fullDs;

        std::vector<Operon::Scalar> lowFlat;
        std::unique_ptr<Operon::Dataset> lowDs;

        std::vector<Operon::Scalar> midFlat;
        std::unique_ptr<Operon::Dataset> midDs;

        std::unique_ptr<Operon::Problem> prob;
        std::vector<Operon::Hash> inputs;
        std::unique_ptr<Operon::BalancedTreeCreator> creator;
        Operon::NormalCoefficientInitializer coeffInit;
        Sf2JeffreysCoeffInit sf2CoeffInit;
        Operon::CoefficientInitializerBase const* coeffInitP{&coeffInit};
        int sf2GuidedLevel{-1};
        std::unique_ptr<Operon::SubtreeCrossover> crossover;
        std::unique_ptr<MutationSet> mut;

        Operon::RandomGenerator rng;
        std::mt19937 sel;

        bool lowRes{false};
        bool stagedMid{false};
        int annealSwitchGen{-1};
        std::size_t switchGen1{std::numeric_limits<std::size_t>::max()};
        std::size_t switchGen2{std::numeric_limits<std::size_t>::max()};

        Operon::Dataset const* activeDs{nullptr};
        std::size_t activeM{0};
        Operon::Range range{0, 0};

        std::vector<PairIndividual> pop;
        std::vector<PairIndividual> archive;
        std::size_t archiveCap{0};
        std::vector<PairIndividual> hallOfFame;

        std::vector<std::vector<PairIndividual>> layers;

        std::size_t currentGen{0};
        int lastGen{0};
        bool truncated{false};

        bool inSimplify{false};
        std::size_t simplifyLeft{0};
        std::size_t simplifyGensRan{0};

        double semBound{0.0};

        std::chrono::steady_clock::time_point fitStart;

        explicit SearchState(std::size_t seed)
            : rng(seed)
            , sel(static_cast<std::uint32_t>(seed) ^ 0x9e3779b9u)
        {}
    };

    std::unique_ptr<SearchState> state_;

    void reevalAll(SearchState& s, std::vector<PairIndividual>& v) const
    {
        auto const nWorkers = workerCount();
        if (oneResGpuRoute()) {
            if (oneResDiagEnabled()) { oneResDiagCompare("reevalAll", v, *s.activeDs, s.range); }
            std::vector<PairIndividual*> ptrs(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) { ptrs[i] = &v[i]; }
            rescoreBatched(ptrs, *s.activeDs, s.range, nWorkers);
            return;
        }
        if (oneResDiagEnabled()) { oneResDiagCompare("reevalAll", v, *s.activeDs, s.range); }
        parallelApply(v.size(), nWorkers, [&](std::size_t i) {
            evaluateGuarded(v[i], *s.activeDs, s.range);
        });
    }

    auto oneResGpuRoute() const -> bool
    {
        return oneResOn() && !oneResBit(16) && config_.useBatchedGpuPointSeg && cisFitnessMode() && kPointSegGpuCompiledIn;
    }

    void rescoreBatched(std::vector<PairIndividual*> const& inds, Operon::Dataset const& ds,
                        Operon::Range range, std::size_t nWorkers) const
    {
        if (inds.empty()) { return; }
        std::vector<std::vector<Operon::Scalar>> gX, gY;
        CisPointSegBatch b, b48;
        bool const ok = evaluateBatchedPointSegPrepassPtrs(inds, ds, range, gX, gY, b, &b48);
        parallelApply(inds.size(), nWorkers, [&](std::size_t k) {
            if (ok && k < gX.size()) {
                evaluateGuardedFromGrid(*inds[k], ds, range, gX[k], gY[k], &b, k, &b48);
            } else {
                evaluateGuarded(*inds[k], ds, range);
            }
        });
    }

    static auto oneResDiagEnabled() -> bool
    {
        static bool const v = std::getenv("SKETCH_T11_DIAG") != nullptr;
        return v;
    }
    void oneResDiagCompare(char const* site, std::vector<PairIndividual> const& v,
                           Operon::Dataset const& ds, Operon::Range range) const
    {
        std::vector<PairIndividual const*> src;
        for (auto const& ind : v) { src.push_back(&ind); }
        oneResDiagComparePtrs(site, src, ds, range);
    }
    void oneResDiagComparePtrs(char const* site, std::vector<PairIndividual const*> const& src,
                               Operon::Dataset const& ds, Operon::Range range) const
    {
        if (src.empty() || !config_.useBatchedGpuPointSeg || !cisFitnessMode()) { return; }
        std::size_t const P = src.size();
        std::vector<PairIndividual> gpuC(P), cpuLegacy(P), cpuOne(P);
        for (std::size_t i = 0; i < P; ++i) { gpuC[i] = *src[i]; cpuLegacy[i] = *src[i]; cpuOne[i] = *src[i]; }
        std::vector<std::vector<Operon::Scalar>> gX, gY;
        CisPointSegBatch b;
        bool const ok = evaluateBatchedPointSegPrepass(gpuC, ds, range, gX, gY, b);
        if (!ok) { std::fprintf(stderr, "[T11DIAG] %s gpu prepass failed\n", site); return; }
        std::size_t gridMismatch = 0; double gridMaxAbs = 0;
        {
            std::vector<Operon::Scalar> cx, cy;
            for (std::size_t i = 0; i < P && i < gX.size(); ++i) {
                try {
                    interpInto(src[i]->x, src[i]->x.GetCoefficients(), ds, range, cx);
                    interpInto(src[i]->y, src[i]->y.GetCoefficients(), ds, range, cy);
                } catch (...) { continue; }
                bool mm = cx.size() != gX[i].size() || cy.size() != gY[i].size();
                for (std::size_t j = 0; !mm && j < cx.size(); ++j) {
                    double const a = std::abs(static_cast<double>(cx[j]) - gX[i][j]), c = std::abs(static_cast<double>(cy[j]) - gY[i][j]);
                    if (std::isfinite(a) && std::isfinite(c)) { gridMaxAbs = std::max({gridMaxAbs, a, c}); }
                    if (cx[j] != gX[i][j] || cy[j] != gY[i][j]) { mm = true; }
                }
                if (mm) { ++gridMismatch; }
            }
        }
        std::fprintf(stderr, "[T11DIAGGRID] site=%s gridMismatch=%zu/%zu gridMaxAbs=%.6g\n", site, gridMismatch, P, gridMaxAbs);
        int const saveMode = config_.oneResidualMode;
        auto& cfgMut = const_cast<MultitreeConfig&>(config_);
        for (std::size_t i = 0; i < P; ++i) { evaluateGuardedFromGrid(gpuC[i], ds, range, gX[i], gY[i], &b, i); }
        cfgMut.oneResidualMode = saveMode & ~1;
        for (std::size_t i = 0; i < P; ++i) { evaluateGuarded(cpuLegacy[i], ds, range); }
        cfgMut.oneResidualMode = saveMode | 1;
        for (std::size_t i = 0; i < P; ++i) { evaluateGuarded(cpuOne[i], ds, range); }
        cfgMut.oneResidualMode = saveMode;
        auto stats = [&](std::vector<PairIndividual> const& a, std::vector<PairIndividual> const& c,
                         double& maxAbs, double& maxRel, long long& inv, long long& pairs, std::size_t& nDiff) {
            maxAbs = 0; maxRel = 0; inv = 0; pairs = 0; nDiff = 0;
            for (std::size_t i = 0; i < P; ++i) {
                double const x = a[i].shapeErr, y = c[i].shapeErr;
                if (!std::isfinite(x) || !std::isfinite(y) || x >= 1e30 || y >= 1e30) { continue; }
                double const d = std::abs(x - y);
                if (d > 0) { ++nDiff; }
                maxAbs = std::max(maxAbs, d);
                maxRel = std::max(maxRel, d / std::max(std::abs(x), 1e-30));
            }
            for (std::size_t i = 0; i < P; ++i) {
                for (std::size_t j = i + 1; j < P; ++j) {
                    double const ai = a[i].shapeErr, aj = a[j].shapeErr, ci = c[i].shapeErr, cj = c[j].shapeErr;
                    if (!std::isfinite(ai) || !std::isfinite(aj) || !std::isfinite(ci) || !std::isfinite(cj)) { continue; }
                    if (ai >= 1e30 || aj >= 1e30) { continue; }
                    ++pairs;
                    if ((ai < aj && ci > cj) || (ai > aj && ci < cj)) { ++inv; }
                }
            }
        };
        double mA1, mR1, mA2, mR2; long long i1, p1, i2, p2; std::size_t d1, d2;
        stats(gpuC, cpuLegacy, mA1, mR1, i1, p1, d1);
        stats(gpuC, cpuOne, mA2, mR2, i2, p2, d2);
        std::vector<std::size_t> ord(P);
        for (std::size_t i = 0; i < P; ++i) { ord[i] = i; }
        std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t c) { return gpuC[a].shapeErr < gpuC[c].shapeErr; });
        std::size_t const top = std::max<std::size_t>(2, P / 10);
        auto topInv = [&](std::vector<PairIndividual> const& c) {
            long long inv = 0;
            for (std::size_t a = 0; a < top && a < P; ++a) {
                for (std::size_t b = a + 1; b < top && b < P; ++b) {
                    std::size_t const i = ord[a], j = ord[b];
                    if (gpuC[i].shapeErr < gpuC[j].shapeErr && c[i].shapeErr > c[j].shapeErr) { ++inv; }
                }
            }
            return inv;
        };
        std::size_t w = 0; double wr = -1;
        for (std::size_t i = 0; i < P; ++i) {
            double const x = gpuC[i].shapeErr, y = cpuOne[i].shapeErr;
            if (!std::isfinite(x) || !std::isfinite(y) || x >= 1e30 || y >= 1e30) { continue; }
            double const r = std::abs(x - y) / std::max(std::abs(x), 1e-30);
            if (r > wr) { wr = r; w = i; }
        }
        std::size_t oneSided1 = 0, oneSided2 = 0;
        auto noCliff = [&](std::vector<PairIndividual> const& c, double& mA, double& mR, std::size_t& nc, std::size_t& os) {
            mA = 0; mR = 0; nc = 0; os = 0;
            for (std::size_t i = 0; i < P; ++i) {
                double const x = gpuC[i].shapeErr, y = c[i].shapeErr;
                bool const bx = !std::isfinite(x) || x >= 1e30, by = !std::isfinite(y) || y >= 1e30;
                if (bx != by) { ++os; }
                if (!std::isfinite(x) || !std::isfinite(y) || x >= 1e30 || y >= 1e30) { continue; }
                double const d = std::abs(x - y);
                if (d > 1e3) { ++nc; continue; }
                mA = std::max(mA, d); mR = std::max(mR, d / std::max(std::abs(x), 1e-30));
            }
        };
        double cA1, cR1, cA2, cR2; std::size_t nc1, nc2;
        noCliff(cpuLegacy, cA1, cR1, nc1, oneSided1);
        noCliff(cpuOne, cA2, cR2, nc2, oneSided2);
        std::fprintf(stderr, "[T11DIAGNC] site=%s legacy: cliff=%zu maxAbs=%.6g maxRel=%.6g | one: cliff=%zu maxAbs=%.6g maxRel=%.6g"
                     " | oneSidedFail legacy=%zu one=%zu law=%d\n",
                     site, nc1, cA1, cR1, nc2, cA2, cR2, oneSided1, oneSided2, cisExactModeDefault());
        std::fprintf(stderr,
            "[T11DIAG] site=%s P=%zu N=%zu | GPU vs CPU-legacy: nDiff=%zu maxAbs=%.6g maxRel=%.6g inversions=%lld/%lld"
            " | GPU vs CPU-one: nDiff=%zu maxAbs=%.6g maxRel=%.6g inversions=%lld/%lld"
            " | topDecile inv legacy=%lld one=%lld | worst gpu=%.9g legacy=%.9g one=%.9g cplx=%g\n",
            site, P, dataX_.size(), d1, mA1, mR1, i1, p1, d2, mA2, mR2, i2, p2, topInv(cpuLegacy), topInv(cpuOne),
            static_cast<double>(gpuC[w].shapeErr), static_cast<double>(cpuLegacy[w].shapeErr),
            static_cast<double>(cpuOne[w].shapeErr), static_cast<double>(gpuC[w].complexity));
    }

    auto trySpectralSplice(SearchState& s, Operon::Tree const& tree, PairIndividual const& parent,
                            bool isX) const -> std::optional<Operon::Tree>
    {
        if (s.activeDs == nullptr) { return std::nullopt; }
        auto const& ds = *s.activeDs;
        auto const range = s.range;
        std::size_t const m = range.Size();
        if (m < 8) { return std::nullopt; }

        auto cx = parent.x.GetCoefficients();
        auto cy = parent.y.GetCoefficients();
        std::vector<Operon::Scalar> fgx, fgy;
        interpInto(parent.x, cx, ds, range, fgx);
        interpInto(parent.y, cy, ds, range, fgy);
        if (fgx.size() != m || fgy.size() != m) { return std::nullopt; }

        auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
        if (tValues.size() != m) { return std::nullopt; }
        std::vector<double> res(m), tD(m);
        for (std::size_t i = 0; i < m; ++i) {
            double const rx = static_cast<double>(dataX_[i]) - static_cast<double>(fgx[i]);
            double const ry = static_cast<double>(dataY_[i]) - static_cast<double>(fgy[i]);
            if (!std::isfinite(rx) || !std::isfinite(ry)) { return std::nullopt; }
            res[i] = isX ? rx : ry;
            double const tv = static_cast<double>(tValues[i]);
            if (!std::isfinite(tv)) { return std::nullopt; }
            tD[i] = tv;
        }
        double const span = tD.back() - tD.front();
        double const useSpan = (span > 0.0 && std::isfinite(span)) ? span : 1.0;

        std::vector<double> seedSigma;
        if (config_.useNoiseWeightedFrequencySeed) {
            std::vector<Operon::Scalar> const* extSig = nullptr;
            std::size_t Ksig = 0;
            Operon::Scalar constSig2 = Operon::Scalar{0};
            if (config_.useExternalPointSigma2 && !externalPointSigma2_.empty()) {
                extSig = &externalPointSigma2_; Ksig = externalPointSigma2_.size();
            } else if (noise_.size() > 0) {
                extSig = &noise_.s2; Ksig = noise_.size();
            } else {
                constSig2 = noiseMeanS2_;
            }
            seedSigma.resize(m);
            for (std::size_t i = 0; i < m; ++i) {
                Operon::Scalar const s2 = (extSig != nullptr)
                    ? (*extSig)[noiseIndexAtArc(i, m, Ksig)] : constSig2;
                seedSigma[i] = std::sqrt(std::max(static_cast<double>(s2), 0.0));
            }
        }
        freqdisc::FrequencyProposal prop;
        if (config_.useMatrixPencilFrequencySeed) {
            freqdisc::FrequencyProposal unusedY;
            freqdisc::proposeFrequencyPencilXY(res, {}, useSpan, prop, unusedY);
        } else {
            prop = config_.useNoiseWeightedFrequencySeed
                ? freqdisc::proposeFrequencyWeighted(res, useSpan, seedSigma)
                : freqdisc::proposeFrequency(res, useSpan);
        }
        if (!prop.valid || !std::isfinite(prop.omega)) { return std::nullopt; }

        auto const qf = freqdisc::quadratureFit(tD, res, prop.omega);
        if (!qf.valid || !std::isfinite(qf.a) || !std::isfinite(qf.b)) { return std::nullopt; }

        auto const& oldNodes = tree.Nodes();
        if (oldNodes.empty()) { return std::nullopt; }

        using BO = Operon::BuiltinOp;
        auto varNode = [&]() {
            Operon::Node n(Operon::NodeType::Variable);
            n.HashValue = s.inputs[0];
            n.Value = Operon::Scalar{1};
            return n;
        };

        std::vector<Operon::Node> out;
        out.reserve(oldNodes.size() + 13);
        out.insert(out.end(), oldNodes.begin(), oldNodes.end());
        out.push_back(Operon::Node::Constant(prop.omega));
        out.push_back(varNode());
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Mul), 2));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Sin), 1));
        out.push_back(Operon::Node::Constant(qf.a));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Mul), 2));
        out.push_back(Operon::Node::Constant(prop.omega));
        out.push_back(varNode());
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Mul), 2));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Cos), 1));
        out.push_back(Operon::Node::Constant(qf.b));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Mul), 2));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Add), 2));
        out.push_back(Operon::Node::Function(static_cast<Operon::Hash>(BO::Add), 2));

        Operon::Tree child(out);
        child.UpdateNodes();
        return child;
    }

    auto jacobianLeafPick(SearchState& s, Operon::Tree& tree, Operon::Tree const& sibling,
                          bool isX, Operon::Dataset const& ds, Operon::Range range) const
        -> Operon::Node*
    {
        auto& nodes = tree.Nodes();
        std::vector<std::size_t> leafIdx;
        leafIdx.reserve(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].IsLeaf()) { leafIdx.push_back(i); }
        }
        if (leafIdx.empty()) { return nullptr; }
        auto uniformPick = [&]() -> Operon::Node* {
            std::uniform_int_distribution<std::size_t> ud(0, leafIdx.size() - 1);
            return &nodes[leafIdx[ud(s.rng)]];
        };

        std::vector<std::size_t> constIdx;
        constIdx.reserve(leafIdx.size());
        for (auto li : leafIdx) {
            if (nodes[li].Type == Operon::NodeType::Constant) { constIdx.push_back(li); }
        }
        if (constIdx.size() < 2) { return uniformPick(); }

        auto const N = dataX_.size();
        if (N == 0) { return uniformPick(); }

        auto buildCandPts = [&](Operon::Tree const& xt, Operon::Tree const& yt,
                                 std::vector<Vec2>& out) -> bool {
            auto cx = xt.GetCoefficients();
            auto cy = yt.GetCoefficients();
            auto gx = interp(xt, cx, ds, range);
            auto gy = interp(yt, cy, ds, range);
            if (gx.size() != gy.size() || gx.empty()) { return false; }
            for (std::size_t j = 0; j < gx.size(); ++j) {
                if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return false; }
            }
            std::vector<Operon::Scalar> nx(gx), ny(gy);
            normalizeInPlace(nx, ny);
            out.resize(nx.size());
            for (std::size_t j = 0; j < nx.size(); ++j) { out[j] = Vec2{nx[j], ny[j]}; }
            return true;
        };

        std::vector<Vec2> strokePts(N);
        for (std::size_t i = 0; i < N; ++i) { strokePts[i] = Vec2{dataX_[i], dataY_[i]}; }

        std::vector<Vec2> candBase;
        bool const baseOk = isX ? buildCandPts(tree, sibling, candBase)
                                 : buildCandPts(sibling, tree, candBase);
        if (!baseOk) { return uniformPick(); }
        auto const psBase = pointSegOp(strokePts, candBase, false, false);
        if (!psBase.ok || psBase.residual.size() != N) { return uniformPick(); }
        std::vector<Operon::Scalar> rBase(psBase.residual);
        for (std::size_t i = 0; i < N; ++i) {
            if (!std::isfinite(rBase[i])) { return uniformPick(); }
        }

        std::vector<double> s2(N, 1.0);
        if (config_.useExternalPointSigma2 && !externalPointSigma2_.empty()) {
            auto const& rlut = noiseRegionLUT(N, externalPointSigma2_.size());
            for (std::size_t i = 0; i < N; ++i) {
                s2[i] = static_cast<double>(externalPointSigma2_[rlut[originalArcIndex(i)]]);
            }
        } else if (noise_.size() > 0U) {
            auto const& rlut = noiseRegionLUT(N, noise_.size());
            for (std::size_t i = 0; i < N; ++i) {
                s2[i] = static_cast<double>(noise_.s2[rlut[originalArcIndex(i)]]);
            }
        }
        for (std::size_t i = 0; i < N; ++i) {
            if (!std::isfinite(s2[i]) || !(s2[i] > 0.0)) { s2[i] = 1.0; }
        }

        std::vector<double> W(constIdx.size(), 0.0);
        for (std::size_t k = 0; k < constIdx.size(); ++k) {
            auto const nj = constIdx[k];
            Operon::Scalar const orig = nodes[nj].Value;
            Operon::Scalar const delta = std::max(Operon::Scalar{1e-3},
                Operon::Scalar{1e-3} * std::abs(orig));
            nodes[nj].Value = orig + delta;
            std::vector<Vec2> candPert;
            bool const pertOk = isX ? buildCandPts(tree, sibling, candPert)
                                     : buildCandPts(sibling, tree, candPert);
            PointSegmentResult psPert;
            bool okPert = false;
            if (pertOk) {
                psPert = pointSegOp(strokePts, candPert, false, false);
                okPert = psPert.ok && psPert.residual.size() == N;
            }
            nodes[nj].Value = orig;
            if (!okPert) { continue; }
            double wj = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                double const dr = (static_cast<double>(psPert.residual[i])
                                    - static_cast<double>(rBase[i])) / static_cast<double>(delta);
                double const term = std::abs(dr)
                    * (static_cast<double>(rBase[i]) * static_cast<double>(rBase[i])) / s2[i];
                if (std::isfinite(term)) { wj += term; }
            }
            if (std::isfinite(wj) && wj >= 0.0) { W[k] = wj; }
        }

        double wsum = 0.0;
        for (auto v : W) { wsum += v; }
        if (!(std::isfinite(wsum) && wsum > 0.0)) { return uniformPick(); }

        std::uniform_real_distribution<double> ud(0.0, wsum);
        double const r = ud(s.rng);
        double acc = 0.0;
        for (std::size_t k = 0; k < constIdx.size(); ++k) {
            acc += W[k];
            if (r <= acc) { return &nodes[constIdx[k]]; }
        }
        return &nodes[constIdx.back()];
    }

    struct GmLin {
        std::size_t N{0}, M{0};
        double invS{0};
        std::size_t iUmin{0}, iUmax{0}, iVmin{0}, iVmax{0};
        int sBranch{0};
        std::vector<double> ax, ay;
        std::vector<int> jb;
        std::vector<double> tau, nhx, nhy, sd;
        std::vector<double> mhx, mhy, sc;
        std::vector<double> r;
        std::vector<double> Q;
        std::size_t nq{0};
        double rr{0};
        double rPerp{0};
    };

    auto gmBuildLin(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range,
                    GmLin& L, std::vector<double>* jOut = nullptr) const -> bool
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        std::size_t const nx = cx.size();
        std::vector<Operon::Scalar> v = cx;
        v.insert(v.end(), cy.begin(), cy.end());
        std::size_t const n = v.size();
        std::vector<double> J;
        std::vector<int> corr;
        if (!lmBuild(ind, v, nx, ds, range, L.r, J, &corr)) { return false; }
        std::vector<Operon::Scalar> U, V;
        interpInto(ind.x, cx, ds, range, U);
        interpInto(ind.y, cy, ds, range, V);
        std::size_t const M = U.size();
        std::size_t const N = dataX_.size();
        std::size_t const m = N + M;
        if (M < 2 || V.size() != M || L.r.size() != m || corr.size() != m + 4) { return false; }
        L.N = N; L.M = M;
        L.iUmin = static_cast<std::size_t>(corr[m + 0]); L.iUmax = static_cast<std::size_t>(corr[m + 1]);
        L.iVmin = static_cast<std::size_t>(corr[m + 2]); L.iVmax = static_cast<std::size_t>(corr[m + 3]);
        double const Cx = 0.5 * (static_cast<double>(U[L.iUmax]) + static_cast<double>(U[L.iUmin]));
        double const Cy = 0.5 * (static_cast<double>(V[L.iVmax]) + static_cast<double>(V[L.iVmin]));
        double const spanU = static_cast<double>(U[L.iUmax]) - static_cast<double>(U[L.iUmin]);
        double const spanV = static_cast<double>(V[L.iVmax]) - static_cast<double>(V[L.iVmin]);
        double const S = std::max({spanU, spanV, 1e-9});
        if (!(S > 0.0) || !std::isfinite(S)) { return false; }
        L.sBranch = (spanU >= spanV && spanU > 1e-9) ? 1 : ((spanV > 1e-9) ? 2 : 0);
        L.invS = 1.0 / S;
        L.ax.resize(M); L.ay.resize(M);
        std::vector<double> qfx(M), qfy(M);
        for (std::size_t j = 0; j < M; ++j) {
            L.ax[j] = (static_cast<double>(U[j]) - Cx) * L.invS;
            L.ay[j] = (static_cast<double>(V[j]) - Cy) * L.invS;
            qfx[j] = static_cast<double>(static_cast<Operon::Scalar>(L.ax[j]));
            qfy[j] = static_cast<double>(static_cast<Operon::Scalar>(L.ay[j]));
        }
        auto const K = noise_.size();
        auto regionOf = [&](std::size_t i) -> std::size_t { return noiseIndexAtArc(originalArcIndex(i), N, K); };
        auto weightOf = [&](double u, std::size_t k) -> double {
            if (K == 0) { return 1.0; }
            double const s2 = static_cast<double>(noise_.s2[k]);
            double const nu = static_cast<double>(noise_.nu[k]);
            if (!(s2 > 0.0)) { return 1.0; }
            if (nu <= 0.0) { return 1.0 / s2; }
            return 0.5 * (nu + 2.0) / (nu * s2 + u);
        };
        constexpr double EPS = 1e-12;
        L.jb.assign(N, 0); L.tau.assign(N, 0.0); L.nhx.assign(N, 0.0); L.nhy.assign(N, 0.0); L.sd.assign(N, 0.0);
        if (oneResBit(4)) {
            thread_local OneResLmCore c;
            if (!oneResLmCoreEval(L.ax.data(), L.ay.data(), M, dataX_, dataY_, dataBlocks_, nullptr, c)) { return false; }
            for (std::size_t i = 0; i < N; ++i) {
                double const r = c.ps.residual[i];
                double tau = std::clamp(c.ps.tau[i], 0.0, static_cast<double>(M - 1));
                auto s0 = static_cast<std::size_t>(std::floor(tau));
                if (s0 > M - 2) { s0 = M - 2; }
                L.jb[i] = static_cast<int>(s0);
                L.tau[i] = std::clamp(tau - static_cast<double>(s0), 0.0, 1.0);
                L.nhx[i] = c.ps.normalX[i]; L.nhy[i] = c.ps.normalY[i];
                double const w = weightOf(r * r, c.dataRegion[i]);
                L.sd[i] = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(N)));
            }
            L.mhx.assign(M, 0.0); L.mhy.assign(M, 0.0); L.sc.assign(M, 0.0);
            for (std::size_t j = 0; j < M; ++j) {
                double const e = std::sqrt(c.e2[j]);
                if (e <= EPS) { continue; }
                double const w = weightOf(c.e2[j], c.curveRegion[j]);
                L.sc[j] = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(M)));
                L.mhx[j] = (L.ax[j] - c.zx[j]) / e; L.mhy[j] = (L.ay[j] - c.zy[j]) / e;
            }
        } else {
        for (std::size_t i = 0; i < N; ++i) {
            int const j = corr[i];
            if (j < 0 || static_cast<std::size_t>(j) + 1 >= M) { continue; }
            double const px = static_cast<double>(dataX_[i]), py = static_cast<double>(dataY_[i]);
            double const ax = qfx[j], ay = qfy[j], bx = qfx[j + 1], by = qfy[j + 1];
            double const ex = bx - ax, ey = by - ay;
            double const len2 = ex * ex + ey * ey;
            double t = 0.0;
            if (len2 > 1e-12) { t = std::clamp(((px - ax) * ex + (py - ay) * ey) / len2, 0.0, 1.0); }
            double const zx = ax + t * ex, zy = ay + t * ey;
            double const d2 = (px - zx) * (px - zx) + (py - zy) * (py - zy);
            double const d = std::sqrt(d2);
            double const w = weightOf(d2, regionOf(i));
            L.jb[i] = j; L.tau[i] = t;
            if (d <= EPS) { continue; }
            L.sd[i] = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(N)));
            L.nhx[i] = (px - zx) / d; L.nhy[i] = (py - zy) / d;
        }
        L.mhx.assign(M, 0.0); L.mhy.assign(M, 0.0); L.sc.assign(M, 0.0);
        for (std::size_t j = 0; j < M; ++j) {
            int const kb = corr[N + j];
            if (kb < 0 || static_cast<std::size_t>(kb) >= N) { continue; }
            double const ddx = qfx[j] - static_cast<double>(dataX_[kb]);
            double const ddy = qfy[j] - static_cast<double>(dataY_[kb]);
            double const e2 = ddx * ddx + ddy * ddy;
            double const e = std::sqrt(e2);
            if (e <= EPS) { continue; }
            double const w = weightOf(e2, regionOf(static_cast<std::size_t>(kb)));
            L.sc[j] = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(M)));
            L.mhx[j] = ddx / e; L.mhy[j] = ddy / e;
        }
        }
        double const rankTol = std::sqrt(std::numeric_limits<double>::epsilon());
        L.Q.clear(); L.nq = 0;
        std::vector<double> col(m);
        for (std::size_t k = 0; k < n; ++k) {
            double c0 = 0.0;
            for (std::size_t i = 0; i < m; ++i) { col[i] = J[i * n + k]; c0 += col[i] * col[i]; }
            if (!(c0 > 0.0) || !std::isfinite(c0)) { continue; }
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t q = 0; q < L.nq; ++q) {
                    double const* Qq = L.Q.data() + q * m;
                    double dot = 0.0;
                    for (std::size_t i = 0; i < m; ++i) { dot += Qq[i] * col[i]; }
                    for (std::size_t i = 0; i < m; ++i) { col[i] -= dot * Qq[i]; }
                }
            }
            double c1 = 0.0;
            for (std::size_t i = 0; i < m; ++i) { c1 += col[i] * col[i]; }
            if (!(c1 > rankTol * rankTol * c0)) { continue; }
            double const inv = 1.0 / std::sqrt(c1);
            L.Q.resize((L.nq + 1) * m);
            for (std::size_t i = 0; i < m; ++i) { L.Q[L.nq * m + i] = col[i] * inv; }
            ++L.nq;
        }
        L.rr = 0.0;
        for (auto z : L.r) { L.rr += z * z; }
        std::vector<double> e;
        gmPerp(L, L.r.data(), e);
        L.rPerp = 0.0;
        for (auto z : e) { L.rPerp += z * z; }
        if (jOut != nullptr) { *jOut = std::move(J); }
        return std::isfinite(L.rr) && std::isfinite(L.rPerp);
    }

    static void gmPerp(GmLin const& L, double const* v, std::vector<double>& e)
    {
        std::size_t const m = L.N + L.M;
        e.assign(v, v + m);
        for (std::size_t q = 0; q < L.nq; ++q) {
            double const* Qq = L.Q.data() + q * m;
            double dot = 0.0;
            for (std::size_t i = 0; i < m; ++i) { dot += Qq[i] * v[i]; }
            for (std::size_t i = 0; i < m; ++i) { e[i] -= dot * Qq[i]; }
        }
    }

    static void gmApplyJU(GmLin const& L, bool isX, double const* d, std::vector<double>& out)
    {
        std::size_t const N = L.N, M = L.M;
        out.assign(N + M, 0.0);
        std::size_t const iMax = isX ? L.iUmax : L.iVmax;
        std::size_t const iMin = isX ? L.iUmin : L.iVmin;
        double const dC = 0.5 * (d[iMax] + d[iMin]);
        double dS = 0.0;
        if (isX && L.sBranch == 1) { dS = d[L.iUmax] - d[L.iUmin]; }
        if (!isX && L.sBranch == 2) { dS = d[L.iVmax] - d[L.iVmin]; }
        thread_local std::vector<double> dqx, dqy;
        dqx.resize(M); dqy.resize(M);
        double const invS = L.invS;
        for (std::size_t j = 0; j < M; ++j) {
            if (isX) {
                dqx[j] = (d[j] - dC) * invS - L.ax[j] * dS * invS;
                dqy[j] = -L.ay[j] * dS * invS;
            } else {
                dqx[j] = -L.ax[j] * dS * invS;
                dqy[j] = (d[j] - dC) * invS - L.ay[j] * dS * invS;
            }
        }
        for (std::size_t i = 0; i < N; ++i) {
            if (L.sd[i] == 0.0) { continue; }
            std::size_t const j = static_cast<std::size_t>(L.jb[i]);
            double const t = L.tau[i];
            double const g0 = L.nhx[i] * dqx[j] + L.nhy[i] * dqy[j];
            double const g1 = L.nhx[i] * dqx[j + 1] + L.nhy[i] * dqy[j + 1];
            out[i] = -L.sd[i] * ((1.0 - t) * g0 + t * g1);
        }
        for (std::size_t j = 0; j < M; ++j) {
            if (L.sc[j] == 0.0) { continue; }
            out[N + j] = L.sc[j] * (L.mhx[j] * dqx[j] + L.mhy[j] * dqy[j]);
        }
    }

    static auto gmUnary(RegOp op, double x) -> double
    {
        switch (op) {
            case RegOp::Sin: return std::sin(x);
            case RegOp::Cos: return std::cos(x);
            case RegOp::Exp: return std::exp(x);
            case RegOp::Log: return std::log(x);
            case RegOp::Sqrt: return std::sqrt(x);
            case RegOp::Square: return x * x;
            case RegOp::Sub: return -x;
            case RegOp::Div: return 1.0 / x;
            default: return x;
        }
    }
    static auto gmBinary(RegOp op, double x1, double x2) -> double
    {
        switch (op) {
            case RegOp::Add: return x1 + x2;
            case RegOp::Sub: return x2 - x1;
            case RegOp::Mul: return x1 * x2;
            case RegOp::Div: return x2 / x1;
            default: return x1;
        }
    }
    static auto gmHash(RegOp op) -> Operon::Hash
    {
        using BO = Operon::BuiltinOp;
        switch (op) {
            case RegOp::Add: return static_cast<Operon::Hash>(BO::Add);
            case RegOp::Sub: return static_cast<Operon::Hash>(BO::Sub);
            case RegOp::Mul: return static_cast<Operon::Hash>(BO::Mul);
            case RegOp::Div: return static_cast<Operon::Hash>(BO::Div);
            case RegOp::Sin: return static_cast<Operon::Hash>(BO::Sin);
            case RegOp::Cos: return static_cast<Operon::Hash>(BO::Cos);
            case RegOp::Exp: return static_cast<Operon::Hash>(BO::Exp);
            case RegOp::Log: return static_cast<Operon::Hash>(BO::Log);
            case RegOp::Sqrt: return static_cast<Operon::Hash>(BO::Sqrt);
            default: return static_cast<Operon::Hash>(BO::Square);
        }
    }

    struct GmTape {
        RegisterProgram prog;
        std::size_t R{0}, M{0};
        std::vector<double> t, val, adj;
        std::vector<char> ok;
    };
    auto gmTape(Operon::Tree const& tree, Operon::Dataset const& ds, Operon::Range range, GmTape& T) const -> bool
    {
        auto const& nodes = tree.Nodes();
        for (auto const& nd : nodes) { if (nd.Arity > 2 || nd.Type == Operon::NodeType::Ref) { return false; } }
        auto c = tree.GetCoefficients();
        T.prog = RegisterProgram::Compile(tree, {c.data(), c.size()});
        if (!T.prog.valid || T.prog.insts.size() != nodes.size()) { return false; }
        auto tv = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
        std::size_t const M = tv.size();
        std::size_t const R = T.prog.insts.size();
        T.R = R; T.M = M;
        T.t.resize(M);
        for (std::size_t j = 0; j < M; ++j) { T.t[j] = static_cast<double>(tv[j]); }
        T.val.assign(R * M, 0.0);
        T.adj.assign(R * M, 0.0);
        for (std::size_t k = 0; k < R; ++k) {
            auto const& in = T.prog.insts[k];
            double const w = static_cast<double>(in.weight);
            double* o = T.val.data() + k * M;
            double const* v1 = T.val.data() + static_cast<std::size_t>(in.src1) * M;
            double const* v2 = T.val.data() + static_cast<std::size_t>(in.src2) * M;
            for (std::size_t j = 0; j < M; ++j) {
                switch (in.op) {
                    case RegOp::LoadConst: o[j] = static_cast<double>(in.imm); break;
                    case RegOp::LoadVar: o[j] = T.t[j] * w; break;
                    case RegOp::Identity: o[j] = v1[j] * w; break;
                    case RegOp::Add: case RegOp::Mul:
                        o[j] = gmBinary(in.op, v1[j], v2[j]) * w; break;
                    case RegOp::Sub: case RegOp::Div:
                        o[j] = (in.arity == 1 ? gmUnary(in.op, v1[j]) : gmBinary(in.op, v1[j], v2[j])) * w; break;
                    default: o[j] = gmUnary(in.op, v1[j]) * w; break;
                }
            }
        }
        for (std::size_t j = 0; j < M; ++j) { T.adj[(R - 1) * M + j] = 1.0; }
        for (std::size_t kk = R; kk-- > 0;) {
            auto const& in = T.prog.insts[kk];
            if (in.op == RegOp::LoadConst || in.op == RegOp::LoadVar) { continue; }
            double const w = static_cast<double>(in.weight);
            double const* a = T.adj.data() + kk * M;
            std::size_t const s1 = in.src1, s2 = in.src2;
            double const* v1 = T.val.data() + s1 * M;
            double const* v2 = T.val.data() + s2 * M;
            double* a1 = T.adj.data() + s1 * M;
            double* a2 = T.adj.data() + s2 * M;
            for (std::size_t j = 0; j < M; ++j) {
                double const g = a[j] * w;
                if (g == 0.0) { continue; }
                switch (in.op) {
                    case RegOp::Identity: a1[j] += g; break;
                    case RegOp::Add: a1[j] += g; a2[j] += g; break;
                    case RegOp::Sub:
                        if (in.arity == 1) { a1[j] -= g; } else { a2[j] += g; a1[j] -= g; }
                        break;
                    case RegOp::Mul: a1[j] += g * v2[j]; a2[j] += g * v1[j]; break;
                    case RegOp::Div:
                        if (in.arity == 1) { a1[j] -= g / (v1[j] * v1[j]); }
                        else { a2[j] += g / v1[j]; a1[j] -= g * v2[j] / (v1[j] * v1[j]); }
                        break;
                    case RegOp::Sin: a1[j] += g * std::cos(v1[j]); break;
                    case RegOp::Cos: a1[j] -= g * std::sin(v1[j]); break;
                    case RegOp::Exp: a1[j] += g * std::exp(v1[j]); break;
                    case RegOp::Log: a1[j] += g / v1[j]; break;
                    case RegOp::Sqrt: a1[j] += g * 0.5 / std::sqrt(v1[j]); break;
                    case RegOp::Square: a1[j] += 2.0 * g * v1[j]; break;
                    default: break;
                }
            }
        }
        T.ok.assign(R, 1);
        for (std::size_t k = 0; k < R; ++k) {
            for (std::size_t j = 0; j < M; ++j) {
                if (!std::isfinite(T.val[k * M + j]) || !std::isfinite(T.adj[k * M + j])) { T.ok[k] = 0; break; }
            }
        }
        return T.ok[R - 1] != 0;
    }

    auto gmLoss(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range) const -> double
    {
        std::vector<Operon::Scalar> gx, gy;
        interpInto(ind.x, ind.x.GetCoefficients(), ds, range, gx);
        interpInto(ind.y, ind.y.GetCoefficients(), ds, range, gy);
        if (gx.empty() || gx.size() != gy.size()) { return std::numeric_limits<double>::infinity(); }
        for (std::size_t j = 0; j < gx.size(); ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return std::numeric_limits<double>::infinity(); }
        }
        auto const sd = oneResBit(4)
            ? oneResLmLoss(gx, gy, dataX_, dataY_, dataBlocks_, nullptr)
            : shapeDistance(gx, gy, dataX_, dataY_, noise_, nullptr, 4, &dataBlocks_,
                                      config_.useGridSdf ? &gridSdf_ : nullptr,
                                      config_.useAmxChamfer ? &amxData_ : nullptr);
        return static_cast<double>(sd);
    }

    auto gmNatsAt(double delta, double level) const -> double
    {
        double const N = static_cast<double>(dataX_.size());
        if (!noise_.empty()) { return N * delta; }
        if (!(level > 0.0) || !std::isfinite(level)) { return N * delta; }
        double const nl = std::max(level + delta, level * 1e-300);
        return 0.5 * N * std::log(nl / level);
    }

    auto gmFitNats(double L) const -> double
    {
        double const N = static_cast<double>(dataX_.size());
        if (!noise_.empty()) { return N * L; }
        return 0.5 * N * std::log(std::max(L, 1e-300));
    }
    auto gmParamPerK() const -> double
    {
        return 0.5 * std::log(static_cast<double>(dataX_.size()) / (2.0 * Sketch::nml::kPi));
    }
    auto gmNmlOf(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range) const -> double
    {
        double const L = gmLoss(ind, ds, range);
        if (!std::isfinite(L)) { return std::numeric_limits<double>::infinity(); }
        GmLin lin;
        std::size_t kid = ind.x.GetCoefficients().size() + ind.y.GetCoefficients().size();
        if (gmBuildLin(ind, ds, range, lin)) { kid = lin.nq; }
        return gmFitNats(L)
            + Sketch::nml::structCodeCanonNats(static_cast<int>(ind.x.Length()))
            + Sketch::nml::structCodeCanonNats(static_cast<int>(ind.y.Length()))
            + static_cast<double>(kid) * gmParamPerK();
    }

    struct GmCand {
        int node{0};
        int kind{0};
        RegOp op{RegOp::Identity};
        int leaf{0};
        int child{-1};
        int newLen{0};
        int dK{0};
        double predFixed{0}, predRefit{0}, beta{0}, nats{0};
        double natsMut{0};
        double dStruct{0}, dParam{0};
        double pred{0};
        double loc{0};
    };

    auto gmApplyCand(SearchState const& s, Operon::Tree const& tree, GmCand const& c) const -> Operon::Tree
    {
        auto const& nodes = tree.Nodes();
        auto const k = static_cast<std::size_t>(c.node);
        std::size_t const start = k - nodes[k].Length;
        auto leafNode = [&](int leaf, double beta) -> Operon::Node {
            if (leaf == 2) {
                Operon::Node n(Operon::NodeType::Variable);
                n.HashValue = s.inputs[0];
                n.CalculatedHashValue = s.inputs[0];
                n.Value = static_cast<Operon::Scalar>(beta);
                return n;
            }
            return Operon::Node::Constant(beta);
        };
        std::vector<Operon::Node> out;
        out.reserve(nodes.size() + 2);
        out.insert(out.end(), nodes.begin(), nodes.begin() + static_cast<long>(start));
        switch (c.kind) {
            case 0: {
                out.insert(out.end(), nodes.begin() + static_cast<long>(start), nodes.begin() + static_cast<long>(k) + 1);
                out.back().Value = static_cast<Operon::Scalar>(static_cast<double>(out.back().Value) + c.beta);
                break;
            }
            case 1: {
                out.insert(out.end(), nodes.begin() + static_cast<long>(start), nodes.begin() + static_cast<long>(k));
                auto f = Operon::Node::Function(gmHash(c.op), nodes[k].Arity);
                f.Value = nodes[k].Value;
                out.push_back(f);
                break;
            }
            case 2: case 4: out.push_back(leafNode(c.leaf, c.beta)); break;
            case 3: {
                auto const ch = static_cast<std::size_t>(c.child);
                std::size_t const cs = ch - nodes[ch].Length;
                out.insert(out.end(), nodes.begin() + static_cast<long>(cs), nodes.begin() + static_cast<long>(ch) + 1);
                break;
            }
            case 5: {
                out.insert(out.end(), nodes.begin() + static_cast<long>(start), nodes.begin() + static_cast<long>(k) + 1);
                out.push_back(Operon::Node::Function(gmHash(c.op), 1));
                break;
            }
            default: {
                out.insert(out.end(), nodes.begin() + static_cast<long>(start), nodes.begin() + static_cast<long>(k) + 1);
                out.push_back(leafNode(c.leaf, c.beta));
                out.push_back(Operon::Node::Function(gmHash(c.op), 2));
                break;
            }
        }
        out.insert(out.end(), nodes.begin() + static_cast<long>(k) + 1, nodes.end());
        Operon::Tree t(out);
        t.UpdateNodes();
        return t;
    }

    static void gmApplyJUT(GmLin const& L, bool isX, double const* w, std::vector<double>& out)
    {
        std::size_t const N = L.N, M = L.M;
        thread_local std::vector<double> gx, gy;
        gx.assign(M, 0.0); gy.assign(M, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            if (L.sd[i] == 0.0) { continue; }
            std::size_t const j = static_cast<std::size_t>(L.jb[i]);
            double const t = L.tau[i];
            double const c = -L.sd[i] * w[i];
            gx[j] += c * (1.0 - t) * L.nhx[i]; gy[j] += c * (1.0 - t) * L.nhy[i];
            gx[j + 1] += c * t * L.nhx[i];     gy[j + 1] += c * t * L.nhy[i];
        }
        for (std::size_t j = 0; j < M; ++j) {
            if (L.sc[j] == 0.0) { continue; }
            double const c = L.sc[j] * w[N + j];
            gx[j] += c * L.mhx[j]; gy[j] += c * L.mhy[j];
        }
        out.assign(M, 0.0);
        double const invS = L.invS;
        auto const& g = isX ? gx : gy;
        double sumG = 0.0, sumS = 0.0;
        for (std::size_t j = 0; j < M; ++j) {
            out[j] += g[j] * invS;
            sumG += g[j];
            sumS += gx[j] * L.ax[j] + gy[j] * L.ay[j];
        }
        std::size_t const iMax = isX ? L.iUmax : L.iVmax;
        std::size_t const iMin = isX ? L.iUmin : L.iVmin;
        out[iMax] += -0.5 * invS * sumG;
        out[iMin] += -0.5 * invS * sumG;
        if ((isX && L.sBranch == 1) || (!isX && L.sBranch == 2)) {
            out[iMax] += -invS * sumS;
            out[iMin] -= -invS * sumS;
        }
    }

    struct GmEvent {
        int applied{0};
        int kind{-1};
        bool coin{false};
        std::size_t cands{0};
        double msLin{0}, msTape{0}, msCand{0}, msSel{0}, msLog{0};
        double predChosen{std::numeric_limits<double>::quiet_NaN()};
    };

    auto childTreeFiniteOnGrid(SearchState const& s, Operon::Tree const& t) const -> bool
    {
        thread_local std::vector<Operon::Scalar> buf;
        interpInto(t, t.GetCoefficients(), *s.activeDs, s.range, buf);
        for (auto v : buf) { if (!std::isfinite(static_cast<double>(v))) { return false; } }
        return true;
    }

    auto mutMultiValid(SearchState const& s, Operon::RandomGenerator& rng, Operon::Tree const& base,
                       PairIndividual const& parent, int* drawsOut = nullptr) const -> Operon::Tree
    {
        if (s.activeDs == nullptr) { return s.mut->multi(rng, base); }
        int const nOps = s.mut->insertIncluded ? 5 : 4;
        std::uniform_int_distribution<int> pick(0, nOps - 1);
        int const idx = pick(rng);
        int draws = 0;
        for (;;) {
            Operon::Tree child = gmArmMutate(s, rng, idx, base, parent);
            ++draws;
            bool const finite = childTreeFiniteOnGrid(s, child);
            if (op2DiagEnabled()) {
                auto& c = op2ValidCounters();
                c.draws[static_cast<std::size_t>(idx)].fetch_add(1, std::memory_order_relaxed);
                if (finite) { c.events[static_cast<std::size_t>(idx)].fetch_add(1, std::memory_order_relaxed); }
            }
            if (finite) { if (drawsOut != nullptr) { *drawsOut = draws; } return child; }
        }
    }
    auto mutRandomOp(SearchState const& s, Operon::RandomGenerator& rng, Operon::Tree t,
                     PairIndividual const& parent) const -> Operon::Tree
    {
        if (!op2(1)) { return s.mut->multi(rng, std::move(t)); }
        return mutMultiValid(s, rng, t, parent);
    }

    auto gmRandomMutate(SearchState const& s, Operon::RandomGenerator& rng, Operon::Tree t,
                        PairIndividual const& parent) const -> Operon::Tree
    {
        if (!config_.useNoiseGuidedMutation || parent.perPointResid_.empty()) { return mutRandomOp(s, rng, t, parent); }
        std::uniform_real_distribution<double> coin(0.0, 1.0);
        if (coin(rng) < 0.5) { return mutRandomOp(s, rng, t, parent); }
        auto* node = noiseGuidedLeafPick(rng, t, parent.perPointResid_);
        if (node == nullptr) { return mutRandomOp(s, rng, t, parent); }
        std::normal_distribution<Operon::Scalar> dist(Operon::Scalar{0}, Operon::Scalar{3});
        node->Value += dist(rng);
        return t;
    }

    auto gmErrModel(int kind) const -> GradMutErrModel const&
    {
        if (config_.gradientMutationInRunErr) { return gmIrModel_[kind < 0 ? 7 : static_cast<std::size_t>(kind)]; }
        return (kind < 0) ? config_.gradientMutationRandomArm : config_.gradientMutationErr[static_cast<std::size_t>(kind)];
    }
    auto gmErrKnown(int kind) const -> bool
    {
        return !config_.gradientMutationInRunErr || gmIrFit_[kind < 0 ? 7 : static_cast<std::size_t>(kind)] != 0;
    }
    auto gmG6eActive() const -> bool
    {
        static bool const envOn = [] { char const* e = std::getenv("SKETCH_G6E"); return e == nullptr || std::atoi(e) != 0; }();
        return config_.gradientMutationG6e && envOn && config_.gradientMutationMode == 1
            && config_.gradientMutationInRunErr;
    }
    auto gmPImprove(int kind, double pred) const -> double
    {
        auto const& e = gmErrModel(kind);
        if (!gmErrKnown(kind)) { return 0.5 * (1.0 - e.pi0); }
        double const a = std::abs(pred);
        double const loc = pred + e.mu + e.b * a;
        double const sc = std::max(e.s0 + e.s1 * a, std::numeric_limits<double>::min());
        return (1.0 - e.pi0) * gmstat::tCdf(-loc / sc, e.nu);
    }
