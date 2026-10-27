
    auto Finalize() -> MultitreeResult
    {
#ifdef __APPLE__
        pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
        if (!state_) {
            InitSearchState();
        }
        auto& s = *state_;
        if (sf2_ != 0 && std::getenv("SKETCH_SF2_DIAG") != nullptr) {
            std::fprintf(stderr, "[SF2END] mask=%d tour=%zu cross=%.3f annealRes=%.3f lowN=%zu lowM=%zu closed=%d k16 mean=%.1f min=%ld max=%ld last=%ld null=%.6g",
                         sf2_, config_.tournamentSize, static_cast<double>(config_.crossTreeProb), config_.annealResolution,
                         config_.annealLowN, config_.annealLowM, dataClosed_ ? 1 : 0,
                         sf2AdmN_ > 0 ? sf2AdmSum_ / static_cast<double>(sf2AdmN_) : -1.0, sf2AdmMin_, sf2AdmMax_, sf2Admitted_,
                         static_cast<double>(sf2NullLoss()));
            for (std::size_t L = 0; L < sf2LvTrials_.size(); ++L) {
                if (sf2LvTrials_[L] > 0) { std::fprintf(stderr, " L%zu=%.0f/%.0f", L, sf2LvSucc_[L], sf2LvTrials_[L]); }
            }
            std::fprintf(stderr, "\n");
        }

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

        if (config_.alpsLayers > 0 && !s.layers.empty()) {
            s.pop.clear();
            for (auto& lay : s.layers) {
                for (auto& ind : lay) { s.pop.push_back(std::move(ind)); }
            }
            rankAndCrowd(s.pop);
        }

        if (s.lowRes) { switchToFull(s, &s.pop, nullptr, s.lastGen + 1); }

        errorClampValue() = config_.floorClamp;

        auto* best = bestAccuracy(s.pop);
        PairIndividual prePolish = *best;
        Operon::Range const fullRange{0, s.M};
        optConstsGuarded(*best, *s.fullDs, fullRange, config_.polishFullObjective);
        evaluateGuarded(*best, *s.fullDs, fullRange);
        if (polishWorsened(best->shapeErr, prePolish.shapeErr)) { *best = std::move(prePolish); }
        rankAndCrowd(s.pop);

        if (config_.useArchive) {
            std::vector<PairIndividual> polished;
            polished.push_back(*best);
            archiveUpdate(s.archive, polished, s.archiveCap);
        }
        if (config_.useComplexityHallOfFame) {
            std::vector<PairIndividual> polished;
            polished.push_back(*best);
            hallOfFameUpdate(s.hallOfFame, polished, kHallOfFameK);
        }

        std::vector<PairIndividual*> f0;
        if (config_.useArchive) {
            for (auto& ind : s.archive) { f0.push_back(&ind); }
        } else {
            for (auto& ind : s.pop) { if (ind.rank == 0) { f0.push_back(&ind); } }
        }

        MultitreeResult res;
        res.front = buildFront(f0, *s.fullDs, fullRange);

        if (config_.useComplexityHallOfFame && !s.hallOfFame.empty()) {
            std::vector<PairIndividual*> hofPtrs;
            hofPtrs.reserve(s.hallOfFame.size());
            for (auto& ind : s.hallOfFame) { hofPtrs.push_back(&ind); }
            res.hallOfFame = buildFront(hofPtrs, *s.fullDs, fullRange);
        }

        PairIndividual* ref = nullptr;
        if (f0.empty()) {
            ref = bestAccuracy(s.pop);
        } else if (config_.legacyFrontRecommendation) {
            ref = f0.back();
        } else {
            ref = siCisRecommend(f0, *s.fullDs, fullRange);
        }
        res.x = ref->x; res.y = ref->y;
        res.shapeError = ref->shapeErr;
        if (config_.computeStrictGateTelemetry) {
            res.siCisRecommendedShapeNll = siCisRecommendationScore(*ref, *s.fullDs, fullRange);
            PairIndividual lineNull = buildLineNull(s);
            optConstsGuarded(lineNull, *s.fullDs, fullRange, config_.polishFullObjective);
            res.nullLineStrictShapeNll = siCisRecommendationScore(lineNull, *s.fullDs, fullRange);
        }

        if (config_.useChowAbstention) {
            constexpr double kLog2E = 1.4426950408889634;
            double best = std::numeric_limits<double>::infinity();
            double second = std::numeric_limits<double>::infinity();
            for (auto* p : f0) {
                double const sc = static_cast<double>(siCisRecommendationScore(*p, *s.fullDs, fullRange));
                if (!std::isfinite(sc)) { continue; }
                if (sc < best) { second = best; best = sc; }
                else if (sc < second) { second = sc; }
            }
            res.chowAbstainThresholdBits =
                static_cast<double>(Sketch::chowAbstainThresholdBits(config_.chowAbstentionCost));
            if (std::isfinite(best) && std::isfinite(second)) {
                double const leadNats = second - best;
                double const leadBits = (leadNats > 0.0 ? leadNats : 0.0) * kLog2E;
                res.chowAbstainLeadBits = leadBits;
                res.chowShouldAnswer =
                    Sketch::chowShouldAnswer(static_cast<float>(leadBits), config_.chowAbstentionCost);
            } else {
                res.chowAbstainLeadBits = std::numeric_limits<double>::infinity();
                res.chowShouldAnswer = true;
            }
        }

        int const t13Mask = Sketch::sketchT13Mask();
        bool const t13Bit0 = Sketch::sketchT13On(t13Mask, 0);
        bool const t13Bit1 = Sketch::sketchT13On(t13Mask, 1);
        if ((config_.reportNmlPick || t13Bit0 || t13Bit1) && !f0.empty()) {
            auto const tNml0 = std::chrono::steady_clock::now();
            auto ranges = nmlRanges();
            ranges.gridPeriod = 2.0 * Sketch::nml::kPi * static_cast<double>(s.M > 1 ? s.M - 1 : 1);
            res.nmlFreqRange = ranges.freq;
            res.nmlExpRateRange = ranges.expRate;
            res.nmlExtent = ranges.extent;
            res.nmlGridPeriod = ranges.gridPeriod;
            res.nmlSigmaBar = ranges.sigmaBar;
            res.nmlFront.reserve(f0.size());
            struct CounterSnap {
                std::size_t d, l, it, ee, nf, lg, sq, dv, ex, none, rec, rej;
            } const snap{nConstDescents_.load(), nConstLossCalls_.load(), nConstIters_.load(), nEvalErrors_.load(),
                         nNonFiniteFitness_.load(), nNonFiniteWithLog_.load(), nNonFiniteWithSqrt_.load(),
                         nNonFiniteWithDiv_.load(), nNonFiniteWithExp_.load(), nNonFiniteWithNoneOfThese_.load(),
                         nLmRecoveredFromNonFiniteStart_.load(), nRejectedNonFiniteBeforeConstOpt_.load()};
            double bestTotal = std::numeric_limits<double>::infinity();
            double bestTotalAbs = std::numeric_limits<double>::infinity();
            bool const dedupeOn = Sketch::nmlDedupeEnabledFromEnv();
            std::unordered_map<std::string, std::size_t> dedupeIndex;
            if (dedupeOn) { dedupeIndex.reserve(f0.size()); }
            std::size_t nDistinct = 0;
            for (std::size_t i = 0; i < f0.size(); ++i) {
                std::size_t reuse = std::string::npos;
                std::string key;
                if (dedupeOn) {
                    auto const cxr = f0[i]->x.GetCoefficients();
                    auto const cyr = f0[i]->y.GetCoefficients();
                    key.reserve(res.front[i].exprX.size() + res.front[i].exprY.size()
                                + (cxr.size() + cyr.size()) * sizeof(Operon::Scalar) + 2);
                    key = res.front[i].exprX;
                    key += '\x1e';
                    key += res.front[i].exprY;
                    key += '\x1e';
                    key.append(reinterpret_cast<char const*>(cxr.data()), cxr.size() * sizeof(Operon::Scalar));
                    key.append(reinterpret_cast<char const*>(cyr.data()), cyr.size() * sizeof(Operon::Scalar));
                    auto const it = dedupeIndex.find(key);
                    if (it != dedupeIndex.end()) { reuse = it->second; }
                }
                if (reuse != std::string::npos) {
                    res.nmlFront.push_back(res.nmlFront[reuse]);
                    if (config_.nmlPolishMembers) {
                        res.nmlFrontRaw.push_back(res.nmlFrontRaw[reuse]);
                        res.nmlPolishedExprX.push_back(res.nmlPolishedExprX[reuse]);
                        res.nmlPolishedExprY.push_back(res.nmlPolishedExprY[reuse]);
                    }
                } else {
                    if (dedupeOn) { dedupeIndex.emplace(std::move(key), i); }
                    ++nDistinct;
                    auto raw = nmlTermsFor(*f0[i], *s.fullDs, fullRange, ranges);
                    if (config_.nmlPolishMembers) {
                        auto const tp0 = std::chrono::steady_clock::now();
                        PairIndividual cp = *f0[i];
                        evaluateGuarded(cp, *s.fullDs, fullRange);
                        PairIndividual pre = cp;
                        optConstsGuarded(cp, *s.fullDs, fullRange, config_.polishFullObjective);
                        evaluateGuarded(cp, *s.fullDs, fullRange);
                        bool const rolled = polishWorsened(cp.shapeErr, pre.shapeErr);
                        double const after = static_cast<double>(cp.shapeErr);
                        double const before = static_cast<double>(pre.shapeErr);
                        auto const nllNoRb = siCisRecommendationScore(cp, *s.fullDs, fullRange);
                        if (rolled) { cp = std::move(pre); }
                        double const pms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tp0).count();
                        res.nmlPolishMs += pms;
                        auto pol = nmlTermsFor(cp, *s.fullDs, fullRange, ranges);
                        pol.polished = true;
                        pol.polishRolledBack = rolled;
                        pol.polishMs = pms;
                        pol.shapeErrBefore = before;
                        pol.nllPolishNoRollback = (std::isfinite(nllNoRb) && nllNoRb < std::numeric_limits<Operon::Scalar>::max())
                            ? static_cast<double>(nllNoRb) : std::numeric_limits<double>::quiet_NaN();
                        pol.shapeErrAfter = after;
                        res.nmlPolishedExprX.push_back(Operon::InfixFormatter::Format(cp.x, *s.fullDs, 4));
                        res.nmlPolishedExprY.push_back(Operon::InfixFormatter::Format(cp.y, *s.fullDs, 4));
                        res.nmlFrontRaw.push_back(std::move(raw));
                        res.nmlFront.push_back(std::move(pol));
                    } else {
                        res.nmlFront.push_back(std::move(raw));
                    }
                }
                auto const& t = res.nmlFront.back();
                if (t.ok && t.total < bestTotal) { bestTotal = t.total; res.nmlPickIndex = static_cast<int>(i); }
                if (t.absOk && t.absTotal < bestTotalAbs) { bestTotalAbs = t.absTotal; res.nmlPickIndexAbs = static_cast<int>(i); }
                if (f0[i] == ref) { res.nmlRecommendedIndex = static_cast<int>(i); }
                if (t13Bit0) {
                    res.front[i].nmlAbsTotal = t.absOk ? t.absTotal : std::numeric_limits<double>::quiet_NaN();
                }
            }
            res.nmlNMembers = f0.size();
            res.nmlNDistinct = nDistinct;
            nConstDescents_ = snap.d; nConstLossCalls_ = snap.l; nConstIters_ = snap.it; nEvalErrors_ = snap.ee;
            nNonFiniteFitness_ = snap.nf; nNonFiniteWithLog_ = snap.lg; nNonFiniteWithSqrt_ = snap.sq;
            nNonFiniteWithDiv_ = snap.dv; nNonFiniteWithExp_ = snap.ex; nNonFiniteWithNoneOfThese_ = snap.none;
            nLmRecoveredFromNonFiniteStart_ = snap.rec; nRejectedNonFiniteBeforeConstOpt_ = snap.rej;
            res.nmlComputeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tNml0).count();
            if (t13Bit1 && res.nmlPickIndexAbs >= 0 && static_cast<std::size_t>(res.nmlPickIndexAbs) < f0.size()) {
                auto* recAbs = f0[static_cast<std::size_t>(res.nmlPickIndexAbs)];
                res.x = recAbs->x; res.y = recAbs->y; res.shapeError = recAbs->shapeErr;
            }
        }

        res.dataTurning = dataTurningF_;
        res.rotationGaugeAngle = rotationGaugeAngle_;
        res.parametrizationGauge = parametrizationGaugeRecord_;
        res.parametrizationPerm = parametrizationPerm_;
        res.generations = s.lastGen + 1;
        res.truncated = s.truncated;
        res.nEvalErrors = nEvalErrors_;
        res.nConstDescents = nConstDescents_;
        res.nConstLossCalls = nConstLossCalls_;
        res.nConstIters = nConstIters_;
        res.nLmStopUnresolvableDamped = nLmStopUnresolvableDamped_;
        res.nLmStopRealizedSmall = nLmStopRealizedSmall_;
        res.nLmStopFixedPoint = nLmStopFixedPoint_;
        if (config_.useConcurrentCpuGpuLM) {
            res.cclmStats = cclmStats_;
            res.cclmParity = cclmParity_;
            res.cclmGpuAllocEvents = cclmGpu_ ? cclmAllocEvents(*cclmGpu_) : 0;
        }
        res.nNonFiniteFitness = nNonFiniteFitness_;
        res.nMutationOversizeRetries = nMutationOversizeRetries_;
        res.nMutationFallbackToParent = nMutationFallbackToParent_;
        res.nRejectedNonFiniteBeforeConstOpt = nRejectedNonFiniteBeforeConstOpt_;
        res.nNonFiniteWithLog = nNonFiniteWithLog_;
        res.nNonFiniteWithSqrt = nNonFiniteWithSqrt_;
        res.nNonFiniteWithDiv = nNonFiniteWithDiv_;
        res.nNonFiniteWithExp = nNonFiniteWithExp_;
        res.nNonFiniteWithNoneOfThese = nNonFiniteWithNoneOfThese_;
        res.nLmRecoveredFromNonFiniteStart = nLmRecoveredFromNonFiniteStart_;
        res.nCascadeCullGens = nCascadeCullGens_;
        res.nCascadeCullSeen = nCascadeCullSeen_;
        res.nCascadeCullSurvived = nCascadeCullSurvived_;
        res.nLiveCoarseSieveSeen = nLiveCoarseSieveSeen_;
        res.nLiveCoarseSieveCulled = nLiveCoarseSieveCulled_;
        auto const nLiveFinite = nLiveCoarseSieveFinite_.load(std::memory_order_relaxed);
        res.liveCoarseSieveBoundMin = liveCoarseSieveBoundMin_.load(std::memory_order_relaxed);
        res.liveCoarseSieveBoundMax = liveCoarseSieveBoundMax_.load(std::memory_order_relaxed);
        if (nLiveFinite > 0) {
            res.liveCoarseSieveBoundMean = liveCoarseSieveBoundSum_.load(std::memory_order_relaxed) /
                static_cast<double>(nLiveFinite);
        }
        res.simplifyGensRan = s.simplifyGensRan;
        res.annealSwitchGen = s.annealSwitchGen;
        res.precomputeTimeMs = precomputeTimeMs_;
        if (config_.profileTiming) {
            res.profEvalMs = static_cast<double>(profEvalNs_.load()) / 1e6;
            res.profLmMs = static_cast<double>(profLmNs_.load()) / 1e6;
            res.profScoreMs = static_cast<double>(profScoreNs_.load()) / 1e6;
            res.profBreedMs = static_cast<double>(profBreedNs_.load()) / 1e6;
            res.profSelectMs = static_cast<double>(profSelectNs_.load()) / 1e6;
            res.profSelectCombineMs = static_cast<double>(profSelectCombineNs_.load()) / 1e6;
            res.profSelectNdsMs = static_cast<double>(profSelectNdsNs_.load()) / 1e6;
            res.profSelectCrowdMs = static_cast<double>(profSelectCrowdNs_.load()) / 1e6;
            res.profGenMs = static_cast<double>(profGenNs_.load()) / 1e6;
            res.profScoreDistMs = static_cast<double>(profScoreDistNs_.load()) / 1e6;
            res.profScoreCisMs = static_cast<double>(profScoreCisNs_.load()) / 1e6;
            res.profScoreOrderMs = static_cast<double>(profScoreOrderNs_.load()) / 1e6;
            res.profScoreCurvMs = static_cast<double>(profScoreCurvNs_.load()) / 1e6;
            res.profLmBuildMs = static_cast<double>(profLmBuildNs_.load()) / 1e6;
            res.profLmLossMs = static_cast<double>(profLmLossNs_.load()) / 1e6;
            res.profCisPointSegMs = static_cast<double>(profCisPointSegNs_.load()) / 1e6;
            res.profCisNllMs = static_cast<double>(profCisNllNs_.load()) / 1e6;
            res.profCisReverseMs = static_cast<double>(profCisReverseNs_.load()) / 1e6;
            res.profCisGrid48Ms = static_cast<double>(profCisGrid48Ns_.load()) / 1e6;
            res.profGpuPrepassMs = static_cast<double>(profGpuPrepassNs_.load()) / 1e6;
            res.profParallelWallMs = static_cast<double>(profParallelWallNs_.load()) / 1e6;
        }
        if (config_.verbose && (nEvalErrors_ > 0 || nInvalidChildren_ > 0)) {
            std::cout << "[robustness] eval exceptions absorbed: " << nEvalErrors_
                      << ", corrupted children discarded: " << nInvalidChildren_ << "\n";
        }
        if (config_.verbose && config_.semanticCrossover) {
            std::cout << "[ssc] semantic crossovers: " << nSemCross_
                      << ", fell back to random subtree crossover: " << nSemFallback_ << "\n";
        }
#ifdef SKETCH_DUP_STATS
        dupStats_.print("multitree");
#endif
        return res;
    }

    struct GradMutTrajPoint { int gen{0}; double wallMs{0}; double bestNml{0}; int N{0}; };
    struct GradMutStats {
        std::size_t events{0};
        std::size_t applied{0};
        std::size_t randomModel{0};
        std::size_t fallback{0};
        std::size_t banditRandomArm{0};
        std::size_t candidates{0};
        std::size_t kindApplied[7]{};
        double msLin{0}, msTape{0}, msCand{0}, msSel{0}, msLog{0};
        double passWallMs{0};
        double banditMs{0};
        std::size_t armPulls[7]{};
        std::size_t armObs[7]{};
        double armMean[7]{};
        std::size_t irCoin{0};
        std::size_t irStubs{0};
        std::size_t irObs{0};
        std::size_t irDup{0};
        std::size_t irN[8]{}, irN0[8]{};
        bool irFit[8]{};
        GradMutErrModel irModel[8]{};
        double irMs{0};
        double irFitMs{0};
        std::vector<GradMutTrajPoint> traj;
    };
    [[nodiscard]] auto gradMutStats() const -> GradMutStats const& { return gmStats_; }
    [[nodiscard]] auto mutOpDiag() const -> MutOpDiagCounters const& { return mutDiag_; }
    [[nodiscard]] auto gradMutEff() const -> std::vector<GradMutEff> const& { return gmEff_; }
    [[nodiscard]] auto gradMutLog() const -> std::vector<GradMutRecord> const& { return gmLog_; }

    auto Run() -> MultitreeResult
    {
        if (inputRejected_) {
            return MultitreeResult{};
        }
        InitSearchState();
        RunTo(config_.generations);
        return Finalize();
    }

    double BestShapeError() const
    {
        if (!state_) { return std::numeric_limits<double>::infinity(); }
        auto const& s = *state_;
        if (config_.alpsLayers > 0 && !s.layers.empty()) {
            double best = std::numeric_limits<double>::infinity();
            for (auto const& lay : s.layers) {
                auto const* b = bestAccuracy(lay);
                if (b != nullptr && b->shapeErr < best) {
                    best = b->shapeErr;
                }
            }
            return best;
        }
        if (config_.useArchive && !s.archive.empty()) {
            auto const* b = bestAccuracy(s.archive);
            return b != nullptr ? static_cast<double>(b->shapeErr) : std::numeric_limits<double>::infinity();
        }
        if (!s.pop.empty()) {
            auto const* b = bestAccuracy(s.pop);
            return b != nullptr ? static_cast<double>(b->shapeErr) : std::numeric_limits<double>::infinity();
        }
        return std::numeric_limits<double>::infinity();
    }

private:
    auto buildEngineAmplitudeContext(PairIndividual const& candidate,
                                     std::vector<PairAmplitudeNodeSlot> const& slots) const
        -> std::optional<BracketNewtonRuntimeContext> {
        if (!state_ || state_->activeDs == nullptr || slots.empty() ||
            state_->range.Size() == 0 || state_->range.Size() > dataX_.size() ||
            state_->range.Size() > dataY_.size()) {
            return std::nullopt;
        }
        auto const tValues = state_->activeDs->GetValues(static_cast<int64_t>(0))
                                 .subspan(state_->range.Start(), state_->range.Size());
        std::vector<double> sampleT(tValues.begin(), tValues.end());
        std::vector<double> targetX(dataX_.begin(), dataX_.begin() + state_->range.Size());
        std::vector<double> targetY(dataY_.begin(), dataY_.begin() + state_->range.Size());
        auto evaluate = [&](PairIndividual const& trial) {
            std::vector<Operon::Scalar> coeffs = trial.x.GetCoefficients();
            std::vector<Operon::Scalar> values(state_->range.Size());
            try {
                interpInto(trial.x, coeffs, *state_->activeDs, state_->range, values);
            } catch (...) {
                values.clear();
            }
            return std::vector<double>(values.begin(), values.end());
        };
        auto context = buildPerturbedAmplitudeContext(candidate, slots, std::move(sampleT),
                                                      std::move(targetX), std::move(targetY), evaluate);
        if (!context.has_value()) return std::nullopt;
        std::vector<Operon::Scalar> yCoeffs = candidate.y.GetCoefficients();
        std::vector<Operon::Scalar> yValues(state_->range.Size());
        try {
            interpInto(candidate.y, yCoeffs, *state_->activeDs, state_->range, yValues);
        } catch (...) {
            return std::nullopt;
        }
        context->candidateY.assign(yValues.begin(), yValues.end());
        if (!validRuntimeContext(*context, slots.size())) return std::nullopt;
        return context;
    }

    auto originalArcIndex(std::size_t idx) const -> std::size_t {
        if ((config_.useParametrizationGauge || config_.useQuotientProjection)
            && !parametrizationPerm_.empty() && parametrizationPerm_.size() == dataX_.size()
            && idx < parametrizationPerm_.size()) {
            return static_cast<std::size_t>(parametrizationPerm_[idx]);
        }
        return idx;
    }

    auto newtonPointSigma2ForRange(Operon::Range range) const -> std::vector<double> {
        std::size_t const N = dataX_.size();
        std::size_t const start = range.Start();
        std::size_t const count = range.Size();
        if (N == 0U || count == 0U || start + count > N) { return {}; }
        std::vector<double> s2(count, 0.0);
        if (config_.useExternalPointSigma2 && !externalPointSigma2_.empty()) {
            auto const& rlut = noiseRegionLUT(N, externalPointSigma2_.size());
            for (std::size_t i = 0; i < count; ++i) {
                s2[i] = static_cast<double>(externalPointSigma2_[rlut[originalArcIndex(start + i)]]);
            }
        } else if (noise_.size() > 0U) {
            auto const& rlut = noiseRegionLUT(N, noise_.size());
            for (std::size_t i = 0; i < count; ++i) {
                s2[i] = static_cast<double>(noise_.s2[rlut[originalArcIndex(start + i)]]);
            }
        } else {
            return {};
        }
        bool constant = true;
        for (std::size_t i = 1; i < count; ++i) {
            if (s2[i] != s2[0]) { constant = false; break; }
        }
        if (constant) { return {}; }
        for (double v : s2) { if (!std::isfinite(v) || v <= 0.0) { return {}; } }
        return s2;
    }

    DTable const* dtable_;
#if defined(__APPLE__)
    std::optional<DTable> accelDTable_;
#endif
    CurveBlocks dataBlocks_;
    GridSdf gridSdf_;
    GridSdf fullGridSdf_;
    AmxData amxData_;
    AmxData fullAmxData_;
    AmxData midAmxData_;
    double precomputeTimeMs_{0.0};
    bool inputRejected_{false};
    std::vector<Operon::Scalar> dataX_;
    std::vector<Operon::Scalar> dataY_;
    std::vector<Operon::Scalar> dataArcX_;
    std::vector<Operon::Scalar> dataArcY_;
    std::vector<Vec2> cisGrid48Stroke_;
    bool cisGrid48Ok_{false};
    std::vector<Vec2> cisGrid48StrokeData_;
    double cisGrid48Scale_{1.0};
    std::vector<double> cisGrid48Zi5LnE_;
    std::vector<double> cisGrid48Zi5Sd_;
    std::vector<std::vector<double>> cisGrid48Zi5MemberLnE_;
    bool cisGrid48Zi5Ok_{false};
    std::vector<Operon::Scalar> dataXDs_;
    std::vector<Operon::Scalar> dataYDs_;
    CurveBlocks dataBlocksDs_;
    std::vector<std::size_t> dsOrigIdx_;
    CurveBlocks fullBlocks_;
    std::vector<Operon::Scalar> fullX_;
    std::vector<Operon::Scalar> fullY_;
    bool annealActive_{false};
    CurveBlocks midBlocks_;
    std::vector<Operon::Scalar> midX_;
    std::vector<Operon::Scalar> midY_;
    bool annealMidActive_{false};
    double turningNumber_ = 0.0;
    double dataTurningF_ = 0.0;
    double curvatureFloor_ = 0.0;
    bool dataClosed_{false};
    mutable float topoClosedness_ = -1.0f;
    double rotationGaugeAngle_ = 0.0;
    ParametrizationGauge parametrizationGaugeRecord_{};
    std::vector<int> parametrizationPerm_{};
    struct QuotientFiniteMinBundle {
        std::vector<Operon::Scalar> x;
        std::vector<Operon::Scalar> y;
    };
    std::vector<QuotientFiniteMinBundle> quotientFiniteMinBundles_{};
    static constexpr double kQfmEigenGapRelTol = 1e-2;
    static constexpr double kQfmThetaTol = 0.05;
    static constexpr double kQfmRadiusRelTol = 1e-6;
    static constexpr int kQfmMaxCyclicAlts = 4;
    static constexpr std::size_t kQfmMaxBundles = 16;

    int sf2_{sf2Mask()};
    auto sf2(int bit) const -> bool { return sf2On(sf2_, bit); }
    int op2_{op2Mask()};
    auto op2(int bit) const -> bool { return op2On(op2_, bit); }
    mutable double sf2Sigma2_{-1.0};
    auto sf2NoiseVar() const -> double
    {
        if (sf2Sigma2_ >= 0.0) { return sf2Sigma2_; }
        double v = 0.0;
        if (!noise_.empty()) {
            double acc = 0.0;
            for (auto e : noise_.s2) { acc += static_cast<double>(e); }
            v = acc / static_cast<double>(noise_.size());
        } else if (dataX_.size() >= 3) {
            std::vector<double> e;
            e.reserve(2 * dataX_.size());
            double const inv = 1.0 / std::sqrt(6.0);
            for (std::size_t j = 1; j + 1 < dataX_.size(); ++j) {
                e.push_back(std::abs(static_cast<double>(dataX_[j - 1]) - 2.0 * dataX_[j] + dataX_[j + 1]) * inv);
                e.push_back(std::abs(static_cast<double>(dataY_[j - 1]) - 2.0 * dataY_[j] + dataY_[j + 1]) * inv);
            }
            auto mid = e.begin() + static_cast<std::ptrdiff_t>(e.size() / 2);
            std::nth_element(e.begin(), mid, e.end());
            double const sd = 1.4826 * *mid;
            v = sd * sd;
        }
        sf2Sigma2_ = (std::isfinite(v) && v > 0.0) ? v : 0.0;
        return sf2Sigma2_;
    }
    mutable std::size_t sf2KMax_{0};
    auto sf2NyquistCount() const -> std::size_t
    {
        auto const n = dataX_.size();
        if (n < 4) { return n; }
        double const s2 = sf2NoiseVar();
        double const thr = 2.0 * s2 * std::log(static_cast<double>(n));
        std::vector<double> zr(n), zi(n);
        double const x0 = dataX_.front(), y0 = dataY_.front();
        double const dxc = static_cast<double>(dataX_.back()) - x0, dyc = static_cast<double>(dataY_.back()) - y0;
        for (std::size_t j = 0; j < n; ++j) {
            double const f = static_cast<double>(j) / static_cast<double>(n - 1);
            zr[j] = static_cast<double>(dataX_[j]) - (x0 + f * dxc);
            zi[j] = static_cast<double>(dataY_[j]) - (y0 + f * dyc);
        }
        double const norm = 1.0 / std::sqrt(static_cast<double>(n));
        std::size_t kMax = 0;
        long const half = static_cast<long>(n / 2);
        for (long k = -half; k <= half; ++k) {
            if (k == 0) { continue; }
            double re = 0.0, im = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                double const a = -2.0 * M_PI * static_cast<double>(k) * static_cast<double>(j) / static_cast<double>(n);
                double const c = std::cos(a), sn = std::sin(a);
                re += zr[j] * c - zi[j] * sn;
                im += zr[j] * sn + zi[j] * c;
            }
            double const p = (re * re + im * im) * norm * norm;
            if (p > thr) { kMax = std::max<std::size_t>(kMax, static_cast<std::size_t>(std::labs(k))); }
        }
        sf2KMax_ = kMax;
        return 2 * kMax + 1;
    }
    void sf2ApplyBudgetRules()
    {
        if (sf2_ == 0) { return; }
        double const P = static_cast<double>(std::max<std::size_t>(3, config_.populationSize));
        double const G = static_cast<double>(std::max<std::size_t>(1, config_.generations));
        if (sf2(0)) {
            double const lnP = std::log(P);
            double const sOpt = std::exp((lnP + std::log(lnP)) / G);
            config_.tournamentSize = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(sOpt)));
        }
        if (sf2(1) && config_.crossTreeProb > Operon::Scalar{0}) {
            config_.crossTreeProb = Operon::Scalar{0.5};
        }
        if (sf2(4) && config_.annealResolution > 0.0) {
            double const rho = 1.0 - std::ceil(std::log2(P)) / G;
            config_.annealResolution = std::max(0.0, rho);
        }
    }
    static auto sf2OrderNats(std::vector<Operon::Scalar> const& tau) -> Operon::Scalar
    {
        auto const n = tau.size();
        if (n < 2) { return Operon::Scalar{0}; }
        thread_local std::vector<double> st;
        st.resize(n - 1);
        double lo = tau[0], hi = tau[0];
        for (std::size_t i = 0; i + 1 < n; ++i) {
            st[i] = std::abs(static_cast<double>(tau[i + 1]) - static_cast<double>(tau[i]));
            lo = std::min(lo, static_cast<double>(tau[i + 1]));
            hi = std::max(hi, static_cast<double>(tau[i + 1]));
        }
        auto mid = st.begin() + static_cast<std::ptrdiff_t>(st.size() / 2);
        std::nth_element(st.begin(), mid, st.end());
        double const sFloor = static_cast<double>(std::numeric_limits<float>::epsilon()) * std::max(1.0, hi - lo);
        double const sc = std::max(*mid, sFloor);
        double acc = 0.0;
        for (std::size_t i = 0; i + 1 < n; ++i) {
            double const b = static_cast<double>(tau[i]) - static_cast<double>(tau[i + 1]);
            if (b > 0.0) { acc += b * b; }
        }
        double const v = acc / (2.0 * sc * sc);
        return std::isfinite(v) ? static_cast<Operon::Scalar>(std::min(v, static_cast<double>(std::numeric_limits<Operon::Scalar>::max())))
                                : std::numeric_limits<Operon::Scalar>::max();
    }
    auto sf2OrderValue(std::vector<Operon::Scalar> const& tau) const -> Operon::Scalar
    {
        if (sf2(2)) { return sf2OrderNats(tau); }
        return monotonicityPenalty(tau, monotonicityEpsJitter(tau));
    }
    auto sf2NullLoss() const -> Operon::Scalar
    {
        std::vector<Operon::Scalar> const cx{Operon::Scalar{0}, Operon::Scalar{0}}, cy{Operon::Scalar{0}, Operon::Scalar{0}};
        auto const v = shapeDistance(cx, cy, dataX_, dataY_, noise_, nullptr, 4, &dataBlocks_,
                                     config_.useGridSdf ? &gridSdf_ : nullptr, config_.useAmxChamfer ? &amxData_ : nullptr);
        return std::isfinite(v) ? v : std::numeric_limits<Operon::Scalar>::infinity();
    }
    void sf2LmRules(LevenbergMarquardt& lm) const
    {
        if (sf2(19)) { lm.enableScreening = false; }
        if (sf2(20) || sf2(22)) {
            lm.maxIter = std::numeric_limits<int>::max();
            lm.maxTrials = static_cast<int>(std::ceil(std::log(lm.lambdaMax / lm.lambdaMin) / std::log(lm.lambdaUp)));
        }
        if (sf2(21)) { lm.floorThreshold = 0.0; }
        if (sf2(22)) { lm.resolutionStop = true; lm.dataCount = static_cast<double>(dataX_.size()); }
        if (op2(2)) { lm.gaugeProject = true; lm.jacobianEps = static_cast<double>(std::numeric_limits<float>::epsilon()); }
    }
    static auto sf2TreeKey(Operon::Tree const& t, std::uint64_t h) -> std::uint64_t
    {
        for (auto const& nd : t.Nodes()) {
            std::uint64_t v = 0;
            auto const val = static_cast<double>(nd.Value);
            std::memcpy(&v, &val, sizeof v);
            for (std::uint64_t w : {static_cast<std::uint64_t>(nd.HashValue), v, static_cast<std::uint64_t>(nd.Arity)}) {
                h ^= w; h *= 1099511628211ULL;
            }
        }
        return h;
    }
    static auto sf2SameTree(Operon::Tree const& a, Operon::Tree const& b) -> bool
    {
        auto const& na = a.Nodes(); auto const& nb = b.Nodes();
        if (na.size() != nb.size()) { return false; }
        for (std::size_t i = 0; i < na.size(); ++i) {
            if (na[i].HashValue != nb[i].HashValue || na[i].Arity != nb[i].Arity) { return false; }
            if (std::memcmp(&na[i].Value, &nb[i].Value, sizeof na[i].Value) != 0) { return false; }
        }
        return true;
    }
    static void sf2DedupeOrder(std::vector<PairIndividual> const& off, std::vector<std::size_t>& order,
                               std::vector<std::pair<std::size_t, std::size_t>>& dups)
    {
        std::unordered_map<std::uint64_t, std::vector<std::size_t>> reps;
        std::vector<std::size_t> kept;
        kept.reserve(order.size());
        for (auto i : order) {
            auto const key = sf2TreeKey(off[i].y, sf2TreeKey(off[i].x, 1469598103934665603ULL));
            auto& bucket = reps[key];
            std::size_t rep = off.size();
            for (auto r : bucket) {
                if (sf2SameTree(off[r].x, off[i].x) && sf2SameTree(off[r].y, off[i].y)) { rep = r; break; }
            }
            if (rep == off.size()) { bucket.push_back(i); kept.push_back(i); }
            else { dups.emplace_back(i, rep); }
        }
        order.swap(kept);
    }
    mutable std::vector<double> sf2LvSucc_, sf2LvTrials_;
    mutable long sf2Admitted_{-1};
    mutable double sf2AdmSum_{0.0};
    mutable long sf2AdmN_{0}, sf2AdmMin_{-1}, sf2AdmMax_{-1};
    mutable std::vector<Operon::Scalar> sf2LexEps_;
    void sf2ComputeLexEps(std::vector<PairIndividual> const& pop) const
    {
        sf2LexEps_.clear();
        std::size_t R = 0;
        for (auto const& p : pop) { if (!p.regionErrCache.empty()) { R = p.regionErrCache.size(); break; } }
        if (R == 0) { return; }
        sf2LexEps_.assign(R, Operon::Scalar{0});
        std::vector<double> col;
        col.reserve(pop.size());
        for (std::size_t k = 0; k < R; ++k) {
            col.clear();
            for (auto const& p : pop) {
                if (k < p.regionErrCache.size() && std::isfinite(p.regionErrCache[k])) { col.push_back(p.regionErrCache[k]); }
            }
            if (col.empty()) { continue; }
            auto m1 = col.begin() + static_cast<std::ptrdiff_t>(col.size() / 2);
            std::nth_element(col.begin(), m1, col.end());
            double const med = *m1;
            for (auto& v : col) { v = std::abs(v - med); }
            auto m2 = col.begin() + static_cast<std::ptrdiff_t>(col.size() / 2);
            std::nth_element(col.begin(), m2, col.end());
            sf2LexEps_[k] = static_cast<Operon::Scalar>(*m2);
        }
    }
    MultitreeConfig config_;
    NoiseModel noise_{};
    Operon::Scalar noiseMeanS2_{1};
    std::vector<Operon::Scalar> externalPointSigma2_{};
    std::vector<PairIndividual> eliteSeeds_{};
    mutable std::atomic<std::size_t> nConstDescents_{0};
    mutable std::atomic<std::size_t> nConstLossCalls_{0};
    mutable std::atomic<std::size_t> nConstIters_{0};
    mutable std::atomic<std::size_t> nLmStopUnresolvableDamped_{0};
    mutable std::atomic<std::size_t> nLmStopRealizedSmall_{0};
    mutable std::atomic<std::size_t> nLmStopFixedPoint_{0};
    mutable std::array<std::atomic<std::size_t>, 64> lmIterHist_{};
    std::shared_ptr<CclmGpuRuntime> cclmGpu_{};
    bool cclmGpuTried_{false};
    double cclmRateCpu_{0.0};
    double cclmRateGpu_{0.0};
    std::vector<ConcurrentLmGenStat> cclmStats_{};
    std::vector<ConcurrentLmParityGen> cclmParity_{};
    CclmFitData cclmFit_{};
    std::vector<CclmJob> cclmJobs_{};
    std::vector<std::size_t> cclmGpuPos_{}, cclmCpuPos_{};
    std::vector<CclmOut> cclmOuts_{};
    bool cclmPrevUncollected_{false};
    mutable std::atomic<std::size_t> nEvalErrors_{0};
    mutable std::atomic<std::size_t> nNonFiniteFitness_{0};
    mutable std::atomic<std::size_t> nInvalidChildren_{0};
    mutable std::atomic<std::size_t> nNonFiniteMutatedConst_{0};
    mutable std::atomic<std::size_t> nNonFiniteLmConst_{0};
    mutable std::atomic<std::size_t> nMutationOversizeRetries_{0};
    mutable std::atomic<std::size_t> nMutationFallbackToParent_{0};
    mutable std::atomic<std::size_t> nRejectedNonFiniteBeforeConstOpt_{0};
    mutable std::atomic<std::size_t> nLmRecoveredFromNonFiniteStart_{0};
    mutable std::atomic<std::size_t> nNonFiniteWithLog_{0};
    mutable std::atomic<std::size_t> nNonFiniteWithSqrt_{0};
    mutable std::atomic<std::size_t> nNonFiniteWithDiv_{0};
    mutable std::atomic<std::size_t> nNonFiniteWithExp_{0};
    mutable std::atomic<std::size_t> nNonFiniteWithNoneOfThese_{0};
    mutable std::atomic<std::size_t> nCascadeCullGens_{0};
    mutable std::atomic<std::size_t> nCascadeCullSeen_{0};
    mutable std::atomic<std::size_t> nCascadeCullSurvived_{0};
    mutable std::atomic<std::size_t> nLiveCoarseSieveSeen_{0};
    mutable std::atomic<std::size_t> nLiveCoarseSieveCulled_{0};
    mutable std::atomic<std::size_t> nLiveCoarseSieveFinite_{0};
    mutable std::atomic<double> liveCoarseSieveBoundMin_{std::numeric_limits<double>::infinity()};
    mutable std::atomic<double> liveCoarseSieveBoundMax_{-std::numeric_limits<double>::infinity()};
    mutable std::atomic<double> liveCoarseSieveBoundSum_{0.0};
#ifdef SKETCH_DUP_STATS
    mutable DupStats dupStats_;
#endif
    mutable std::atomic<std::size_t> nSemCross_{0};
    mutable std::atomic<std::size_t> nSemFallback_{0};

    mutable std::atomic<int64_t> profEvalNs_{0};
    mutable std::atomic<int64_t> profLmNs_{0};
    mutable std::atomic<int64_t> profScoreNs_{0};
    mutable std::atomic<int64_t> profBreedNs_{0};
    mutable std::atomic<int64_t> profSelectNs_{0};
    mutable std::atomic<int64_t> profGenNs_{0};
    mutable std::atomic<int64_t> profSelectCombineNs_{0};
    mutable std::atomic<int64_t> profSelectNdsNs_{0};
    mutable std::atomic<int64_t> profSelectCrowdNs_{0};

    mutable std::atomic<int64_t> profScoreDistNs_{0};
    mutable std::atomic<int64_t> profScoreCisNs_{0};
    mutable std::atomic<int64_t> profScoreOrderNs_{0};
    mutable std::atomic<int64_t> profScoreCurvNs_{0};
    mutable GradMutStats gmStats_{};
    mutable MutOpDiagCounters mutDiag_{};
    mutable std::vector<GradMutRecord> gmLog_;
    mutable std::vector<GradMutEff> gmEff_;
    mutable std::mutex gmLogMu_;
    mutable std::atomic<int> gmEventId_{0};
    mutable bool gmTrajStarted_{false};
    mutable std::chrono::steady_clock::time_point gmTrajT0_{};
    mutable double gmTrajHookMs_{0.0};
    mutable double gmArmN_[7]{}, gmArmSum_[7]{}, gmArmSumSq_[7]{};
    mutable std::array<std::vector<double>, 8> gmIrP_{}, gmIrY_{}, gmIrD_{};
    mutable std::array<std::size_t, 8> gmIrN_{}, gmIrN0_{};
    mutable std::array<char, 8> gmIrFit_{};
    mutable std::unordered_set<std::uint64_t> gmIrSeen_;
    mutable std::array<GradMutErrModel, 8> gmIrModel_ = [] {
        std::array<GradMutErrModel, 8> a{};
        for (auto& e : a) { e.pi0 = 0.5; }
        return a;
    }();
    mutable std::atomic<int64_t> profLmBuildNs_{0};
    mutable std::atomic<int64_t> profLmLossNs_{0};

    mutable std::atomic<int64_t> profCisPointSegNs_{0};
    mutable std::atomic<int64_t> profCisNllNs_{0};
    mutable std::atomic<int64_t> profCisReverseNs_{0};
    mutable std::atomic<int64_t> profCisGrid48Ns_{0};

    mutable std::atomic<int64_t> profGpuPrepassNs_{0};
    mutable std::atomic<int64_t> profParallelWallNs_{0};

    template <typename F>
    void parallelApply(std::size_t n, std::size_t workers, F&& f) const
    {
        if (workers <= 1 || n <= 1) {
            for (std::size_t i = 0; i < n; ++i) { f(i); }
            return;
        }
        workers = std::min(workers, n);
        bool const credit = perTreeCreditMode();
        Operon::Scalar const clamp = errorClampValue();
        bool const floorMode = floorOnChamferMode();
        bool const ptClamp = pointClampMode();
        Operon::Scalar const arcGuard = arcLengthGuardRatio();
        bool const constMag = constMagnitudeComplexityMode();
        bool const cisFit = cisFitnessMode();
        constexpr std::size_t kChunk = 4;
        std::atomic<std::size_t> next{0};
        std::function<void()> body = [&f, &next, credit, clamp, floorMode, ptClamp, arcGuard,
                               constMag, cisFit, n]() {
#ifdef __APPLE__
                pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
                CreditModeGuard cg(credit);
                FloorClampGuard fg(clamp);
                FloorModeGuard fm(floorMode);
                PointClampGuard pg(ptClamp);
                ArcLengthGuardRatioGuard ag(arcGuard);
                ConstMagnitudeComplexityGuard cmg(constMag);
                CisFitnessModeGuard cfg(cisFit);
                for (;;) {
                    std::size_t const lo = next.fetch_add(kChunk, std::memory_order_relaxed);
                    if (lo >= n) { break; }
                    std::size_t const hi = std::min(lo + kChunk, n);
                    for (std::size_t i = lo; i < hi; ++i) { f(i); }
                }
            };

#ifndef _WIN32
        if (persistentPoolEnabledFromEnv()) {
            auto& pool = PersistentWorkerPool::instance();
            auto dispatcherLock = pool.tryAcquireDispatcher();
            if (dispatcherLock.owns_lock()) {
                pool.dispatch(workers, body);
                return;
            }
        }
        constexpr std::size_t kWorkerStackBytes = std::size_t{16} << 20;
        std::vector<std::unique_ptr<std::function<void()>>> ctxs;
        std::vector<pthread_t> tids;
        ctxs.reserve(workers); tids.reserve(workers);
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, kWorkerStackBytes);
        for (std::size_t w = 0; w < workers; ++w) {
            ctxs.push_back(std::make_unique<std::function<void()>>(body));
            pthread_t tid{};
            if (pthread_create(&tid, &attr, &MultitreeGP::workerTrampoline, ctxs.back().get()) == 0) {
                tids.push_back(tid);
            } else {
                (*ctxs.back())();
            }
        }
        pthread_attr_destroy(&attr);
        for (pthread_t t : tids) { pthread_join(t, nullptr); }
#else
        std::vector<std::thread> pool;
        pool.reserve(workers);
        for (std::size_t w = 0; w < workers; ++w) { pool.emplace_back(body); }
        for (auto& th : pool) { th.join(); }
#endif
    }

#ifndef _WIN32
    static void* workerTrampoline(void* p) {
        (*static_cast<std::function<void()>*>(p))();
        return nullptr;
    }
#endif

    struct EvaluatorArena {
        std::vector<Operon::Scalar> gx;
        std::vector<Operon::Scalar> gy;
        std::vector<Operon::Scalar> scratch_x;
        std::vector<Operon::Scalar> scratch_y;
    };

    static auto getEvaluatorArena() -> EvaluatorArena&
    {
        thread_local EvaluatorArena arena;
        return arena;
    }

    public:
    void interpInto(Operon::Tree const& tree, std::vector<Operon::Scalar> const& coeffs,
                    Operon::Dataset const& ds, Operon::Range range,
                    std::vector<Operon::Scalar>& out) const
    {
        out.resize(range.Size());
        if (config_.useRegisterInterp) {
            auto bc = BytecodeProgram::Compile(tree, {coeffs.data(), coeffs.size()});
            if (bc.valid) {
                auto rp = Sketch::RegisterProgram::FromBytecode(bc);
                if (rp.valid) {
                    auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
                    Sketch::RegisterInterpreter::EvaluateSingle(rp, tValues, {out.data(), out.size()});
                    return;
                }
            }
        }
        if (config_.useBytecodeEvaluator) {
            auto prog = BytecodeProgram::Compile(tree, {coeffs.data(), coeffs.size()});
            if (prog.valid) {
                auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
                prog.Evaluate(tValues, out);
                return;
            }
        }
        Operon::Interpreter<Operon::Scalar, DTable> in(dtable_, &ds, &tree);
        in.Evaluate({coeffs.data(), coeffs.size()}, range, {out.data(), out.size()});
    }

    auto interp(Operon::Tree const& tree, std::vector<Operon::Scalar> const& coeffs,
                Operon::Dataset const& ds, Operon::Range range) const
        -> std::vector<Operon::Scalar>
    {
        std::vector<Operon::Scalar> out;
        interpInto(tree, coeffs, ds, range, out);
        return out;
    }
    private:

    static inline thread_local bool mgTauCacheValid_ = false;
    static inline thread_local std::vector<Operon::Scalar> mgTauCache_;
    static inline thread_local std::size_t mgTauCacheN_ = 0, mgTauCacheM_ = 0;
    static inline thread_local Operon::Scalar mgTauCacheGx0_ = 0, mgTauCacheGy0_ = 0,
                                                mgTauCacheGxL_ = 0, mgTauCacheGyL_ = 0;

    auto monotonicityOrderTerm(std::vector<Operon::Scalar> const& gx,
                                std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto const M = gx.size();
        auto const N = dataX_.size();
        if (M < 2 || N < 2) { return BIG; }

        if (oneResOn()) {
            thread_local std::vector<Operon::Scalar> oRes, oTau;
            if (!oneResidualOf(gx, gy, true, oRes, &oTau) || oTau.size() != N) { return BIG; }
            if (sf2(2)) { return sf2OrderValue(oTau); }
            Operon::Scalar const epsJ = monotonicityEpsJitter(oTau);
            return monotonicityPenalty(oTau, epsJ);
        }

        if (mgTauCacheValid_ && mgTauCacheN_ == N && mgTauCacheM_ == M
            && mgTauCacheGx0_ == gx[0] && mgTauCacheGy0_ == gy[0]
            && mgTauCacheGxL_ == gx[M - 1] && mgTauCacheGyL_ == gy[M - 1]) {
            mgTauCacheValid_ = false;
            if (sf2(2)) { return sf2OrderValue(mgTauCache_); }
            Operon::Scalar const epsJitter = monotonicityEpsJitter(mgTauCache_);
            return monotonicityPenalty(mgTauCache_, epsJitter);
        }

        thread_local std::vector<Operon::Scalar> mgCx, mgCy;
        mgCx.resize(M);
        mgCy.resize(M);
        Operon::Scalar gx0 = gx[0], gy0 = gy[0];
        if (!std::isfinite(gx0) || !std::isfinite(gy0)) { return BIG; }
        Operon::Scalar xlo = gx0, xhi = gx0, ylo = gy0, yhi = gy0;
        mgCx[0] = gx0; mgCy[0] = gy0;
        for (std::size_t j = 1; j < M; ++j) {
            Operon::Scalar const vx = gx[j], vy = gy[j];
            if (!std::isfinite(vx) || !std::isfinite(vy)) { return BIG; }
            mgCx[j] = vx; mgCy[j] = vy;
            if (vx < xlo) { xlo = vx; } else if (vx > xhi) { xhi = vx; }
            if (vy < ylo) { ylo = vy; } else if (vy > yhi) { yhi = vy; }
        }
        Operon::Scalar const cxc = (xhi + xlo) / 2, cyc = (yhi + ylo) / 2;
        Operon::Scalar const sp = std::max({xhi - xlo, yhi - ylo, Operon::Scalar{1e-9}});
        for (std::size_t j = 0; j < M; ++j) {
            mgCx[j] = (mgCx[j] - cxc) / sp;
            mgCy[j] = (mgCy[j] - cyc) / sp;
        }

        thread_local std::vector<Vec2> mgStrokePts, mgCandPts;
        mgStrokePts.resize(N);
        mgCandPts.resize(M);
        for (std::size_t i = 0; i < N; ++i) { mgStrokePts[i] = Vec2{dataX_[i], dataY_[i]}; }
        for (std::size_t j = 0; j < M; ++j) { mgCandPts[j] = Vec2{mgCx[j], mgCy[j]}; }

        PointSegmentConfig psCfg;
        psCfg.computeTau = true;
        auto const psRes = evaluatePointSegment(mgStrokePts, mgCandPts, psCfg);
        if (!psRes.ok || psRes.tau.size() != N) { return BIG; }

        if (sf2(2)) { return sf2OrderValue(psRes.tau); }
        Operon::Scalar const epsJitter = monotonicityEpsJitter(psRes.tau);
        return monotonicityPenalty(psRes.tau, epsJitter);
    }

    auto orderTerm(std::vector<Operon::Scalar> const& gx,
                   std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        auto const rawOrder = config_.useMonotonicityOrderGuard
            ? monotonicityOrderTerm(gx, gy)
            : (config_.useBandedOrderDistance
                ? orderDistanceBanded(gx, gy, dataX_, dataY_, dataClosed_, config_.orderBandWidth, 64, &dataArcX_, &dataArcY_)
                : (config_.useArcOrderDistance
                    ? orderDistanceArc(gx, gy, dataX_, dataY_, dataClosed_, 64, &dataArcX_, &dataArcY_)
                    : orderDistance(gx, gy, dataX_, dataY_, dataClosed_, 64, &dataArcX_, &dataArcY_)));
        if (sf2(2) && config_.useMonotonicityOrderGuard) { return rawOrder; }
        auto od = config_.orderWeight * rawOrder;
        if (config_.whitenOrderTerm && !noise_.empty()) { od /= noiseMeanS2_; }
        return od;
    }

    auto curvatureTerm(std::vector<Operon::Scalar> const& gx,
                       std::vector<Operon::Scalar> const& gy) const -> Operon::Scalar
    {
        double const tc = turningFourier(gx, gy, dataClosed_, config_.curvatureHarmonics);
        double const excess = std::abs(tc - dataTurningF_) - curvatureFloor_;
        if (!(excess > 0.0)) { return Operon::Scalar{0}; }
        auto term = config_.curvatureWeight * static_cast<Operon::Scalar>(excess * excess);
        if (config_.curvatureMultiplicative) { return term; }
        if (config_.whitenOrderTerm && !noise_.empty()) { term /= noiseMeanS2_; }
        return term;
    }

    auto applyCurvature(Operon::Scalar sd, Operon::Scalar term) const -> Operon::Scalar
    {
        if (!std::isfinite(sd)) { return sd; }
        return config_.curvatureMultiplicative ? sd * (Operon::Scalar{1} + term) : sd + term;
    }

    auto shapeOnly(PairIndividual const& ind, Operon::Dataset const& ds,
                   Operon::Range range) const -> Operon::Scalar
    {
        return shapeOnlyParts(ind, ds, range, nullptr);
    }
