
    static auto dataSeedTurnCounts(double T, std::size_t maxCount) -> std::vector<long>
    {
        std::vector<long> out;
        if (maxCount == 0 || !(T >= 1.0)) { return out; }
        long const klo = std::max(1L, static_cast<long>(std::floor(0.5 * T)));
        long const khi = std::max(klo, static_cast<long>(std::llround(T)));
        long kc = static_cast<long>(std::llround(0.5 * (static_cast<double>(klo) + static_cast<double>(khi))));
        kc = std::min(std::max(kc, klo), khi);
        out.push_back(kc);
        for (long d = 1; (kc - d) >= klo || (kc + d) <= khi; ++d) {
            if (out.size() >= maxCount) { break; }
            if (kc + d <= khi) { out.push_back(kc + d); }
            if (out.size() >= maxCount) { break; }
            if (kc - d >= klo) { out.push_back(kc - d); }
        }
        if (out.size() > maxCount) { out.resize(maxCount); }
        return out;
    }

    static auto freqRestartValue(Operon::Scalar v, std::size_t j) -> Operon::Scalar
    {
        constexpr double TWO_PI = 6.283185307179586;
        if (j == 1) { return -v; }
        auto const s = (v < Operon::Scalar{0}) ? -1.0 : 1.0;
        auto const k0 = static_cast<long long>(std::llround(std::abs(static_cast<double>(v)) / TWO_PI));
        auto const m = static_cast<long long>((j - 2) / 2 + 1);
        auto k = ((j - 2) % 2 == 0) ? k0 + m : k0 - m;
        if (k < 1) { k = k0 + m + 1; }
        return static_cast<Operon::Scalar>(s * static_cast<double>(k) * TWO_PI);
    }

    auto treeJacobian(Operon::Tree const& tree, std::vector<Operon::Scalar> const& c,
                      Operon::Dataset const& ds, Operon::Range range, std::size_t M,
                      std::vector<double>& out) const -> bool
    {
        auto const nc = c.size();
        out.assign(M * nc, 0.0);
        if (nc == 0) { return true; }
        Operon::Interpreter<Operon::Scalar, DTable> in(dtable_, &ds, &tree);
        auto J = in.JacFwd({c.data(), c.size()}, range);
        if (static_cast<std::size_t>(J.rows()) != M || static_cast<std::size_t>(J.cols()) != nc) { return false; }
        static thread_local std::vector<char> bad_buf;
        bad_buf.assign(nc, 0);
        auto& bad = bad_buf;
        if (J.allFinite()) {
            for (std::size_t k = 0; k < nc; ++k) {
                for (std::size_t j = 0; j < M; ++j) {
                    out[j * nc + k] = static_cast<double>(J(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)));
                }
            }
            return true;
        }
        bool any = false;
        for (std::size_t k = 0; k < nc; ++k) {
            for (std::size_t j = 0; j < M; ++j) {
                auto const g = static_cast<double>(J(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)));
                if (std::isfinite(g)) { out[j * nc + k] = g; }
                else { out[j * nc + k] = std::numeric_limits<double>::quiet_NaN(); bad[k] = 1; any = true; }
            }
        }
        if (!any) { return true; }
        static thread_local std::vector<Operon::Scalar> cp_buf, cm_buf, fp_buf, fm_buf;
        for (std::size_t k = 0; k < nc; ++k) {
            if (bad[k] == 0) { continue; }
            double const h = 1e-4 * std::max(std::abs(static_cast<double>(c[k])), 1.0);
            cp_buf = c;
            cm_buf = c;
            cp_buf[k] = static_cast<Operon::Scalar>(static_cast<double>(c[k]) + h);
            cm_buf[k] = static_cast<Operon::Scalar>(static_cast<double>(c[k]) - h);
            interpInto(tree, cp_buf, ds, range, fp_buf);
            interpInto(tree, cm_buf, ds, range, fm_buf);
            if (fp_buf.size() != M || fm_buf.size() != M) { return false; }
            for (std::size_t j = 0; j < M; ++j) {
                if (std::isfinite(out[j * nc + k])) { continue; }
                double const g = (static_cast<double>(fp_buf[j]) - static_cast<double>(fm_buf[j])) / (2 * h);
                out[j * nc + k] = std::isfinite(g) ? g : 0.0;
            }
        }
        return true;
    }

    struct OneResLmCore {
        PointSegOneResult<double> ps;
        std::vector<std::size_t> dataRegion;
        std::vector<double> e2, zx, zy;
        std::vector<std::size_t> curveRegion;
        std::vector<int> curveCorr;
    };
    auto oneResLmCoreEval(double const* qx, double const* qy, std::size_t M,
                          std::vector<Operon::Scalar> const& dataXr, std::vector<Operon::Scalar> const& dataYr,
                          CurveBlocks const& dataBlocksR, std::vector<std::size_t> const* origIdxOv,
                          OneResLmCore& c) const -> bool
    {
        auto const N = dataXr.size();
        auto const K = noise_.size();
        auto const NFull = dataX_.size();
        if (N < 2 || M < 2) { return false; }
        thread_local std::vector<OneResP2d> dataD, candD, drawnD;
        dataD.resize(N);
        candD.resize(M);
        for (std::size_t i = 0; i < N; ++i) { dataD[i] = OneResP2d{static_cast<double>(dataXr[i]), static_cast<double>(dataYr[i])}; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(qx[j]) || !std::isfinite(qy[j])) { return false; }
            candD[j] = OneResP2d{qx[j], qy[j]};
        }
        evaluatePointSegmentOneInto<double>(dataD.data(), N, candD.data(), M,true,true, c.ps, oneResUnsigned());
        if (!c.ps.ok) { return false; }
        c.dataRegion.resize(N);
        for (std::size_t i = 0; i < N; ++i) {
            c.dataRegion[i] = (origIdxOv != nullptr)
                ? noiseIndexAtArc(originalArcIndex((*origIdxOv)[i]), NFull, K)
                : noiseIndexAtArc(originalArcIndex(i), N, K);
        }
        c.e2.resize(M); c.zx.resize(M); c.zy.resize(M); c.curveRegion.resize(M); c.curveCorr.resize(M);
        bool const segRev = oneResBit(2) && origIdxOv == nullptr;
        if (segRev) {
            drawnD.resize(N);
            for (std::size_t i = 0; i < N; ++i) { drawnD[originalArcIndex(i)] = dataD[i]; }
        }
        int ptHint = -1;
        for (std::size_t j = 0; j < M; ++j) {
            if (segRev) {
                auto const h = pointToPolylineSegOne(qx[j], qy[j], drawnD.data(), N);
                c.e2[j] = h.d2; c.zx[j] = h.zx; c.zy[j] = h.zy;
                c.curveRegion[j] = noiseIndexAtArc(h.nearestVertex, N, K);
                c.curveCorr[j] = -1;
            } else {
                auto const a = nearestPointD(qx[j], qy[j], dataXr, dataYr, dataBlocksR, &ptHint);
                std::size_t const kb = a.i;
                double const ddx = qx[j] - static_cast<double>(dataXr[kb]);
                double const ddy = qy[j] - static_cast<double>(dataYr[kb]);
                c.e2[j] = ddx * ddx + ddy * ddy;
                c.zx[j] = static_cast<double>(dataXr[kb]); c.zy[j] = static_cast<double>(dataYr[kb]);
                c.curveRegion[j] = c.dataRegion[kb];
                c.curveCorr[j] = static_cast<int>(kb);
            }
        }
        return true;
    }

    auto oneResLmWeight(double u, std::size_t k) const -> double
    {
        if (noise_.size() == 0) { return 1.0; }
        double const s2 = static_cast<double>(noise_.s2[k]);
        double const nu = static_cast<double>(noise_.nu[k]);
        if (!(s2 > 0.0)) { return 1.0; }
        if (nu <= 0.0) { return 1.0 / s2; }
        return 0.5 * (nu + 2.0) / (nu * s2 + u);
    }

    auto oneResLmLoss(std::vector<Operon::Scalar> const& gx, std::vector<Operon::Scalar> const& gy,
                      std::vector<Operon::Scalar> const& dataXr, std::vector<Operon::Scalar> const& dataYr,
                      CurveBlocks const& dataBlocksR, std::vector<std::size_t> const* origIdxOv) const -> Operon::Scalar
    {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto const M = gx.size();
        auto const N = dataXr.size();
        if (M < 2 || N < 2 || gy.size() != M) { return BIG; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(gx[j]) || !std::isfinite(gy[j])) { return BIG; }
        }
        std::size_t iUmin = 0, iUmax = 0, iVmin = 0, iVmax = 0;
        for (std::size_t j = 1; j < M; ++j) {
            if (gx[j] < gx[iUmin]) { iUmin = j; }
            if (gx[j] > gx[iUmax]) { iUmax = j; }
            if (gy[j] < gy[iVmin]) { iVmin = j; }
            if (gy[j] > gy[iVmax]) { iVmax = j; }
        }
        double const Cx = 0.5 * (static_cast<double>(gx[iUmax]) + static_cast<double>(gx[iUmin]));
        double const Cy = 0.5 * (static_cast<double>(gy[iVmax]) + static_cast<double>(gy[iVmin]));
        double const spanU = static_cast<double>(gx[iUmax]) - static_cast<double>(gx[iUmin]);
        double const spanV = static_cast<double>(gy[iVmax]) - static_cast<double>(gy[iVmin]);
        double const S = std::max({spanU, spanV, 1e-9});
        if (!(S > 0.0) || !std::isfinite(S)) { return BIG; }
        double const invS = 1.0 / S;
        thread_local std::vector<double> qdx, qdy;
        qdx.resize(M); qdy.resize(M);
        for (std::size_t j = 0; j < M; ++j) {
            qdx[j] = (static_cast<double>(gx[j]) - Cx) * invS;
            qdy[j] = (static_cast<double>(gy[j]) - Cy) * invS;
        }
        thread_local OneResLmCore c;
        if (!oneResLmCoreEval(qdx.data(), qdy.data(), M, dataXr, dataYr, dataBlocksR, origIdxOv, c)) { return BIG; }
        auto const K = noise_.size();
        auto rho = [&](double u, std::size_t k) -> double {
            if (K == 0) { return u; }
            return static_cast<double>(pointLoss(static_cast<Operon::Scalar>(u), noise_.s2[k], noise_.nu[k]));
        };
        double accD = 0.0;
        for (std::size_t i = 0; i < N; ++i) { double const r = c.ps.residual[i]; accD += rho(r * r, c.dataRegion[i]); }
        double accC = 0.0;
        for (std::size_t j = 0; j < M; ++j) { accC += rho(c.e2[j], c.curveRegion[j]); }
        double const L = 0.5 * (accD / static_cast<double>(N) + accC / static_cast<double>(M));
        if (!std::isfinite(L)) { return BIG; }
        return static_cast<Operon::Scalar>(L);
    }

    auto lmBuild(PairIndividual const& ind, std::vector<Operon::Scalar> const& v, std::size_t nx,
                 Operon::Dataset const& ds, Operon::Range range,
                 std::vector<double>& resid, std::vector<double>& jac,
                 std::vector<int>* corr = nullptr,
                 BytecodeProgram* progXCache = nullptr, BytecodeProgram* progYCache = nullptr,
                 std::vector<Operon::Scalar> const* dataXOv = nullptr,
                 std::vector<Operon::Scalar> const* dataYOv = nullptr,
                 CurveBlocks const* dataBlocksOv = nullptr,
                 std::vector<std::size_t> const* origIdxOv = nullptr) const -> bool
    {
        auto const& dataXr = (dataXOv != nullptr) ? *dataXOv : dataX_;
        auto const& dataYr = (dataYOv != nullptr) ? *dataYOv : dataY_;
        auto const& dataBlocksR = (dataBlocksOv != nullptr) ? *dataBlocksOv : dataBlocks_;
        auto const N = dataXr.size();
        auto const n = v.size();
        auto const ny = n - nx;
        static thread_local std::vector<Operon::Scalar> vx_buf, vy_buf, U_buf, V_buf, qx_buf, qy_buf;
        static thread_local std::vector<double> Ju_buf, Jv_buf, dCx_buf, dCy_buf, dS_buf, dqx_buf, dqy_buf;

        vx_buf.assign(v.begin(), v.begin() + static_cast<long>(nx));
        vy_buf.assign(v.begin() + static_cast<long>(nx), v.end());
        bool const useCacheX = progXCache != nullptr && progXCache->valid;
        bool const useCacheY = progYCache != nullptr && progYCache->valid;
        if (useCacheX) {
            U_buf.resize(range.Size());
            progXCache->UpdateCoefficients({vx_buf.data(), vx_buf.size()});
            auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
            progXCache->Evaluate(tValues, U_buf);
        } else {
            interpInto(ind.x, vx_buf, ds, range, U_buf);
        }
        if (useCacheY) {
            V_buf.resize(range.Size());
            progYCache->UpdateCoefficients({vy_buf.data(), vy_buf.size()});
            auto tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());
            progYCache->Evaluate(tValues, V_buf);
        } else {
            interpInto(ind.y, vy_buf, ds, range, V_buf);
        }
        auto const& U = U_buf;
        auto const& V = V_buf;
        auto const M = U.size();
        if (M < 2 || N < 2 || V.size() != M) { return false; }
        for (std::size_t j = 0; j < M; ++j) {
            if (!std::isfinite(U[j]) || !std::isfinite(V[j])) { return false; }
        }

        if (!treeJacobian(ind.x, vx_buf, ds, range, M, Ju_buf)) { return false; }
        if (!treeJacobian(ind.y, vy_buf, ds, range, M, Jv_buf)) { return false; }
        auto const& Ju = Ju_buf;
        auto const& Jv = Jv_buf;
        auto ju = [&](std::size_t j, std::size_t k) -> double { return Ju[j * nx + k]; };
        auto jv = [&](std::size_t j, std::size_t k) -> double { return Jv[j * ny + k]; };

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
        if (!(S > 0.0) || !std::isfinite(S)) { return false; }

        dCx_buf.assign(n, 0.0);
        dCy_buf.assign(n, 0.0);
        dS_buf.assign(n, 0.0);
        auto& dCx = dCx_buf;
        auto& dCy = dCy_buf;
        auto& dS = dS_buf;
        for (std::size_t k = 0; k < nx; ++k) { dCx[k] = 0.5 * (ju(iUmax, k) + ju(iUmin, k)); }
        for (std::size_t k = 0; k < ny; ++k) { dCy[nx + k] = 0.5 * (jv(iVmax, k) + jv(iVmin, k)); }
        if (spanU >= spanV && spanU > 1e-9) {
            for (std::size_t k = 0; k < nx; ++k) { dS[k] = ju(iUmax, k) - ju(iUmin, k); }
        } else if (spanV > 1e-9) {
            for (std::size_t k = 0; k < ny; ++k) { dS[nx + k] = jv(iVmax, k) - jv(iVmin, k); }
        }

        qx_buf.resize(M);
        qy_buf.resize(M);
        dqx_buf.resize(M * n);
        dqy_buf.resize(M * n);
        auto& qx = qx_buf;
        auto& qy = qy_buf;
        auto& dqx = dqx_buf;
        auto& dqy = dqy_buf;
        double const invS = 1.0 / S;
        bool const t11Lm = oneResBit(4);
        thread_local std::vector<double> qdx_buf, qdy_buf;
        if (t11Lm) { qdx_buf.resize(M); qdy_buf.resize(M); }
        for (std::size_t j = 0; j < M; ++j) {
            double const ax = (static_cast<double>(U[j]) - Cx) * invS;
            double const ay = (static_cast<double>(V[j]) - Cy) * invS;
            qx[j] = static_cast<Operon::Scalar>(ax);
            qy[j] = static_cast<Operon::Scalar>(ay);
            if (t11Lm) { qdx_buf[j] = ax; qdy_buf[j] = ay; }
            double* rx = dqx.data() + j * n;
            double* ry = dqy.data() + j * n;
            for (std::size_t k = 0; k < n; ++k) {
                double const gu = (k < nx) ? ju(j, k) : 0.0;
                double const gv = (k >= nx) ? jv(j, k - nx) : 0.0;
                rx[k] = (gu - dCx[k]) * invS - ax * dS[k] * invS;
                ry[k] = (gv - dCy[k]) * invS - ay * dS[k] * invS;
            }
        }

        auto const K = noise_.size();
        auto const NFull = dataX_.size();
        auto regionOf = [&](std::size_t i) -> std::size_t {
            if (origIdxOv != nullptr) { return noiseIndexAtArc(originalArcIndex((*origIdxOv)[i]), NFull, K); }
            return noiseIndexAtArc(originalArcIndex(i), N, K);
        };
        auto weightOf = [&](double u, std::size_t k) -> double {
            if (K == 0) { return 1.0; }
            double const s2 = static_cast<double>(noise_.s2[k]);
            double const nu = static_cast<double>(noise_.nu[k]);
            if (!(s2 > 0.0)) { return 1.0; }
            if (nu <= 0.0) { return 1.0 / s2; }
            return 0.5 * (nu + 2.0) / (nu * s2 + u);
        };

        auto const m = N + M;
        resid.assign(m, 0.0);
        jac.assign(m * n, 0.0);
        if (corr != nullptr) {
            corr->assign(m + 4, -1);
            (*corr)[m + 0] = static_cast<int>(iUmin); (*corr)[m + 1] = static_cast<int>(iUmax);
            (*corr)[m + 2] = static_cast<int>(iVmin); (*corr)[m + 3] = static_cast<int>(iVmax);
        }
        constexpr double EPS = 1e-12;

        if (t11Lm) {
            thread_local OneResLmCore c;
            if (!oneResLmCoreEval(qdx_buf.data(), qdy_buf.data(), M, dataXr, dataYr, dataBlocksR, origIdxOv, c)) { return false; }
            for (std::size_t i = 0; i < N; ++i) {
                double const r = c.ps.residual[i];
                double const w = weightOf(r * r, c.dataRegion[i]);
                double const scale = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(N)));
                resid[i] = scale * r;
                double tau = c.ps.tau[i];
                tau = std::clamp(tau, 0.0, static_cast<double>(M - 1));
                auto s0 = static_cast<std::size_t>(std::floor(tau));
                if (s0 > M - 2) { s0 = M - 2; }
                double const f = std::clamp(tau - static_cast<double>(s0), 0.0, 1.0);
                if (corr != nullptr) { (*corr)[i] = static_cast<int>(s0); }
                double const nhx = c.ps.normalX[i], nhy = c.ps.normalY[i];
                double* row = jac.data() + i * n;
                double const* r0x = dqx.data() + s0 * n;       double const* r0y = dqy.data() + s0 * n;
                double const* r1x = dqx.data() + (s0 + 1) * n; double const* r1y = dqy.data() + (s0 + 1) * n;
                for (std::size_t k = 0; k < n; ++k) {
                    double const g0 = nhx * r0x[k] + nhy * r0y[k];
                    double const g1 = nhx * r1x[k] + nhy * r1y[k];
                    row[k] = -scale * ((1.0 - f) * g0 + f * g1);
                }
            }
            for (std::size_t j = 0; j < M; ++j) {
                double const e2 = c.e2[j];
                double const w = weightOf(e2, c.curveRegion[j]);
                double const scale = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(M)));
                double const e = std::sqrt(e2);
                resid[N + j] = scale * e;
                if (corr != nullptr) { (*corr)[N + j] = c.curveCorr[j]; }
                if (e <= EPS) { continue; }
                double const mhx = (qdx_buf[j] - c.zx[j]) / e;
                double const mhy = (qdy_buf[j] - c.zy[j]) / e;
                double* row = jac.data() + (N + j) * n;
                double const* rx = dqx.data() + j * n;
                double const* ry = dqy.data() + j * n;
                for (std::size_t k = 0; k < n; ++k) { row[k] = scale * (mhx * rx[k] + mhy * ry[k]); }
            }
            for (auto const& z : resid) { if (!std::isfinite(z)) { return false; } }
            for (auto const& z : jac) { if (!std::isfinite(z)) { return false; } }
            return true;
        }

        bool const useBlk = (config_.lmBlockNn != 0);
        thread_local CurveBlocks qBlk;
        if (useBlk) { buildBlocks(qx, qy,true, qBlk); }

        int segHint = -1;
        for (std::size_t i = 0; i < N; ++i) {
            double const px = static_cast<double>(dataXr[i]), py = static_cast<double>(dataYr[i]);
            std::size_t jb = 0; double tau = 0.0, best = std::numeric_limits<double>::max();
            if (useBlk) {
                auto const a = nearestSegmentD(px, py, qx, qy, qBlk, &segHint);
                best = a.d2; jb = a.j; tau = a.tau;
            } else {
            for (std::size_t j = 0; j + 1 < M; ++j) {
                double const ax = static_cast<double>(qx[j]), ay = static_cast<double>(qy[j]);
                double const bx = static_cast<double>(qx[j + 1]), by = static_cast<double>(qy[j + 1]);
                double const ex = bx - ax, ey = by - ay;
                double const len2 = ex * ex + ey * ey;
                double t = 0.0;
                if (len2 > 1e-12) { t = std::clamp(((px - ax) * ex + (py - ay) * ey) / len2, 0.0, 1.0); }
                double const zx = ax + t * ex, zy = ay + t * ey;
                double const d2 = (px - zx) * (px - zx) + (py - zy) * (py - zy);
                if (d2 < best) { best = d2; jb = j; tau = t; }
            }
            }
            double const w = weightOf(best, regionOf(i));
            double const scale = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(N)));
            double const d = std::sqrt(best);
            resid[i] = scale * d;
            if (corr != nullptr) { (*corr)[i] = static_cast<int>(jb); }
            if (d <= EPS) { continue; }
            double const ax = static_cast<double>(qx[jb]), ay = static_cast<double>(qy[jb]);
            double const bx = static_cast<double>(qx[jb + 1]), by = static_cast<double>(qy[jb + 1]);
            double const zx = ax + tau * (bx - ax), zy = ay + tau * (by - ay);
            double const nhx = (px - zx) / d, nhy = (py - zy) / d;
            double* row = jac.data() + i * n;
            double const* r0x = dqx.data() + jb * n;       double const* r0y = dqy.data() + jb * n;
            double const* r1x = dqx.data() + (jb + 1) * n; double const* r1y = dqy.data() + (jb + 1) * n;
            for (std::size_t k = 0; k < n; ++k) {
                double const g0 = nhx * r0x[k] + nhy * r0y[k];
                double const g1 = nhx * r1x[k] + nhy * r1y[k];
                row[k] = -scale * ((1.0 - tau) * g0 + tau * g1);
            }
        }

        int ptHint = -1;
        for (std::size_t j = 0; j < M; ++j) {
            double const cxj = static_cast<double>(qx[j]), cyj = static_cast<double>(qy[j]);
            std::size_t kb = 0; double best = std::numeric_limits<double>::max();
            if (useBlk) {
                auto const a = nearestPointD(cxj, cyj, dataXr, dataYr, dataBlocksR, &ptHint);
                best = a.d2; kb = a.i;
            } else {
            for (std::size_t i = 0; i < N; ++i) {
                double const ddx = cxj - static_cast<double>(dataXr[i]);
                double const ddy = cyj - static_cast<double>(dataYr[i]);
                double const d2 = ddx * ddx + ddy * ddy;
                if (d2 < best) { best = d2; kb = i; }
            }
            }
            double const w = weightOf(best, regionOf(kb));
            double const scale = std::sqrt(std::max(w, 0.0) / (2.0 * static_cast<double>(M)));
            double const e = std::sqrt(best);
            resid[N + j] = scale * e;
            if (corr != nullptr) { (*corr)[N + j] = static_cast<int>(kb); }
            if (e <= EPS) { continue; }
            double const mhx = (cxj - static_cast<double>(dataXr[kb])) / e;
            double const mhy = (cyj - static_cast<double>(dataYr[kb])) / e;
            double* row = jac.data() + (N + j) * n;
            double const* rx = dqx.data() + j * n;
            double const* ry = dqy.data() + j * n;
            for (std::size_t k = 0; k < n; ++k) { row[k] = scale * (mhx * rx[k] + mhy * ry[k]); }
        }
        for (auto const& z : resid) { if (!std::isfinite(z)) { return false; } }
        for (auto const& z : jac) { if (!std::isfinite(z)) { return false; } }
        return true;
    }

    void optConstsGn(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range,
                     bool fullObjective) const
    {
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        std::size_t const nx = cx.size();
        std::vector<Operon::Scalar> c0 = cx;
        c0.insert(c0.end(), cy.begin(), cy.end());
        if (c0.empty()) { return; }

        BytecodeProgram progX, progY;
        bool const useBc = config_.useBytecodeEvaluator;
        if (useBc) {
            progX = BytecodeProgram::Compile(ind.x, {cx.data(), cx.size()});
            progY = BytecodeProgram::Compile(ind.y, {cy.data(), cy.size()});
        }
        auto const tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());

        bool const useDs = config_.lmDownsampleK > 0 && !dataXDs_.empty() && dataXDs_.size() < dataX_.size();

        auto lossOn = [&](std::vector<Operon::Scalar> const& v, bool ds_) -> Operon::Scalar {
            thread_local std::vector<Operon::Scalar> gx, gy;
            if (useBc && progX.valid && progY.valid) {
                gx.resize(range.Size());
                gy.resize(range.Size());
                progX.UpdateCoefficients({v.data(), nx});
                progY.UpdateCoefficients({v.data() + nx, v.size() - nx});
                progX.Evaluate(tValues, gx);
                progY.Evaluate(tValues, gy);
            } else {
                std::vector<Operon::Scalar> vx(v.begin(), v.begin() + static_cast<long>(nx));
                std::vector<Operon::Scalar> vy(v.begin() + static_cast<long>(nx), v.end());
                interpInto(ind.x, vx, ds, range, gx);
                interpInto(ind.y, vy, ds, range, gy);
            }
            auto const& dxr = (ds_ && useDs) ? dataXDs_ : dataX_;
            auto const& dyr = (ds_ && useDs) ? dataYDs_ : dataY_;
            auto const* blkr = (ds_ && useDs) ? &dataBlocksDs_ : &dataBlocks_;
            auto sd = oneResBit(4)
                ? oneResLmLoss(gx, gy, dxr, dyr, *blkr, (ds_ && useDs) ? &dsOrigIdx_ : nullptr)
                : shapeDistance(gx, gy, dxr, dyr, noise_, nullptr, 4, blkr, config_.useGridSdf ? &gridSdf_ : nullptr, config_.useAmxChamfer ? &amxData_ : nullptr);
            if (fullObjective && config_.orderWeight > Operon::Scalar{0}) { sd += orderTerm(gx, gy); }
            if (fullObjective && config_.curvatureWeight > Operon::Scalar{0}) { sd = applyCurvature(sd, curvatureTerm(gx, gy)); }
            return sd;
        };
        auto loss = [&](std::vector<Operon::Scalar> const& v) -> Operon::Scalar { return lossOn(v, false); };
        bool const profT10 = config_.profileTiming;
        auto lossTimedOn = [&](std::vector<Operon::Scalar> const& v, bool ds_) -> Operon::Scalar {
            if (!profT10) { return lossOn(v, ds_); }
            auto t0 = std::chrono::steady_clock::now();
            auto sd = lossOn(v, ds_);
            auto t1 = std::chrono::steady_clock::now();
            profLmLossNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count(), std::memory_order_relaxed);
            return sd;
        };
        auto lossTimed = [&](std::vector<Operon::Scalar> const& v) -> Operon::Scalar { return lossTimedOn(v, false); };
        auto buildOn = [&](std::vector<Operon::Scalar> const& v,
                           std::vector<double>& r, std::vector<double>& J, bool ds_) -> bool {
            bool const cacheOk = config_.useLmBytecodeCache && useBc && progX.valid && progY.valid;
            BytecodeProgram* pxc = cacheOk ? &progX : nullptr;
            BytecodeProgram* pyc = cacheOk ? &progY : nullptr;
            std::vector<Operon::Scalar> const* dxOv = (ds_ && useDs) ? &dataXDs_ : nullptr;
            std::vector<Operon::Scalar> const* dyOv = (ds_ && useDs) ? &dataYDs_ : nullptr;
            CurveBlocks const* blkOv = (ds_ && useDs) ? &dataBlocksDs_ : nullptr;
            std::vector<std::size_t> const* idxOv = (ds_ && useDs) ? &dsOrigIdx_ : nullptr;
            if (!profT10) { return lmBuild(ind, v, nx, ds, range, r, J, nullptr, pxc, pyc, dxOv, dyOv, blkOv, idxOv); }
            auto t0 = std::chrono::steady_clock::now();
            auto ok = lmBuild(ind, v, nx, ds, range, r, J, nullptr, pxc, pyc, dxOv, dyOv, blkOv, idxOv);
            auto t1 = std::chrono::steady_clock::now();
            profLmBuildNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count(), std::memory_order_relaxed);
            return ok;
        };
        auto build = [&](std::vector<Operon::Scalar> const& v,
                         std::vector<double>& r, std::vector<double>& J) -> bool { return buildOn(v, r, J, false); };

        LevenbergMarquardt lm;
        if (config_.floorClamp > 0.0f) {
            lm.floorThreshold = static_cast<double>(config_.floorClamp);
        } else if (config_.floorStop > 0.0f) {
            lm.floorThreshold = static_cast<double>(config_.floorStop);
        }
        if (sf2(11)) { lm.screenLoss = static_cast<double>(sf2NullLoss()); }
        sf2LmRules(lm);
        auto recordLmStop = [&](LevenbergMarquardt const& l) {
            if (!l.resolutionStop) { return; }
            if (l.stopUnresolvableDamped) { ++nLmStopUnresolvableDamped_; }
            if (l.stopRealizedSmall) { ++nLmStopRealizedSmall_; }
            if (l.stopFixedPoint) { ++nLmStopFixedPoint_; }
            auto const bucket = static_cast<std::size_t>(std::clamp(l.itersRun, 0, 63));
            ++lmIterHist_[bucket];
        };
        auto recordGaugeDiag = [&](LevenbergMarquardt const& l) {
            if (!l.gaugeProject) { return; }
            for (int d : l.gaugeNullDimHistory) { op2GaugeRecordNullDim(d); }
        };
        Operon::Scalar bestv = 0;
        std::vector<Operon::Scalar> best;
        if (useDs) {
            if (lm.resolutionStop) { lm.dataCount = static_cast<double>(dataXDs_.size()); }
            auto buildDs = [&](std::vector<Operon::Scalar> const& v, std::vector<double>& r, std::vector<double>& J) -> bool { return buildOn(v, r, J, true); };
            auto lossDs = [&](std::vector<Operon::Scalar> const& v) -> Operon::Scalar { return lossTimedOn(v, true); };
            auto mid = lm.Minimize(c0, buildDs, lossDs, bestv);
            ++nConstDescents_;
            nConstLossCalls_ += static_cast<std::size_t>(lm.lossCalls);
            nConstIters_ += static_cast<std::size_t>(lm.itersRun);
            recordLmStop(lm);
            recordGaugeDiag(lm);
            LevenbergMarquardt lmPolish = lm;
            lmPolish.maxIter = std::max(0, config_.lmDownsamplePolishIters);
            if (lmPolish.resolutionStop) { lmPolish.dataCount = static_cast<double>(dataX_.size()); }
            Operon::Scalar bestv2 = bestv;
            best = (lmPolish.maxIter > 0) ? lmPolish.Minimize(mid, build, lossTimed, bestv2) : mid;
            nConstLossCalls_ += static_cast<std::size_t>(lmPolish.lossCalls);
            nConstIters_ += static_cast<std::size_t>(lmPolish.itersRun);
            if (lmPolish.maxIter > 0) { recordLmStop(lmPolish); recordGaugeDiag(lmPolish); }
        } else if (config_.useVarProLinear && useBc && progX.valid && progY.valid && c0.size() <= 64) {
            auto const maskX = computeAffineCoeffMask(progX, nx);
            auto const maskY = computeAffineCoeffMask(progY, c0.size() - nx);
            auto const linMask = maskX | (maskY << nx);
            best = lm.MinimizeVarPro(c0, build, lossTimed, bestv, linMask);
            ++nConstDescents_;
            nConstLossCalls_ += static_cast<std::size_t>(lm.lossCalls);
            nConstIters_ += static_cast<std::size_t>(lm.itersRun);
            recordLmStop(lm);
            recordGaugeDiag(lm);
        } else {
            best = lm.Minimize(c0, build, lossTimed, bestv);
            ++nConstDescents_;
            nConstLossCalls_ += static_cast<std::size_t>(lm.lossCalls);
            nConstIters_ += static_cast<std::size_t>(lm.itersRun);
            recordLmStop(lm);
            recordGaugeDiag(lm);
        }
        if (best.size() != c0.size()) { return; }
        for (auto const v : best) {
            if (!std::isfinite(static_cast<double>(v))) { ++nNonFiniteLmConst_; nSafetyCounters().nonFiniteLmConst.fetch_add(1, std::memory_order_relaxed); return; }
        }
        ind.x.SetCoefficients(std::vector<Operon::Scalar>(best.begin(), best.begin() + static_cast<long>(nx)));
        ind.y.SetCoefficients(std::vector<Operon::Scalar>(best.begin() + static_cast<long>(nx), best.end()));
        if (op2(2) || op2DiagEnabled()) {
            for (auto const v : best) { op2GaugeRecordAbsConst(std::abs(static_cast<double>(v))); }
        }
    }

    void optConstsGnBatchedPopulation(std::vector<PairIndividual*> const& inds,
                                      Operon::Dataset const& ds, Operon::Range range,
                                      bool fullObjective) const
    {
        std::size_t const P = inds.size();
        if (P == 0) { return; }
        auto const tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());

        struct Cand {
            PairIndividual* ind{nullptr};
            std::size_t nx{0};
            std::vector<Operon::Scalar> x;
            BytecodeProgram progX, progY;
            Sketch::RegisterProgram regX, regY;
            double lambda = 1e-3;
            Operon::Scalar f = 0;
            int outerIt = 0;
            bool haveResult = false;
        };

        std::vector<Cand> cs(P);
        std::vector<std::size_t> active;
        active.reserve(P);
        for (std::size_t i = 0; i < P; ++i) {
            auto& c = cs[i];
            c.ind = inds[i];
            auto cx = c.ind->x.GetCoefficients();
            auto cy = c.ind->y.GetCoefficients();
            c.nx = cx.size();
            c.x = cx;
            c.x.insert(c.x.end(), cy.begin(), cy.end());
            if (c.x.empty()) { continue; }
            c.progX = BytecodeProgram::Compile(c.ind->x, {cx.data(), cx.size()});
            c.progY = BytecodeProgram::Compile(c.ind->y, {cy.data(), cy.size()});
            if (!c.progX.valid || !c.progY.valid) { continue; }
            c.regX = Sketch::RegisterProgram::FromBytecode(c.progX);
            c.regY = Sketch::RegisterProgram::FromBytecode(c.progY);
            if (!c.regX.valid || !c.regY.valid) { continue; }
            active.push_back(i);
        }
        if (active.empty()) { return; }

        LevenbergMarquardt lmCfg;
        if (config_.floorClamp > 0.0f) { lmCfg.floorThreshold = static_cast<double>(config_.floorClamp); }
        else if (config_.floorStop > 0.0f) { lmCfg.floorThreshold = static_cast<double>(config_.floorStop); }
        if (sf2(11)) { lmCfg.screenLoss = static_cast<double>(sf2NullLoss()); }
        sf2LmRules(lmCfg);

        std::vector<Sketch::RegisterProgram> bx, by;
        std::vector<std::vector<Operon::Scalar>> gxs, gys;
        auto batchedLoss = [&](std::vector<std::size_t> const& idxs,
                                std::vector<std::vector<Operon::Scalar>> const& xs) -> std::vector<Operon::Scalar> {
            std::size_t const K = idxs.size();
            bx.resize(K); by.resize(K);
            for (std::size_t k = 0; k < K; ++k) {
                auto& c = cs[idxs[k]];
                c.regX.UpdateCoefficients({xs[k].data(), c.nx});
                c.regY.UpdateCoefficients({xs[k].data() + c.nx, xs[k].size() - c.nx});
                bx[k] = c.regX;
                by[k] = c.regY;
            }
            Sketch::BatchedRegisterInterp::EvaluatePopulation(bx, tValues, gxs);
            Sketch::BatchedRegisterInterp::EvaluatePopulation(by, tValues, gys);
            std::vector<Operon::Scalar> out(K);
            for (std::size_t k = 0; k < K; ++k) {
                auto sd = shapeDistance(gxs[k], gys[k], dataX_, dataY_, noise_, nullptr, 4,
                                        &dataBlocks_, config_.useGridSdf ? &gridSdf_ : nullptr,
                                        config_.useAmxChamfer ? &amxData_ : nullptr);
                if (fullObjective && config_.orderWeight > Operon::Scalar{0}) { sd += orderTerm(gxs[k], gys[k]); }
                if (fullObjective && config_.curvatureWeight > Operon::Scalar{0}) { sd = applyCurvature(sd, curvatureTerm(gxs[k], gys[k])); }
                out[k] = sd;
            }
            return out;
        };

        {
            std::vector<std::size_t> idxs(active);
            std::vector<std::vector<Operon::Scalar>> xs(idxs.size());
            for (std::size_t k = 0; k < idxs.size(); ++k) { xs[k] = cs[idxs[k]].x; }
            auto fs = batchedLoss(idxs, xs);
            std::vector<std::size_t> next;
            next.reserve(active.size());
            for (std::size_t k = 0; k < idxs.size(); ++k) {
                auto& c = cs[idxs[k]];
                c.f = fs[k];
                c.haveResult = true;
                if (!std::isfinite(static_cast<double>(c.f))) { continue; }
                if (lmCfg.floorThreshold > 0.0 && static_cast<double>(c.f) <= lmCfg.floorThreshold) { continue; }
                next.push_back(idxs[k]);
            }
            active.swap(next);
        }

        for (int it = 0; it < lmCfg.maxIter && !active.empty(); ++it) {
            struct IterState {
                std::size_t ci;
                std::vector<double> H, g, diag;
                bool accepted = false;
                bool unresolvable = false;
                double rel = 0.0;
                double fOld = 0.0;
                std::vector<Operon::Scalar> xBefore;
                std::size_t m = 0;
            };
            std::vector<IterState> iters;
            iters.reserve(active.size());
            std::vector<double> r, J;
            for (auto ci : active) {
                auto& c = cs[ci];
                r.clear(); J.clear();
                BytecodeProgram* pxc = (config_.useLmBytecodeCache && c.progX.valid) ? &c.progX : nullptr;
                BytecodeProgram* pyc = (config_.useLmBytecodeCache && c.progY.valid) ? &c.progY : nullptr;
                bool buildOk;
                if (config_.profileTiming) {
                    auto tb0 = std::chrono::steady_clock::now();
                    buildOk = lmBuild(*c.ind, c.x, c.nx, ds, range, r, J, nullptr, pxc, pyc, nullptr, nullptr, nullptr, nullptr);
                    auto tb1 = std::chrono::steady_clock::now();
                    profLmBuildNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(tb1 - tb0).count(), std::memory_order_relaxed);
                } else {
                    buildOk = lmBuild(*c.ind, c.x, c.nx, ds, range, r, J, nullptr, pxc, pyc, nullptr, nullptr, nullptr, nullptr);
                }
                if (!buildOk) { continue; }
                auto const n = c.x.size();
                auto const m = r.size();
                if (m == 0 || J.size() != m * n) { continue; }
                IterState st; st.ci = ci; st.m = m;
                st.H.assign(n * n, 0.0); st.g.assign(n, 0.0); st.diag.assign(n, 0.0);
                for (std::size_t i = 0; i < m; ++i) {
                    double const* row = J.data() + i * n;
                    double const ri = r[i];
                    for (std::size_t a = 0; a < n; ++a) {
                        double const ra = row[a];
                        st.g[a] += ra * ri;
                        double* Ha = st.H.data() + a * n;
                        for (std::size_t b = a; b < n; ++b) { Ha[b] += ra * row[b]; }
                    }
                }
                for (std::size_t a = 0; a < n; ++a) {
                    for (std::size_t b = 0; b < a; ++b) { st.H[a * n + b] = st.H[b * n + a]; }
                    st.diag[a] = st.H[a * n + a];
                    if (!(st.diag[a] > 0.0)) { st.diag[a] = 1.0; }
                }
                st.fOld = static_cast<double>(c.f);
                st.xBefore = c.x;
                iters.push_back(std::move(st));
            }

            std::vector<char> trying(iters.size(), 1);
            static thread_local std::vector<double> L_scratch, y_scratch, dx_scratch;
            for (int trial = 0; lmCfg.resolutionStop || trial < lmCfg.maxTrials; ++trial) {
                if (std::none_of(trying.begin(), trying.end(), [](char t) { return t != 0; })) { break; }
                std::vector<std::size_t> trialCi;
                std::vector<std::vector<Operon::Scalar>> trialX;
                for (std::size_t ii = 0; ii < iters.size(); ++ii) {
                    if (!trying[ii]) { continue; }
                    auto& c = cs[iters[ii].ci];
                    auto const n = c.x.size();
                    if (lmCfg.resolutionStop && !std::isfinite(c.lambda)) {
                        iters[ii].unresolvable = true; trying[ii] = 0; continue;
                    }
                    if (!LevenbergMarquardt::solveDamped(iters[ii].H, iters[ii].g, iters[ii].diag, c.lambda, n,
                                           dx_scratch, L_scratch, y_scratch)) {
                        c.lambda = lmCfg.resolutionStop ? (c.lambda * lmCfg.lambdaUp)
                                                         : std::min(c.lambda * lmCfg.lambdaUp, lmCfg.lambdaMax);
                        continue;
                    }
                    if (lmCfg.gaugeProject) {
                        op2GaugeRecordNullDim(LevenbergMarquardt::projectGaugeNull(
                            iters[ii].H, n, lmCfg.jacobianEps, iters[ii].m, dx_scratch));
                    }
                    if (lmCfg.resolutionStop) {
                        double gdx = 0.0;
                        for (std::size_t a = 0; a < n; ++a) { gdx += iters[ii].g[a] * dx_scratch[a]; }
                        double dHd = 0.0;
                        for (std::size_t a = 0; a < n; ++a) {
                            double const* Ha = iters[ii].H.data() + a * n;
                            double hd = 0.0;
                            for (std::size_t b = 0; b < n; ++b) { hd += Ha[b] * dx_scratch[b]; }
                            dHd += dx_scratch[a] * hd;
                        }
                        double const pred = -2.0 * gdx - dHd;
                        double const fCur = static_cast<double>(c.f);
                        bool resolvable;
                        if (pred <= 0.0) { resolvable = false; }
                        else if (fCur - pred <= 0.0) { resolvable = true; }
                        else { resolvable = ((lmCfg.dataCount / 2.0) * std::log(fCur / (fCur - pred)) >= 0.5); }
                        if (!resolvable) { iters[ii].unresolvable = true; trying[ii] = 0; continue; }
                    }
                    std::vector<Operon::Scalar> xn = c.x;
                    bool finite = true;
                    for (std::size_t a = 0; a < n; ++a) {
                        auto const v = static_cast<double>(c.x[a]) + dx_scratch[a];
                        if (!std::isfinite(v)) { finite = false; break; }
                        xn[a] = static_cast<Operon::Scalar>(v);
                    }
                    if (!finite) {
                        c.lambda = lmCfg.resolutionStop ? (c.lambda * lmCfg.lambdaUp)
                                                         : std::min(c.lambda * lmCfg.lambdaUp, lmCfg.lambdaMax);
                        continue;
                    }
                    trialCi.push_back(ii);
                    trialX.push_back(std::move(xn));
                }
                if (trialCi.empty()) { continue; }
                std::vector<std::size_t> candIdxs(trialCi.size());
                for (std::size_t t = 0; t < trialCi.size(); ++t) { candIdxs[t] = iters[trialCi[t]].ci; }
                auto fns = batchedLoss(candIdxs, trialX);
                for (std::size_t t = 0; t < trialCi.size(); ++t) {
                    auto ii = trialCi[t];
                    auto& c = cs[iters[ii].ci];
                    auto const fn = fns[t];
                    if (std::isfinite(static_cast<double>(fn)) && fn < c.f) {
                        iters[ii].rel = (static_cast<double>(c.f) - static_cast<double>(fn))
                            / std::max(std::abs(static_cast<double>(c.f)), 1e-30);
                        c.x = std::move(trialX[t]);
                        c.f = fn;
                        c.lambda = std::max(c.lambda * lmCfg.lambdaDown, lmCfg.lambdaMin);
                        iters[ii].accepted = true;
                        trying[ii] = 0;
                        continue;
                    }
                    c.lambda = lmCfg.resolutionStop ? (c.lambda * lmCfg.lambdaUp)
                                                     : std::min(c.lambda * lmCfg.lambdaUp, lmCfg.lambdaMax);
                    if (!lmCfg.resolutionStop && lmCfg.enableScreening && it == 0 && trial == 0
                        && static_cast<double>(c.f) > lmCfg.screenLoss) {
                        trying[ii] = 0;
                    }
                }
            }

            std::vector<std::size_t> next;
            next.reserve(active.size());
            for (auto const& st : iters) {
                if (st.unresolvable) { ++nLmStopUnresolvableDamped_; continue; }
                if (!st.accepted) {
                    if (lmCfg.resolutionStop) { ++nLmStopUnresolvableDamped_; }
                    continue;
                }
                if (lmCfg.resolutionStop) {
                    auto const& c = cs[st.ci];
                    bool same = (c.x.size() == st.xBefore.size());
                    for (std::size_t a = 0; same && a < c.x.size(); ++a) { same = same && (c.x[a] == st.xBefore[a]); }
                    if (same) { ++nLmStopFixedPoint_; continue; }
                    auto const& gain = (lmCfg.dataCount / 2.0) * std::log(st.fOld / static_cast<double>(c.f));
                    if (gain < 0.5) { ++nLmStopRealizedSmall_; continue; }
                } else if (st.rel < lmCfg.ftol) { continue; }
                next.push_back(st.ci);
            }
            active.swap(next);
        }

        for (auto& c : cs) {
            if (!c.haveResult || c.x.empty() || c.x.size() < c.nx) { continue; }
            c.ind->x.SetCoefficients(std::vector<Operon::Scalar>(c.x.begin(), c.x.begin() + static_cast<long>(c.nx)));
            c.ind->y.SetCoefficients(std::vector<Operon::Scalar>(c.x.begin() + static_cast<long>(c.nx), c.x.end()));
            ++nConstDescents_;
        }
    }

    [[nodiscard]] auto frequencySeedingTopologyAllowed() const -> bool
    {
        if (!config_.useTopologyGatedFrequency) { return true; }
        if (topoClosedness_ < 0.0f) {
            std::size_t const n = dataX_.size();
            if (n >= 2 && dataY_.size() == n) {
                std::vector<float> fx(n), fy(n);
                for (std::size_t i = 0; i < n; ++i) {
                    fx[i] = static_cast<float>(dataX_[i]);
                    fy[i] = static_cast<float>(dataY_[i]);
                }
                topoClosedness_ = Sketch::strokeClosedness(fx.data(), fy.data(), static_cast<int>(n));
            } else {
                topoClosedness_ = 0.0f;
            }
        }
        return topoClosedness_ > (1.0f - Sketch::kClosedGapFraction);
    }

    void optConsts(PairIndividual& ind, Operon::Dataset const& ds, Operon::Range range,
                   bool fullObjective = false) const
    {
#ifdef SKETCH_DUP_STATS
        {
            std::uint64_t hs = DupStats::kFnvBasis;
            DupStats::mixTree(hs, ind.x); DupStats::mixTree(hs, ind.y);
            std::uint64_t hx = hs;
            auto dcx = ind.x.GetCoefficients(); auto dcy = ind.y.GetCoefficients();
            DupStats::mixCoefs(hx, dcx); DupStats::mixCoefs(hx, dcy);
            dupStats_.record(hs, hx, !(dcx.empty() && dcy.empty()));
        }
#endif
        if (config_.useLm) { optConstsLm(ind, ds, range); return; }
        if (config_.constOptimizer == 1) { optConstsGn(ind, ds, range, fullObjective); return; }
        auto cx = ind.x.GetCoefficients();
        auto cy = ind.y.GetCoefficients();
        std::size_t nx = cx.size();
        std::vector<Operon::Scalar> c0 = cx;
        c0.insert(c0.end(), cy.begin(), cy.end());
        if (c0.empty()) { return; }
        BytecodeProgram progX, progY;
        bool const useBc = config_.useBytecodeEvaluator;
        if (useBc) {
            progX = BytecodeProgram::Compile(ind.x, {cx.data(), cx.size()});
            progY = BytecodeProgram::Compile(ind.y, {cy.data(), cy.size()});
        }
        auto const tValues = ds.GetValues(static_cast<int64_t>(0)).subspan(range.Start(), range.Size());

        auto loss = [&](std::vector<Operon::Scalar> const& v) -> Operon::Scalar {
            thread_local std::vector<Operon::Scalar> gx, gy;
            if (useBc && progX.valid && progY.valid) {
                gx.resize(range.Size());
                gy.resize(range.Size());
                progX.UpdateCoefficients({v.data(), nx});
                progY.UpdateCoefficients({v.data() + nx, v.size() - nx});
                progX.Evaluate(tValues, gx);
                progY.Evaluate(tValues, gy);
            } else {
                std::vector<Operon::Scalar> vx(v.begin(), v.begin() + nx);
                std::vector<Operon::Scalar> vy(v.begin() + nx, v.end());
                interpInto(ind.x, vx, ds, range, gx);
                interpInto(ind.y, vy, ds, range, gy);
            }
            auto sd = shapeDistance(gx, gy, dataX_, dataY_, noise_, nullptr, 4, &dataBlocks_, config_.useGridSdf ? &gridSdf_ : nullptr, config_.useAmxChamfer ? &amxData_ : nullptr);
            if (fullObjective && config_.orderWeight > Operon::Scalar{0}) { sd += orderTerm(gx, gy); }
            if (fullObjective && config_.curvatureWeight > Operon::Scalar{0}) { sd = applyCurvature(sd, curvatureTerm(gx, gy)); }
            ++nConstLossCalls_;
            return sd;
        };
        ++nConstDescents_;
        NelderMead nm;
        if (fullObjective) {
            nm.maxIter = (config_.polishConstIter > 0) ? config_.polishConstIter : static_cast<int>(30 * c0.size());
        } else {
            nm.maxIter = (config_.inLoopConstIter > 0) ? config_.inLoopConstIter : static_cast<int>(30 * c0.size());
        }
        Operon::Scalar bestv;
        auto best = nm.Minimize(c0, loss, bestv);

        if (config_.constMultiStart > 1) {
            auto mx = frequencyParamMask(ind.x);
            auto my = frequencyParamMask(ind.y);
            if (mx.size() == nx && my.size() == cy.size()) {
                std::vector<char> freq;
                freq.reserve(c0.size());
                freq.insert(freq.end(), mx.begin(), mx.end());
                freq.insert(freq.end(), my.begin(), my.end());
                auto const nFreq = std::count(freq.begin(), freq.end(), char{1});
                if (nFreq > 0) {
                    constexpr double TWO_PI = 6.283185307179586;
                    std::vector<long> turns;
                    if (config_.constMultiStartDataSeed) {
                        turns = dataSeedTurnCounts(turningNumber_, config_.constMultiStart - 1);
                    }
                    for (std::size_t j = 1; j < config_.constMultiStart; ++j) {
                        auto start = c0;
                        bool seeded = false;
                        if (config_.constMultiStartDataSeed && (j - 1) < turns.size()) {
                            double const omega = TWO_PI * static_cast<double>(turns[j - 1]);
                            for (std::size_t k = 0; k < start.size(); ++k) {
                                if (freq[k] != 0) {
                                    double const s = (c0[k] < Operon::Scalar{0}) ? -1.0 : 1.0;
                                    start[k] = static_cast<Operon::Scalar>(s * omega);
                                }
                            }
                            seeded = true;
                        }
                        if (!seeded) {
                            for (std::size_t k = 0; k < start.size(); ++k) {
                                if (freq[k] != 0) { start[k] = freqRestartValue(c0[k], j); }
                            }
                        }
                        Operon::Scalar v = 0;
                        auto cand = nm.Minimize(start, loss, v);
                        if (std::isfinite(v) && v < bestv) { bestv = v; best = std::move(cand); }
                    }

                    if (config_.useFrequencyDiscovery && frequencySeedingTopologyAllowed()) {
                        std::vector<Operon::Scalar> bx(best.begin(), best.begin() + static_cast<long>(nx));
                        std::vector<Operon::Scalar> by(best.begin() + static_cast<long>(nx), best.end());
                        std::vector<Operon::Scalar> fgx, fgy;
                        interpInto(ind.x, bx, ds, range, fgx);
                        interpInto(ind.y, by, ds, range, fgy);
                        std::size_t const m = range.Size();
                        std::vector<double> resX(m), resY(m), tD(m);
                        bool finiteAll = (fgx.size() == m && fgy.size() == m);
                        for (std::size_t i = 0; finiteAll && i < m; ++i) {
                            double const rx = static_cast<double>(dataX_[i]) - static_cast<double>(fgx[i]);
                            double const ry = static_cast<double>(dataY_[i]) - static_cast<double>(fgy[i]);
                            if (!std::isfinite(rx) || !std::isfinite(ry)) { finiteAll = false; break; }
                            resX[i] = rx; resY[i] = ry;
                            tD[i] = static_cast<double>(tValues[i]);
                        }
                        if (finiteAll && m > 0) {
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
                            freqdisc::FrequencyProposal propX, propY;
                            if (config_.useMatrixPencilFrequencySeed) {
                                freqdisc::proposeFrequencyPencilXY(resX, resY, useSpan, propX, propY);
                            } else {
                                propX = config_.useNoiseWeightedFrequencySeed
                                    ? freqdisc::proposeFrequencyWeighted(resX, useSpan, seedSigma)
                                    : freqdisc::proposeFrequency(resX, useSpan);
                                propY = config_.useNoiseWeightedFrequencySeed
                                    ? freqdisc::proposeFrequencyWeighted(resY, useSpan, seedSigma)
                                    : freqdisc::proposeFrequency(resY, useSpan);
                            }
                            freqdisc::FrequencyProposal const* prop = nullptr;
                            bool useX = true;
                            if (propX.valid && propY.valid) {
                                useX = propX.magnitude >= propY.magnitude;
                                prop = useX ? &propX : &propY;
                            } else if (propX.valid) { prop = &propX; useX = true; }
                            else if (propY.valid) { prop = &propY; useX = false; }
                            if (prop != nullptr) {
                                auto const& targetRes = useX ? resX : resY;
                                double seedOmega = prop->omega;
                                auto qf = freqdisc::quadratureFit(tD, targetRes, seedOmega);
                                if (qf.valid && config_.useVarProFrequencyGn) {
                                    std::vector<double> variance;
                                    if (config_.useWeightedVarProFrequencyGn && !noise_.empty()) {
                                        variance.resize(m);
                                        for (std::size_t i = 0; i < m; ++i) {
                                            auto const region = noiseIndexAtArc(i, m, noise_.size());
                                            variance[i] = static_cast<double>(noise_.s2[region]);
                                        }
                                    }
                                    auto const refined = freqdisc::refineQuadratureFrequencyOneStep(
                                        tD, targetRes, seedOmega, variance);
                                    if (refined.valid) {
                                        seedOmega = refined.omega;
                                        qf = refined.fit;
                                    }
                                }
                                if (qf.valid) {
                                    auto start = c0;
                                    for (std::size_t k = 0; k < start.size(); ++k) {
                                        if (freq[k] != 0) {
                                            double const s = (c0[k] < Operon::Scalar{0}) ? -1.0 : 1.0;
                                            start[k] = static_cast<Operon::Scalar>(s * seedOmega);
                                        }
                                    }
                                    Operon::Scalar v2 = 0;
                                    auto cand2 = nm.Minimize(start, loss, v2);
                                    if (std::isfinite(v2) && v2 < bestv) { bestv = v2; best = std::move(cand2); }
                                }
                            }
                        }
                    }
                }
            }
        }

        ind.x.SetCoefficients(std::vector<Operon::Scalar>(best.begin(), best.begin() + nx));
        ind.y.SetCoefficients(std::vector<Operon::Scalar>(best.begin() + nx, best.end()));
    }

    auto buildFront(std::vector<PairIndividual*>& members, Operon::Dataset const& ds,
                    Operon::Range range) const -> std::vector<FrontEntry>
    {
        std::sort(members.begin(), members.end(), [](auto* a, auto* b) {
            return a->complexity < b->complexity
                || (a->complexity == b->complexity && a->shapeErr < b->shapeErr);
        });
        std::vector<FrontEntry> out;
        out.reserve(members.size());
        for (auto* p : members) {
            FrontEntry e;
            try {
                e.exprX = Operon::InfixFormatter::Format(p->x, ds, 4);
                e.exprY = Operon::InfixFormatter::Format(p->y, ds, 4);
                e.shapeErr = p->shapeErr;
                if (config_.directFinalistEval) {
                    DirectFinalistEvalGuard dfe(true);
                    Operon::Scalar chamfer{0};
                    e.shapeErr = shapeOnlyParts(*p, ds, range, &chamfer);
                }
                e.nodes = static_cast<int>(p->totalLength());
                e.errX = p->errX; e.errY = p->errY;
                e.regionErr = regionResid(*p, ds, range);
                auto [scx, scy] = sampledCurve(*p, ds, range);
                e.curveX = std::move(scx); e.curveY = std::move(scy);
                if (!e.curveX.empty()) {
                    e.orderErr = orderDistance(e.curveX, e.curveY, dataX_, dataY_, dataClosed_, 64, &dataArcX_, &dataArcY_);
                    e.turning = turningFourier(e.curveX, e.curveY, dataClosed_, config_.curvatureHarmonics);
                }
            } catch (std::exception const&) {
                ++nEvalErrors_;
            } catch (...) {
                ++nEvalErrors_;
            }
            out.push_back(std::move(e));
        }
        return out;
    }

    static constexpr std::size_t kArchiveCapMultiple = 2;
    static constexpr std::size_t kArchiveReinjectDivisor = 20;

    using ArchiveKey = std::array<Operon::Scalar, 3>;
    static auto archiveKey(PairIndividual const& p) -> ArchiveKey {
        if (perTreeCreditMode()) { return {p.errX, p.errY, p.complexity}; }
        return {p.shapeErr, Operon::Scalar{0}, p.complexity};
    }

    static auto archiveEligible(PairIndividual const& p) -> bool {
        constexpr auto BIG = std::numeric_limits<Operon::Scalar>::max();
        auto ok = [](Operon::Scalar v) { return std::isfinite(v) && v < BIG; };
        if (!ok(p.shapeErr)) { return false; }
        if (perTreeCreditMode()) { return ok(p.errX) && ok(p.errY); }
        return true;
    }

    static void archiveUpdate(std::vector<PairIndividual>& archive,
                              std::vector<PairIndividual> const& candidates,
                              std::size_t cap)
    {
        std::set<ArchiveKey> seen;
        for (auto const& p : archive) { seen.insert(archiveKey(p)); }
        for (auto const& c : candidates) {
            if (!archiveEligible(c)) { continue; }
            if (!seen.insert(archiveKey(c)).second) { continue; }
            archive.push_back(c);
        }
        if (archive.empty()) { return; }

        nonDominatedSort(archive);
        archive.erase(std::remove_if(archive.begin(), archive.end(),
                                     [](PairIndividual const& p) { return p.rank != 0; }),
                      archive.end());

        if (archive.size() <= cap) { return; }
        std::vector<std::size_t> idx(archive.size());
        std::iota(idx.begin(), idx.end(), std::size_t{0});
        assignCrowding(archive, idx);
        std::vector<std::size_t> order(archive.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return archive[a].crowding > archive[b].crowding;
        });
        order.resize(cap);
        std::sort(order.begin(), order.end());
        std::vector<PairIndividual> kept;
        kept.reserve(cap);
        for (auto i : order) { kept.push_back(std::move(archive[i])); }
        archive = std::move(kept);
    }

    static constexpr std::size_t kHallOfFameK = 50;

    static auto hofBucket(PairIndividual const& p) -> long long {
        return static_cast<long long>(p.totalLength());
    }

    static void hallOfFameUpdate(std::vector<PairIndividual>& hof,
                                 std::vector<PairIndividual> const& candidates,
                                 std::size_t kPerBucket)
    {
        std::set<ArchiveKey> seen;
        for (auto const& p : hof) { seen.insert(archiveKey(p)); }
        for (auto const& c : candidates) {
            if (!archiveEligible(c)) { continue; }
            if (!seen.insert(archiveKey(c)).second) { continue; }
            hof.push_back(c);
        }
        if (hof.empty()) { return; }

        std::stable_sort(hof.begin(), hof.end(), [](PairIndividual const& a, PairIndividual const& b) {
            auto const ba = hofBucket(a), bb = hofBucket(b);
            if (ba != bb) { return ba < bb; }
            return a.shapeErr < b.shapeErr;
        });
        std::vector<PairIndividual> kept;
        kept.reserve(std::min(hof.size(), kPerBucket * 64));
        long long curBucket = std::numeric_limits<long long>::min();
        std::size_t inBucket = 0;
        for (auto& p : hof) {
            auto const b = hofBucket(p);
            if (b != curBucket) { curBucket = b; inBucket = 0; }
            if (inBucket < kPerBucket) { kept.push_back(std::move(p)); ++inBucket; }
        }
        hof = std::move(kept);
    }

    static void archiveReinject(std::vector<PairIndividual>& pop,
                                std::vector<PairIndividual> const& archive,
                                std::mt19937& sel)
    {
        if (pop.empty() || archive.empty()) { return; }
        std::set<ArchiveKey> inPop;
        for (auto const& p : pop) { inPop.insert(archiveKey(p)); }
        std::vector<std::size_t> cand;
        for (std::size_t i = 0; i < archive.size(); ++i) {
            if (inPop.find(archiveKey(archive[i])) == inPop.end()) { cand.push_back(i); }
        }
        if (cand.empty()) { return; }
        std::size_t const K = std::max<std::size_t>(1, pop.size() / kArchiveReinjectDivisor);
        std::size_t const n = std::min({K, cand.size(), pop.size()});
        std::shuffle(cand.begin(), cand.end(), sel);
        std::vector<std::size_t> order(pop.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return crowdedLess(pop[a], pop[b]);
        });
        for (std::size_t k = 0; k < n; ++k) {
            pop[order[pop.size() - 1 - k]] = archive[cand[k]];
        }
    }

    static void fastNonDominatedSort2Obj(std::vector<PairIndividual>& pop) {
        auto const n = pop.size();
        if (n == 0) { return; }
        std::vector<std::pair<Operon::Scalar, Operon::Scalar>> key(n);
        for (std::size_t i = 0; i < n; ++i) {
            key[i] = {clampedShapeErr(pop[i]), effectiveComplexity(pop[i])};
        }
        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            if (key[a].first != key[b].first) { return key[a].first < key[b].first; }
            return key[a].second < key[b].second;
        });
        std::vector<Operon::Scalar> compVals(n);
        for (std::size_t i = 0; i < n; ++i) { compVals[i] = key[i].second; }
        std::sort(compVals.begin(), compVals.end());
        compVals.erase(std::unique(compVals.begin(), compVals.end()), compVals.end());
        auto compIndex = [&](Operon::Scalar c) -> std::size_t {
            return static_cast<std::size_t>(
                std::lower_bound(compVals.begin(), compVals.end(), c) - compVals.begin());
        };
        std::size_t const m = compVals.size();
        std::vector<int> fenwick(m + 1, 0);
        auto fenwickQueryPrefixMax = [&](std::size_t idx1) -> int {
            int best = 0;
            for (auto i = idx1; i > 0; i -= (i & (~i + 1))) { best = std::max(best, fenwick[i]); }
            return best;
        };
        auto fenwickUpdateMax = [&](std::size_t idx1, int value) {
            for (auto i = idx1; i <= m; i += (i & (~i + 1))) {
                fenwick[i] = std::max(fenwick[i], value);
            }
        };
        std::size_t i = 0;
        while (i < n) {
            std::size_t j = i;
            while (j < n && key[order[j]].first == key[order[i]].first
                          && key[order[j]].second == key[order[i]].second) {
                ++j;
            }
            int groupRank = 0;
            for (std::size_t k = i; k < j; ++k) {
                auto const ci = compIndex(key[order[k]].second);
                groupRank = std::max(groupRank, fenwickQueryPrefixMax(ci + 1));
            }
            for (std::size_t k = i; k < j; ++k) {
                pop[order[k]].rank = groupRank;
            }
            for (std::size_t k = i; k < j; ++k) {
                auto const ci = compIndex(key[order[k]].second);
                fenwickUpdateMax(ci + 1, groupRank + 1);
            }
            i = j;
        }
    }

    static void nonDominatedSort(std::vector<PairIndividual>& pop) {
        if (useExactNlogNSortMode() && !perTreeCreditMode()) {
            fastNonDominatedSort2Obj(pop);
            return;
        }
        auto const n = pop.size();
        std::vector<std::vector<std::size_t>> dominated(n);
        std::vector<int> domCount(n, 0);
        std::vector<std::size_t> front;
        if (fastNonDominatedSortMode()) {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = i + 1; j < n; ++j) {
                    if (dominates(pop[i], pop[j])) {
                        dominated[i].push_back(j);
                        ++domCount[j];
                    } else if (dominates(pop[j], pop[i])) {
                        dominated[j].push_back(i);
                        ++domCount[i];
                    }
                }
            }
            for (std::size_t i = 0; i < n; ++i) {
                if (domCount[i] == 0) { pop[i].rank = 0; front.push_back(i); }
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    if (i == j) { continue; }
                    if (dominates(pop[i], pop[j])) { dominated[i].push_back(j); }
                    else if (dominates(pop[j], pop[i])) { ++domCount[i]; }
                }
                if (domCount[i] == 0) { pop[i].rank = 0; front.push_back(i); }
            }
        }
        int r = 0;
        while (!front.empty()) {
            std::vector<std::size_t> next;
            for (auto i : front) {
                for (auto j : dominated[i]) {
                    if (--domCount[j] == 0) { pop[j].rank = r + 1; next.push_back(j); }
                }
            }
            ++r;
            front = std::move(next);
        }
    }

    static void assignCrowding(std::vector<PairIndividual>& pop,
                               std::vector<std::size_t>& idx) {
        constexpr auto INF = std::numeric_limits<Operon::Scalar>::max();
        for (auto i : idx) { pop[i].crowding = 0; }
        if (idx.size() <= 2) { for (auto i : idx) { pop[i].crowding = INF; } return; }
        auto crowd1 = [&](auto keyFn) {
            std::sort(idx.begin(), idx.end(),
                      [&](std::size_t a, std::size_t b) { return keyFn(a) < keyFn(b); });
            pop[idx.front()].crowding = INF;
            pop[idx.back()].crowding = INF;
            auto span = keyFn(idx.back()) - keyFn(idx.front());
            if (!(span > 0)) { return; }
            for (std::size_t k = 1; k + 1 < idx.size(); ++k) {
                if (pop[idx[k]].crowding == INF) { continue; }
                pop[idx[k]].crowding += (keyFn(idx[k+1]) - keyFn(idx[k-1])) / span;
            }
        };
        if (perTreeCreditMode()) {
            crowd1([&](std::size_t i) { return clampedErr(pop[i].errX); });
            crowd1([&](std::size_t i) { return clampedErr(pop[i].errY); });
        } else {
            crowd1([&](std::size_t i) { return clampedShapeErr(pop[i]); });
        }
        crowd1([&](std::size_t i) { return effectiveComplexity(pop[i]); });
    }

    static void rankAndCrowd(std::vector<PairIndividual>& pop) {
        nonDominatedSort(pop);
        int maxRank = 0;
        for (auto const& p : pop) { maxRank = std::max(maxRank, p.rank); }
        for (int r = 0; r <= maxRank; ++r) {
            std::vector<std::size_t> idx;
            for (std::size_t i = 0; i < pop.size(); ++i) { if (pop[i].rank == r) { idx.push_back(i); } }
            assignCrowding(pop, idx);
        }
    }

    static void selectNextInto(std::vector<PairIndividual>& combined, std::size_t N,
                               std::vector<PairIndividual>& out,
                               int64_t* ndsNsOut = nullptr, int64_t* crowdNsOut = nullptr)
    {
        if (ndsNsOut != nullptr) {
            auto const t0 = std::chrono::steady_clock::now();
            nonDominatedSort(combined);
            *ndsNsOut += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t0).count();
        } else {
            nonDominatedSort(combined);
        }
        int maxRank = 0;
        for (auto const& p : combined) { maxRank = std::max(maxRank, p.rank); }
        out.clear();
        out.reserve(N);
        std::vector<std::size_t> idx;
        std::chrono::steady_clock::time_point crowdT0;
        for (int r = 0; r <= maxRank && out.size() < N; ++r) {
            idx.clear();
            for (std::size_t i = 0; i < combined.size(); ++i) {
                if (combined[i].rank == r) { idx.push_back(i); }
            }
            if (crowdNsOut != nullptr) { crowdT0 = std::chrono::steady_clock::now(); }
            assignCrowding(combined, idx);
            if (crowdNsOut != nullptr) {
                *crowdNsOut += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - crowdT0).count();
            }
            if (out.size() + idx.size() <= N) {
                for (auto i : idx) { out.push_back(std::move(combined[i])); }
            } else {
                if (crowdNsOut != nullptr) { crowdT0 = std::chrono::steady_clock::now(); }
                std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
                    return combined[a].crowding > combined[b].crowding;
                });
                if (crowdNsOut != nullptr) {
                    *crowdNsOut += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - crowdT0).count();
                }
                for (std::size_t k = 0; out.size() < N; ++k) {
                    out.push_back(std::move(combined[idx[k]]));
                }
            }
        }
    }

    static auto selectNext(std::vector<PairIndividual>& combined, std::size_t N)
        -> std::vector<PairIndividual>
    {
        std::vector<PairIndividual> out;
        selectNextInto(combined, N, out);
        return out;
    }

public:
    static void nonDominatedSortForTest(std::vector<PairIndividual>& pop) { nonDominatedSort(pop); }

    static auto islandRanges(std::size_t N, std::size_t k)
        -> std::vector<std::pair<std::size_t, std::size_t>>
    {
        if (k < 1) { k = 1; }
        if (N > 0) { k = std::min(k, N); }
        std::vector<std::pair<std::size_t, std::size_t>> r;
        r.reserve(k);
        std::size_t const base = N / k;
        std::size_t const rem = N % k;
        std::size_t lo = 0;
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t const sz = base + (i < rem ? 1 : 0);
            r.emplace_back(lo, lo + sz);
            lo += sz;
        }
        return r;
    }

    static void migrateRing(std::vector<PairIndividual>& pop,
                            std::vector<std::pair<std::size_t, std::size_t>> const& ranges,
                            std::size_t mCount)
    {
        std::size_t const k = ranges.size();
        if (k < 2 || mCount == 0) { return; }
        std::vector<std::vector<PairIndividual>> emigrants(k);
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t const lo = ranges[i].first, hi = ranges[i].second;
            std::size_t const sz = hi - lo;
            std::size_t const m = std::min(mCount, sz);
            std::vector<std::size_t> idx(sz);
            for (std::size_t j = 0; j < sz; ++j) { idx[j] = lo + j; }
            std::partial_sort(idx.begin(), idx.begin() + static_cast<std::ptrdiff_t>(m), idx.end(),
                [&](std::size_t a, std::size_t b) { return crowdedLess(pop[a], pop[b]); });
            emigrants[i].reserve(m);
            for (std::size_t j = 0; j < m; ++j) { emigrants[i].push_back(pop[idx[j]]); }
        }
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t const dst = (i + 1) % k;
            std::size_t const lo = ranges[dst].first, hi = ranges[dst].second;
            std::size_t const sz = hi - lo;
            std::size_t const m = std::min(emigrants[i].size(), sz);
            std::vector<std::size_t> idx(sz);
            for (std::size_t j = 0; j < sz; ++j) { idx[j] = lo + j; }
            std::partial_sort(idx.begin(), idx.begin() + static_cast<std::ptrdiff_t>(m), idx.end(),
                [&](std::size_t a, std::size_t b) { return crowdedLess(pop[b], pop[a]); });
            for (std::size_t j = 0; j < m; ++j) { pop[idx[j]] = emigrants[i][j]; }
        }
    }

private:
    auto tournament(std::vector<PairIndividual> const& pop, std::mt19937& sel) const
        -> PairIndividual const&
    {
        std::uniform_int_distribution<std::size_t> pick(0, pop.size() - 1);
        std::size_t k = std::max<std::size_t>(2, config_.tournamentSize);
        std::size_t best = pick(sel);
        for (std::size_t i = 1; i < k; ++i) {
            auto c = pick(sel);
            if (crowdedLess(pop[c], pop[best])) { best = c; }
        }
        return pop[best];
    }

    auto tournament(std::vector<PairIndividual> const& pop, std::mt19937& sel,
                    std::size_t lo, std::size_t hi) const
        -> PairIndividual const&
    {
        std::uniform_int_distribution<std::size_t> pick(lo, hi - 1);
        std::size_t k = std::max<std::size_t>(2, config_.tournamentSize);
        std::size_t best = pick(sel);
        for (std::size_t i = 1; i < k; ++i) {
            auto c = pick(sel);
            if (crowdedLess(pop[c], pop[best])) { best = c; }
        }
        return pop[best];
    }

    auto epsilonLexicaseSelect(std::vector<PairIndividual> const& pop, std::mt19937& sel) const
        -> PairIndividual const&
    {
        return epsilonLexicaseSelect(pop, sel, 0, pop.size());
    }

    auto epsilonLexicaseSelect(std::vector<PairIndividual> const& pop, std::mt19937& sel,
                               std::size_t lo, std::size_t hi) const
        -> PairIndividual const&
    {
        constexpr auto INF = std::numeric_limits<Operon::Scalar>::max();
        constexpr std::size_t kPoolSize = 12;
        bool const sf2Lex = sf2(5) && !sf2LexEps_.empty();

        std::uniform_int_distribution<std::size_t> pick(lo, hi - 1);
        std::vector<std::size_t> alive;
        if (sf2Lex) {
            alive.resize(hi - lo);
            std::iota(alive.begin(), alive.end(), lo);
        } else {
            alive.reserve(kPoolSize);
            for (std::size_t i = 0; i < kPoolSize; ++i) { alive.push_back(pick(sel)); }
        }

        auto pickRandom = [&]() -> PairIndividual const& {
            std::uniform_int_distribution<std::size_t> among(0, alive.size() - 1);
            return pop[alive[among(sel)]];
        };

        std::size_t R = 0;
        for (auto idx : alive) {
            if (!pop[idx].regionErrCache.empty()) { R = pop[idx].regionErrCache.size(); break; }
        }
        if (R == 0) { return pickRandom(); }

        std::vector<std::size_t> order(R);
        for (std::size_t k = 0; k < R; ++k) { order[k] = k; }
        std::shuffle(order.begin(), order.end(), sel);

        auto caseErr = [&](std::size_t idx, std::size_t k) -> Operon::Scalar {
            auto const& c = pop[idx].regionErrCache;
            return (k < c.size() && std::isfinite(c[k])) ? c[k] : INF;
        };

        for (auto k : order) {
            if (alive.size() <= 1) { break; }
            Operon::Scalar minVal = INF;
            for (auto idx : alive) { minVal = std::min(minVal, caseErr(idx, k)); }
            if (!std::isfinite(minVal)) { continue; }
            Operon::Scalar epsilon = sf2Lex ? (k < sf2LexEps_.size() ? sf2LexEps_[k] : Operon::Scalar{0})
                                            : Operon::Scalar{0.1} * minVal + Operon::Scalar{1e-9};
            std::vector<std::size_t> next;
            for (auto idx : alive) {
                if (caseErr(idx, k) <= minVal + epsilon) { next.push_back(idx); }
            }
            if (!next.empty()) { alive = std::move(next); }
            if (alive.size() == 1) { return pop[alive.front()]; }
        }
        return pickRandom();
    }

    static auto bestAccuracy(std::vector<PairIndividual>& pop) -> PairIndividual* {
        if (pop.empty()) { return nullptr; }
        return &*std::min_element(pop.begin(), pop.end(),
            [](auto const& a, auto const& b) { return a.shapeErr < b.shapeErr; });
    }
    static auto bestAccuracy(std::vector<PairIndividual> const& pop) -> PairIndividual const* {
        if (pop.empty()) { return nullptr; }
        return &*std::min_element(pop.begin(), pop.end(),
            [](auto const& a, auto const& b) { return a.shapeErr < b.shapeErr; });
    }
    auto floorReached(std::vector<PairIndividual>& pop) const -> bool {
        if (pop.empty()) { return false; }
        if (!config_.floorOnChamfer) {
            auto* best = bestAccuracy(pop);
            return best != nullptr && best->shapeErr <= config_.floorStop;
        }
        auto const& bc = *std::min_element(pop.begin(), pop.end(),
            [](auto const& a, auto const& b) { return a.chamferErr < b.chamferErr; });
        return bc.chamferErr <= config_.floorStop;
    }
    static auto countFront0(std::vector<PairIndividual> const& pop) -> std::size_t {
        return static_cast<std::size_t>(
            std::count_if(pop.begin(), pop.end(), [](auto const& p) { return p.rank == 0; }));
    }
