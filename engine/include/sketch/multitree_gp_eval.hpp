
    auto gmRandomBeats(Operon::RandomGenerator& rng, int kind, double pred, bool& coin) const -> bool
    {
        coin = false;
        if (config_.gradientMutationInRunErr && (!gmErrKnown(-1) || !gmErrKnown(kind))) {
            coin = true;
            std::uniform_real_distribution<double> u01(0.0, 1.0);
            return u01(rng) < 0.5;
        }
        return gmDraw(rng, -1, 0.0) < gmDraw(rng, kind, pred);
    }
    auto gmDraw(Operon::RandomGenerator& rng, int kind, double pred) const -> double
    {
        auto const& e = gmErrModel(kind);
        double const a = std::abs(pred);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        if (config_.gradientMutationInRunErr && !gmErrKnown(kind)) {
            if (e.pi0 > 0.0 && u01(rng) < e.pi0) { return 0.0; }
            return u01(rng) < 0.5 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
        }
        if (e.pi0 > 0.0 && u01(rng) < e.pi0) { return 0.0; }
        double const loc = pred + e.mu + e.b * a;
        double const sc = std::max(e.s0 + e.s1 * a, 0.0);
        std::student_t_distribution<double> td(std::max(e.nu, 1e-3));
        return loc + sc * td(rng);
    }
    auto gmLoc(int kind, double pred) const -> double
    {
        if (!config_.gradientMutationUseErrModel) { return pred; }
        auto const& e = gmErrModel(kind);
        return (1.0 - e.pi0) * (pred + e.mu + e.b * std::abs(pred));
    }
    auto gmLinRefitNml(PairIndividual const& ind, Operon::Dataset const& ds, Operon::Range range) const -> double
    {
        double const L = gmLoss(ind, ds, range);
        if (!std::isfinite(L)) { return std::numeric_limits<double>::quiet_NaN(); }
        GmLin lin;
        if (!gmBuildLin(ind, ds, range, lin)) { return std::numeric_limits<double>::quiet_NaN(); }
        return gmFitNats(L) + gmNatsAt(lin.rPerp - lin.rr, lin.rr)
            + Sketch::nml::structCodeCanonNats(static_cast<int>(ind.x.Length()))
            + Sketch::nml::structCodeCanonNats(static_cast<int>(ind.y.Length()))
            + static_cast<double>(lin.nq) * gmParamPerK();
    }

    void gradientMutate(SearchState const& s, Operon::RandomGenerator& rng, std::uint64_t seed,
                        Operon::Tree const& tree, Operon::Tree const& sibling, bool isX,
                        Operon::Tree const& parentTree, PairIndividual const& parent,
                        Operon::Tree& out, GmEvent& ev) const
    {
        using clk = std::chrono::steady_clock;
        auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        auto fallback = [&]() { ev.applied = 3; out = gmRandomMutate(s, rng, tree, parent); };
        if (s.activeDs == nullptr) { fallback(); return; }
        auto const& ds = *s.activeDs;
        auto const range = s.range;
        auto const t0 = clk::now();
        PairIndividual base;
        base.x = isX ? tree : sibling;
        base.y = isX ? sibling : tree;
        GmLin L;
        if (!gmBuildLin(base, ds, range, L)) { fallback(); return; }
        std::size_t const M = L.M;
        std::size_t const m = L.N + L.M;
        std::size_t const nq = L.nq;
        std::vector<double> G(M * nq), Qr(nq, 0.0), tmp;
        for (std::size_t q = 0; q < nq; ++q) {
            double const* Qq = L.Q.data() + q * m;
            gmApplyJUT(L, isX, Qq, tmp);
            std::copy(tmp.begin(), tmp.end(), G.begin() + static_cast<long>(q * M));
            double d = 0.0;
            for (std::size_t i = 0; i < m; ++i) { d += Qq[i] * L.r[i]; }
            Qr[q] = d;
        }
        auto const t1 = clk::now();
        ev.msLin = ms(t0, t1);
        GmTape T;
        if (!gmTape(tree, ds, range, T) || T.M != M) { ev.msTape = ms(t1, clk::now()); fallback(); return; }
        auto const t2 = clk::now();
        ev.msTape = ms(t1, t2);

        auto const& nodes = tree.Nodes();
        std::size_t const len = nodes.size();
        auto const& pset = s.prob->GetPrimitiveSet();
        auto enabled = [&](RegOp op) { return pset.IsEnabled(gmHash(op)); };
        static constexpr RegOp kUnary[] = {RegOp::Sin, RegOp::Cos, RegOp::Exp, RegOp::Log, RegOp::Sqrt, RegOp::Square};
        static constexpr RegOp kBinary[] = {RegOp::Add, RegOp::Sub, RegOp::Mul, RegOp::Div};
        static constexpr RegOp kInsert[] = {RegOp::Add, RegOp::Mul, RegOp::Div};
        std::vector<int> coefIn(len, 0);
        for (std::size_t k = 0; k < len; ++k) {
            std::size_t const st = k - nodes[k].Length;
            int c = 0;
            for (std::size_t q = st; q <= k; ++q) { c += nodes[q].Optimize ? 1 : 0; }
            coefIn[k] = c;
        }
        double const structBase = Sketch::nml::structCodeCanonNats(static_cast<int>(len));
        double const perK = gmParamPerK();
        double const refLevel = L.rPerp;
        double const refShift = L.rPerp - L.rr;

        std::vector<double> dvec(M), bvec(M), dr0, db, qv(nq), qb(nq), vres(m);
        std::vector<GmCand> cands;
        cands.reserve(len * 24);
        double const f32 = static_cast<double>(std::numeric_limits<float>::epsilon());
        auto score = [&](GmCand c, bool hasD, bool hasB) {
            double vv = 0.0;
            if (hasD) {
                gmApplyJU(L, isX, dvec.data(), dr0);
                for (std::size_t i = 0; i < m; ++i) { vres[i] = L.r[i] + dr0[i]; vv += vres[i] * vres[i]; }
                for (std::size_t q = 0; q < nq; ++q) {
                    double const* Gq = G.data() + q * M;
                    double d = 0.0;
                    for (std::size_t j = 0; j < M; ++j) { d += Gq[j] * dvec[j]; }
                    qv[q] = Qr[q] + d;
                }
            } else {
                for (std::size_t i = 0; i < m; ++i) { vres[i] = L.r[i]; }
                vv = L.rr;
                qv = Qr;
            }
            double qvv = 0.0;
            for (std::size_t q = 0; q < nq; ++q) { qvv += qv[q] * qv[q]; }
            double const ee = vv - qvv;
            double fixed = vv - L.rr, refit = ee - L.rr, beta = 0.0;
            bool newIdent = false;
            if (hasB) {
                gmApplyJU(L, isX, bvec.data(), db);
                double bb = 0.0, bv = 0.0;
                for (std::size_t i = 0; i < m; ++i) { bb += db[i] * db[i]; bv += db[i] * vres[i]; }
                if (bb > 0.0 && std::isfinite(bb)) {
                    beta = -bv / bb;
                    fixed = vv - bv * bv / bb - L.rr;
                    double qbb = 0.0, qbv = 0.0;
                    for (std::size_t q = 0; q < nq; ++q) {
                        double const* Gq = G.data() + q * M;
                        double d = 0.0;
                        for (std::size_t j = 0; j < M; ++j) { d += Gq[j] * bvec[j]; }
                        qb[q] = d; qbb += d * d; qbv += d * qv[q];
                    }
                    double const ff = bb - qbb, fe = bv - qbv;
                    if (c.kind == 0) { refit = L.rPerp - L.rr; }
                    else if (ff > f32 * bb) { refit = ee - fe * fe / ff - L.rr; newIdent = true; }
                }
            }
            if (!std::isfinite(fixed) || !std::isfinite(refit) || !std::isfinite(beta)) { return; }
            c.predFixed = fixed; c.predRefit = refit; c.beta = beta;
            c.nats = gmNatsAt(refit, L.rr);
            c.natsMut = gmNatsAt(refit - refShift, refLevel);
            c.dK += newIdent ? 1 : 0;
            c.dStruct = Sketch::nml::structCodeCanonNats(c.newLen) - structBase;
            c.dParam = static_cast<double>(c.dK) * perK;
            c.pred = c.natsMut + c.dStruct + c.dParam;
            if (!std::isfinite(c.nats) || !std::isfinite(c.pred)) { return; }
            c.loc = gmLoc(c.kind, c.pred);
            cands.push_back(c);
        };
        auto allFinite = [&](std::vector<double> const& v) {
            for (auto z : v) { if (!std::isfinite(z)) { return false; } }
            return true;
        };
        int const L0 = static_cast<int>(len);
        for (std::size_t k = 0; k < len; ++k) {
            if (T.ok[k] == 0) { continue; }
            double const* u = T.val.data() + k * M;
            double const* a = T.adj.data() + k * M;
            auto const& nd = nodes[k];
            int const span = static_cast<int>(nd.Length) + 1;
            bool const isConst = nd.IsConstant();
            bool const isVar = nd.IsVariable();
            GmCand c; c.node = static_cast<int>(k);
            if (isConst || isVar) {
                c.newLen = L0;
                if (nd.Optimize) {
                    for (std::size_t j = 0; j < M; ++j) { bvec[j] = a[j] * (isConst ? 1.0 : T.t[j]); }
                    c.kind = 0; c.dK = 0; score(c, false, true);
                }
                for (std::size_t j = 0; j < M; ++j) { dvec[j] = -a[j] * u[j]; bvec[j] = a[j] * (isConst ? T.t[j] : 1.0); }
                c.kind = 2; c.leaf = isConst ? 2 : 1; c.dK = -(nd.Optimize ? 1 : 0); score(c, true, true); c.leaf = 0;
            } else {
                auto const& in = T.prog.insts[k];
                double const w = static_cast<double>(in.weight);
                double const* x1 = T.val.data() + static_cast<std::size_t>(in.src1) * M;
                double const* x2 = T.val.data() + static_cast<std::size_t>(in.src2) * M;
                c.newLen = L0; c.dK = 0;
                if (in.arity == 1) {
                    for (auto g : kUnary) {
                        if (g == in.op || !enabled(g)) { continue; }
                        for (std::size_t j = 0; j < M; ++j) { dvec[j] = a[j] * (gmUnary(g, x1[j]) * w - u[j]); }
                        if (!allFinite(dvec)) { continue; }
                        c.kind = 1; c.op = g; score(c, true, false);
                    }
                } else {
                    for (auto g : kBinary) {
                        if (g == in.op || !enabled(g)) { continue; }
                        for (std::size_t j = 0; j < M; ++j) { dvec[j] = a[j] * (gmBinary(g, x1[j], x2[j]) * w - u[j]); }
                        if (!allFinite(dvec)) { continue; }
                        c.kind = 1; c.op = g; score(c, true, false);
                    }
                }
                c.op = RegOp::Identity;
                for (int q = 0; q < static_cast<int>(in.arity); ++q) {
                    std::size_t const ch = (q == 0) ? in.src1 : in.src2;
                    if (T.ok[ch] == 0) { continue; }
                    double const* uc = T.val.data() + ch * M;
                    for (std::size_t j = 0; j < M; ++j) { dvec[j] = a[j] * (uc[j] - u[j]); }
                    c.kind = 3; c.child = static_cast<int>(ch);
                    c.newLen = L0 - span + static_cast<int>(nodes[ch].Length) + 1;
                    c.dK = -(coefIn[k] - coefIn[ch]);
                    score(c, true, false);
                }
                c.child = -1;
                for (int leaf = 1; leaf <= 2; ++leaf) {
                    for (std::size_t j = 0; j < M; ++j) { dvec[j] = -a[j] * u[j]; bvec[j] = a[j] * (leaf == 1 ? 1.0 : T.t[j]); }
                    c.kind = 4; c.leaf = leaf; c.newLen = L0 - span + 1; c.dK = -coefIn[k];
                    score(c, true, true);
                }
                c.leaf = 0;
            }
            if (isConst) { continue; }
            if (len + 1 <= config_.maxLength) {
                c.newLen = L0 + 1; c.dK = 0;
                for (auto g : kUnary) {
                    if (!enabled(g)) { continue; }
                    for (std::size_t j = 0; j < M; ++j) { dvec[j] = a[j] * (gmUnary(g, u[j]) - u[j]); }
                    if (!allFinite(dvec)) { continue; }
                    c.kind = 5; c.op = g; score(c, true, false);
                }
            }
            if (len + 2 <= config_.maxLength) {
                c.newLen = L0 + 2;
                for (auto g : kInsert) {
                    if (!enabled(g)) { continue; }
                    for (int leaf = 1; leaf <= 2; ++leaf) {
                        for (std::size_t j = 0; j < M; ++j) {
                            double const l = (leaf == 1) ? 1.0 : T.t[j];
                            if (g == RegOp::Add) { dvec[j] = 0.0; bvec[j] = a[j] * l; }
                            else if (g == RegOp::Mul) { dvec[j] = -a[j] * u[j]; bvec[j] = a[j] * l * u[j]; }
                            else { dvec[j] = -a[j] * u[j]; bvec[j] = a[j] * l / u[j]; }
                        }
                        if (!allFinite(bvec)) { continue; }
                        c.kind = 6; c.op = g; c.leaf = leaf; c.dK = 0;
                        score(c, g != RegOp::Add, true);
                    }
                }
                c.leaf = 0; c.op = RegOp::Identity;
            }
        }
        auto const t3 = clk::now();
        ev.msCand = ms(t2, t3);
        ev.cands = cands.size();
        if (cands.empty()) { fallback(); ev.msSel = ms(t3, clk::now()); return; }

        bool const em = config_.gradientMutationUseErrModel;
        std::vector<char> alive(cands.size(), 1);
        std::size_t chosen = cands.size();
        bool pickRandom = false;
        for (std::size_t tries = 0; tries < cands.size(); ++tries) {
            std::size_t pick = cands.size();
            double best = std::numeric_limits<double>::infinity();
            if ((config_.gradientMutationMode == 2 || config_.gradientMutationMode == 3) && em
                && config_.gradientMutationErrcalPImprove) {
                double bestP = -1.0;
                for (std::size_t i = 0; i < cands.size(); ++i) {
                    if (alive[i] == 0) { continue; }
                    double const pi = gmPImprove(cands[i].kind, cands[i].pred);
                    if (pi > bestP || (pi == bestP && cands[i].loc < best)) { bestP = pi; best = cands[i].loc; pick = i; }
                }
                bool coin = false;
                if (tries == 0 && config_.gradientMutationMode == 2 && pick < cands.size()
                    && gmRandomBeats(rng, cands[pick].kind, cands[pick].pred, coin)) {
                    ev.coin = coin; pickRandom = true; break;
                }
                if (tries == 0) { ev.coin = coin; }
            } else if ((config_.gradientMutationMode == 2 || config_.gradientMutationMode == 3) && em) {
                bool const ir = config_.gradientMutationInRunErr;
                for (std::size_t i = 0; i < cands.size(); ++i) {
                    if (alive[i] == 0) { continue; }
                    double const y = gmDraw(rng, cands[i].kind, cands[i].pred);
                    if (y < best || (ir && y == best && pick < cands.size() && cands[i].loc < cands[pick].loc)) { best = y; pick = i; }
                }
                if (tries == 0 && config_.gradientMutationMode == 2 && gmDraw(rng, -1, 0.0) < best) { pickRandom = true; break; }
            } else {
                for (std::size_t i = 0; i < cands.size(); ++i) {
                    if (alive[i] != 0 && cands[i].loc < best) { best = cands[i].loc; pick = i; }
                }
                if (config_.gradientMutationInRunErr) {
                    bool coin = false;
                    if (tries == 0 && em && config_.gradientMutationMode == 1 && pick < cands.size()
                        && gmRandomBeats(rng, cands[pick].kind, cands[pick].pred, coin)) {
                        ev.coin = coin; pickRandom = true; break;
                    }
                    if (tries == 0) { ev.coin = coin; }
                } else if (tries == 0 && em && config_.gradientMutationMode == 1 && pick < cands.size()
                    && gmDraw(rng, -1, 0.0) < gmDraw(rng, cands[pick].kind, cands[pick].pred)) {
                    pickRandom = true; break;
                }
            }
            if (pick >= cands.size()) { break; }
            auto t = gmApplyCand(s, tree, cands[pick]);
            if (t.Length() <= config_.maxLength && t.Depth() <= config_.maxDepth && !t.Nodes().empty()) {
                out = std::move(t);
                chosen = pick;
                break;
            }
            alive[pick] = 0;
        }
        auto const t4 = clk::now();
        ev.msSel = ms(t3, t4);
        if (pickRandom) { ev.applied = 2; out = gmRandomMutate(s, rng, tree, parent); }
        else if (chosen >= cands.size()) { fallback(); }
        else { ev.applied = 1; ev.kind = cands[chosen].kind; ev.predChosen = cands[chosen].pred; }

        if (config_.gradientMutationLog) {
            int const evId = gmEventId_.fetch_add(1);
            bool clone = parentTree.Length() == tree.Length();
            if (clone) {
                auto const& pn = parentTree.Nodes();
                for (std::size_t q = 0; q < len && clone; ++q) {
                    clone = pn[q].HashValue == nodes[q].HashValue && pn[q].Arity == nodes[q].Arity && pn[q].Value == nodes[q].Value;
                }
            }
            double const L0loss = gmLoss(base, ds, range);
            PairIndividual ref = base;
            optConstsGuarded(ref, ds, range);
            double const Lref0 = gmLoss(ref, ds, range);
            double const nmlRef = gmNmlOf(ref, ds, range);
            Operon::RandomGenerator rngShadow(seed ^ 0x5DEECE66DULL);
            double randomNml = std::numeric_limits<double>::quiet_NaN();
            int randomLen = 0;
            {
                auto rt = gmRandomMutate(s, rngShadow, tree, parent);
                randomLen = static_cast<int>(rt.Length());
                PairIndividual pr = base;
                (isX ? pr.x : pr.y) = rt;
                optConstsGuarded(pr, ds, range);
                randomNml = gmNmlOf(pr, ds, range) - nmlRef;
            }
            std::vector<std::size_t> order(cands.size());
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q) { return cands[p].pred < cands[q].pred; });
            std::vector<std::size_t> idx;
            std::size_t const K = std::min<std::size_t>(config_.gradientMutationLogTopK, order.size());
            for (std::size_t q = 0; q < K; ++q) { idx.push_back(order[q]); }
            if (chosen < cands.size() && std::find(idx.begin(), idx.end(), chosen) == idx.end()) { idx.push_back(chosen); }
            std::vector<GradMutRecord> recs;
            for (auto ci : idx) {
                auto const& c = cands[ci];
                auto t = gmApplyCand(s, tree, c);
                if (!(t.Length() <= config_.maxLength && t.Depth() <= config_.maxDepth)) { continue; }
                GradMutRecord rec;
                rec.gen = static_cast<int>(s.currentGen);
                rec.event = evId;
                rec.rankPred = static_cast<int>(std::find(order.begin(), order.end(), ci) - order.begin());
                rec.kind = c.kind; rec.node = c.node; rec.treeLen = static_cast<int>(len);
                rec.nCands = static_cast<int>(cands.size());
                rec.isX = isX; rec.chosen = (ci == chosen && ev.applied == 1); rec.childIsParentClone = clone;
                rec.randomChosen = (ev.applied == 2);
                rec.predFixed = c.predFixed; rec.predRefitTotal = c.predRefit;
                rec.predRefitMut = c.predRefit - refShift;
                rec.predRefitMutNats = c.natsMut;
                rec.predNml = c.pred; rec.dStruct = c.dStruct; rec.dParam = c.dParam; rec.dK = c.dK;
                rec.L0 = L0loss; rec.Lref0 = Lref0;
                PairIndividual pm = base;
                (isX ? pm.x : pm.y) = t;
                double const Lfix = gmLoss(pm, ds, range);
                PairIndividual pr = pm;
                optConstsGuarded(pr, ds, range);
                double const Lr = gmLoss(pr, ds, range);
                rec.realFixed = Lfix - L0loss;
                rec.realRefitTotal = Lr - L0loss;
                rec.realRefitMut = Lr - Lref0;
                rec.realRefitMutNats = gmNatsAt(rec.realRefitMut, Lref0);
                rec.realFixedNats = gmNatsAt(rec.realFixed, L0loss);
                rec.predFixedNats = gmNatsAt(c.predFixed, L.rr);
                rec.realNmlMut = gmNmlOf(pr, ds, range) - nmlRef;
                rec.randomRealNmlMut = randomNml;
                rec.randomLen = randomLen;
                recs.push_back(rec);
            }
            std::lock_guard<std::mutex> lk(gmLogMu_);
            gmLog_.insert(gmLog_.end(), recs.begin(), recs.end());
            ev.msLog = ms(t4, clk::now());
        }
    }

    void gmEffRecord(SearchState const& s, Operon::Tree const& pre, Operon::Tree const& post,
                     Operon::Tree const& sibling, bool isX, PairIndividual const& parent, int type, int kind,
                     double pred, int arm = -1) const
    {
        if (s.activeDs == nullptr) { return; }
        auto const& ds = *s.activeDs;
        auto const range = s.range;
        auto refitNml = [&](PairIndividual ind) {
            optConstsGuarded(ind, ds, range);
            return gmNmlOf(ind, ds, range);
        };
        PairIndividual a, b;
        a.x = isX ? pre : sibling;  a.y = isX ? sibling : pre;
        b.x = isX ? post : sibling; b.y = isX ? sibling : post;
        PairIndividual p;
        p.x = parent.x; p.y = parent.y;
        GradMutEff e;
        e.gen = static_cast<int>(s.currentGen);
        e.type = type; e.kind = kind; e.pred = pred; e.arm = arm;
        double const nPost = refitNml(b);
        e.dVsPre = nPost - refitNml(a);
        e.dVsParent = nPost - refitNml(p);
        e.lenDelta = static_cast<int>(post.Length()) - static_cast<int>(pre.Length());
        std::lock_guard<std::mutex> lk(gmLogMu_);
        gmEff_.push_back(e);
    }

    struct GmPendingEvent { bool isX{true}; std::uint64_t seed{0}; int arm{-1}; };
    struct GmPendingChild { std::size_t idx{0}; PairIndividual const* parent{nullptr}; std::vector<GmPendingEvent> ev; };
    struct GmBanditRec { std::size_t idx{0}; int arm{0}; double ms{0}; PairIndividual pre; };
    struct GmIrStub { std::size_t idx{0}; int model{0}; double pred{0}; PairIndividual pre; bool noop{false}; std::uint64_t key{0}; };
    auto gmNatsResolution(double L) const -> double
    {
        double const Lf = static_cast<double>(static_cast<float>(L));
        double const up = static_cast<double>(std::nextafter(static_cast<float>(L), std::numeric_limits<float>::infinity()));
        double const d = std::abs(gmFitNats(up) - gmFitNats(Lf));
        return (d > 0.0 && std::isfinite(d)) ? d : std::numeric_limits<double>::min();
    }
    static auto gmTreeHash(Operon::Tree const& t, std::uint64_t h) -> std::uint64_t
    {
        auto mix = [&h](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
        mix(t.Nodes().size());
        for (auto const& nd : t.Nodes()) {
            mix(static_cast<std::uint64_t>(nd.HashValue));
            mix(static_cast<std::uint64_t>(nd.Arity));
            std::uint32_t vb = 0;
            static_assert(sizeof(nd.Value) == sizeof(vb));
            std::memcpy(&vb, &nd.Value, sizeof(vb));
            mix(vb);
        }
        return h;
    }

    static auto parallelBreedSeed(std::uint64_t runSeed, std::size_t gen, std::size_t idx, int tag) -> std::uint64_t
    {
        auto splitmix64 = [](std::uint64_t x) -> std::uint64_t {
            x += 0x9E3779B97F4A7C15ULL;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
            return x ^ (x >> 31);
        };
        std::uint64_t h = splitmix64(runSeed);
        h = splitmix64(h ^ static_cast<std::uint64_t>(gen));
        h = splitmix64(h ^ static_cast<std::uint64_t>(idx));
        h = splitmix64(h ^ static_cast<std::uint64_t>(tag));
        return h;
    }

    static constexpr int kGmArms = 7;
    auto gmBanditPick(Operon::RandomGenerator& rng) const -> int
    {
        if (config_.gradientMutationMode != 3) { return -1; }
        int nLow = 0;
        for (int a = 0; a < kGmArms; ++a) { if (gmArmN_[a] < 2.0) { ++nLow; } }
        int arm = 0;
        if (nLow > 0) {
            std::uniform_int_distribution<int> ud(0, nLow - 1);
            int k = ud(rng);
            for (int a = 0; a < kGmArms; ++a) {
                if (gmArmN_[a] < 2.0) { if (k == 0) { arm = a; break; } --k; }
            }
        } else {
            double best = -std::numeric_limits<double>::infinity();
            for (int a = 0; a < kGmArms; ++a) {
                double const n = gmArmN_[a];
                double const mean = gmArmSum_[a] / n;
                double const var = std::max((gmArmSumSq_[a] - n * mean * mean) / (n - 1.0), 0.0);
                std::student_t_distribution<double> td(n - 1.0);
                double const draw = mean + std::sqrt(var / n) * td(rng);
                if (draw > best) { best = draw; arm = a; }
            }
        }
        ++gmStats_.armPulls[arm];
        return arm;
    }
    auto gmArmMutate(SearchState const& s, Operon::RandomGenerator& rng, int arm, Operon::Tree t,
                     PairIndividual const& parent) const -> Operon::Tree
    {
        switch (arm) {
            case 0: return s.mut->onePoint(rng, std::move(t));
            case 1: return (*s.mut->changeFunc)(rng, std::move(t));
            case 2: return (*s.mut->replaceSub)(rng, std::move(t));
            case 3: return (*s.mut->removeSub)(rng, std::move(t));
            case 4: return (*s.mut->insertOp)(rng, std::move(t));
            default: {
                auto* node = parent.perPointResid_.empty() ? nullptr : noiseGuidedLeafPick(rng, t, parent.perPointResid_);
                if (node == nullptr) { return s.mut->onePoint(rng, std::move(t)); }
                std::normal_distribution<Operon::Scalar> dist(Operon::Scalar{0}, Operon::Scalar{3});
                node->Value += dist(rng);
                return t;
            }
        }
    }
    void gmBanditUpdate(SearchState const& s, std::vector<PairIndividual> const& offspring,
                        std::vector<GmBanditRec> const& recs, double scoreWallMs, std::size_t workers,
                        std::vector<char> const* refitMask = nullptr) const
    {
        bool const refitRw = config_.gradientMutationBanditRefitReward;
        auto const t0 = std::chrono::steady_clock::now();
        if (s.activeDs == nullptr || offspring.empty()) { return; }
        double const perChildMs = scoreWallMs * static_cast<double>(std::max<std::size_t>(workers, 1))
                                / static_cast<double>(offspring.size());
        std::vector<double> r(recs.size(), std::numeric_limits<double>::quiet_NaN());
        parallelApply(recs.size(), workers, [&](std::size_t i) {
            auto const& rc = recs[i];
            if (refitRw && (refitMask == nullptr || rc.idx >= refitMask->size() || (*refitMask)[rc.idx] == 0)) { return; }
            double const post = gmNmlOf(offspring[rc.idx], *s.activeDs, s.range);
            double const pre = refitRw ? gmLinRefitNml(rc.pre, *s.activeDs, s.range)
                                       : gmNmlOf(rc.pre, *s.activeDs, s.range);
            double const cost = rc.ms + perChildMs;
            if (std::isfinite(post) && std::isfinite(pre) && cost > 0.0) { r[i] = -(post - pre) / cost; }
        });
        for (std::size_t i = 0; i < recs.size(); ++i) {
            if (!std::isfinite(r[i])) { continue; }
            int const a = recs[i].arm;
            gmArmN_[a] += 1.0; gmArmSum_[a] += r[i]; gmArmSumSq_[a] += r[i] * r[i];
        }
        for (int a = 0; a < kGmArms; ++a) {
            gmStats_.armObs[a] = static_cast<std::size_t>(gmArmN_[a]);
            gmStats_.armMean[a] = gmArmN_[a] > 0 ? gmArmSum_[a] / gmArmN_[a] : 0.0;
        }
        gmStats_.banditMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }

    void gmIrObserve(SearchState const& s, std::vector<PairIndividual> const& offspring,
                     std::vector<GmIrStub> const& stubs, std::vector<char> const& refitMask, std::size_t workers) const
    {
        auto const t0 = std::chrono::steady_clock::now();
        if (s.activeDs == nullptr) { return; }
        gmStats_.irStubs += stubs.size();
        std::vector<double> y(stubs.size(), std::numeric_limits<double>::quiet_NaN()), dl(stubs.size(), 0.0);
        parallelApply(stubs.size(), workers, [&](std::size_t i) {
            auto const& st = stubs[i];
            if (st.idx >= refitMask.size() || refitMask[st.idx] == 0) { return; }
            if (st.noop) { y[i] = 0.0; return; }
            double const post = gmNmlOf(offspring[st.idx], *s.activeDs, s.range);
            double const ref = gmLinRefitNml(st.pre, *s.activeDs, s.range);
            if (std::isfinite(post) && std::isfinite(ref)) {
                y[i] = post - ref;
                dl[i] = gmNatsResolution(gmLoss(offspring[st.idx], *s.activeDs, s.range))
                      + gmNatsResolution(gmLoss(st.pre, *s.activeDs, s.range));
            }
        });
        std::array<char, 8> grew{};
        for (std::size_t i = 0; i < stubs.size(); ++i) {
            if (!std::isfinite(y[i]) || (stubs[i].model < 7 && !std::isfinite(stubs[i].pred))) { continue; }
            if (!gmIrSeen_.insert(stubs[i].key).second) { ++gmStats_.irDup; continue; }
            auto const m = static_cast<std::size_t>(stubs[i].model);
            ++gmIrN_[m];
            ++gmStats_.irObs;
            if (stubs[i].noop || y[i] == 0.0) { ++gmIrN0_[m]; continue; }
            gmIrP_[m].push_back(stubs[i].pred);
            gmIrY_[m].push_back(y[i]);
            gmIrD_[m].push_back(dl[i]);
            static bool const g6dDebug = std::getenv("SKETCH_G6D_DEBUG") != nullptr;
            if (g6dDebug) {
                std::fprintf(stderr, "G6DOBS %zu %.9g %.9g %.3g lx=%zu ly=%zu %s | %s\n", m, stubs[i].pred, y[i], dl[i],
                             offspring[stubs[i].idx].x.Length(), offspring[stubs[i].idx].y.Length(),
                             Operon::InfixFormatter::Format(offspring[stubs[i].idx].x, *s.activeDs, 3).c_str(),
                             Operon::InfixFormatter::Format(offspring[stubs[i].idx].y, *s.activeDs, 3).c_str());
            }
            grew[m] = 1;
        }
        std::vector<std::size_t> todo;
        for (std::size_t m = 0; m < 8; ++m) {
            gmIrModel_[m].pi0 = (static_cast<double>(gmIrN0_[m]) + 0.5) / (static_cast<double>(gmIrN_[m]) + 1.0);
            std::size_t const dim = (m == 7) ? 3 : 5;
            if (grew[m] != 0 && gmIrY_[m].size() > dim) { todo.push_back(m); }
        }
        auto const f0 = std::chrono::steady_clock::now();
        std::vector<char> ok(todo.size(), 0);
        std::vector<GradMutErrModel> nm(todo.size());
        parallelApply(todo.size(), workers, [&](std::size_t j) {
            auto const m = todo[j];
            nm[j] = gmIrModel_[m];
            ok[j] = gmstat::fitErrMap(gmIrP_[m], gmIrY_[m], gmIrD_[m], m != 7, nm[j], gmIrFit_[m] != 0) ? 1 : 0;
        });
        for (std::size_t j = 0; j < todo.size(); ++j) {
            auto const m = todo[j];
            double const pi0 = gmIrModel_[m].pi0;
            gmIrModel_[m] = (ok[j] != 0) ? nm[j] : GradMutErrModel{};
            gmIrModel_[m].pi0 = pi0;
            gmIrFit_[m] = ok[j] != 0 ? 1 : 0;
        }
        auto const t1 = std::chrono::steady_clock::now();
        for (std::size_t m = 0; m < 8; ++m) {
            gmStats_.irN[m] = gmIrN_[m]; gmStats_.irN0[m] = gmIrN0_[m];
            gmStats_.irFit[m] = gmIrFit_[m] != 0; gmStats_.irModel[m] = gmIrModel_[m];
        }
        gmStats_.irFitMs += std::chrono::duration<double, std::milli>(t1 - f0).count();
        gmStats_.irMs += std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (gmG6eActive()) { gmIrPoolG6e(); }
    }

    void gmIrPoolG6e() const
    {
        std::vector<double> allP, allY;
        for (std::size_t m = 0; m < 7; ++m) {
            allP.insert(allP.end(), gmIrP_[m].begin(), gmIrP_[m].end());
            allY.insert(allY.end(), gmIrY_[m].begin(), gmIrY_[m].end());
        }
        double poolMed = 0.0, poolMad = 1.0;
        if (!allY.empty()) { gmstat::medmad(allP, allY, true, poolMed, poolMad); }
        poolMad = std::max(poolMad, 1e-9);

        std::array<double, 8> nAll{}, x0{};
        int kindsWithData = 0;
        for (std::size_t m = 0; m < 8; ++m) {
            nAll[m] = static_cast<double>(gmIrN_[m]);
            x0[m] = static_cast<double>(gmIrN0_[m]);
            if (gmIrN_[m] > 0) { ++kindsWithData; }
        }
        auto const [a, b] = gmstat::betaBinomialTypeII(
            std::vector<double>(nAll.begin(), nAll.end()), std::vector<double>(x0.begin(), x0.end()));

        std::array<double, 8> muHat{}, lnsHat{}, muV{}, lnsV{};
        std::vector<char> has(8, 0);
        for (std::size_t m = 0; m < 8; ++m) {
            std::size_t const nk = gmIrY_[m].size();
            if (nk == 0) { continue; }
            has[m] = 1;
            double med = 0.0, mad = 0.0;
            gmstat::medmad(gmIrP_[m], gmIrY_[m], m != 7, med, mad);
            double const sk = (nk >= 2) ? std::max(mad, 1e-9) : poolMad;
            muHat[m] = med;
            muV[m] = std::max(sk * sk / static_cast<double>(nk), 1e-15);
            lnsHat[m] = std::log(sk);
            lnsV[m] = (nk >= 2) ? 1.0 / (2.0 * (static_cast<double>(nk) - 1.0)) : std::numeric_limits<double>::infinity();
        }
        for (std::size_t m = 0; m < 8; ++m) { if (!std::isfinite(lnsV[m])) { lnsV[m] = 1e18; } }
        auto const [muOut, muPool] = gmstat::dlPool(
            std::vector<double>(muHat.begin(), muHat.end()), std::vector<double>(muV.begin(), muV.end()), has);
        double const lnsPoolDefault = std::log(poolMad);
        auto [lnsOut, lnsPool] = gmstat::dlPool(
            std::vector<double>(lnsHat.begin(), lnsHat.end()), std::vector<double>(lnsV.begin(), lnsV.end()), has);
        if (kindsWithData == 0) { lnsPool = lnsPoolDefault; }

        double poolB = 0.0, poolS1 = 0.0, poolNu = 1.0;
        {
            std::vector<double> pp, yy, dd;
            for (std::size_t m = 0; m < 7; ++m) {
                double const mu = muOut.empty() ? 0.0 : muOut[m];
                for (std::size_t i = 0; i < gmIrY_[m].size(); ++i) {
                    pp.push_back(gmIrP_[m][i]);
                    yy.push_back(gmIrY_[m][i] - mu);
                    dd.push_back(gmIrD_[m][i]);
                }
            }
            GradMutErrModel gm{};
            if (pp.size() > 5 && gmstat::fitErrMap(pp, yy, dd, true, gm, false)) {
                poolB = gm.b; poolS1 = gm.s1; poolNu = gm.nu;
            }
        }
        double nu7 = (gmIrFit_[7] != 0) ? gmIrModel_[7].nu : ((poolNu > 0.0) ? poolNu : 1.0);

        bool const anyData = kindsWithData > 0;
        for (std::size_t m = 0; m < 8; ++m) {
            GradMutErrModel gm{};
            gm.pi0 = (x0[m] + a) / (nAll[m] + a + b);
            gm.mu = muOut.empty() ? 0.0 : muOut[m];
            gm.s0 = std::exp(lnsOut.empty() ? lnsPool : lnsOut[m]);
            gm.b = (m < 7) ? poolB : 0.0;
            gm.s1 = (m < 7) ? poolS1 : 0.0;
            gm.nu = (m < 7) ? poolNu : nu7;
            gmIrModel_[m] = gm;
            gmIrFit_[m] = anyData ? 1 : 0;
            gmStats_.irFit[m] = gmIrFit_[m] != 0;
            gmStats_.irModel[m] = gmIrModel_[m];
        }
    }

    void gmDeferredPass(SearchState const& s, std::vector<PairIndividual>& offspring,
                        std::vector<GmPendingChild> const& pend, std::vector<GmBanditRec>* bandit,
                        std::vector<GmIrStub>* irObs = nullptr) const
    {
        auto const w0 = std::chrono::steady_clock::now();
        bool const shadow = config_.gradientMutationLog || config_.gradientMutationEffLog;
        auto const d0 = nConstDescents_.load(), l0 = nConstLossCalls_.load(), i0 = nConstIters_.load(), ee0 = nEvalErrors_.load();
        std::vector<std::vector<GmEvent>> res(pend.size());
        std::vector<double> evMs(pend.size(), 0.0);
        std::vector<PairIndividual> preChild(bandit != nullptr ? pend.size() : 0);
        parallelApply(pend.size(), workerCount(), [&](std::size_t i) {
            auto const& pc = pend[i];
            auto& child = offspring[pc.idx];
            if (bandit != nullptr && pc.ev.size() == 1) { preChild[i].x = child.x; preChild[i].y = child.y; }
            for (auto const& pe : pc.ev) {
                auto const e0 = std::chrono::steady_clock::now();
                Operon::RandomGenerator rng(pe.seed);
                auto& tr = pe.isX ? child.x : child.y;
                auto const& sib = pe.isX ? child.y : child.x;
                Operon::Tree const pre = tr;
                Operon::Tree out;
                GmEvent ev;
                if (pe.arm >= 0 && pe.arm != kGmArms - 1) {
                    out = gmArmMutate(s, rng, pe.arm, pre, *pc.parent);
                    ev.applied = 4;
                } else {
                    gradientMutate(s, rng, pe.seed, pre, sib, pe.isX, pe.isX ? pc.parent->x : pc.parent->y, *pc.parent, out, ev);
                }
                evMs[i] += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - e0).count() - ev.msLog;
                if (config_.gradientMutationEffLog) {
                    GmEvent evc = ev;
                    gmEffRecord(s, pre, out, sib, pe.isX, *pc.parent, evc.applied, evc.kind, evc.predChosen, pe.arm);
                }
                tr = std::move(out);
                res[i].push_back(ev);
            }
        });
        if (shadow) { nConstDescents_ = d0; nConstLossCalls_ = l0; nConstIters_ = i0; nEvalErrors_ = ee0; }
        for (std::size_t i = 0; i < pend.size(); ++i) {
            auto const& pc = pend[i];
            auto& child = offspring[pc.idx];
            bool grad = false;
            for (auto const& ev : res[i]) {
                ++gmStats_.events;
                if (ev.applied == 1) { ++gmStats_.applied; ++gmStats_.kindApplied[ev.kind]; grad = true; }
                else if (ev.applied == 2) { ++gmStats_.randomModel; }
                else if (ev.applied == 4) { ++gmStats_.banditRandomArm; }
                else { ++gmStats_.fallback; }
                if (ev.coin) { ++gmStats_.irCoin; }
                gmStats_.candidates += ev.cands;
                gmStats_.msLin += ev.msLin; gmStats_.msTape += ev.msTape; gmStats_.msCand += ev.msCand;
                gmStats_.msSel += ev.msSel; gmStats_.msLog += ev.msLog;
            }
            bool const tooBig = child.x.Length() > 2 * config_.maxLength || child.y.Length() > 2 * config_.maxLength;
            bool const invalid = !tooBig && !structurallyValid(child);
            if (tooBig || invalid) {
                child = *pc.parent;
                if (invalid) { ++nInvalidChildren_; }
                if (!config_.bracketNewtonConfig.warmStart) {
                    child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                    child.warmCompletedSweeps = 0;
                }
                continue;
            }
            if (grad) {
                child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                child.warmCompletedSweeps = 0;
            }
            if (irObs != nullptr && pc.ev.size() == 1 && (res[i][0].applied == 1 || res[i][0].applied == 2)) {
                int const model = res[i][0].applied == 1 ? res[i][0].kind : 7;
                auto const& pre = preChild[i];
                std::uint64_t const hPre = gmTreeHash(pre.y, gmTreeHash(pre.x, 1469598103934665603ULL));
                std::uint64_t const hPost = gmTreeHash(child.y, gmTreeHash(child.x, 1469598103934665603ULL));
                std::uint64_t key = hPre * 1099511628211ULL ^ hPost;
                key = key * 1099511628211ULL ^ static_cast<std::uint64_t>(model);
                irObs->push_back({pc.idx, model, res[i][0].applied == 1 ? res[i][0].predChosen : 0.0, pre, hPre == hPost, key});
            }
            if (bandit != nullptr && pc.ev.size() == 1 && pc.ev[0].arm >= 0) {
                bandit->push_back({pc.idx, pc.ev[0].arm, evMs[i], std::move(preChild[i])});
            }
        }
        gmStats_.passWallMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - w0).count();
    }

    void gmTrajectoryPoint(std::vector<PairIndividual> const& pop, SearchState const& s) const
    {
        auto const h0 = std::chrono::steady_clock::now();
        if (!gmTrajStarted_) { gmTrajStarted_ = true; gmTrajT0_ = h0; gmTrajHookMs_ = 0.0; }
        std::vector<double> v(pop.size(), std::numeric_limits<double>::infinity());
        if (s.activeDs != nullptr) {
            parallelApply(pop.size(), workerCount(), [&](std::size_t i) { v[i] = gmNmlOf(pop[i], *s.activeDs, s.range); });
        }
        double best = std::numeric_limits<double>::infinity();
        for (auto z : v) { if (std::isfinite(z)) { best = std::min(best, z); } }
        double const wall = std::chrono::duration<double, std::milli>(h0 - gmTrajT0_).count() - gmTrajHookMs_;
        gmStats_.traj.push_back({static_cast<int>(s.currentGen), wall, best, static_cast<int>(dataX_.size())});
        gmTrajHookMs_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - h0).count();
    }

    auto mutateChildTree(SearchState& s, Operon::Tree tree, PairIndividual const& parent,
                          bool isX, bool& structuralChange, Operon::Tree const* childSibling = nullptr) const -> Operon::Tree
    {
        structuralChange = false;
        if (config_.gradientMutationMode != 0) {
            std::uint64_t const seed = s.rng();
            Operon::RandomGenerator rng(seed);
            Operon::Tree const& sib = (childSibling != nullptr) ? *childSibling : (isX ? parent.y : parent.x);
            Operon::Tree out;
            GmEvent ev;
            gradientMutate(s, rng, seed, tree, sib, isX, isX ? parent.x : parent.y, parent, out, ev);
            ++gmStats_.events;
            if (ev.applied == 1) { ++gmStats_.applied; ++gmStats_.kindApplied[ev.kind]; }
            else if (ev.applied == 2) { ++gmStats_.randomModel; }
            else { ++gmStats_.fallback; }
            if (ev.coin) { ++gmStats_.irCoin; }
            gmStats_.candidates += ev.cands;
            structuralChange = ev.applied == 1;
            return out;
        }
        auto normalMutate = [&](Operon::Tree t) -> Operon::Tree {
            if (config_.useJacobianNodeSensitivity) {
                std::uniform_real_distribution<double> coinJ(0.0, 1.0);
                if (coinJ(s.rng) < 0.5) {
                    return mutRandomOp(s, s.rng, t, parent);
                }
                Operon::Tree const& sibling = isX ? parent.y : parent.x;
                auto* node = (s.activeDs != nullptr)
                    ? jacobianLeafPick(s, t, sibling, isX, *s.activeDs, s.range)
                    : nullptr;
                if (node == nullptr) { return mutRandomOp(s, s.rng, t, parent); }
                std::normal_distribution<Operon::Scalar> distJ(Operon::Scalar{0}, Operon::Scalar{3});
                auto const beforeJ = node->Value;
                node->Value += distJ(s.rng);
                if (opDiagEnabled()) { opDiagRecord(OpMutNoiseGuided, node->Value == beforeJ, false); }
                return t;
            }
            if (!config_.useNoiseGuidedMutation || parent.perPointResid_.empty()) {
                return mutRandomOp(s, s.rng, t, parent);
            }
            std::uniform_real_distribution<double> coin(0.0, 1.0);
            if (coin(s.rng) < 0.5) {
                return mutRandomOp(s, s.rng, t, parent);
            }
            auto* node = sf2(7) ? sf2LevelLeafPick(s.rng, t) : noiseGuidedLeafPick(s.rng, t, parent.perPointResid_);
            if (node == nullptr) { return mutRandomOp(s, s.rng, t, parent); }
            if (sf2(7)) { s.sf2GuidedLevel = static_cast<int>(node->Level); }
            std::normal_distribution<Operon::Scalar> dist(Operon::Scalar{0}, Operon::Scalar{3});
            auto const before = node->Value;
            node->Value += dist(s.rng);
            if (opDiagEnabled()) { opDiagRecord(OpMutNoiseGuided, node->Value == before, false); }
            return t;
        };
        if (!config_.useSpectralStructuralMutation) {
            return normalMutate(std::move(tree));
        }
        std::uniform_real_distribution<double> coinF7(0.0, 1.0);
        if (coinF7(s.rng) >= 0.5) {
            return normalMutate(std::move(tree));
        }
        auto spliced = trySpectralSplice(s, tree, parent, isX);
        if (spliced.has_value()) {
            structuralChange = true;
            return std::move(*spliced);
        }
        return normalMutate(std::move(tree));
    }

    static constexpr int kMaxMutationRetries = 3;
    auto mutateChildTreeCapped(SearchState& s, Operon::Tree base, Operon::Tree const& parentTree,
                                PairIndividual const& parent, bool isX, bool& structuralChange,
                                Operon::Tree const* childSibling = nullptr) const
        -> Operon::Tree
    {
        if (config_.mutationOpDiag) { tlLastMutOp = 5; mutDiag_.mutateCalls.fetch_add(1, std::memory_order_relaxed); }
        auto t = mutateChildTree(s, base, parent, isX, structuralChange, childSibling);
        if (config_.gradientMutationEffLog && config_.gradientMutationMode == 0) {
            auto const d0 = nConstDescents_.load(), l0 = nConstLossCalls_.load(), i0 = nConstIters_.load(), ee0 = nEvalErrors_.load();
            Operon::Tree const& sib = (childSibling != nullptr) ? *childSibling : (isX ? parent.y : parent.x);
            gmEffRecord(s, base, t, sib, isX, parent, 0, -1, std::numeric_limits<double>::quiet_NaN(),
                        config_.mutationOpDiag ? tlLastMutOp : -1);
            nConstDescents_ = d0; nConstLossCalls_ = l0; nConstIters_ = i0; nEvalErrors_ = ee0;
        }
        if (!allConstantsFinite(t)) {
            ++nNonFiniteMutatedConst_; nSafetyCounters().nonFiniteMutatedConst.fetch_add(1, std::memory_order_relaxed);
            structuralChange = false;
            t = base;
        }
        if (!config_.enforceStrictMaxLength) { return t; }
        int tries = 0;
        while ((t.Length() > config_.maxLength || t.Depth() > config_.maxDepth)
               && tries < kMaxMutationRetries) {
            ++nMutationOversizeRetries_;
            t = mutateChildTree(s, base, parent, isX, structuralChange, childSibling);
            ++tries;
        }
        if (t.Length() > config_.maxLength || t.Depth() > config_.maxDepth) {
            ++nMutationFallbackToParent_;
            structuralChange = false;
            return parentTree;
        }
        return t;
    }

    void switchToMid(SearchState& s, std::vector<PairIndividual>* popPtr,
                     std::vector<std::vector<PairIndividual>>* layersPtr, int atGen)
    {
        if (!s.lowRes || !s.stagedMid) { return; }
        if (annealMidActive_) {
            dataX_ = midX_;
            dataY_ = midY_;
            dataBlocks_ = midBlocks_;
#if defined(__APPLE__)
            if (config_.useAmxChamfer) {
                amxData_ = midAmxData_;
            }
#endif
        }
        s.stagedMid = false;
        s.annealSwitchGen = atGen;
        if (s.midDs != nullptr) {
            s.activeDs = s.midDs.get();
            s.activeM = s.midM;
            s.range = Operon::Range{0, s.activeM};
        }
        if (popPtr != nullptr && !popPtr->empty()) { reevalAll(s, *popPtr); rankAndCrowd(*popPtr); }
        if (layersPtr != nullptr) {
            for (auto& lay : *layersPtr) {
                if (!lay.empty()) { reevalAll(s, lay); rankAndCrowd(lay); }
            }
        }
        if (config_.useArchive && !s.archive.empty()) {
            std::vector<PairIndividual> old = std::move(s.archive);
            s.archive.clear();
            reevalAll(s, old);
            archiveUpdate(s.archive, old, s.archiveCap);
        }
        if (config_.verbose) {
            std::cout << "  [anneal: switched to mid resolution N=" << dataX_.size()
                      << ", M=" << s.activeM << " at gen " << atGen << "]\n";
        }
    }

    void switchToFull(SearchState& s, std::vector<PairIndividual>* popPtr,
                      std::vector<std::vector<PairIndividual>>* layersPtr, int atGen)
    {
        if (!s.lowRes) { return; }
        if (annealActive_) {
            dataX_ = fullX_;
            dataY_ = fullY_;
            dataBlocks_ = fullBlocks_;
            if (config_.useGridSdf) {
                gridSdf_ = std::move(fullGridSdf_);
            }
#if defined(__APPLE__)
            if (config_.useAmxChamfer) {
                amxData_ = std::move(fullAmxData_);
            }
#endif
        }
        s.lowRes = false;
        s.stagedMid = false;
        s.annealSwitchGen = atGen;
        s.activeDs = s.fullDs.get();
        s.activeM = s.M;
        s.range = Operon::Range{0, s.activeM};
        if (popPtr != nullptr && !popPtr->empty()) { reevalAll(s, *popPtr); rankAndCrowd(*popPtr); }
        if (layersPtr != nullptr) {
            for (auto& lay : *layersPtr) {
                if (!lay.empty()) { reevalAll(s, lay); rankAndCrowd(lay); }
            }
        }
        if (config_.useArchive && !s.archive.empty()) {
            std::vector<PairIndividual> old = std::move(s.archive);
            s.archive.clear();
            reevalAll(s, old);
            archiveUpdate(s.archive, old, s.archiveCap);
        }
        if (config_.verbose) {
            std::cout << "  [anneal: switched to full resolution N=" << dataX_.size()
                      << ", M=" << s.activeM << " at gen " << atGen << "]\n";
        }
    }

public:
    void SetDispatchTable(DTable const* dt)
    {
        dtable_ = dt;
#if defined(__APPLE__)
        if (config_.useAccelerateTranscendentals) {
            accelDTable_ = *dtable_;
            RegisterAccelerateTranscendentals(*accelDTable_);
            dtable_ = &(*accelDTable_);
        }
#endif
    }

    void InitSearchState()
    {
#ifdef __APPLE__
        pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
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

        nEvalErrors_ = 0;
#ifdef SKETCH_DUP_STATS
        dupStats_.reset();
#endif
        nConstDescents_ = 0; nConstLossCalls_ = 0; nConstIters_ = 0;
        nLmStopUnresolvableDamped_ = 0; nLmStopRealizedSmall_ = 0; nLmStopFixedPoint_ = 0;
        for (auto& h : lmIterHist_) { h = 0; }
        if (config_.useConcurrentCpuGpuLM) {
            cclmStats_.clear(); cclmParity_.clear(); cclmRateCpu_ = 0.0; cclmRateGpu_ = 0.0;
        }
        nNonFiniteFitness_ = 0;
        nInvalidChildren_ = 0;
        nNonFiniteMutatedConst_ = 0;
        nNonFiniteLmConst_ = 0;
        nMutationOversizeRetries_ = 0;
        nMutationFallbackToParent_ = 0;
        nRejectedNonFiniteBeforeConstOpt_ = 0;
        nLmRecoveredFromNonFiniteStart_ = 0;
        nNonFiniteWithLog_ = 0; nNonFiniteWithSqrt_ = 0; nNonFiniteWithDiv_ = 0;
        nNonFiniteWithExp_ = 0; nNonFiniteWithNoneOfThese_ = 0;
        nSemCross_ = 0; nSemFallback_ = 0;
        nCascadeCullGens_ = 0; nCascadeCullSeen_ = 0; nCascadeCullSurvived_ = 0;
        nLiveCoarseSieveSeen_ = 0; nLiveCoarseSieveCulled_ = 0;
        liveCoarseSieveBoundMin_ = std::numeric_limits<double>::infinity();
        liveCoarseSieveBoundMax_ = -std::numeric_limits<double>::infinity();
        liveCoarseSieveBoundSum_ = 0.0;
        nLiveCoarseSieveFinite_ = 0;

        state_ = std::make_unique<SearchState>(config_.seed);
        auto& s = *state_;
        s.fitStart = std::chrono::steady_clock::now();

        auto const M = config_.gridSize;
        auto const N = config_.populationSize;
        s.M = M;
        s.N = N;

        s.lowM = (config_.annealLowM >= 2 && config_.annealLowM < M) ? config_.annealLowM : M;
        s.midM = (config_.annealMidM >= 2 && config_.annealMidM < M) ? config_.annealMidM : M;

        s.fullFlat.resize(M * 2);
        for (std::size_t i = 0; i < M; ++i) {
            s.fullFlat[i] = static_cast<Operon::Scalar>(i) / static_cast<Operon::Scalar>(M - 1);
            s.fullFlat[M + i] = 0;
        }
        s.fullDs = std::make_unique<Operon::Dataset>(gsl::not_null{s.fullFlat.data()}, static_cast<int>(M), 2);

        if (s.lowM < M) {
            s.lowFlat.resize(s.lowM * 2);
            for (std::size_t i = 0; i < s.lowM; ++i) {
                s.lowFlat[i] = static_cast<Operon::Scalar>(i) / static_cast<Operon::Scalar>(s.lowM - 1);
                s.lowFlat[s.lowM + i] = 0;
            }
            s.lowDs = std::make_unique<Operon::Dataset>(gsl::not_null{s.lowFlat.data()}, static_cast<int>(s.lowM), 2);
        }

        if (s.midM < M && s.midM != s.lowM) {
            s.midFlat.resize(s.midM * 2);
            for (std::size_t i = 0; i < s.midM; ++i) {
                s.midFlat[i] = static_cast<Operon::Scalar>(i) / static_cast<Operon::Scalar>(s.midM - 1);
                s.midFlat[s.midM + i] = 0;
            }
            s.midDs = std::make_unique<Operon::Dataset>(gsl::not_null{s.midFlat.data()}, static_cast<int>(s.midM), 2);
        }

        auto vars = s.fullDs->GetVariables();

        s.prob = std::make_unique<Operon::Problem>(s.fullDs.get());
        s.prob->SetTrainingRange(0, static_cast<int>(M));
        s.prob->SetTarget(vars[1].Hash);
        s.prob->SetInputs(std::vector{vars[0].Hash});

        using NT = Operon::NodeType;
        using BO = Operon::BuiltinOp;
        auto mask = NT::Constant | NT::Variable
                  | BO::Add | BO::Sub | BO::Mul | BO::Div
                  | BO::Sin | BO::Cos | BO::Exp | BO::Log | BO::Sqrt | BO::Square;
        s.prob->ConfigurePrimitiveSet(mask);

        auto& pset = s.prob->GetPrimitiveSet();
        s.inputs = {vars[0].Hash};

        if (config_.useOperatorFamilyBias && !config_.usePriorMinimalComplexity) {
            applyOperatorFamilyBias(pset, familyWeightsFromProbe(familyProbe(dataX_, dataY_), sf2(8) ? dataX_.size() : 0));
        }

        s.creator = std::make_unique<Operon::BalancedTreeCreator>(&pset, s.inputs, 0.0, config_.maxLength);
        s.coeffInit.ParameterizeDistribution(Operon::Scalar{0},
                                           static_cast<Operon::Scalar>(config_.coeffSpread));
        s.crossover = std::make_unique<Operon::SubtreeCrossover>(0.9, config_.maxDepth, config_.maxLength);
        if (sf2(12)) {
            double const Mm1 = static_cast<double>(std::max<std::size_t>(2, config_.gridSize) - 1);
            s.sf2CoeffInit.logLo = std::log(1.0 / Mm1);
            s.sf2CoeffInit.logHi = std::log(M_PI * Mm1);
            s.coeffInitP = &s.sf2CoeffInit;
        }
        if (opDiagEnabled()) { config_.mutationOpDiag = true; }
        s.mut = std::make_unique<MutationSet>(pset, *s.creator, *s.coeffInitP, config_.maxDepth, config_.maxLength,
                                              config_.binaryInsertSubtree,
                                              config_.mutationOpDiag ? &mutDiag_ : nullptr,
                                              op2(0));

        s.lowRes = annealActive_ || (s.lowDs != nullptr);
        s.stagedMid = annealMidActive_ || (s.midDs != nullptr);
        s.annealSwitchGen = -1;
        s.switchGen1 = s.lowRes
            ? static_cast<std::size_t>(std::floor(config_.annealResolution
                                                  * static_cast<double>(config_.generations)))
            : std::numeric_limits<std::size_t>::max();
        s.switchGen2 = (s.lowRes && s.stagedMid)
            ? static_cast<std::size_t>(std::floor(config_.annealStage2Resolution
                                                  * static_cast<double>(config_.generations)))
            : s.switchGen1;

        s.activeDs = (s.lowDs != nullptr && s.lowRes) ? s.lowDs.get() : s.fullDs.get();
        s.activeM = (s.lowDs != nullptr && s.lowRes) ? s.lowM : M;
        s.range = Operon::Range{0, s.activeM};

        auto makeTree = [&]() {
            auto t = (*s.creator)(s.rng, 5, 1, config_.maxDepth);
            (*s.coeffInitP)(s.rng, t);
            return t;
        };

        auto const nWorkers = workerCount();
        s.pop.resize(N);
        for (auto& ind : s.pop) {
            ind.x = makeTree();
            ind.y = makeTree();
        }
        std::vector<std::vector<Operon::Scalar>> initGridsX, initGridsY;
        if (config_.useBatchedRegisterEval) {
            evaluateBatchedPrepass(s.pop, *s.activeDs, s.range, initGridsX, initGridsY);
        }
        std::vector<std::vector<Operon::Scalar>> initPsGridsX, initPsGridsY;
        CisPointSegBatch initPsBatch;
        CisPointSegBatch initPsBatch48;
        PsStreamGen initPsGen;
        bool const initPsBatchOk = (psStreamWanted() && psStreamBegin(s.pop.size(), s.range, initPsGen))
            ? [&] {
                  for (std::size_t i = 0; i < s.pop.size(); ++i) {
                      psStreamAdd(s.pop[i], i, *s.activeDs, s.range, initPsGridsX, initPsGridsY, initPsGen);
                      psStreamPump(initPsGen);
                  }
                  return psStreamFinish(s.pop.size(), initPsGen, initPsBatch);
              }()
            : (config_.useBatchedGpuPointSeg && cisFitnessMode())
                ? evaluateBatchedPointSegPrepass(s.pop, *s.activeDs, s.range, initPsGridsX, initPsGridsY, initPsBatch, &initPsBatch48)
                : false;
        parallelApply(s.pop.size(), nWorkers, [&](std::size_t i) {
            if (initPsBatchOk && i < initPsGridsX.size()) {
                evaluateGuardedFromGrid(s.pop[i], *s.activeDs, s.range, initPsGridsX[i], initPsGridsY[i], &initPsBatch, i, &initPsBatch48);
            } else if (config_.useBatchedRegisterEval && i < initGridsX.size()) {
                evaluateGuardedFromGrid(s.pop[i], *s.activeDs, s.range, initGridsX[i], initGridsY[i]);
            } else {
                evaluateGuarded(s.pop[i], *s.activeDs, s.range);
            }
        });
        rankAndCrowd(s.pop);

        if (!eliteSeeds_.empty()) {
            std::size_t const nSeed = std::min(eliteSeeds_.size(), s.pop.size());
            for (std::size_t i = 0; i < nSeed; ++i) {
                s.pop[s.pop.size() - 1 - i] = eliteSeeds_[i];
            }
            parallelApply(nSeed, nWorkers, [&](std::size_t i) {
                evaluateGuarded(s.pop[s.pop.size() - 1 - i], *s.activeDs, s.range);
            });
            rankAndCrowd(s.pop);
        }

        s.archiveCap = kArchiveCapMultiple * N;
        if (config_.useArchive) { archiveUpdate(s.archive, s.pop, s.archiveCap); }
        if (config_.useComplexityHallOfFame) { hallOfFameUpdate(s.hallOfFame, s.pop, kHallOfFameK); }

        if (config_.alpsLayers > 0) {
            std::size_t const Lay = config_.alpsLayers;
            std::size_t const cap = std::max<std::size_t>(N / Lay, 8);
            s.layers.resize(Lay);
            auto seeded = selectNext(s.pop, std::min(cap, s.pop.size()));
            for (auto& ind : seeded) { ind.age = 0; s.layers[0].push_back(std::move(ind)); }
        }
    }

    void SetThreads(std::size_t t) { config_.nThreads = t; }

    void SetEliteSeeds(std::vector<PairIndividual> seeds) { eliteSeeds_ = std::move(seeds); }

    void SetExternalPointSigma2(std::vector<Operon::Scalar> sigma2) { externalPointSigma2_ = std::move(sigma2); }

    [[nodiscard]] auto ElitePopulation(std::size_t k) const -> std::vector<PairIndividual>
    {
        std::vector<PairIndividual> out;
        if (!state_) { return out; }
        auto const& s = *state_;
        std::vector<PairIndividual const*> ptrs;
        ptrs.reserve(s.pop.size());
        for (auto const& ind : s.pop) { if (ind.rank == 0) { ptrs.push_back(&ind); } }
        std::stable_sort(ptrs.begin(), ptrs.end(), [](auto const* a, auto const* b) {
            return a->shapeErr < b->shapeErr;
        });
        std::size_t const n = std::min(k, ptrs.size());
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i) { out.push_back(*ptrs[i]); }
        return out;
    }

    void RunTo(std::size_t targetGen)
    {
#ifdef __APPLE__
        pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
        if (!state_) {
            InitSearchState();
        }
        auto& s = *state_;
        if (s.truncated) { return; }

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
        if (s.inSimplify) {
            errorClampValue() = std::max(config_.floorClamp, config_.floorStop);
        }

        auto overBudget = [&]() -> bool {
            if (config_.maxFitSeconds <= 0) { return false; }
            std::chrono::duration<double> const el = std::chrono::steady_clock::now() - s.fitStart;
            return el.count() > config_.maxFitSeconds;
        };

        auto const M = s.M;
        auto const N = s.N;
        auto const nWorkers = workerCount();

        std::uniform_real_distribution<double> uni(0.0, 1.0);

        auto makeTree = [&]() {
            auto t = (*s.creator)(s.rng, 5, 1, config_.maxDepth);
            (*s.coeffInitP)(s.rng, t);
            return t;
        };

        constexpr std::size_t kSemCandidates = 8;
        constexpr std::size_t kSemSpreadSample = 16;
        constexpr std::size_t kSemSamples = 32;
        std::size_t const semStride = std::max<std::size_t>(1, M / kSemSamples);

        auto semanticsOf = [&](Operon::Tree const& t) {
            std::vector<Operon::Scalar> out;
            auto const coeffs = t.GetCoefficients();
            Operon::Interpreter<Operon::Scalar, DTable> in(dtable_, s.activeDs, &t);
            auto v = in.Evaluate({coeffs.data(), coeffs.size()}, s.range);
            for (std::size_t k = 0; k < v.size(); k += semStride) { out.push_back(v[k]); }
            return out;
        };
        auto subtreeOf = [](Operon::Tree const& t, std::size_t i) {
            auto const& nodes = t.Nodes();
            auto const len = static_cast<std::ptrdiff_t>(nodes[i].Length + 1);
            auto const end = static_cast<std::ptrdiff_t>(i) + 1;
            Operon::Tree st(std::vector<Operon::Node>(nodes.begin() + (end - len), nodes.begin() + end));
            st.UpdateNodes();
            return st;
        };
        auto semDist = [](std::vector<Operon::Scalar> const& a,
                          std::vector<Operon::Scalar> const& b) -> double {
            constexpr auto INF = std::numeric_limits<double>::infinity();
            if (a.size() != b.size() || a.empty()) { return INF; }
            double acc = 0.0;
            for (std::size_t k = 0; k < a.size(); ++k) {
                double const d = static_cast<double>(a[k]) - static_cast<double>(b[k]);
                if (!std::isfinite(d)) { return INF; }
                acc += d * d;
            }
            return std::sqrt(acc / static_cast<double>(a.size()));
        };
        auto updateSemBound = [&](std::vector<PairIndividual> const& p) {
            std::vector<std::vector<Operon::Scalar>> S;
            for (std::size_t i = 0; i < p.size() && i < kSemSpreadSample; ++i) {
                S.push_back(semanticsOf(p[i].x));
                S.push_back(semanticsOf(p[i].y));
            }
            std::vector<double> d;
            for (std::size_t i = 0; i < S.size(); ++i) {
                for (std::size_t j = i + 1; j < S.size(); ++j) {
                    double const v = semDist(S[i], S[j]);
                    if (std::isfinite(v)) { d.push_back(v); }
                }
            }
            if (d.empty()) { s.semBound = 0.0; return; }
            auto mid = d.begin() + static_cast<std::ptrdiff_t>(d.size() / 2);
            std::nth_element(d.begin(), mid, d.end());
            s.semBound = *mid;
        };
        auto ssCrossover = [&](Operon::RandomGenerator& r, Operon::Tree const& A,
                               Operon::Tree const& B) -> Operon::Tree {
            auto const& an = A.Nodes();
            auto const& bn = B.Nodes();
            if (an.size() < 2 || bn.empty() || !(s.semBound > 0.0)) { return (*s.crossover)(r, A, B); }
            ++nSemCross_;
            std::uniform_int_distribution<std::size_t> pa(0, an.size() - 2);
            std::size_t const i = pa(r);
            std::size_t const iLen = an[i].Length + 1;
            auto const target = semanticsOf(subtreeOf(A, i));
            std::uniform_int_distribution<std::size_t> pb(0, bn.size() - 1);
            std::size_t bestJ = bn.size();
            double bestD = std::numeric_limits<double>::infinity();
            for (std::size_t c = 0; c < kSemCandidates; ++c) {
                std::size_t const j = pb(r);
                std::size_t const jLen = bn[j].Length + 1;
                if (an.size() - iLen + jLen > config_.maxLength) { continue; }
                double const d = semDist(target, semanticsOf(subtreeOf(B, j)));
                if (d > 0.0 && d <= s.semBound && d < bestD) { bestD = d; bestJ = j; }
            }
            if (bestJ == bn.size()) { ++nSemFallback_; return (*s.crossover)(r, A, B); }
            auto const jLen = static_cast<std::ptrdiff_t>(bn[bestJ].Length + 1);
            auto const jEnd = static_cast<std::ptrdiff_t>(bestJ) + 1;
            auto const iEnd = static_cast<std::ptrdiff_t>(i) + 1;
            std::vector<Operon::Node> out;
            out.reserve(an.size() + static_cast<std::size_t>(jLen));
            out.insert(out.end(), an.begin(), an.begin() + (iEnd - static_cast<std::ptrdiff_t>(iLen)));
            out.insert(out.end(), bn.begin() + (jEnd - jLen), bn.begin() + jEnd);
            out.insert(out.end(), an.begin() + iEnd, an.end());
            Operon::Tree child(out);
            child.UpdateNodes();
            if (child.Depth() > config_.maxDepth) { ++nSemFallback_; return (*s.crossover)(r, A, B); }
            return child;
        };
        auto xover = [&](Operon::RandomGenerator& r, Operon::Tree const& A, Operon::Tree const& B) {
            return config_.semanticCrossover ? ssCrossover(r, A, B) : (*s.crossover)(r, A, B);
        };

        auto enterSimplify = [&](std::vector<PairIndividual>& p) {
            s.inSimplify = true;
            s.simplifyLeft = config_.floorSimplifyGens;
            errorClampValue() = std::max(config_.floorClamp, config_.floorStop);
            rankAndCrowd(p);
            if (config_.verbose) {
                auto* b = bestAccuracy(p);
                std::cout << "  [floor reached: best_err=" << (b != nullptr ? b->shapeErr : 0)
                          << " <= " << config_.floorStop << "] simplifying for "
                          << config_.floorSimplifyGens << " gens at clamp "
                          << errorClampValue() << "\n";
            }
        };

        targetGen = std::min(targetGen, config_.generations);

        if (config_.alpsLayers == 0) {
            std::vector<PairIndividual> offspring;
            offspring.reserve(N);
            std::vector<char> wantOpt;
            wantOpt.reserve(N);
            std::vector<PairIndividual> combined;
            combined.reserve(2 * N);
            std::vector<PairIndividual> nextPop;
            nextPop.reserve(N);

            bool const islandsActive =
                config_.useIslands && config_.numIslands > 1 && N > 1;
            std::vector<std::pair<std::size_t, std::size_t>> islRanges;
            std::vector<std::size_t> islId;
            if (islandsActive) {
                std::size_t const k =
                    std::min<std::size_t>(static_cast<std::size_t>(config_.numIslands), N);
                islRanges = islandRanges(N, k);
                islId.assign(N, 0);
                for (std::size_t i = 0; i < islRanges.size(); ++i) {
                    for (std::size_t j = islRanges[i].first; j < islRanges[i].second; ++j) {
                        islId[j] = i;
                    }
                }
            }
            std::vector<PairIndividual> islCombined;
            std::vector<PairIndividual> islOut;

            if (config_.gradientMutationTrajectory && s.currentGen < targetGen) { gmTrajectoryPoint(s.pop, s); }
            for (; s.currentGen < targetGen; ++s.currentGen) {
                std::size_t gen = s.currentGen;
                s.lastGen = static_cast<int>(gen);
                std::chrono::steady_clock::time_point profGenT0;
                if (config_.profileTiming) { profGenT0 = std::chrono::steady_clock::now(); }

                if (s.lowRes && s.stagedMid && gen >= s.switchGen1) { switchToMid(s, &s.pop, nullptr, static_cast<int>(gen)); }
                else if (s.lowRes && gen >= s.switchGen2) { switchToFull(s, &s.pop, nullptr, static_cast<int>(gen)); }

                if (config_.semanticCrossover) { updateSemBound(s.pop); }

                std::chrono::steady_clock::time_point profBreedT0;
                if (config_.profileTiming) { profBreedT0 = std::chrono::steady_clock::now(); }
                offspring.clear();
                wantOpt.clear();
                std::vector<std::vector<Operon::Scalar>> psGridsX, psGridsY;
                PsStreamGen psGen;
                if (config_.gpuPointSegOverlap) { (void)psStreamBegin(N, s.range, psGen); }
                PsOverlapHelper psHelper;
                if (psGen.active && config_.gpuPointSegOverlapHelper && offspring.capacity() >= N) {
                    psGridsX.resize(N); psGridsY.resize(N);
                    psHelper.start(offspring.data(), [this, &s, &psGridsX, &psGridsY, &psGen](PairIndividual const& ind, std::size_t i) {
                        psStreamAdd(ind, i, *s.activeDs, s.range, psGridsX, psGridsY, psGen);
                        psStreamPump(psGen);
                    });
                }
                static bool const measureTreeSizes = std::getenv("SKETCH_MEASURE_TREE_SIZES") != nullptr;
                struct TreeSizeGenStats {
                    std::size_t total{0}, overMaxLen{0};
                    std::size_t overFromMutation{0}, overFromCrossover{0};
                    std::size_t maxSeen{0};
                } tsStats;
                bool const gmDefer = config_.gradientMutationMode != 0 && !psGen.active;
                std::vector<GmPendingChild> gmPend;
                std::vector<GmBanditRec> gmBandit;
                bool const gmIr = gmDefer && config_.gradientMutationInRunErr;
                bool const gmRefitRw = gmDefer && config_.gradientMutationMode == 3 && config_.gradientMutationBanditRefitReward;
                std::vector<GmIrStub> gmIrStubs;
                std::vector<char> gmRefitMask;
                if (sf2(5)) { sf2ComputeLexEps(s.pop); }
                bool const diagOnGen = opDiagEnabled();
                bool const parallelBreedOk = config_.useParallelBreed
                    && !diagOnGen && !measureTreeSizes && gmDefer && config_.gradientMutationMode == 1
                    && !config_.semanticCrossover && !psGen.active && config_.lmTargeted
                    && !psHelper.running();
                if (parallelBreedOk) {
                    offspring.resize(N);
                    std::vector<GmPendingChild> gmPendSlots(N);
                    bool const lexParentsCap = config_.useLexicase && !config_.useCoherentSelection;
                    parallelApply(N, workerCount(), [&](std::size_t i) {
                        bool const lexParents = lexParentsCap;
                        auto const seedRng = parallelBreedSeed(config_.seed, gen, i, 0);
                        auto const seedSel = parallelBreedSeed(config_.seed, gen, i, 1);
                        Operon::RandomGenerator rngLocal(seedRng);
                        std::mt19937 selLocal(static_cast<std::uint32_t>(seedSel));
                        std::uniform_real_distribution<double> uniLocal(0.0, 1.0);
                        PairIndividual const* ap;
                        PairIndividual const* bp;
                        if (islandsActive) {
                            auto const& rg = islRanges[islId[i]];
                            ap = lexParents ? &epsilonLexicaseSelect(s.pop, selLocal, rg.first, rg.second)
                                            : &tournament(s.pop, selLocal, rg.first, rg.second);
                            bp = lexParents ? &epsilonLexicaseSelect(s.pop, selLocal, rg.first, rg.second)
                                            : &tournament(s.pop, selLocal, rg.first, rg.second);
                        } else {
                            ap = lexParents ? &epsilonLexicaseSelect(s.pop, selLocal) : &tournament(s.pop, selLocal);
                            bp = lexParents ? &epsilonLexicaseSelect(s.pop, selLocal) : &tournament(s.pop, selLocal);
                        }
                        auto const& a = *ap;
                        auto const& b = *bp;
                        auto xoverLocal = [&](Operon::Tree const& A, Operon::Tree const& B) {
                            return (*s.crossover)(rngLocal, A, B);
                        };
                        PairIndividual child;
                        auto breedOnce = [&]() {
                            if (config_.crossTreeProb > Operon::Scalar{0}) {
                                auto const p = config_.crossTreeProb;
                                auto const& dx = (uniLocal(selLocal) < p) ? b.y : b.x;
                                auto const& dy = (uniLocal(selLocal) < p) ? b.x : b.y;
                                bool const cx = uniLocal(selLocal) < 0.9;
                                bool const cy = uniLocal(selLocal) < 0.9;
                                child.x = cx ? xoverLocal(a.x, dx) : a.x;
                                child.y = cy ? xoverLocal(a.y, dy) : a.y;
                                bool const doExch = uniLocal(selLocal) < p;
                                if (doExch) {
                                    bool const exchY = uniLocal(selLocal) < 0.5;
                                    if (exchY) { child.y = xoverLocal(child.y, child.x); }
                                    else       { child.x = xoverLocal(child.x, child.y); }
                                }
                            } else {
                                bool const cx = uniLocal(selLocal) < 0.9;
                                bool const cy = uniLocal(selLocal) < 0.9;
                                child.x = cx ? xoverLocal(a.x, b.x) : a.x;
                                child.y = cy ? xoverLocal(a.y, b.y) : a.y;
                            }
                        };
                        breedOnce();
                        if (op2(1)) {
                            std::vector<Operon::Scalar> localBufX, localBufY;
                            auto finiteXY = [&]() -> bool {
                                if (s.activeDs == nullptr) { return true; }
                                interpInto(child.x, child.x.GetCoefficients(), *s.activeDs, s.range, localBufX);
                                interpInto(child.y, child.y.GetCoefficients(), *s.activeDs, s.range, localBufY);
                                for (auto v : localBufX) { if (!std::isfinite(static_cast<double>(v))) { return false; } }
                                for (auto v : localBufY) { if (!std::isfinite(static_cast<double>(v))) { return false; } }
                                return true;
                            };
                            bool finite = finiteXY();
                            while (!finite) { breedOnce(); finite = finiteXY(); }
                        }
                        GmPendingChild gmPc;
                        bool const doMutX = uniLocal(selLocal) < 0.25;
                        if (doMutX) { gmPc.ev.push_back({true, rngLocal(), gmBanditPick(rngLocal)}); }
                        bool const doMutY = uniLocal(selLocal) < 0.25;
                        if (doMutY) { gmPc.ev.push_back({false, rngLocal(), gmBanditPick(rngLocal)}); }
                        if (child.x.Length() > 2 * config_.maxLength ||
                            child.y.Length() > 2 * config_.maxLength) { child = a; }
                        if (!structurallyValid(child)) { child = a; ++nInvalidChildren_; }
                        if (sf2_ != 0) {
                            child.sf2LvX = -1; child.sf2LvY = -1;
                            child.sf2ParentErr = a.shapeErr;
                            child.sf2Child = true;
                        }
                        if (config_.bracketNewtonConfig.warmStart) {
                            child.warmSegmentIndex = a.warmSegmentIndex;
                            child.warmSegmentT = a.warmSegmentT;
                            child.warmCompletedSweeps = a.warmCompletedSweeps;
                        } else {
                            child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                            child.warmCompletedSweeps = 0;
                        }
                        if (!gmPc.ev.empty()) {
                            gmPc.idx = i;
                            gmPc.parent = &a;
                            gmPendSlots[i] = std::move(gmPc);
                        }
                        offspring[i] = std::move(child);
                    });
                    for (auto& pc : gmPendSlots) { if (!pc.ev.empty()) { gmPend.push_back(std::move(pc)); } }
                } else {
                while (offspring.size() < N) {
                    bool const lexParents = config_.useLexicase && !config_.useCoherentSelection;
                    PairIndividual const* ap;
                    PairIndividual const* bp;
                    if (islandsActive) {
                        auto const& rg = islRanges[islId[offspring.size()]];
                        ap = lexParents ? &epsilonLexicaseSelect(s.pop, s.sel, rg.first, rg.second)
                                        : &tournament(s.pop, s.sel, rg.first, rg.second);
                        bp = lexParents ? &epsilonLexicaseSelect(s.pop, s.sel, rg.first, rg.second)
                                        : &tournament(s.pop, s.sel, rg.first, rg.second);
                    } else {
                        ap = lexParents ? &epsilonLexicaseSelect(s.pop, s.sel) : &tournament(s.pop, s.sel);
                        bp = lexParents ? &epsilonLexicaseSelect(s.pop, s.sel) : &tournament(s.pop, s.sel);
                    }
                    auto const& a = *ap;
                    auto const& b = *bp;
                    bool const diagOn = opDiagEnabled();
                    Operon::Tree diagPreCrossX, diagPreCrossY;
                    Operon::Tree diagPostCrossX, diagPostCrossY;
                    bool diagCrossXApplied = false, diagCrossYApplied = false;
                    bool diagExchApplied = false, diagExchTouchedX = false;
                    PairIndividual child;
                    if (config_.crossTreeProb > Operon::Scalar{0}) {
                        auto const p = config_.crossTreeProb;
                        auto const& dx = (uni(s.sel) < p) ? b.y : b.x;
                        auto const& dy = (uni(s.sel) < p) ? b.x : b.y;
                        bool const cx = uni(s.sel) < 0.9;
                        bool const cy = uni(s.sel) < 0.9;
                        if (diagOn) { diagPreCrossX = a.x; diagPreCrossY = a.y; diagCrossXApplied = cx; diagCrossYApplied = cy; }
                        child.x = cx ? xover(s.rng, a.x, dx) : a.x;
                        child.y = cy ? xover(s.rng, a.y, dy) : a.y;
                        if (diagOn) { diagPostCrossX = child.x; diagPostCrossY = child.y; }
                        bool const doExch = uni(s.sel) < p;
                        if (doExch) {
                            bool const exchY = uni(s.sel) < 0.5;
                            if (exchY) { child.y = xover(s.rng, child.y, child.x); }
                            else       { child.x = xover(s.rng, child.x, child.y); }
                            if (diagOn) { diagExchApplied = true; diagExchTouchedX = !exchY; }
                        }
                    } else {
                        bool const cx = uni(s.sel) < 0.9;
                        bool const cy = uni(s.sel) < 0.9;
                        if (diagOn) { diagPreCrossX = a.x; diagPreCrossY = a.y; diagCrossXApplied = cx; diagCrossYApplied = cy; }
                        child.x = cx ? xover(s.rng, a.x, b.x) : a.x;
                        child.y = cy ? xover(s.rng, a.y, b.y) : a.y;
                        if (diagOn) { diagPostCrossX = child.x; diagPostCrossY = child.y; }
                    }
                    if (op2(1)) {
                        auto redrawCross = [&]() {
                            if (config_.crossTreeProb > Operon::Scalar{0}) {
                                auto const p = config_.crossTreeProb;
                                auto const& dx = (uni(s.sel) < p) ? b.y : b.x;
                                auto const& dy = (uni(s.sel) < p) ? b.x : b.y;
                                bool const cx = uni(s.sel) < 0.9;
                                bool const cy = uni(s.sel) < 0.9;
                                child.x = cx ? xover(s.rng, a.x, dx) : a.x;
                                child.y = cy ? xover(s.rng, a.y, dy) : a.y;
                                bool const doExch = uni(s.sel) < p;
                                if (doExch) {
                                    bool const exchY = uni(s.sel) < 0.5;
                                    if (exchY) { child.y = xover(s.rng, child.y, child.x); }
                                    else       { child.x = xover(s.rng, child.x, child.y); }
                                }
                            } else {
                                bool const cx = uni(s.sel) < 0.9;
                                bool const cy = uni(s.sel) < 0.9;
                                child.x = cx ? xover(s.rng, a.x, b.x) : a.x;
                                child.y = cy ? xover(s.rng, a.y, b.y) : a.y;
                            }
                        };
                        thread_local std::vector<Operon::Scalar> op2BufX, op2BufY;
                        auto finiteXY = [&]() -> bool {
                            if (s.activeDs == nullptr) { return true; }
                            interpInto(child.x, child.x.GetCoefficients(), *s.activeDs, s.range, op2BufX);
                            interpInto(child.y, child.y.GetCoefficients(), *s.activeDs, s.range, op2BufY);
                            for (auto v : op2BufX) { if (!std::isfinite(static_cast<double>(v))) { return false; } }
                            for (auto v : op2BufY) { if (!std::isfinite(static_cast<double>(v))) { return false; } }
                            return true;
                        };
                        bool const diagEv = op2DiagEnabled();
                        auto recordCrossDraw = [&](bool finite) {
                            if (!diagEv) { return; }
                            auto& oc = op2ValidCounters();
                            oc.draws[5].fetch_add(1, std::memory_order_relaxed);
                            oc.draws[6].fetch_add(1, std::memory_order_relaxed);
                            if (finite) {
                                oc.events[5].fetch_add(1, std::memory_order_relaxed);
                                oc.events[6].fetch_add(1, std::memory_order_relaxed);
                            }
                        };
                        bool finite = finiteXY();
                        recordCrossDraw(finite);
                        while (!finite) {
                            redrawCross();
                            finite = finiteXY();
                            recordCrossDraw(finite);
                        }
                    }
                    bool f7X = false, f7Y = false;
                    GmPendingChild gmPc;
                    int sf2LvX = -1, sf2LvY = -1;
                    int diagMutXSlot = -1, diagMutYSlot = -1;
                    bool const doMutX = uni(s.sel) < 0.25;
                    if (doMutX) {
                        if (gmDefer) { gmPc.ev.push_back({true, s.rng(), gmBanditPick(s.rng)}); }
                        else {
                            s.sf2GuidedLevel = -1;
                            child.x = mutateChildTreeCapped(s, child.x, a.x, a, true, f7X, &child.y);
                            sf2LvX = s.sf2GuidedLevel;
                            if (diagOn) { diagMutXSlot = config_.gradientMutationMode != 0 ? -2 : tlLastMutOp; }
                        }
                    }
                    bool const doMutY = uni(s.sel) < 0.25;
                    if (doMutY) {
                        if (gmDefer) { gmPc.ev.push_back({false, s.rng(), gmBanditPick(s.rng)}); }
                        else {
                            s.sf2GuidedLevel = -1;
                            child.y = mutateChildTreeCapped(s, child.y, a.y, a, false, f7Y, &child.x);
                            sf2LvY = s.sf2GuidedLevel;
                            if (diagOn) { diagMutYSlot = config_.gradientMutationMode != 0 ? -2 : tlLastMutOp; }
                        }
                    }
                    if (diagOn) {
                        if (diagCrossXApplied) { opDiagRecord(OpCrossoverX, sameTreeStructure(diagPreCrossX, diagPostCrossX), false); }
                        if (diagCrossYApplied) { opDiagRecord(OpCrossoverY, sameTreeStructure(diagPreCrossY, diagPostCrossY), false); }
                        if (diagExchApplied) {
                            bool const exchNoop = diagExchTouchedX
                                ? sameTreeStructure(diagPostCrossX, child.x)
                                : sameTreeStructure(diagPostCrossY, child.y);
                            opDiagRecord(OpCrossTreeExchange, exchNoop, false);
                        }
                        std::vector<Operon::Scalar> diagGridX, diagGridY;
                        bool const xFinite = [&] {
                            interpInto(child.x, child.x.GetCoefficients(), *s.activeDs, s.range, diagGridX);
                            for (auto v : diagGridX) { if (!std::isfinite(v)) { return false; } }
                            return true;
                        }();
                        bool const yFinite = [&] {
                            interpInto(child.y, child.y.GetCoefficients(), *s.activeDs, s.range, diagGridY);
                            for (auto v : diagGridY) { if (!std::isfinite(v)) { return false; } }
                            return true;
                        }();
                        auto mutOpIdFor = [](int slot) -> int {
                            if (slot >= 0 && slot <= 4) { return kMutSlotToOpDiagId[static_cast<std::size_t>(slot)]; }
                            if (slot == 5) { return OpMutNoiseGuided; }
                            return OpMutGradient;
                        };
                        if (!xFinite || !yFinite) {
                            if (doMutX && diagMutXSlot != -1) { opDiagCounters().invalid[static_cast<std::size_t>(mutOpIdFor(diagMutXSlot))].fetch_add(!xFinite ? 1 : 0, std::memory_order_relaxed); }
                            else if (diagExchApplied && diagExchTouchedX) { opDiagCounters().invalid[OpCrossTreeExchange].fetch_add(!xFinite ? 1 : 0, std::memory_order_relaxed); }
                            else if (diagCrossXApplied) { opDiagCounters().invalid[OpCrossoverX].fetch_add(!xFinite ? 1 : 0, std::memory_order_relaxed); }
                            if (doMutY && diagMutYSlot != -1) { opDiagCounters().invalid[static_cast<std::size_t>(mutOpIdFor(diagMutYSlot))].fetch_add(!yFinite ? 1 : 0, std::memory_order_relaxed); }
                            else if (diagExchApplied && !diagExchTouchedX) { opDiagCounters().invalid[OpCrossTreeExchange].fetch_add(!yFinite ? 1 : 0, std::memory_order_relaxed); }
                            else if (diagCrossYApplied) { opDiagCounters().invalid[OpCrossoverY].fetch_add(!yFinite ? 1 : 0, std::memory_order_relaxed); }
                        }
                    }
                    if (measureTreeSizes) {
                        auto const lx = child.x.Length(), ly = child.y.Length();
                        ++tsStats.total;
                        tsStats.maxSeen = std::max({tsStats.maxSeen, lx, ly});
                        if (lx > config_.maxLength) {
                            ++tsStats.overMaxLen;
                            (doMutX ? tsStats.overFromMutation : tsStats.overFromCrossover)++;
                        }
                        if (ly > config_.maxLength) {
                            ++tsStats.overMaxLen;
                            (doMutY ? tsStats.overFromMutation : tsStats.overFromCrossover)++;
                        }
                    }
                    if (child.x.Length() > 2 * config_.maxLength ||
                        child.y.Length() > 2 * config_.maxLength) { child = a; f7X = f7Y = false; sf2LvX = sf2LvY = -1; }
                    if (!structurallyValid(child)) { child = a; ++nInvalidChildren_; f7X = f7Y = false; sf2LvX = sf2LvY = -1; }
                    if (sf2_ != 0) {
                        child.sf2LvX = sf2LvX; child.sf2LvY = sf2LvY;
                        child.sf2ParentErr = a.shapeErr;
                        child.sf2Child = true;
                    }
                    if (config_.bracketNewtonConfig.warmStart) {
                        child.warmSegmentIndex = a.warmSegmentIndex;
                        child.warmSegmentT = a.warmSegmentT;
                        child.warmCompletedSweeps = a.warmCompletedSweeps;
                    } else {
                        child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                        child.warmCompletedSweeps = 0;
                    }
                    if (f7X || f7Y) {
                        child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                        child.warmCompletedSweeps = 0;
                    }
                    if (!config_.lmTargeted) {
                        wantOpt.push_back(uni(s.sel) < config_.optFraction ? char{1} : char{0});
                    }
                    if (gmDefer && !gmPc.ev.empty()) {
                        gmPc.idx = offspring.size();
                        gmPc.parent = &a;
                        gmPend.push_back(std::move(gmPc));
                    }
                    offspring.push_back(std::move(child));
                    if (psHelper.running()) {
                        psHelper.publish(offspring.size());
                    } else if (psGen.active) {
                        psStreamAdd(offspring.back(), offspring.size() - 1, *s.activeDs, s.range, psGridsX, psGridsY, psGen);
                        psStreamPump(psGen);
                    }
                }
                }
                if (psHelper.running() && !psHelper.finish()) { psGen.active = false; }
                if (gmDefer && !gmPend.empty()) { gmDeferredPass(s, offspring, gmPend, &gmBandit, gmIr ? &gmIrStubs : nullptr); }
                if (gmIr || gmRefitRw) { gmRefitMask.assign(offspring.size(), 0); }
                std::chrono::steady_clock::time_point gmScoreT0;
                if (config_.gradientMutationMode == 3) { gmScoreT0 = std::chrono::steady_clock::now(); }
                if (measureTreeSizes) {
                    double const frac = tsStats.total > 0
                        ? static_cast<double>(tsStats.overMaxLen) / static_cast<double>(2 * tsStats.total) : 0.0;
                    std::cout << "[measure_tree_sizes] gen=" << gen
                              << " nCoordsSeen=" << (2 * tsStats.total)
                              << " overMaxLen(24)=" << tsStats.overMaxLen
                              << " frac=" << frac
                              << " maxSeen=" << tsStats.maxSeen
                              << " overFromMutation=" << tsStats.overFromMutation
                              << " overFromCrossover=" << tsStats.overFromCrossover << "\n";
                }
                if (config_.profileTiming) {
                    profBreedNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - profBreedT0).count(), std::memory_order_relaxed);
                }

                auto const curWorkers = workerCount();
                std::vector<std::vector<Operon::Scalar>> genGridsX, genGridsY;
                if (config_.useBatchedRegisterEval) {
                    evaluateBatchedPrepass(offspring, *s.activeDs, s.range, genGridsX, genGridsY);
                }
                CisPointSegBatch psBatch;
                CisPointSegBatch psBatch48;
                std::chrono::steady_clock::time_point profGpuPrepassT0;
                if (config_.profileTiming) { profGpuPrepassT0 = std::chrono::steady_clock::now(); }
                bool const psBatchOk = psGen.active ? psStreamFinish(offspring.size(), psGen, psBatch)
                    : (psStreamWanted() && psStreamBegin(offspring.size(), s.range, psGen))
                        ? [&] {
                              for (std::size_t i = 0; i < offspring.size(); ++i) {
                                  psStreamAdd(offspring[i], i, *s.activeDs, s.range, psGridsX, psGridsY, psGen);
                                  psStreamPump(psGen);
                              }
                              return psStreamFinish(offspring.size(), psGen, psBatch);
                          }()
                    : (config_.useBatchedGpuPointSeg && cisFitnessMode())
                        ? evaluateBatchedPointSegPrepass(offspring, *s.activeDs, s.range, psGridsX, psGridsY, psBatch, &psBatch48)
                        : false;
                if (config_.profileTiming) {
                    profGpuPrepassNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - profGpuPrepassT0).count(), std::memory_order_relaxed);
                }
                std::chrono::steady_clock::time_point profParallelWallT0;
                if (config_.profileTiming) { profParallelWallT0 = std::chrono::steady_clock::now(); }
                bool const batchedLmOk = config_.useBatchedLmInnerLoop && !config_.lmTargeted
                    && config_.constOptimizer == 1 && !config_.useLm && !config_.guardConstOpt
                    && config_.lmDownsampleK == 0 && config_.useBytecodeEvaluator;
                bool const t11Route = oneResGpuRoute();
                if (batchedLmOk) {
                    std::vector<PairIndividual*> wantPtrs;
                    for (std::size_t i = 0; i < offspring.size(); ++i) {
                        if (wantOpt[i] != 0) { wantPtrs.push_back(&offspring[i]); }
                    }
                    if (config_.profileTiming) {
                        auto t0 = std::chrono::steady_clock::now();
                        optConstsGnBatchedPopulation(wantPtrs, *s.activeDs, s.range, config_.polishFullObjective);
                        auto t1 = std::chrono::steady_clock::now();
                        profLmNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count(), std::memory_order_relaxed);
                    } else {
                        optConstsGnBatchedPopulation(wantPtrs, *s.activeDs, s.range, config_.polishFullObjective);
                    }
                    parallelApply(offspring.size(), curWorkers, [&](std::size_t i) {
                        if (wantOpt[i] != 0) {
                            if (t11Route) { return; }
                            evaluateGuarded(offspring[i], *s.activeDs, s.range);
                        } else if (psBatchOk && i < psGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, psGridsX[i], psGridsY[i], &psBatch, i, &psBatch48);
                        } else if (config_.useBatchedRegisterEval && i < genGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, genGridsX[i], genGridsY[i]);
                        } else {
                            evaluateGuarded(offspring[i], *s.activeDs, s.range);
                        }
                    });
                    if (t11Route) { rescoreBatched(wantPtrs, *s.activeDs, s.range, curWorkers); }
                } else if (!config_.lmTargeted) {
                    parallelApply(offspring.size(), curWorkers, [&](std::size_t i) {
                        if (wantOpt[i] != 0) {
                            if (t11Route) { optConstsInLoopDeferEval(offspring[i], *s.activeDs, s.range); }
                            else { optConstsInLoopGuardedFinite(offspring[i], *s.activeDs, s.range); }
                        } else if (psBatchOk && i < psGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, psGridsX[i], psGridsY[i], &psBatch, i, &psBatch48);
                        } else if (config_.useBatchedRegisterEval && i < genGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, genGridsX[i], genGridsY[i]);
                        } else {
                            evaluateGuarded(offspring[i], *s.activeDs, s.range);
                        }
                    });
                    if (t11Route) {
                        std::vector<PairIndividual*> wantPtrs;
                        for (std::size_t i = 0; i < offspring.size(); ++i) {
                            if (wantOpt[i] != 0) { wantPtrs.push_back(&offspring[i]); }
                        }
                        rescoreBatched(wantPtrs, *s.activeDs, s.range, curWorkers);
                    }
                } else {
                    bool const umaDiag = pointSegSchedDiagEnabled();
                    auto const fp0 = umaDiag ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
                    parallelApply(offspring.size(), curWorkers, [&](std::size_t i) {
                        if (psBatchOk && i < psGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, psGridsX[i], psGridsY[i], &psBatch, i, &psBatch48);
                        } else if (config_.useBatchedRegisterEval && i < genGridsX.size()) {
                            evaluateGuardedFromGrid(offspring[i], *s.activeDs, s.range, genGridsX[i], genGridsY[i]);
                        } else {
                            evaluateGuarded(offspring[i], *s.activeDs, s.range);
                        }
                    });
                    if (umaDiag) {
                        pointSegSchedDiag().firstPassNs.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - fp0).count(), std::memory_order_relaxed);
                    }
                    bool const sf2FitAll = sf2(17) || sf2(18);
                    std::size_t const k = sf2FitAll ? offspring.size() : sf2(16)
                        ? std::min(offspring.size(), sf2Admitted_ >= 0 ? static_cast<std::size_t>(sf2Admitted_) : offspring.size() / 2)
                        : std::min(
                        offspring.size(),
                        static_cast<std::size_t>(std::lround(config_.optFraction * static_cast<double>(offspring.size()))));
                    auto order = lmRefinementOrder(offspring, k, false);
                    std::vector<std::pair<std::size_t, std::size_t>> sf2Dups;
                    if (sf2(18)) { sf2DedupeOrder(offspring, order, sf2Dups); }
                    std::size_t kEff = std::min(k, order.size());
                    if (config_.useCascadeCull) {
                        Operon::Scalar const floor = populationFloor(s.pop);
                        ++nCascadeCullGens_;
                        nCascadeCullSeen_ += k;
                        kEff = 0;
                        while (kEff < k && offspring[order[kEff]].shapeErr <= floor) { ++kEff; }
                        nCascadeCullSurvived_ += kEff;
                    }
                    bool const t11Track = t11Route || oneResDiagEnabled();
                    std::vector<char> lmRan(t11Track ? kEff : 0, char{0});
                    if (!gmRefitMask.empty() && config_.useConcurrentCpuGpuLM && cclmStageOk(s)) {
                        for (std::size_t j = 0; j < kEff; ++j) { gmRefitMask[order[j]] = 1; }
                    }
                    if (config_.useConcurrentCpuGpuLM && cclmStageOk(s)) {
                        runConcurrentLmStage(offspring, order, kEff, s, curWorkers, static_cast<int>(gen));
                    } else {
                    parallelApply(kEff, curWorkers, [&](std::size_t j) {
                        std::size_t const i = order[j];
                        if (config_.useLiveCoarseSieve) {
                            nLiveCoarseSieveSeen_.fetch_add(1, std::memory_order_relaxed);
                            auto const bound = liveCoarseBound(offspring[i], *s.activeDs, s.range);
                            if (std::isfinite(bound)) {
                                auto const value = static_cast<double>(bound);
                                liveCoarseSieveBoundSum_.fetch_add(value, std::memory_order_relaxed);
                                nLiveCoarseSieveFinite_.fetch_add(1, std::memory_order_relaxed);
                                auto oldMin = liveCoarseSieveBoundMin_.load(std::memory_order_relaxed);
                                while (value < oldMin && !liveCoarseSieveBoundMin_.compare_exchange_weak(
                                    oldMin, value, std::memory_order_relaxed)) {}
                                auto oldMax = liveCoarseSieveBoundMax_.load(std::memory_order_relaxed);
                                while (value > oldMax && !liveCoarseSieveBoundMax_.compare_exchange_weak(
                                    oldMax, value, std::memory_order_relaxed)) {}
                            }
                            if (!std::isfinite(bound) ||
                                bound > config_.liveCoarseSieveThreshold) {
                                nLiveCoarseSieveCulled_.fetch_add(1, std::memory_order_relaxed);
                                return;
                            }
                        }
                        if (!gmRefitMask.empty()) { gmRefitMask[i] = 1; }
                        if (t11Route) { optConstsInLoopDeferEval(offspring[i], *s.activeDs, s.range); }
                        else { optConstsInLoopGuardedFinite(offspring[i], *s.activeDs, s.range); }
                        if (t11Track) { lmRan[j] = 1; }
                    });
                    if (t11Track) {
                        std::vector<PairIndividual*> refined;
                        for (std::size_t j = 0; j < kEff; ++j) {
                            if (lmRan[j] != 0) { refined.push_back(&offspring[order[j]]); }
                        }
                        if (oneResDiagEnabled()) {
                            std::vector<PairIndividual const*> src(refined.begin(), refined.end());
                            oneResDiagComparePtrs("lmRefined", src, *s.activeDs, s.range);
                        }
                        if (t11Route) { rescoreBatched(refined, *s.activeDs, s.range, curWorkers); }
                    }
                    }
                    for (auto const& [dup, rep] : sf2Dups) {
                        auto const lvX = offspring[dup].sf2LvX, lvY = offspring[dup].sf2LvY;
                        auto const pe = offspring[dup].sf2ParentErr;
                        offspring[dup] = offspring[rep];
                        offspring[dup].sf2LvX = lvX; offspring[dup].sf2LvY = lvY; offspring[dup].sf2ParentErr = pe;
                    }
                }
                if (!gmRefitMask.empty() && !config_.lmTargeted) {
                    for (std::size_t i = 0; i < offspring.size() && i < wantOpt.size(); ++i) { gmRefitMask[i] = wantOpt[i] != 0 ? 1 : 0; }
                }
                if (config_.profileTiming) {
                    profParallelWallNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - profParallelWallT0).count(), std::memory_order_relaxed);
                }
                if (config_.gradientMutationMode == 3 && !gmBandit.empty()) {
                    double const scoreMs = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - gmScoreT0).count();
                    gmBanditUpdate(s, offspring, gmBandit, scoreMs, curWorkers, gmRefitRw ? &gmRefitMask : nullptr);
                }
                if (gmIr && !gmIrStubs.empty()) { gmIrObserve(s, offspring, gmIrStubs, gmRefitMask, curWorkers); }

                if (config_.bracketNewtonRuntime != nullptr &&
                    config_.bracketNewtonGenerationBuilder && config_.bracketNewtonResultConsumer) {
                    BracketNewtonParentPacingMap parentMap;
                    bool const parentMapReady = !config_.bracketNewtonConfig.warmStart ||
                        (config_.bracketNewtonParentMapBuilder &&
                         (parentMap = config_.bracketNewtonParentMapBuilder(offspring, s.pop),
                          validParentPacingMap(parentMap, offspring.size(), s.pop.size())));
                    bool const useAsyncUma = config_.bracketNewtonUmaRunner != nullptr &&
                        config_.bracketNewtonAsyncSubmit && config_.bracketNewtonAsyncDrain;
                    if (!useAsyncUma) {
                        for (auto& candidate : offspring) {
                            auto generation = config_.bracketNewtonGenerationBuilder(candidate);
                            std::optional<std::size_t> parentIndex;
                            std::vector<double> parentPacing;
                            if (config_.bracketNewtonConfig.warmStart && parentMapReady) {
                                auto const index = &candidate - offspring.data();
                                parentIndex = parentMap.parentByOffspring[index];
                                parentPacing = parentMap.pacingStateByParent[*parentIndex];
                                generation.parentIndex = parentIndex;
                                generation.parentPacingState = parentPacing;
                                if (!parentMap.footpointStateByParent.empty())
                                    generation.footpoints = parentMap.footpointStateByParent[*parentIndex];
                            }
                            auto result = submitGenerationWithFallback(
                                *config_.bracketNewtonRuntime, config_.bracketNewtonConfig,
                                std::move(generation.residual), std::move(generation.gram),
                                std::move(generation.rhs), std::move(generation.parentPacingState),
                                generation.parentIndex, std::move(generation.footpoints));
                            config_.bracketNewtonResultConsumer(candidate, result);
                        }
                    } else {
                        struct PendingUma {
                            PairIndividual* candidate{nullptr};
                            std::shared_ptr<std::atomic<bool>> done;
                        };
                        std::array<std::shared_ptr<PendingUma>, UmaPingPongContract::slotCount> pending{};
                        auto drain = [&]() {
                            config_.bracketNewtonAsyncDrain();
                            for (std::size_t slot = 0; slot < pending.size(); ++slot) {
                                if (pending[slot] == nullptr || pending[slot]->done->load()) continue;
                                auto result = config_.bracketNewtonUmaRunner->complete(slot);
                                if (pending[slot]->candidate != nullptr)
                                    config_.bracketNewtonResultConsumer(*pending[slot]->candidate, result.transaction);
                                pending[slot].reset();
                            }
                        };
                        for (auto& candidate : offspring) {
                            std::size_t const slot = (&candidate - offspring.data()) % pending.size();
                            if (pending[slot] != nullptr) drain();
                            auto generation = config_.bracketNewtonGenerationBuilder(candidate);
                            if (config_.bracketNewtonConfig.warmStart && parentMapReady) {
                                auto const index = &candidate - offspring.data();
                                generation.parentIndex = parentMap.parentByOffspring[index];
                                generation.parentPacingState =
                                    parentMap.pacingStateByParent[*generation.parentIndex];
                                if (!parentMap.footpointStateByParent.empty())
                                    generation.footpoints =
                                        parentMap.footpointStateByParent[*generation.parentIndex];
                            }
                            if (!validGeneration(config_.bracketNewtonConfig, generation)) {
                                if (config_.bracketNewtonAsyncSkip) config_.bracketNewtonAsyncSkip();
                                continue;
                            }
                            auto state = std::make_shared<PendingUma>();
                            state->candidate = &candidate;
                            state->done = std::make_shared<std::atomic<bool>>(false);
                            bool const queued = config_.bracketNewtonUmaRunner->submitAsync(
                                slot, config_.bracketNewtonConfig, generation.residual,
                                generation.gram, generation.rhs, generation.parentPacingState,
                                generation.parentIndex, generation.footpoints);
                            if (!queued) {
                                BracketNewtonResult failed;
                                failed.usedFallback = true;
                                config_.bracketNewtonResultConsumer(candidate, failed);
                                continue;
                            }
                            pending[slot] = state;
                            bool const dispatched = config_.bracketNewtonAsyncSubmit(
                                config_.bracketNewtonConfig, generation,
                                [this, state, slot](BracketNewtonResult result) {
                                    auto completed = config_.bracketNewtonUmaRunner->completeFromBackend(
                                        slot, std::move(result));
                                    state->done->store(true);
                                    if (state->candidate != nullptr)
                                        config_.bracketNewtonResultConsumer(*state->candidate,
                                                                            completed.transaction);
                                });
                            if (!dispatched) drain();
                        }
                        drain();
                    }
                }

                if (config_.useEngineAmplitudeRuntime && config_.bracketNewtonRuntime != nullptr) {
                    for (auto& candidate : offspring) {
                        std::vector<PairAmplitudeNodeSlot> slots;
                        for (std::size_t node = 0; node < candidate.x.Nodes().size() && slots.size() < 1U; ++node) {
                            if (candidate.x.Nodes()[node].IsConstant()) {
                                slots.push_back({PairAmplitudeAxis::x, node,
                                                  BracketNewtonParameterKind::amplitude});
                            }
                        }
                        auto context = buildEngineAmplitudeContext(candidate, slots);
                        if (!context.has_value()) { continue; }
                        auto runtimeConfig = config_.bracketNewtonConfig;
                        runtimeConfig.population = 1U;
                        runtimeConfig.samples = context->sampleT.size();
                        runtimeConfig.parameters = slots.size();
                        auto const newtonSigma2 = runtimeConfig.useNoiseWeightedNewton
                            ? newtonPointSigma2ForRange(s.range) : std::vector<double>{};
                        auto generation = buildLinearAmplitudeGeneration(
                            *context, runtimeConfig,
                            newtonSigma2.empty() ? nullptr : &newtonSigma2);
                        if (!generation.has_value()) { continue; }
                        if (runtimeConfig.warmStart) {
                            generation->footpoints.segmentIndex = candidate.warmSegmentIndex;
                            generation->footpoints.segmentT = candidate.warmSegmentT;
                            generation->footpoints.completedSweeps = candidate.warmCompletedSweeps;
                            generation->parentIndex = 0U;
                            generation->parentPacingState = {
                                static_cast<double>(candidate.warmCompletedSweeps)};
                        }
                        if (runtimeConfig.useFp16Storage) {
                            generation->fp16Basis = context->linearBasisX;
                            generation->fp16Residual = context->targetX;
                        }
                        auto result = submitGenerationWithFallback(
                            *config_.bracketNewtonRuntime, runtimeConfig, std::move(*generation));
                        if (!result.ok || result.update.size() != slots.size()) { continue; }
                        if (runtimeConfig.warmStart && !result.nextFootpoints.segmentIndex.empty()) {
                            candidate.warmSegmentIndex = result.nextFootpoints.segmentIndex;
                            candidate.warmSegmentT = result.nextFootpoints.segmentT;
                            candidate.warmCompletedSweeps = result.nextFootpoints.completedSweeps;
                        }
                        PairIndividual original = candidate;
                        auto linearScore = [&](PairIndividual const& trial) {
                            std::vector<double> values;
                            for (auto const& slot : slots) {
                                values.push_back(static_cast<double>(trial.x.Nodes()[slot.nodeIndex].Value));
                            }
                            return linearAmplitudeSse(*context, values);
                        };
                        if (!applyMappedAmplitudeNodes(candidate, slots, result.update, linearScore)) {
                            continue;
                        }
                        auto const before = shapeOnlyGuarded(original, *s.activeDs, s.range);
                        evaluateGuarded(candidate, *s.activeDs, s.range);
                        if (!std::isfinite(candidate.shapeErr) || candidate.shapeErr > before) {
                            candidate = std::move(original);
                        }
                    }
                }

                if (sf2(7)) {
                    for (auto& o : offspring) {
                        bool const imp = std::isfinite(o.shapeErr) && o.shapeErr < o.sf2ParentErr;
                        for (int lv : {o.sf2LvX, o.sf2LvY}) {
                            if (lv < 0) { continue; }
                            auto const L = static_cast<std::size_t>(lv);
                            if (sf2LvTrials_.size() <= L) { sf2LvTrials_.resize(L + 1, 0.0); sf2LvSucc_.resize(L + 1, 0.0); }
                            sf2LvTrials_[L] += 1.0;
                            if (imp) { sf2LvSucc_[L] += 1.0; }
                        }
                        o.sf2LvX = o.sf2LvY = -1;
                    }
                }
                std::chrono::steady_clock::time_point profSelectT0;
                if (config_.profileTiming) { profSelectT0 = std::chrono::steady_clock::now(); }
                if (config_.useArchive) { archiveUpdate(s.archive, offspring, s.archiveCap); }
                if (config_.useComplexityHallOfFame) {
                    hallOfFameUpdate(s.hallOfFame, offspring, kHallOfFameK);
                }

                std::chrono::steady_clock::time_point combineT0;
                if (config_.profileTiming) { combineT0 = std::chrono::steady_clock::now(); }
                if (islandsActive) {
                    nextPop.clear();
                    for (auto const& rg : islRanges) {
                        std::size_t const lo = rg.first, hi = rg.second, sz = hi - lo;
                        islCombined.clear();
                        for (std::size_t j = lo; j < hi; ++j) { islCombined.push_back(std::move(s.pop[j])); }
                        for (std::size_t j = lo; j < hi; ++j) { islCombined.push_back(std::move(offspring[j])); }
                        selectNextInto(islCombined, sz, islOut);
                        for (auto& p : islOut) { nextPop.push_back(std::move(p)); }
                    }
                    if (config_.profileTiming) {
                        profSelectCombineNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - combineT0).count(), std::memory_order_relaxed);
                    }
                } else {
                    combined.clear();
                    for (auto& p : s.pop) { combined.push_back(std::move(p)); }
                    for (auto& p : offspring) { combined.push_back(std::move(p)); }
                    if (config_.profileTiming) {
                        profSelectCombineNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - combineT0).count(), std::memory_order_relaxed);
                        int64_t ndsNs = 0, crowdNs = 0;
                        selectNextInto(combined, N, nextPop, &ndsNs, &crowdNs);
                        profSelectNdsNs_.fetch_add(ndsNs, std::memory_order_relaxed);
                        profSelectCrowdNs_.fetch_add(crowdNs, std::memory_order_relaxed);
                    } else {
                        selectNextInto(combined, N, nextPop);
                    }
                }
                std::swap(s.pop, nextPop);
                if (sf2_ != 0) {
                    long adm = 0;
                    for (auto& p : s.pop) { if (p.sf2Child) { ++adm; } p.sf2Child = false; p.sf2LvX = p.sf2LvY = -1; }
                    sf2Admitted_ = adm;
                    sf2AdmSum_ += static_cast<double>(adm); ++sf2AdmN_;
                    sf2AdmMin_ = (sf2AdmMin_ < 0) ? adm : std::min(sf2AdmMin_, adm);
                    sf2AdmMax_ = std::max(sf2AdmMax_, adm);
                }

                if (islandsActive && config_.migrationInterval > 0 &&
                    ((gen + 1) % static_cast<std::size_t>(config_.migrationInterval) == 0)) {
                    migrateRing(s.pop, islRanges,
                                static_cast<std::size_t>(std::max(1, config_.migrationCount)));
                }

                if (config_.useArchive) {
                    archiveReinject(s.pop, s.archive, s.sel);
                    rankAndCrowd(s.pop);
                }
                if (config_.profileTiming) {
                    profSelectNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - profSelectT0).count(), std::memory_order_relaxed);
                }

                if (config_.floorStop > Operon::Scalar{0} && !noise_.empty()) {
                    if (s.inSimplify) {
                        ++s.simplifyGensRan;
                        if (--s.simplifyLeft == 0) { s.lastGen = static_cast<int>(gen); ++s.currentGen; break; }
                    } else {
                        if (floorReached(s.pop)) {
                            s.lastGen = static_cast<int>(gen);
                            if (config_.floorSimplifyGens == 0) { ++s.currentGen; break; }
                            enterSimplify(s.pop);
                        }
                    }
                }

                if (config_.verbose && (gen % 5 == 0 || gen == config_.generations - 1)) {
                    auto* best = bestAccuracy(s.pop);
                    std::cout << "  gen " << gen << ": best_err=" << best->shapeErr
                              << " nodes=" << best->totalLength()
                              << " front=" << countFront0(s.pop);
                    if (config_.useArchive) { std::cout << " archive=" << s.archive.size(); }
                    std::cout << "\n";
                }

                if (config_.profileTiming) {
                    profGenNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - profGenT0).count(), std::memory_order_relaxed);
                }
                if (config_.gradientMutationTrajectory) { gmTrajectoryPoint(s.pop, s); }
                if (overBudget()) { s.truncated = true; ++s.currentGen; break; }
            }
        } else {
            std::size_t const Lay = config_.alpsLayers;
            std::size_t const cap = std::max<std::size_t>(N / Lay, 8);
            auto ageLimit = [&](std::size_t i) -> int {
                return (i + 1 == Lay) ? std::numeric_limits<int>::max()
                                      : static_cast<int>(config_.ageGap * (i + 1));
            };
            auto tourPtr = [&](std::vector<PairIndividual*>& poolv) -> PairIndividual* {
                std::uniform_int_distribution<std::size_t> pick(0, poolv.size() - 1);
                auto* a = poolv[pick(s.sel)]; auto* b = poolv[pick(s.sel)];
                return crowdedLess(*a, *b) ? a : b;
            };

            for (; s.currentGen < targetGen; ++s.currentGen) {
                std::size_t gen = s.currentGen;
                s.lastGen = static_cast<int>(gen);

                if (s.lowRes && s.stagedMid && gen >= s.switchGen1) { switchToMid(s, nullptr, &s.layers, static_cast<int>(gen)); }
                else if (s.lowRes && gen >= s.switchGen2) { switchToFull(s, nullptr, &s.layers, static_cast<int>(gen)); }

                for (auto& lay : s.layers) { for (auto& ind : lay) { ++ind.age; } }

                while (s.layers[0].size() < cap) {
                    PairIndividual ind; ind.x = makeTree(); ind.y = makeTree(); ind.age = 0;
                    evaluateGuarded(ind, *s.activeDs, s.range); s.layers[0].push_back(std::move(ind));
                }

                std::vector<std::vector<PairIndividual>> off(Lay);
                for (std::size_t i = 0; i < Lay; ++i) {
                    std::vector<PairIndividual*> poolv;
                    for (auto& ind : s.layers[i]) { poolv.push_back(&ind); }
                    if (i > 0) { for (auto& ind : s.layers[i - 1]) { poolv.push_back(&ind); } }
                    if (poolv.empty()) { continue; }
                    std::vector<char> wantOpt;
                    wantOpt.reserve(cap);
                    for (std::size_t k = 0; k < cap; ++k) {
                        auto* a = tourPtr(poolv); auto* b = tourPtr(poolv);
                        PairIndividual child;
                        child.x = (uni(s.sel) < 0.9) ? (*s.crossover)(s.rng, a->x, b->x) : a->x;
                        child.y = (uni(s.sel) < 0.9) ? (*s.crossover)(s.rng, a->y, b->y) : a->y;
                        bool f7X = false, f7Y = false;
                        if (uni(s.sel) < 0.25) { child.x = mutateChildTreeCapped(s, child.x, a->x, *a, true, f7X, &child.y); }
                        if (uni(s.sel) < 0.25) { child.y = mutateChildTreeCapped(s, child.y, a->y, *a, false, f7Y, &child.x); }
                        if (child.x.Length() > 2 * config_.maxLength ||
                            child.y.Length() > 2 * config_.maxLength) { child = *a; f7X = f7Y = false; }
                        if (!structurallyValid(child)) { child = *a; ++nInvalidChildren_; f7X = f7Y = false; }
                        if (f7X || f7Y) {
                            child.warmSegmentIndex.clear(); child.warmSegmentT.clear();
                            child.warmCompletedSweeps = 0;
                        }
                        child.age = std::max(a->age, b->age);
                        wantOpt.push_back(uni(s.sel) < config_.optFraction ? char{1} : char{0});
                        off[i].push_back(std::move(child));
                    }
                    parallelApply(off[i].size(), nWorkers, [&](std::size_t k) {
                        if (wantOpt[k] != 0) {
                            optConstsInLoopGuardedFinite(off[i][k], *s.activeDs, s.range);
                        } else {
                            evaluateGuarded(off[i][k], *s.activeDs, s.range);
                        }
                    });
                }

                if (config_.useArchive) {
                    for (auto const& o : off) { archiveUpdate(s.archive, o, s.archiveCap); }
                }
                if (config_.useComplexityHallOfFame) {
                    for (auto const& o : off) { hallOfFameUpdate(s.hallOfFame, o, kHallOfFameK); }
                }

                std::vector<PairIndividual> allp;
                for (auto& lay : s.layers) { for (auto& ind : lay) { allp.push_back(std::move(ind)); } }
                for (auto& o : off) { for (auto& ind : o) { allp.push_back(std::move(ind)); } }
                for (auto& lay : s.layers) { lay.clear(); }
                for (auto& ind : allp) {
                    std::size_t tgt = Lay - 1;
                    for (std::size_t i = 0; i < Lay; ++i) {
                        if (ind.age <= ageLimit(i)) { tgt = i; break; }
                    }
                    s.layers[tgt].push_back(std::move(ind));
                }
                for (auto& lay : s.layers) {
                    if (lay.size() > cap) { lay = selectNext(lay, cap); }
                    else if (!lay.empty()) { rankAndCrowd(lay); }
                }

                if (config_.floorStop > Operon::Scalar{0} && !noise_.empty()) {
                    if (s.inSimplify) {
                        ++s.simplifyGensRan;
                        if (--s.simplifyLeft == 0) { s.lastGen = static_cast<int>(gen); ++s.currentGen; break; }
                    } else {
                        bool reached = false;
                        for (auto& lay : s.layers) {
                            if (lay.empty()) { continue; }
                            if (floorReached(lay)) { reached = true; break; }
                        }
                        if (reached) {
                            s.lastGen = static_cast<int>(gen);
                            if (config_.floorSimplifyGens == 0) { ++s.currentGen; break; }
                            s.inSimplify = true;
                            s.simplifyLeft = config_.floorSimplifyGens;
                            errorClampValue() = std::max(config_.floorClamp, config_.floorStop);
                            for (auto& lay : s.layers) { if (!lay.empty()) { rankAndCrowd(lay); } }
                        }
                    }
                }

                if (overBudget()) { s.truncated = true; ++s.currentGen; break; }
            }
        }
    }

