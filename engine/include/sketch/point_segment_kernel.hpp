#pragma once


#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <operon/operon/core/types.hpp>

#if defined(__APPLE__) && defined(SKETCH_POINTSEG_ENABLE_VECTOR_EXP)
#include <Accelerate/Accelerate.h>
#define SKETCH_POINTSEG_HAVE_ACCELERATE 1
#endif

namespace Sketch {

struct PointSegmentConfig {
    bool usePointSegment{false};

    bool useVectorizedExp{false};

    bool useBandedSearch{false};
    std::size_t bandWindow{20};

    bool useSpatialHash{false};

    bool computeTau{false};

    bool computeNormal{false};
};

struct Vec2 {
    Operon::Scalar x{0};
    Operon::Scalar y{0};
};

struct PointSegmentResult {
    std::vector<Operon::Scalar> residual;
    std::vector<Operon::Scalar> tau;
    std::vector<Vec2> normalAtClosest;
    Operon::Scalar fitness{0};
    bool ok{true};
};

namespace detail {

inline auto psSafeSqrt(Operon::Scalar v) -> Operon::Scalar
{
    return std::sqrt(std::max(v, Operon::Scalar{0}));
}

}

inline auto evaluatePointSegment(
    std::vector<Vec2> const& strokePoints,
    std::vector<Vec2> const& candidatePolyline,
    PointSegmentConfig const& cfg = {}) -> PointSegmentResult
{
    PointSegmentResult out;
    auto const numPoints = strokePoints.size();
    auto const numVerts = candidatePolyline.size();

    if (numPoints == 0 || numVerts < 2) {
        out.ok = false;
        out.fitness = std::numeric_limits<Operon::Scalar>::max();
        return out;
    }

    auto const numSegs = numVerts - 1;

    constexpr auto kEps = std::numeric_limits<Operon::Scalar>::epsilon();
    Operon::Scalar const lenSqEps = kEps * Operon::Scalar{100};

    thread_local std::vector<Vec2> segV;
    thread_local std::vector<Operon::Scalar> segLenSq;
    thread_local std::vector<Vec2> segNormal;
    segV.resize(numSegs);
    segLenSq.resize(numSegs);
    segNormal.resize(numSegs);
    Operon::Scalar sumSegLen{0};
    std::size_t finiteSegCount{0};

    for (std::size_t j = 0; j < numSegs; ++j) {
        Vec2 const& a = candidatePolyline[j];
        Vec2 const& b = candidatePolyline[j + 1];
        Vec2 v{b.x - a.x, b.y - a.y};
        Operon::Scalar lenSq = v.x * v.x + v.y * v.y;
        if (!std::isfinite(lenSq)) {
            lenSq = 0;
            out.ok = false;
        }
        segV[j] = v;
        segLenSq[j] = lenSq;

        if (lenSq > lenSqEps) {
            Operon::Scalar const invLen = Operon::Scalar{1} / detail::psSafeSqrt(lenSq);
            Vec2 const tangent{v.x * invLen, v.y * invLen};
            segNormal[j] = Vec2{-tangent.y, tangent.x};
            sumSegLen += detail::psSafeSqrt(lenSq);
            ++finiteSegCount;
        } else {
            segNormal[j] = Vec2{0, 1};
        }
    }

    Operon::Scalar const meanSegLen = (finiteSegCount > 0)
        ? (sumSegLen / static_cast<Operon::Scalar>(finiteSegCount))
        : Operon::Scalar{0};
    Operon::Scalar const softMinScale = std::max(meanSegLen, static_cast<Operon::Scalar>(kEps));

    Operon::Scalar const gatherC = std::log(Operon::Scalar{1} / kEps) + Operon::Scalar{1};

    auto evalSegAtPoint = [&](std::size_t j, Vec2 const& p) {
        Vec2 const& a = candidatePolyline[j];
        Vec2 const& v = segV[j];
        Operon::Scalar const lenSq = segLenSq[j];
        Vec2 const pa{p.x - a.x, p.y - a.y};
        Operon::Scalar const paSq = pa.x * pa.x + pa.y * pa.y;
        Operon::Scalar const paDotV = pa.x * v.x + pa.y * v.y;
        Operon::Scalar tC;
        if (lenSq > lenSqEps) {
            Operon::Scalar const tRaw = paDotV / lenSq;
            tC = std::min(std::max(tRaw, Operon::Scalar{0}), Operon::Scalar{1});
        } else {
            tC = 0;
        }
        Operon::Scalar d2 = paSq - Operon::Scalar{2} * tC * paDotV + tC * tC * lenSq;
        bool segOk = true;
        if (!std::isfinite(d2) || d2 < 0) {
            d2 = std::isfinite(d2) ? Operon::Scalar{0} : std::numeric_limits<Operon::Scalar>::max();
            if (!std::isfinite(paSq) || !std::isfinite(paDotV)) { segOk = false; }
        }
        Operon::Scalar const dVal = detail::psSafeSqrt(d2);
        Vec2 const diff{pa.x - tC * v.x, pa.y - tC * v.y};
        Operon::Scalar signedVal = diff.x * segNormal[j].x + diff.y * segNormal[j].y;
        if (!std::isfinite(signedVal)) { signedVal = 0; segOk = false; }
        struct SegEval { Operon::Scalar d2, dVal, signedVal, tC; bool ok; };
        return SegEval{d2, dVal, signedVal, tC, segOk};
    };

    thread_local std::vector<std::vector<std::size_t>> psGridCells;
    thread_local std::vector<std::size_t> psTouchedCells;
    thread_local std::vector<std::size_t> psSegStamp;
    thread_local std::size_t psStampCounter{0};
    thread_local std::vector<std::size_t> psGatherIdx;
    thread_local std::vector<Operon::Scalar> psGatherD;
    thread_local std::vector<Operon::Scalar> psGatherSigned;
    thread_local std::vector<Operon::Scalar> psGatherTc;
    thread_local std::vector<Operon::Scalar> psGatherExp;

    Operon::Scalar gridOriginX{0}, gridOriginY{0}, gridCellSize{1};
    std::size_t gridNumCellsX{1}, gridNumCellsY{1};

    if (cfg.useSpatialHash) {
        Operon::Scalar minX = candidatePolyline[0].x, maxX = candidatePolyline[0].x;
        Operon::Scalar minY = candidatePolyline[0].y, maxY = candidatePolyline[0].y;
        for (auto const& v : candidatePolyline) {
            minX = std::min(minX, v.x); maxX = std::max(maxX, v.x);
            minY = std::min(minY, v.y); maxY = std::max(maxY, v.y);
        }
        gridOriginX = minX;
        gridOriginY = minY;
        gridCellSize = softMinScale;
        constexpr std::size_t kMaxGridCells = 1u << 20;
        auto computeDims = [&]() {
            Operon::Scalar const w = std::max(maxX - minX, Operon::Scalar{0});
            Operon::Scalar const h = std::max(maxY - minY, Operon::Scalar{0});
            gridNumCellsX = static_cast<std::size_t>(std::ceil(static_cast<double>(w) / static_cast<double>(gridCellSize))) + 1;
            gridNumCellsY = static_cast<std::size_t>(std::ceil(static_cast<double>(h) / static_cast<double>(gridCellSize))) + 1;
        };
        computeDims();
        for (int guard = 0; guard < 64 && gridNumCellsX * gridNumCellsY > kMaxGridCells; ++guard) {
            gridCellSize *= Operon::Scalar{2};
            computeDims();
        }

        for (auto idx : psTouchedCells) {
            if (idx < psGridCells.size()) { psGridCells[idx].clear(); }
        }
        psTouchedCells.clear();
        psGridCells.resize(gridNumCellsX * gridNumCellsY);
        psSegStamp.assign(numSegs, 0);

        auto cellCoord = [&](Operon::Scalar x, Operon::Scalar y) {
            auto cx = static_cast<std::ptrdiff_t>(std::floor(static_cast<double>((x - gridOriginX) / gridCellSize)));
            auto cy = static_cast<std::ptrdiff_t>(std::floor(static_cast<double>((y - gridOriginY) / gridCellSize)));
            return std::pair<std::ptrdiff_t, std::ptrdiff_t>{cx, cy};
        };
        auto clampDim = [](std::ptrdiff_t v, std::size_t dim) {
            std::ptrdiff_t const hi = static_cast<std::ptrdiff_t>(dim) - 1;
            if (v < 0) { return std::ptrdiff_t{0}; }
            if (v > hi) { return hi; }
            return v;
        };

        for (std::size_t j = 0; j < numSegs; ++j) {
            Vec2 const& a = candidatePolyline[j];
            Vec2 const& b = candidatePolyline[j + 1];
            Operon::Scalar const segMinX = std::min(a.x, b.x), segMaxX = std::max(a.x, b.x);
            Operon::Scalar const segMinY = std::min(a.y, b.y), segMaxY = std::max(a.y, b.y);
            auto [cxLo0, cyLo0] = cellCoord(segMinX, segMinY);
            auto [cxHi0, cyHi0] = cellCoord(segMaxX, segMaxY);
            std::ptrdiff_t const cxLo = clampDim(cxLo0, gridNumCellsX);
            std::ptrdiff_t const cxHi = clampDim(cxHi0, gridNumCellsX);
            std::ptrdiff_t const cyLo = clampDim(cyLo0, gridNumCellsY);
            std::ptrdiff_t const cyHi = clampDim(cyHi0, gridNumCellsY);
            for (std::ptrdiff_t cy = cyLo; cy <= cyHi; ++cy) {
                for (std::ptrdiff_t cx = cxLo; cx <= cxHi; ++cx) {
                    std::size_t const cellIdx = static_cast<std::size_t>(cy) * gridNumCellsX + static_cast<std::size_t>(cx);
                    psGridCells[cellIdx].push_back(j);
                    psTouchedCells.push_back(cellIdx);
                }
            }
        }
    }

    out.residual.assign(numPoints, Operon::Scalar{0});
    if (cfg.computeTau) { out.tau.assign(numPoints, Operon::Scalar{0}); } else { out.tau.clear(); }
    if (cfg.computeNormal) { out.normalAtClosest.assign(numPoints, Vec2{0, 0}); } else { out.normalAtClosest.clear(); }
    out.fitness = 0;

    thread_local std::vector<Operon::Scalar> dJ;
    thread_local std::vector<Operon::Scalar> signedJ;
    thread_local std::vector<Operon::Scalar> expJ;
    thread_local std::vector<Operon::Scalar> argJ;
    thread_local std::vector<Operon::Scalar> tcJ;
    dJ.resize(numSegs);
    signedJ.resize(numSegs);
    expJ.resize(numSegs);
    argJ.resize(numSegs);
    if (cfg.computeTau) { tcJ.resize(numSegs); }

    for (std::size_t i = 0; i < numPoints; ++i) {
        Vec2 const& p = strokePoints[i];
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) {
            out.ok = false;
            out.residual[i] = 0;
            continue;
        }

        Operon::Scalar hardMinD2 = std::numeric_limits<Operon::Scalar>::max();

        if (cfg.useSpatialHash) {
            auto [pcx0, pcy0] = [&]() {
                auto cx = static_cast<std::ptrdiff_t>(std::floor(static_cast<double>((p.x - gridOriginX) / gridCellSize)));
                auto cy = static_cast<std::ptrdiff_t>(std::floor(static_cast<double>((p.y - gridOriginY) / gridCellSize)));
                return std::pair<std::ptrdiff_t, std::ptrdiff_t>{cx, cy};
            }();
            std::size_t argMinIdx = 0;
            bool anyFound = false;
            std::size_t const stampA = ++psStampCounter;
            std::size_t const maxRing = gridNumCellsX + gridNumCellsY + 1;
            for (std::size_t ring = 0; ring <= maxRing; ++ring) {
                std::ptrdiff_t const r = static_cast<std::ptrdiff_t>(ring);
                for (std::ptrdiff_t cy = pcy0 - r; cy <= pcy0 + r; ++cy) {
                    if (cy < 0 || cy >= static_cast<std::ptrdiff_t>(gridNumCellsY)) { continue; }
                    bool const yEdge = (cy == pcy0 - r || cy == pcy0 + r);
                    for (std::ptrdiff_t cx = pcx0 - r; cx <= pcx0 + r; ++cx) {
                        if (cx < 0 || cx >= static_cast<std::ptrdiff_t>(gridNumCellsX)) { continue; }
                        if (!yEdge && cx != pcx0 - r && cx != pcx0 + r) { continue; }
                        std::size_t const cellIdx = static_cast<std::size_t>(cy) * gridNumCellsX + static_cast<std::size_t>(cx);
                        for (auto segIdx : psGridCells[cellIdx]) {
                            if (psSegStamp[segIdx] == stampA) { continue; }
                            psSegStamp[segIdx] = stampA;
                            auto ev = evalSegAtPoint(segIdx, p);
                            if (!ev.ok) { out.ok = false; }
                            if (ev.d2 < hardMinD2) { hardMinD2 = ev.d2; argMinIdx = segIdx; anyFound = true; }
                        }
                    }
                }
                if (anyFound) {
                    Operon::Scalar const bound = static_cast<Operon::Scalar>(ring) * gridCellSize;
                    if (bound * bound >= hardMinD2) { break; }
                }
            }
            if (!anyFound) {
                out.ok = false;
                hardMinD2 = std::numeric_limits<Operon::Scalar>::max();
            }

            if (std::isfinite(hardMinD2)) {
                out.fitness += hardMinD2;
            } else {
                out.ok = false;
            }

            Operon::Scalar const m = detail::psSafeSqrt(hardMinD2);
            Operon::Scalar const rad = m + gatherC * softMinScale;
            auto const cellRadius = static_cast<std::ptrdiff_t>(std::ceil(static_cast<double>(rad) / static_cast<double>(gridCellSize))) + 1;

            psGatherIdx.clear();
            std::size_t const stampB = ++psStampCounter;
            if (anyFound) { psSegStamp[argMinIdx] = stampB; psGatherIdx.push_back(argMinIdx); }
            std::ptrdiff_t const cyLo = std::max<std::ptrdiff_t>(0, pcy0 - cellRadius);
            std::ptrdiff_t const cyHi = std::min<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(gridNumCellsY) - 1, pcy0 + cellRadius);
            std::ptrdiff_t const cxLo = std::max<std::ptrdiff_t>(0, pcx0 - cellRadius);
            std::ptrdiff_t const cxHi = std::min<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(gridNumCellsX) - 1, pcx0 + cellRadius);
            for (std::ptrdiff_t cy = cyLo; cy <= cyHi; ++cy) {
                for (std::ptrdiff_t cx = cxLo; cx <= cxHi; ++cx) {
                    std::size_t const cellIdx = static_cast<std::size_t>(cy) * gridNumCellsX + static_cast<std::size_t>(cx);
                    for (auto segIdx : psGridCells[cellIdx]) {
                        if (psSegStamp[segIdx] == stampB) { continue; }
                        psSegStamp[segIdx] = stampB;
                        psGatherIdx.push_back(segIdx);
                    }
                }
            }

            std::size_t const kCount = psGatherIdx.size();
            psGatherD.resize(kCount);
            psGatherSigned.resize(kCount);
            psGatherTc.resize(kCount);
            psGatherExp.resize(kCount);
            for (std::size_t k = 0; k < kCount; ++k) {
                auto ev = evalSegAtPoint(psGatherIdx[k], p);
                if (!ev.ok) { out.ok = false; }
                psGatherD[k] = ev.dVal;
                psGatherSigned[k] = ev.signedVal;
                psGatherTc[k] = ev.tC;
            }
            Operon::Scalar sumExp{0};
            for (std::size_t k = 0; k < kCount; ++k) {
                Operon::Scalar const e = std::exp(-(psGatherD[k] - m) / softMinScale);
                psGatherExp[k] = e;
                sumExp += e;
            }
            Operon::Scalar blended{0};
            for (std::size_t k = 0; k < kCount; ++k) {
                blended += (psGatherExp[k] / sumExp) * psGatherSigned[k];
            }
            if (!std::isfinite(blended)) { blended = 0; out.ok = false; }
            out.residual[i] = blended;

            if (cfg.computeTau) {
                Operon::Scalar tauBlend{0};
                for (std::size_t k = 0; k < kCount; ++k) {
                    Operon::Scalar const paramJ = static_cast<Operon::Scalar>(psGatherIdx[k]) + psGatherTc[k];
                    tauBlend += (psGatherExp[k] / sumExp) * paramJ;
                }
                if (!std::isfinite(tauBlend)) { tauBlend = 0; out.ok = false; }
                out.tau[i] = tauBlend;
            }

            if (cfg.computeNormal) {
                Vec2 normalBlend{0, 0};
                for (std::size_t k = 0; k < kCount; ++k) {
                    Operon::Scalar const w = psGatherExp[k] / sumExp;
                    Vec2 const& segN = segNormal[psGatherIdx[k]];
                    normalBlend.x += w * segN.x;
                    normalBlend.y += w * segN.y;
                }
                if (!std::isfinite(normalBlend.x) || !std::isfinite(normalBlend.y)) {
                    normalBlend = Vec2{0, 0};
                    out.ok = false;
                }
                out.normalAtClosest[i] = normalBlend;
            }

            continue;
        }

        std::size_t jLo = 0, jHi = numSegs;
        if (cfg.useBandedSearch && numPoints > 1) {
            double const frac = static_cast<double>(i) / static_cast<double>(numPoints - 1);
            auto const center = static_cast<std::ptrdiff_t>(frac * static_cast<double>(numSegs - 1) + 0.5);
            auto const win = static_cast<std::ptrdiff_t>(cfg.bandWindow);
            std::ptrdiff_t lo = center - win;
            std::ptrdiff_t hi = center + win + 1;
            if (lo < 0) { lo = 0; }
            if (hi > static_cast<std::ptrdiff_t>(numSegs)) { hi = static_cast<std::ptrdiff_t>(numSegs); }
            jLo = static_cast<std::size_t>(lo);
            jHi = static_cast<std::size_t>(hi);
        }

        for (std::size_t j = jLo; j < jHi; ++j) {
            Vec2 const& a = candidatePolyline[j];
            Vec2 const& v = segV[j];
            Operon::Scalar const lenSq = segLenSq[j];

            Vec2 const pa{p.x - a.x, p.y - a.y};
            Operon::Scalar const paSq = pa.x * pa.x + pa.y * pa.y;
            Operon::Scalar const paDotV = pa.x * v.x + pa.y * v.y;

            Operon::Scalar tC;
            if (lenSq > lenSqEps) {
                Operon::Scalar const tRaw = paDotV / lenSq;
                tC = std::min(std::max(tRaw, Operon::Scalar{0}), Operon::Scalar{1});
            } else {
                tC = 0;
            }

            Operon::Scalar d2 = paSq - Operon::Scalar{2} * tC * paDotV + tC * tC * lenSq;
            if (!std::isfinite(d2) || d2 < 0) {
                d2 = std::isfinite(d2) ? Operon::Scalar{0} : std::numeric_limits<Operon::Scalar>::max();
                if (!std::isfinite(paSq) || !std::isfinite(paDotV)) { out.ok = false; }
            }

            Operon::Scalar const dVal = detail::psSafeSqrt(d2);

            Vec2 const diff{pa.x - tC * v.x, pa.y - tC * v.y};
            Operon::Scalar signedVal = diff.x * segNormal[j].x + diff.y * segNormal[j].y;
            if (!std::isfinite(signedVal)) { signedVal = 0; out.ok = false; }

            dJ[j] = dVal;
            signedJ[j] = signedVal;
            if (cfg.computeTau) { tcJ[j] = tC; }
            if (d2 < hardMinD2) { hardMinD2 = d2; }
        }

        if (std::isfinite(hardMinD2)) {
            out.fitness += hardMinD2;
        } else {
            out.ok = false;
        }

        Operon::Scalar const m = detail::psSafeSqrt(hardMinD2);
        Operon::Scalar sumExp{0};
        bool tookVectorizedPath = false;
#if defined(SKETCH_POINTSEG_HAVE_ACCELERATE)
        if constexpr (std::is_same_v<Operon::Scalar, float>) {
            if (cfg.useVectorizedExp) {
                std::size_t const winLen = jHi - jLo;
                for (std::size_t j = jLo; j < jHi; ++j) {
                    argJ[j] = -(dJ[j] - m) / softMinScale;
                }
                int const n = static_cast<int>(winLen);
                vvexpf(expJ.data() + jLo, argJ.data() + jLo, &n);
                for (std::size_t j = jLo; j < jHi; ++j) {
                    sumExp += expJ[j];
                }
                tookVectorizedPath = true;
            }
        }
#endif
        if (!tookVectorizedPath) {
            for (std::size_t j = jLo; j < jHi; ++j) {
                Operon::Scalar const e = std::exp(-(dJ[j] - m) / softMinScale);
                expJ[j] = e;
                sumExp += e;
            }
        }
        Operon::Scalar blended{0};
        for (std::size_t j = jLo; j < jHi; ++j) {
            blended += (expJ[j] / sumExp) * signedJ[j];
        }
        if (!std::isfinite(blended)) { blended = 0; out.ok = false; }
        out.residual[i] = blended;

        if (cfg.computeTau) {
            Operon::Scalar tauBlend{0};
            for (std::size_t j = jLo; j < jHi; ++j) {
                Operon::Scalar const paramJ = static_cast<Operon::Scalar>(j) + tcJ[j];
                tauBlend += (expJ[j] / sumExp) * paramJ;
            }
            if (!std::isfinite(tauBlend)) { tauBlend = 0; out.ok = false; }
            out.tau[i] = tauBlend;
        }

        if (cfg.computeNormal) {
            Vec2 normalBlend{0, 0};
            for (std::size_t j = jLo; j < jHi; ++j) {
                Operon::Scalar const w = expJ[j] / sumExp;
                normalBlend.x += w * segNormal[j].x;
                normalBlend.y += w * segNormal[j].y;
            }
            if (!std::isfinite(normalBlend.x) || !std::isfinite(normalBlend.y)) {
                normalBlend = Vec2{0, 0};
                out.ok = false;
            }
            out.normalAtClosest[i] = normalBlend;
        }
    }

    if (!std::isfinite(out.fitness)) {
        out.fitness = std::numeric_limits<Operon::Scalar>::max();
        out.ok = false;
    }

    return out;
}

inline auto monotonicityEpsJitter(std::vector<Operon::Scalar> const& tau) -> Operon::Scalar
{
    auto const n = tau.size();
    if (n < 2) { return Operon::Scalar{0}; }
    thread_local std::vector<Operon::Scalar> steps;
    steps.resize(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        steps[i] = std::abs(tau[i + 1] - tau[i]);
    }
    thread_local std::vector<Operon::Scalar> sorted;
    sorted.assign(steps.begin(), steps.end());
    auto const mid = sorted.begin() + static_cast<std::ptrdiff_t>(sorted.size() / 2);
    std::nth_element(sorted.begin(), mid, sorted.end());
    Operon::Scalar const medianStep = sorted.empty() ? Operon::Scalar{0} : *mid;
    return Operon::Scalar{2} * medianStep;
}

inline auto monotonicityPenalty(std::vector<Operon::Scalar> const& tau, Operon::Scalar epsJitter) -> Operon::Scalar
{
    auto const n = tau.size();
    if (n < 2) { return Operon::Scalar{0}; }
    Operon::Scalar acc{0};
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Operon::Scalar const backward = (tau[i] - tau[i + 1]) - epsJitter;
        if (backward > Operon::Scalar{0}) { acc += backward * backward; }
    }
    return acc / static_cast<Operon::Scalar>(n - 1);
}

}
