#pragma once


#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Sketch {

struct SymPoint {
    double x{0.0};
    double y{0.0};
};

struct RotationCandidate {
    int order{0};
    double residual{0.0};
};

struct ReflectionAxis {
    double angle{0.0};
    double residual{0.0};
};

struct SymmetryReport {
    double cx{0.0};
    double cy{0.0};
    double meanSpacing{1.0};

    int rotationalOrder{0};
    double rotationalResidual{std::numeric_limits<double>::infinity()};
    std::vector<RotationCandidate> rotationCandidates;

    std::vector<ReflectionAxis> reflectionAxes;

};

namespace detail {

inline double meanSpacing(const std::vector<SymPoint>& pts) {
    if (pts.size() < 2) return 1.0;
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
        double dx = pts[i + 1].x - pts[i].x;
        double dy = pts[i + 1].y - pts[i].y;
        total += std::sqrt(dx * dx + dy * dy);
    }
    double avg = total / static_cast<double>(pts.size() - 1);
    return (avg > 1e-12) ? avg : 1e-12;
}

inline void centroid(const std::vector<SymPoint>& pts, double& cx, double& cy) {
    cx = 0.0;
    cy = 0.0;
    if (pts.empty()) return;
    for (const auto& p : pts) {
        cx += p.x;
        cy += p.y;
    }
    cx /= static_cast<double>(pts.size());
    cy /= static_cast<double>(pts.size());
}

inline double symmetrizedResidual(const std::vector<SymPoint>& a, const std::vector<SymPoint>& b) {
    if (a.empty() || b.empty()) return std::numeric_limits<double>::infinity();
    auto oneWay = [](const std::vector<SymPoint>& from, const std::vector<SymPoint>& to) {
        double sum = 0.0;
        for (const auto& p : from) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& q : to) {
                double dx = p.x - q.x;
                double dy = p.y - q.y;
                double d2 = dx * dx + dy * dy;
                if (d2 < best) best = d2;
            }
            sum += std::sqrt(best);
        }
        return sum / static_cast<double>(from.size());
    };
    return 0.5 * (oneWay(a, b) + oneWay(b, a));
}

inline std::vector<SymPoint> rotateAbout(const std::vector<SymPoint>& pts, double cx, double cy, double theta) {
    std::vector<SymPoint> out;
    out.reserve(pts.size());
    double c = std::cos(theta), s = std::sin(theta);
    for (const auto& p : pts) {
        double dx = p.x - cx, dy = p.y - cy;
        out.push_back({cx + dx * c - dy * s, cy + dx * s + dy * c});
    }
    return out;
}

inline std::vector<SymPoint> reflectAbout(const std::vector<SymPoint>& pts, double cx, double cy, double angle) {
    std::vector<SymPoint> out;
    out.reserve(pts.size());
    double c = std::cos(2.0 * angle), s = std::sin(2.0 * angle);
    for (const auto& p : pts) {
        double dx = p.x - cx, dy = p.y - cy;
        out.push_back({cx + dx * c + dy * s, cy + dx * s - dy * c});
    }
    return out;
}

}

struct SymmetryDetectOptions {
    int maxOrder{8};
    double tolGaps{0.75};
    int axisSteps{180};
    int axisPolishSteps{20};
    double axisDedupRad{0.15};
};

inline SymmetryReport detectSymmetry(const std::vector<SymPoint>& pointsIn,
                                      const SymmetryDetectOptions& opt = {}) {
    SymmetryReport report;
    if (pointsIn.size() < 4) {
        return report;
    }

    std::vector<SymPoint> pts = pointsIn;
    detail::centroid(pts, report.cx, report.cy);
    report.meanSpacing = detail::meanSpacing(pts);
    const double tol = opt.tolGaps * report.meanSpacing;

    report.rotationCandidates.reserve(std::max(0, opt.maxOrder - 1));
    for (int m = 2; m <= opt.maxOrder; ++m) {
        double theta = 2.0 * M_PI / static_cast<double>(m);
        auto rotated = detail::rotateAbout(pts, report.cx, report.cy, theta);
        double res = detail::symmetrizedResidual(pts, rotated);
        report.rotationCandidates.push_back({m, res});
        if (res < tol && res < report.rotationalResidual) {
            report.rotationalResidual = res;
            report.rotationalOrder = m;
        }
    }
    report.rotationalOrder = 0;
    report.rotationalResidual = std::numeric_limits<double>::infinity();
    for (const auto& c : report.rotationCandidates) {
        if (c.residual < tol && c.order > report.rotationalOrder) {
            report.rotationalOrder = c.order;
            report.rotationalResidual = c.residual;
        }
    }

    std::vector<ReflectionAxis> coarse;
    coarse.reserve(opt.axisSteps);
    for (int i = 0; i < opt.axisSteps; ++i) {
        double angle = M_PI * static_cast<double>(i) / static_cast<double>(opt.axisSteps);
        auto reflected = detail::reflectAbout(pts, report.cx, report.cy, angle);
        double res = detail::symmetrizedResidual(pts, reflected);
        coarse.push_back({angle, res});
    }
    std::vector<ReflectionAxis> minima;
    int n = static_cast<int>(coarse.size());
    for (int i = 0; i < n; ++i) {
        double res = coarse[i].residual;
        double prev = coarse[(i - 1 + n) % n].residual;
        double next = coarse[(i + 1) % n].residual;
        if (res <= prev && res <= next && res < tol) {
            minima.push_back(coarse[i]);
        }
    }
    double step0 = M_PI / static_cast<double>(opt.axisSteps);
    for (auto& cand : minima) {
        double bestAngle = cand.angle;
        double bestRes = cand.residual;
        double halfWidth = step0;
        for (int iter = 0; iter < 4; ++iter) {
            double lo = bestAngle - halfWidth, hi = bestAngle + halfWidth;
            for (int k = 0; k <= opt.axisPolishSteps; ++k) {
                double a = lo + (hi - lo) * static_cast<double>(k) / static_cast<double>(opt.axisPolishSteps);
                auto reflected = detail::reflectAbout(pts, report.cx, report.cy, a);
                double res = detail::symmetrizedResidual(pts, reflected);
                if (res < bestRes) {
                    bestRes = res;
                    bestAngle = a;
                }
            }
            halfWidth *= 0.4;
        }
        cand.angle = bestAngle;
        cand.residual = bestRes;
    }
    for (auto& cand : minima) {
        cand.angle = std::fmod(cand.angle, M_PI);
        if (cand.angle < 0) cand.angle += M_PI;
    }
    std::sort(minima.begin(), minima.end(), [](const ReflectionAxis& a, const ReflectionAxis& b) {
        return a.angle < b.angle;
    });
    for (const auto& cand : minima) {
        bool merged = false;
        for (auto& kept : report.reflectionAxes) {
            double d = std::fabs(cand.angle - kept.angle);
            d = std::min(d, M_PI - d);
            if (d < opt.axisDedupRad) {
                if (cand.residual < kept.residual) kept = cand;
                merged = true;
                break;
            }
        }
        if (!merged) report.reflectionAxes.push_back(cand);
    }
    std::sort(report.reflectionAxes.begin(), report.reflectionAxes.end(),
              [](const ReflectionAxis& a, const ReflectionAxis& b) { return a.angle < b.angle; });

    return report;
}

}
