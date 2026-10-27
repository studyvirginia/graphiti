#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace Sketch {

struct CandidatePopulationRecord {
    std::string dataset;
    std::string candidateId;
    std::size_t rank{0};
    bool selected{false};
    double fullScore{0.0};
    double coarseBound{0.0};
    std::vector<double> curveX;
    std::vector<double> curveY;
};

inline auto writeCandidatePopulationJsonl(std::string const& path,
                                          std::vector<CandidatePopulationRecord> const& rows)
    -> bool {
    std::ofstream out(path);
    if (!out) { return false; }
    out << std::setprecision(17);
    for (auto const& row : rows) {
        if (row.curveX.size() != row.curveY.size() || row.curveX.empty() ||
            !std::isfinite(row.fullScore) || !std::isfinite(row.coarseBound)) { return false; }
        for (double value : row.curveX) { if (!std::isfinite(value)) { return false; } }
        for (double value : row.curveY) { if (!std::isfinite(value)) { return false; } }
        out << "{\"schema\":\"graphiti.candidate.v1\",\"dataset\":\""
            << row.dataset << "\",\"candidate_id\":\"" << row.candidateId
            << "\",\"rank\":" << row.rank << ",\"selected\":"
            << (row.selected ? "true" : "false") << ",\"full_score\":" << row.fullScore
            << ",\"coarse_bound\":" << row.coarseBound << ",\"curve_x\":[";
        for (std::size_t i = 0; i < row.curveX.size(); ++i) {
            if (i) { out << ','; }
            out << row.curveX[i];
        }
        out << "],\"curve_y\":[";
        for (std::size_t i = 0; i < row.curveY.size(); ++i) {
            if (i) { out << ','; }
            out << row.curveY[i];
        }
        out << "]}\n";
    }
    return static_cast<bool>(out);
}

}
