#pragma once


#include <cmath>
#include <limits>

namespace Sketch {

inline auto chowAbstainThresholdBits(float abstentionCost) -> float
{
    double const ca = static_cast<double>(abstentionCost);
    if (!(ca > 0.0)) { return std::numeric_limits<float>::infinity(); }
    if (ca >= 1.0) { return -std::numeric_limits<float>::infinity(); }
    double const tau = std::log2((1.0 - ca) / ca);
    return static_cast<float>(tau);
}

inline auto chowShouldAnswer(float leadBits, float abstentionCost) -> bool
{
    float const tau = chowAbstainThresholdBits(abstentionCost);
    double lead = static_cast<double>(leadBits);
    if (!(lead > 0.0) || !std::isfinite(lead)) { lead = 0.0; }
    return lead > static_cast<double>(tau);
}

inline auto chowLeadBitsFromNegLog2(float bestNegLog2Score, float secondNegLog2Score) -> float
{
    double const best = static_cast<double>(bestNegLog2Score);
    double const second = static_cast<double>(secondNegLog2Score);
    if (!std::isfinite(best) || !std::isfinite(second)) { return 0.0f; }
    double const lead = second - best;
    return static_cast<float>(lead > 0.0 ? lead : 0.0);
}

inline auto chowShouldAnswerFromNegLog2(float bestNegLog2Score, float secondNegLog2Score, float abstentionCost) -> bool
{
    return chowShouldAnswer(chowLeadBitsFromNegLog2(bestNegLog2Score, secondNegLog2Score), abstentionCost);
}

}
