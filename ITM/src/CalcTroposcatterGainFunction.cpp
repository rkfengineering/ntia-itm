#include "ITM/itm.h"

namespace ITM::NTIA::ItmHelpers {
    namespace {
        // values from [Algorithm, 6.13]
        double constexpr aList[] = { 25.0, 80.0, 177.0, 395.0, 705.0 };
        double constexpr bList[] = { 24.0, 45.0, 68.0, 80.0, 105.0 };
    }

    double calcTropoFreqGainCurveFit_dB(const std::size_t arrayInd, const double &rTerm) {
        const double inv_rTerm = 1.0 / rTerm;
        const double inv_rTermSqrd = inv_rTerm * inv_rTerm;
        return 10.0 * std::log10(1.0 + aList[arrayInd] * inv_rTermSqrd * inv_rTermSqrd + 
                        bList[arrayInd] * inv_rTermSqrd);   // related to TN101v2, Eqn III.49, but from [Algorithm, 6.13]
    }

    double calcTropoFreqGain_dB(const double& rParam, double& scatterEfficiency) {
        // Force scatterEfficiency term to fall in between 1 <= eta_s <= 5
        scatterEfficiency = std::min({std::max({scatterEfficiency, 1.0}), 5.0});

        const std::size_t scatterInd = static_cast<std::size_t>(scatterEfficiency);
        const double scatterEffRemainder = scatterEfficiency - static_cast<double>(scatterInd);

        const double tropoGain_dB = calcTropoFreqGainCurveFit_dB(scatterInd - 1u, rParam);
        
        // If the scatter efficiency term is not an exact integer, interpolate
        if (scatterEffRemainder != 0.0) {
            return (1.0 - scatterEffRemainder) * tropoGain_dB + 
                        scatterEffRemainder * calcTropoFreqGainCurveFit_dB(scatterInd, rParam);
        }

        return tropoGain_dB;
    }
}
