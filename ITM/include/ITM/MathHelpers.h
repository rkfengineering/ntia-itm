#ifndef ITM_MATH_HELPERS_H
#define ITM_MATH_HELPERS_H

#define _USE_MATH_DEFINES
#include <cmath>

namespace NTIA::ITM::MathHelpers {
    namespace {
        double constexpr kC0 {2.515516};
        double constexpr kC1 {0.802853};
        double constexpr kC2 {0.010328};
        double constexpr kD1 {1.432788};
        double constexpr kD2 {0.189269};
        double constexpr kD3 {0.001308};
    }

    struct TerrainFitResults {
        double m_y1Value;
        double m_y2Value;
    };
    
    /// @brief This function computes the inverse complementary cumulative distribution function approximation as
    /// described in Formula 26.2.23 in Abramowitz & Stegun. This approximation has an error of abs(epsilon(p)) < 4.5e-4
    /// @param q Quantile fraction (0.0 < q < 1.0)
    /// @return Inverse complementary cumulative distribution function, Q(q)^-1
    double calcInvComplCumulDistribFunc(const double& q) {
        const double xVal = (q > 0.5) ? 1.0 - q : q;

        const double T_x = std::sqrt(-2.0 * std::log(xVal));

        const double zetaNumer = (kC2 * T_x + kC1) * T_x + kC0;
        const double zetaDenom = ((kD3 * T_x + kD2) * T_x + kD1) * T_x + 1.0;
        const double zeta_x = zetaNumer / zetaDenom;

        const double Q_q = T_x - zeta_x;
        return (q > 0.5) ? -Q_q : Q_q;
    }
    /*=============================================================================
    |
    |  Description:  Perform a linear least squares fit to the terrain data
    |
    |        Input:  pfl[2]         - Input data array, in pfl format
    |                d_start        - Start distance
    |                d_end          - End distance
    |
    |      Outputs:  fit_y1         - Fitted y1 value
    |                fit_y2         - Fitted y2 value
    |
    |      Returns:  [None]
    |
    *===========================================================================*/
    TerrainFitResults fitTerrainProfile_linearLeastSquares(const TerrainProfile& terrainProfile, 
                const double& distToStart_m, const double& distToEnd_m) {
        // For ease of reference in the code
        const std::size_t& numPointsMinusTx = terrainProfile.m_numPointsMinusTx;
        const double& sampleResolution_m = terrainProfile.m_sampleResolution_m;
        const auto& terrainHeightList_m = terrainProfile.m_terrainHeightList_m;

        int startInd = std::abs(distToStart_m / sampleResolution_m);
        int endInd = numPointsMinusTx - std::abs(numPointsMinusTx - distToEnd_m / sampleResolution_m);

        // TODO(vmartin): Figure out with Alex what this means? Why is it done this way?
        if (endInd <= startInd) {
            startInd = std::abs(startInd - 1);
            const int endIndDiff = static_cast<int>(numPointsMinusTx) - (endInd + 1);
            endInd = numPointsMinusTx - std::abs(endIndDiff);
        }

        const std::size_t xLength = endInd - startInd;

        double middleShiftedInd_double = -0.5 * static_cast<double>(xLength);
        double middleShiftedEndInd_double = endInd + middleShiftedInd_double;

        double sumOfY = 0.5 * (terrainHeightList_m[startInd] + terrainHeightList_m[endInd]);
        double scaledSumOfY = 0.5 * (terrainHeightList_m[startInd] - terrainHeightList_m[endInd]) * middleShiftedInd_double;

        for (std::size_t profileInd = startInd; profileInd <= endInd; profileInd++) {
            startInd++;
            middleShiftedInd_double++;

            sumOfY += terrainHeightList_m[startInd];
            scaledSumOfY += terrainHeightList_m[startInd] * middleShiftedInd_double;
        }

        sumOfY /= xLength;
        const double sumOfYScale = 12.0 / ((xLength * xLength + 2.0) * xLength);
        scaledSumOfY *= sumOfYScale;

        TerrainFitResults results;
        results.m_y1Value = sumOfY - scaledSumOfY * middleShiftedEndInd_double;
        results.m_y2Value = sumOfY + scaledSumOfY * (static_cast<double>(numPointsMinusTx) - middleShiftedEndInd_double);

        return results;
    }
} // end namespace MathHelpers

#endif // ITM_MATH_HELPERS_H
