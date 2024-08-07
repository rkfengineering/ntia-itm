#include <ITM/ItmCommonCalculator.h>
#include <ITM/MathHelpers.h>

#include <algorithm>

/*=============================================================================
 |
 |  Description:  Compute the terrain irregularity parameter, delta_h
 |
 |        Input:  pfl[]          - Terrain data
 |                d_start__meter - Distance into the terrain profile to start 
 |                                 considering data, in meters
 |                d_end__meter   - Distance into the terrain profile to end 
 |                                 considering data, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  delta_h__meter - Terrain irregularity parameter, in meters
 |
 *===========================================================================*/

namespace NTIA::ITM {
    double ItmCommonCalculator::calcTerrainIrreg_m(const double& distToStart_m, const double& distToEnd_m) {
        const std::size_t& numPointsMinusTx = m_itmResults.m_intermResults.m_terrainProfile.m_numPointsMinusTx;
        const double& sampleResolution_m = m_itmResults.m_intermResults.m_terrainProfile.m_sampleResolution_m;

        double xStart = distToStart_m / sampleResolution_m;    // index to start considering terrain points
        double xEnd = distToEnd_m /sampleResolution_m;         // index to stop considering terrain points

        // Not enough data
        if (xEnd - xStart < 2.0) {
            return 0.0;
        }

        // TODO(vmartin): Work with Alex to figure out what this *intends* to do and maybe rewrite entirely?
        std::size_t tenPercentInd = static_cast<std::size_t>(0.1 * (xEnd - xStart + 8.0));
        tenPercentInd = std::min({std::max({4u, tenPercentInd}), 25u});

        std::size_t maxInd = 10u * tenPercentInd - 5u;
        std::size_t ninetyPercentInd = maxInd - tenPercentInd;

        TerrainProfile adjustedProfile;
        adjustedProfile.m_numPointsMinusTx = maxInd - 1u;
        adjustedProfile.m_sampleResolution_m = 1.0;

        xEnd = (xEnd - xStart) / static_cast<double>(adjustedProfile.m_numPointsMinusTx);
        std::size_t xInd = static_cast<std::size_t>(xStart);
        xStart -= static_cast<double>(xInd) + 1.0;

        const auto& terrainHeightList_m = m_itmResults.m_intermResults.m_terrainProfile.m_terrainHeightList_m;
        for (std::size_t profileInd = 0; profileInd < maxInd; profileInd++)
        {
            while (xStart > 0.0 && (xInd + 1u) < numPointsMinusTx) {
                xStart--;
                xInd++;
            }

            const double adjustedTerrainHeight_m = terrainHeightList_m[xInd + 1u] + (terrainHeightList_m[xInd + 1u] - terrainHeightList_m[xInd]) * xStart;
            adjustedProfile.m_terrainHeightList_m.push_back(adjustedTerrainHeight_m);

            xStart += xEnd;
        }

        auto fitResults = MathHelpers::fitTerrainProfile_linearLeastSquares(adjustedProfile, distToStart_m, distToEnd_m);

        fitResults.m_y2Value = (fitResults.m_y2Value - fitResults.m_y1Value) / static_cast<double>(numPointsMinusTx);

        std::vector<double> fittedDiffList;
        // Calculate the difference between fitted line and actual data
        for (std::size_t adjustedProfileInd = 0; adjustedProfileInd < maxInd; adjustedProfileInd++) {
            fittedDiffList.push_back(adjustedProfile.m_terrainHeightList_m[adjustedProfileInd] - fitResults.m_y1Value);

            fitResults.m_y1Value += fitResults.m_y2Value;
        }

        std::nth_element(fittedDiffList.begin(), fittedDiffList.begin() + tenPercentInd - 1u, fittedDiffList.end(), std::greater<double>());
        const double q10 = fittedDiffList[tenPercentInd - 1u];

        std::nth_element(fittedDiffList.begin(), fittedDiffList.begin() + ninetyPercentInd, fittedDiffList.end(), std::greater<double>());
        const double q90 = fittedDiffList[ninetyPercentInd];

        double terrainIrreg_m = q10 - q90;

        // [ERL 79-ITS 67, Eqn 3], inverted
        terrainIrreg_m /= (1.0 - 0.8 * exp(-(distToEnd_m - distToStart_m) / 50.0e3));

        return terrainIrreg_m;
    }
} // end namespace