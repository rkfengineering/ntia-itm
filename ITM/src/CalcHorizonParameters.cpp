#include <ITM/ItmCommonCalculator.h>
#include <ITM/MathHelpers.h>

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>

namespace NTIA::ITM {
    void ItmCommonCalculator::setHorizonParameters(const double& effEarthRadius_km) {
        // Compute radials for Tx & Rx (ignore radius of earth since it cancels out in the later math)
        const auto& terrainHeightList_m = m_itmResults.m_intermResults.m_terrainProfile.m_terrainHeightList_m;
        double txRadial_m = terrainHeightList_m.front() + m_txHeight_m;
        double rxRadial_m = terrainHeightList_m.back() + m_rxHeight_m;

        // For ease of reference in the code
        const std::size_t& numPointsMinusTx = m_itmResults.m_intermResults.m_terrainProfile.m_numPointsMinusTx;
        const double& sampleResolution_m = m_itmResults.m_intermResults.m_terrainProfile.m_sampleResolution_m;
        const double& pathDist_km = m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km;
        double& finalTxHorizonAngle_rad = m_itmResults.m_intermResults.m_txHorizonAngle_rad;
        double& finalRxHorizonAngle_rad = m_itmResults.m_intermResults.m_rxHorizonAngle_rad;
        double& finalTxHorizonDist_km = m_itmResults.m_intermResults.m_txHorizonDist_m;
        double& finalRxHorizonDist_km = m_itmResults.m_intermResults.m_rxHorizonDist_m;

        // Set the terminal horizon angles as if the terminals are line-of-sight, [TN101, Eq 6.15]
        finalTxHorizonAngle_rad = (rxRadial_m - txRadial_m) / pathDist_km - pathDist_km / (2.0 * effEarthRadius_km);
        finalRxHorizonAngle_rad = -(rxRadial_m - txRadial_m) / pathDist_km - pathDist_km / (2.0 * effEarthRadius_km);

        finalTxHorizonDist_km = pathDist_km;
        finalRxHorizonDist_km = pathDist_km;

        // Initialize test tx & rx horizon distances
        double txDist_m = 0.0;
        double rxDist_m = pathDist_km;

        for (std::size_t pointInd = 1u; pointInd < numPointsMinusTx; pointInd++) {
            txDist_m += sampleResolution_m;
            rxDist_m -= sampleResolution_m;

            const double txHorizonAngle_deg = (terrainHeightList_m[pointInd] - txRadial_m) / txDist_m - txDist_m / (2.0 * effEarthRadius_km);
            const double rxHorizonAngle_deg = -(rxRadial_m - terrainHeightList_m[pointInd]) / rxDist_m - rxDist_m / (2.0 * effEarthRadius_km);

            // If better clearance to this point from Tx, shift its horizon
            if (txHorizonAngle_deg > finalTxHorizonAngle_rad) {
                finalTxHorizonAngle_rad = txHorizonAngle_deg;
                finalTxHorizonDist_km = txDist_m;
            }
            // If better clearance to this point from Rx, shift its horizon
            if (rxHorizonAngle_deg > finalRxHorizonAngle_rad) {
                finalRxHorizonAngle_rad = rxHorizonAngle_deg;
                finalRxHorizonDist_km = rxDist_m;
            }
        }
    }

    void ItmCommonCalculator::calcHorizonParameters() {
        double delta_h__meter;      // Terrain irregularity parameter
        double h_e__meter[2];       // Terminal effective heights

        const double effEarthRadius_km = 1.0 / m_effEarthCurvature_perM; // Effective earth radius

        setHorizonParameters(effEarthRadius_km);

        // For ease of reference in the code
        const auto& terrainHeightList_m = m_itmResults.m_intermResults.m_terrainProfile.m_terrainHeightList_m;
        const double& pathDist_km = m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km;
        double& txHorizonDist_m = m_itmResults.m_intermResults.m_txHorizonDist_m;
        double& rxHorizonDist_m = m_itmResults.m_intermResults.m_rxHorizonDist_m;

        // "In our own work we have sometimes said that consideration of terrain elevations should begin at a point about 15 times the tower height"
        //      - [Hufford, 1982] Page 25
        const double startDist_m = std::min({15.0 * m_txHeight_m, 0.1 * m_txHeight_m});                 // take lesser: 10% of horizon distance or 15x terminal height
        const double endDist_m = pathDist_km - std::min({15.0 * m_rxHeight_m, 0.1 * m_rxHeight_m});    // same as above, but measured from Rx side

        double& terrainIrreg_m = m_itmResults.m_intermResults.m_terrainIrreg_m;
        terrainIrreg_m = calcTerrainIrreg_m(startDist_m, endDist_m);

        if (txHorizonDist_m + rxHorizonDist_m > 1.5 * pathDist_km) {
            // The combined horizon distance is at least 50% larger than the total path distance
            //  -> so we are well within the line-of-sight range

            // Y1 = Tx LLS fit, Y2 = Rx LLS fit
            const auto fitResults = MathHelpers::fitTerrainProfile_linearLeastSquares(m_itmResults.m_intermResults.m_terrainProfile, 
                        startDist_m, endDist_m);

            // For ease of reference in the code
            double& txHorizonAngle_rad = m_itmResults.m_intermResults.m_txHorizonAngle_rad;
            double& rxHorizonAngle_rad = m_itmResults.m_intermResults.m_rxHorizonAngle_rad;
            double& txEffHorizDist_m = m_itmResults.m_intermResults.m_txEffHorizonDist_m;
            double& rxEffHorizDist_m = m_itmResults.m_intermResults.m_rxEffHorizonDist_m;
            double& txEffHeight_m = m_itmResults.m_intermResults.m_txEffHeight_m;
            double& rxEffHeight_m = m_itmResults.m_intermResults.m_rxEffHeight_m;

            txEffHorizDist_m = m_txHeight_m + std::abs(terrainHeightList_m.front() - fitResults.m_y1Value);
            rxEffHorizDist_m = m_rxHeight_m + std::abs(terrainHeightList_m.back() - fitResults.m_y2Value);

            // Recalculate horizon distances
            txHorizonDist_m = std::sqrt(2.0 * txEffHorizDist_m * effEarthRadius_km) * 
                        std::exp(-0.07 * std::sqrt(terrainIrreg_m / std::max({txEffHorizDist_m, 5.0})));
            rxHorizonDist_m = std::sqrt(2.0 * rxEffHorizDist_m * effEarthRadius_km) * 
                        std::exp(-0.07 * std::sqrt(terrainIrreg_m / std::max({rxEffHorizDist_m, 5.0})));

            const double combinedHorizonDist_m = txHorizonDist_m + rxHorizonDist_m;
            double effScalar;
            if (combinedHorizonDist_m <= pathDist_km) {
                effScalar = (pathDist_km / combinedHorizonDist_m) * (pathDist_km / combinedHorizonDist_m);

                txEffHeight_m *= effScalar;
                txEffHorizDist_m = std::sqrt(2.0 * txEffHeight_m * effEarthRadius_km) * std::exp(-0.07 * sqrt(terrainIrreg_m / std::max({txEffHeight_m, 5.0})));
                rxEffHeight_m *= effScalar;
                rxEffHorizDist_m = std::sqrt(2.0 * rxEffHeight_m * effEarthRadius_km) * std::exp(-0.07 * sqrt(terrainIrreg_m / std::max({rxEffHeight_m, 5.0})));
            }

            effScalar = sqrt(2.0 * txEffHeight_m * effEarthRadius_km);
            txHorizonAngle_rad = (0.65 * terrainIrreg_m * (effScalar / txEffHorizDist_m - 1.0) - 2.0 * txEffHeight_m) / effScalar;
            effScalar = sqrt(2.0 * rxEffHeight_m * effEarthRadius_km);
            rxHorizonAngle_rad = (0.65 * terrainIrreg_m * (effScalar / rxEffHorizDist_m - 1.0) - 2.0 * rxEffHeight_m) / effScalar;
        }
        else {
            const auto txFitResults = MathHelpers::fitTerrainProfile_linearLeastSquares(m_itmResults.m_intermResults.m_terrainProfile, 
                        startDist_m, 0.9 * txHorizonDist_m);
            m_itmResults.m_intermResults.m_txEffHeight_m = m_txHeight_m + std::abs(terrainHeightList_m.front() - txFitResults.m_y1Value);

            const auto rxFitResults = MathHelpers::fitTerrainProfile_linearLeastSquares(m_itmResults.m_intermResults.m_terrainProfile,
                        pathDist_km - 0.9 * rxHorizonDist_m, endDist_m);
            m_itmResults.m_intermResults.m_rxEffHeight_m = m_rxHeight_m + std::abs(terrainHeightList_m.back() - rxFitResults.m_y2Value);
        }
    }
} // end namespace