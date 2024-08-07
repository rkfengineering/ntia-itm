#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

/*=============================================================================
 |
 |  Description:  Compute the diffraction loss at a specified distance
 |
 |        Input:  d__meter          - Path distance, in meters
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                Z_g            - Complex ground impedance
 |                a_e__meter     - Effective earth radius, in meters
 |                delta_h__meter - Terrain irregularity parameter, in meters
 |                h__meter[2]       - Terminal heights, in meters
 |                mode           - Area or Point-to-Point mode flag
 |                theta_los      - Angular distance of line-of-sight region
 |                d_sML__meter   - Maximum line-of-sight distance for 
 |                                 a smooth earth, in meters
 |                f__mhz         - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_d__db        - Diffraction loss, in dB
 |
 *===========================================================================*/

namespace NTIA::ITM {
    double ItmCommonCalculator::calcDiffractLoss_dB(const double& inputDist_m, const double& effEarthRadius_m, const bool isP2P, 
                const double& angularDist_LoS_rad, const double& maxDistSmoothEarth_LoS_m) {
        const double attenKnifeEdge_dB = calcKnifeEdgeDiffractLoss_dB(inputDist_m, effEarthRadius_m, angularDist_LoS_rad);

        const double attenSmoothEarth_dB = calcSmoothEarthDiffractLoss_dB(inputDist_m, effEarthRadius_m, angularDist_LoS_rad);

        // Terrain roughness term, using d_sML__meter, per [ERL 79-ITS 67, page 3-13]
        const double tempTerrainIrreg_m = ItmHelpers::calcTerrainRoughness_m(maxDistSmoothEarth_LoS_m, m_itmResults.m_intermResults.m_terrainIrreg_m);

        const double sigmaH_m = ItmHelpers::calcSigmaH_m(tempTerrainIrreg_m);

        // Clutter factor
        // [ERL 79-ITS 67, Eqn 3.38c]
        const double attenClutterFactor_dB = std::min({15.0, 5.0 * std::log10(1.0 + 1.0e-5 * m_txHeight_m * m_rxHeight_m * m_freq_MHz * sigmaH_m)});

        // compute the weighting factor in the following calculations
        double temp2_terrainIrreg_m = ItmHelpers::calcTerrainRoughness_m(inputDist_m, m_itmResults.m_intermResults.m_terrainIrreg_m);

        double q = m_txHeight_m * m_rxHeight_m;
        const double qSubK = m_itmResults.m_intermResults.m_txEffHeight_m * m_itmResults.m_intermResults.m_rxEffHeight_m - q;

        // For low antennas with known path parameters, C ~= 10 [ERL 79-ITS 67, page 3-8]
        if (isP2P) {
            q += 10.0;
        }
        const double term1 = sqrt(1.0 + qSubK / q);                              // square root term in [ERL 79-ITS 67, Eqn 2.23]

        const double maxDist_LoS_m = m_itmResults.m_intermResults.m_txHorizonDist_m + 
                        m_itmResults.m_intermResults.m_rxHorizonDist_m;         // Maximum line-of-sight distance for actual path
        q = (term1 + (-angularDist_LoS_rad * effEarthRadius_m + maxDist_LoS_m) / inputDist_m) * 
                    std::min({temp2_terrainIrreg_m * m_freq_MHz / ItmHelpers::kWaveToMHzFreqTerm, 6283.2});

        // weighting factor [ERL 17-ITS 67, Eqn 3.23]
        double weightFactor = 25.1 / (25.1 + sqrt(q));

        return weightFactor * attenSmoothEarth_dB + (1.0 - weightFactor) * attenKnifeEdge_dB + attenClutterFactor_dB;
    }
}