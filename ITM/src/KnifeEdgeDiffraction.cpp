#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

/*=============================================================================
 |
 |  Description:  Compute the knife-edge diffraction loss
 |
 |        Input:  d__meter          - Distance of interest, in meters
 |                f__mhz            - Frequency, in MHz
 |                a_e__meter        - Effective earth radius, in meters
 |                theta_los         - Angular distance of line-of-sight region
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_k__db        - Knife-edge diffraction loss, in dB
 |
 *===========================================================================*/

namespace NTIA::ITM {
    double ItmCommonCalculator::calcKnifeEdgeDiffractLoss_dB(const double& inputDist_m, const double& effEarthRadius_km, const double& angularDist_LoS_rad) {
        const double& txHorizonDist_m = m_itmResults.m_intermResults.m_txHorizonDist_m;
        const double& rxHorizonDist_m = m_itmResults.m_intermResults.m_rxHorizonDist_m;

        const double maxDist_LoS_m = txHorizonDist_m + rxHorizonDist_m;                         // Maximum line-of-sight distance for actual path
        const double angularDist_nLoS_rad = inputDist_m / effEarthRadius_km - angularDist_LoS_rad;    // Angular distance of diffraction region [Algorithm, Eqn 4.12]

        const double diffractDist_nLoS_m = inputDist_m - maxDist_LoS_m;                                  // Diffraction distance, in meters

        // 1 / (4 pi) = 0.0795775
        // [TN101, Eqn I.7]
        const double angularDistSqrd = angularDist_nLoS_rad * angularDist_nLoS_rad;
        const double nuCommonTerm = 0.0795775 * (m_freq_MHz / ItmHelpers::kWaveToMHzFreqTerm) * angularDistSqrd * diffractDist_nLoS_m;
        const double nu1 = nuCommonTerm * txHorizonDist_m / (diffractDist_nLoS_m + txHorizonDist_m);
        const double nu2 = nuCommonTerm * rxHorizonDist_m / (diffractDist_nLoS_m + rxHorizonDist_m);

        double A_k__db = ItmHelpers::calcFresnelIntegral(nu1) + ItmHelpers::calcFresnelIntegral(nu2);                   // [TN101, Eqn I.1]

        return A_k__db;
    }
}