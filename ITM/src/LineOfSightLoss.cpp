#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

/*=============================================================================
 |
 |  Description:  Compute the loss in the line-of-sight region
 |
 |        Input:  inputDist_m          - Path distance, in meters
 |                h_e__meter[2]     - Terminal effective heights, in meters
 |                m_groundImpedance               - Complex surface transfer impedance
 |                delta_h__meter    - Terrain irregularity parameter
 |                diffractSlope               - Diffraction slope
 |                diffractLineIntercept              - Diffraction intercept
 |                d_sML__meter      - Maximum line-of-sight distance for
 |                                    a smooth earth, in meters
 |                m_freq_MHz            - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_los__db         - Loss, in dB
 |
 *===========================================================================*/

namespace NTIA::ITM {
    double ItmCommonCalculator::calcLineOfSightLoss_dB(const double& inputDist_m, 
                const double& diffractSlope, const double& diffractLineIntercept, const double& maxDistSmoothEarth_LoS_m) {
        const double tempTerrainIrreg_m = ItmHelpers::calcTerrainRoughness_m(inputDist_m, m_itmResults.m_intermResults.m_terrainIrreg_m);
        const double tempSigmaH = ItmHelpers::calcSigmaH_m(tempTerrainIrreg_m);

        // Angular wavenumber, k
        const double waveNumber = m_freq_MHz / ItmHelpers::kWaveToMHzFreqTerm;

        // [Algorithm, Eqn 4.46]
        const double& txEffHeight_m = m_itmResults.m_intermResults.m_txEffHeight_m;
        const double& rxEffHeight_m = m_itmResults.m_intermResults.m_rxEffHeight_m;

        const double effHeightSum_m = txEffHeight_m + rxEffHeight_m;
        const double sinOfPsi = effHeightSum_m / std::sqrt(inputDist_m * inputDist_m + effHeightSum_m * effHeightSum_m);

        // [Algorithm, Eqn 4.47]
        std::complex<double> reflCoeff_e = (sinOfPsi - m_groundImpedance) / (sinOfPsi + m_groundImpedance) * 
                    std::exp(-std::min({10.0, waveNumber * tempSigmaH * sinOfPsi}));

        // |R_e| = Magnitude of R_e', [Algorithm, Eqn 4.48]
        const double reflCoeff_mag = reflCoeff_e.real() * reflCoeff_e.real() + reflCoeff_e.imag() * reflCoeff_e.imag();
        if (reflCoeff_mag < 0.25 || reflCoeff_mag < sinOfPsi) {
            reflCoeff_e *= std::sqrt(sinOfPsi / reflCoeff_mag);
        }
        // phase difference between rays, [Algorithm, Eqn 4.49]
        double rayPhaseDiff_rad = waveNumber * 2.0 * txEffHeight_m * rxEffHeight_m / inputDist_m;

        // [Algorithm, Eqn 4.50]
        if (rayPhaseDiff_rad > M_PI / 2.0)
            rayPhaseDiff_rad = M_PI - (M_PI / 2.0) * (M_PI / 2.0) / rayPhaseDiff_rad;

        // Two-ray attenuation
        std::complex<double> twoRayReflCoeff = std::complex<double>(std::cos(rayPhaseDiff_rad), -std::sin(rayPhaseDiff_rad)) + reflCoeff_e;
        const double attenTwoRay_dB = -10.0 * std::log10(twoRayReflCoeff.real() * twoRayReflCoeff.real() + twoRayReflCoeff.imag() * twoRayReflCoeff.imag());

        // Extended diffraction attenuation
        const double diffractLoss_dB = diffractSlope * inputDist_m + diffractLineIntercept;

        // weighting factor
        const double w = 1.0 / (1.0 + m_freq_MHz * m_itmResults.m_intermResults.m_terrainIrreg_m / std::max({10.0e3, maxDistSmoothEarth_LoS_m}));

        return w * attenTwoRay_dB + (1.0 - w) * diffractLoss_dB;
    }
}