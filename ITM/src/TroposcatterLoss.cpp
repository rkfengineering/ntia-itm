#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

/*=============================================================================
 |
 |  Description:  The attenuation function, F(th * d)
 |
 |        Input:  td             - angularDist_rad * distance
 |
 |      Outputs:  [None]
 |
 |      Returns:  F()            - in dB
 |
 *===========================================================================*/

namespace NTIA::ITM::ItmHelpers {
    namespace {
        // Constants from [Algorithm, 6.9]
        double constexpr aAttenArray[3] = { 133.4, 104.6, 71.8 };
        double constexpr bAttenArray[3] = { 0.332e-3, 0.212e-3, 0.157e-3 };
        double constexpr cAttenArray[3] = { -10.0, -2.5, 5.0 };
    }

    double calcTropoAttenFunction_dB(const double& inputDist_m) {
        std::size_t arrayInd;

        // select the set of values to use
        if (inputDist_m <= 10.0e3)          // <= 10 km
            arrayInd = 0;
        else if (inputDist_m <= 70.0e3)     // 10 km to 70 km 
            arrayInd = 1;
        else                                // > 70 km
            arrayInd = 2;

        const double attenFunc_dB = aAttenArray[arrayInd] + bAttenArray[arrayInd] * inputDist_m 
                    + cAttenArray[arrayInd] * std::log10(inputDist_m); // [Algorithm, 6.9]

        return attenFunc_dB;
    }
}

/*=============================================================================
 |
 |  Description:  Troposcatter loss
 |
 |        Input:  tropoPathLength_m - Path distance, in meters
 |                theta_hzn[2]      - Terminal horizon angles
 |                d_hzn__meter[2]   - Terminal horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                earthEffRadius_m        - Effective earth radius, in meters
 |                N_s               - Surface refractivity, in N-Units
 |                f__mhz            - Frequency, in MHz
 |                theta_los         - Angular distance of LOS region
 |
 |      Outputs:  initialH0_dB                - H0_m() value
 |
 |      Returns:  F()               - in dB
 |
 *===========================================================================*/

namespace NTIA::ITM {
    double ItmCommonCalculator::calcTroposcatterLoss_dB(const double& tropoPathLength_m, const double& earthEffRadius_m, 
                const double& angularDist_LoS_rad, double& initialH0_dB) {
        double finalH0_dB = initialH0_dB;

        // Calculate angular wavelength
        const double waveLength_m = ItmHelpers::kSpeedOfLight_mPerS * 1.0e-6 / m_freq_MHz;
        const double waveNumber_radPerM = 2.0 * M_PI / waveLength_m;

        // If initialH0_dB is already set to a value > 15, no need to perform these calculations
        if (initialH0_dB <= 15.0) {
            const double& txHorizonDist_m = m_itmResults.m_intermResults.m_txHorizonDist_m;
            const double& rxHorizonDist_m = m_itmResults.m_intermResults.m_rxHorizonDist_m;
            const double& txEffHeight_m = m_itmResults.m_intermResults.m_txEffHeight_m;
            const double& rxEffHeight_m = m_itmResults.m_intermResults.m_rxEffHeight_m;
            const double& txHorizAngle_rad = m_itmResults.m_intermResults.m_txHorizonAngle_rad;
            const double& rxHorizAngle_rad = m_itmResults.m_intermResults.m_rxHorizonAngle_rad;

            double horizonDistDelta_m = txHorizonDist_m - rxHorizonDist_m;
            double effHeightRatio = rxEffHeight_m / txEffHeight_m;

            if (horizonDistDelta_m < 0.0)       // ensure correct frame of reference
            {
                horizonDistDelta_m = -horizonDistDelta_m;
                effHeightRatio = 1.0 / effHeightRatio;
            }

            const double angularDist_rad = txHorizAngle_rad + rxHorizAngle_rad + tropoPathLength_m / earthEffRadius_m;    // angular distance, in radians

            // [TN101, Eqn 9.4a]
            double r1_radSqrd = 2.0 * waveNumber_radPerM * angularDist_rad * txEffHeight_m;
            double r2_radSqrd = 2.0 * waveNumber_radPerM * angularDist_rad * rxEffHeight_m;

            // "If both r_1 and r_2 are less than 0.2 the function A_scat is not defined (or is infinite)" [Algorithm, page 11]
            if (r1_radSqrd < 0.2 && r2_radSqrd < 0.2) {
                return kDefaultMaxLoss_dB;
            }

            double asymmetryParam = (tropoPathLength_m - horizonDistDelta_m) / (tropoPathLength_m + horizonDistDelta_m);       // asymmetry parameter

            // "In all of this, we truncate the values of asymmetryParam and q at 0.1 and 10" [Algorithm, page 16]
            double q = std::min({std::max({0.1, effHeightRatio / asymmetryParam}), 10.0});      // TN101, Eqn 9.5
            asymmetryParam = std::max({0.1, asymmetryParam});                                   // TN101, Eqn 9.5

            double h_0__meter = (tropoPathLength_m - horizonDistDelta_m) * (tropoPathLength_m + horizonDistDelta_m) * angularDist_rad * 0.25 / tropoPathLength_m;   // height of cross-over, [Algorithm, 4.66] [TN101v1, 9.3b]

            double Z_0__meter = 1.7556e3;       // Scale height, [Algorithm, 4.67]
            double Z_1__meter = 8.0e3;          // [Algorithm, 4.67]
            double scatterEffTerm = (h_0__meter / Z_0__meter) * (1.0 + (0.031 - N_s * 2.32e-3 + N_s * N_s * 5.67e-6) * exp(-pow(std::min({1.7, h_0__meter / Z_1__meter}), 6)));     // Scattering efficiency factor, scatterEffTerm [TN101 Eqn 9.3a]

            const double tropoGain_r1 = ItmHelpers::calcTropoFreqGain_dB(r1_radSqrd, scatterEffTerm);
            const double tropoGain_r2 = ItmHelpers::calcTropoFreqGain_dB(r2_radSqrd, scatterEffTerm);
            const double avgTropoGain_dB = 0.5 * (tropoGain_r1 + tropoGain_r2);                        // First term in TN101v1, Eqn 9.5
            const double deltaH_minTerm = 6.0 * (0.6 - std::log10(std::max({scatterEffTerm, 1.0}))) * std::log10(asymmetryParam) * log10(q);
            const double deltaH_dB = std::min({avgTropoGain_dB, deltaH_minTerm});

            finalH0_dB = avgTropoGain_dB + deltaH_dB;       // TN101, Eqn 9.5
            finalH0_dB = std::max({finalH0_dB, 0.0});       // "If Delta_H_0 would make finalH0_dB negative, use finalH0_dB = 0" [TN101v1, p9.4] 

            if (scatterEffTerm < 1.0) { 
                const double sqrtTwo = std::sqrt(2.0);
                const double logTerm_sqrTerm = (1.0 + sqrtTwo / r1_radSqrd) * (1.0 + sqrtTwo / r2_radSqrd);
                const double logTerm_scalar = (r1_radSqrd + r2_radSqrd) / (r1_radSqrd + r2_radSqrd + 2.0 * sqrtTwo);
                const double finalH0_logTerm = std::log10(logTerm_sqrTerm * logTerm_sqrTerm * logTerm_scalar);
                finalH0_dB = scatterEffTerm * finalH0_dB + (1.0 - scatterEffTerm) * 10.0 * finalH0_logTerm;
            } // if <=1, interpolate with the special case of scatterEffTerm = 0

            // TODO(vmartin): Conditions here seem to be at odds with the if statement containing it, which requires (h0 <= 15)
            // "If, at d_5, calculations show that finalH0_dB will exceed 15 dB, they are replaced by the value it has at d_6" [Algorithm, page 12]
            if (finalH0_dB > 15.0 && initialH0_dB >= 0.0) {
                finalH0_dB = initialH0_dB;
            }

            initialH0_dB = finalH0_dB;
        }

        const double thConst = tropoPathLength_m / earthEffRadius_m - angularDist_LoS_rad;

        const double kD0_m = 40.0e3;   // [Algorithm, 6.8]
        const double logTerm = waveNumber_radPerM * ItmHelpers::kWaveToMHzFreqTerm * thConst * thConst * thConst * thConst;
        return ItmHelpers::calcTropoAttenFunction_dB(thConst * tropoPathLength_m) + 
                    10.0 * log10(logTerm) - 
                    0.1 * (m_surfaceRefractivity_N - 301.0) * exp(-thConst * tropoPathLength_m / kD0_m) +
                    finalH0_dB;    // [Algorithm, 4.63]
    }
}