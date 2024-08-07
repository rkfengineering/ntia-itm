#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

namespace NTIA::ITM {
    namespace {
        double constexpr kOneThird { 1.0 / 3.0 };
    }

    /*=============================================================================
    |
    |  Description:  Compute the smooth earth diffraction loss using the 
    |                Vogler 3-radii method
    |
    |        Input:  diffractPathLength_m          - Path distance, in meters
    |                f__mhz            - Frequency, in MHz
    |                effEarthRadius_km        - Effective earth radius, in meters
    |                angularDist_LoS_rad         - Angular distance of line-of-sight region
    |                d_hzn__meter[2]   - Horizon distances, in meters
    |                h_e__meter[2]     - Effective terminal heights, in meters
    |                Z_g               - Complex ground impedance
    |
    |      Outputs:  [None]
    |
    |      Returns:  A_r__db           - Smooth-earth diffraction loss, in dB
    |
    *===========================================================================*/
    double ItmCommonCalculator::calcSmoothEarthDiffractLoss_dB(const double& diffractPathLength_m, const double& effEarthRadius_km, const double& angularDist_LoS_rad) {
        const double& txHorizonDist_m = m_itmResults.m_intermResults.m_txHorizonDist_m;
        const double& rxHorizonDist_m = m_itmResults.m_intermResults.m_rxHorizonDist_m;
        const double& txEffHeight_m = m_itmResults.m_intermResults.m_txEffHeight_m;
        const double& rxEffHeight_m = m_itmResults.m_intermResults.m_rxEffHeight_m;

        const double angularDist_nonLoS_rad = diffractPathLength_m / effEarthRadius_km - angularDist_LoS_rad;   // [Algorithm, Eqn 4.12]
        const double actualDist_maxLoS_m = txHorizonDist_m + rxHorizonDist_m;                                   // Maximum line-of-sight distance for actual path

        // compute 3 radii
        double adjEffEarthRadiusList_km[3];
        // which is effEarthRadius_km when angularDist_LoS_rad = d_ML__meter / effEarthRadius_km
        adjEffEarthRadiusList_km[0] = (diffractPathLength_m - actualDist_maxLoS_m) / (diffractPathLength_m / effEarthRadius_km - angularDist_LoS_rad);
        // Compute the radius of the effective earth for terminal j using[Volger 1964, Eqn 3] re - arranged
        adjEffEarthRadiusList_km[1] = 0.5 * txHorizonDist_m * txHorizonDist_m / txEffHeight_m;
        adjEffEarthRadiusList_km[2] = 0.5 * rxHorizonDist_m * rxHorizonDist_m / rxEffHeight_m;

        double kValueList[3];
        double b0List[3];
        double earthRadiusConstList[3];
        for (std::size_t arrayInd = 0; arrayInd < 3; arrayInd++)
        {
            const double earthRadius_km = 1.0 / kActualEarthCurvature_perMeter;
            // C_0 is the ratio of the 4/3 earth to effective earth (technically Vogler 1964 ratio is 4/3 to effective earth k value), all raised to the (1/3) power.
            // C_0 = (4 / 3k) ^ (1 / 3) [Vogler 1964, Eqn 2]
            earthRadiusConstList[arrayInd] = std::pow((4.0 / 3.0) * earthRadius_km / adjEffEarthRadiusList_km[arrayInd], kOneThird);

            // [Vogler 1964, Eqn 6a / 7a]
            kValueList[arrayInd] = 0.017778 * earthRadiusConstList[arrayInd] * std::pow(m_freq_MHz, -kOneThird) / std::abs(m_groundImpedance);

            // Compute B_0 for each radius
            // [Vogler 1964, Fig 4]
            b0List[arrayInd] = 1.607 - kValueList[arrayInd];
        }

        double diffractDistList_km[3];
        diffractDistList_km[0] = (adjEffEarthRadiusList_km[0] * angularDist_nonLoS_rad) * 1.0e-3;   // angular distance of the "diffraction path"
        diffractDistList_km[1] = txHorizonDist_m * 1.0e-3;
        diffractDistList_km[2] = rxHorizonDist_m * 1.0e-3;

        double inputDistList_km[3];
        // Compute inputDistList_km for each radius [Vogler 1964, Eqn 2]
        const double freqPowerTerm = std::pow(m_freq_MHz, kOneThird);
        inputDistList_km[1] = b0List[1] * earthRadiusConstList[1] * earthRadiusConstList[1] * freqPowerTerm * diffractDistList_km[1];
        inputDistList_km[2] = b0List[2] * earthRadiusConstList[2] * earthRadiusConstList[2] * freqPowerTerm * diffractDistList_km[2];
        inputDistList_km[0] = b0List[0] * earthRadiusConstList[0] * earthRadiusConstList[0] * freqPowerTerm * diffractDistList_km[0] + 
                    inputDistList_km[1] + inputDistList_km[2];

        // Compute height gain functions for Tx & Rx
        const double& txGainHeight_dB = ItmHelpers::calcSmoothEarthGainHeight_dB(inputDistList_km[1], kValueList[1]);
        const double& rxGainHeight_dB = ItmHelpers::calcSmoothEarthGainHeight_dB(inputDistList_km[2], kValueList[2]);

        // Compute distance gain function
        const double gainDist_dB = 0.05751 * inputDistList_km[0] - 10.0 * std::log10(inputDistList_km[0]);  // [TN101, Eqn 8.4] & [Volger 1964, Eqn 13]

        return gainDist_dB - txGainHeight_dB - rxGainHeight_dB - 20.0;                                      // [Algorithm, Eqn 4.20] & [Volger 1964]
    }

    namespace ItmHelpers {
        /*=============================================================================
        |
        |  Description:  Height Function, F(x, K) for smooth earth diffraction
        |
        |        Input:  inputDist_km          - Normalized distance, in meters
        |                K              - K value
        |
        |      Outputs:  [None]
        |
        |      Returns:  F(x, K)        - in dB
        |
        *===========================================================================*/
        double calcSmoothEarthGainHeight_dB(const double& inputDist_km, const double& kValue) {
            if (inputDist_km < 200.0) {
                // TODO(vmartin): Is this supposed to be a log10 instead of log?
                const double w = -std::log(kValue);

                if (kValue < 1e-5 || inputDist_km * w * w * w > 5495.0) {
                    return (inputDist_km > 1.0) ? 17.372 * std::log(inputDist_km) - 117.0 : -117.0;
                }
                else {
                    return 2.5e-5 * inputDist_km * inputDist_km / kValue - 8.686 * w - 15.0;
                }
            }
            else {
                const double intermResult = 0.05751 * inputDist_km - 4.343 * std::log(inputDist_km);

                if (inputDist_km < 2.0e3) {
                    const double w = 0.0134 * inputDist_km * exp(-0.005 * inputDist_km);
                    return (1.0 - w) * intermResult + w * (17.372 * std::log(inputDist_km) - 117.0);
                }

                return intermResult;
            }
        }
    }
}