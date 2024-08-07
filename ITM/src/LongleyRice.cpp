#include <ITM/ItmCommonCalculator.h>

#include <algorithm>

/*=============================================================================
 |
 |  Description:  Compute the reference attenuation, using the 
 |                Longley-Rice method
 |
 |        Input:  theta_hzn[2]      - Terminal horizon angles
 |                m_freq_MHz            - Frequency, in MHz
 |                Z_g               - Complex surface transfer impedance
 |                d_hzn__meter[2]   - Terminal horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                gamma_e           - Curvature of the effective earth
 |                N_s               - Surface refractivity, in N-Units
 |                delta_h__meter    - Terrain irregularity parameter
 |                h__meter[2]       - Terminal structural heights, in meters
 |                pathDist_m          - Path distance, in meters
 |                mode              - Mode of operation (P2P or Area)
 |
 |      Outputs:  A_ref__db         - Reference attenuation, in dB
 |                warnings          - Warning flags
 |                propmode          - Mode of propagation value
 |
 |      Returns:  error             - Error code
 |
 *===========================================================================*/
namespace NTIA::ITM {
    double ItmCommonCalculator::calcLongleyRiceLoss_dB(PropagationMode& propMode) {
        const double effEarthRadius_m = 1.0 / m_effEarthCurvature_perM;

        // Terrestrial smooth earth horizon distance approximation
        const double txSmoothEarthHorizonDist_m = std::sqrt(2.0 * m_itmResults.m_intermResults.m_txEffHeight_m * effEarthRadius_m);
        const double rxSmoothEarthHorizonDist_m = std::sqrt(2.0 * m_itmResults.m_intermResults.m_rxEffHeight_m * effEarthRadius_m);

        // Maximum line-of-sight distance for smooth earth
        double smoothEarthDist_maxLoS_m = txSmoothEarthHorizonDist_m + rxSmoothEarthHorizonDist_m;

        // Maximum line-of-sight distance for actual path
        double actualDist_maxLoS_m = m_itmResults.m_intermResults.m_txHorizonDist_m + m_itmResults.m_intermResults.m_rxHorizonDist_m;

        // Angular distance of line-of-sight region
        const double angularDistInLoS_rad = -std::max({m_itmResults.m_intermResults.m_txHorizonAngle_rad + m_itmResults.m_intermResults.m_rxHorizonAngle_rad, 
                    -actualDist_maxLoS_m / effEarthRadius_m});

        // Select two distances far in the diffraction region
        const double diffractDist3_m = std::max({smoothEarthDist_maxLoS_m, actualDist_maxLoS_m + 5.0 * pow(pow(effEarthRadius_m, 2) / m_freq_MHz, 1.0 / 3.0)});
        const double diffractDist4_m = diffractDist3_m + 10.0 * pow(pow(effEarthRadius_m, 2) / m_freq_MHz, 1.0 / 3.0);

        // Compute the diffraction loss at the two distances
        double A_3__db = DiffractionLoss(diffractDist3_m, d_hzn__meter, h_e__meter, Z_g, effEarthRadius_m, delta_h__meter, h__meter, mode, angularDistInLoS_rad, smoothEarthDist_maxLoS_m, m_freq_MHz);
        double A_4__db = DiffractionLoss(diffractDist4_m, d_hzn__meter, h_e__meter, Z_g, effEarthRadius_m, delta_h__meter, h__meter, mode, angularDistInLoS_rad, smoothEarthDist_maxLoS_m, m_freq_MHz);

        // Compute the slope and intercept of the diffraction line
        const double diffractLineSlope = (A_4__db - A_3__db) / (diffractDist4_m - diffractDist3_m);
        const double diffractLineIntercept_dB = A_3__db - diffractLineSlope * diffractDist3_m;

        const double pathDist_m = m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km * 1.0e3;

        // if the path distance is less than the maximum smooth earth line of sight distance...
        if (pathDist_m < smoothEarthDist_maxLoS_m)
        {
            // Compute the diffraction loss at the maximum smooth earth line of sight distance
            const double diffractLoss_smoothEarth_maxLoS_dB = smoothEarthDist_maxLoS_m * diffractLineSlope + diffractLineIntercept_dB;

            // [ERL 79-ITS 67, Eqn 3.16a], in meters instead of km and with MIN() part below
            double diffractDist0_m = 0.04 * m_freq_MHz * h_e__meter[0] * h_e__meter[1];
            double diffractDist1_m = 0.0;
            if (diffractLineIntercept_dB >= 0.0)
            {
                diffractDist0_m = std::min({diffractDist0_m, 0.5 * actualDist_maxLoS_m});               // other part of [ERL 79-ITS 67, Eqn 3.16a]
                diffractDist1_m = diffractDist0_m + 0.25 * (actualDist_maxLoS_m - diffractDist0_m);     // [ERL 79-ITS 67, Eqn 3.16d]
            }
            else
                diffractDist1_m = std::max({-diffractLineIntercept_dB / diffractLineSlope, 0.25 * actualDist_maxLoS_m});

            const double losLoss1_dB = LineOfSightLoss(diffractDist1_m, h_e__meter, Z_g, delta_h__meter, diffractLineSlope, diffractLineIntercept_dB, smoothEarthDist_maxLoS_m, m_freq_MHz);

            bool foundPositiveValues = false;

            double kHat1_dBPerM = 0.0, kHat2_dBPerM = 0.0;

            if (diffractDist0_m < diffractDist1_m) {
                const double losLoss0_dB = LineOfSightLoss(diffractDist0_m, h_e__meter, Z_g, delta_h__meter, diffractLineSlope, diffractLineIntercept_dB, smoothEarthDist_maxLoS_m, m_freq_MHz);

                // TODO(vmartin): Is this log supposed to be a log10??
                const double q = std::log(smoothEarthDist_maxLoS_m / diffractDist0_m);

                // [ERL 79-ITS 67, Eqn 3.20]
                const double kHat2_part2_numer = (smoothEarthDist_maxLoS_m - diffractDist0_m) * (losLoss1_dB - losLoss0_dB) - (diffractDist1_m - diffractDist0_m) * (diffractLoss_smoothEarth_maxLoS_dB - losLoss0_dB);
                // TODO(vmartin): Is this log supposed to be a log10??
                const double kHat2_part2_denom = (smoothEarthDist_maxLoS_m - diffractDist0_m) * std::log(diffractDist1_m / diffractDist0_m) - (diffractDist1_m - diffractDist0_m) * q;
                kHat2_dBPerM = std::max({0.0, kHat2_part2_numer / kHat2_part2_denom });

                foundPositiveValues = diffractLineIntercept_dB > 0.0 || kHat2_dBPerM > 0.0;
                if (foundPositiveValues) {
                    // [ERL 79-ITS 67, Eqn 3.21]
                    kHat1_dBPerM = (diffractLoss_smoothEarth_maxLoS_dB - losLoss0_dB - kHat2_dBPerM * q) / (smoothEarthDist_maxLoS_m - diffractDist0_m);

                    if (kHat1_dBPerM < 0.0) {
                        kHat1_dBPerM = 0.0;
                        kHat2_dBPerM = std::abs(diffractLoss_smoothEarth_maxLoS_dB - losLoss0_dB) / q;

                        if (kHat2_dBPerM == 0.0) {
                            kHat1_dBPerM = diffractLineSlope;
                        }
                    }
                }
            }

            if (!foundPositiveValues) {
                kHat1_dBPerM = std::abs(diffractLoss_smoothEarth_maxLoS_dB - losLoss1_dB) / (smoothEarthDist_maxLoS_m - diffractDist1_m);
                kHat2_dBPerM = 0.0;

                if (kHat1_dBPerM == 0.0)
                    kHat1_dBPerM = diffractLineSlope;
            }

            // TODO(vmartin): Is this log supposed to be a log10??
            const double intermAtten_dB = diffractLoss_smoothEarth_maxLoS_dB - kHat1_dBPerM * smoothEarthDist_maxLoS_m - kHat2_dBPerM * log(smoothEarthDist_maxLoS_m);

            // [ERL 79-ITS 67, Eqn 3.19]
            // TODO(vmartin): Is this log supposed to be a log10??
            return intermAtten_dB + kHat1_dBPerM * pathDist_m + kHat2_dBPerM * log(pathDist_m);
            propMode = PropagationMode::LineOfSight;
        }
        else {
            // this is a trans-horizon path
            // select to points far into the troposcatter region
            double d_5__meter = actualDist_maxLoS_m + 200.0e3;
            double d_6__meter = actualDist_maxLoS_m + 400.0e3;

            // Compute the troposcatter loss at the two distances
            double h0 = -1.0;
            double A_6__db = TroposcatterLoss(d_6__meter, theta_hzn, d_hzn__meter, h_e__meter, effEarthRadius_m, N_s, m_freq_MHz, angularDistInLoS_rad, &h0);
            double A_5__db = TroposcatterLoss(d_5__meter, theta_hzn, d_hzn__meter, h_e__meter, effEarthRadius_m, N_s, m_freq_MHz, angularDistInLoS_rad, &h0);

            double M_s, A_s0__db, d_x__meter;

            // if we got a reasonable prediction value back...
            if (A_5__db < 1000.0)
            {
                // Compute the slope of the troposcatter line
                M_s = (A_6__db - A_5__db) / 200e3;

                // Find the diffraction-troposcatter transition distance
                d_x__meter = MAX(MAX(smoothEarthDist_maxLoS_m, actualDist_maxLoS_m + 1.088 * pow(pow(effEarthRadius_m, 2) / m_freq_MHz, 1.0 / 3.0) * log(m_freq_MHz)), (A_5__db - diffractLineIntercept_dB - M_s * d_5__meter) / (diffractLineSlope - M_s));

                // Compute the intercept of the troposcatter line
                A_s0__db = (diffractLineSlope - M_s) * d_x__meter + diffractLineIntercept_dB;
            }
            else
            {
                // troposcatter gives no real results - so use diffraction line parameters for tropo line
                M_s = diffractLineSlope;
                A_s0__db = diffractLineIntercept_dB;
                d_x__meter = 10e6;
            }

            // Determine if its diffraction or troposcatter and compute the loss
            if (pathDist_m > d_x__meter)
            {
                *A_ref__db = M_s * pathDist_m + A_s0__db;
                *propmode = MODE__TROPOSCATTER;
            }
            else
            {
                *A_ref__db = diffractLineSlope * pathDist_m + diffractLineIntercept_dB;
                *propmode = MODE__DIFFRACTION;
            }
        }

        // Don't allow a negative loss
        *A_ref__db = MAX(*A_ref__db, 0.0);

        return SUCCESS;
    }
}