#include <ITM/ItmCommonCalculator.h>

#include <algorithm>

namespace NTIA::ITM {
    double ItmCommonCalculator::calcLongleyRiceLoss_dB(PropagationMode& propMode, const bool isP2P) {
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
        const double attenDiffract3_dB = calcDiffractLoss_dB(diffractDist3_m, effEarthRadius_m, isP2P, angularDistInLoS_rad, smoothEarthDist_maxLoS_m);
        const double attenDiffract4_dB = calcDiffractLoss_dB(diffractDist4_m, effEarthRadius_m, isP2P, angularDistInLoS_rad, smoothEarthDist_maxLoS_m);

        // Compute the slope and intercept of the diffraction line
        const double diffractLineSlope = (attenDiffract4_dB - attenDiffract3_dB) / (diffractDist4_m - diffractDist3_m);
        const double diffractLineIntercept_dB = attenDiffract3_dB - diffractLineSlope * diffractDist3_m;

        const double pathDist_m = m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km * 1.0e3;

        // if the path distance is less than the maximum smooth earth line of sight distance...
        if (pathDist_m < smoothEarthDist_maxLoS_m)
        {
            const double& txEffHeight_m = m_itmResults.m_intermResults.m_txEffHeight_m;
            const double& rxEffHeight_m = m_itmResults.m_intermResults.m_rxEffHeight_m;
            // Compute the diffraction loss at the maximum smooth earth line of sight distance
            const double diffractLoss_smoothEarth_maxLoS_dB = smoothEarthDist_maxLoS_m * diffractLineSlope + diffractLineIntercept_dB;

            // [ERL 79-ITS 67, Eqn 3.16a], in meters instead of km and with MIN() part below
            double diffractDist0_m = 0.04 * m_freq_MHz * txEffHeight_m * rxEffHeight_m;
            double diffractDist1_m = 0.0;
            if (diffractLineIntercept_dB >= 0.0)
            {
                diffractDist0_m = std::min({diffractDist0_m, 0.5 * actualDist_maxLoS_m});               // other part of [ERL 79-ITS 67, Eqn 3.16a]
                diffractDist1_m = diffractDist0_m + 0.25 * (actualDist_maxLoS_m - diffractDist0_m);     // [ERL 79-ITS 67, Eqn 3.16d]
            }
            else
                diffractDist1_m = std::max({-diffractLineIntercept_dB / diffractLineSlope, 0.25 * actualDist_maxLoS_m});

            const double losLoss1_dB = calcLineOfSightLoss_dB(diffractDist1_m, diffractLineSlope, diffractLineIntercept_dB, smoothEarthDist_maxLoS_m);

            bool foundPositiveValues = false;

            double kHat1_dBPerM = 0.0, kHat2_dBPerM = 0.0;

            if (diffractDist0_m < diffractDist1_m) {
                const double losLoss0_dB = calcLineOfSightLoss_dB(diffractDist0_m, diffractLineSlope, diffractLineIntercept_dB, smoothEarthDist_maxLoS_m);

                // TODO(vmartin): Is this log supposed to be a log10??
                const double q = std::log(smoothEarthDist_maxLoS_m / diffractDist0_m);

                // [ERL 79-ITS 67, Eqn 3.20]
                const double kHat2_part2_numer = (smoothEarthDist_maxLoS_m - diffractDist0_m) * (losLoss1_dB - losLoss0_dB) - 
                            (diffractDist1_m - diffractDist0_m) * (diffractLoss_smoothEarth_maxLoS_dB - losLoss0_dB);
                // TODO(vmartin): Is this log supposed to be a log10??
                const double kHat2_part2_denom = (smoothEarthDist_maxLoS_m - diffractDist0_m) * std::log(diffractDist1_m / diffractDist0_m) - 
                            (diffractDist1_m - diffractDist0_m) * q;
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
            const double finalLoss_dB = intermAtten_dB + kHat1_dBPerM * pathDist_m + kHat2_dBPerM * log(pathDist_m);
            propMode = PropagationMode::LineOfSight;

            // Don't allow a negative loss
            return std::max({finalLoss_dB, 0.0});
        }
        else {
            // this is a trans-horizon path
            // select to points far into the troposcatter region
            double tropoDist5_m = actualDist_maxLoS_m + 200.0e3;
            double tropoDist6_m = actualDist_maxLoS_m + 400.0e3;

            // Compute the troposcatter loss at the two distances
            double currentH0_dB = -1.0;
            double attenTropo6_dB = calcTroposcatterLoss_dB(tropoDist6_m, effEarthRadius_m, angularDistInLoS_rad, &currentH0_dB);
            double attenTropo5_dB = calcTroposcatterLoss_dB(tropoDist5_m, effEarthRadius_m, angularDistInLoS_rad, &currentH0_dB);

            double tropoLineSlope, tropoLineIntercept_dB, diffractTropoTransitionDist_m;

            // if we got a reasonable prediction value back...
            if (attenTropo5_dB < 1.0e3) {
                // Compute the slope of the troposcatter line
                tropoLineSlope = (attenTropo6_dB - attenTropo5_dB) / 200.0e3;

                // Find the diffraction-troposcatter transition distance
                diffractTropoTransitionDist_m = std::max({std::max({smoothEarthDist_maxLoS_m, 
                            actualDist_maxLoS_m + 1.088 * pow(pow(effEarthRadius_m, 2) / m_freq_MHz, 1.0 / 3.0) * log(m_freq_MHz)}), 
                            (attenTropo5_dB - diffractLineIntercept_dB - tropoLineSlope * tropoDist5_m) / (diffractLineSlope - tropoLineSlope)});

                // Compute the intercept of the troposcatter line
                tropoLineIntercept_dB = (diffractLineSlope - tropoLineSlope) * diffractTropoTransitionDist_m + diffractLineIntercept_dB;
            }
            else {
                // troposcatter gives no real results - so use diffraction line parameters for tropo line
                tropoLineSlope = diffractLineSlope;
                tropoLineIntercept_dB = diffractLineIntercept_dB;
                diffractTropoTransitionDist_m = 10e6;
            }

            double finalLoss_dB = 0.0;
            // Determine if its diffraction or troposcatter and compute the loss
            if (pathDist_m > diffractTropoTransitionDist_m)
            {
                finalLoss_dB = tropoLineSlope * pathDist_m + tropoLineIntercept_dB;
                propMode = Troposcatter;
            }
            else
            {
                finalLoss_dB = diffractLineSlope * pathDist_m + diffractLineIntercept_dB;
                propMode = Diffraction;
            }
        
            // Don't allow a negative loss
            return std::max({finalLoss_dB, 0.0});
        }
    }
}