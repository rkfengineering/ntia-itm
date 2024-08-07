#include <ITM/ItmCommonCalculator.h>
#include <complex>

namespace NTIA::ITM {
    ItmResults ItmCommonCalculator::calcItmLoss_P2P_dB(const std::vector<double>& terrainHeightList_m, 
                const double& terrainSampleResolution_m) {
        // Zero out / reset ITM results object
        m_itmResults = ItmResults();

        // Populate terrainProfile
        m_itmResults.m_intermResults.m_terrainProfile.m_sampleResolution_m = terrainSampleResolution_m;
        m_itmResults.m_intermResults.m_terrainProfile.m_terrainHeightList_m = terrainHeightList_m;
        m_itmResults.m_intermResults.m_terrainProfile.m_numPointsMinusTx = terrainHeightList_m.size() - 1u;
        
        // For ease of reference in the code
        const std::size_t& numPointsMinusTx = m_itmResults.m_intermResults.m_terrainProfile.m_numPointsMinusTx;
        const double numPointsMinusTx_double = static_cast<double>(numPointsMinusTx);
        
        m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km = numPointsMinusTx_double * terrainSampleResolution_m;
        m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km *= 1.0e-3;

        // Calculate average path height, ignoring first & last 10% of the path
        const std::size_t oneTenthNumPoints = 0.1 * numPointsMinusTx_double;
        double avgPathHeightAmsl_m = 0;
        for (std::size_t pointInd = oneTenthNumPoints; pointInd <= numPointsMinusTx - oneTenthNumPoints; pointInd++) {
            avgPathHeightAmsl_m += terrainHeightList_m[pointInd];
        }
        avgPathHeightAmsl_m /= static_cast<double>(numPointsMinusTx - 2u * oneTenthNumPoints + 1u);

        initialize_P2P(avgPathHeightAmsl_m);
        calcHorizonParameters();

        // Reference attenuation, in dB
        double A_ref__db = 0;
        int propmode = MODE__NOT_SET;
        rtn = LongleyRice(theta_hzn, f__mhz, Z_g, d_hzn__meter, h_e__meter, gamma_e, N_s, delta_h__meter, h__meter, d__meter, MODE__P2P, 
            &A_ref__db, warnings, &propmode);
        if (rtn != SUCCESS)
            return rtn;

        double A_fs__db = FreeSpaceLoss(d__meter, f__mhz);

        // switch from percentages to ratios
        const double timeFrac = m_timePercent / 100.0;
        const double locationFrac = m_locationPercent / 100.0;
        const double situationFrac = m_situationPercent / 100.0;

        *A__db = Variability(timeFrac, locationFrac, situationFrac, h_e__meter, delta_h__meter, f__mhz, d__meter, A_ref__db, climate, mdvar, warnings) + A_fs__db;

        // Save off intermediate values
        interValues->A_ref__db = A_ref__db;
        interValues->A_fs__db = A_fs__db;
        interValues->delta_h__meter = delta_h__meter;
        interValues->d_hzn__meter[0] = d_hzn__meter[0];
        interValues->d_hzn__meter[1] = d_hzn__meter[1];
        interValues->h_e__meter[0] = h_e__meter[0];
        interValues->h_e__meter[1] = h_e__meter[1];
        interValues->N_s = N_s;
        interValues->theta_hzn[0] = theta_hzn[0];
        interValues->theta_hzn[1] = theta_hzn[1];
        interValues->mode = propmode;

        if (*warnings != NO_WARNINGS)
            return SUCCESS_WITH_WARNINGS;

        return SUCCESS;
    }
} // end namespace