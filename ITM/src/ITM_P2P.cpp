#include <ITM/ItmCommonCalculator.h>
#include <ITM/ItmHelpers.h>

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
        PropagationMode propMode = NotSet;
        const double finalLoss_dB = calcLongleyRiceLoss_dB(propMode, true);
        
        m_itmResults.m_intermResults.m_fsplAtten_dB = ItmHelpers::calcFSPL_dB(m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km * 1.0e3, m_freq_MHz);

        // switch from percentages to ratios
        const double timeFrac = m_timePercent / 100.0;
        const double locationFrac = m_locationPercent / 100.0;
        const double situationFrac = m_situationPercent / 100.0;

        m_itmResults.m_atten_dB = Variability(timeFrac, locationFrac, situationFrac, d__meter, finalLoss_dB) + m_itmResults.m_intermResults.m_fsplAtten_dB;

        return m_itmResults;
    }
} // end namespace