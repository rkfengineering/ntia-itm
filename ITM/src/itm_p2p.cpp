#include <ITM/ItmCommonCalculator.h>
#include <complex>

namespace NTIA::ITM {
    ItmResults ItmCommonCalculator::calcItmLoss_P2P_dB(const std::vector<double>& terrainHeightList_m, 
                const double& terrainSampleResolution_m) {
        double N_s;                 // Surface refractivity, in N-Units
        double gamma_e;             // Curvature of the effective earth
        double delta_h__meter;      // Terrain irregularity parameter
        double d__meter;            // Path distance, in meters
        std::complex<double> Z_g;   // Ground impedance
        double theta_hzn[2];        // Terminal horizon angles
        double d_hzn__meter[2];     // Terminal horizon distances
        double h_e__meter[2];       // Terminal effective heights

        IntermResults intermResults;

        const std::size_t numPoints = terrainHeightList_m.size() - 1u;
        intermResults.m_pathDist_km = static_cast<double>(numPoints) * terrainSampleResolution_m;
        intermResults.m_pathDist_km *= 1.0e-3;

        // Calculate average path height, ignoring first & last 10% of the path
        const std::size_t oneTenthNumPoints = 0.1 * static_cast<double>(numPoints);
        double sysHeightAmsl_m = 0;
        for (std::size_t pointInd = oneTenthNumPoints; pointInd <= numPoints - oneTenthNumPoints; pointInd++)
            sysHeightAmsl_m += terrainHeightList_m[pointInd];
        sysHeightAmsl_m /= static_cast<double>(numPoints - 2u * oneTenthNumPoints + 1u);

        InitializePointToPoint(f__mhz, h_sys__meter, N_0, pol, epsilon, sigma, &Z_g, &gamma_e, &N_s);

        double h__meter[2] = { h_tx__meter, h_rx__meter };
        QuickPfl(pfl, gamma_e, h__meter, theta_hzn, d_hzn__meter, h_e__meter, &delta_h__meter, &d__meter);

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