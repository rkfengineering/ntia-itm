#ifndef ITM_CORE_HELPERS_H
#define ITM_CORE_HELPERS_H

#define _USE_MATH_DEFINES
#include <cmath>

namespace NTIA::ITM::ItmHelpers {
    double constexpr kSpeedOfLight_mPerS { 299792458.0 };
    double constexpr kWaveToMHzFreqTerm { kSpeedOfLight_mPerS * 1.0e-6 / (2.0 * M_PI) };

    double calcSmoothEarthGainHeight_dB(const double& inputDist_km, const double& kValue);
    double calcSigmaH_m(const double& terrainIrreg_m);
    double calcTerrainRoughness_m(const double& pathDist_m, const double& terrainIrreg_m);

    /// @brief Approximation of the Fresnel integral, as defined in "6. Addenda - Numerical Approximations" from ITM Algorithm Whitepaper
    /// @param nu Nu is the input the Frensel integral
    /// @return Frensel integration result from nu --> infinity    
    double calcFresnelIntegral(const double& nu);

    /// @brief Curve fit helper for calculating troposcatter frequency gain function, H_0()
    /// @param arrayInd Index of array defined in algorithm document (a & b)
    /// @param rTerm Input parameter defined in algorithm document (r_1 or r_2)
    /// @return Curve fit value from the defined troposcatter frequency gain function's curve (dB)
    double calcTropoFreqGainCurveFit_dB(const std::size_t arrayInd, const double &rTerm);

    /// @brief Troposcatter frequency gain function, H_0(), from [TN101v1, Ch 9.2]
    /// @param rParam Input parameter defined in algorithm document (r_1 or r_2)
    /// @param scatterEfficiency Scatter efficiency found in algorithm document (eta_s)
    /// @return Troposcatter frequency gain (dB)
    double calcTropoFreqGain_dB(const double& rParam, double& scatterEfficiency);

    double calcTropoAttenFunction_dB(const double& inputDist_m);
} // end namespace

#endif // ITM_CORE_HELPERS_H