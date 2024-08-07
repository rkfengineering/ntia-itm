#include <ITM/ItmCommonCalculator.h>

/*=============================================================================
 |
 |  Description:  Initialize parameters for point-to-point mode
 |
 |        Input:  f__mhz            - Frequency, in MHz
 |                h_sys__meter      - Average height of the path above
 |                                    mean sea level, in meters
 |                N_0               - Refractivity, in N-Units
 |                pol               - Polarization
 |                                      + 0 : POLARIZATION__HORIZONTAL
 |                                      + 1 : POLARIZATION__VERTICAL
 |                epsilon           - Relative permittivity
 |                sigma             - Conductivity
 |
 |      Outputs:  Z_g               - Complex ground impedance
 |                gamma_e           - Curvature of the effective earth
 |                N_s               - Surface refractivity, in N-Units
 |
 |      Returns:  [None]
 |
 *===========================================================================*/

namespace NTIA::ITM {
    void ItmCommonCalculator::initialize_P2P(const double& avgPathHeightAmsl_m) {
        // Scale local refractivity into a surface refractivity based on the path's average elevation AMSL
        m_surfaceRefractivity_N = (avgPathHeightAmsl_m <= 0.0) 
                    ? m_refractivity_N 
                    : m_refractivity_N * std::exp(-avgPathHeightAmsl_m / 9460.0);   // [TN101, Eq 4.3]
        
        const double effEarthCurvatureScaleTerm = 1.0 - 0.04665 * std::exp(m_surfaceRefractivity_N / 179.3);
        m_effEarthCurvature_perM = kActualEarthCurvature_perMeter * effEarthCurvatureScaleTerm;   // [TN101, Eq 4.4], reworked

        std::complex<double> complexRelPermittivity(m_relPermittivity, 18.0e3 * m_conductivity / m_freq_MHz);

        // Ground impedance for horizontal polarization
        m_groundImpedance = std::sqrt(complexRelPermittivity - 1.0);
        if (m_isTxHorizPolariz) {
            // Adjust for vertical polarization
            m_groundImpedance /= complexRelPermittivity;
        }
    }

}