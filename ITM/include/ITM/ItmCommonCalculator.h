#ifndef ITM_COMMON_CALCULATOR_H
#define ITM_COMMON_CALCULATOR_H

#include <ITM/ItmConstructs.h>

#include <complex>
#include <iostream>
#include <sstream>
#include <vector>

namespace NTIA::ITM {
    namespace {
        // NOTE: WGS-84 mean Earth radius is 6371008.7714 meters
        double constexpr kActualEarthCurvature_perMeter { 1.0 / 6371008.7714 };
        double constexpr kDefaultMaxLoss_dB { 999.0 };
    }

    class ItmCommonCalculator {
    public:
        /// @brief Construct generic ITM calculator for calling the model in either point-to-point or area mode
        /// @param txHeight_m Structural height of Tx (meters)
        /// @param rxHeight_m Structural height of Rx (meters)
        /// @param climateCode Radio climate
        /// @param refractivity_N Refractivity (N-units)
        /// @param freq_MHz Frequency (MHz)
        /// @param isTxHorizPolariz Indicates transmitter antenna polarization (true = horizontal, false = vertical)
        /// @param relPermittivity Relative permittivity
        /// @param conductivity Conductivity
        /// @param varMode Mode of variability
        /// @param timePercent Time percentage (0 < time < 100%)
        /// @param locationPercent Location percentage (0 < location < 100%)
        /// @param situationPercent Situation percentage (0 < situation < 100%)
        /// @param performValidation Optional parameter indicating whether validation should be performed (toggle off to improve speed)
        ItmCommonCalculator(const double& txHeight_m, const double& rxHeight_m, 
                const RadioClimate& climateCode, const double& refractivity_N, const double& freq_MHz,
                const bool isTxHorizPolariz, const double& relPermittivity, const double& conductivity, 
                const VariabilityMode& varMode, const double& timePercent, const double& locationPercent, const double& situationPercent,
                const bool performValidation = true) : 
                    m_txHeight_m(txHeight_m), m_rxHeight_m(rxHeight_m), 
                    m_radioClimate(climateCode), m_refractivity_N(refractivity_N), m_freq_MHz(freq_MHz) {
            if (performValidation) {
                validateInputs();
            }
        }

        /// @brief The ITS Irregular Terrain Model (ITM).
        /// This function exposes point-to-point mode functionality, 
        /// with variability specified with time/location/situation (TLS)
        /// @param terrainHeightList_m List of terrain heights along path between Tx --> Rx (meters)
        /// @param terrainSampleResolution_m Sample resolution between successive terrain height values in terrainHeightList_m (meters)
        /// @return Results struct containing ITM basic transmission loss (dB) and various intermediate calculated values
        ItmResults calcItmLoss_P2P_dB(const std::vector<double>& terrainHeightList_m, 
                    const double& terrainSampleResolution_m);
        
        /// @brief The ITS Irregular Terrain Model (ITM).
        /// This function exposes area mode functionality, 
        /// with variability specified with time/location/situation (TLS)
        /// @param txSitingCriteria Tx siting criteria (indicating how well the Tx was sited to communicate with the Rx)
        /// @param rxSitingCriteria Rx siting criteria (indicating how well the Rx was sited to communicate with the Tx)
        /// @param dist_km Path length (km)
        /// @param terrainIrregularityParam_m Parameter indicating how much the regional terrain fluctuates over space (meters)
        /// @return Results struct containing ITM basic transmission loss (dB) and various intermediate calculated values
        ItmResults calcItmLoss_area_dB(const SitingCriteria& txSitingCriteria, const SitingCriteria& rxSitingCriteria, const double& dist_km,
                const double& terrainIrregularityParam_m);
    private:
        void validateInputs() {
            std::ostringstream oStrStream;
            if (m_txHeight_m < 1.0 || m_txHeight_m > 1.0e3) {
                std::cerr << "WARNING: ItuCommonCalculator::validateInputs(): " 
                            << "ITM was only designed to support transmitter heights between 1 m < txHeight_m < 1 km (txHeight_m = " 
                            << m_txHeight_m << ")";
            }
            if (m_rxHeight_m < 1.0 || m_rxHeight_m > 1.0e3) {
                std::cerr << "WARNING: ItuCommonCalculator::validateInputs(): " 
                            << "ITM was only designed to support receiver heights between 1 m < rxHeight_m < 1 km (rxHeight_m = " 
                            << m_rxHeight_m << ")";
            }
            if (m_txHeight_m < 0.5 || m_txHeight_m > 3.0e3) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support transmitter heights outside of the range 0.5 m < txHeight_m < 3 km (txHeight_m = " 
                            << m_txHeight_m << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_rxHeight_m < 0.5 || m_rxHeight_m > 3.0e3) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support receiver heights outside of the range 0.5 m < rxHeight_m < 3 km (rxHeight_m = " 
                            << m_rxHeight_m << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_refractivity_N < 250.0 || m_refractivity_N > 400.0) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support refractivity values outside of the range 250 < N < 400 (N = " 
                            << m_refractivity_N << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_freq_MHz < 40.0 || m_freq_MHz > 10.0e3) {
                std::cerr << "WARNING: ItuCommonCalculator::validateInputs(): " 
                            << "ITM was only designed to support frequencies between 40 MHz < freq_MHz < 10 GHz (freq_MHz = " 
                            << m_freq_MHz << ")";
            }
            if (m_freq_MHz < 20.0|| m_freq_MHz > 20.0e3) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support frequencies outside of the range 20 MHz < freq_MHz < 20 GHz (freq_MHz = " 
                            << m_freq_MHz << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_relPermittivity < 1.0) {
                std::ostringstream oStrStream;
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support relative permittivity values < 1 (relPermittivity = " 
                            << m_relPermittivity << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_conductivity <= 0.0) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support conductivity values <= 0 (conductivity = " 
                            << m_conductivity << ")";
                throw std::domain_error(oStrStream.str());
            }
            if (m_timePercent <= 0.0 || m_timePercent >= 100.0) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support time percentages outside of the range 0 < timePercent < 100 (timePercent = " 
                            << m_timePercent << ")";
                throw std::domain_error(oStrStream.str());   
            }
            if (m_locationPercent <= 0.0 || m_locationPercent >= 100.0) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support location percentages outside of the range 0 < locationPercent < 100 (locationPercent = " 
                            << m_locationPercent << ")";
                throw std::domain_error(oStrStream.str());   
            }
            if (m_situationPercent <= 0.0 || m_situationPercent >= 100.0) {
                oStrStream << "ERROR: ItuCommonCalculator::validateInputs(): " 
                            << "ITM does not support situation percentages outside of the range 0 < situationPercent < 100 (situationPercent = " 
                            << m_situationPercent << ")";
                throw std::domain_error(oStrStream.str());   
            }
        }

        void validateIntermValues() {
            /*
            // Check validity of small angle approximation
            if (abs(theta_hzn[0]) > 200e-3)
                *warnings |= WARN__TX_HORIZON_ANGLE;
            if (abs(theta_hzn[1]) > 200e-3)
                *warnings |= WARN__RX_HORIZON_ANGLE;

            // Checks that the actual horizon distance can't be less than 1/10 of the smooth earth horizon distance
            if (d_hzn__meter[0] < 0.1 * d_hzn_s__meter[0])
                *warnings |= WARN__TX_HORIZON_DISTANCE_1;
            if (d_hzn__meter[1] < 0.1 * d_hzn_s__meter[1])
                *warnings |= WARN__RX_HORIZON_DISTANCE_1;

            // Checks that the actual horizon distance can't be greater than 3 times the smooth earth horizon distance
            if (d_hzn__meter[0] > 3.0 * d_hzn_s__meter[0])
                *warnings |= WARN__TX_HORIZON_DISTANCE_2;
            if (d_hzn__meter[1] > 3.0 * d_hzn_s__meter[1])
                *warnings |= WARN__RX_HORIZON_DISTANCE_2;

            // Check the surface refractivity
            if (N_s < 150)
                return ERROR__SURFACE_REFRACTIVITY_SMALL;
            if (N_s > 400)
                return ERROR__SURFACE_REFRACTIVITY_LARGE;
            if (N_s < 250) // 150 <= N_s < 250
                *warnings |= WARN__SURFACE_REFRACTIVITY;

            // Check effective earth size
            if (a_e__meter < 4000000 || a_e__meter > 13333333)
                return ERROR__EFFECTIVE_EARTH;

            // Check ground impedance
            if (Z_g.real() <= abs(Z_g.imag()))
                return ERROR__GROUND_IMPEDANCE;

            const double minPathDist_m = std::abs(m_itmResults.m_intermResults.m_txEffHeight_m - m_itmResults.m_intermResults.m_rxEffHeight_m) / 0.2;

            const double pathDist_m = m_itmResults.m_intermResults.m_terrainProfile.m_pathDist_km * 1.0e3;

            if (pathDist_m < minPathDist_m)
                *warnings |= WARN__PATH_DISTANCE_TOO_SMALL_1;
            if (pathDist_m < 1.0e3)
                *warnings |= WARN__PATH_DISTANCE_TOO_SMALL_2;
            if (pathDist_m > 1.0e6)
                *warnings |= WARN__PATH_DISTANCE_TOO_BIG_1;
            if (pathDist_m > 2.0e6)
                *warnings |= WARN__PATH_DISTANCE_TOO_BIG_2;
            */
        }

        void initialize_P2P(const double& avgPathHeightAmsl_m);
        void setHorizonParameters(const double& effEarthRadius_km);
        void calcHorizonParameters();
        double calcTerrainIrreg_m(const double& distToStart_m, const double& distToEnd_m);
        double calcLongleyRiceLoss_dB(PropagationMode& propMode);
        double calcSmoothEarthDiffractLoss_dB(const double& diffractPathLength_m, const double& effEarthRadius_km, const double& angularDist_LoS_rad);
        double calcKnifeEdgeDiffractLoss_dB(const double& inputDist_m, const double& effEarthRadius_km, const double& angularDist_LoS_rad);
        double calcTroposcatterLoss_dB(const double& tropoPathLength_m, const double& earthEffRadius_m, 
                const double& angularDist_LoS_rad, double& initialH0_dB);

        // Initial parameters
        double m_txHeight_m;
        double m_rxHeight_m;
        RadioClimate m_radioClimate;
        double m_refractivity_N;
        double m_freq_MHz;
        double m_isTxHorizPolariz;
        double m_relPermittivity;
        double m_conductivity;
        double m_varMode;
        double m_timePercent;
        double m_locationPercent;
        double m_situationPercent;

        // Intermediate parameters
        std::complex<double> m_groundImpedance;
        double m_surfaceRefractivity_N;     // Surface refractivity, in N-Units
        double m_effEarthCurvature_perM;    // Curvature of the effective earth

        // Output parameters (updated by each member function)
        ItmResults m_itmResults;
    };
} // end namespace

#endif // ITM_COMMON_CALCULATOR_H