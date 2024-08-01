#ifndef ITM_COMMON_CALCULATOR_H
#define ITM_COMMON_CALCULATOR_H

#include <ITM/ItmConstructs.h>

#include <iostream>
#include <sstream>
#include <vector>

namespace NTIA::ITM {
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
    }
}

#endif // ITM_COMMON_CALCULATOR_H