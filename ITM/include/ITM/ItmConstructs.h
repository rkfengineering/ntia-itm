#ifndef NTIA_ITM_CONSTRUCTS_H
#define NTIA_ITM_CONSTRUCTS_H

namespace NTIA::ITM {
    struct TerrainProfile {
        // Default construct will zero out all values
        TerrainProfile() = default;
        
        std::vector<double> m_terrainHeightList_m;  // Terrain heights along the path (first ind = Tx --> last ind = Rx)
        std::size_t m_numPointsMinusTx;             // Number of points in the path, not including the Tx
        double m_pathDist_km;                       // Path distance, in km
        double m_sampleResolution_m;                // Sampling resolution between terrain heights, in meters
    };

    struct IntermResults
    {   
        // Default constructor will zero out all values
        IntermResults() = default;

        double m_txHorizonAngle_rad;        // Terminal horizon angles, in degrees
        double m_rxHorizonAngle_rad;    
        double m_txHorizonDist_m;           // Terminal horizon distances, in meters
        double m_rxHorizonDist_m;   
        double m_txEffHorizonDist_m;        // Terminal effective horizon distances, in meters
        double m_rxEffHorizonDist_m;    
        double m_txEffHeight_m;             // Terminal effective heights, in meters
        double m_rxEffHeight_m;
        double m_surfRefract_N;             // Surface refractivity, in N-Units
        double m_terrainIrreg_m;            // Terrain irregularity parameter, in meters
        double m_refAtten_dB;               // Reference attenuation, in dB
        double m_fsplAtten_dB;              // Free space basic transmission loss, in dB
        TerrainProfile m_terrainProfile;    // Terrain profile along Tx --> Rx path
        PropagationMode m_propMode;         // Mode of propagation value
    };

    struct ItmResults {
        // Default constructor will zero out all values
        ItmResults() = default;

        double m_atten_dB;
        IntermResults m_intermResults;
    };
    
    /// @brief Tx & Rx siting criteria required as an input to area-mode ITM calculations
    enum SitingCriteria {
        Random,
        Careful,
        VeryCareful
    };

    enum VariabilityMode {
        SingleMessageMode,
        AccidentalMode,
        MobileMode,
        BroadcastMode
    };

    enum PropagationMode {
        NotSet,
        LineOfSight,
        Diffraction,
        Troposcatter
    };

    enum RadioClimate {
        Equatorial,
        ContinentalSubtropical,
        MaritimeSubtropical,
        Desert,
        Temperate,
        MaritimeTemperateOverLand,
        MaritimeTemperateOverSea
    };
}

#endif // NTIA_ITM_CONSTRUCTS_H