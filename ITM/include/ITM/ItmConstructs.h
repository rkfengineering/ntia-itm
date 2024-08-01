#ifndef NTIA_ITM_CONSTRUCTS_H
#define NTIA_ITM_CONSTRUCTS_H

namespace NTIA::ITM {
    struct IntermResults
    {
        double m_txHorizonAngle_deg;    // Terminal horizon angles, in degrees
        double m_rxHorizonAngle_deg;
        double m_txHorizonDist_m;       // Terminal horizon distances, in meters
        double m_rxHorizonDist_m;
        double m_txEffHorizonDist_m;    // Terminal effective heights, in meters
        double m_rxEffHorizonDist_m;
        double m_surfRefract_N;         // Surface refractivity, in N-Units
        double m_terrainIrreg_m;        // Terrain irregularity parameter, in meters
        double m_refAtten_dB;           // Reference attenuation, in dB
        double m_fsplAtten_dB;          // Free space basic transmission loss, in dB
        double m_pathDist_km;           // Path distance, in km
        PropagationMode m_propMode;     // Mode of propagation value
    };

    struct ItmResults {
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