#include <ITM/ItmHelpers.h>

namespace ITM::NTIA::ItmHelpers {
    double calcSigmaH_m(const double& terrainIrreg_m) {
        // "RMS deviation of terrain and terrain clutter within the limits of the first Fresnel zone in the dominant reflecting plane"
        // [ERL 79-ITS 67, Eqn 3.6a]
        return 0.78 * terrainIrreg_m * std::exp(-0.5 * std::pow(terrainIrreg_m, 0.25));
    }
}