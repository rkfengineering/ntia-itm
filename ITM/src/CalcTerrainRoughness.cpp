#include <ITM/ItmHelpers.h>

namespace NTIA::ITM::ItmHelpers {
    double calcTerrainRoughness_m(const double& pathDist_m, const double& terrainIrreg_m) {
        // [ERL 79-ITS 67, Eqn 3], with distance in meters instead of kilometers
        return terrainIrreg_m * (1.0 - 0.8 * std::exp(-pathDist_m / 50.0e3));
    }
}