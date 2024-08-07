#include "ITM/itm.h"

namespace NTIA::ITM::ItmHelpers {
    double calcFSPL_dB(const double& dist_m, const double& freq_MHz) {
        return 32.45 + 20.0 * std::log10(freq_MHz) + 20.0 * std::log10(dist_m * 1.0e-3);
    }
}