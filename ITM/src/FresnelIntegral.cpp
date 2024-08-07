#include <ITM/ItmHelpers.h>

namespace NTIA::ITM::ItmHelpers {
    double calcFresnelIntegral(const double& nu)
    {
        if (nu < 2.4)
            return 6.02 + 9.11 * nu - 1.27 * nu * nu;   // [TN101v2, Eqn III.24b] and [ERL 79-ITS 67, Eqn 3.27a & 3.27b]
        else
            return 12.953 + 20.0 * std::log10(nu);      // [TN101v2, Eqn III.24c] and [ERL 79-ITS 67, Eqn 3.27a & 3.27b]
    }
}