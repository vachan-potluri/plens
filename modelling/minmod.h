/**
 * @file minmod.h
 * @brief Class for Minmod slope limiter
 */

// header guard including namespace because minmod is a function in utilities namespace.
#ifndef SLOPE_LIMITERS_MINMOD_H
#define SLOPE_LIMITERS_MINMOD_H

#include "slope_limiter.h"

namespace slope_limiters
{

/**
 * @class Minmod
 * Minmod slope limiter. OpenFOAM's implementation is followed. See the documentation of
 * SubcellInterpolator and BTP-2 report Appendix.
 */
class Minmod: public SlopeLimiter
{
    public:
    Minmod(): SlopeLimiter() {};
    virtual double value(const double slope_ratio) const override
    {
        if(slope_ratio < 0) return 0;
        else if(slope_ratio < 1) return slope_ratio;
        else return 1;
    }
};

} // namespace slope_limiters

#endif
