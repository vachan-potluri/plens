/**
 * @file none.h
 * @brief Class for "None" slope limiter: the slope of interpolation is always 0
 */

#ifndef NONE_H
#define NONE_H

#include "slope_limiter.h"

namespace slope_limiters
{

/**
 * @class None
 * The value of limiter returned here is always 0 so that the net slope is always 0. This results
 * in a 1st order subcell FV algorithm.
 */
class None: public SlopeLimiter
{
    public:
    None(): SlopeLimiter() {};
    virtual double value(const double slope_ratio) const override
    {
        return 0;
    }
};

} // namespace slope_limiters

#endif