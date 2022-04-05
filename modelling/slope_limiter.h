/**
 * @file slope_limiter.h
 * @brief Base class for slope limiter
 */

#ifndef SLOPE_LIMITER_H
#define SLOPE_LIMITER_H

namespace slope_limiters
{

/**
 * @class SlopeLimiter
 * Base class for slope limiter objects. The function value() is overloaded in derived types.
 */
class SlopeLimiter
{
    public:
    SlopeLimiter() = default;
    virtual double value(const double slope_ratio) const = 0;
};

} // namespace slope_limiters

#endif
