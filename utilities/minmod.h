/**
 * @file minmod.h
 * Implements the minmod function
 */

#ifndef MINMOD_H
#define MINMOD_H

namespace utilities{

/**
 * Minmod function for two arguments
 * @f[
    \text{minmod}(a,b) = \left{
        \begin{array}{cc}
            0 & \text{if}\ ab < 0
            \sgn (a) \min (|a|, |b|) & \text{otherwise}
        \end{array}
    \right.
 * @f]
 * where @f$\sgn(a)@f$ is the signum function.
 */
inline double minmod(const double a, const double b)
{
    if(a*b <= 0) return 0;
    else{
        const double a_abs = fabs(a), b_abs = fabs(b);
        const float sign  = (a>0 ? 1 : -1);
        return std::min(a_abs, b_abs)*sign;
    }
}

} // namespace utilities

#endif
