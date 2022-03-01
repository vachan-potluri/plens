/**
 * @file minmod.h
 * Implements the minmod function
 */

#ifndef MINMOD_H
#define MINMOD_H

#include <initializer_list>
#include <vector>

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



/**
 * Minmod function for multiple arguments
 * @f[
    \text{minmod}(a_1,a_2,\ldots a_n) = \left{
        \begin{array}{cc}
            \sgn (a_1) \min (|a_1|, |a_2|, \ldots |a_n|) &
                \text{if}\ a_1,a_2,\ldots,a_n\ \text{have same sign}
            0 & \text{otherwise}
        \end{array}
    \right.
 * @f]
 * where @f$\sgn(x)@f$ is the signum function.
 */
inline double minmod(std::initializer_list<double> arr)
{
    std::vector<float> signs(arr.size());
    float sign_product(1); // product of all signs
    int i=0;
    for(const double a : arr){
        signs[i] = (a>0 ? 1.0 : -1.0);
        sign_product *= signs[i];
        i++;
    }

    if(sign_product > 0){
        std::vector<double> arr_abs(arr.size());
        i = 0;
        for(const double a : arr){
            arr_abs[i] = fabs(a);
            i++;
        }
        return *std::min_element(arr_abs.begin(), arr_abs.end());
    }
    else return 0;
}

} // namespace utilities

#endif
