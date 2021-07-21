/**
 * @file rk4_stage5_register3.h
 * @brief A class for 5 stage, 3 register RK4 time integration.
 */

#ifndef RK4_STAGE5_REGISTER3_H
#define RK4_STAGE5_REGISTER3_H

#include <array>

/**
 * @class RK4Stage5Register3
 * This class provides Butcher coefficients for 5 stage, 3 register RK4 update. The coefficients
 * are taken from Kennedy, Carpenter & Lewis (2000). See WJ-09-Jun-2021 for details.
 *
 * Currently, only the "N" version of this scheme is used. In future, this class can be extended
 * to provide "M" or "C" version also by taking a constructor argument. And hence the data is
 * stored in a class, and not in some namespace. The coefficients are stored as public variables
 * and are calculated in the constructor. Thus a `const` instance to this class is recommended for
 * use.
 *
 * For a 5 stage 3 register RK4 scheme, the Butcher coefficients take the form
 * $$
 * \begin{array}{r|ccccc}
 * c_0=0 & \\\\
 * c_1 & a_{10} \\\\
 * c_2 & a_{20} & a_{21} \\\\
 * c_3 & b_0 & a_{31} & a_{32} \\\\
 * c_4 & b_0 & b_1 & a_{42} & a_{43} \\\\
 * \hline
 * & b_0 & b_1 & b_2 & b_3 & b_4
 * \end{array}
 * $$
 * Notice that the indexing of coefficients starts from 0, unlike in Kennedy et al (2000) where it
 * starts from 1.
 *
 * The update steps thus are
 * $$
 * \begin{align}
 * Q^{(0)} &= Q^{(n)} \\\\
 * Q^{(1)} &= Q^{(0)} + \Delta t a_{10} R^{(0)} \\\\
 * Q^{(2)} &= Q^{(0)} + \Delta t \left[ a_{20}R^{(0)} + a_{21} R^{(1)} \right] \\\\
 *         &= Q^{(1)} + \Delta t \left[ (a_{20}-a_{10})R^{(0)} + a_{21} R^{(1)} \right] \\\\
 * Q^{(3)} &= Q^{(2)} + \Delta t \left[
 *              (b_0-a_{20})R^{(0)} + (a_{31}-a_{21})R^{(1)} + a_{32}R^{(2)}
 *            \right] \\\\
 * Q^{(4)} &= Q^{(3)} + \Delta t \left[
 *              (b_1-a_{31})R^{(1)} + (a_{42}-a_{32})R^{(2)} + a_{43}R^{(3)}
 *            \right] \\\\
 * Q^{(5)} &= Q^{(4)} + \Delta t \left[
 *              (b_2-a_{42})R^{(2)} + (b_3-a_{43})R^{(3)} + b_4 R^{(4)}
 *            \right] \\\\
 * Q^{(n+1)} &= Q^{(5)}
 * \end{align}
 * $$
 *
 * The coefficients are stored as
 * - `b` (size 5) such that `b[i]` gives @f$b_i,\ i=0,\ldots,4@f$
 * - `a_outer` (size 5) such that `a_outer[i]` gives @f$a_{(i+1)i},\ i=0,\ldots,3@f$ and
 *   `a_outer[4]` gives @f$a_{54}:=b_4@f$
 * - `a_inner` (size 4) such that `a_inner[i]` gives @f$a_{(i+2)i},\ i=0,1,2@f$ and `a_inner[3]`
 *   gives @f$a_{53}:=b_3@f$
 *
 * The introduction of @f$a_{54}@f$ and @f$a_{53}@f$ is to facilitate a loop over stages 3-5 which
 * calculate @f$Q^{(3-5)}@f$
 */
class RK4Stage5Register3
{
    public:

    /**
     * The "b" array. See class documentation for details
     */
    std::array<double, 5> b;

    /**
     * Outer line of "a" array. See class documentation for details
     */
    std::array<double, 5> a_outer;

    /**
     * Inner line of "a" array. See class documentation for details
     */
    std::array<double, 4> a_inner;

    RK4Stage5Register3();
};

#endif
