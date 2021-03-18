/**
 * @file outflow.cc
 * @brief Outflow BC
 */

#include "outflow.h"

using namespace BCs;

/**
 * @brief Sets the ghost conservative state to inner conservative state obtained using
 * BC::get_state().
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    cons_gh = cons_in;
}



/**
 * Sets the ghost state based on algo described in the class documentation. Useful formula:
 * @f[
 * (\rho E)_{gh} = (\rho E)_{in} + 2\frac{p_{pr}-p_{in}}{\gamma-1}
 * @f]
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    double M = ns_ptr_->get_M(cons_in);
    if(M >= 1) cons_gh = cons_in;
    else{
        // subsonic
        cons_gh = cons_in;
        double p_in = ns_ptr_->get_p(cons_in);
        cons_gh[4] += 2*(p_pr_ - p_in)/(ns_ptr_->get_gma()-1);
    }
}



/**
 * @brief Sets the ghost cavars to inner cavars obtained using BC::get_cavars().
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons;
    Avars av;
    CAvars cav_in(&cons, &av);
    get_cavars(ldd, cav_in);
    cav_gh = cav_in;
}

