/**
 * @file uniform_inflow.cc
 * @brief Spatially uniform iutflow boundary condition
 */

#include "uniform_inflow.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo described in class documentation. Useful
 * formula:
 * @f[
 * \vec{Q}_{gh} = 2\vec{Q}_{pr} - \vec{Q}_{in}\\
 * \text{where } \vec{Q} = \{\rho, \rho u, \rho v, \rho w, \rho E \}
 * @f]
 * @note The ghost state thus computed might in itself be a physically invalid. However, that is not
 * of concern because it is going to be passed to the auxiliary flux function for the values on
 * surface. Here, the BR1 flux (average) of ghost and inner states is the prescribed state.
 *
 * @note @p normal is unused
 */
void UniformInflow::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr_[var] - cons_in[var];
}



/**
 * Sets the ghost state equal to prescribed state
 *
 * @note @p normal and @p ldd are unused
 */
void UniformInflow::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    cons_gh = cons_pr_;
}



/**
 * Sets the ghost cavars according to the algo described in class documentation. Ghost auxililary
 * variables equal inner values and ghost conservative state is calculated like in stage 1. With
 * this, the BR1 viscous flux of ghost and inner cavars is equivalent to theoretical viscous flux
 * calculated using inner avars and prescribed conservative state.
 */
void UniformInflow::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    
    State& cons_gh = cav_gh.get_state();
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr_[var] - cons_in[var];
    
    Avars& av_gh = cav_gh.get_avars();
    get_avars(ldd, av_gh);
}

