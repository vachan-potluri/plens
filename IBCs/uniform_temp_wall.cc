/**
 * @file uniform_temp_wall.h
 * @brief Boundary condition for a moving wall with spatially uniform prescribed temperature
 */

#include "uniform_temp_wall.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo mentioned in the class documentation. The
 * approach taken is very similar to that used in UniformInflow::get_ghost_stage1()
 *
 * @note @p normal is unused
 */
void UniformTempWall::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in, cons_pr; // cons_pr will be defined here
    get_state(ldd, cons_in);
    
    double p_in = ns_ptr_->get_p(cons_in);
    double rho_pr = p_in/(ns_ptr_->get_R()*T_pr_);
    Tensor<1,dim> vel; // initialise to 0
    ns_ptr_->prim_to_cons(rho_pr, vel, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons_in[var];
}



/**
 * Ghost state has reversed velocity with all other variables unchanged.
 *
 * @note @p normal is unused
 */
void UniformTempWall::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    cons_gh[0] = cons_in[0];
    cons_gh[4] = cons_in[4];
    for(int d=0; d<dim; d++) cons_gh[1+d] = -cons_in[1+d]; // reverse velocity
}

