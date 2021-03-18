/**
 * @file uniform_temp_wall.h
 * @brief Boundary condition for a moving wall with spatially uniform prescribed temperature
 */

#include "uniform_inflow.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo mentioned in the class documentation. The
 * approach taken is very similar to that used in UniformInflow::get_ghost_stage1()
 *
 * @note @p normal is unused
 */
void UniformInflow::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in, cons_pr; // cons_pr will be defined here
    get_state(ldd, cons_in);
    
    double p_in = ns_ptr_->get_p(cons_in);
    cons_pr[0] = p_in/(ns_ptr_->get_R()*T_pr_);
    cons_pr[4] = p_in/(ns_ptr_->get-gma()-1); // initialise
    for(int d=0; d<dim; d++){
        cons_pr[1+d] = cons_pr[0]*vel_pr_[d];
        cons_pr[4] += 0.5*cons_pr[0]*vel_pr_[d]*vel_pr_[d];
    }
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons_in[var];
}

