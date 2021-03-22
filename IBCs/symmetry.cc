/**
 * @file symmetry.cc
 * @brief Symmetry boundary condition
 */

#include "symmetry.h"

using namespace BCs;

/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Symmetry::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in, cons_pr; // cons_pr will be set in this fn
    get_state(ldd, cons_in);
    
    Tensor<1,dim> vel_in, vel_pr;
    for(int d=0; d<dim; d++) vel_in[d] = cons_in[1+d]/cons_in[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_pr[d] = vel_in[d] - normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons_in);
    ns_ptr_->prim_to_cons(cons_in[0], vel_pr, p_in, cons_pr); // sets cons_pr
    
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons_in[var];
}

