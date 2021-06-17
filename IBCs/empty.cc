/**
 * @file empty.cc
 * @brief Empty boundary condition
 */

#include "empty.h"

using namespace BCs;

/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Empty::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_pr; // cons_pr will be set in this fn
    
    Tensor<1,dim> vel_in, vel_pr;
    for(int d=0; d<dim; d++) vel_in[d] = cons[1+d]/cons[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_pr[d] = vel_in[d] - normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr, p_in, cons_pr); // sets cons_pr
    
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
}



/**
 * Sets the ghost state according to the algo described in class documentation. Here, the normal
 * velocity is nulled and tangential velocity is retained. `cons_gh` here is like the prescribed
 * conservative state calculated within stage 1 getter.
 *
 * @pre @p normal has to be a unit vector
 */
void Empty::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    Tensor<1,dim> vel_in, vel_gh;
    for(int d=0; d<dim; d++) vel_in[d] = cons[1+d]/cons[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_gh[d] = vel_in[d] - normal_vel*normal[d];

    double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_gh, p_in, cons_gh); // sets cons_pr
}



/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Empty::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons_pr; // cons_pr will be set in this fn
    State& cons_gh = cav_gh.get_state();
    const State& cons = cav.get_state();
    
    Tensor<1,dim> vel_in, vel_pr;
    for(int d=0; d<dim; d++) vel_in[d] = cons[1+d]/cons[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_pr[d] = vel_in[d] - normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr, p_in, cons_pr); // sets cons_pr
    
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
    
    Avars& av_gh = cav_gh.get_avars();
    const Avars &av = cav.get_avars();
    av_gh = av;
}
