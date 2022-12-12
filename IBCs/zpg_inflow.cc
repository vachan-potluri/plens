/**
 * @file zpg_inflow.cc
 */

#include "zpg_inflow.h"

using namespace BCs;

void ZPGInflow::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    State cons_pr; // cons_pr will be defined here
    
    const double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
}



void ZPGInflow::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    const double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr_, p_in, cons_gh);
}



void ZPGInflow::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
)
{
    State cons_pr;
    const State& cons = cav.get_state();
    State& cons_gh = cav_gh.get_state();
    const double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];

    const Avars& av = cav.get_avars();
    Avars& av_gh = cav_gh.get_avars();
    av_gh = av;
}
