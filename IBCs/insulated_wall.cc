/**
 * @file insulated_wall.cc
 * @brief A BC for insulated (moving) wall
 */

#include "insulated_wall.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo mentioned in the class documentation.
 * The approach taken is very similar to UniformTempWall::get_ghost_stage1(), except that here
 * inner temperature is used. See eq (52) of Mengaldo et al (2014).
 *
 * @note @p normal is unused
 */
void InsulatedWall::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    State cons_pr; // cons_pr will be defined here
    
    double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
}



/**
 * Compressible NS: Ghost state has reversed velocity with all other variables unchanged, see eq.
 * (44) of Mengaldo et al (2014).
 * Euler: only normal component reversed, see eq. (33) of Megaldo et al (2014).
 *
 * This function is identical to UniformTempWall::get_ghost_stage2().
 *
 * @note Eq. (34) of Mengaldo et al. (2014) rightly shows the total energy to be unchanged. That's
 * because reversing an orthogonal component doesn't change the velocity magnitude.
 * @note The equations given by Mengaldo et al (2014) are for stationary wall. To apply here,
 * replace @f$\vec{v}@f$ in the paper by @f$\vec{v}-\vec{v}_w@f$ which gives velocity relative to
 * the wall velocity (InsulatedWall::vel_pr_).
 */
void InsulatedWall::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    // total energy doesn't change in both cases the since velocity magnitude remains same
    cons_gh[0] = cons[0];
    cons_gh[4] = cons[4];
    if(ns_ptr_->is_inviscid()){
        Tensor<1,dim> smom_in, smom_in_rel; //  absolute and relative specific momentum
        for(int d=0; d<dim; d++){
            smom_in[d] = cons[1+d];
            smom_in_rel[d] = cons[1+d] - cons[0]*vel_pr_[d];
        }
        double normal_smom_rel = scalar_product(smom_in_rel, normal); // smom_in_rel dot normal
        for(int d=0; d<dim; d++) cons_gh[1+d] = smom_in[d] - 2*normal_smom_rel*normal[d];
    }
    else{
        // reverse relative velocity
        for(int d=0; d<dim; d++) cons_gh[1+d] = 2*cons[0]*vel_pr_[d] - cons[1+d];
    }
}



/**
 * Calculate ghost velocity like in get_ghost_stage1() and auxiliary variables such that shear
 * stresses remain constant and heat flux has reversed normal component. See the note in detailed
 * documentation to see why this works. The documentation of UniformTempWall::get_ghost_stage3()
 * might also be useful.
 */
void InsulatedWall::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
)
{
    State cons_pr;
    const State& cons = cav.get_state();
    State& cons_gh = cav_gh.get_state();
    
    double p_in = ns_ptr_->get_p(cons);
    ns_ptr_->prim_to_cons(cons[0], vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];

    const Avars& av = cav.get_avars();
    Avars& av_gh = cav_gh.get_avars();
    for(usi i=0; i<6; i++) av_gh[i] = av[i]; // shear stresses

    Tensor<1,dim> heat_flux_in;
    for(usi i=6; i<9; i++) heat_flux_in[i] = av[i];
    double normal_heat_flux = scalar_product(heat_flux_in, normal);
    for(usi i=6; i<9; i++) av_gh[i] = heat_flux_in[i-6] - 2*normal_heat_flux*normal[i-6];
}
