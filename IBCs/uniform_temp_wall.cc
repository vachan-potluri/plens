/**
 * @file uniform_temp_wall.h
 * @brief Boundary condition for a moving wall with spatially uniform prescribed temperature
 */

#include "uniform_temp_wall.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo mentioned in the class documentation. The
 * approach taken is very similar to that used in UniformInflow::get_ghost_stage1(). See eq. (50)
 * of Mengaldo et al (2014). Note: that equation assumes zero wall velocity.
 *
 * @note @p normal is unused
 */
void UniformTempWall::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_pr; // cons_pr will be defined here
    
    double p_in = ns_ptr_->get_p(cons);
    double rho_pr = p_in/(ns_ptr_->get_R()*T_pr_);
    ns_ptr_->prim_to_cons(rho_pr, vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
}



/**
 * Compressible NS: Ghost state has reversed velocity with all other variables unchanged, see eq.
 * (44) of Mengaldo et al (2014).
 * Euler: only normal component reversed, see eq. (33) of Megaldo et al (2014).
 *
 * @note Eq. (34) of Mengaldo et al. (2014) rightly shows the total energy to be unchanged. That's
 * because reversing an orthogonal component doesn't change the velocity magnitude.
 * @note The equations given by Mengaldo et al (2014) are for stationary wall. To apply here,
 * replace @f$\vec{v}@f$ in the paper by @f$\vec{v}-\vec{v}_w@f$ which gives velocity relative to
 * the wall velocity (UniformTempWall::vel_pr_).
 * @note The wall temperature is irrelevant for this function, and hence also for purely inviscid
 * simulations.
 */
void UniformTempWall::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
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
 * Ghost aux vars are set to inner aux vars and ghost cons state is set as in stage 1. Note here
 * again that since auxiliary variables are the set same, the energy flux becomes linear in
 * velocity. And as a result, taking an average of separately calculated inner and ghost diffusive
 * fluxes gives the same flux as Menglado et al (2014) suggest using weak-prescribed approach in
 * eq. (51). The momentum flux depends solely on auxiliary variables and has nothing to do with
 * wall velocity or temperature. The prescribed wall temperature is irrelevant here, even for
 * energy flux.
 *
 * @note @p normal is unused
 */
void UniformTempWall::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons_pr;
    const State& cons = cav.get_state();
    State& cons_gh = cav_gh.get_state();
    
    double p_in = ns_ptr_->get_p(cons);
    double rho_pr = p_in/(ns_ptr_->get_R()*T_pr_);
    ns_ptr_->prim_to_cons(rho_pr, vel_pr_, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
    
    const Avars& av = cav.get_avars();
    Avars& av_gh = cav_gh.get_avars();
    av_gh = av;
}



#ifdef DEBUG
void UniformTempWall::test()
{
    utilities::Testing t("UniformTempWall", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
    
    NavierStokes ns("air");
    Tensor<1,dim> vel({0.5,0.5,0.5});
    
    std::unique_ptr<BC> bc_p =
        std::make_unique<UniformTempWall>(
            bctd.dof_handler,
            bctd.g_cvars,
            bctd.g_avars,
            1.5,
            vel,
            &ns);
    
    {
        t.new_block("testing ghost getters");
        Tensor<1,dim> normal; // immaterial
        FaceLocalDoFData ldd(1, 1, 3); // cell id, face id, face dof id
        
        // modify bctd.g_cvars to get subsonic/supersonic state
        // State cons({1,1,1,2,53}); // p=20, subsonic
        State cons({1,1,1,2,13}), cons_gh; // p=4, supersonic
        psize gdof_id = bc_p->get_global_dof_id(ldd);
        for(cvar var: cvar_list) bctd.g_cvars[var][gdof_id] = cons[var];
        
        Avars av, av_gh;
        bc_p->get_avars(ldd, av);
        CAvars cav(&cons, &av), cav_gh(&cons_gh, &av_gh);
        
        for(int i=0; i<3; i++){
            std::cout << "Stage " << i << "\n";
            bc_p->get_ghost_wrappers[i](ldd, cav, normal, cav_gh);
            utilities::print_state(cav_gh.get_state());
            utilities::print_avars(cav_gh.get_avars());
        }
    }
}
#endif

