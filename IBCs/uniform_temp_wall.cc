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
    const FaceLocalDoFData &ldd,
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
    const FaceLocalDoFData &ldd,
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



/**
 * Ghost aux vars are set to inner aux vars and ghost cons state is set as in stage 1.
 *
 * @note @p normal is unused
 */
void UniformTempWall::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &ca_gh
) const
{
    State cons_in, cons_pr;
    State& cons_gh = ca_gh.get_state();
    get_state(ldd, cons_in);
    
    double p_in = ns_ptr_->get_p(cons_in);
    double rho_pr = p_in/(ns_ptr_->get_R()*T_pr_);
    Tensor<1,dim> vel; // initialise to 0
    ns_ptr_->prim_to_cons(rho_pr, vel, p_in, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons_in[var];
    
    Avars& av_gh = ca_gh.get_avars();
    get_avars(ldd, av_gh);
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
        State cons({1,1,1,2,13}); // p=4, supersonic
        psize gdof_id = bc_p->get_global_dof_id(ldd);
        for(cvar var: cvar_list) bctd.g_cvars[var][gdof_id] = cons[var];
        
        Avars av;
        CAvars cav(&cons, &av);
        
        for(int i=0; i<3; i++){
            std::cout << "Stage " << i << "\n";
            bc_p->get_ghost_wrappers[i](ldd, normal, cav);
            utilities::print_state(cav.get_state());
            utilities::print_avars(cav.get_avars());
        }
    }
}
#endif

