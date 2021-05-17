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
 * velocity is reversed and tangential velocity is retained. Which means the kinetic energy doesn't
 * change. Hence, only components 1-3 of the conservative state need to be changed.
 *
 * @pre @p normal has to be a unit vector
 */
void Symmetry::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    cons_gh[0] = cons[0];
    cons_gh[4] = cons[4];
    
    Tensor<1,dim> smom_in, smom_gh; // specific momentum
    for(int d=0; d<dim; d++) smom_in[d] = cons[1+d];
    double normal_smom = scalar_product(smom_in, normal); // smom_in dot normal
    for(int d=0; d<dim; d++) smom_gh[d] = smom_in[d] - 2*normal_smom*normal[d];
}



/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Symmetry::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons_pr; // cons_pr will be set in this fn
    State& cons_gh = cav_gh.get_state();
    
    Tensor<1,dim> vel_in, vel_pr;
    for(int d=0; d<dim; d++) vel_in[d] = cons[1+d]/cons[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_pr[d] = vel_in[d] - normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons_in);
    ns_ptr_->prim_to_cons(cons[0], vel_pr, p_in, cons_pr); // sets cons_pr
    
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
    
    Avars& av_gh = cav_gh.get_avars();
    av_gh = cav.get_avars();
}



#ifdef DEBUG
void Symmetry::test()
{
    utilities::Testing t("Symmetry", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
    
    NavierStokes ns("air");
    
    std::unique_ptr<BC> bc_p =
        std::make_unique<Symmetry>(bctd.dof_handler, bctd.g_cvars, bctd.g_avars, &ns);
    
    {
        t.new_block("testing ghost getters");
        Tensor<1,dim> normal({0,1,0});
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

