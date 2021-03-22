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



/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Symmetry::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    
    Tensor<1,dim> vel_in, vel_gh;
    for(int d=0; d<dim; d++) vel_in[d] = cons_in[1+d]/cons_in[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_gh[d] = vel_in[d] - 2*normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons_in);
    ns_ptr_->prim_to_cons(cons_in[0], vel_gh, p_in, cons_gh); // sets cons_gh
}



/**
 * Sets the ghost state according to the algo described in class documentation
 *
 * @pre @p normal has to be a unit vector
 */
void Symmetry::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons_in, cons_pr; // cons_pr will be set in this fn
    State& cons_gh = cav_gh.get_state();
    get_state(ldd, cons_in);
    
    Tensor<1,dim> vel_in, vel_pr;
    for(int d=0; d<dim; d++) vel_in[d] = cons_in[1+d]/cons_in[0];
    double normal_vel = scalar_product(vel_in, normal); // vel_in dot normal
    for(int d=0; d<dim; d++) vel_pr[d] = vel_in[d] - normal_vel*normal[d];
    
    double p_in = ns_ptr_->get_p(cons_in);
    ns_ptr_->prim_to_cons(cons_in[0], vel_pr, p_in, cons_pr); // sets cons_pr
    
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons_in[var];
    
    Avars& av_gh = cav_gh.get_avars();
    get_avars(ldd, av_gh);
}



#ifdef DEBUG
void Symmetry::test()
{
    utilities::Testing t("Symmetry", "class");
    utilities::BCTestData bctd(2,2); // refinement and degree
    
    NavierStokes ns("air");
    
    std::unique_ptr<BC> bc_p =
        std::make_unique<Symmetry>(bctd.dof_handler, bctd.g_cvars, bctd.g_avars, &ns);
    
    {
        t.new_block("testing ghost getters");
        Tensor<1,dim> normal; // immaterial
        LocalDoFData ldd(1, 1, 3); // cell id, face id, face dof id
        
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

