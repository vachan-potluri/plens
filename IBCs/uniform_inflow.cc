/**
 * @file uniform_inflow.cc
 * @brief Spatially uniform iutflow boundary condition
 */

#include "uniform_inflow.h"

using namespace BCs;

/**
 * Sets the ghost conservative state according to the algo described in class documentation. Useful
 * formula:
 * @f[
 * \vec{Q}_{gh} = 2\vec{Q}_{pr} - \vec{Q}_{in}\\
 * \text{where } \vec{Q} = \{\rho, \rho u, \rho v, \rho w, \rho E \}
 * @f]
 * @note The ghost state thus computed might in itself be a physically invalid. However, that is not
 * of concern because it is going to be passed to the auxiliary flux function for the values on
 * surface. Here, the BR1 flux (average) of ghost and inner states is the prescribed state.
 *
 * @note @p normal is unused
 */
void UniformInflow::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr_[var] - cons[var];
}



/**
 * Sets the ghost state equal to prescribed state
 *
 * @note @p normal and @p ldd are unused
 */
void UniformInflow::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    cons_gh = cons_pr_;
}



/**
 * Sets the ghost cavars according to the algo described in class documentation. Ghost auxililary
 * variables equal inner values and ghost conservative state is calculated like in stage 1. With
 * this, the BR1 viscous flux of ghost and inner cavars is equivalent to theoretical viscous flux
 * calculated using inner avars and prescribed conservative state.
 *
 * @note The above statement is not true in general, but is true here because the gradients are the
 * same. This makes the diffusive flux (of energy alone, to be precise, momentum fluxes depend only
 * on the gradients) linear in velocity. Thus calculating the flux and taking an average become
 * commutable in this special case.
 */
void UniformInflow::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    const State& cons = cav.get_state();
    State& cons_gh = cav_gh.get_state();
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr_[var] - cons[var];
    
    const Avars& av = cav.get_avars();
    Avars& av_gh = cav_gh.get_avars();
    av_gh = av;
}



#ifdef DEBUG
void UniformInflow::test()
{
    utilities::Testing t("UniformInflow", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
    
    State cons_pr({2,1,4,2,20});
    
    std::unique_ptr<BC> bc_p =
        std::make_unique<UniformInflow>(bctd.dof_handler, bctd.g_cvars, bctd.g_avars, cons_pr);
    
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

