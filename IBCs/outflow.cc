/**
 * @file outflow.cc
 * @brief Outflow BC
 */

#include "outflow.h"

using namespace BCs;

/**
 * @brief Sets the ghost conservative state to inner conservative state obtained using
 * BC::get_state().
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    get_state(ldd, cons_gh);
}



/**
 * Sets the ghost state based on algo described in the class documentation. Useful formula:
 * @f[
 * (\rho E)_{gh} = (\rho E)_{in} + 2\frac{p_{pr}-p_{in}}{\gamma-1}
 * @f]
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    double M = ns_ptr_->get_M(cons_in);
    if(M >= 1) cons_gh = cons_in;
    else{
        // subsonic
        cons_gh = cons_in;
        double p_in = ns_ptr_->get_p(cons_in);
        cons_gh[4] += 2*(p_pr_ - p_in)/(ns_ptr_->get_gma()-1);
    }
}



/**
 * @brief Sets the ghost cavars to inner cavars obtained using BC::get_cavars().
 *
 * @note @p normal is unused
 */
void Outflow::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    get_cavars(ldd, cav_gh);
}



#ifdef DEBUG
void Outflow::test()
{
    utilities::Testing t("Outflow", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
    NavierStokes ns("air");
    
    std::unique_ptr<BC> bc_p =
        std::make_unique<Outflow>(bctd.dof_handler, bctd.g_cvars, bctd.g_avars, 1, &ns);
    
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

