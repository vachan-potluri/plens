/**
 * @file free.cc
 * @brief Free BC
 */

#include "free.h"

using namespace BCs;

/**
 * @brief Sets the ghost conservative state to inner conservative state obtained using
 * BC::get_state().
 *
 * @note @p normal is unused
 */
void Free::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    cons_gh = cons_in;
}



/**
 * @brief Sets the ghost conservative state to inner conservative state obtained using
 * BC::get_state().
 *
 * @note @p normal is unused
 */
void Free::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    State cons_in;
    get_state(ldd, cons_in);
    cons_gh = cons_in;
}



/**
 * @brief Sets the ghost cavars to inner cavars obtained using BC::get_cavars().
 *
 * @note @p normal is unused
 */
void Free::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    State cons;
    Avars av;
    CAvars cav_in(&cons, &av);
    get_cavars(ldd, cav_in);
    cav_gh = cav_in;
}



#ifdef DEBUG
void Free::test()
{
    utilities::Testing t("Free", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
    
    // Free fbc(bctd.dof_handler, bctd.g_cvars, bctd.g_avars);
    std::unique_ptr<BC> fbc_p =
        std::make_unique<Free>(bctd.dof_handler, bctd.g_cvars, bctd.g_avars);
    
    {
        t.new_block("testing ghost getters");
        Tensor<1,dim> normal; // immaterial
        LocalDoFData ldd(1, 1, 3); // cell id, face id, face dof id
        
        State cons;
        Avars av;
        CAvars cav(&cons, &av);
        
        for(int i=0; i<3; i++){
            std::cout << "Stage " << i << "\n";
            fbc_p->get_ghost_wrappers[i](ldd, normal, cav);
            utilities::print_state(cav.get_state());
            utilities::print_avars(cav.get_avars());
        }
    }
}
#endif

