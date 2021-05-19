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
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    cons_gh = cons;
}



/**
 * @brief Sets the ghost conservative state to inner conservative state obtained using
 * BC::get_state().
 *
 * @note @p normal is unused
 */
void Free::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    cons_gh = cons;
}



/**
 * @brief Sets the ghost cavars to inner cavars obtained using BC::get_cavars().
 *
 * @note @p normal is unused
 */
void Free::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    cav_gh = cav;
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
        FaceLocalDoFData ldd(1, 1, 3); // cell id, face id, face dof id
        
        State cons, cons_gh;
        Avars av, av_gh;
        CAvars cav(&cons, &av), cav_gh(&cons_gh, &av_gh);
        fbc_p->get_cavars(ldd, cav); // set cav
        
        for(int i=0; i<3; i++){
            std::cout << "Stage " << i << "\n";
            fbc_p->get_ghost_wrappers[i](ldd, cav, normal, cav_gh);
            utilities::print_state(cav_gh.get_state());
            utilities::print_avars(cav_gh.get_avars());
        }
    }
}
#endif

