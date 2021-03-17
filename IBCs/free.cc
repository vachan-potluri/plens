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

