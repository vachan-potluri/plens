/**
 * @file free.h
 * @brief Free boundary condition
 */

#ifndef FREE_H
#define FREE_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/local_dof_data.h>
#include "BC.h"

#include <array>

#ifdef DEBUG
#include <iostream>
#include <memory>
#endif

namespace BCs{

/**
 * @class Free
 * @brief Free boundary condition
 *
 * In all stages, ghost values are set to inner values. The purpose of this BC is to enable running
 * test cases like shock tube and Riemann 2D problems. The boundary condition for such cases are
 * free. This is different from outflow BC where boundary pressure can affect the boundary flux if
 * the outflow is subsonic.
 */
class Free: public BC
{
    public:
    /**
     * @brief Constructor. Nothing special, just calls the base constructor and exits.
     */
    Free(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav
    ): BC(dh, gcv, gav) {}
    
    virtual void get_ghost_stage1(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage2(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage3(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        CAvars &ca_gh
    ) const override;
    
    #ifdef DEBUG
    static void test();
    #endif
};

}

#endif

