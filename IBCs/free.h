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
#include <dgsem/face_local_dof_data.h>
#include "BC.h"

#include <array>

#ifdef DEBUG
#include <iostream>
#include <memory>
#endif

namespace BCs{

/**
 * @class Free
 * @brief Free boundary condition. See the class documentation of BCs::BC
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
     * @brief Constructor. Just calls the base constructor and exits.
     */
    Free(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav
    ): BC("free", dh, gcv, gav) {}
    
    virtual void get_ghost_stage1(
        const FaceLocalDoFData &ldd,
        const State &cons,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) override;
    
    virtual void get_ghost_stage2(
        const FaceLocalDoFData &ldd,
        const State &cons,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) override;
    
    virtual void get_ghost_stage3(
        const FaceLocalDoFData &ldd,
        const CAvars &cav,
        const Tensor<1,dim> &normal,
        CAvars &ca_gh
    ) override;
    
    #ifdef DEBUG
    static void test();
    #endif
};

} // namespace BCs

#endif

