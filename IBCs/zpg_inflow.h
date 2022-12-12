/**
 * @file zpg_inflow.h
 * @brief Inflow BC with prescribed inflow velocity and zero pressure gradient
 */

#ifndef ZPG_INFLOW_H
#define ZPG_INFLOW_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/tensor.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <modelling/var_enums.h>
#include <modelling/navier_stokes.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/face_local_dof_data.h>
#include "BC.h"

namespace BCs{

/**
 * @class ZPGInflow
 * @brief This class is similar to UniformInflow, just that it now uses inner variables for pressure
 * and density, and specifies the velocity.
 */
class ZPGInflow: public BC
{
    private:
    /**
     * The inflow velocity.
     */
    const Tensor<1,dim> vel_pr_;

    /**
     * Pointer to a NavierStokes instance. Required for stage 1. This variable is kept private
     * because it is a raw pointer.
     */
    const NavierStokes* ns_ptr_;

    public:
    ZPGInflow(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const Tensor<1,dim>& vel_pr,
        const NavierStokes* ns_ptr
    ): BC("ZPG inflow", dh, gcv, gav), vel_pr_(vel_pr), ns_ptr_(ns_ptr) {}

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
};
}
#endif
