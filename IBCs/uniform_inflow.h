/**
 * @file uniform_inflow.h
 * @brief Spatially uniform iutflow boundary condition
 */

#ifndef UNIFORM_INFLOW_H
#define UNIFORM_INFLOW_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <modelling/var_enums.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/local_dof_data.h>
#include "BC.h"

namespace BCs
{

/**
 * @class UniformInflow
 * @brief Spatially uniform inflow boundary condition. See the class documentation of BCs::BC
 *
 * This BC assumes that all getters are called with velocity component along normal being negative.
 * The algorithm for 3 stages is as follows.
 * 1. The ghost conservative state is set such that the BR1 flux of auxiliary variables equals
 * the prescribed conservative state taken in constructor.
 * 2. Ghost conservative state equals prescribed conservative state.
 * 3. Ghost auxiliary variables equal inner auxiliary variables. Ghost conservative state equals
 * that used in stage 1. In this way, the BR1 viscous flux gives the theoretical viscous flux
 * obtained using inner auxiliary variables and prescribed velocity.
 *
 * @note Note the following documentation for this class from pens2D project
 * For this BC, neither of Mengaldo et al. (2014) or Bassi & Rebay (1997) fully describe the
 * implementation. For step 2, it is clear from the former reference that ghost state is to be set
 * to prescribed state. The former doesn't describe for steps 1 and 3, while the latter does.
 * However the former and latter don't exactly use the same wording for how ghost state is to be
 * calculated for step 2. This class may therefore be modified in future. However, since diffusive
 * contribution tends to be low at inflow, the detail for step 1 and 3 should not matter much.
 *
 * Temporal and spatial variation is not supported here.
 */
class UniformInflow: public BC
{
    private:
    /**
     * 'Pr'escribed conservative state: the conservative state specified by us at the inlet. This
     * is stored as a value and not by reference because the state passed to ctor is generally
     * temporarily declared.
     */
    const State cons_pr_;
    
    public:
    /**
     * Constructor. Calls the base constructor and sets UniformInflow::cons_pr_.
     */
    UniformInflow(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const State& cons_pr
    ): BC(dh, gcv, gav), cons_pr_(cons_pr) {}
    
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

} // namespace BCs

#endif

