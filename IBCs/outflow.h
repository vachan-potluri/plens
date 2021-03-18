/**
 * @file outflow.h
 * @brief Outflow boundary condition
 */

#ifndef OUTFLOW_H
#define OUTFLOW_H

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

namespace BCs
{
/**
 * @class Outflow
 * @brief Outflow boundary condition. See the class documentation of BCs::BC
 *
 * This BC assumes that all getters are called with velocity component along normal being positive.
 * This BC behaves like BCs::Free if the outflow is supersonic. If it is subsonic, then the outlet
 * pressure shows up in the ghost values. The algorithm for ghost values of 3 stages is as follows.
 * 1. Ghost conservative state equals inner conservative state
 * 2. For supersonic case, ghost conservative state equals inner conservative state. For subsonic
 * case, ghost values of @f$\rho@f$, @f$\vec{u}@f$ equal inner values, and ghost value of pressure
 * is calculated as @f$p_{gh}=2p_{pr}-p_{in}@f$ where @f$p_{pr}@f$ is the prescribed pressure,
 * passed in ctor. Temperature is calculated using density and pressure.
 * 3. Ghost cavars equal inner cavars
 *
 * For BR1 algorithm of stages 1 and 3, the aforementioned strategy is similar to what Mengaldo et
 * al (2014) suggest. See the documentation of outflow BC in pens2D for more info.
 *
 * Although not generally required, this BC accounts for spatially varying outlet pressure. Time
 * variation is not supported.
 */
class Outflow: public BC
{
    private:
    /**
     * Outlet pressure. Will be termed 'pr'escribed. This is of use only when the outflow is
     * subsonic.
     */
    double p_pr_;
    
    public:
    /**
     * Constructor. Calls the base constructor and sets Outflow::p_pr_.
     */
    Outflow(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const double p_pr
    ): BC(dh, gcv, gav), p_pr_(p_pr) {}
    
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

