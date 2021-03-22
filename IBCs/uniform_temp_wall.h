/**
 * @file uniform_temp_wall.h
 * @brief Boundary condition for a moving wall with spatially uniform prescribed temperature
 */

#ifndef UNIFORM_TEMP_WALL_H
#define UNIFORM_TEMP_WALL_H

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
#include <dgsem/local_dof_data.h>
#include "BC.h"

namespace BCs
{

/**
 * @class UniformTempWall
 * Boundary condition class for a moving wall with spatially uniform prescribed temperature. See
 * the class documentation of BCs::BC
 *
 * The algorithm for the 3 stages is as follows.
 * 1. The ghost conservative state is calculated such that the BR1 auxiliary flux of ghost and inner
 * conservative states gives a conservative state where
 *    - Interface velocity equals wall velocity
 *    - Interface temperature equals prescribed temperature
 *    - Interface density is calculated using prescribed temperature and inner pressure
 * 2. Ghost conservative state is obtained by reversing inner velocity, keeping other variables
 * unchanged.
 * 3. Ghost auxiliary variables equal inner auxiliary variables and ghost velocity (as part of
 * conservative state) is set such that the interface velocity becomes zero (with BR1 viscous flux).
 *
 * Because stage 1 requires calculating density from pressure and temperature, this class requires
 * a NavierStokes object pointer.
 *
 * @note Note the following comments from the class documentation from pens2D:
 * For this BC, Mengaldo et al. (2014) and Bassi & Rebay (1997) disagree only on the way inviscid
 * boundary flux is set. The former uses Weak-Riemann approach while the latter uses Weak-Prescribed
 * approach.
 *
 * Temporal and spatial variation of prescribed temperature and wall velocity is not supported.
 *
 * @warning This class uses a raw pointer to a NavierStokes instance.
 */
class UniformTempWall: public BC
{
    private:
    /**
     * 'Pr'escribed temperature: the temperature which we wish to specify at the wall
     */
    const double T_pr_;
    
    /**
     * The velocity of the wall. This is stored directly by value and not by reference.
     */
    const Tensor<1,dim> vel_pr_;
    
    /**
     * Pointer to a NavierStokes instance. Required for stage 1. This variable is kept private
     * because it is a raw pointer.
     */
    const NavierStokes* ns_ptr_;
    
    public:
    /**
     * Constructor. Calls the base constructor and sets UniformTempWall::T_pr_ and
     * UniformTempWall::vel_pr_.
     */
    UniformTempWall(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const double T_pr,
        const Tensor<1,dim>& vel_pr,
        const NavierStokes* ns_ptr
    ): BC(dh, gcv, gav), T_pr_(T_pr), vel_pr_(vel_pr), ns_ptr_(ns_ptr) {}
    
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

