/**
 * @file insulated_wall.h
 * @brief A BC for insulated (moving) wall
 */

#ifndef INSULATE_WALL_H
#define INSULATE_WALL_H

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
#include "BC.h"

using namespace BCs;

namespace BCs{

/**
 * @class InsulatedWall
 * Boundary condition for a moving insulated wall. See the class documentation of BCs::BC.
 *
 * The algorithm for the 3 stages is as follows (Mengaldo et al (2014)).
 * 1. Ghost conservative state is set such that BR1 flux of ghost and inner state gives and
 *    interface state with
 *    - Velocity equalling wall velocity
 *    - Temperature and pressure equalling inner values
 * 2. Ghost state has reversed inner velocity with same pressure and density. For Euler equations,
 *    only the normal component is reversed. This part is exactly same as in BCs::UniformTempWall.
 * 3. Ghost shear stresses equal inner shear stresses. Ghost heat flux has reversed wall normal
 *    component. Ghost velocity is set like in stage 1. This gives an effect as described by
 *    Mengaldo, in eqs. (53-54) if BR1 viscous flux is used. See the following note.
 *
 * @note Some notes on stage 3 flux. Since the ghost shear stresses are set equal to inner values,
 * the BR1 viscous flux (which averages diffusive flux based on inner and ghost conservative and
 * auxiliary variables) is constant for momentum (i.e.; equals the flux based on inner stresses
 * alone), and the flux for energy becomes linear in velocity and heat flux. Hence, taking the
 * average of diffusive fluxes amounts to taking the average of velocity and heat flux first, and
 * then compute the flux using averaged velocity and heat flux. The linearity makes both these
 * operations commutative. In summary:
 * `diffusive flux(average velocity, average heat flux, inner shear stress)` equals
 * `average(diffusive fluxes of ghost and inner states)`.
 * Mengaldo's paper provides the former in some sense.
 *
 * Because stage 1 requires calculating conservative state from pressure, density and velocity,
 * this class requires a NavierStokes instance.
 *
 * @warning This class uses a raw pointer to a NavierStokes instance.
 */
class InsulatedWall: public BC
{
    private:
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
     * Constructor. Calls the base constructor and sets InsulatedWall::vel_pr_.
     */
    InsulatedWall(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const Tensor<1,dim>& vel_pr,
        const NavierStokes* ns_ptr
    ): BC("insulated wall", dh, gcv, gav), vel_pr_(vel_pr), ns_ptr_(ns_ptr) {}

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

} // namespace BCs

#endif
