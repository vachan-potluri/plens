/**
 * @file symmetry.h
 * @brief Symmetry boundary condition
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

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
 * @class Symmetry
 * @brief Implements a symmetric BC. See the class documentation of BCs::BC
 *
 * @note The following note is from the class documentation of pens2D:
 * For this BC, neither Mengaldo et al. (2014) or Bassi & Rebay (1997) have given details. I am
 * following the general logic that face values for step 1 (which are subsequently used in step 3)
 * are the physical values that would occur at a face and the ghost variabels for step 2 correspond
 * to a ficticious cell beyond the boundary
 *
 * @remark pens2D follows weak-Prescribed approach for "steps" ("stages" here) 1 and 3 while here,
 * we follow weak-Riemann approach
 *
 * The algorithm for the 3 stages is as follows.
 * 1. The ghost conservative state is set such that BR1 auxiliary flux of ghost and inner state
 * gives a state with
 *    - Pressure and density unchanged
 *    - Tangential velocity to face unchanged
 *    - Normal velocity component becomes 0
 * 2. The ghost conservative state has inner pressure, density and mirrored inner velocity
 * 3. Ghost avars equal inner avars and ghost cons state is calculated as in stage 1.
 *
 * Mathematically, the ghost velocity for stage 2 satisfies
 * @f[
 * (\vec{u}_g\cdot\vec{n})\vec{n} = -(\vec{u}_i\cdot\vec{n})\vec{n}\\
 * \vec{u}_g - (\vec{u}_g\cdot\vec{n})\vec{n} = \vec{u}_i - (\vec{u}_i\cdot\vec{n})\vec{n}\\
 * \implies \vec{u}_g = \vec{u}_i - 2(\vec{u}_i\cdot\vec{n})\vec{n}
 * @f]
 *
 * For stage 1 & 3, the ghost velocity is just the tangential component of inner velocity
 * @f[
 * \vec{u}_g = \vec{u}_i - (\vec{u}_i\cdot\vec{n})\vec{n}
 * @f]
 *
 * This class requires NS pointer for all stages.
 *
 * @warning This class uses a raw pointer to a NavierStokes instance.
 */
class Symmetry: public BC
{
    private:
    const NavierStokes* ns_ptr_;

    public:
    /**
     * @brief Constructor. Calls the base constructor and sets NS pointer.
     */
    Symmetry(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const NavierStokes* ns_ptr
    ): BC(dh, gcv, gav), ns_ptr_(ns_ptr) {}
    
    virtual void get_ghost_stage1(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage2(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage3(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        CAvars &cav_gh
    ) const override;
    
    #ifdef DEBUG
    static void test();
    #endif
};

} // namespace BCs

#endif

