/**
 * @file empty.h
 * @brief Empty boundary condition
 */

#ifndef EMPTY_H
#define EMPTY_H

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

#endif

namespace BCs{

/**
 * @class Empty
 * @brief Implements an empty BC. This is mainly to be used for transverse boundaries in 1d/2d
 * simulations. As such, this boundary only makes sense when the flow is parallel to the boundary
 * face in these special cases. See also the class documentation of BCs::BC
 *
 * Generally, this viscous algorithm implemented here doesn't matter. The use cases of this BC are
 * mostly for inviscid Riemann problems. However, to be complete, for stages 1 and 3, BCs::Symmetry
 * will be followed.
 *
 * The algorithm for the 3 stages is as follows.
 * 1. The ghost conservative state is set such that BR1 auxiliary flux of ghost and inner state
 * gives a state with
 *    - Pressure and density unchanged
 *    - Tangential velocity to face unchanged
 *    - Normal velocity component becomes 0
 * 2. The ghost conservative state has inner pressure, density and tangential component of inner
 * velocity
 * 3. Ghost avars equal inner avars and ghost cons state is calculated as in stage 1.
 *
 * For stages 1 and 3, the resultant velocity (after taking BR1 flux) becomes
 * @f[
 * \vec{u} = \vec{u}_i - (\vec{u}_i \cdot \vec{n}) \vec{n}
 * @f]
 * For stage 2, the same formula is used to set the ghost state. So the ghost velocity in stage 2
 * differs from BCs::Symmetry by a factor of 2 in the normal component.
 *
 * This class requires NS pointer for all stages.
 *
 * @warning This class uses a raw pointer to a NavierStokes instance.
 *
 * @remark This class was introduced as an experiment. It was found that the results it gives are
 * not good. See WJ-17-Jun-2021.
 */
class Empty: public BC
{
    private:
    const NavierStokes* ns_ptr_;

    public:
    /**
     * @brief Constructor. Calls the base constructor and sets NS pointer.
     */
    Empty(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const NavierStokes* ns_ptr
    ): BC("empty", dh, gcv, gav), ns_ptr_(ns_ptr) {}

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
        CAvars &cav_gh
    ) override;
};

} // namespace BCs
