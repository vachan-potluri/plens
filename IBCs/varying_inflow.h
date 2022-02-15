/**
 * @file varying_inflow.h
 * @brief Spatially and temporally varying inflow boundary condition
 */

#ifndef VARYING_INFLOW_H
#define VARYING_INFLOW_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <modelling/var_enums.h>
#include <modelling/navier_stokes.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/face_local_dof_data.h>
#include "BC.h"

#include <string>
#include <map>

namespace BCs{

/**
 * @class VaryingInflow
 * @brief Spatially and temporally varying inflow boundary condition. See the class documentation of
 * BCs::BC.
 *
 * This BC is almost the same as BCs::UniformInflow, except that the prescribed properties can be
 * spatially and temporally varying. Since it can be temporally varying, the simulation time has to
 * be set for this BC for the getters to give accurate values. Since this setting of time is not
 * anything specific to this class, the time setter is inherited from the base class BCs::BC.
 *
 * For the spatially varying part, BCs::BC::get_global_dof_id() is used along with
 * VaryingInflow::dof_locations to evaluate the function parsers.
 */
class VaryingInflow: public BC
{
    private:
    /**
     * Pointer to a NavierStokes instance. Required for stage getting conservative variables from
     * prescribed states.
     */
    const NavierStokes* ns_ptr_;

    /**
     * Prescribed pressure function.
     */
    FunctionParser<dim> p_pr_fn_;

    /**
     * Prescribed temperature function.
     */
    FunctionParser<dim> T_pr_fn_;

    /**
     * Prescribed velocity function. Unlike VaryingInflow::p_pr_fn_ and VaryingInflow::T_pr_fn_,
     * this is a vector valued function.
     */
    FunctionParser<dim> velocity_pr_fn_;

    void calculate_prescribed_state(const FaceLocalDoFData &ldd, State &s);

    public:
    /**
     * Constant reference to dof locations. Set in the constructor.
     */
    const std::map<psize, Point<dim>>& dof_locations;

    VaryingInflow(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const std::map<psize, Point<dim>> &dl,
        const std::string& p_expr,
        const std::string& T_expr,
        const std::string& velocity_expr,
        const NavierStokes* ns_ptr
    );

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
