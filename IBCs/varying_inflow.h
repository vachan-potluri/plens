/**
 * @file varying_inflow.h
 * @brief Spatially and temporally varying inflow boundary condition
 */

#ifndef VARYING_INFLOW_H
#define VARYING_INFLOW_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/function_parser.h>

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
 * @class VaryingInflow
 * @brief Spatially and temporally varying inflow boundary condition. See the class documentation of
 * BCs::BC.
 *
 * This BC is almost the same as BCs::UniformInflow, except that the prescribed properties can be
 * spatially and temporally varying. Since it can be temporally varying, the simulation time has to
 * be set for this BC for the getters to give accurate values. Since this setting of time is not
 * anything specific to this class, the time setter is inherited from the base class BCs::BC.
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
    FunctionParser<dim> p_pr_;

    /**
     * Prescribed temperature function.
     */
    FunctionParser<dim> T_pr_;

    /**
     * Prescribed velocity function. Unlike VaryingInflow::p_pr_ and VaryingInflow::T_pr_, this is
     * a vector valued function.
     */
    FunctionParser<dim> velocity_pr_;
};

} // namespace BCs

#endif
