/**
 * @file double_mach_reflection.h
 * @brief Sets initial condition for the double Mach reflection (DMR) case.
 */

#ifndef DOUBLE_MACH_REFLECTION_H
#define DOUBLE_MACH_REFLECTION_H

#include <deal.II/base/point.h>
#include <deal.II/base/function_parser.h>

#include "IC.h"
#include <dgsem/dtype_aliases.h>
#include <dgsem/LA.h>

using namespace dealii;

namespace ICs{

/**
 * @class DoubleMachReflection
 * @brief Sets initial condition for the double Mach reflection (DMR) case. See Woodward & Colella
 * (1984), Kemm (2016) and Hennemann et al (2021).
 *
 * This class takes 3 inputs:
 * 1. Wedge LE location
 * 2. Wedge angle
 * 3. Offset (from LE) for initial shock position
 *
 * In a conventional setup on rectangular domain with cartesian grid, both the wedge angle and
 * offset are zero, and the wedge LE location is @f$(1/6,0,0)@f$. The offset is taken as a
 * coordinate vector such that the vector addition of wedge LE and the offset gives a point on the
 * initial shock. These parameters are used to form the shock position function object which is
 * evaluated in DoubleMachReflection::set() to set the IC.
 *
 * The initial condition is set here dof-wise using ICs::IC::dof_locations.
 *
 * @note It is assumed that the fluid is calorically perfect gas air.
 */
class DoubleMachReflection: public IC
{
    private:

    /**
     * Analytical initial shock position as a function. If this evaluates to positive (negative)
     * value, it means that the point lies in pre (post) shock flow.
     */
    FunctionParser<dim> shock_pos_function_;

    public:
    DoubleMachReflection(
        const DoFHandler<dim> &dh,
        const std::map<psize, Point<dim>> &dl,
        std::array<LA::MPI::Vector, 5> &gcv,
        const Point<dim>& wedge_le_loc,
        const double wedge_angle,
        const Tensor<1,dim>& shock_offset
    );
    virtual void set() override;
};

}

#endif
