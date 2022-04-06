/**
 * @file subcell_interpolator.h
 * A class for performing subcell interpolation
 */

#ifndef SUBCELL_INTERPOLATOR_H
#define SUBCELL_INTERPOLATOR_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/dofs/dof_handler.h>

#include <modelling/minmod.h>
#include <modelling/navier_stokes.h>
#include <modelling/state.h>

#include "cell_dof_info.h"
#include "dtype_aliases.h"
#include "LA.h"

#include <array>
#include <cmath>
#include <vector>

using namespace dealii;
using namespace slope_limiters;

/**
 * @class SubcellInterpolator
 * A class for performing subcell interpolation. Fortunately, in Hennemann (2021)'s algorithm,
 * the interpolation has to be done only on reference cell and hence the implementation of this
 * class is highly simplified.
 *
 * This class is intended for use in calculating the low order inviscid residual using a 2nd order
 * Finite Volume algorithm. It takes the solution vectors during construction and calculates
 * the linear slopes when reinit() is called. The getter get_left_right_states() then uses these
 * linear slopes along with the limiter to return the left and right states. The approach taken
 * is exactly how it is done in OpenFOAM. See BTP-2 report, Appendix.
 */
class SubcellInterpolator
{
    public:
    static constexpr int dim = 3;

    /**
     * FE degree
     */
    const usi degree;

    /**
     * Dofs per cell
     */
    const usi dofs_per_cell;

    /**
     * Cell dof infor object
     */
    CellDoFInfo cdi;

    /**
     * 1d LGL weights. These will give the subcell interface locations.
     */
    std::vector<double> w_1d;

    /**
     * 1d LGL nodal locations. These will be use to calculate the slopes.
     */
    std::vector<double> node_loc_1d;

    /**
     * 1d subcell interface locations. These will also be used to calculate the slopes.
     */
    std::vector<double> subcell_face_loc_1d;

    /**
     * Const ref to conervative variable vectors
     */
    const std::array<LA::MPI::Vector, 5>& gcrk_cvars;

    /**
     * Const ref to blender value
     */
    const LA::MPI::Vector& gcrk_alpha;

    /**
     * Const ref to slope limiter object
     */
    const SlopeLimiter& slope_lim;

    /**
     * Pointer to NavierStokes object
     */
    const NavierStokes* ns_ptr;

    /**
     * Conservative variable states. Will be set in reinit(). These are public and can also be used
     * outside the class.
     */
    std::vector<State> cell_states;

    /**
     * Conservative variable linear (non-limited) slope. Access:
     * `slopes[dof][dir][cvar]`.
     * Will be set in reinit(). The word "linear" is used to denote linear interpolation, as opposed
     * to limited linear interpolation.
     */
    std::vector<std::array<State, dim>> cell_cvar_linslopes;

    /**
     * Blender value in the current cell
     */
    double alpha;

    SubcellInterpolator(
        const usi d,
        const std::array<LA::MPI::Vector, 5>& cvar_vecs,
        const LA::MPI::Vector& alpha_vec,
        const SlopeLimiter& s,
        const NavierStokes* ns_p
    );

    void reinit(const DoFHandler<dim>::active_cell_iterator& cell);

    void get_left_right_states(
        const TableIndices<dim>& ti,
        const usi dir,
        State& cl,
        State& cr
    );
};

#endif