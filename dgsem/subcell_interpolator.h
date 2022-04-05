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
#include <modelling/state.h>

#include "cell_dof_info.h"
#include "dtype_aliases.h"
#include "LA.h"

#include <vector>
#include <array>

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
 * the limited slopes when reinit() is called.
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
     * Const ref to ghosted conervative variable vectors
     */
    const std::array<LA::MPI::Vector, 5>& gh_gcrk_cvars;

    /**
     * Const ref to slope limiter object
     */
    const SlopeLimiter& slope_lim;

    /**
     * Conservative variable slope. Access:
     * `slopes[dof][dir][cvar]`
     */
    std::vector<std::array<State, dim>> cell_cvar_slopes;

    SubcellInterpolator(
        const usi d,
        const std::array<LA::MPI::Vector, 5>& vecs,
        const std::array<LA::MPI::Vector, 5>& gh_vecs,
        const SlopeLimiter& s
    );

    void reinit(const DoFHandler<dim>::active_cell_iterator& cell);

    void get_left_right_states(
        const TableIndices<dim>& ti,
        const usi dir,
        State& cons_left,
        State& cons_right
    );
};

#endif