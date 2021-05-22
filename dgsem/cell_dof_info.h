/**
 * @file cell_dof_info.h
 * Provides maps to go back and forth cell-local dof ordering and tensorial dof ordering.
 */

#ifndef CELL_DOF_INFO_H
#define CELL_DOF_INFO_H

#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>

#include <cmath>

#include "dtype_aliases.h"

using namespace dealii;

/**
 * @class CellDoFInfo
 * This class gives maps to go back and forth cell-local dof ordering and tensorial dof ordering.
 * This switching is commonly required when calculating contribution of volumetric terms. The class
 * takes polynomial degree for construction and provides two maps: cell-local to tensorial and
 * vice-versa.
 */
class CellDoFInfo
{
    public:
    static constexpr usi dim = 3;

    /**
     * Degree of polynomial basis. Of course, assuming that the `FiniteElement` is `FE_DGQ`.
     */
    const usi degree;

    /**
     * Cell-local to tensorial map. Access:
     * `local_to_tensorial[cell-local dof id][dir]`
     * gives the tensor index for direction `dir`. The access can also be done using `operator()`.
     */
    Table<2, usi> local_to_tensorial;

    /**
     * Tensorial indies to cell-local map. Access:
     * `tensorial_to_local[i][j][k]`
     * gives the cell-local dof id corresponding to tensor indices i, j and k. The access can also
     * be done using `operator()`.
     */
    Table<3, usi> tensorial_to_local;



    /**
     * Constructor. Sets the degree and forms the maps
     */
    CellDoFInfo(const usi deg)
    : degree(deg)
    {
        // first set sizes
        const usi n_dofs_per_cell = std::pow(degree+1, dim);
        local_to_tensorial.reinit(TableIndices<2>(n_dofs_per_cell, dim));
        tensorial_to_local.reinit(TableIndices<3>(degree+1, degree+1, degree+1));
    }
};

#endif
