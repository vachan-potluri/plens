/**
 * @file cell_dof_info.h
 * Provides maps to go back and forth cell-local dof ordering and tensorial dof ordering.
 */

#ifndef CELL_DOF_INFO_H
#define CELL_DOF_INFO_H

#include "dtype_aliases.h"

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
    /**
     * Degree of polynomial basis. Of course, assuming that the `FiniteElement` is `FE_DGQ`.
     */
    const usi degree;

    /**
     * Cell-local to tensorial map. Access:
     * `local_to_tensorial[cell-local dof id][dir]`
     * gives the tensor index for direction `dir`.
     */
    std::vector<std::array<usi, 3>> local_to_tensorial;

    /**
     * Tensorial indies to cell-local map. Access:
     * `tensorial_to_local[i][j][k]`
     * gives the cell-local dof id corresponding to tensor indices i, j and k.
     */
    std::vector<
        std::vector<
            std::vector<usi>
        >
    > tensorial_to_local;



    /**
     * Constructor. Sets the degree and forms the maps
     */
    CellDoFInfo(const usi deg)
    : degree(deg)
    {}
};

#endif
