/**
 * @file cell_dof_info.h
 * @brief Provides maps to go back and forth cell-local dof ordering and tensorial dof ordering.
 */

#ifndef CELL_DOF_INFO_H
#define CELL_DOF_INFO_H

#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>

#include <cmath>

#include "dtype_aliases.h"

#ifdef DEBUG
#include <iostream>
#include <utilities/testing.h>
#endif

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

        // local to tensorial
        for(usi i=0; i<n_dofs_per_cell; i++){
            // get the tensorial indices
            std::array<usi, dim> ids = {
                i%(degree+1),
                std::floor( (i%((degree+1)*(degree+1)))/(degree+1) ),
                std::floor(i/((degree+1)*(degree+1)))
            };

            for(usi dir=0; dir<dim; dir++) local_to_tensorial[i][dir] = ids[dir];
        }

        // tensorial to local
        for(usi i=0; i<=degree; i++){
            for(usi j=0; j<=degree; j++){
                for(usi k=0; k<=degree; k++){
                    tensorial_to_local[i][j][k] = i + j*(degree+1) + k*(degree+1)*(degree+1);
                }
            }
        }
    }

    #ifdef DEBUG
    static void test()
    {
        utilities::Testing t("CellDoFInfo", "class");
        t.new_block("Testing maps construction");

        const usi degree = 2;
        const usi n_dofs_per_cell = pow(degree+1,3);

        CellDoFInfo cdi(degree);
        std::cout << "Degree " << degree << " maps:\n";

        std::cout << "Local to tensorial:\n";
        for(usi i=0; i<n_dofs_per_cell; i++){
            std::cout << "\tDof " << i << " : ";
            for(usi dir=0; dir<3; dir++){
                std::cout << cdi.local_to_tensorial[i][dir] << " ";
            }
            std::cout << "\n";
        }

        std::cout << "Tensorial to local:\n";
        for(usi i=0; i<=degree; i++){
            for(usi j=0; j<=degree; j++){
                for(usi k=0; k<=degree; k++){
                    std::cout << "\tTensor indices: " << i << " " << j << " " << k << " : "
                        << cdi.tensorial_to_local[i][j][k] << "\n";
                }
            }
        }
    }
    #endif
};

#endif
