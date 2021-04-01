/**
 * @file piecewise_function.h
 * @brief Sets initial condition based on piecewise defined function.
 */

#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>

#include "IC.h"
#include <dgsem/dtype_aliases.h>
#include <dgsem/LA.h>
#include <utilities/split_string.h>
#include <modelling/var_enums.h>

#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <map>

using namespace dealii;
namespace ICs
{

/**
 * @class PiecewiseFunction
 * @brief Sets initial condition based on piecewise defined function.
 *
 * This class is a direct extension of `piecewise_function` class of pens2D project. This class
 * sets IC based on function defined on pieces. The pieces are cartesian cuboids formed by tensor
 * products of interface locations provided. The interface locations and the functions are to be
 * provided in a file in the following syntax.
@verbatim
<p/c>
<number of x pieces (nx)>
<number of y pieces (ny)>
<number of y pieces (nz)>
<(inline) list of position of x interfaces in ascending order separated by space>
<(inline) list of position of y interfaces in ascending order separated by space>
<(inline) list of position of z interfaces in ascending order separated by space>
<variable 1 function of piece 0>
<variable 2 function of piece 0>
<variable 3 function of piece 0>
<variable 4 function of piece 0>
<variable 5 function of piece 0>
.
.
.
<variable 1 function of piece nx*ny*nz-1>
<variable 2 function of piece nx*ny*nz-1>
<variable 3 function of piece nx*ny*nz-1>
<variable 4 function of piece nx*ny*nz-1>
<variable 5 function of piece nx*ny*nz-1>
@endverbatim
 * The first line must be a character "p" or "c" denoting primitive and conservative respectively.
 * The variables 1-5 will be treated as primitive/conservative variables accoring to this
 * character. The number of pieces must be at least 1 and the number of interfaces provided
 * subsequently must be 1 less than the number of pieces. For instance, if nx=1, then the list of
 * x-interfaces will be blank line.
 *
 * The ordering of pieces follows dealii's convention: index increments fastest in x-direction,
 * followed by y and z directions. An example is shown here for nx=3, ny=2, nz=1.
@verbatim
                         ___ ___ ___
                        |   |   |   |
                        | 3 | 4 | 5 |
                y       |___|___|___|
                ^       |   |   |   |
                |       | 0 | 1 | 2 |
                |       |___|___|___|
                --> x
@endverbatim
 * For the above example, if nz>1, then the index of pieces will continue on the next z level in
 * the same manner.
 *
 * This class sets the IC cell-wise (equivalent to taking `cell_based=true` in piecewise_function
 * class of pens2D project). Meaning, all dofs within a cell are assumed to be in a single piece.
 * This is ensured by using the cell center to determine the piece. The exact location of the
 * remaining dofs doesn't matter. In case the cell center itself lies exactly on an interface, then
 * an exception is thrown. As long as the cell centers are "inside" pieces, the interfaces can
 * occur within cells too. For smooth transitions between pieces, this could be an undesired
 * behaviour. So this function best suits discontinuous piecewise functions. Examples: riemann 2D
 * problems, Shu-Osher test.
 *
 * With all this said, it is however useful if the mesh is constructed in such a way that the
 * interfaces for pieces coincide with some cells' faces. In such special case, even smooth
 * transition between pieces would be well captured.
 */
class PiecewiseFunction: public IC
{
    private:
    /**
     * Number of pieces in each direction
     */
    std::array<usi, dim> np_;

    /**
     * Interface locations in each direction. Access: `ilocs_[dir][interface id]`
     */
    std::array<std::vector<double>, dim> ilocs_;

    /**
     * Pointers to function parser objects. The usage of pointers is to enable resizing the
     * container by using copy and move semantics which are disabled for FunctionParser. See
     * https://groups.google.com/g/dealii/c/gznY81h1Jmw for details. Access:
     * `fpps_[piece id][variable id]`
     */
    std::vector< std::array< std::unique_ptr<FunctionParser<dim>>, 5 > > fpps_

    /**
     * A boolen variable to indicate whether the functions provided in file are primitive
     * variable functions. If @p false, then the functions are conservative variable functions.
     * This value is set based on the character in the first line of provided file. See class
     * documentation for more details.
     */
    bool prim_fns_;

    usi get_piece_id(const Point<dim> &p);

    public:
    static constexpr int dim;
    PiecewiseFunction(
        const DoFHandler<dim> &dh,
        std::array<LA::MPI::Vector, 5> &gcv,
        const std::string &filename
    );
    virtual void set() override;
};

} // namespace ICs

#endif
