/**
 * @file piecewise_function.h
 * @brief Sets initial condition based on piecewise defined function.
 */

#include "piecewise_function.h"

using namespace ICs;

/**
 * Constructor. Calls the base constructor and parses all the functions in the file @p filename.
 */
PiecewiseFunction::PiecewiseFunction
    const DoFHandler<dim> &dh,
    std::array<LA::MPI::Vector, 5> &gcv,
    const std::string &filename
)
: IC(dh, gcv)
{
    std::ifstream fn_file(filename);
    AssertThrow(
        fn_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file provided for reading IC functions."
        )
    );
}
