/**
 * @file periodic.cc
 * @brief Implements periodic BC
 */

#include "periodic.h"

using namespace BCs;

/**
 * Constructor. Calls the base constructor and sets Periodic::per_paris and Periodic::pid. If pid
 * is not 0 or 1, exception is raised.
 */
Periodic::Periodic(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav,
    const std::vector<GridTools::PeriodicFacePair<DoFHandler<dim>::cell_iterator>>& pairs,
    const usi id
): BC(dh, gcv, gav), per_pairs(pairs), pid(id)
{
    AssertThrow(
        pid == 0 || pid == 1,
        StandardExceptions::ExcMessage(
            "The id used for constructing periodic BC must be 0 or 1"
        )
    );
}