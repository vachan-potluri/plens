/**
 * @file BC.h
 * @brief Base class for BC
 */

#include "BC.h"

using namespace BCs;

/**
 * @brief Constructor. Takes and sets the relevant references.
 */
BC::BC(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav
):
dof_handler(dh),
g_cvars(gcv),
g_avars(gav),
degree(dh.get_fe().degree),
fdi(dh.get_fe().degree)
{}

/**
 * @brief Default destructor.
 */
BC::~BC() = default;

