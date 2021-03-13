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
    const DoFHandler<dim>& dof_handler,
    const std::array<LA::MPI::Vector, 5>& g_cvars,
    const std::array<LA::MPI::Vector, 9>& g_avars
):
dof_handler_(dof_handler),
g_cvars_(g_cvars),
g_avars_(g_avars),
degree_(dof_handler.get_fe().degree),
fdi_(dof_handler.get_fe().degree)
{}

/**
 * @brief Default destructor.
 */
BC::~BC() = default;

