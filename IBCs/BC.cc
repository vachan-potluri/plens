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
{
    form_cell_map();
}



/**
 * @brief Default destructor.
 */
BC::~BC() = default;



/**
 * @brief Populates BC::cell_map_
 *
 * The function loops over active cell iterators and adds the cell index to the map if it is owned.
 */
void BC::form_cell_map()
{
    for(auto cell: dof_handler.active_cell_iterators()){
        if(cell->is_locally_owned()) cell_map_[cell->index()] = cell;
    }
}

