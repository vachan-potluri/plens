/**
 * @file periodic.h
 * @brief Implements periodic BC
 */

#ifndef PERIODIC_H
#define PERIODIC_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/exceptions.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/local_dof_data.h>
#include "BC.h"

#include <array>
#include <vector>

#ifdef DEBUG
#include <iostream>
#include <memory>
#endif

namespace BCs
{

/**
 * @class Periodic
 * @brief Implements periodic BCs. See the documentation of BCs::BC.
 * See the note 'pens2D to plens'. The class additionally takes a vector
 * of PeriodicFacePair for construction. This can be generated using
 * GridTools::collect_periodic_faces(). Further, an integer argument is taken for construction
 * which defines which face in the face pair is the face in consideration for this BC object.
 *
 * For all the 3 stages, the ghost values are set to corresponding values on the periodic face
 * pair.
 *
 * It is mandatory that BC::g_cvars and BC::g_avars are ghosted vectors. This would enable access
 * to periodic cells/faces in case they are not owned by this process. Of course, it is also
 * mandatory that the relevant dofs of BC::g_cvars and BC::g_avars are set appropriately to include
 * these periodic ghosted cells. `add_periodicity()` function can be used for this purpose.
 */
class Periodic: public BC
{
    public:
    /**
     * The vector containing periodic face pairs. See the class documentation.
     */
    const std::vector<GridTools::PeriodicFacePair<DoFHandler<dim>::cell_iterator>>& per_pairs;
    /**
     * The 'p'air id. Must be 0 or 1. This indicates which set of Periodic::per_pairs is in
     * consideration for this BC. See the class documentation.
     */
    const usi pid;

    Periodic(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const std::vector<GridTools::PeriodicFacePair<DoFHandler<dim>::cell_iterator>>& pairs,
        const usi id
    );
};

} // namespace BCs

#endif
