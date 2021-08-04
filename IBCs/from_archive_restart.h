/**
 * @file from_archive_restart
 * @brief A class to set IC from archive, in a special case
 */

#ifndef FROM_ARCHIVE_RESTART_H
#define FROM_ARCHIVE_RESTART_H

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q_generic.h>

#include "IC.h"
#include <dgsem/dtype_aliases.h>
#include <dgsem/LA.h>
#include <modelling/var_enums.h>
#include <manifolds/manifold_description.h>
#include <manifolds/cylinder.h>
#include <manifolds/nose_cylinder.h>

#include <string>
#include <fstream>
#include <algorithm>

using namespace dealii;
namespace ICs{

/**
 * @class FromArchiveRestart
 * @brief A class to set IC from archive, in a special case
 *
 * This class is a simplification of ICs::FromArchive in a special case. For the class
 * ICs::FromArchive to work, the only requirement is that the archive's triangulation must be
 * encompassing the problems' triangulation. For this class, the requirements are more strict:
 * 1. The archive triangulation and the problem's triangulation must be the same. However, a
 *    separate triangulation has to be constructed to be able to call the `load()` function.
 * 2. The finite element of the archive's and problem's dof handler must be the same. Both the
 *    type and the degree have to match. This essentially makes the dof handlers also match.
 * 3. The number of processors used for archive's solution must be exactly same as the number of
 *    processors being used by the current problem. Otherwise, using the loaded solution directly
 *    would not be possible and we would have to evaluate the solution at points which are not
 *    owned (like in ICs::FromArchive).
 *
 * For more details, see https://groups.google.com/g/dealii/c/p5CE15yChrI.
 *
 * The job is done mostly in the constructor which is very similar to ICs::FromArchive's
 * constructor, except that ghosted vectors are not used. The set function then simply equates the
 * problem's solution vectors to the archive's vectors dof-wise.
 */
class FromArchiveRestart: public IC
{
    private:

    /**
     * The triangulation relevant for the solution vectors stored in archive. This must exactly
     * match with the problem's triangulation. See detailed documentation. This is constructed
     * from the problem's dof handler provided in contsructor.
     */
    parallel::distributed::Triangulation<dim> ar_triang_;

    /**
     * DoF handler for archived solution. Again, constructed from the problem's dof handler
     * provided in contsructor.
     */
    DoFHandler<dim> ar_dof_handler_;

    /**
     * Un-ghosted vectors of archived solution.
     */
    std::array<LA::MPI::Vector, 5> ar_gcvars_;

    public:
    static constexpr int dim = 3;

    const ConditionalOStream pcout;

    FromArchiveRestart(
        const DoFHandler<dim> &dh,
        const std::map<psize, Point<dim>> &dl,
        std::array<LA::MPI::Vector, 5> &gcv,
        const std::string &ar_filename
    );
};

} // namespace ICs

#endif
