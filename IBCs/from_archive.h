/**
 * @file from_archive.h
 * @brief An IC that is set using archives written by solution transfer
 */

#ifndef FROM_ARCHIVE_H
#define FROM_ARCHIVE_H

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include "IC.h"
#include <dgsem/dtype_aliases.h>
#include <dgsem/LA.h>
#include <modelling/var_enums.h>

#include <string>
#include <fstream>
#include <map>

using namespace dealii;
namespace ICs
{

/**
 * @class FromArchive
 * An IC class that is used to set IC from archives saved by solution transfer. This is either used
 * for restarting a simulation, or setting a solution as IC to a probably more finer grid. The
 * algorithm used here is as in pens2D. First, the archive triangulation is read. The archive dof
 * handler is constructed using the constructor parameter provided. Once that is done, the
 * solution vectors stored in the archive are read. Now for these solution vectors, their ghosted
 * versions are created by specifying all dofs as relevant. This way, the ghosted vectors can be
 * used to set the IC of the actual problem regardless of the newer partition.
 *
 * Once the archive triang, dof handler and solution vectors are ready, the only information of the
 * actual problem we require is the dof locations which is provided by the ctor.
 */
class FromArchive: public IC
{
    private:

    /**
     * The triangulation relevant for the solution vectors stored in archive. Preferably, this has
     * to be the exactly same triangulation on which the simulation was done before performing
     * solution transfer. The triangulation of the actual problem passed through dof handler in the
     * constructor has nothing to do with this. This is constructed using the file
     * `ar_mesh_filename` provided in the constructor.
     */
    parallel::distributed::Triangulation<dim> ar_triang_;

    /**
     * DoF handler for archived solution.
     */
    DoFHandler<dim> ar_dof_handler_;

    public:
    static constexpr int dim = 3;

    FromArchive(
        const DoFHandler<dim> &dh,
        const std::map<psize, Point<dim>> &dl,
        std::array<LA::MPI::Vector, 5> &gcv,
        const std::string &ar_filename,
        const std::string &ar_mesh_filename,
        const usi ar_fe_degree
    );
};

} // namespace ICs

#endif
