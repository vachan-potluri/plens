/**
 * @file from_archive.h
 * @brief An IC that is set using archives written by solution transfer
 */

#ifndef FROM_ARCHIVE_H
#define FROM_ARCHIVE_H

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/numerics/fe_field_function.h>

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
namespace ICs
{

/**
 * @class FromArchive
 * An IC class that is used to set IC from archives saved by solution transfer. This is either used
 * for restarting a simulation, or setting a solution as IC to a probably more finer grid. The
 * algorithm used here is as in pens2D.
 *
 * This class does its job partly in the constructor and partly in set(). The job in the
 * constructor is to reproduce the solution saved in the archive. This involves the following steps
 * 1. Read the archived triangulation. This need not be the same as the triangulation on which the
 *    IC is being set. However, this class uses FEFieldFunction to set the IC and hence the
 *    archived triangulation must be a geometrical superset of the problem's triangulation. Once
 *    read, manifold is applied to the archived triangulation. This is done using the pointer
 *    passed in the constructor, when it is not a `nullptr`.
 * 2. Form the archived dof handler. This is easy, just use the archived triangulation and the fe
 *    degree.
 * 3. Read the archived solution. To read the archived solution, solution vectors must be
 *    constructed. To allow for difference in partitioning, these solution vectors are constructed
 *    with all dofs as relevant.
 *
 * Once this job is done, the rest is made easy by dealii's Functions::FEFieldFunction. The archive
 * dof handler and solution vectors are used to generate such field functions. These field
 * functions can then be asked for values at any point inside the archive's dof handler domain.
 * Using this, the
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
     * Mapping used for archived triangulation
     */
    MappingQGeneric<dim> ar_mapping_;

    /**
     * DoF handler for archived solution.
     */
    DoFHandler<dim> ar_dof_handler_;

    /**
     * Un-ghosted vectors of archived solution.
     */
    std::array<LA::MPI::Vector, 5> ar_gcvars_;

    /**
     * Ghosted vectors of archived solution. All dofs are specified as relevant for these vectors.
     */
    std::array<LA::MPI::Vector, 5> ar_gh_gcvars_;

    public:
    static constexpr int dim = 3;

    FromArchive(
        const DoFHandler<dim> &dh,
        const std::map<psize, Point<dim>> &dl,
        std::array<LA::MPI::Vector, 5> &gcv,
        const MPI_Comm &mpi_comm,
        const std::string &ar_mesh_filename,
        const std::unique_ptr<ManifoldDescriptions::ManifoldDescription> ar_mfld_desc_ptr,
        const std::unique_ptr<MappingQGeneric<dim>> ar_mapping_ptr,
        const usi ar_fe_degree,
        const std::string &ar_filename
    );

    virtual void set() override;
};

} // namespace ICs

#endif
