/**
 * @file from_archive.h
 * @brief An IC that is set using archives written by solution transfer
 */

#ifndef FROM_ARCHIVE_H
#define FROM_ARCHIVE_H

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
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
#include <memory>

using namespace dealii;
namespace ICs
{

/**
 * @class FromArchive
 * An IC class that is used to set IC from archives saved by solution transfer. This is generally
 * used for transferring solution between different grids. This class will work as long as the
 * current problem's domain lies completely within the archived solution's domain. The archives are
 * written during every data output operation in PLENS (see PLENS::write() and
 * PLENS::do_solution_transfer()).
 *
 * This class does its job partly in the constructor and partly in set(). The job in the
 * constructor is to reproduce the solution saved in the archive. These solution vectors
 * (ar_gcvars_) are later used to create FEFieldFunction objects (in set()), which are finally used
 * to set the IC dof-wise (i.e. set the values in IC::g_cvars). To ensure that this class works
 * irrespective of the current partitioning, each processor reads its own copy of the archived
 * triangulation (ar_triang_ptr_) and generates its own solution vectors (ar_gcvars_). The catch
 * here is that for solution transfer to work, the triangulation which calls `load()` must also be
 * a parallel distributed one and as a result, ar_triang_ptr_ has to point to a distributed
 * triangulation; not a serial triangulation. For this purpose, the global mpi communicator
 * (`MPI_COMM_WORLD`) is split into individual communicators such that each processor gets a local
 * communicator (local_mpi_comm_ptr_) that only contains itself. The archived triangulation and
 * solution vectors are loaded using this local communicator.
 *
 * Once this job is done, the rest is made easy by dealii's Functions::FEFieldFunction. Though in
 * principle, FEFieldFunction requires a ghosted vector, since we know before hand that archived
 * triangulation is completely held by ar_triang_ptr_ (there is no 'partitioning' because a local
 * communicator containing only current process is being used), a non-ghosted vector (e.g. any of
 * ar_gcvars_) does the job as well.
 *
 * @note There are better ways to do the solution evaluation part. See the thread
 * https://groups.google.com/g/dealii/c/p5CE15yChrI and also WJ-04-Aug-2021.
 * Wolfgang and Peter have given some alternatives. However, I don't think any changes are
 * currently required because this class is not used much.
 * @note If the aim is to restart the same simulation (on the same mesh, using same FE and same
 * number of processors), then use ICs::FromArchiveRestart instead. That doesn't involve evaluating
 * a solution and hence is exact. With this class, there are bound to be some distortions if a
 * point belongs to multiple cells (i.e.; lies on their boundaries).
 */
class FromArchive: public IC
{
    private:

    std::unique_ptr<MPI_Comm> local_mpi_comm_ptr_;

    /**
     * A pointer to the triangulation relevant for the solution vectors stored in archive. This can
     * simply be triangulation on which the simulation was done and solution transfer was
     * performed. The triangulation of the actual problem passed through dof handler in the
     * constructor has nothing to do with this. This is constructed using the file
     * `ar_mesh_filename` provided in the constructor. The domain of this must completely encompass
     * the actual problem's triangulation domain (i.e. the triangulation linked to
     * IC::dof_handler).
     */
    std::unique_ptr<parallel::distributed::Triangulation<dim>> ar_triang_ptr_;

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

    public:
    static constexpr int dim = 3;

    const ConditionalOStream pcout;

    FromArchive(
        const DoFHandler<dim> &dh,
        const std::map<psize, Point<dim>> &dl,
        std::array<LA::MPI::Vector, 5> &gcv,
        const MPI_Comm &mpi_comm,
        const std::string &ar_mesh_filename,
        const std::string &ar_mesh_format,
        const std::unique_ptr<ManifoldDescriptions::ManifoldDescription> &ar_mfld_desc_ptr,
        const std::unique_ptr<MappingQGeneric<dim>> &ar_mapping_ptr,
        const usi ar_fe_degree,
        const std::string &ar_filename
    );

    virtual void set() override;
};

} // namespace ICs

#endif
