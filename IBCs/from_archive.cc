/**
 * @file from_archive.cc
 * @brief An IC that is set using archives written by solution transfer
 */

#include "from_archive.h"

using namespace ICs;

/**
 * Constructor. Does the tasks mentioned in detailed documentation.
 * 1. Set FromArchive::ar_triang_ along with its manifold using `ar_mesh_filename` and
 *    `mfld_desc_ptr`. Assumes mesh in msh format.
 * 2. Set FromArchive::ar_dof_handler_ using FromArchive::ar_triang_ and `ar_fe_degree`.
 * 3. Set FromArchive::ar_gcvars_ and FromArchive::ar_gh_gcvars_ using all dofs as relevant
 */
FromArchive::FromArchive(
    const DoFHandler<dim> &dh,
    const std::map<psize, Point<dim>> &dl,
    std::array<LA::MPI::Vector, 5> &gcv,
    const MPI_Comm &mpi_comm,
    const std::string &ar_mesh_filename,
    const std::unique_ptr<ManifoldDescriptions::ManifoldDescription> mfld_desc_ptr,
    const std::unique_ptr<MappingQGeneric<dim>> mapping_ptr,
    const usi ar_fe_degree,
    const std::string &ar_filename
):
IC(dh, dl, gcv),
ar_triang_(mpi_comm),
ar_mapping_(*mapping_ptr),
ar_dof_handler_(ar_triang_)
{
    // 1. Read and set manifolds to triangulation
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(ar_triang_);

    std::ifstream ar_mesh_file(ar_mesh_filename);
    AssertThrow(
        ar_mesh_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open mesh file specified for archived solution."
        )
    );
    grid_in.read_msh(ar_mesh_file);
    ar_mesh_file.close();

    if(mfld_desc_ptr != nullptr) mfld_desc_ptr->set(ar_triang_); 

    // 2. Set dof handler
    FE_DGQ<dim> ar_fe(ar_fe_degree);
    ar_dof_handler_.distribute_dofs(ar_fe);

    // 3. Set and read solution vectors, make the ghosted vectors ready-to-use
    IndexSet ar_locally_owned_dofs = ar_dof_handler_.locally_owned_dofs();
    IndexSet ar_all_dofs(ar_dof_handler_.n_dofs());
    ar_all_dofs.add_range(0, ar_dof_handler_.n_dofs());

    std::vector<LA::MPI::Vector*> ar_gcvar_ptrs;
    for(cvar var: cvar_list){
        ar_gcvars_[var].reinit(ar_locally_owned_dofs, mpi_comm);
        ar_gh_gcvars_[var].reinit(ar_locally_owned_dofs, ar_all_dofs, mpi_comm);
        ar_gcvar_ptrs.emplace_back(&ar_gcvars_[var]);
    }

    ar_triang_.load(ar_filename);
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> sol_trans(ar_dof_handler_);
    sol_trans.deserialize(ar_gcvar_ptrs);
    for(cvar var: cvar_list){
        // form ghosted vectors from serial vectors
        ar_gh_gcvars_[var] = ar_gcvars_[var];
    }
}
