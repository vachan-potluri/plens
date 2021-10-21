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
 *
 * @param[in] dh Dof handler of the problem for which IC is to be set (as in class IC)
 * @param[in] dl Dof locations corresponding to `dh` (as in class IC)
 * @param[out] gcv Conservative variable vectors of the problem which are to be set (as in class
 *                 IC). Will be set in set().
 * @param[in] mpi_comm MPI communicator of the problem, generally taken directly from PLENS
 * @param[in] ar_mesh_filename The filename of archived solution's mesh, `msh` format is assumed
 * @param[in] ar_mfld_desc_ptr Pointer to manifold description object to apply on archived
 *                             triangulation. If this is `nullptr`, nothing is done and dealii's
 *                             default settings apply. Generally taken directly from PLENS and
 *                             matches with that of the actual problem's
 * @param[in] ar_mapping_ptr Pointer to mapping object to be used for archived triangulation.
 *                           Generally taken directly from PLENS and matches with that of actual
 *                           problem's
 * @param[in] ar_fe_degree Archived solution's fe degree
 * @param[in] ar_filename The actual archive filename, the one which is generated by solution
 *                        transfer serialising.
 *
 * @note `std::unique_ptr` can only be passed by reference. See
 * https://stackoverflow.com/questions/30905487/how-can-i-pass-stdunique-ptr-into-a-function
 */
FromArchive::FromArchive(
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
):
IC(dh, dl, gcv),
pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm)==0)),
ar_triang_(mpi_comm),
ar_mapping_(*ar_mapping_ptr),
ar_dof_handler_(ar_triang_)
{
    // 1. Read and set manifolds to triangulation
    pcout << "\nSetting IC from archive\nReading archive mesh and applying manifolds\n";
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(ar_triang_);

    std::ifstream ar_mesh_file(ar_mesh_filename);
    AssertThrow(
        ar_mesh_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open mesh file specified for archived solution."
        )
    );

    // set format
    GridIn<dim>::Format fmt;
    if(ar_mesh_format == "msh") fmt = GridIn<dim>::Format::msh;
    else if(ar_mesh_format == "vtk") fmt = GridIn<dim>::Format::vtk;
    else if(ar_mesh_format == "unv") fmt = GridIn<dim>::Format::unv;
    else{
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "Currently ICs::FromArchive can only read archive meshes in these formats: "
                "'msh|vtk|unv'."
            )
        );
    }
    grid_in.read(ar_mesh_file, fmt);
    ar_mesh_file.close();

    if(ar_mfld_desc_ptr != nullptr) ar_mfld_desc_ptr->set(ar_triang_); 

    // 2. Set dof handler
    pcout << "Forming dof handler for archive\n";
    FE_DGQ<dim> ar_fe(ar_fe_degree);
    ar_dof_handler_.distribute_dofs(ar_fe);

    // 3. Set and read solution vectors, make the ghosted vectors ready-to-use
    pcout << "Reading and setting archive solution vectors, this may take some time ...";
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
    pcout << " completed\n";
}



/**
 * Sets FromArchive::g_cvars. The algorithm is simple:
 * - For every conservative variable
 *   - Construct an fe field function for archived solution FromArchive::ar_gh_gcvars_
 *   - Loop over locally owned dofs of the present problem FromArchive::locally_owned_dofs_
 *     - Set the conservative variable value using archived solutions's field function
 */
void FromArchive::set()
{
    pcout << "Setting IC on all owned dofs. This may take long time\n";
    for(cvar var: cvar_list){
        Functions::FEFieldFunction<dim, DoFHandler<dim>, LA::MPI::Vector> ar_cvar_fn(
            ar_dof_handler_,
            ar_gh_gcvars_[var],
            ar_mapping_
        );

        for(auto cur_dof: locally_owned_dofs_){
            try{
                g_cvars[var][cur_dof] = ar_cvar_fn.value(dof_locations[cur_dof]);
            }
            catch(...){
                AssertThrow(
                    false,
                    StandardExceptions::ExcMessage(
                        "Unable to probe for archived solution at one dof location of the "
                        "problem. This commonly occurs if the domain of the problem has portions "
                        "that lie out of the domain of the archived solution domain."
                    )
                );
            }
        } // loop over owned dofs
    } // loop over cvars
}
