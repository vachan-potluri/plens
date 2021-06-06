/**
 * @file plens_test.cc
 * @brief A class to test PLENS class.
 */

#include "plens_test.h"

/**
 * Constructor. Sets the testing object
 */
plens_test::plens_test()
:
t("PLENS", "class")
{
    // read_mesh_test();
    // set_NS_test();
    // set_IC_test();
    // collect_periodic_faces_test();
    // set_BC_test();
    // face_dof_matching_test();
    // calc_surf_flux_test();
    // calc_cell_cons_grad_test();
    // calc_aux_vars_test();
    // calc_cell_ho_residual_test();
    // plens_test::mapping_ho_metrics_test();
    // plens_test::calc_cell_lo_inv_residual_test();
    // plens_test::calc_blender_test();
    plens_test::calc_rhs_test();
}



/**
 * Desctructor
 */
plens_test::~plens_test() = default;



/**
 * Tests PLENS::read_mesh(). Also outputs the read mesh in the file "read_mesh_test.vtk". Although
 * only triangulation is to be printed, by default grid out functions don't add additional points
 * (based on manifold) for internal cells. See
 * https://groups.google.com/g/dealii/c/UOIvkNV5va4.
 * So to really see if manifold is applied on internal cells, DataOut has to be used. For this
 * purpose a temporary dof handler is constructed here. VtkFlags are not required to be set here
 * because any solution data is not being written.
 */
void plens_test::read_mesh_test() const
{
    t.new_block("testing read_mesh()");
    PLENS problem;
    problem.read_mesh();

    // construct temporary dof handler
    DoFHandler<PLENS::dim> dof_handler(problem.triang);

    MappingQGeneric<PLENS::dim> mapping(4);
    DataOut<PLENS::dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.build_patches(
        mapping, mapping.get_degree(), DataOut<PLENS::dim>::CurvedCellRegion::curved_inner_cells
    );
    std::ofstream file("read_mesh_test.vtk");
    AssertThrow(
        file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file 'read_mesh_test.vtk' for outputting the triangulation"
        )
    );
    data_out.write_vtk(file);
    file.close();
    std::cout << "Written the triangulation into 'read_mesh_test.vtk'. Go ahead and check!\n";
}



/**
 * Tests the construction of NS object
 */
void plens_test::set_NS_test() const
{
    t.new_block("testing set_NS()");
    PLENS problem;
    problem.set_NS();

    #ifdef DEBUG
    problem.ns_ptr->print_modelling_params();
    #endif
}



/**
 * Tests setting of IC. Run this preferrably in solo and in parallel.
 */
void plens_test::set_IC_test() const
{
    t.new_block("testing set_IC()");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();

    DataOut<PLENS::dim> data_out;
    data_out.attach_dof_handler(problem.dof_handler);
    for(cvar var: cvar_list) data_out.add_data_vector(problem.g_cvars[var], cvar_names[var]);

    Vector<float> subdom(problem.triang.n_active_cells());
    for(float &x: subdom) x = problem.triang.locally_owned_subdomain();
    data_out.add_data_vector(subdom, "Subdomain");

    data_out.build_patches(
        *(problem.mapping_ptr),
        problem.mapping_ho_degree,
        DataOut<PLENS::dim>::CurvedCellRegion::curved_inner_cells
    );

    // individual processor files
    std::ofstream proc_file(
        "set_IC_test" +
        Utilities::int_to_string(Utilities::MPI::this_mpi_process(problem.mpi_comm), 2) +
        ".vtu"
    );
    AssertThrow(
        proc_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open processor file"
        )
    );
    data_out.write_vtu(proc_file);

    // master file
    if(Utilities::MPI::this_mpi_process(problem.mpi_comm) == 0){
        std::vector<std::string> filenames;
        for(psize i=0; i<Utilities::MPI::n_mpi_processes(problem.mpi_comm); i++){
            filenames.emplace_back(
                "set_IC_test" +
                Utilities::int_to_string(i,2) +
                ".vtu"
            );
        }
        std::ofstream master_file("set_IC_test.pvtu");
        AssertThrow(
            master_file.good(),
            StandardExceptions::ExcMessage(
                "Unable to open master file"
            )
        );
        data_out.write_pvtu_record(master_file, filenames);

        std::cout << "Written pvtu and vtu files named 'set_IC_test'.\n";
    }
}



/**
 * See WJ-04-May-2021. This function mimics the function at
 * https://www.dealii.org/current/doxygen/deal.II/grid__tools__dof__handlers_8cc_source.html#l02115
 * to see what is causing the assertion to fail. Apparently, the faces are not at all getting
 * detected.
 */
void plens_test::collect_periodic_faces_test() const
{
    t.new_block("testing collect_periodic_faces() by mimicing");
    PLENS problem(2,2);
    problem.read_mesh();
    for(auto cell: problem.triang.active_cell_iterators()){
        const auto face_1 = cell->face(4);
        const auto face_2 = cell->face(5);
        if(face_1->at_boundary()){
            std::cout << "Found face1 of local id 4 at boundary with bid " << face_1->boundary_id()
                << "\n";
        }
        if(face_2->at_boundary()){
            std::cout << "Found face2 of local id 5 at boundary with bid " << face_2->boundary_id()
                << "\n";
        }
    }
}



void plens_test::set_BC_test() const
{
    t.new_block("testing set_BC()");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();
}



/**
 * See WJ-08-May-2021 and WJ-10-May-2021. This function loops over all internal faces and sees if
 * the ordering of dofs used in FaceDoFInfo would be consistent across cells for 'all' meshes. Test
 * this function with all sorts of wierd meshes. Don't run this in parallel
 */
void plens_test::face_dof_matching_test() const
{
    t.new_block("testing face_dof_matching_test() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_dof_handler();

    const FaceDoFInfo fdi(problem.dof_handler.get_fe().degree);
    const usi n_dofs_per_face = (fdi.degree+1)*(fdi.degree+1);
    const usi n_dofs_per_cell = n_dofs_per_face*(fdi.degree+1);

    std::vector<unsigned int> dof_ids(n_dofs_per_cell);
    std::vector<unsigned int> dof_ids_neighbor(n_dofs_per_cell);
    
    for(const auto &cell: problem.dof_handler.active_cell_iterators()){
        for(usi fid=0; fid<problem.n_faces_per_cell; fid++){
            if(cell->face(fid)->at_boundary()) continue;

            usi fid_neighbor = cell->neighbor_of_neighbor(fid); // face id wrt neighbor
            const auto &neighbor = cell->neighbor(fid);
            std::cout << "Cell id: " << cell->index()
                << "\n\tFace id: " << fid
                << "\n\tNeighbor cell index: " << neighbor->index()
                << "\n\t Face id wrt neighbor: " << fid_neighbor << "\n";
            
            cell->get_dof_indices(dof_ids);
            neighbor->get_dof_indices(dof_ids_neighbor);

            std::cout << "\tLooping over face dofs\n";
            for(usi i=0; i<n_dofs_per_face; i++){
                usi dof_id = dof_ids[fdi.maps[fid].at(i)];
                usi dof_id_neighbor = dof_ids_neighbor[fdi.maps[fid_neighbor].at(i)];
                std::cout << "\t\tLocation wrt cell: " << problem.dof_locations[dof_id]
                    << "\n\t\tLocation wrt neighbor: " << problem.dof_locations[dof_id_neighbor]
                    << "\n";
            } // loop over face dofs
        } // loop over internal faces
    } // loop over active cells
}



/**
 * Simply runs PLENS::calc_surf_flux() for all stages. Just to see that there are no run-time
 * issues. Preferable to run this in parallel.
 */
void plens_test::calc_surf_flux_test() const
{
    t.new_block("testing calc_surf_flux() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    PLENS::locly_ord_surf_flux_term_t<double> surf_flux_term;
    for(usi stage=1; stage<=3; stage++){
        std::cout << "Stage " << stage << "\n";
        problem.calc_surf_flux(stage, surf_flux_term);
    }
}



/**
 * Tests calculation of conervative gradients
 */
void plens_test::calc_cell_cons_grad_test() const
{
    t.new_block("testing calc_cell_cons_grad() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    PLENS::locly_ord_surf_flux_term_t<double> s1_surf_flux;
    problem.calc_surf_flux(1, s1_surf_flux); // stage 1 flux

    auto cell = problem.dof_handler.begin_active();
    for(auto c: problem.dof_handler.active_cell_iterators()){
        if(!(c->at_boundary())){
            cell = c;
            break;
        }
    }
    std::vector<std::array<State, 3>> cons_grad(problem.fe.dofs_per_cell);
    problem.calc_cell_cons_grad(cell, s1_surf_flux, cons_grad);

    std::cout << "Cell " << cell->index() << " conservative gradients:\n";
    for(usi i=0; i<problem.fe.dofs_per_cell; i++){
        std::cout << "\tDoF " << i << "\n";
        for(usi dir=0; dir<PLENS::dim; dir++){
            std::cout << "\t\tDirection " << dir << ": ";
            for(cvar var: cvar_list){
                std::cout << cons_grad[i][dir][var] << " ";
            }
            std::cout << "\n";
        }
    }
}



/**
 * Tests the calculation of auxiliary variables.
 */
void plens_test::calc_aux_vars_test() const
{
    t.new_block("testing calc_aux_vars() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    problem.calc_aux_vars();

    std::vector<psize> dof_ids(problem.fe.dofs_per_cell);
    for(auto cell: problem.dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        cell->get_dof_indices(dof_ids);

        std::cout << "Cell: " << cell->index() << "\n";
        for(usi i=0; i<problem.fe.dofs_per_cell; i++){
            std::cout << "\tDoF: " << i << "\n";
            for(avar var: avar_list){
                std::cout << "\t\tAvar " << avar_names[var] << ": "
                    << problem.gcrk_avars[var][dof_ids[i]] << "\n";
            }
        }
    }
}



/**
 * Tests the high order residual calculation. See WJ-31-May-2021. The residual is basically some
 * form of approximating derivatives of flux terms.
 */
void plens_test::calc_cell_ho_residual_test() const
{
    t.new_block("testing calc_cell_ho_residual() function");
    PLENS problem(2,4);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    // get an internal cell
    auto cell = problem.dof_handler.begin_active();
    for(auto c: problem.dof_handler.active_cell_iterators()){
        if(!c->is_locally_owned()) continue;
        if(!(c->at_boundary())){
            cell = c;
            /*
            const usi direction = 2;
            std::vector<psize> dof_ids(problem.fe.dofs_per_cell), dof_ids_nei(problem.fe.dofs_per_cell);
            for(usi lr_id=0; lr_id<2; lr_id++){
                usi face_id = 2*direction + lr_id;
                std::cout << "Face id: " << face_id
                    << "Neighbor cell index: " << cell->neighbor_index(face_id) << "\n";
                for(usi face_dof_id=0; face_dof_id<problem.fe_face.dofs_per_face; face_dof_id++){
                    const auto &neighbor = cell->neighbor(face_id);
                    usi face_id_nei = cell->neighbor_of_neighbor(face_id);
                    cell->get_dof_indices(dof_ids);
                    neighbor->get_dof_indices(dof_ids_nei);

                    psize gdof_id = dof_ids[problem.fdi.maps[face_id].at(face_dof_id)];
                    usi face_dof_nei = problem.nei_face_matching_dofs.at(cell->index())[face_id][face_dof_id];
                    psize gdof_id_nei = dof_ids_nei[problem.fdi.maps[face_id_nei].at(face_dof_nei)];

                    std::cout << "\tOwner dof: " << gdof_id << "\n"
                        << "\tOwner dof location: " << problem.dof_locations[gdof_id]
                        << "\tNeighbor dof location: " << problem.dof_locations[gdof_id_nei]
                        << "\n";
                    
                    State cons, cons_nei;
                    for(cvar var: cvar_list){
                        cons[var] = problem.gcrk_cvars[var][gdof_id];
                        cons_nei[var] = problem.gh_gcrk_cvars[var][gdof_id_nei];
                    }

                    std::cout << "\tOwner state: ";
                    utilities::print_state(cons);
                    std::cout << "\tNeighbor state: ";
                    utilities::print_state(cons_nei);
                }
            }
            */
            break;
        }
    }

    problem.pcout << "Cell: " << cell->index() << "\n";
    problem.pcout << "Cell center: " << cell->center() << "\n";
    problem.pcout << "Cell diameter: " << cell->diameter() << "\n";
    problem.pcout << "Cell jacobian (at dof 13): "
        << problem.metrics.at(cell->index()).detJ[13] << "\n";
    
    std::vector<psize> dof_ids(problem.fe.dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    // calculate auxiliary variables
    problem.calc_aux_vars();
    problem.pcout << "Auxiliary variables calculated:\n";
    for(usi i=0; i<problem.fe.dofs_per_cell; i++){
        problem.pcout << "\tDoF: " << i << "\n";
        for(avar var: avar_list){
            problem.pcout << "\t\tAvar " << avar_names[var] << ": "
                << problem.gcrk_avars[var][dof_ids[i]] << "\n";
        }
    }

    for(usi stage=2; stage<=3; stage++){
        problem.pcout << "Stage " << stage << "\n";
        PLENS::locly_ord_surf_flux_term_t<double> s_surf_flux;
        problem.calc_surf_flux(stage, s_surf_flux); // stage flux

        std::vector<State> residual(problem.fe.dofs_per_cell);
        problem.calc_cell_ho_residual(stage, cell, s_surf_flux, residual);
        
        for(usi i=0; i<problem.fe.dofs_per_cell; i++){
            problem.pcout << "\tDoF: " << i << "\n";
            for(cvar var: cvar_list){
                problem.pcout << "\t\t" << cvar_names[var] << ": " << residual[i][var] << "\n";
            }
        }
    }
}



/**
 * Tests the metric terms using high order mapping on a single cell. For using this function, set
 * mesh in input.prm to one_cylindrical_cell.msh and type to "curved" and subtype to
 * "cylinder flare" (this will set cylindrical manifold). Also set all BCs to free.
 */
void plens_test::mapping_ho_metrics_test() const
{
    t.new_block("testing metric terms using high order mapping");
    PLENS problem(1,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    const auto& cell = problem.dof_handler.begin_active();
    std::cout << "Face connectors (to match physical and references axes):\n";
    for(usi dir=0; dir<3; dir++){
        std::cout << cell->face(2*dir+1)->center() - cell->face(2*dir)->center() << "\n";
    }

    std::cout << "Contravariant vectors\n";
    for(usi i=0; i<problem.fe.dofs_per_cell; i++){
        std::cout << "\tDoF " << i << ": Jacobian: "
            << problem.metrics.at(0).detJ[i] << "\n";
        for(usi dir=0; dir<3; dir++){
            std::cout << "\t\tDirection " << dir << ": "
                << problem.metrics.at(0).JxContra_vecs[i][dir] << "\n";
        }
    }

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    PLENS::locly_ord_surf_flux_term_t<double> s1_surf_flux;
    problem.calc_surf_flux(1, s1_surf_flux); // stage 1 flux
    std::vector<std::array<State, 3>> cons_grad(problem.fe.dofs_per_cell);
    problem.calc_cell_cons_grad(cell, s1_surf_flux, cons_grad);

    std::vector<psize> dof_ids(problem.fe.dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    std::cout << "Cell " << cell->index() << " conservative gradients:\n";
    for(usi i=0; i<problem.fe.dofs_per_cell; i++){
        std::cout << "\tDoF " << i << " location: "
            << problem.dof_locations[dof_ids[i]] << "\n";
        for(usi dir=0; dir<PLENS::dim; dir++){
            std::cout << "\t\tDirection " << dir << ": ";
            for(cvar var: cvar_list){
                std::cout << cons_grad[i][dir][var] << " ";
            }
            std::cout << "\n";
        }
    }
}



/**
 * Tests the calculation of low order inviscid residual
 */
void plens_test::calc_cell_lo_inv_residual_test() const
{
    t.new_block("testing calc_cell_lo_inv_residual() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    // get an internal cell
    auto cell = problem.dof_handler.begin_active();
    for(auto c: problem.dof_handler.active_cell_iterators()){
        if(!c->is_locally_owned()) continue;
        if(!(c->at_boundary())){
            cell = c;
            break;
        }
    }

    PLENS::locly_ord_surf_flux_term_t<double> s2_surf_flux;
    problem.calc_surf_flux(2, s2_surf_flux); // stage flux

    std::vector<State> ho_residual(problem.fe.dofs_per_cell),
        lo_residual(problem.fe.dofs_per_cell);
    problem.calc_cell_ho_residual(2, cell, s2_surf_flux, ho_residual);
    problem.calc_cell_lo_inv_residual(cell, s2_surf_flux, lo_residual);

    problem.pcout << "Cell: " << cell->index() << "\n";
    problem.pcout << "center: " << cell->center() << "\n";

    problem.pcout << "Face connectors:\n";
    for(usi dir=0; dir<3; dir++){
        problem.pcout << "Direction " << dir << ": "
            << cell->face(2*dir+1)->center() - cell->face(2*dir)->center() << "\n";
    }

    problem.pcout << "Contravariant vectors:\n";
    for(usi dir=0; dir<3; dir++){
        problem.pcout << "\tDirection: " << dir << "\n";
        for(usi i=0; i<problem.fe.dofs_per_cell; i++){
            problem.pcout << "\t\tDoF " << i << " : "
                << problem.metrics.at(cell->index()).JxContra_vecs[i][dir] << "\n";
        }
    }
    problem.pcout << "Subcell normals:\n";
    for(usi dir=0; dir<3; dir++){
        problem.pcout << "\tDirection: " << dir << "\n";

        usi dir1 = (dir+1)%3;
        usi dir2 = (dir+2)%3;
        for(usi id=0; id<=problem.fe.degree+1; id++){
            for(usi id1=0; id1<=problem.fe.degree; id1++){
                for(usi id2=0; id2<=problem.fe.degree; id2++){
                    TableIndices<3> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;

                    problem.pcout << "\t\t" << ti << ": "
                        << problem.metrics.at(cell->index()).subcell_normals[dir](ti) << "\n";
                }
            }
        }
    }
    for(usi i=0; i<problem.fe.dofs_per_cell; i++){
        problem.pcout << "\tDoF: " << i << "\n";
        for(cvar var: cvar_list){
            problem.pcout << "\t\t" << cvar_names[var] << ":\n"
                << "\t\t\tHigh order: " << ho_residual[i][var] << "\n"
                << "\t\t\tLow order: " << lo_residual[i][var] << "\n";
        }
    }
}



/**
 * Tests the calculation of blender
 */
void plens_test::calc_blender_test() const
{
    t.new_block("testing calc_blender() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    problem.calc_blender();

    for(auto cell: problem.dof_handler.active_cell_iterators()){
        if(!cell->is_locally_owned()) continue;

        problem.pcout << "Global cell index: " << cell->global_active_cell_index()
            << ", alpha: " << problem.gcrk_alpha[cell->global_active_cell_index()] << "\n";
    }
}



/**
 * Tests PLENS::calc_rhs()
 */
void plens_test::calc_rhs_test() const
{
    t.new_block("testing calc_rhs() function");
    PLENS problem(2,2);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();

    // set gcrk_cvars to g_cvars
    // since time loop is not started, this is manually done
    for(cvar var: cvar_list){
        // std::cout << "\tVar: " << cvar_names[var] << "\n";
        for(auto i: problem.locally_owned_dofs){
            problem.gcrk_cvars[var][i] = problem.g_cvars[var][i];
            // std::cout << "\t\tDoF " << i << ": " << problem.gcrk_cvars[var][i] << "\n";
        }
    }
    for(cvar var: cvar_list){
        problem.gcrk_cvars[var].compress(VectorOperation::insert);
        problem.gh_gcrk_cvars[var] = problem.gcrk_cvars[var];
    }

    problem.calc_aux_vars();
    problem.calc_blender();
    problem.calc_rhs();

    std::vector<psize> dof_ids(problem.fe.dofs_per_cell);
    for(const auto& cell: problem.dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        cell->get_dof_indices(dof_ids);
        problem.pcout << "Cell " << cell->index() << "\n";
        for(usi i=0; i<problem.fe.dofs_per_cell; i++){
            problem.pcout << "\tDoF: " << i << "\n";
            for(cvar var: cvar_list){
                problem.pcout << "\t\t" << cvar_names[var] << " "
                    << problem.gcrk_rhs[var][dof_ids[i]] << "\n";
            }
        }
    } // loop over owned cells
}
