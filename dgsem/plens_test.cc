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
    set_IC_test();
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

    problem.ns_ptr->print_modelling_params();
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
    }
}
