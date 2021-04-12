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
    read_mesh_test();
    set_NS_test();
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
    t.new_block("testing construction of NS object");
    PLENS problem;
    problem.set_NS();

    problem.ns_ptr->print_modelling_params();
}
