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
}



/**
 * Desctructor
 */
plens_test::~plens_test() = default;



/**
 * Tests PLENS::read_mesh(). Also outputs the read mesh in the file "read_mesh_test.eps". Eps is
 * used instead of vtk because it supports Q2 mapping.
 */
void plens_test::read_mesh_test() const
{
    t.new_block("testing read_mesh()");
    PLENS problem;
    problem.read_mesh();

    std::ofstream file("read_mesh_test.gpl");
    AssertThrow(
        file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file 'read_mesh_test.gpl' for outputting the triangulation"
        )
    );
    MappingQGeneric<PLENS::dim> mapping(2);
    GridOut grid_out;
    grid_out.write_gnu_plot(problem.triang, file, &mapping);
    file.close();
    std::cout << "Written the triangulation into 'read_mesh_test.gpl'. Go ahead and check!\n";
}
