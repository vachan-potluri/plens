/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * Constructor.
 */
PLENS::PLENS()
:
mpi_comm(MPI_COMM_WORLD),
pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm)==0)),
triang(mpi_comm)
{
    declare_parameters();
}



/**
 * Destructor.
 */
PLENS::~PLENS()
{}



/**
 * Declares all parameters.
 */
void PLENS::declare_parameters()
{
    pcout << "\nDeclaring parameters\n";

    prm.enter_subsection("mesh");
    {
        prm.declare_entry(
            "type",
            "straight",
            Patterns::Selection("curved|straight"),
            "Type of mesh. Options: curved|straight"
        );

        prm.declare_entry(
            "mesh file format",
            "msh",
            Patterns::Selection("msh|vtk"),
            "Format of mesh file. Only those that permit setting boundary id are of use. "
            "Options: msh|vtk"
        );

        prm.declare_entry(
            "mesh file name",
            "",
            Patterns::Anything(),
            "Mesh file name. No checks are done on the format"
        );
    } // subsection mesh
    prm.leave_subsection();

    std::ofstream sample_file("sample_input_file.prm");
    AssertThrow(
        sample_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file for writing sample input parameters settings."
        )
    );
    prm.print_parameters(sample_file, ParameterHandler::Text);
    sample_file.close();
    pcout << "Sample input file named 'sample_input_file.prm' written\n";
}



/**
 * Reads mesh based on the settings in prm file.
 *
 * For straight edged meshes, the procedure is simple.
 */
void PLENS::read_mesh()
{}



#ifdef DEBUG
void PLENS::test()
{
    utilities::Testing t("PLENS", "class");

    PLENS problem;
}
#endif
