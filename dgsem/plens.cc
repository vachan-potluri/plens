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
}



/**
 * Reads mesh based on the settings in prm file.
 *
 * For straight edged meshes, the procedure is simple.
 */
void PLENS::read_mesh()
{}
