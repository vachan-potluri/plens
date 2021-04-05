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
{}



/**
 * Reads mesh based on the settings in prm file.
 *
 * For straight edged meshes, the procedure is simple.
 */
void PLENS::read_mesh()
{}
