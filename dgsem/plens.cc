/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * @brief Constructor.
 */
PLENS::PLENS()
:
mpi_comm(MPI_COMM_WORLD),
triang(mpi_comm)
{}



/**
 * @brief Destructor.
 */
PLENS::~PLENS()
{}



/**
 * Reads mesh based on the settings in prm file.
 *
 * For straight edged meshes, the procedure is simple.
 */
void PLENS::read_mesh()
{}
