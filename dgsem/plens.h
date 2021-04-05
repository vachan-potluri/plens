/**
 * @file plens.h
 * @brief The main class of this project.
 */

#ifndef PLENS_H
#define PLENS_H

#include <iostream>
#include <fstream>
#include <string>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/grid_in.h>

#include "dgsem/dtype_aliases.h"

#ifdef DEBUG
#include "utilities/testing.h"
#endif

using namespace dealii;

/**
 * @class PLENS
 * @brief A class for perfect-gas laminar explicit Navier-Stokes computations.
 *
 * The algorithm to be implemented is discussed in entries after WJ-15-Feb-2021. The relevant
 * papers are
 * 1. [dgsem1] Hennemann, Ramirez & Hindenlang, JCP, 2021.
 * 2. [dgsem2] Gassner, Winters & Kopriva, JCP, 2016.
 * 3. [dgsem3] Fisher & Carpenter, JCP, 2013.
 *
 * The documentation for this class is being written on-the-fly as the code is being built module
 * by module. Many settings required will be done through dealii::ParameterHandler. Explicit
 * documentation for this will not be provided here, but in the comments of sample prm file. The
 * prm file should be named `input.prm`.
 *
 *
 *
 * @section meshing Meshing: forming the Triangulation object
 *
 * We explicitly distinguish between curved and straight-edge triangulations. The latter are easy
 * to read and construct using dealii's built-in functions.
 *
 * Curved meshes are significantly complicated because dealii doesn't provide a support to directly
 * read higher order meshes. See the note __pens2D to plens__ for more details. These are not yet
 * supported.
 */
class PLENS
{
    private:
    /**
     * The dimension
     */
    static constexpr int dim = 3;

    /**
     * The MPI communicator. Set to MPI_COMM_WORLD in the constructor
     */
    MPI_Comm mpi_comm;

    /**
     * Output stream for parallel processes
     */
    ConditionalOStream pcout;

    /**
     * The ParameterHandler object. The parameters are declared in declare_parameters(). They will
     * be parsed in the constructor
     */
    ParameterHandler prm;

    /**
     * The triangulation object. Completely set in read_mesh()
     */
    parallel::distributed::Triangulation<dim> triang;

    public:
    PLENS();
    ~PLENS();
    void declare_parameters();
    void read_mesh();



    #ifdef DEBUG
    static void test();
    #endif
};

#endif
