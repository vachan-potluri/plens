/**
 * @file plens.h
 * @brief The main class of this project.
 */

#ifndef PLENS_H
#define PLENS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/mapping_q_generic.h>

#include "dtype_aliases.h"
#include <utilities/split_string.h>
#include <modelling/navier_stokes.h>

#ifdef DEBUG
#include <utilities/testing.h>
#endif

using namespace dealii;

// see https://stackoverflow.com/questions/4757565/what-are-forward-declarations-in-c
class plens_test; // forward declaration

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
 * read higher order meshes. See the note __pens2D to plens__ for more details. Also see entries
 * around WJ-05-Apr-2021. To truly have curved meshes, i.e.; curvature even in the internal cell
 * edges, manifold must be set to all the internal features too. Because specifying region wise
 * manifold through prm file is difficult, currently only "cylinder flare" and
 * "blunted double cone" geometries are supported. See read_mesh() for more details.
 */
class PLENS
{
    friend plens_test; // for testing

    private:
    /**
     * The dimension
     */
    static constexpr int dim = 3;
    static constexpr usi n_bc_max = 12;

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

    /**
     * The NavierStokes object. Set in set_NS(). A unique pointer is used. When raw pointer is
     * required, use get(). Until set_NS() is called, it remains a nullptr as set in constructor
     */
    std::unique_ptr<NavierStokes> ns_ptr;

    /**
     * Pointer to the mapping object. This will be set in read_mesh() where the mesh type
     * (straight/curved) is used to set the degree of mapping. For 'curved' mesh, the mapping
     * degree is set equal to the value taken through constructor.
     */

    public:
    PLENS();
    ~PLENS();
    void declare_parameters();
    void read_mesh();
    void set_NS();
    void set_dof_handler();



    #ifdef DEBUG
    static void test();
    #endif
};

#endif
