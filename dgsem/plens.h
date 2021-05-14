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
#include <map>
#include <set>

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
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/index_set.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/geometry_info.h>

#include "dtype_aliases.h"
#include <utilities/split_string.h>
#include <modelling/navier_stokes.h>
#include "LA.h"
#include <modelling/var_enums.h>
#include <IBCs/IC.h>
#include <IBCs/piecewise_function.h>
#include <IBCs/BC.h>
#include <IBCs/free.h>
#include <IBCs/outflow.h>
#include <IBCs/uniform_inflow.h>
#include <IBCs/uniform_temp_wall.h>
#include <IBCs/symmetry.h>
#include <IBCs/periodic.h>
#include "face_dof_info.h"

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
 *
 * If there are any periodic boundary conditions, then the mesh has to be in 'standard orientation'
 * as described in dealii documentation. The best example of such meshes is a cartesian mesh.
 *
 * @section dof_handler DoFHandler and solution vectors
 *
 * Because periodic boundary conditions are supported, the dof handler being constructed must take
 * this into account. This is done in set_dof_handler(). Note that there will be two loops required
 * over BCs. Once for setting dof handler and once for setting BCs. These two steps cannot be
 * combined because BCs require dof handler and solution vectors for construction.
 *
 * Once dof handler is constructed, the solution vectors are constructed using the owned and
 * relevant dofs. Relevant dofs here are combination of owned dofs and dofs in neighboring cells
 * and the dofs in cells connected by periodicity (if any). Unlike in pens2D, here all the dofs
 * are taken as relevant rather than just the ones lying on the common face. This may be modified
 * in the future.
 *
 * @section modelling Fluid and flow modelling
 *
 * This is done through the NavierStokes class. The fluid is assumed to be a perfect gas and the
 * flow is assumed laminar. Viscosity is obtained from Sutherland's law and thermal conductivity
 * is obtained using constant Prandtl number.
 *
 * Additionally, the NavierStokes class also provides functions for flux calculations and some
 * helper functions too. The flux functions are extensively used in updating the solution.
 *
 * @section ICs Initial conditions
 *
 * Currently, the code supports only one initial condition: ICs::PiecewiseFunction. See the class
 * documentation for more details. The domain is divided into pieces using cartesian bifurcators
 * and in each piece, a function is set.
 *
 *
 * @section BCs Boundary conditions
 *
 * Currently, the code supports 6 boundary conditions.
 * - BCs::Free
 * - BCs::Outflow
 * - BCs::UniformInflow
 * - BCs::UniformTempWall
 * - BCs::Symmetry
 * - BCs::Periodic
 *
 * See the class documentations for details. The periodic boundary conditions is a very fragile one
 * to handle. For periodic BC to be applied, the mesh must be in standard orientation (see
 * dealii's documentation for the meaning of this). Further, there must be two entries per pair
 * for a periodic BC and each of those entries must be consistent with each other. By consistent,
 * we mean
 * 1. The 'periodic direction' must be given same for both entries
 * 1. The 'other surface boundary id' must be complementary to each other
 * 1. The 'periodic orientation' must be complementary to each other
 *
 * These are absolutely essential and no checks are done on these. If these are not followed, the
 * code may produce unexpected behaviour.
 *
 * @section face_assem Face term assembly
 *
 * The residual for temporal update consists of a surface term and a volume term. The surface term
 * requires calculation of numerical flux. For this purpose, it is required to obtain a matching
 * dof from neighbor for every dof on a face.
 *
 * For a 2D mesh, this matching is straight forward: face-local dofs from owner and neighbor side
 * have same index for matching. However, this is not so for 3D. See
 * https://groups.google.com/g/dealii/c/u8e2mLq3qeQ
 * and plens_test::face_dof_matching_test(). This sort of approach may work for most meshes, but is
 * not guaranteed to work according to dealii. So instead, the approach taken is as suggested by
 * Wolfgang in the above question: loop over neighbor side dofs on a face and see which of them
 * matches. This procedure is employed in .
 */
class PLENS
{
    friend plens_test; // for testing

    public:
    /**
     * The dimension
     */
    static constexpr int dim = 3;

    /**
     * The maximum number of boundary conditions. Use to declare parameters and set BCs
     */
    static constexpr usi n_bc_max = 12;
    
    /**
     * Number of faces per cell
     */
    static const usi n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;

    /**
     * Locally order surface term type. "Locally ordered" is to emphasize that the access must
     * happen through `cell id --> cell-local face id --> face-local dof id`. Commonly, this type
     * is used to store surface fluxes on faces of owned cells.
     *
     * @note It is very easy to encounter seg fault with this type if the size of inner vectors is
     * not set before accessing them.
     */
    template <class T>
    using locly_ord_surf_term_t = std::map<
        psize,
        std::array<
            std::vector<T>,
            n_faces_per_cell
        >
    >;

    private:
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
     * Degree used to set maping_ptr in case read_mesh() detects a curved mesh. 'HO' abbreviates
     * high order here. Thus this variable makes sense only when >1. If it is equal to 1, then it
     * means straight edges will be used even if mesh subsection in prm file mentions curved mesh.
     */
    const usi mapping_ho_degree;

    /**
     * Pointer to the mapping object. This will be set in read_mesh() where the mesh type
     * (straight/curved) is used to set the degree of mapping. For 'curved' mesh, the mapping
     * degree is set equal to the value taken through constructor. Until read_mesh() is called,
     * it remains a nullptr as set in constructor
     */
    std::unique_ptr<MappingQGeneric<dim>> mapping_ptr;

    /**
     * Finite element object (for volume). Set in constructor.
     */
    FE_DGQ<dim> fe;

    /**
     * Finite element object for face. Generally used to get number of dofs per face. Set in
     * constructor.
     */
    FE_FaceQ<dim> fe_face;

    /**
     * Dof handler object. Set in set_dof_handler()
     */
    DoFHandler<dim> dof_handler;

    /**
     * All relevant dof locations. Set in set_dof_handler()
     */
    std::map<unsigned int, Point<dim>> dof_locations;

    /**
     * Locally owned dofs. Set in set_sol_vecs()
     */
    IndexSet locally_owned_dofs;

    /**
     * Locally relevant dofs. Set in set_sol_vecs()
     */
    IndexSet locally_relevant_dofs;

    /**
     * Global conservative variable vectors. These are mostly intended to serve as buffer to
     * gcrk_cvars
     */
    std::array<LA::MPI::Vector, 5> g_cvars;

    /**
     * Global conservative variable vectors of previous time step. These are required in RK
     * updates
     */
    std::array<LA::MPI::Vector, 5> gold_cvars;

    /**
     * Global conservative variable vectors of 'c'urrent 'RK' solution. These are the main vectors
     * which are updated in every time step
     */
    std::array<LA::MPI::Vector, 5> gcrk_cvars;

    /**
     * Ghosted version of gcrk_cvars
     */
    std::array<LA::MPI::Vector, 5> gh_gcrk_cvars;

    /**
     * Like PLENS::gcrk_cvars, but for auxiliary variables
     */
    std::array<LA::MPI::Vector, 9> gcrk_avars;

    /**
     * Like PLENS::gh_gcrk_avars, but for auxiliary variables
     */
    std::array<LA::MPI::Vector, 9> gh_gcrk_avars;

    /**
     * A list of boundary ids of the decomposed mesh held by this process. This list can be empty
     * too. This list is populated by looping over all faces held by this process. This is then
     * used to form the boundary condition objects. This list enables construction of BC object
     * only for those boundaries held by this mpi process.
     */
    std::set<usi> bid_list;

    /**
     * A map of pointers to BC objects for the boundaries held by this mpi process. The keys to
     * this map are from PLENS::bid_list
     */
    std::map<usi, BCs::BC*> bc_list;

    /**
     * FaceDoFInfo object. Used for looping over faces for surface assembly.
     */
    FaceDoFInfo fdi;

    /**
     * A map giving matching neighbor side dof for a dof on a surface. See @ref face_assem for more
     * details. Access:
     * `nei_face_matching_dofs[cell id][face id][face dof id]=j`
     * such that `dof_ids_neighbor[fdi[neighbor face id][j]]` and
     * `dof_ids[fdi[face id][face dof id]]` lie at the same location. Here `neighbor face id` can
     * be obtained using `cell->neighbor_of_neighbor(face id)`. For faces on boundary, the data
     * held by this object is garbage, unused. Before using, its size must be set.
     */
    // std::map<psize, std::array<std::vector<usi>, n_faces_per_cell> > nei_face_matching_dofs;
    locly_ord_surf_term_t<usi> nei_face_matching_dofs;



    void form_neighbor_face_matchings(const double tol = 1e-4);

    public:
    PLENS(const usi mhod = 2, const usi fe_degree = 1);
    ~PLENS();
    void declare_parameters();
    void read_mesh();
    void set_NS();
    void set_dof_handler();
    void set_sol_vecs();
    void set_IC();
    void set_BC();



    #ifdef DEBUG
    static void test();
    #endif
};

#endif
