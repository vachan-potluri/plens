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
#include <algorithm>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "dtype_aliases.h"
#include "set_manifold.h"
#include "LA.h"
#include "face_dof_info.h"
#include "metric_terms.h"
#include "cell_dof_info.h"
#include "change_of_basis_matrix.h"
#include "blender_calculator.h"
#include "rk4_stage5_register3.h"
#include <utilities/split_string.h>
#include <manifolds/manifold_description.h>
#include <manifolds/cylinder.h>
#include <manifolds/nose_cylinder.h>
#include <modelling/navier_stokes.h>
#include <modelling/var_enums.h>
#include <IBCs/IC.h>
#include <IBCs/piecewise_function.h>
#include <IBCs/from_archive.h>
#include <IBCs/from_archive_restart.h>
#include <IBCs/BC.h>
#include <IBCs/free.h>
#include <IBCs/outflow.h>
#include <IBCs/uniform_inflow.h>
#include <IBCs/uniform_temp_wall.h>
#include <IBCs/symmetry.h>
#include <IBCs/empty.h>
#include <IBCs/periodic.h>

#ifdef DEBUG
#include <utilities/printing.h>
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
 * 1. [1] Hennemann, Ramirez & Hindenlang, JCP, 2021.
 * 2. [2] Gassner, Winters & Kopriva, JCP, 2016.
 * 3. [3] Fisher & Carpenter, JCP, 2013.
 *
 * All solver settings required will be done taken dealii::ParameterHandler. Explicit documentation
 * for this will not be provided here, but will be auto-generated in the comments of sample prm
 * file which is print when a simulation is run. The prm file should be named `input.prm`, no other
 * variations to this are possible. In case multiple prm files are to be stored, they can be kept
 * with particular names and input.prm can be symbolically linked to them as required.
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
 * edges, a manifold must be set to all the internal features too. Because specifying region wise
 * manifold through prm file is difficult, currently only "cylinder flare" and
 * "blunted double cone" geometries are supported. See read_mesh() and SetManifold for more
 * details.
 *
 * If there are any periodic boundary conditions, then the mesh has to be in 'standard orientation'
 * as described in dealii documentation. The best example of such meshes is a cartesian mesh.
 * Periodic boundary conditions will be taken care in @ref dof_handler.
 *
 * @section dof_handler DoFHandler and solution vectors
 *
 * Because periodic boundary conditions are supported, the dof handler being constructed must take
 * this into account. This is done in set_dof_handler(). Note that there will be two loops required
 * over BCs. Once for setting dof handler (viz. here) and once for setting BCs, in set_BC(). These
 * two steps cannot be combined because BCs require dof handler and solution vectors for
 * construction. Periodicity relations are added to PLENS::triang in set_dof_handler() using
 * standard dealii's functions.
 *
 * To create the index sets for owned and relevant dofs, two different approaches are used
 * depending on the scenario.
 * - If the simulation doesn't have periodic boundary conditions, then only the dofs of the
 *   adjacent cell, lying on the common face are considered relevant. These are really all that is
 *   required. Such relevant dofs are formed in form_neighbor_face_matchings(). This function
 *   serves another major purpose too, see @ref face_assem.
 * - If the simulation has periodic boundary conditions, then simply dealii's in-built function
 *   is used to get the relevant dofs. In this case, all dofs of adjacent cells are considered
 *   relevant.
 *
 * The vaiable PLENS::has_periodic_bc indicates whether or not the simulation has periodic boundary
 * conditions. Once dof handler is constructed, the solution vectors are constructed using the
 * owned and relevant dofs.
 *
 * @subsection cell_indices A note on cell indices
 *
 * There are many map-like objects in this project that use "cell index" or "cell id" as the key
 * for mapping. In most cases, these ids are based on `cell->index()`. In many places of the
 * documentation, the word "cell id" is encountered and it always means `cell->index()`. Although,
 * keep in mind that its name may be misleading because `DoFCellAccessor` has another function
 * named `id()` and this is not what is being referred to.
 *
 * However, the problem with these indices is that they may not be contiguous within a processor.
 * The cell indices are generally generated by the triangulation without considerations for
 * partitioning. Moreover, once a triangulation is partitioned, each processor have completely
 * different sets of indices, mostly unrelated to each other. A same index might occur in different
 * processor for a different cell. This is not a problem for data that is only relevant locally in
 * a process (e.g. metric terms). However, there are some variables which take only single value
 * in a cell and also require communication. Presently, the only such data object in this project
 * is the blender value PLENS::gcrk_alpha. Hennemann et al (2021)'s algorithm requires some
 * diffusion of @f$\alpha@f$ which obviously needs communication across processes.
 *
 * In such cases, the index `cell->global_active_cell_index()` can be used. For its explanation,
 * see Wolfgang's reply here:
 * https://groups.google.com/g/dealii/c/_lmP3VCLBsw/m/f8dAEgdPAAAJ.
 * This is introduced in version 9.3.0. It provides a contiguous cell indexing so that parallel
 * vectors may be constructed using these indices. Only for PLENS::gcrk_alpha, the index returned
 * by this function is used as the "key" (actually, it is not a map, so using the term "key" is
 * not strictly correct).
 *
 * And finally to close, using `cell->index()` is valid only when
 * - You want to distinguish between only cells owned by a process
 * - There is no mesh refinement
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
 * and in each piece, a function is set. Shortly, an IC to restart a simulation will be added.
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
 * @remark BCs::Empty was introduced on trial basis and found to not work like expected, and hence
 * is not mentioned here.
 *
 * The philosophy of the boundary condition handling is described in BCs::BC. Basically all that is
 * required are the ghost variables. Now what exactly are these ghost variables depends on what is
 * being calculated. Read about the "stages" described in @ref face_assem and also in NavierStokes.
 *
 * See the individual class documentations for algo-related details. The periodic boundary
 * condition is a very fragile one to handle. For periodic BC to be applied, the mesh must be in
 * standard orientation (see dealii's documentation for the meaning of this). Further, there must
 * be two entries per pair for a periodic BC and each of those entries must be consistent with each
 * other. By consistent, we mean
 * 1. The 'periodic direction' must be given same for both entries
 * 2. The 'other surface boundary id' must be complementary to each other
 * 3. The 'periodic orientation' must be complementary to each other
 *
 * See the documentation of parameters to know what these mean. These are absolutely essential and
 * no checks are done on these. If these are not followed, the code may produce unexpected
 * behaviour.
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
 * matches. This procedure is employed in form_neighbor_face_matchings().
 *
 * Once these matchings are available, it is simple to loop over all faces of actively owned cells
 * and calculate the flux, provided distributed vectors are available. The function
 * calc_surf_flux() does this. It takes an argument for the 'stage' of flux computation. The
 * surface fluxes are required in 3 stages. The information about these stages is given in detail
 * in NavierStokes class.
 * 1. Auxiliary variable calculation. Here, surface values of conservative variables are required
 *    to compute @f$\nabla\vec{v}@f$ and @f$\nabla T@f$.
 * 2. Inviscid conservative flux.
 * 3. Diffusive or viscous conservative flux.
 *
 * The functions in NavierStokes class and all classes inherited from BCs::BC have wrappers to the
 * flux getters of these 3 stages so that they can be indexed by a stage variable. The only
 * distinction between these 3 stages is that while inviscid and viscous fluxes are calculated
 * as normal components at a face dof, the auxiliary flux is left as is. This is because different
 * components of @f$\nabla\vec{v}@f$ and @f$\nabla T@f$ require different normal vector components
 * and hence, as such, there is no "normal" component for auxiliary flux.
 *
 * The presence of flux calculation wrappers in NavierStokes and BCs allows for a single function
 * calc_surf_flux() to do the job for all stages. This tremendously improves code maintainability
 * and significantly assists easy development.
 *
 * @section vol_contrib Volumetric contribution
 *
 * This is the most important, most involved part of the code. See Hennemann et al (2021) [1] for
 * the algorithm.
 *
 * Took all relevant notes in short and attached them to WJ-22-Feb-2021. Also, took detailed notes
 * and attached them to WJ-20-May-2021. These notes are also present physically in TW1 book.
 *
 * Explaining the algo here is impossible. However, some very important notes are mentioned here.
 * 1. Like in refs [1-3], dealii also uses @f$[0,1]^3@f$ as the reference cell. Had it been
 *    @f$[-1,1]^3@f$, things would have got slightly complicated.
 * 2. The volumetric contribution is always calculated by transforming the physical cell to
 *    reference space. Thus, every cell's calculation is completely isolated from other cells.
 * 3. The subcell normal vectors obtained in eq. (B.53) of [1]
 *    - May not be unit vectors
 *    - Do not match (in direction) with the physical normals at faces (local indices) 0, 2 and 4.
 *      Both these facts are emphasised sufficiently in MetricTerms::subcell_normals.
 *
 * The metric terms are calculated using the class MetricTerms and stored as a map. Read the class
 * documentation and also that of MetricTerms::reinit() to get an idea of what is being done.
 * Specifically, note the comments made in these documentations about subcell normals. They don't
 * match with cell normals on all faces.
 *
 * Unlike for surface flux calculation, volumetric contribution cannot be unified for conservative
 * and auxiliary variable calculation. This is because auxiliary variables require gradients of
 * conservative variables in all 3 directions. This means a total of 5*3=15 equations need to be
 * set up for these gradients (see NavierStokes class for reference) whereas inviscid/diffusive
 * contribution for conservative variable residual requires only 5 equations. Hence a
 * straight-forward unified approach is not possible. However, the inviscid and diffusive high
 * order residual contributions can be calculated in a unified manner. The inviscid subcell
 * contribution is then added separately.
 *
 * @section final_residual_calc Final residual calculation
 *
 * The final residual calculation happens in calc_rhs(). It acts like an orchestrator for all other
 * functions required to calculate the rhs, i.e.; the residual at every dof. And then, this
 * function is invoked in update() function which does an RK update of the current solution. The
 * function calc_rhs() already takes care of limiting by blending with a low order solution
 * appropriately.
 *
 * @section time_integration Time integration
 *
 * This is done in the update() function. Currently, only 5 stage, 3 storage RK4 methods are
 * supported. There are many variants of this specific method and any of those could be used
 * by changing the coefficients in RK4Stage5Register3 class. While the original plan was to take
 * RK order as a parameter, it was realised later that methods higher than RK3 cannot be put in
 * a generic algorithm.
 *
 * The algorithm currently used is from Kennedy, Carpenter & Lewis (2000).
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
     * is used to store surface fluxes on faces of owned cells. More commonly, this type is wrapped
     * in an array to be used to store conservative variable flux. See locly_ord_surf_flux_term_t
     *
     * @note It is very easy to encounter seg fault with this type if the size of inner vectors is
     * not set before accessing them. The functions defined here which modify such variables also
     * set the size.
     */
    template <class T>
    using locly_ord_surf_term_t = std::map<
        psize,
        std::array<
            std::vector<T>,
            n_faces_per_cell
        >
    >;

    /**
     * A wrapper to locly_ord_surf_term_t for storing conservative flux in local ordering. Access:
     * `locly_ord_flux_term_t[cvar][cell id][cell-local face id][face-local dof id] = type T`
     */
    template <class T>
    using locly_ord_surf_flux_term_t = std::array<locly_ord_surf_term_t<T>, 5>;

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
     * Pointer to ManifoldDescription object. This is set to null in constructor and will be set
     * if the mesh type is specified curved. Set in read_mesh(). Normally, this (manifold
     * description) information is not required beyond this function. However, if the initial
     * condition used is ICs::FromArchive, then this pointer is passed so that the triangulation
     * for the solution being read can be assigned correct manifold.
     */
    std::unique_ptr<ManifoldDescriptions::ManifoldDescription> mfld_desc_ptr;

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
    std::map<psize, Point<dim>> dof_locations;

    /**
     * A boolean variable indicating if the problem has periodic BC(s). This is required to set
     * the relevant dofs correctly.
     */
    bool has_periodic_bc;

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
     * Global conservative variable vectors of 'c'urrent 'RK' solution. These are the main vectors
     * which are updated in every (sub) time step
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
     * The vector containing "right hand side" or residual for all dofs globally. Will be used to
     * update PLENS::gcrk_cvars.
     */
    std::array<LA::MPI::Vector, 5> gcrk_rhs;

    /**
     * Like PLENS::gcrk_rhs, but for previous stage RK solution. Used in the last three stages of
     * RK4 update.
     */
    std::array<LA::MPI::Vector, 5> gprk_rhs;

    /**
     * Like PLENS::gprk_rhs, but for previous-of-previous stage RK solution. Used in the last three
     * stages of RK4 update.
     */
    std::array<LA::MPI::Vector, 5> gpprk_rhs;

    /**
     * Vector holding dof-wise viscosity for current RK solution.
     */
    LA::MPI::Vector gcrk_mu;

    /**
     * Vector holding dof-wise thermal conductivity for current RK solution.
     */
    LA::MPI::Vector gcrk_k;

    /**
     * Variable used for calculating blender (@f$\alpha@f$) value. This need not be ghosted. Its
     * value will be set in PLENS::calc_blender() based on the parameters provided.
     */
    LA::MPI::Vector gcrk_blender_var;

    /**
     * An array that stores the value of blender itself ($\alpha$). The elements of this vector
     * are to be accessed using `CellAccessor::global_active_cell_index()` which is introduced in
     * dealii-9.3.0. The construction of this vector is done using
     * `Triangulation::global_active_cell_index_partitioner()` which returns (pointer to) an object
     * of type `Utilities::MPI::Partitioner`. The functions
     * `Utilities::MPI::Partitioner::locally_owned_range()` and
     * `Utilities::MPI::Partitioner::ghost_indices()` are used to "reinit" this vector and its
     * ghosted version (see gcrk_gh_alpha). See @ref cell_indices.
     */
    LA::MPI::Vector gcrk_alpha;

    /**
     * Ghosted version of PLENS::gcrk_alpha.
     */
    LA::MPI::Vector gh_gcrk_alpha;

    /**
     * Pressure. Used in PLENS::write(), set in PLENS::post_process().
     */
    LA::MPI::Vector gcrk_p;

    /**
     * Temperature. Used in PLENS::write(), set in PLENS::post_process().
     */
    LA::MPI::Vector gcrk_T;

    /**
     * Velocity. Used in PLENS::write(), set in PLENS::post_process().
     */
    std::array<LA::MPI::Vector, dim> gcrk_vel;

    /**
     * Old solution of $\rho E$. Used for steady state error calculation
     */
    LA::MPI::Vector rhoE_old;

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

    /**
     * 1D weights corresponding to LGL quadrature
     */
    std::vector<double> w_1d;

    /**
     * 1D differentiation matrix in reference cell. Note that the algo uses strong form of DG and
     * hence, the components are (from [2])
     * @f[
     * D_{ij} = \frac{\partial l_j}{\partial \xi}(\xi_i)
     * @f]
     * @note This involves gradient wrt reference cell, and not wrt physical cell.
     */
    FullMatrix<double> ref_D_1d;

    /**
     * Q matrix
     * @f[
     * Q = \text{diag}(w_0, \ldots, w_N) D
     * @f]
     */
    FullMatrix<double> ref_Q_1d;

    /**
     * A map of metric terms. The map will be formed in PLENS::calc_metric_terms().
     * @warning The data stored in metric terms is all public and hence modifiable externally.
     * Don't accidentally modify this data.
     */
    std::map<psize, MetricTerms<dim>> metrics;

    /**
     * A CellDoFInfo object. Useful for calculating residuals.
     */
    CellDoFInfo cdi;

    /**
     * The blender calculator object
     */
    BlenderCalculator blender_calc;

    /**
     * The class instance containing RK4 coefficients.
     */
    const RK4Stage5Register3 rk4_coeffs;

    /**
     * Current simulation time.
     */
    double cur_time;

    /**
     * Simulation end time.
     */
    double end_time;

    /**
     * Current time step.
     */
    double time_step;

    /**
     * Courant number for the simulation.
     */
    double Co;

    /**
     * Output counter.
     */
    unsigned int output_counter;

    /**
     * Number of time steps taken to reach current time from start time.
     */
    unsigned int n_time_steps;

    /**
     * The write frequency
     */
    usi write_freq;

    /**
     * A timer for wall time calculation in the simulation. Used in print statements
     */
    Timer clk;

    /**
     * A timer used to time certain important sections of the code.
     */
    TimerOutput timer;

    /**
     * This is used for generating a pvd file containing all the data files and time value for the
     * entire simulation (of course, only when write() is called). See
     * https://www.dealii.org/current/doxygen/deal.II/namespaceDataOutBase.html#a6f1c052ba49fd44cd8e3f35ba871aebd
     */
    std::vector<std::pair<double, std::string>> times_and_names;



    void form_neighbor_face_matchings(
        IndexSet& loc_rel_dofs,
        const double tol = 1e-8
    );
    void calc_metric_terms();
    void calc_surf_flux(
        const usi stage,
        locly_ord_surf_flux_term_t<double> &surf_flux_term
    ) const;
    void calc_cell_cons_grad(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s1_surf_flux,
        std::vector<std::array<State, 3>>& cons_grad
    ) const;
    void assert_positivity() const;
    void calc_aux_vars();
    void calc_blender();
    void calc_cell_ho_residual(
        const usi stage,
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_cell_lo_inv_residual(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_rhs();
    void calc_time_step();
    void post_process();
    double calc_ss_error(Vector<double>& cell_ss_error) const;
    void do_solution_transfer(const std::string& filename);
    void write();
    void update();

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
    void read_time_settings();
    void run();



    #ifdef DEBUG
    static void test();
    #endif
};

#endif
