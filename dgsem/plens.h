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
#include <cstdlib>

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
#include <deal.II/base/function_parser.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
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
#include "subcell_interpolator.h"
#include <utilities/split_string.h>
#include <manifolds/manifold_description.h>
#include <manifolds/cylinder.h>
#include <manifolds/nose_cylinder.h>
#include <modelling/navier_stokes.h>
#include <modelling/var_enums.h>
#include <modelling/minmod.h>
#include <modelling/none.h>
#include <IBCs/IC.h>
#include <IBCs/piecewise_function.h>
#include <IBCs/from_archive.h>
#include <IBCs/from_archive_restart.h>
#include <IBCs/double_mach_reflection.h>
#include <IBCs/BC.h>
#include <IBCs/free.h>
#include <IBCs/outflow.h>
#include <IBCs/uniform_inflow.h>
#include <IBCs/uniform_temp_wall.h>
#include <IBCs/insulated_wall.h>
#include <IBCs/symmetry.h>
#include <IBCs/empty.h>
#include <IBCs/periodic.h>
#include <IBCs/varying_inflow.h>
#include <IBCs/zpg_inflow.h>
#include <utilities/minmod.h>

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
 *
 * @subsection local_time_stepping Local time stepping
 *
 * To accelerate convergence to steady state, local time stepping is useful. Generally, there are
 * different classifications for this.
 *
 * 1. Time accurate
 * 2. Time inaccurate
 *
 * 1. Local update
 * 2. Global update
 *
 * A time accurate local time stepping is generally done for space-time DG methods where solution
 * has high order variation in time as well. This is used to communicate the conservative state
 * data and thus the calculated flux is accurate in time. A time inaccurate stepping doesn't care
 * for this.
 *
 * In a local update, only certain cells which have time deficit are updated. This is generally
 * done for time accurate algorithms. For a global update, all cells are updated based on their
 * individual stable time step.
 *
 * Here, a time inaccurate global update algorithm is used. PLENS::loc_time_steps stores the stable
 * time step for every cell. However, if local time stepping is not requirested for (in prm file)
 * or if the criterion for activation is not met, then all elements of PLENS::loc_time_steps are
 * overr-written with the global time step (PLENS::time_step). The criterion is that steady state
 * error must be smaller than PLENS::local_stepping_threshold times PLENS::time_step. All this is
 * done in calc_time_step().
 *
 * To facilitate using multiple time steps in a domain, the function multiply_time_step_to_rhs()
 * is written. Once this function is called, PLENS::gcrk_rhs gets multiplied dof-wise by
 * PLENS::loc_time_steps.
 *
 * Once local time stepping is activated, PLENS::cur_time loses most of its significance. It
 * becomes the accumulation of smallest time steps over time, which may not correspond to the time
 * in any of the cells. Nevertheless, it is still compared with PLENS::end_time to end the
 * simulation.
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
     * This variable type is for storing variables on the surfaces of a cell. Access:
     * `data[cell-local face id][face-local dof id]`. This is used to define
     * PLENS::locly_ord_surf_term_t which is for storing such data for all cells.
     */
    template <class T>
    using cell_surf_term_t = std::array<
        std::vector<T>,
        n_faces_per_cell
    >;

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
        cell_surf_term_t<T>
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
     * degree is set equal to the value taken through constructor (subsequently stored in
     * PLENS::mapping_ho_degree). Until read_mesh() is called, it remains a nullptr as set in
     * constructor.
     */
    std::unique_ptr<MappingQGeneric<dim>> mapping_ptr;

    /**
     * Pointer to ManifoldDescriptions::ManifoldDescription object. This is set to null in
     * constructor and will be set if the mesh type is specified curved. Set in read_mesh().
     * Normally, this (manifold description) information is not required beyond this function.
     * However, if the initial condition used is ICs::FromArchive, then this pointer is passed so
     * that the triangulation for the solution being read can be assigned correct manifold.
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
     * the relevant dofs correctly. See set_dof_handler() and set_sol_vecs().
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
     * which are updated in every (sub) time step in update().
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
     * update PLENS::gcrk_cvars. Will be calculated in calc_rhs().
     */
    std::array<LA::MPI::Vector, 5> gcrk_rhs;

    /**
     * Like PLENS::gcrk_rhs, but for previous stage RK solution. Used in the last four stages of
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
     * value will be set in PLENS::calc_blender() based on the parameters provided. Generally, this
     * is one of @f$p,\ \rho,\ p\rho@f$.
     */
    LA::MPI::Vector gcrk_blender_var;

    /**
     * An array that stores the value of blender itself (@f$\alpha@f$). The elements of this vector
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
     * Old solution of @f$\rho E@f$. Used for steady state error calculation
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
     * The BlenderCcalculator object
     */
    BlenderCalculator blender_calc;

    /**
     * A function of simulation time which specifies the limit of blender value (PLENS::gcrk_alpha)
     * on wall BCs. This is evaluated and used in calc_blender(). For cells adjacent to wall, the
     * blender value is set to be the maximum of this function value and the value returned by
     * BlenderCalculator::get_blender().
     */
    FunctionParser<dim> wall_blender_limit_function;

    /**
     * The class instance containing RK4 coefficients.
     */
    const RK4Stage5Register3 rk4_coeffs;

    /**
     * RK update order/degree. Takes value 3 or 4. Depending on this value, update() calls
     * update_rk3() and update_rk4() respectively. The value is set in read_time_settings().
     */
    usi rk_order;

    /**
     * Current simulation time. If local time stepping gets activated, this variable would lose its
     * significance. It would have no physical meaning then. But algorithm would still compare this
     * with PLENS::end_time to end the simulation. See @ref local_time_stepping.
     */
    double cur_time;

    /**
     * Simulation end time. Like PLENS::cur_time, this variable would lose significance once local
     * time stepping gets activated. See @ref local_time_stepping.
     */
    double end_time;

    /**
     * Current (global) time step. If local time step has been activated, then this would be the
     * minimum of PLENS::loc_time_steps over all cells. In that case, this variable will lose its
     * relevance. If local time stepping gets activated, this variable would store the minimum time
     * step over all cells. See @ref local_time_stepping.
     */
    double time_step;

    /**
     * Local time steps for all cells. See @ref local_time_stepping. If local stepping has not
     * been asked for, or has not been activated, then all values in this map are set to
     * PLENS::time_step.
     */
    LA::MPI::Vector loc_time_steps;

    /**
     * Ghosted version of PLENS::loc_time_steps
     */
    LA::MPI::Vector gh_loc_time_steps;

    /**
     * The function that dynamically evaluates Courant number (PLENS::Co) based on the simulation
     * time (PLENS::cur_time). This will be properly initialised in PLENS::read_time_settings().
     * This function will be used in PLENS::calc_time_step() to set the Courant number.
     */
    FunctionParser<dim> courant_function;

    /**
     * Courant number for the simulation at the current time. This is obtained by
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
     * A boolean that indicates whether or not local stepping was requested. This does not indicate
     * whether local stepping is actually used in this time step. That would depend on one more
     * criterion, see @ref local_time_stepping.
     */
    bool requested_local_stepping;

    /**
     * The simulation time after which local stepping is switched on if requested.
     */
    double local_stepping_start_time;

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

    /**
     * Cell-averaged values of conservative variables. Used in calc_cell_cons_grad_fv_gl() for
     * calculating conservative gradients in FV sense.
     */
    std::array<LA::MPI::Vector, 5> gcrk_cvar_avg;

    /**
     * Ghosted version of PLENS::gcrk_cvar_avg
     */
    std::array<LA::MPI::Vector, 5> gh_gcrk_cvar_avg;

    /**
     * Auxiliary variables calculated in FV sense. They have a single value for a given cell.
     */
    std::array<LA::MPI::Vector, 9> gcrk_avar_fv;

    /**
     * Ghosted version of PLENS::gcrk_avar_fv.
     */
    std::array<LA::MPI::Vector, 9> gh_gcrk_avar_fv;

    /**
     * A function parser to calculate error in some variable. This will be treated as the exact
     * solution. This is generally useful for doing convergence studies. The error calculation
     * itself is done in . The error variable will be defined in the input file.
     */
    FunctionParser<dim> conv_exact_function;

    /**
     * Pointer to slope limiter object, required for subcell interpolation
     */
    // std::unique_ptr<slope_limiters::SlopeLimiter> slope_lim_ptr;

    /**
     * Pointer to subcell interpolator object
     */
    // std::unique_ptr<SubcellInterpolator> subcell_interp_ptr;

    /**
     * The blender value for diffusive residual scaling. See WJ notes around 08-Jun-2022.
     * $1-\alpha_d$ will be the scaling of the diffusive residual.
     */
    LA::MPI::Vector gcrk_alpha_d;

    /**
     * A function of space and time which, if evaluates to a +ve quantity for a cell center and for
     * the current time, determines if $\alpha_d=\alpha/\alpha_\text{max}$ is to be used. If it
     * evaluates to a negative value, then $\alpha_d=0$ is used.
     */
    FunctionParser<dim> vis_blending_region;



    void form_neighbor_face_matchings(
        IndexSet& loc_rel_dofs,
        const double tol = 1e-8
    );
    void calc_metric_terms();
    void calc_surf_flux(
        const usi stage,
        locly_ord_surf_flux_term_t<double> &surf_flux_term
    ) const;
    // void calc_surf_flux(
    //     const usi stage,
    //     locly_ord_surf_term_t<State> &surf_flux_term
    // ) const;
    void calc_cell_cons_grad(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s1_surf_flux,
        std::vector<std::array<State, 3>>& cons_grad
    ) const;
    void calc_cell_evar_grad(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& evar_surf,
        std::vector<std::array<State, 3>>& evar_grad
    ) const;
    void calc_cell_cons_grad_fv_gl(
        const DoFHandler<dim>::active_cell_iterator& cell,
        std::array<State, 3>& cons_grad
    ) const;
    void calc_cvar_avg();
    void assert_positivity() const;
    void calc_aux_vars();
    void calc_aux_vars_fv();
    void calc_blender(const bool print_wall_blender_limit = false);
    void calc_cell_ho_residual(
        const usi stage,
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_cell_ho_inv_residual(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s2_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_cell_ho_dif_residual(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s3_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_cell_lo_inv_residual(
        const DoFHandler<dim>::active_cell_iterator& cell,
        const locly_ord_surf_flux_term_t<double>& s_surf_flux,
        std::vector<State>& residual
    ) const;
    void calc_cell_dif_residual_fv_gl(
        const DoFHandler<dim>::active_cell_iterator& cell,
        State& residual
    ) const;
    void calc_rhs(
        const bool print_wall_blender_limit = false
    );
    void calc_time_step();
    void multiply_time_step_to_rhs();
    void post_process();
    double calc_ss_error(Vector<double>& cell_ss_error) const;
    double calc_convergence_error();
    void do_solution_transfer(const std::string& filename);
    void write();
    void update_rk4();
    void update_rk3();
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
    void set_blender();
    void read_time_settings();
    void run();



    #ifdef DEBUG
    static void test();
    #endif
};

#endif
