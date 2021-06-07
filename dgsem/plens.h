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
#include <deal.II/distributed/tria.h>
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

#include "dtype_aliases.h"
#include "LA.h"
#include "face_dof_info.h"
#include "metric_terms.h"
#include "cell_dof_info.h"
#include "change_of_basis_matrix.h"
#include "blender_calculator.h"
#include <utilities/split_string.h>
#include <modelling/navier_stokes.h>
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
#include "rk_coeffs.h"

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
 * matches. This procedure is employed in form_neighbor_face_matchings().
 *
 * Once this matchings are available, it is simple to loop over all faces of actively owned cells
 * and calculate the flux. The function calc_surf_flux() does this. It takes an argument for the
 * 'stage' of flux computation. The surface fluxes are required in 3 stages. The information about
 * these stages is given in detail in NavierStokes class.
 * 1. Auxiliary variable calculation. Here, surface values of conservative variables are required
 * to compute @f$\nabla\vec{v}@f$ and @f$\nabla T@f$.
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
 * @section vol_contrib Volumetric contribution
 *
 * This is the most important, most involved part of the code. The relevant papers are
 *
 * [1] Hennemann, Ramirez, Hindenlang et al, JCP, 2021.
 *
 * [2] Gassner, Winters & Kopriva, JCP, 2016.
 *
 * [3] Fischer & Carpenter, JCP, 2013.
 *
 * Took all relevant notes in short and attached them to WJ-22-Feb-2021. Also, took detailed notes
 * and attached them to WJ-20-May-2021. These notes are also present physically in TW1 book.
 * Further notes taken will be mentioned as and when done.
 * Explaining the algo here is impossible. However, some very important notes are mentioned here.
 * 1. Like in [1-3], dealii also uses @f$[0,1]^3@f$ as the reference cell. Had it been
 * @f$[-1,1]^3@f$, things would have got slightly complicated.
 * 1. The volumetric contribution is always calculated by transforming the physical cell to
 * reference space. Thus, every cell's calculation is completely isolated from other cells.
 * 1. The subcell normal vectors obtained in (1-B.53)
 *    - May not be unit vectors
 *    - Do not match (in direction) with the physical normals at faces (local indices) 0, 2 and 4.
 *
 * The metric terms are calculated using the class MetricTerms and stored as a map. Read the class
 * documentation and also that of MetricTerms::reinit() to get an idea of what is being done.
 * Specifically, note the comments made in these documentations about subcell normals. They don't
 * match with cell normals on all faces.
 *
 * Unlike for surface flux calculation, volumetric contribution cannot be unified for conservative
 * and auxiliary variable calculation. This is because auxiliary variables require gradients of
 * conservative variables in all 3 directions. This means a total of 5*3=15 equations need to be
 * set up for these gradients whereas inviscid/diffusive contribution for conservative variable
 * residual requires only 5 equations. Hence a straight-forward unified approach is not possible.
 * However, the inviscid and diffusive high order residual contributions can be calculated in a
 * single loop. The inviscid subcell contribution can then be added separately.
 *
 * @section final_residual_calc Final residual calculation
 *
 * For residual contribution calculation, there are some functions which calculate residual in a
 * given cell and then some outer functions which invoke these functions cell-by-cell. The term
 * "residual" will generally be used for a single cell and the term "rhs" will be used for a vector
 * holding the residual of all cells/dofs. Unlike in @ref face_assem, different "stages" cannot be
 * combined into a single loop even though they use the same formula (eq. (B.14) of Hennemann et al
 * (2021)). Instead, the calculation of conservative variable gradients is kept separate and the
 * calculation of high order inviscid and diffusive contributions is combined in a single function.
 * Low order inviscid contribution is done through a separate function. This separation of inviscid
 * and diffusive contributions allows for any future changes in the algorithm being used for
 * incorporating diffusive terms.
 *
 * This is generally done in the following steps
 * - PLENS::assert_positivity()
 * - PLENS::calc_aux_vars()
 *   - Calculates the auxiliary variables (using eq. (B.14) of Hennemann et al (2021))
 *   - Invokes PLENS::calc_surf_flux() and PLENS::calc_cell_cons_grad() cell-by-cell
 * - PLENS::calc_blender()
 *   - Calculates the value of @f$\alpha@f$. This will subsequently be used for calculating
 *     inviscid contribution
 *   - This function may even be called before PLENS::calc_aux_vars()
 * - PLENS::calc_rhs()
 *   - Calculates the complete residual.
 *   - Internally invokes PLENS::calc_surf_flux(), PLENS::calc_cell_ho_residual() and
 *     PLENS::calc_cell_lo_inv_residual()
 *
 * Then, all these steps are put in a time loop to complete the simulation.
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
     * Vector holding dof-wise viscosity for current RK solution.
     */
    LA::MPI::Vector gcrk_mu;

    /**
     * Vector holding dof-wise thermal conductivity for current RK solution.
     */
    LA::MPI::Vector gcrk_k;

    /**
     * Variable used for calculating blender (@f$\apha@f$) value. This need not be ghosted. Its
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



    void form_neighbor_face_matchings(const double tol = 1e-4);
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
