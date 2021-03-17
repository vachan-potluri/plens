/**
 * @file BC.h
 * @brief Base class for BC
 */

#ifndef BC_H
#define BC_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/local_dof_data.h>
#include <dgsem/face_dof_info.h>

#include <array>
#include <map>
#include <functional>

#ifdef DEBUG
#include <iostream>
#include <vector>
#include <utilities/testing.h>
#include <utilities/printing.h>
#include "bc_test_data.h"
#endif

namespace BCs{

using namespace dealii;

/**
 * @class BC
 * @brief Base class for weak Riemann imposition of BCs.
 *
 * This class is a base class for a set of BC objects that impose BCs in weak-Riemann manner. The
 * terminology "weak-Riemann" as well as the implementation of different kinds of BCs is taken from
 * Mengaldo et al (2014) and Bassi & Rebay (1997).
 *
 * Boundary conditions are required at 3 stages in the code.
 * 1. For calculation of auxiliary variables @f$\tau_{ij}@f$ and @f$q_j@f$. In this project, they
 * are calculated using conservative variables, as described by Bassi & Rebay (1997). The BC objects
 * give ghost values of conservative variables @f$\{\rho, \rho u, \rho v, \rho w, \rho E\}@f$
 * which can be used along with NavierStokes::surf_flux_wrappers[0] to get the interface values
 * @f$\{\rho^*, (\rho u)^*, (\rho v)^*, (\rho w)^*, (\rho E)^*\}@f$
 *
 * 2. For calculation of inviscid flux on the boundary surface. Here again, the BC objects give
 * ghost values of conservative variables. NavierStokes::surf_flux_wrappers[1] can then be used to
 * get the inviscid surface flux.
 *
 * 3. For calculation of viscous (diffusive) flux on boundary surface. Here, the BC objects give
 * ghost values of auxiliary variables, along with ghost values of velocity (in conservative state,
 * the remaining variables are not relevant). NavierStokes::surf_flux_wrappers[2] can then be used
 * to obtain the diffusive surface flux.
 *
 * Note that in all cases, the ghost values are given by the BC object and hence the approach is
 * weak-Riemann type. Getter functions are provided for all three stages which are overwridden in
 * classes dedicated to specific boundary condition. And then, similar to the wrappers in
 * NavierStokes, wrappers to the getters are provided here also to enable a uniform call signature
 * which will be useful for assembling surface flux. See NavierStokes class documentation for more
 * details on the usage of wrappers.
 *
 * Keeping in mind periodic BC and spatially varying BC, this base class takes const references to
 * dof handler object, conservative and auxiliary variables. The periodic BC class will additionally
 * have some more entities to be set which will be done by that specific class implementation. In
 * each stage getter, the local information about dof will be taken through LocalDoFData. This
 * will be inevitable for periodic BCs and for spatially varying BCs. Additionally, the local unit
 * normal will also be required. Although the normal could in principle be calculated using dof
 * handler and the local data provided, it would be expensive. Imagine initialising an fe values
 * object for each function call.
 *
 * These BC objects will work for boundaries involving curved manifolds also provided the normal
 * used in the getters is correct.
 *
 * @note For auxiliary variable and viscous fluxes, weak-Riemann and weak-Prescribed approaches are
 * equivalent, it is possible to establish a simple algebraic relation between the weak-Prescribed
 * fluxes and weak-Riemann ghost variables. However, underlying this relation, the algorithms used
 * for auxiliary and viscous flux evaluation are sort of assumed. Currently, all derived classes
 * assume BR1 algorithm (simple avg, see NavierStokes) for these two fluxes.
 */
class BC
{
    public:
    static constexpr int dim = 3; // dimension
    
    const DoFHandler<dim>& dof_handler;
    const FaceDoFInfo fdi;
    
    // Refs to cvars and avars. These are assumed to hold the values of cvars and avars at the time
    // of calling getter functions and thus used in the same to finally give the ghost values. One
    // way to set these correectly is to assign current RK step solutions
    const std::array<LA::MPI::Vector, 5>& g_cvars;
    const std::array<LA::MPI::Vector, 9>& g_avars;
    
    /**
     * @brief Wrappers to get_ghost_stagex for x=1,2,3
     *
     * Having a unified call signature necessitates the usage of CAvars, like in NavierStokes class.
     * These are set in BC::set_wrappers() which is called from BC::BC(). Hence, these are
     * automatically set in all derived classes if the derived ctor calls the base ctor. To know if
     * inheritance of `this` works as expected, see the small code snippet in WJ-17-Mar-2021.
     */
    std::array< std::function<
        void(const LocalDoFData&, const Tensor<1,dim>&, CAvars&)
    >, 3 > get_ghost_wrappers;
    
    protected:
    /**
     * A map between cell id and cell iterator for owned cells.
     *
     * This will work only for unrefined grids. The relation satisfied by a pair @p p in this map
     * is: `p.second->index() = p.first`. For hp grids, a cell index cannot uniquely determine a
     * cell iterator.
     */
    std::map<psize, DoFHandler<dim>::active_cell_iterator> cell_map_;
    
    public:
    BC(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav
    );
    virtual ~BC();
    
    /**
     * @brief Get ghost values of conservative variables for auxiliary variable calculation
     *
     * See the class documentation for more details. This needs to be overridden in derived classes.
     */
    virtual void get_ghost_stage1(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const {}
    
    /**
     * @brief Get ghost values of conservative variables for inviscid flux calculation
     *
     * See the class documentation for more details. This needs to be overridden in derived classes.
     */
    virtual void get_ghost_stage2(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const {}
    
    /**
     * @brief Get ghost values of conservative and auxiliary variables for calculation of viscous
     * flux.
     *
     * See the class documentation for more details. This needs to be overridden in derived classes.
     */
    virtual void get_ghost_stage3(
        const LocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        CAvars &cav_gh
    ) const {}
    
    protected:
    void form_cell_map();
    psize get_global_dof_id(const LocalDoFData &ldd) const;
    void get_state(const LocalDoFData &ldd, State &s) const;
    void get_avars(const LocalDoFData &ldd, Avars &a) const;
    void get_cavars(const LocalDoFData &ldd, CAvars &ca) const;
    void set_wrappers();
    
    public:
    #ifdef DEBUG
    static void test();
    #endif
};

}

#endif

