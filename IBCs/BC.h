/**
 * @file BC.h
 * @brief Base class for BC
 */

#ifndef BC_H
#define BC_H

namespace BCs{

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
 * dof handler object. The periodic BC class will additionally have some more entities to be set
 * which will be done by that specific class implementation. In each stage getter, the local
 * information about dof will be taken through ldof_data. This will be inevitable for periodic BCs
 * and for spatially varying BCs. Additionally, the local unit normal will also be required.
 */

}

#endif

