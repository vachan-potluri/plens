/**
 * @file NavierStokes.h
 * @brief Class for navier stokes solver
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/physics/transformations.h>

#include <cmath> // for pow
#include <array>

#include "var_enums.h"
#include "state.h"

#ifdef DEBUG
#include <iostream>
#include <utilities/printing.h>
#endif

#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H

/**
 * @class NavierStokes
 * @brief Provides functionality for 3D perfect gas laminar simulations of NS equations
 * 
 * The class uses calorically perfect gas model along with Sutherland's law. The class can be
 * constructed by specifying molecular weight, @f$\gamma@f$, Pr, and Sutherland's constants. The
 * relevant equations are (see APS-1 report)
 * @f[
 * p = \rho RT, \quad R=R^0/M\\
 * e = \frac{R}{\gamma-1}T\\
 * a = \sqrt{\gamma RT}
 * \mu = \mu_0 \left( \frac{T}{T_0} \right)^{3/2} \frac{T_0+S}{T+S}
 * k = \frac{\mu c_p}{\text{Pr}}
 * @f]
 * For air and N2, a second constructor is provided to set all these values.
 *
 * Notes about the DGSEM limiter algorithm are linked in WJ-22-Feb-2021.
 * This class provides functionality for flux calculation in 3 stages of the DGSEM algorithm:
 * 1. Auxiliary variable calculation (@f$\nabla\vec{v}@f$ and @f$\nabla T@f$)
 * 2. Inviscid flux calculation
 * 3. Viscous flux calculation
 *
 * Functions required in each of these stages are identified by 'aux', 'inv' and 'vis'. Further,
 * these fluxes are required in two types: 1. surface normal flux and 2. two-point volume flux. The
 * two-point volume flux is unique to the DGSEM algorithm. This distinction in functions is made
 * using 'surf' and 'vol'. Thus in total, there are 6 kinds of fluxes. For all these types of
 * fluxes, the arguments are two appropriate conservative states. The auxiliary variables could
 * also directly be calculated from primitive state, but using conservative state was suggested in
 * Bassi & Rebay (1997). For viscous fluxes, additionally, auxiliary variables
 * (@f$\mathbf{\tau}@f$ and @f$\vec{q}^{''}@f$) are also required.
 *
 * If any identifier from {surf, vol} is missing, then that function is a theoretical function.
 * E.g. NavierStokes::get_inv_flux() takes a conservative state and direction to give the
 * theoretical inviscid flux in that direction.
 *
 * Prefixes 'get' and 'set' are used only for public functions. Private functions don't need them.
 *
 * Strategy for calculating surface normal (numerical) inviscid flux:
 * The relevant function is NavierStokes::get_inv_surf_flux(). This function returns surface normal
 * flux on a face between two states. It internally uses NavierStokes::inv_surf_xflux by rotating
 * coordinate system to align with normal provided. See the function description for more details.
 */
class NavierStokes
{
    public:
    static constexpr int dim = 3; // dimension
    // Choices for inviscid (surface) flux scheme
    enum class inv_surf_flux_scheme{
        hllc,
        rusanov
    };
    
    private:
    double gma_, M_, Pr_, mu0_, T0_, S_;
    std::function< void (const state&, const state&, state&) > inv_surf_xflux;
    void hllc_xflux(const state &lcs, const state &rcs, state &f) const;
    void rusanov_xflux(const state &lcs, const state &rcs, state &f) const;
    
    public:
    NavierStokes(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S,
        const inv_surf_flux_scheme isfs
    );
    NavierStokes(const std::string gas_name, const inv_surf_flux_scheme ifs);
    void set_modelling_params(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S
    );
    void set_inv_surf_flux_scheme(const inv_surf_flux_scheme ifs);
    
    static void assert_positivity(const state &cons);
    static double get_e(const state &cons);
    double get_p(const state &cons) const;
    double get_a(const state &cons) const;
    
    void get_inv_flux(const state &cons, const dealii::Tensor<1,dim> &dir, state &f) const;
    
    void get_inv_surf_flux(
        const state &ocs, const state &ncs, const dealii::Tensor<1,dim> &normal, state &f
    ) const;
    
    #ifdef DEBUG
    static void test();
    void print_modelling_params() const;
    #endif
};

#endif

