/**
 * @file navier_stokes.h
 * @brief Class for navier stokes solver
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/physics/transformations.h>

#include <cmath> // for pow
#include <array>

#include "var_enums.h"
#include "state.h"
#include "avars.h"
#include "cavars.h"

#ifdef DEBUG
#include <iostream>
#include <utilities/printing.h>
#include <utilities/testing.h>
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
 * For air and N2, a second constructor is provided to set all these values without explicitly
 * passing them in constructor.
 *
 * Notes about the DGSEM limiter algorithm are linked in WJ-22-Feb-2021 entry.
 * This class provides functionality for flux calculation in 3 stages of the DGSEM algorithm:
 * 1. Auxiliary variable calculation (@f$\nabla\vec{v}@f$ and @f$\nabla T@f$). To be mathematically
 * precise, @f$\tau@f$ and @f$\vec{q}^{''}@f$ are the auxiliary variables. Since these can be
 * obtained through linear operations of @f$\nabla\vec{v}@f$ and @f$\nabla T@f$ (appropriately 
 * using @f$\mu@f$ and @f$k@f$), the term is used to refer to either set of variables depending on
 * the context.
 * 2. Inviscid flux calculation
 * 3. Viscous flux calculation
 *
 * Functions required in each of these stages are identified by 'aux', 'inv' and 'vis'. Further,
 * these fluxes are required in two types: 1. surface normal flux and 2. two-point volume flux. The
 * two-point volume flux is unique to the DGSEM algorithm. This distinction in functions is made
 * using phrases 'surf' and 'vol'. Thus in total, there are 3*2 = 6 kinds of fluxes. For the first
 * two kinds of fluxes (aux and inv), the arguments are two appropriate conservative states. The
 * auxiliary variables are calculated indirectly using the conservative variables, as was suggested
 * in Bassi & Rebay (1997). For example,
 * @f[
 * \nabla u = \frac{1}{\rho} \left( \nabla (\rho u) - u \nabla\rho \right)
 * @f]
 * This would make the calculation of @f$\nabla T@f$ very cumbersome, but some helping functions
 * will probably be provided here when the code reaches that stage. Different algorithms differ in
 * how the auxiliary variable flux is calculated at the surface (and also in volume specifically for
 * the DGSEM method). Currently, only a simple average, as used by Bassi & Rebay (1997) is
 * implemented. This is labelled BR1. As a result, the 'aux' functions return
 * @f$\{\rho^*, (\rho u)^*, (\rho v)^*, (\rho w)^*, (\rho E)^*\}@f$ for surface flux (note usage of
 * conservative variables).
 *
 * One important note here. Unlike inv fluxes, the aux fluxes need not vary with direction. For
 * illustration, consider solving
 * @f[
 * W - \frac{\partial \rho}{\partial x} = 0
 * @f]
 * where @f$W=\partial \rho/\partial x@f$ is the auxiliary/additional variable we wish to solve. The
 * surface term in weak formulation would be
 * @f[
 * \int_{\partial \Omega_e} \rho^* \vec{n}\,dA
 * @f]
 * and the volume contribution will be calculated using two-point fluxes as described in Appendix B
 * of Gassner, Winters & Kopriva (2016). Note that for this specific case, @f$G@f$ and @f$H@f$ used
 * are in eqs (B.5) and (B.6) identically zero. And moreover, if the gradient direction is changed
 * such that we now require @f$W=\partial \rho/\partial y@f$, then @f$F@f$ and @f$H@f$ are
 * identically zero along with a special equality: @f$F(W)@f$ of the first case is exactly equal to
 * @f$G(W)@f$ of the second case. Thus, the direction is really not required. Only the metric term
 * factors change for different directions and this is taken care elsewhere. That said, there can
 * however be formulations where both surface and volume fluxes are direction dependent. Hence the
 * variables NavierStokes::get_aux_surf_flux and NavierStokes::get_aux_vol_flux have direction.
 * However, BR1 does not require direction.
 *
 * For calculating surface normal (numerical) inviscid flux, the relevant function is
 * NavierStokes::get_inv_surf_flux(). This function returns surface normal flux on a face between
 * two states. It internally uses NavierStokes::inv_surf_xflux by rotating coordinate system to
 * align with normal provided. See the function description for more details.
 * Ideally, it would be preferred to link inviscid surface flux to the inviscid volume flux as
 * described in section 3.4 of Gassner, Winters & Kopriva (2016). However, this is not strictly
 * necessary. Currently the linking is not being done.
 * For volume flux, Chandrashekhar flux described in eqs. (3.16, 3.18-20) of Gassner, Winters &
 * Kopriva (2016) is used. However, the framework itself is kept open to addition of any further
 * options. The volume fluxes can handle direction directly and there is no need for the two step
 * approach taken for inviscid surface fluxes. Also, unlike inv surf flux which is a member
 * function, the inv vol flux is a variable of type <tt>std::function</tt>.
 * Moreover, volume fluxes are already formulated in a directional manner and don't need the
 * two-function approach used here for inviscid surface flux.
 * For both surface and volume fluxes, the options are private functions and the getters (which will
 * be used directly eventually) are public functions/function variables.
 *
 * For viscous fluxes, additionally, auxiliary variables
 * (@f$\mathbf{\tau}@f$ and @f$\vec{q}^{''}@f$) are also required. For this purpose, a new class
 * cavars is introduced. It is just a container like class to allow storing conservative and
 * auxiliary variables in a single address. Currently, only BR1 type viscous/diffusive flux is
 * implemented. It just returns an average of diffusive fluxes based on variables of both sides.
 *
 * Finally, this class provides wrappers: NavierStokes::surf_flux_wrappers which is an array of surf
 * flux functions in the order aux, inv, dif. This is useful for assembling surface flux terms of
 * all cells. It can avoid code repetition. The repetition is otherwise inevitable because diffusive
 * fluxes can be calculated only after the auxiliary variables are calculated using auxiliary
 * fluxes. This means, at least two decoupled loops will be required: either aux, inv+dif or
 * aux+inv, dif. The wrappers provide unified interface to all types of fluxes thus enabling
 * automation of assembly.
 *
 * @note If any identifier from {surf, vol} is missing, then that function is a theoretical
 * function. E.g. NavierStokes::get_inv_flux() takes a conservative state and direction to give the
 * theoretical inviscid flux in that direction.
 *
 * @note Prefixes 'get' and 'set' are used only for public functions. Private functions don't have
 * them.
 *
 * 
 */
class NavierStokes
{
    public:
    static constexpr int dim = 3; // dimension
    
    /**
     * Universal gas constant in SI units. Accurate upto 10 significant digits. Source: Wikipedia
     */
    static constexpr double R0 = 8.314462618;
    
    // Choices for auxiliary variable surface and volume flux scheme
    enum class aux_surf_flux_scheme{
        BR1
    };
    enum class aux_vol_flux_scheme{
        BR1
    };
    
    // Choices for inviscid surface and volume flux scheme
    enum class inv_surf_flux_scheme{
        hllc,
        rusanov
    };
    enum class inv_vol_flux_scheme{
        chandrashekhar
    };
    
    // Choices for diffusive surface and volume flux scheme
    enum class dif_surf_flux_scheme{
        BR1
    };
    enum class dif_vol_flux_scheme{
        BR1
    };
    
    private:
    double gma_, M_, Pr_, mu0_, T0_, S_;
    
    std::function< void (const State&, const State&, State&) > inv_surf_xflux;
    
    // inv surf fluxes
    void hllc_xflux(const State &lcs, const State &rcs, State &f) const;
    void rusanov_xflux(const State &lcs, const State &rcs, State &f) const;
    
    void chandrashekhar_flux(
        const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
    ) const; // inv vol flux
    
    static void br1_flux(const State &cs1, const State &cs2, State &f); // aux surf & vol flux
    
    static void br1_flux(
        const CAvars &cav1, const CAvars &cav2, const dealii::Tensor<1,dim> &dir, State &f
    ); // dif surf & vol flux
    
    public:
    std::function< void (
        const State&, const State&, const dealii::Tensor<1,dim>&, State&
    ) >
        get_aux_surf_flux,
        get_aux_vol_flux,
        get_inv_vol_flux;
    
    std::function< void (
        const CAvars&, const CAvars&, const dealii::Tensor<1,dim>&, State&
    ) >
        get_dif_surf_flux,
        get_dif_vol_flux;
    
    std::array<
        std::function< void (const CAvars&, const CAvars&, const dealii::Tensor<1,dim>&, State&) >,
        3
    > surf_flux_wrappers, vol_flux_wrappers;

    std::array<
        std::function<void(const CAvars&, const dealii::Tensor<1,dim>&, State&)>,
        3
    > flux_wrappers;
    
    NavierStokes(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S,
        const aux_surf_flux_scheme asfs = aux_surf_flux_scheme::BR1,
        const aux_vol_flux_scheme avfs = aux_vol_flux_scheme::BR1,
        const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
        const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar,
        const dif_surf_flux_scheme dsfs = dif_surf_flux_scheme::BR1,
        const dif_vol_flux_scheme dvfs = dif_vol_flux_scheme::BR1
    );
    NavierStokes(
        const std::string gas_name,
        const bool inviscid = false,
        const aux_surf_flux_scheme asfs = aux_surf_flux_scheme::BR1,
        const aux_vol_flux_scheme avfs = aux_vol_flux_scheme::BR1,
        const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
        const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar,
        const dif_surf_flux_scheme dsfs = dif_surf_flux_scheme::BR1,
        const dif_vol_flux_scheme dvfs = dif_vol_flux_scheme::BR1
    );
    void set_modelling_params(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S
    );
    void set_aux_surf_flux_scheme(const aux_surf_flux_scheme asfs);
    void set_aux_vol_flux_scheme(const aux_vol_flux_scheme avfs);
    void set_inv_surf_flux_scheme(const inv_surf_flux_scheme isfs);
    void set_inv_vol_flux_scheme(const inv_vol_flux_scheme ivfs);
    void set_dif_surf_flux_scheme(const dif_surf_flux_scheme dsfs);
    void set_dif_vol_flux_scheme(const dif_vol_flux_scheme dvfs);
    
    void set_wrappers();
    
    static void assert_positivity(const State &cons);
    static double get_e(const State &cons);
    double get_p(const State &cons) const;
    double get_a(const State &cons) const;
    double get_M(const State &cons) const;
    
    void get_inv_flux(const State &cons, const dealii::Tensor<1,dim> &dir, State &f) const;
    void get_inv_surf_flux(
        const State &ocs, const State &ncs, const dealii::Tensor<1,dim> &normal, State &f
    ) const;
    
    static void get_stress_tensor(const Avars &av, dealii::SymmetricTensor<2,dim> &st);
    static void get_dif_flux(const CAvars &cav, const dealii::Tensor<1,dim> &dir, State &f);
    
    /**
     * Gives @f$\gamma@f$ value held by this instance
     */
    inline double get_gma() const {return gma_;}
    
    /**
     * Gives @f$R=R_0/M@f$ value
     */
    inline double get_R() const {return R0/M_;}

    /**
     * Gives @f$c_v=\frac{R}{\gamma-1}@f$
     */
    inline double get_cv() const {return get_R()/(gma_-1);}

    /**
     * Gives @f$\mu@f$ based on temperature provided. Blindly calculates the value based on
     * Sutherland's formula, without any assertion on `T`.
     */
    inline double get_mu(const double T) const
    {
        return mu0_*std::pow(T/T0_, 1.5)*(T0_+S_)/(T+S_);
    }

    /**
     * Gives @f$k@f$ based on viscosity provided. Assertion on positivity of viscosity is not done.
     */
    inline double get_k(const double mu) const {return mu*get_R()*gma_/((gma_-1)*Pr_);}
    
    /**
     * Gives the conservative state using density, velocity and pressure
     */
    inline void prim_to_cons(
        const double rho,
        const dealii::Tensor<1,dim> &vel,
        const double p,
        State &cons
    ) const
    {
        cons[0] = rho;
        cons[4] = p/(get_gma()-1); // initialise
        for(int d=0; d<dim; d++){
            cons[1+d] = rho*vel[d];
            cons[4] += 0.5*rho*vel[d]*vel[d];
        }
    }
    
    #ifdef DEBUG
    static void test();
    void print_modelling_params() const;
    #endif
};

#endif

