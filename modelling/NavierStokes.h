/**
 * @file NavierStokes.h
 * @brief Class for navier stokes solver
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>

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
 */
class NavierStokes
{
    public:
    static constexpr int dim = 3; // dimension
    
    private:
    double gma_, M_, Pr_, mu0_, T0_, S_;
    
    public:
    NavierStokes(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S
    );
    NavierStokes(const std::string gas_name);
    void set_modelling_params(
        const double gma, const double M, const double Pr,
        const double mu0, const double T0, const double S
    );
    
    void assert_positivity(const state &cons) const;
    double get_p(const state &cons) const;
    
    void get_inv_flux(const state &cons, const dealii::Tensor<1,dim> &dir, state &f) const;
    
    #ifdef DEBUG
    static void test();
    void print_modelling_params() const;
    #endif
};

#endif

