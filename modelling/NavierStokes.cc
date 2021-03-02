/**
 * @file NavierStokes.cc
 * @brief Class for navier stokes solver
 */

#include "NavierStokes.h"



/**
 * @brief Constructor. Set values explicitly.
 * 
 * Calls NavierStokes::set_modelling_params() internally
 */
NavierStokes::NavierStokes(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S
)
{
    set_modelling_params(gma, M, Pr, mu0, T0, S);
}



/**
 * @brief Constructor for air and N2
 */
NavierStokes::NavierStokes(const std::string gas_name)
{
    bool gas_supported = (gas_name=="air" || gas_name=="N2");
    AssertThrow(
        gas_supported,
        dealii::StandardExceptions::ExcMessage(
            "Unsupported gas name. Only 'air' and 'N2' are supported"
        )
    );
    
    double gma=1.4, M, Pr=0.69, mu0, T0=273, S;
    if(gas_name == "air"){
        M = 0.029;
        mu0 = 1.716e-5;
        S = 111;
    }
    else{
        // guaranteed to be nitrogen
        M = 0.028;
        mu0 = 1.663e-5;
        S = 107;
    }
    
    // set the values
    set_modelling_params(gma, M, Pr, mu0, T0, S);
}



/**
 * @brief Sets modelling parameters
 * 
 * All values must be in SI units. No assertions are made on the provided data. They are blindly
 * trusted.
 */
void NavierStokes::set_modelling_params(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S
)
{
    gma_ = gma; M_ = M; Pr_ = Pr; mu0_ = mu0; T0_ = T0; S_ = S;
}



/**
 * @brief Asserts positivity of density and thermal energy of given state
 */
void NavierStokes::assert_positivity(const state &cons)
{
    AssertThrow(
        cons[0] > 0,
        dealii::StandardExceptions::ExcMessage("Negative density encountered!")
    );
    
    AssertThrow(
        get_e(cons) > 0, // safe since density already asserted positive
        dealii::StandardExceptions::ExcMessage("Negative energy encountered!")
    );
}



/**
 * @brief Get thermal energy from given conservative state
 *
 * This is purely a mathematical operation, the returned value might be negative. Also, division
 * by density is involved here, so related errors migh occur.
 */
double NavierStokes::get_e(const state &cons)
{
    double ske=0; // specific kinetic energy
    for(int dir=0; dir<dim; dir++){
        ske += pow(cons[1+dir], 2);
    }
    ske *= 0.5/cons[0];
    
    return (cons[4]-ske)/cons[0];
}



/**
 * @brief Get pressure from given conservative state
 * 
 * This function blindly calculates pressure, based on a formula. The return value can be negative
 * too. There is a division by density involved here, so related errors might occur.
 * 
 * Making this function @p static is not possible because this function requires the value of
 * @f$\gamma@f$.
 */
double NavierStokes::get_p(const state &cons) const
{
    double ske=0; // specific kinetic energy
    for(int dir=0; dir<dim; dir++){
        ske += pow(cons[1+dir], 2);
    }
    ske *= 0.5/cons[0];
    
    return (gma_-1)*(cons[4]-ske);
}



/**
 * @brief Calculates inviscid flux in a given direction based on given conservative state
 * 
 * @pre Parameter @p dir has to be a unit vector. No checks regarding this are done
 * 
 * The positivity of @p cons is not checked in this function. This is just a mathematical operation
 */
void NavierStokes::get_inv_flux(
    const state &cons, const dealii::Tensor<1,dim> &dir, state &f
) const
{
    dealii::Tensor<1,dim> vel; // velocity vector
    for(int d=0; d<dim; d++){
        vel[d] = cons[1+d]/cons[0];
    }
    double p = get_p(cons);
    
    double vel_n = dealii::scalar_product(vel, dir); // normal velocity
    
    f[0] = cons[0]*vel_n;
    for(int d=0; d<dim; d++){
        f[1+d] = vel_n*cons[1+d] + p*dir[d];
    }
    f[4] = (cons[4] + p)*vel_n;
}



/* ------------------------------------------------------------------------------------ */



#ifdef DEBUG
void NavierStokes::test()
{
    std::cout << "\n\n\n\nTesting class NavierStokes\n";
    
    {
        NavierStokes ns("air");
        ns.print_modelling_params();
    }
    
    {
        NavierStokes ns("N2");
        ns.print_modelling_params();
    }
    
    {
        NavierStokes ns(1,2,3,4,5,6);
        ns.print_modelling_params();
    }
    
    {
        NavierStokes ns("air");
        state cons = {2,2,4,6,14.0000000001};
        std::cout << "Pressure " << ns.get_p(cons) << "\n";
        std::cout << "Energy " << ns.get_e(cons) << "\n";
        ns.assert_positivity(cons);
        
        std::array<state, 3> fluxes;
        std::array<dealii::Tensor<1,dim>, 3> dir_vecs; // initialised to 0
        for(int d=0; d<dim; d++) dir_vecs[d][d] = 1;
        for(int d=0; d<dim; d++){
            ns.get_inv_flux(cons, dir_vecs[d], fluxes[d]);
            std::cout << d << " direction flux";
            utilities::print_state(fluxes[d]);
        }
    }
}



void NavierStokes::print_modelling_params() const
{
    std::cout << "Modelling parameters (gamma, M, Pr, mu0, T0, S):\n";
    std::cout << gma_ << ", " << M_ << ", " << Pr_ << ", " << mu0_ << ", " << T0_ << ", " << S_
        << "\n";
}
#endif

