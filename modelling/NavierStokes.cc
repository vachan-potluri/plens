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
        NavierStokes ns("nitrogen");
        ns.print_modelling_params();
    }
}



void NavierStokes::print_modelling_params() const
{
    std::cout << "Modelling parameters (gamma, M, Pr, mu0, T0, S):\n";
    std::cout << gma_ << ", " << M_ << ", " << Pr_ << ", " << mu0_ << ", " << T0_ << ", " << S_
        << "\n";
}
#endif

