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
 * @brief Asserts positivity of density and pressure/thermal energy of given state
 * 
 * Making this function @p static is not possible because this function internally calls
 * NavierStokes::get_p().
 */
void NavierStokes::assert_positivity(const state &cons) const
{
    AssertThrow(
        cons[0] > 0,
        dealii::StandardExceptions::ExcMessage("Negative density encountered!")
    );
    
    AssertThrow(
        get_p(cons) > 0,
        dealii::StandardExceptions::ExcMessage("Negative pressure encountered!")
    );
}



/**
 * @brief Get pressure from given conservative state
 * 
 * This function blindly calculates pressure, based on a formula. The return value can be negative
 * too.
 * 
 * Making this function @p static is not possible because this function requires the value of
 * @f$\gamma@f$. An alternative (to make this static) could be to calculate thermal energy instead
 * of pressure. However, it really adds no value and the non-static nature of this function doesn't
 * add any significant disadvantage.
 */
double NavierStokes::get_p(const state &cons) const
{
    double ske=0; // specific kinetic energy
    for(int dir=0; dir<dim; dir++){
        ske += pow(cons[1+dir], 2);
    }
    ske *= 0.5*cons[0];
    
    return (gma_-1)*(cons[4]-ske);
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
        NavierStokes::state cons = {1,1,1,1,1.50000001};
        std::cout << "Pressure " << ns.get_p(cons) << "\n";
        ns.assert_positivity(cons);
    }
}



void NavierStokes::print_modelling_params() const
{
    std::cout << "Modelling parameters (gamma, M, Pr, mu0, T0, S):\n";
    std::cout << gma_ << ", " << M_ << ", " << Pr_ << ", " << mu0_ << ", " << T0_ << ", " << S_
        << "\n";
}
#endif

