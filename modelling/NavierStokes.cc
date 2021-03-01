/**
 * @file NavierStokes.cc
 * @brief Class for navier stokes solver
 */

#include "NavierStokes.h"

/**
 * @brief Constructor. Set values explicitly. All values must be in SI units
 */
NavierStokes::NavierStokes(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S
): gma_(gma), M_(M), Pr_(Pr), mu0_(mu0), T0_(T0), S_(S)
{}

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
    NavierStokes(gma, M, Pr, mu0, T0, S);
}

