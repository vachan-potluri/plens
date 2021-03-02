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



/**
 * @brief HLLC x-flux function.
 *
 * Prefix 'l' and 'r' for left and right. 'c' for conservative. Assumes the interface between
 * @p lcs and @p rcs has normal in x-direction. See detailed class documentation for refs.
 *
 * Positivity of state is not asserted, this function is a pure mathematical operation.
 */
void NavierStokes::hllc_xflux(const state &lcs, const state &rcs, state &f) const
{
    dealii::Tensor<1,dim> xdir({1,0,0});
    // wave speed estimates
    double sql = sqrt(lcs[0]), sqr = sqrt(rcs[0]); // 'sq'uare roots of densities
    double pl = get_p(lcs), pr = get_p(rcs); // pressures
    double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0];
    
    double ut = (ul*sql + ur*sqr)/(sql + sqr); // u tilde
    double Ht = ( (lcs[4]+pl)/sql + (rcs[4]+pr)/sqr )/(sql + sqr); // H tilde
    double at = sqrt((gma_-1)*(Ht - 0.5*ut*ut)); // a tilde
    
    double sl = ut-at, sr = ut+at; // left and right wave speeds
    // star wave speed
    double s = ( pr-pl + lcs[1]*(sl-ul) - rcs[1]*(sr-ur) )/(lcs[0]*(sl-ul) - rcs[0]*(sr-ur));
    
    // cases
    if(sl>0){
        // left state at interface
        get_inv_flux(lcs, xdir, f);
    }
    else if(s>0){
        // left star state at interface
        double temp = (sl-ul)/(sl-s);
        state lss = {
            temp*lcs[0],
            temp*lcs[0]*s,
            temp*lcs[2],
            temp*lcs[3],
            temp*( lcs[4] + (s-ul)*( lcs[0]*s + pl/(sl-ul) ))
        }; // left star state
        
        state lf; // flux based on lcs
        get_inv_flux(lcs, xdir, lf);
        
        for(cvar var: cvar_list) f[var] = lf[var] + sl*(lss[var] - lcs[var]);
    }
    else if(sr>0){
        // right star state at interface
        double temp = (sr-ur)/(sr-s);
        state rss = {
            temp*rcs[0],
            temp*rcs[0]*s,
            temp*rcs[2],
            temp*rcs[3],
            temp*( rcs[4] + (s-ur)*( rcs[0]*s + pr/(sr-ur) ))
        }; // right star state
        
        state rf; // flux based on rcs
        get_inv_flux(rcs, xdir, rf);
        
        for(cvar var: cvar_list) f[var] = rf[var] + sr*(rss[var] - rcs[var]);
    }
    else{
        // right state at interface
        get_inv_flux(rcs, xdir, f);
    }
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
        state cons = {2,2,4,6,15};
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
    
    {
        NavierStokes ns("air");
        state lcs = {1.5,-3,1.5,4.5,23}, rcs = {2,-2,4,4,34}, f; // from WJ-02-Mar-2021
        ns.hllc_xflux(lcs, rcs, f);
        std::cout << "HLLC flux";
        utilities::print_state(f);
    }
}



void NavierStokes::print_modelling_params() const
{
    std::cout << "Modelling parameters (gamma, M, Pr, mu0, T0, S):\n";
    std::cout << gma_ << ", " << M_ << ", " << Pr_ << ", " << mu0_ << ", " << T0_ << ", " << S_
        << "\n";
}
#endif

