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
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S,
    const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
    const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar
)
{
    set_modelling_params(gma, M, Pr, mu0, T0, S);
    set_inv_surf_flux_scheme(isfs);
    set_inv_vol_flux_scheme(ivfs);
}



/**
 * @brief Constructor for air and N2
 */
NavierStokes::NavierStokes(
    const std::string gas_name,
    const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
    const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar
)
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
    set_inv_surf_flux_scheme(isfs);
    set_inv_vol_flux_scheme(ivfs);
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
 * @brief Sets inviscid surface (numerical) flux function: NavierStokes::inv_surf_xflux
 */
void NavierStokes::set_inv_surf_flux_scheme(const inv_surf_flux_scheme isfs)
{
    if(isfs == inv_surf_flux_scheme::hllc){
        inv_surf_xflux = [=](const state &lcs, const state &rcs, state &f){
            this->hllc_xflux(lcs, rcs, f);
        };
    }
    else{
        inv_surf_xflux = [=](const state &lcs, const state &rcs, state &f){
            this->rusanov_xflux(lcs, rcs, f);
        };
    }
}



/**
 * @brief Sets the invsicid volume flux function NavierStokes::inv_vol_flux
 */
void NavierStokes::set_inv_vol_flux_scheme(const inv_vol_flux_scheme ivfs)
{
    // only one option available currently
    get_inv_vol_flux = [=](
        const state &cs1, const state &cs2, const dealii::Tensor<1,dim> &dir, state &f){
        this->chandrashekhar_flux(cs1, cs2, dir, f);
    };
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
 * @brief Get speed of sound from given conservative state
 *
 * Internally calls NavierStokes::get_p() and returns @f$\sqrt{\gamma p/\rho}@f$. Positivity of
 * density and pressure are not checked.
 */
double NavierStokes::get_a(const state &cons) const
{
    return sqrt(gma_*get_p(cons)/cons[0]);
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
    
    double vel_n = dealii::scalar_product(vel, dir); // velocity in the direction 'dir'
    
    f[0] = cons[0]*vel_n;
    for(int d=0; d<dim; d++){
        f[1+d] = vel_n*cons[1+d] + p*dir[d];
    }
    f[4] = (cons[4] + p)*vel_n;
}



/**
 * @brief Calculates inviscid surface (numerical) flux
 *
 * This function uses NavierStokes::inv_surf_xflux (one of NavierStokes::hllc_xflux() and
 * NavierStokes::rusanov_xflux() based on the argument provided in constructor). The algorithm
 * for this function is described in detail in WJ-23-Feb-2021.
 *
 * 1. Rotate the coordinate system such that x-direction aligned with @p normal. Store the rotation
 * matrix.
 * 2. Rotate the velocites to orient them along new coordinate axes.
 * 3. Pass the rotated states to NavierStokes::inv_surf_xflux and get the normal flux
 * 4. Rotate back the flux momentum components of the above calculated flux
 *
 * @param[in] ocs 'Owner' conservative state
 * @param[in] ncs 'Neighbor' conservative state
 * @param[in] normal The normal vector of the face shared by owner and neighbor. This must be a
 * unit vector and must point from owner side to neighbor side
 * @param[out] f The resultant surface normal flux
 *
 * Positivity of @p ocs and @p ncs is not asserted
 */
void NavierStokes::get_inv_surf_flux(
    const state &ocs, const state &ncs, const dealii::Tensor<1,dim> &normal, state &f
) const
{
    // Step 1: rotate coordinate system
    std::array<dealii::Tensor<1,3>, 3> dir_vecs; // initialised to 0
    for(int d=0; d<dim; d++) dir_vecs[d][d] = 1;
    
    dealii::Tensor<1,dim> m = dealii::cross_product_3d(dir_vecs[0], normal); // m = x cross n
    double M = m.norm(); // magnitude of m
    m /= M; // now m is a unit vector <-- rotation axis
    double theta = asin(M); // <-- rotation angle
    
    // rotation matrix
    dealii::FullMatrix<double> R(dim);
    R.copy_from(
        dealii::Physics::Transformations::Rotations::rotation_matrix_3d(
            dealii::Point<dim>(m), theta
        )// returns tensor
    ); // copies from second order tensor
    
    // Step 2: get rotated states
    dealii::Vector<double> osmom(dim), nsmom(dim), // owner and neighbor specific momentum
        osmom_r(dim), nsmom_r(dim); // rotated specific momentum
    for(int d=0; d<dim; d++){
        osmom[d] = ocs[1+d];
        nsmom[d] = ncs[1+d];
    }
    // get the momentum components wrt rotated coordinate system
    R.Tvmult(osmom_r, osmom); // osmom_r = R^-1 * osmom, R^T = R^{-1}
    R.Tvmult(nsmom_r, nsmom);
    
    state ocs_r, ncs_r; // rotated states
    ocs_r[0] = ocs[0];
    ncs_r[0] = ncs[0];
    for(int d=0; d<dim; d++){
        ocs_r[1+d] = osmom_r[d];
        ncs_r[1+d] = nsmom_r[d];
    }
    ocs_r[4] = ocs[4];
    ncs_r[4] = ncs[4];
    
    // Step 3: get normal flux wrt rotated coordinate system
    state f_r;
    inv_surf_xflux(ocs_r, ncs_r, f_r);
    
    // Step 4: rotate back coordinate system and obtain the appropriate momentum components
    dealii::Vector<double> mom_flux_r(3), mom_flux(3); // momentum fluxes
    for(int d=0; d<dim; d++) mom_flux_r[d] = f_r[d];
    // get momentum flux components w.r.t original coordinate system
    R.vmult(mom_flux, mom_flux_r);
    
    f[0] = f_r[0];
    for(int d=0; d<dim; d++) f[d] = mom_flux[d];
    f[4] = f_r[4];
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



/**
 * @brief Rusanov x-flux function
 *
 * Every other detail same as in NavierStokes::hllc_xflux()
 */
void NavierStokes::rusanov_xflux(const state &lcs, const state &rcs, state &f) const
{
    dealii::Tensor<1,dim> xdir({1,0,0});
    
    state lf, rf; // left and right conservative fluxes
    double al = get_a(lcs), ar = get_a(rcs); // left and right sound speeds
    double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0]; // left & right flow speeds
    double S; // the "single wave" speed
    
    get_inv_flux(lcs, xdir, lf);
    get_inv_flux(rcs, xdir, rf);
    
    S = std::max(fabs(ul)+al, fabs(ur)+ar);
    
    for(cvar var: cvar_list) f[var] = 0.5*(lf[var] + rf[var]) - 0.5*S*(rcs[var] - lcs[var]);
}



/**
 * @brief Chandrashekhar inviscid volume flux.
 *
 * See eqs (3.16, 3.18-3.20) of Gassner, Winters & Kopriva (2016). Positivity of the states
 * provided is not asserted. @p dir has to be a unit vector.
 */
void NavierStokes::chandrashekhar_flux(
    const state &cs1, const state &cs2, const dealii::Tensor<1,dim> &dir, state &f
) const
{
    double p1 = get_p(cs1), p2 = get_p(cs2);
    double beta1 = 0.5*cs1[0]/p1, beta2 = 0.5*cs2[0]/p2;
    double beta_ln = (beta1-beta2)/(log(beta1) - log(beta2));
    
    double p_hat = 0.5*(cs1[0]+cs2[0])/(beta1+beta2);
    double rho_ln = (cs1[0]-cs2[0])/(log(cs1[0]) - log(cs2[0]));
    
    dealii::Tensor<1,dim> vel_avg, vel_sq_avg; // sq for 'sq'uare
    double v1, v2; // temporary quantities
    for(int d=0; d<dim; d++){
        v1 = cs1[1+d]/cs1[0];
        v2 = cs2[1+d]/cs2[0];
        vel_avg[d] = 0.5*(v1+v2);
        vel_sq_avg[d] = 0.5*(v1*v1 + v2*v2);
    }
    
    double H_hat = 0.5/((gma_-1)*beta_ln) + p_hat/rho_ln; // initialise
    for(int d=0; d<dim; d++){
        H_hat += vel_avg[d]*vel_avg[d] - 0.5*vel_sq_avg[d];
    }
    
    double vel_n = dealii::scalar_product(vel_avg, dir); // velocity in the direction 'dir'
    
    f[0] = vel_n*rho_ln;
    for(int d=0; d<dim; d++){
        f[1+d] = rho_ln*vel_n*vel_avg[d] + p_hat*dir[d];
    }
    f[4] = rho_ln*vel_n*H_hat;
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
        state lcs = {1.5,3,1.5,4.5,23}, rcs = {2,2,4,4,34}, f; // from WJ-02-Mar-2021
        ns.hllc_xflux(lcs, rcs, f);
        std::cout << "HLLC x flux";
        utilities::print_state(f);
        
        ns.rusanov_xflux(lcs,rcs,f);
        std::cout << "Rusanov x flux";
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

