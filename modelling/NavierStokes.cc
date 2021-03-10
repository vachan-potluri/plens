/**
 * @file NavierStokes.cc
 * @brief Class for navier stokes solver
 */

#include "NavierStokes.h"



/**
 * @brief Constructor. Set parameter values explicitly.
 *
 * Internally calls NavierStokes::set_modelling_params(), NavierStokes::set_inv_surf_flux_scheme()
 * and NavierStokes::set_inv_vol_flux_scheme().
 */
NavierStokes::NavierStokes(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S,
    const aux_surf_flux_scheme asfs = aux_surf_flux_scheme::BR1,
    const aux_vol_flux_scheme avfs = aux_vol_flux_scheme::BR1,
    const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
    const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar,
    const dif_surf_flux_scheme dsfs = dif_surf_flux_scheme::BR1,
    const dif_vol_flux_scheme dvfs = dif_vol_flux_scheme::BR1
)
{
    set_modelling_params(gma, M, Pr, mu0, T0, S);
    set_aux_surf_flux_scheme(asfs);
    set_aux_vol_flux_scheme(avfs);
    set_inv_surf_flux_scheme(isfs);
    set_inv_vol_flux_scheme(ivfs);
    set_dif_surf_flux_scheme(dsfs);
    set_dif_vol_flux_scheme(dvfs);
    set_wrappers();
}



/**
 * @brief Special constructor for air and N2.
 *
 * Internally calls NavierStokes::set_modelling_params(), NavierStokes::set_inv_surf_flux_scheme()
 * and NavierStokes::set_inv_vol_flux_scheme().
 */
NavierStokes::NavierStokes(
    const std::string gas_name,
    const aux_surf_flux_scheme asfs = aux_surf_flux_scheme::BR1,
    const aux_vol_flux_scheme avfs = aux_vol_flux_scheme::BR1,
    const inv_surf_flux_scheme isfs = inv_surf_flux_scheme::hllc,
    const inv_vol_flux_scheme ivfs = inv_vol_flux_scheme::chandrashekhar,
    const dif_surf_flux_scheme dsfs = dif_surf_flux_scheme::BR1,
    const dif_vol_flux_scheme dvfs = dif_vol_flux_scheme::BR1
)
{
    bool gas_supported = (gas_name=="air" || gas_name=="N2");
    AssertThrow(
        gas_supported,
        dealii::StandardExceptions::ExcMessage(
            "Unsupported gas name. Only 'air' and 'N2' are currently supported"
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
    set_aux_surf_flux_scheme(asfs);
    set_aux_vol_flux_scheme(avfs);
    set_inv_surf_flux_scheme(isfs);
    set_inv_vol_flux_scheme(ivfs);
    set_dif_surf_flux_scheme(dsfs);
    set_dif_vol_flux_scheme(dvfs);
    set_wrappers();
}



/**
 * @brief Sets modelling parameters
 * 
 * All values must be in SI units. No assertions are made on the provided data. They are blindly
 * trusted.
 *
 * @pre All parameters provided must be positive and in SI units
 */
void NavierStokes::set_modelling_params(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S
)
{
    gma_ = gma; M_ = M; Pr_ = Pr; mu0_ = mu0; T0_ = T0; S_ = S;
}



/**
 * @brief Sets the auxiliary surface (numerical) flux function: NavierStokes::get_aux_surf_flux.
 *
 * @note Currently only BR1 flux is supported, @p asfs is unused param
 */
void NavierStokes::set_aux_surf_flux_scheme(const aux_surf_flux_scheme asfs)
{
    get_aux_surf_flux = [=](
        const state &cs1, const state &cs2, const dealii::Tensor<1,dim> &dir, state &f){
        this->br1_flux(cs1, cs2, f); // dir unused for BR1
    };
}



/**
 * @brief Sets the auxiliary volume (numerical) flux function: NavierStokes::get_aux_vol_flux.
 *
 * @note Currently only BR1 flux is supported, @p asfs is unused param
 */
void NavierStokes::set_aux_vol_flux_scheme(const aux_vol_flux_scheme avfs)
{
    get_aux_vol_flux = [=](
        const state &cs1, const state &cs2, const dealii::Tensor<1,dim> &dir, state &f){
        this->br1_flux(cs1, cs2, f); // dir unused for BR1
    };
}



/**
 * @brief Sets the inviscid surface (numerical) uni-directional flux function:
 * NavierStokes::inv_surf_xflux.
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
 *
 * @note Currently only Chandrashekhar volume flux is available. So he parameter @p ivfs is not used
 * at all currently.
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
 * @brief Sets the diffusive surface flux function
 *
 * @note Currently only BR1 flux is supported, @p dsfs is an unused param here
 */
void NavierStokes::set_dif_surf_flux_scheme(const dif_surf_flux_scheme dsfs)
{
    get_dif_surf_flux = [=](
        const cavars &cav1, const cavars cav2, const dealii::Tensor<1,dim> &dir, state &f
    ){
        this->br1_flux(cav1, cav2, dir, f);
    };
}



/**
 * @brief Sets the diffusive volume flux function
 *
 * @note Currently only BR1 flux is supported, @p dvfs is an unused param here
 */
void NavierStokes::set_dif_vol_flux_scheme(const dif_vol_flux_scheme dvfs)
{
    get_dif_vol_flux = [=](
        const cavars &cav1, const cavars cav2, const dealii::Tensor<1,dim> &dir, state &f
    ){
        this->br1_flux(cav1, cav2, dir, f);
    };
}



/**
 * @brief Sets surface flux wrappers
 *
 * See the class documentation for more details
 *
 * @pre All the surf/vol setters must be called before invoking this function
 */
void NavierStokes::set_wrappers()
{
    // surface flux wrappers
    surf_flux_wrappers[0] = [=](
        const cavars &cav1, const cavars cav2, const dealii::Tensor<1,dim> &dir, state &f
    ){
        this->get_aux_surf_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // aux
    
    surf_flux_wrappers[1] = [=](
        const cavars &cav1, const cavars cav2, const dealii::Tensor<1,dim> &dir, state &f
    ){
        this->get_inv_surf_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // inv
    
    surf_flux_wrappers[2] = get_dif_surf_flux; // dif, directly equate function objects
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
 *
 * @pre <tt>cons[0]</tt> must be positive
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
 * @f$\gamma@f$. Doesn't use NavierStokes::get_e(). That approach would involve one additional
 * division and multiplication. Since this function is often called, that would be redundant.
 *
 * @pre <tt>cons[0]</tt> must be positive
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
 *
 * @pre @p cons must pass the test of NavierStokes::assert_positivity()
 */
double NavierStokes::get_a(const state &cons) const
{
    return sqrt(gma_*get_p(cons)/cons[0]);
}



/**
 * @brief Calculates inviscid flux in a given direction based on given conservative state
 * 
 * @pre @p cons must pass the test of NavierStokes::assert_positivity()
 * @pre @p dir has to be a unit vector
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
 * 2. Obtain the velocity components w.r.t. this new coordinate system. The components obtained this
 * way are equivalent to the components of velocity when rotated in a way opposite to what is
 * described in the first point. Thus the new components are @f$R^{-1}\vec{v} = R^T\vec{v}@f$.
 * 3. Pass the 'rotated' states to NavierStokes::inv_surf_xflux and get the 'rotated' normal flux.
 * The adjective rotated here is for the coordinate system, and not the state itself.
 * 4. Rotate back the coordinate system and get the momentum flux components in actual coordinate
 * system using the above calculated flux.
 *
 * @param[in] ocs 'Owner' conservative state
 * @param[in] ncs 'Neighbor' conservative state
 * @param[in] normal The normal vector of the face shared by owner and neighbor. This @b must be a
 * unit vector and @b must point from owner side to neighbor side
 * @param[out] f The resultant surface normal flux
 *
 * @pre @p ocs and @p ncs must pass the test of NavierStokes::assert_positivity()
 * @pre @p normal has to be a unit vector
 */
void NavierStokes::get_inv_surf_flux(
    const state &ocs, const state &ncs, const dealii::Tensor<1,dim> &normal, state &f
) const
{
    // Step 1: rotate coordinate system
    dealii::Tensor<1,dim> xdir({1,0,0});
    
    dealii::Tensor<1,dim> m = dealii::cross_product_3d(xdir, normal); // m = x cross n
    double M = m.norm(); // magnitude of m
    
    dealii::FullMatrix<double> R(dim); // rotation matrix
    if(M > 1e-3){
        // this tolerance corrsponds to an angle of ~0.06 degrees between x and n
        m /= M; // now m is a unit vector <-- rotation axis
        double theta = asin(M); // <-- rotation angle
        
        R.copy_from(
            dealii::Physics::Transformations::Rotations::rotation_matrix_3d(
                dealii::Point<dim>(m), theta
            )// returns tensor
        ); // copies from second order tensor
    }
    else{
        // either x is parallel or anti-parallel to n
        R = dealii::IdentityMatrix(dim);
        if(normal[0] < 0){
            // n is anti-parallel to x
            R *= -1;
        }
    }
    
    // Step 2: get rotated states
    dealii::Vector<double> osmom(dim), nsmom(dim), // owner and neighbor specific momentum
        osmom_r(dim), nsmom_r(dim); // rotated specific momentum
    for(int d=0; d<dim; d++){
        osmom[d] = ocs[1+d];
        nsmom[d] = ncs[1+d];
    }
    // get the momentum components wrt rotated coordinate system
    R.Tvmult(osmom_r, osmom); // osmom_r = R^{-1} * osmom, R^T = R^{-1}
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
    for(int d=0; d<dim; d++) mom_flux_r[d] = f_r[1+d];
    // get momentum flux components w.r.t original coordinate system
    R.vmult(mom_flux, mom_flux_r); // mom_flux = R * mom_flux_r
    
    f[0] = f_r[0];
    for(int d=0; d<dim; d++) f[1+d] = mom_flux[d];
    f[4] = f_r[4];
}



/**
 * @brief Gives the symmetric stress tensor based on avars provided
 */
void NavierStokes::get_stress_tensor(const avars &av, dealii::SymmetricTensor<2,dim> &st)
{
    int i=0;
    for(int row=0; row<dim; row++){
        for(int col=row; col<dim; col++){
            st[row][col] = av[i]; // automatically sets the symmetrical part too
            i++;
        } // loop over cols with col >= row
    } // loop over rows
}



/**
 * @brief Calculates diffusive flux in direction @p dir based on the conservative and auxiliary
 * variables provided by @p cav
 *
 * Internally uses NavierStokes::get_stress_tensor()
 *
 * @pre @p dir has to be a unit vector
 * @pre The conservative state stored in @p cav must have non-zero density
 */
void NavierStokes::get_dif_flux(
    const cavars &cav, const dealii::Tensor<1,dim> &dir, state &f
)
{
    const state &cons = cav.get_state();
    const avars &av = cav.get_avars();
    dealii::SymmetricTensor<2,dim> st;
    get_stress_tensor(av, st);
    
    dealii::Tensor<1,dim> mom_flux = st*dir;
    f[0] = 0; // density flux
    f[4] = 0; // initialise energy flux
    double v; // temporary quantity
    for(int d=0; d<dim; d++){
        v = cons[1+d]/cons[0];
        f[1+d] = mom_flux[d];
        f[4] += av[6+d]*dir[d] + v*mom_flux[d];
    }
    
}



// # # # # # # # # # Private Functions # # # # # # # # # # # #



/**
 * @brief HLLC x-flux function.
 *
 * Prefix 'l' and 'r' for left and right. 'c' for conservative. Assumes the interface between
 * @p lcs and @p rcs has normal in x-direction. See Toro 3rd ed, secs 3.2.4, 10.4 and 10.5. The
 * choice of left and right wave speeds here is as follows
 * @f[
 * S_L = \tilde{u} - \tilde{a}
 * S_L = \tilde{u} + \tilde{a}
 * @f]
 * where @f$\tilde{\cdot}@f$ represents Roe average.
 *
 * @pre @p lcs and @p rcs must pass the test of NavierStokes::assert_positivity()
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
 * Similar to NavierStokes::hllc_xflux(). See Toro 3rd ed, sec 10.5.1.
 *
 * @pre @p lcs and @p rcs must pass the test of NavierStokes::assert_positivity()
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
 * See eqs (3.16, 3.18-3.20) of Gassner, Winters & Kopriva (2016).
 *
 * @pre @p dir has to be a unit vector
 * @pre @p cs1 and @p cs2 must pass the test of NavierStokes::assert_positivity()
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



/**
 * @brief BR1 surface and volume flux for auxiliary variables
 *
 * Simply calculates the average of states. See the class documentation for more details.
 *
 * @note Positivity of @p cs1 and @p cs2 is not checked
 */
void NavierStokes::br1_flux(const state &cs1, const state &cs2, state &f)
{
    for(cvar var: cvar_list) f[var] = 0.5*(cs1[var] + cs2[var]);
}



/**
 * @brief BR1 surface and volume flux for diffusive flux in NS equations
 *
 * Simply returns the average of fluxes based on @p cav1 and @p cav2. Internally uses
 * NavierStokes::get_dif_flux()
 *
 * @pre @p dir must be unit vector
 * @pre Conservative states associated with @p cav1 and @p cav2 must have positive density
 */
void NavierStokes::br1_flux(
    const cavars &cav1, const cavars &cav2, const dealii::Tensor<1,dim> &dir, state &f
)
{
    state f1, f2;
    get_dif_flux(cav1, dir, f1);
    get_dif_flux(cav2, dir, f2);
    for(cvar var: cvar_list) f[var] = 0.5*(f1[var] + f2[var]);
}



/* ------------------------------------------------------------------------------------ */



#ifdef DEBUG
void NavierStokes::test()
{
    utilities::Testing t("NavierStokes", "class");
    
    {
        t.new_block();
        NavierStokes ns("air");
        ns.print_modelling_params();
    }
    
    {
        t.new_block();
        NavierStokes ns("N2");
        ns.print_modelling_params();
    }
    
    {
        t.new_block();
        NavierStokes ns(1,2,3,4,5,6);
        ns.print_modelling_params();
    }
    
    {
        t.new_block("testing get_inv_flux()");
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
        t.new_block("testing hllc_xflux() and rusanov_xflux()");
        NavierStokes ns("air");
        state lcs = {1.5,3,1.5,4.5,23}, rcs = {2,2,4,4,34}, f; // from WJ-02-Mar-2021
        ns.hllc_xflux(lcs, rcs, f);
        std::cout << "HLLC x flux";
        utilities::print_state(f);
        
        ns.rusanov_xflux(lcs,rcs,f);
        std::cout << "Rusanov x flux";
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing get_inv_surf_flux() and get_inv_vol_flux()");
        NavierStokes ns("air");
        
        std::cout << "\nTrivial case: normal along x dir\n";
        state ocs = {1.5,3,1.5,4.5,23}, ncs = {2,2,4,4,34}, f; // from previous block
        dealii::Tensor<1,dim> dir({1,0,0});
        ns.get_inv_surf_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid surface flux";
        utilities::print_state(f); // should be same as HLLC output of previous block
        ns.get_inv_vol_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid volume flux";
        utilities::print_state(f);
        
        std::cout << "\nTrivial case: normal along -x dir\n";
        ocs = {1.5,-3,-1.5,-4.5,23};
        ncs = {2,-2,-4,-4,34};
        dir[0] = -1; dir[1] = 0; dir[2] = 0; // -x dir
        ns.get_inv_surf_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid surface flux";
        utilities::print_state(f); // should be same as HLLC output of previous block, but with
                                   // velocities reversed
        ns.get_inv_vol_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid volume flux";
        utilities::print_state(f);
        
        std::cout << "\nNormal along y dir, with velocities rotated by 90 degrees ccw in xy plane\n";
        ocs = {1.5,-1.5,3,4.5,23};
        ncs = {2,-4,2,4,34}; // from previous case, rotate x & y
        dir[0] = 0; dir[1] = 1; dir[2] = 0; // y dir
        ns.get_inv_surf_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid surface flux";
        utilities::print_state(f); // should be equal to HLLC flux printed in previous block, with
                                   // x & y components rotated
        ns.get_inv_vol_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid volume flux";
        utilities::print_state(f);
        
        std::cout << "\nNormal along z dir, with velocities rotated by 90 degrees cw in xz plane\n";
        ocs = {1.5,-4.5,1.5,3,23};
        ncs = {2,-4,4,2,34}; // from previous case, rotate x & z
        dir[0] = 0; dir[1] = 0; dir[2] = 1; // z dir
        ns.get_inv_surf_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid surface flux";
        utilities::print_state(f); // should be equal to HLLC flux printed in previous block, with
                                   // x & z components rotated
        ns.get_inv_vol_flux(ocs, ncs, dir, f);
        std::cout << "Inviscid volume flux";
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing aux surf and vol functions");
        NavierStokes ns("air");
        state cs1={3,5,7,9,11}, cs2={33,55,77,99,1111}, f;
        dealii::Tensor<1,3> dir({10,20,30}); // invalid, but doesn't matter for BR1 scheme
        ns.get_aux_surf_flux(cs1, cs2, dir, f);
        utilities::print_state(f);
        ns.get_aux_vol_flux(cs1, cs2, dir, f);
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing get_stress_tensor() and get_dif_flux()");
        NavierStokes ns("air");
        state cs = {2,4,6,8,9}, f;
        avars av = {2,3,4,5,6,7,8,9,10};
        cavars cav(&cs, &av);
        dealii::Tensor<1,dim> dir({0,0,1});
        dealii::SymmetricTensor<2,dim> st;
        ns.get_stress_tensor(av, st);
        std::cout << "Stress tensor:\n" << st << "\n";
        ns.get_dif_flux(cav, dir, f);
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing get_dif_surf/vol_flux()");
        NavierStokes ns("air");
        state s1 = {2,4,6,8,9}, s2 = {2,6,4,8,10}, f;
        avars av1 = {2,3,4,5,6,7,8,9,10}, av2 = {12,13,14,15,16,17,18,19,110};
        cavars cav1(&s1, &av1), cav2(&s2, &av2);
        dealii::Tensor<1,dim> dir({0,0,1});
        ns.get_dif_surf_flux(cav1, cav2, dir, f);
        utilities::print_state(f);
        ns.get_dif_vol_flux(cav1, cav2, dir, f);
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

