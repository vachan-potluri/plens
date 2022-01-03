/**
 * @file navier_stokes.cc
 * @brief Class for navier stokes solver
 */

#include "navier_stokes.h"



/**
 * @brief Constructor. Set parameter values explicitly.
 *
 * Internally calls NavierStokes::set_modelling_params(), and other flux scheme setters.
 */
NavierStokes::NavierStokes(
    const double gma, const double M, const double Pr, double mu0, const double T0, const double S,
    const aux_surf_flux_scheme asfs,
    const aux_vol_flux_scheme avfs,
    const inv_surf_flux_scheme isfs,
    const inv_vol_flux_scheme ivfs,
    const dif_surf_flux_scheme dsfs,
    const dif_vol_flux_scheme dvfs
): flux_blender_value(1)
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
 * Internally calls NavierStokes::set_modelling_params(), and other flux scheme setters. If
 * `inviscid==true`, then @f$\mu_0@f$ is set to zero. This results in zero values of @f$\mu@f$
 * and @f$k@f$.
 */
NavierStokes::NavierStokes(
    const std::string gas_name,
    const bool inviscid,
    const aux_surf_flux_scheme asfs,
    const aux_vol_flux_scheme avfs,
    const inv_surf_flux_scheme isfs,
    const inv_vol_flux_scheme ivfs,
    const dif_surf_flux_scheme dsfs,
    const dif_vol_flux_scheme dvfs
): flux_blender_value(1)
{
    bool gas_supported = (gas_name=="air" || gas_name=="N2" || gas_name=="nitrogen");
    AssertThrow(
        gas_supported,
        dealii::StandardExceptions::ExcMessage(
            "Unsupported gas name. Only 'air' and 'N2' (or 'nitrogen') are currently supported"
        )
    );
    
    double gma=1.4, M, Pr=0.69, mu0(0), T0=273, S;
    if(gas_name == "air"){
        M = 0.029;
        if(!inviscid) mu0 = 1.716e-5;
        S = 111;
    }
    else{
        // guaranteed to be nitrogen
        M = 0.028;
        if(!inviscid) mu0 = 1.663e-5;
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
        const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f){
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
        const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f){
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
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->hllc_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::rusanov){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->rusanov_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::ausm_plus_up){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->ausm_plus_up_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::rusanov_hllc_blend){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->rusanov_hllc_blend_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::rusanov_ausm_plus_up_blend){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->rusanov_ausm_plus_up_blend_xflux(lcs, rcs, f);
        };
    }
    else{
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->modified_sw_xflux(lcs, rcs, f);
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
        const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f){
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
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
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
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
    ){
        this->br1_flux(cav1, cav2, dir, f);
    };
}



/**
 * @brief Sets surface, volume and internal (or, exact) flux wrappers. Surface flux wrappers are
 * used in assembling surface flux for the 3 stages. Volume and internal/exact flux wrappers are
 * used for calculating residual (or, RHS) for stages 2 and 3. See PLENS class documentation for
 * details.
 *
 * See the class documentation for more details.
 *
 * @pre All the surf/vol setters must be called before invoking this function
 */
void NavierStokes::set_wrappers()
{
    // surface flux wrappers
    surf_flux_wrappers[0] = [=](
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
    ){
        this->get_aux_surf_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // aux
    
    surf_flux_wrappers[1] = [=](
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
    ){
        this->get_inv_surf_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // inv
    
    surf_flux_wrappers[2] = get_dif_surf_flux; // dif, directly equate function objects

    // volume flux wrappers
    vol_flux_wrappers[0] = [=](
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
    ){
        this->get_aux_vol_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // aux
    
    vol_flux_wrappers[1] = [=](
        const CAvars &cav1, const CAvars cav2, const dealii::Tensor<1,dim> &dir, State &f
    ){
        this->get_inv_vol_flux(cav1.get_state(), cav2.get_state(), dir, f);
    }; // inv
    
    vol_flux_wrappers[2] = get_dif_vol_flux; // dif, directly equate function objects

    // internal flux wrappers
    flux_wrappers[0] = [=](
        const CAvars &cav, const dealii::Tensor<1,dim>& dir, State& flux
    ){
        flux = cav.get_state();
    }; // aux

    flux_wrappers[1] = [=](
        const CAvars &cav, const dealii::Tensor<1,dim>& dir, State& flux
    ){
        this->get_inv_flux(cav.get_state(), dir, flux);
    }; // inv

    flux_wrappers[2] = get_dif_flux; // dif, equate functions directly
}



/**
 * @brief Asserts positivity of density and thermal energy of given state
 */
void NavierStokes::assert_positivity(const State &cons)
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
double NavierStokes::get_e(const State &cons)
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
double NavierStokes::get_p(const State &cons) const
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
double NavierStokes::get_a(const State &cons) const
{
    return sqrt(gma_*get_p(cons)/cons[0]);
}



/**
 * Get the Mach number. Internally calls NavierStokes::get_a()
 *
 * @pre @p cons must pass the test of NavierStokes::assert_positivity()
 */
double NavierStokes::get_M(const State &cons) const
{
    double speed_sq(0);
    for(int d=0; d<dim; d++) speed_sq += (cons[1+d]*cons[1+d])/(cons[0]*cons[0]);
    return sqrt(speed_sq)/get_a(cons);
}



/**
 * @brief Calculates inviscid flux in a given direction based on given conservative state
 * 
 * @pre @p cons must pass the test of NavierStokes::assert_positivity()
 * @pre @p dir has to be a unit vector
 */
void NavierStokes::get_inv_flux(
    const State &cons, const dealii::Tensor<1,dim> &dir, State &f
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
 * for this function is described in detail in WJ-23-Feb-2021 and WJ-01-Jun-2021. Only the outline
 * is noted here.
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
 * See also WJ-01-Jun-2021 and WJ-31-May-2021. Initially, `asin` of the cross product was used to
 * calculate the angle between vectors. However, the return range of `asin` is
 * @f$[-\pi/2,\pi/2]@f$ which is not enough. We require angle to lie in @f$[0,\pi]@f$. The angle
 * cannot be greated than @f$\pi@f$ because in that case, the axis of rotation
 * (@f$\hat{x} \times \hat{n}@f$) would also reverse its direction. So when viewed from the axis of
 * rotation so computed, the rotation is always anti-cloclwise and the rotation angle is always in
 * the range @f$[0,\pi]@f$. So keeping all this in mind, `acos` of the dot product is the best
 * option.
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
    const State &ocs, const State &ncs, const dealii::Tensor<1,dim> &normal, State &f
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
        double theta = acos(
            dealii::scalar_product(xdir, normal)
        ); // <-- rotation angle (both xdir and normal are unit vectors)
        
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
    
    State ocs_r, ncs_r; // rotated states
    ocs_r[0] = ocs[0];
    ncs_r[0] = ncs[0];
    for(int d=0; d<dim; d++){
        ocs_r[1+d] = osmom_r[d];
        ncs_r[1+d] = nsmom_r[d];
    }
    ocs_r[4] = ocs[4];
    ncs_r[4] = ncs[4];
    
    // Step 3: get normal flux wrt rotated coordinate system
    State f_r;
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
void NavierStokes::get_stress_tensor(const Avars &av, dealii::SymmetricTensor<2,dim> &st)
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
    const CAvars &cav, const dealii::Tensor<1,dim> &dir, State &f
)
{
    const State &cons = cav.get_state();
    const Avars &av = cav.get_avars();
    dealii::SymmetricTensor<2,dim> st;
    get_stress_tensor(av, st);
    
    dealii::Tensor<1,dim> mom_flux = st*dir;
    f[0] = 0; // density flux
    f[4] = 0; // initialise energy flux
    double v; // temporary quantity
    for(int d=0; d<dim; d++){
        v = cons[1+d]/cons[0];
        f[1+d] = mom_flux[d];
        f[4] += -av[6+d]*dir[d] + v*mom_flux[d];
    }
    
}



/**
 * Gives the x-directional right eigen vector matrix as @p K. The matrix itself is nothing but a
 * horizontal concatenation of right eigen vectors related to x-directional inviscid flux. See
 * Toro, 3rd ed, sec 3.2.2.
 *
 * @param[in] vel Velocity vector
 * @param[in] a Speed of sound
 * @param[in] H Total enthalpy (@f$ e + \frac{\vec{u} \cdot \vec{u}}{2} + \frac{p}{\rho} @f$)
 * @param[out] K The right eigen vector matrix
 *
 * @pre @p K must be a square matrix of size 5. No assertions on this are made. So the behaviour
 * maybe unexpected when this condition is not met.
 *
 * @note Although in principle this function can be static, it is not being done because get_xKinv
 * cannot be static.
 */
void NavierStokes::get_xK(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    dealii::FullMatrix<double> &K
) const
{
    // 1st col
    K(0,0) = 1;
    for(int d=0; d<dim; d++) K(1+d, 0) = vel[d];
    K(1,0) -= a;
    K(4,0) = H - vel[0]*a;

    // 2nd col
    K(0,1) = 1;
    for(int d=0; d<dim; d++) K(1+d, 1) = vel[d];
    K(4,1) = 0.5*dealii::scalar_product(vel, vel);

    // 3rd col
    K(0,2) = 0;
    K(1,2) = 0;
    K(2,2) = 1;
    K(3,2) = 0;
    K(4,2) = vel[1];

    // 4th col
    K(0,3) = 0;
    K(1,3) = 0;
    K(2,3) = 0;
    K(3,3) = 1;
    K(4,3) = vel[2];

    // 5th col
    K(0,4) = 1;
    for(int d=0; d<dim; d++) K(1+d, 4) = vel[d];
    K(1,4) += a;
    K(4,4) = H + vel[0]*a;
}



/**
 * Gives the inverse of x-directional right eigen vector matrix as @p Kinv. All other details and
 * requirements are exactly as for get_xK(). This function depends on NavierStoes::gma_ and thus
 * cannot be static.
 */
void NavierStokes::get_xKinv(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    dealii::FullMatrix<double> &Kinv
) const
{
    // 1st row
    Kinv(0,0) = H + a*(vel[0]-a)/(gma_-1);
    for(int d=0; d<dim; d++) Kinv(0,1+d) = -vel[d];
    Kinv(0,1) -= a/(gma_-1);
    Kinv(0,4) = 1;

    // 2nd row
    Kinv(1,0) = 2*(2*a*a/(gma_-1) - H);
    for(int d=0; d<dim; d++) Kinv(1,1+d) = 2*vel[d];
    Kinv(1,4) = -2;

    // 3rd row
    Kinv(2,0) = -2*vel[1]*a*a/(gma_-1);
    Kinv(2,1) = 0;
    Kinv(2,2) = 2*a*a/(gma_-1);
    Kinv(2,3) = 0;
    Kinv(2,4) = 0;

    // 4th row
    Kinv(3,0) = -2*vel[2]*a*a/(gma_-1);
    Kinv(3,1) = 0;
    Kinv(3,2) = 0;
    Kinv(3,3) = 2*a*a/(gma_-1);
    Kinv(3,4) = 0;

    // 5th row
    Kinv(4,0) = H - a*(vel[0]+a)/(gma_-1);
    for(int d=0; d<dim; d++) Kinv(4,1+d) = -vel[d];
    Kinv(4,1) += a/(gma_-1);
    Kinv(4,4) = 1;

    const double factor = 0.5*(gma_-1)/(a*a);
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            Kinv(i,j) *= factor;
        }
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
void NavierStokes::hllc_xflux(const State &lcs, const State &rcs, State &f) const
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
        State lss = {
            temp*lcs[0],
            temp*lcs[0]*s,
            temp*lcs[2],
            temp*lcs[3],
            temp*( lcs[4] + (s-ul)*( lcs[0]*s + pl/(sl-ul) ))
        }; // left star state
        
        State lf; // flux based on lcs
        get_inv_flux(lcs, xdir, lf);
        
        for(cvar var: cvar_list) f[var] = lf[var] + sl*(lss[var] - lcs[var]);
    }
    else if(sr>0){
        // right star state at interface
        double temp = (sr-ur)/(sr-s);
        State rss = {
            temp*rcs[0],
            temp*rcs[0]*s,
            temp*rcs[2],
            temp*rcs[3],
            temp*( rcs[4] + (s-ur)*( rcs[0]*s + pr/(sr-ur) ))
        }; // right star state
        
        State rf; // flux based on rcs
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
void NavierStokes::rusanov_xflux(const State &lcs, const State &rcs, State &f) const
{
    dealii::Tensor<1,dim> xdir({1,0,0});
    
    State lf, rf; // left and right conservative fluxes
    double al = get_a(lcs), ar = get_a(rcs); // left and right sound speeds
    double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0]; // left & right flow speeds
    double S; // the "single wave" speed
    
    get_inv_flux(lcs, xdir, lf);
    get_inv_flux(rcs, xdir, rf);
    
    S = std::max(fabs(ul)+al, fabs(ur)+ar);
    
    for(cvar var: cvar_list) f[var] = 0.5*(lf[var] + rf[var]) - 0.5*S*(rcs[var] - lcs[var]);
}



/**
 * @brief AUSM+-up flux function
 *
 * See Liou (2006), JCP, or Anant's thesis.
 *
 * @pre @p lcs and @p rcs must pass the test of NavierStokes::assert_positivity()
 */
void NavierStokes::ausm_plus_up_xflux(const State &lcs, const State &rcs, State &f) const
{
    double pl = get_p(lcs), pr = get_p(rcs); // pressures
    double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0]; // velocities
    double Hl = (lcs[4]+pl)/lcs[0], Hr = (rcs[4]+pr)/rcs[0]; // total/stagnation enthalpies

    double astl = sqrt(2*(gma_-1)/(gma_+1)*Hl),
        astr = sqrt(2*(gma_-1)/(gma_+1)*Hr); // critical speed of sounds ('st'ar)
    double ahl = astl*astl/std::max(astl, ul),
        ahr = astr*astr/std::max(astr, -ur); // 'h'at velocities (eq 30 in Liou 2006)
    double a12 = std::min(ahl, ahr); // a_1/2

    double Ml = ul/a12, Mr = ur/a12; // Mach numbers
    double Mavgsq = 0.5*(Ml*Ml + Mr*Mr);

    double M12 = ausm::mach_split_4_pos(Ml) + ausm::mach_split_4_neg(Mr) -
        ausm::Kp/ausm::fa*std::max(1-ausm::sigma*Mavgsq, 0.0)*2*(pr-pl)/(a12*a12*(lcs[0]+rcs[0]));

    double p_split_pos = ausm::pressure_split_5_pos(Ml),
        p_split_neg = ausm::pressure_split_5_neg(Mr);
    
    double p12 = p_split_pos*pl + p_split_neg*pr -
        ausm::Ku*ausm::fa*p_split_pos*p_split_neg*(lcs[0]+rcs[0])*a12*(ur-ul);
    
    // set the final flux
    if(M12 > 0){
        double m12 = a12*M12*lcs[0]; // mass flow rate
        f[0] = m12;
        f[1] = m12*ul + p12;
        f[2] = m12*lcs[2]/lcs[0];
        f[3] = m12*lcs[3]/lcs[0];
        f[4] = m12*Hl;
    }
    else{
        double m12 = a12*M12*rcs[0]; // mass flow rate
        f[0] = m12;
        f[1] = m12*ur + p12;
        f[2] = m12*rcs[2]/rcs[0];
        f[3] = m12*rcs[3]/rcs[0];
        f[4] = m12*Hr;
    }
}



/**
 * Blended Rusanov-HLLC flux function. Rusanov flux is given weight NavierStokes::flux_blender_value
 * and HLLC is given the complementary weight.
 *
 * @pre The value of NavierStokes::flux_blender_value must be appropriately set using
 * set_flux_blender_value() before using this function.
 */
void NavierStokes::rusanov_hllc_blend_xflux(const State &lcs, const State &rcs, State &f) const
{
    State f_rusanov, f_hllc;
    rusanov_xflux(lcs, rcs, f_rusanov);
    hllc_xflux(lcs, rcs, f_hllc);
    for(cvar var: cvar_list) f[var] = flux_blender_value*f_rusanov[var] +
        (1-flux_blender_value)*f_hllc[var];
}



/**
 * Blended Rusanov-AUSM+-up flux function. Rusanov gets weight NavierStokes::flux_blender_value, and
 * AUSM+-up gets the complement.
 *
 * @pre The value of NavierStokes::flux_blender_value must be appropriately set using
 * set_flux_blender_value() before using this function.
 */
void NavierStokes::rusanov_ausm_plus_up_blend_xflux(const State &lcs, const State &rcs, State &f)
const
{
    State f_rusanov, f_ausm;
    rusanov_xflux(lcs, rcs, f_rusanov);
    ausm_plus_up_xflux(lcs, rcs, f_ausm);
    for(cvar var: cvar_list) f[var] = flux_blender_value*f_rusanov[var] +
        (1-flux_blender_value)*f_ausm[var];
}



/**
 * Modified Steger-Warming flux. See Druguet, Candler & Nompelis (2005). Algo:
 * - Calculate omega based on left and right pressure values
 * - Using omega, calculate pos and neg states based on which Jacobian matrices will be constructed
 * - Calculate pos and neg eigen values with epsilon-based correction
 * - Using pos and neg states, and using pos and neg eigen values, calculate the pos and neg
 *   Jacobian matrices using get_xK() and get_xKinv()
 * - Calculate the final flux
 */
void NavierStokes::modified_sw_xflux(const State &lcs, const State &rcs, State &f) const
{
    // Calculate omega
    const double pl = get_p(lcs), pr = get_p(rcs);
    const double sigma2 = 0.5;
    // const double omega = 0.5/( 1 + std::pow(sigma2*(pr-pl)/std::min(pl,pr), 2) );
    const double omega = 0;

    // pos and neg states which are inputs for calculating pos and neg jacobians
    // these are also used to calculate pos and neg eigenvalues
    dealii::Tensor<1,dim> vel_pos, vel_neg;
    for(int d=0; d<dim; d++){
        vel_pos[d] = (1-omega)*lcs[1+d]/lcs[0] + omega*rcs[1+d]/rcs[0];
        vel_neg[d] = (1-omega)*rcs[1+d]/rcs[0] + omega*lcs[1+d]/lcs[0];
    }
    std::cout << "vel pos: " << vel_pos << "\nvel neg: " << vel_neg << "\n";

    const double rho_pos = (1-omega)*lcs[0] + omega*rcs[0],
        rho_neg = (1-omega)*rcs[0] + omega*lcs[0];
    const double p_pos = (1-omega)*pl + omega*pr,
        p_neg = (1-omega)*pr + omega*pl;
    const double a_pos = std::sqrt(gma_*p_pos/rho_pos),
        a_neg = std::sqrt(gma_*p_neg/rho_neg);
    const double H_pos = gma_*p_pos/((gma_-1)*rho_pos) + dealii::scalar_product(vel_pos, vel_pos),
        H_neg = gma_*p_neg/((gma_-1)*rho_neg) + dealii::scalar_product(vel_neg, vel_neg);
    std::cout << "rho pos and neg: " << rho_pos << " " << rho_neg << "\n";
    std::cout << "p pos and neg: " << p_pos << " " << p_neg << "\n";

    // pos and neg eigenvector matrices
    dealii::FullMatrix<double> K_pos(dim+2), Kinv_pos(dim+2), K_neg(dim+2), Kinv_neg(dim+2);
    get_xK(vel_pos, a_pos, H_pos, K_pos);
    get_xK(vel_neg, a_neg, H_neg, K_neg);
    get_xKinv(vel_pos, a_pos, H_pos, Kinv_pos);
    get_xKinv(vel_neg, a_neg, H_neg, Kinv_neg);

    // pos and neg eigenvalues, and their correction
    std::array<double, dim+2> eig_pos, eig_neg;
    eig_pos[0] = pos(vel_pos[0] - a_pos);
    for(int d=0; d<dim; d++) eig_pos[1+d] = pos(vel_pos[0]);
    eig_pos[4] = pos(vel_pos[0] + a_pos);
    eig_neg[0] = neg(vel_neg[0] - a_neg);
    for(int d=0; d<dim; d++) eig_neg[1+d] = neg(vel_neg[0]);
    eig_neg[4] = neg(vel_neg[0] + a_neg);
    std::cout << "eig pos: ";
    utilities::print_state(eig_pos);
    std::cout << "eig neg: ";
    utilities::print_state(eig_neg);

    // const double eps = 0.3;
    const double eps = 0;
    for(int i=0; i<dim+2; i++){
        eig_pos[i] = 0.5*(eig_pos[i] + sqrt(pow(eig_pos[i], 2) + pow(eps*a_pos, 2)));
        eig_neg[i] = 0.5*(eig_neg[i] - sqrt(pow(eig_neg[i], 2) + pow(eps*a_neg, 2)));
    }

    dealii::FullMatrix<double> A_pos(dim+2), A_neg(dim+2);
    A_pos = 0;
    A_neg = 0;
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            for(int k=0; k<dim+2; k++){
                A_pos(i,j) += K_pos(i,k)*eig_pos[k]*Kinv_pos(k,j);
                A_neg(i,j) += K_neg(i,k)*eig_neg[k]*Kinv_neg(k,j);
            }
        }
    }
    std::cout << "\nK pos:" << "\n";
    K_pos.print_formatted(std::cout);
    std::cout << "\nKinv pos:" << "\n";
    Kinv_pos.print_formatted(std::cout);
    std::cout << "\nK neg:" << "\n";
    K_neg.print_formatted(std::cout);
    std::cout << "\nKinv neg:" << "\n";
    Kinv_neg.print_formatted(std::cout);
    std::cout << "\nA pos:" << "\n";
    A_pos.print_formatted(std::cout);
    std::cout << "\nA neg:" << "\n";
    A_neg.print_formatted(std::cout);

    // calculate the flux
    for(int i=0; i<dim+2; i++){
        f[i] = 0;
        for(int j=0; j<dim+2; j++){
            f[i] += A_pos(i,j)*lcs[j] + A_neg(i,j)*rcs[j];
        }
    }
}



/**
 * @brief Chandrashekhar inviscid volume flux.
 *
 * See eqs (3.16, 3.18-3.20) of Gassner, Winters & Kopriva (2016).
 *
 * On later testing, it was identified that the quantities @f$\rho^{\text{ln}}@f$ and
 * @f$\beta^\text{ln}@f$ can cause problem if both the states have same values
 * (see WJ-31-May-2021). So for such "ln" quantities, if the denominators are less that 1e-8, the
 * "ln" quantities are set directly based on state 1.
 *
 * @pre @p dir has to be a unit vector
 * @pre @p cs1 and @p cs2 must pass the test of NavierStokes::assert_positivity()
 */
void NavierStokes::chandrashekhar_flux(
    const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
) const
{
    double p1 = get_p(cs1), p2 = get_p(cs2);
    double beta1 = 0.5*cs1[0]/p1, beta2 = 0.5*cs2[0]/p2;
    double beta_ln(beta1);
    double denom = log(beta1) - log(beta2);
    if(fabs(denom) > 1e-8) beta_ln = (beta1-beta2)/(denom);
    
    double p_hat = 0.5*(cs1[0]+cs2[0])/(beta1+beta2);
    double rho_ln(cs1[0]);
    denom = log(cs1[0]) - log(cs2[0]);
    if(fabs(denom) > 1e-8) rho_ln = (cs1[0]-cs2[0])/(denom);
    
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
void NavierStokes::br1_flux(const State &cs1, const State &cs2, State &f)
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
    const CAvars &cav1, const CAvars &cav2, const dealii::Tensor<1,dim> &dir, State &f
)
{
    State f1, f2;
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
        const double T = 400;
        const double mu = ns.get_mu(T);
        std::cout << "Viscosity and thermal conductivity at " << T << " K: "
            << mu << ", " << ns.get_k(mu) << "\n";
    }
    
    {
        t.new_block();
        NavierStokes ns("N2");
        ns.print_modelling_params();
    }

    {
        t.new_block("testing inviscid argument");
        NavierStokes ns("N2", true); // inviscid
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
        State cons = {2,2,4,6,15};
        std::cout << "Pressure " << ns.get_p(cons) << "\n";
        std::cout << "Energy " << ns.get_e(cons) << "\n";
        std::cout << "Mach number " << ns.get_M(cons) << "\n";
        ns.assert_positivity(cons);
        
        std::array<State, 3> fluxes;
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
        State lcs = {1.5,3,1.5,4.5,23}, rcs = {2,2,4,4,34}, f; // from WJ-02-Mar-2021
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
        State ocs = {1.5,3,1.5,4.5,23}, ncs = {2,2,4,4,34}, f; // from previous block
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
        State cs1={3,5,7,9,11}, cs2={33,55,77,99,1111}, f;
        dealii::Tensor<1,3> dir({10,20,30}); // invalid, but doesn't matter for BR1 scheme
        ns.get_aux_surf_flux(cs1, cs2, dir, f);
        utilities::print_state(f);
        ns.get_aux_vol_flux(cs1, cs2, dir, f);
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing get_stress_tensor() and get_dif_flux()");
        NavierStokes ns("air");
        State cs = {2,4,6,8,9}, f;
        Avars av = {2,3,4,5,6,7,8,9,10};
        CAvars cav(&cs, &av);
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
        State s1 = {2,4,6,8,9}, s2 = {2,6,4,8,10}, f;
        Avars av1 = {2,3,4,5,6,7,8,9,10}, av2 = {12,13,14,15,16,17,18,19,110};
        CAvars cav1(&s1, &av1), cav2(&s2, &av2);
        dealii::Tensor<1,dim> dir({0,0,1});
        ns.get_dif_surf_flux(cav1, cav2, dir, f);
        utilities::print_state(f);
        ns.get_dif_vol_flux(cav1, cav2, dir, f);
        utilities::print_state(f);
    }
    
    {
        t.new_block("testing surface flux wrappers");
        NavierStokes ns("air");
        State s1 = {1.5,-4.5,1.5,3,23}, s2 = {2,-4,4,2,34}, f;
        Avars av1 = {2,3,4,5,6,7,8,9,10}, av2 = {12,13,14,15,16,17,18,19,110};
        CAvars cav1(&s1, &av1), cav2(&s2, &av2);
        dealii::Tensor<1,dim> dir({0,0,1});
        for(int stage=0; stage<3; stage++){
            ns.surf_flux_wrappers[stage](cav1, cav2, dir, f);
            std::cout << "Stage " << stage << " flux:";
            utilities::print_state(f);
        }
    }

    {
        t.new_block("testing ausm_plus_up_xflux()");
        NavierStokes ns("air");
        State lcs = {1,0,0,0,1000/0.4}, rcs = {1,0,0,0,0.01/0.4}, f;

        ns.ausm_plus_up_xflux(lcs, rcs, f);
        std::cout << "AUSM+-up x flux: ";
        utilities::print_state(f);

        ns.hllc_xflux(lcs, rcs, f);
        std::cout << "HLLC x flux: ";
        utilities::print_state(f);
    }

    {
        t.new_block("Testing modified SW flux");
        NavierStokes ns("air");
        State lcs = {1,0.75,0,0,2.78125}, f;

        ns.modified_sw_xflux(lcs, lcs, f);
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

