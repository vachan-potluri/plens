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
    else if(isfs == inv_surf_flux_scheme::modified_sw){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->modified_sw_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::chandrashekhar){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->chandrashekhar_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::kennedy_gruber){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->kennedy_gruber_xflux(lcs, rcs, f);
        };
    }
    else if(isfs == inv_surf_flux_scheme::rusanov_kennedy_gruber_blend){
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->rusanov_kennedy_gruber_blend_xflux(lcs, rcs, f);
        };
    }
    else{
        inv_surf_xflux = [=](const State &lcs, const State &rcs, State &f){
            this->ismail_roe_xflux(lcs, rcs, f);
        };
    }
}



/**
 * @brief Sets the invsicid volume flux function NavierStokes::inv_vol_flux
 */
void NavierStokes::set_inv_vol_flux_scheme(const inv_vol_flux_scheme ivfs)
{
    if(ivfs == inv_vol_flux_scheme::chandrashekhar){
        get_inv_vol_flux = [=](
            const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
        ){
            this->chandrashekhar_vol_flux(cs1, cs2, dir, f);
        };
    }
    else if(ivfs == inv_vol_flux_scheme::kennedy_gruber){
        get_inv_vol_flux = [=](
            const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
        ){
            this->kennedy_gruber_vol_flux(cs1, cs2, dir, f);
        };
    }
    else{
        get_inv_vol_flux = [=](
            const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
        ){
            this->ismail_roe_vol_flux(cs1, cs2, dir, f);
        };
    }
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
        ske += cons[1+dir]*cons[1+dir];
    }
    const double rho_inv = 1/cons[0];
    ske *= 0.5*rho_inv;
    
    return (cons[4]-ske)*rho_inv;
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
        ske += cons[1+dir]*cons[1+dir];
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
    const double rho_inv = 1/cons[0];
    for(int d=0; d<dim; d++) speed_sq += (cons[1+d]*cons[1+d])*(rho_inv*rho_inv);
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
    const double rhoinv = 1.0/cons[0];
    for(int d=0; d<dim; d++){
        vel[d] = cons[1+d]*rhoinv;
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
    
    dealii::Tensor<2,dim> R; // rotation matrix, all entries are zero
    if(M > 1e-3){
        // this tolerance corrsponds to an angle of ~0.06 degrees between x and n
        m /= M; // now m is a unit vector <-- rotation axis
        const double c = dealii::scalar_product(xdir, normal), // cos of angle
            t = 1-c, // 1 - (cos of angle)
            s = std::sqrt(1-c*c); // sin of angle, will be +ve since angles lies in [0,pi]
        R = dealii::Tensor<2,dim>(
            {{t * m[0] * m[0] + c,
            t * m[0] * m[1] - s * m[2],
            t * m[0] * m[2] + s * m[1]},
            {t * m[0] * m[1] + s * m[2],
            t * m[1] * m[1] + c,
            t * m[1] * m[2] - s * m[0]},
            {t * m[0] * m[2] - s * m[1],
            t * m[1] * m[2] + s * m[0],
            t * m[2] * m[2] + c}}
        ); // taken from <deal.II/physics/transformations.h>
    }
    else{
        // either x is parallel or anti-parallel to n
        if(normal[0] < 0){
            // n is anti-parallel to x
            for(int d=0; d<dim; d++) R[d][d] = -1;
        }
        else{
            // n is parallel to x
            for(int d=0; d<dim; d++) R[d][d] = 1;
        }
    }
    dealii::Tensor<2,dim> RT(dealii::transpose(R)); // transpose of R
    
    // Step 2: get rotated states
    dealii::Tensor<1,dim> osmom, nsmom, // owner and neighbor specific momentum
        osmom_r, nsmom_r; // rotated specific momentum (initialised to 0)
    for(int d=0; d<dim; d++){
        osmom[d] = ocs[1+d];
        nsmom[d] = ncs[1+d];
    }
    // get the momentum components wrt rotated coordinate system
    // osmom_r = R^{-1} * osmom, R^T = R^{-1}
    for(int row=0; row<dim; row++){
        for(int col=0; col<dim; col++){
            const double RT_val = RT[row][col];
            osmom_r[row] += RT_val*osmom[col];
            nsmom_r[row] += RT_val*nsmom[col];
        }
    }
    
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
    dealii::Tensor<1,dim> mom_flux_r, mom_flux; // momentum fluxes (initialised to 0)
    for(int d=0; d<dim; d++) mom_flux_r[d] = f_r[1+d];
    // get momentum flux components w.r.t original coordinate system
    // mom_flux = R * mom_flux_r
    for(int row=0; row<dim; row++){
        for(int col=0; col<dim; col++){
            mom_flux[row] += R[row][col]*mom_flux_r[col];
        }
    }
    
    f[0] = f_r[0];
    for(int d=0; d<dim; d++) f[1+d] = mom_flux[d];
    f[4] = f_r[4];
}



/**
 * @brief Gives the symmetric stress tensor based on avars provided
 */
void NavierStokes::get_stress_tensor(const Avars &av, dealii::Tensor<2,dim> &st)
{
    int i=0;
    for(int row=0; row<dim; row++){
        for(int col=row; col<dim; col++){
            st[row][col] = av[i];
            if(col != row) st[col][row] = av[i];
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
    dealii::Tensor<2,dim> st;
    get_stress_tensor(av, st);
    
    dealii::Tensor<1,dim> mom_flux;
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            mom_flux[i] += st[i][j]*dir[j];
        }
    }
    f[0] = 0; // density flux
    f[4] = 0; // initialise energy flux
    double v; // temporary quantity
    const double rhoinv = 1.0/cons[0];
    for(int d=0; d<dim; d++){
        v = cons[1+d]*rhoinv;
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
 * @note Although in principle this function can be static, it is not being done because get_xKinv
 * cannot be static.
 */
void NavierStokes::get_xK(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    dealii::Tensor<2,dim+2> &K
) const
{
    // 1st col
    K[0][0] = 1;
    for(int d=0; d<dim; d++) K[1+d][0] = vel[d];
    K[1][0] -= a;
    K[4][0] = H - vel[0]*a;

    // 2nd col
    K[0][1] = 1;
    for(int d=0; d<dim; d++) K[1+d][1] = vel[d];
    K[4][1] = 0.5*dealii::scalar_product(vel, vel);

    // 3rd col
    K[0][2] = 0;
    K[1][2] = 0;
    K[2][2] = 1;
    K[3][2] = 0;
    K[4][2] = vel[1];

    // 4th col
    K[0][3] = 0;
    K[1][3] = 0;
    K[2][3] = 0;
    K[3][3] = 1;
    K[4][3] = vel[2];

    // 5th col
    K[0][4] = 1;
    for(int d=0; d<dim; d++) K[1+d][4] = vel[d];
    K[1][4] += a;
    K[4][4] = H + vel[0]*a;
}



/**
 * Gives the right eigenvector matrix in an arbitrary direction @p dir. Similar to get_xK(), but
 * more generic. See Rohde (2001) for the formula.
 *
 * The 4th and 5th columns can get singular which may cause difficulty in inverting the matrix.
 * Thus the form (R-1), (R-2) or (R-3) is used depending on the largest component of @p dir.
 *
 * @pre @p dir must be a unit vector. No checks on this are done.
 *
 * @note Rohde's column ordering is different and some columns also have a sign change. This
 * function is in this sense not consistent with get_xK(). However, as long as this function is
 * used with get_Kinv(), the final results don't matter.
 */
void NavierStokes::get_K(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    const dealii::Tensor<1,dim>& dir,
    dealii::Tensor<2,dim+2> &K
) const
{
    const double vel_n = dealii::scalar_product(vel, dir);

    // 1st col
    K[0][0] = 1;
    for(int d=0; d<dim; d++) K[1+d][0] = vel[d] - a*dir[d];
    K[4][0] = H - vel_n*a;

    // 2nd col
    K[0][1] = 1;
    for(int d=0; d<dim; d++) K[1+d][1] = vel[d];
    K[4][1] = 0.5*dealii::scalar_product(vel, vel);

    // 3rd col
    K[0][2] = 1;
    for(int d=0; d<dim; d++) K[1+d][2] = vel[d] + a*dir[d];
    K[4][2] = H + vel_n*a;

    double dir_abs[3];
    for(int d=0; d<dim; d++) dir_abs[d] = fabs(dir[d]);

    if(dir_abs[0] > dir_abs[1] && dir_abs[0] > dir_abs[2]){
        // (R-1) form
        // 4th col
        K[0][3] = 0;
        K[1][3] = dir[1];
        K[2][3] = -dir[0];
        K[3][3] = 0;
        K[4][3] = vel[0]*dir[1] - vel[1]*dir[0];

        // 5th col
        K[0][4] = 0;
        K[1][4] = -dir[2];
        K[2][4] = 0;
        K[3][4] = dir[0];
        K[4][4] = vel[2]*dir[0] - vel[1]*dir[2];
    }
    else if(dir_abs[1] > dir_abs[0] && dir_abs[1] > dir_abs[2]){
        // (R-2) form
        // 4th col
        K[0][3] = 0;
        K[1][3] = dir[1];
        K[2][3] = -dir[0];
        K[3][3] = 0;
        K[4][3] = vel[0]*dir[1] - vel[1]*dir[0];

        // 5th col
        K[0][4] = 0;
        K[1][4] = 0;
        K[2][4] = dir[2];
        K[3][4] = -dir[1];
        K[4][4] = vel[1]*dir[2] - vel[2]*dir[1];
    }
    else{
        // (R-3) form
        // 4th col
        K[0][4] = 0;
        K[1][4] = -dir[2];
        K[2][4] = 0;
        K[3][4] = dir[0];
        K[4][4] = vel[2]*dir[0] - vel[1]*dir[2];

        // 5th col
        K[0][4] = 0;
        K[1][4] = 0;
        K[2][4] = dir[2];
        K[3][4] = -dir[1];
        K[4][4] = vel[1]*dir[2] - vel[2]*dir[1];
    }
}



/**
 * Gives the inverse of x-directional right eigen vector matrix as @p Kinv. All other details and
 * requirements are exactly as for get_xK(). This function depends on NavierStokes::gma_ and thus
 * cannot be static.
 */
void NavierStokes::get_xKinv(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    dealii::Tensor<2,dim+2> &Kinv
) const
{
    const double temp = 1/(gma_-1);
    // 1st row
    Kinv[0][0] = H + a*(vel[0]-a)*temp;
    for(int d=0; d<dim; d++) Kinv[0][1+d] = -vel[d];
    Kinv[0][1] -= a*temp;
    Kinv[0][4] = 1;

    // 2nd row
    Kinv[1][0] = 2*(2*a*a*temp - H);
    for(int d=0; d<dim; d++) Kinv[1][1+d] = 2*vel[d];
    Kinv[1][4] = -2;

    // 3rd row
    Kinv[2][0] = -2*vel[1]*a*a*temp;
    Kinv[2][1] = 0;
    Kinv[2][2] = 2*a*a*temp;
    Kinv[2][3] = 0;
    Kinv[2][4] = 0;

    // 4th row
    Kinv[3][0] = -2*vel[2]*a*a*temp;
    Kinv[3][1] = 0;
    Kinv[3][2] = 0;
    Kinv[3][3] = 2*a*a*temp;
    Kinv[3][4] = 0;

    // 5th row
    Kinv[4][0] = H - a*(vel[0]+a)*temp;
    for(int d=0; d<dim; d++) Kinv[4][1+d] = -vel[d];
    Kinv[4][1] += a*temp;
    Kinv[4][4] = 1;

    const double factor = 0.5*(gma_-1)/(a*a);
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            Kinv[i][j] *= factor;
        }
    }
}



/**
 * Inverse of get_K(). Depending on the dominant direction, (L-1), (L-2) or (L-3) forms from Rohde
 * (2001) are used.
 */
void NavierStokes::get_Kinv(
    const dealii::Tensor<1,dim> &vel,
    const double a,
    const double H,
    const dealii::Tensor<1,dim>& dir,
    dealii::Tensor<2,dim+2> &Kinv
) const
{
    const double vel_n = dealii::scalar_product(vel, dir),
        ek = dealii::scalar_product(vel, vel),
        temp = 1/(a*a);
    
    // 1st row
    Kinv[0][0] = ((gma_-1)*ek + a*vel_n)*0.5*temp;
    for(int d=0; d<dim; d++) Kinv[0][1+d] = ((1-gma_)*vel[d] - a*dir[d])*0.5*temp;
    Kinv[0][4] = (gma_-1)*0.5*temp;

    // 2nd row
    Kinv[1][0] = 1-(gma_-1)*ek*temp;
    for(int d=0; d<dim; d++) Kinv[1][1+d] = (gma_-1)*vel[d]*temp;
    Kinv[1][4] = (1-gma_)*temp;

    // 3rd row
    Kinv[2][0] = ((gma_-1)*ek - a*vel_n)*0.5*temp;
    for(int d=0; d<dim; d++) Kinv[2][1+d] = ((1-gma_)*vel[d] + a*dir[d])*0.5*temp;
    Kinv[2][4] = (gma_-1)*0.5*temp;

    double dir_abs[3];
    for(int d=0; d<dim; d++) dir_abs[d] = fabs(dir[d]);

    if(dir_abs[0] > dir_abs[1] && dir_abs[0] > dir_abs[2]){
        // (L-1) form
        // 4th row
        Kinv[3][0] = (vel[1] - vel_n*dir[1])/dir[0];
        Kinv[3][1] = dir[1];
        Kinv[3][2] = (dir[1]*dir[1]-1)/dir[0];
        Kinv[3][3] = dir[1]*dir[2]/dir[0];
        Kinv[3][4] = 0;

        // 5th row
        Kinv[4][0] = (vel_n*dir[2] - vel[2])/dir[0];
        Kinv[4][1] = -dir[2];
        Kinv[4][2] = -dir[1]*dir[2]/dir[0];
        Kinv[4][3] = (1-dir[2]*dir[2])/dir[0];
        Kinv[4][4] = 0;
    }
    else if(dir_abs[1] > dir_abs[0] && dir_abs[1] > dir_abs[2]){
        // (L-2) form
        // 4th row
        Kinv[3][0] = (vel_n*dir[0]-vel[0])/dir[1];
        Kinv[3][1] = (1-dir[0]*dir[0])/dir[1];
        Kinv[3][2] = -dir[0];
        Kinv[3][3] = -dir[0]*dir[2]/dir[1];
        Kinv[3][4] = 0;

        // 5th row
        Kinv[4][0] = (vel[2] - vel_n*dir[2])/dir[1];
        Kinv[4][1] = dir[0]*dir[2]/dir[1];
        Kinv[4][2] = dir[2];
        Kinv[4][3] = (dir[2]*dir[2]-1)/dir[1];
        Kinv[4][4] = 0;
    }
    else{
        // (L-3) form
        // 4th row
        Kinv[3][0] = (vel[0]-vel_n*dir[0])/dir[2];
        Kinv[3][1] = (dir[0]*dir[0]-1)/dir[2];
        Kinv[3][2] = dir[0]*dir[1]/dir[2];
        Kinv[3][3] = dir[0];
        Kinv[3][4] = 0;

        // 5th row
        Kinv[4][0] = (vel_n*dir[1]-vel[1])/dir[2];
        Kinv[4][1] = -dir[0]*dir[1]/dir[2];
        Kinv[4][2] = (1-dir[1]*dir[1])/dir[2];
        Kinv[4][3] = -dir[1];
        Kinv[4][4] = 0;
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
 * @note See WJ-16-Feb-2022. A couple of other wave speed estimates are also implemented now. See
 * Toro, 3rd ed, sections 10.5 and 10.6.
 *
 * @pre @p lcs and @p rcs must pass the test of NavierStokes::assert_positivity()
 */
void NavierStokes::hllc_xflux(const State &lcs, const State &rcs, State &f) const
{
    dealii::Tensor<1,dim> xdir({1,0,0});

    // wave speed estimates
    // 1. Using Roe averaged speeds (direct estimation)
    // see Toro section 10.5.1
    /*
    const double sql = sqrt(lcs[0]), sqr = sqrt(rcs[0]); // 'sq'uare roots of densities
    const double pl = get_p(lcs), pr = get_p(rcs); // pressures
    const double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0];
    
    const double ut = (ul*sql + ur*sqr)/(sql + sqr); // u tilde
    const double Ht = ( (lcs[4]+pl)/sql + (rcs[4]+pr)/sqr )/(sql + sqr); // H tilde
    const double at = sqrt((gma_-1)*(Ht - 0.5*ut*ut)); // a tilde
    
    const double sl = ut-at, sr = ut+at; // left and right wave speeds
    */
    // 2. Using adaptive noninterative pressure estimate for star region
    // (pressure-based wave speed estimation)
    // see Toro section 10.6, section 9.5.2 and section 9.4
    // /*
    const double pl = get_p(lcs), pr = get_p(rcs); // pressures
    const double pmin = std::min(pl,pr), pmax = std::max(pl,pr);
    const double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0];
    const double al = std::sqrt(gma_*pl/lcs[0]), ar = std::sqrt(gma_*pr/rcs[0]);
    const double p_pvrs = 0.5*(pl+pr) - 0.5*(ur-ul)*0.5*(lcs[0]+rcs[0])*0.5*(al+ar);
    double ps; // p star
    if(pmax/pmin < 2 && pmin <= p_pvrs && p_pvrs <= pmax){
        // use pvrs estimate
        ps = p_pvrs;
    }
    else if(p_pvrs < pmin){
        // use trrs estimate (two rarefactions)
        const double z = (gma_-1)/(2*gma_);
        ps = std::pow(
            ((al+ar)-(gma_-1)/2*(ur-ul))/(al/std::pow(pl,z)+ar/std::pow(pr,z)),
            1/z
        );
    }
    else{
        // use tsrs estimate (two shocks)
        const double p0 = std::max(0.0, p_pvrs);
        const double gl = std::sqrt(2/((gma_+1)*lcs[0])/(p0+(gma_-1)/(gma_+1)*pl)),
            gr = std::sqrt(2/((gma_+1)*rcs[0])/(p0+(gma_-1)/(gma_+1)*pr));
        ps = (gl*pl + gr*pr - (ur-ul))/(gl+gr);
    }

    const double ql = ( ps <= pl ? 1 : std::sqrt(1+(gma_+1)/(2*gma_)*(ps/pl-1)) ),
        qr = ( ps <= pr ? 1 : std::sqrt(1+(gma_+1)/(2*gma_)*(ps/pr-1)) );
    const double sl = ul-al*ql, sr = ur+ar*qr;
    // */

    // star wave speed
    const double s = ( pr-pl + lcs[1]*(sl-ul) - rcs[1]*(sr-ur) )/(lcs[0]*(sl-ul) - rcs[0]*(sr-ur));
    
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
    const double pl = get_p(lcs), pr = get_p(rcs); // pressures
    const double ul = lcs[1]/lcs[0], ur = rcs[1]/rcs[0]; // velocities
    const double Hl = (lcs[4]+pl)/lcs[0], Hr = (rcs[4]+pr)/rcs[0]; // total/stagnation enthalpies

    const double astl = sqrt(2*(gma_-1)/(gma_+1)*Hl),
        astr = sqrt(2*(gma_-1)/(gma_+1)*Hr); // critical speed of sounds ('st'ar)
    const double ahl = astl*astl/std::max(astl, ul),
        ahr = astr*astr/std::max(astr, ur); // 'h'at velocities (eq 30 in Liou 2006)
    const double a12 = std::min(ahl, ahr); // a_1/2

    const double Ml = ul/a12, Mr = ur/a12; // Mach numbers
    const double Mavgsq = 0.5*(Ml*Ml + Mr*Mr);
    const double M0 = sqrt(std::min(1.0 ,std::max(Mavgsq, ausm::Minfty*ausm::Minfty)));
    const double fa = M0*(2-M0), alpha = 3/16*(-4+5*fa*fa);

    double M12 = ausm::mach_split_4_pos(Ml) + ausm::mach_split_4_neg(Mr) -
        ausm::Kp/fa*std::max(1-ausm::sigma*Mavgsq, 0.0)*2*(pr-pl)/(a12*a12*(lcs[0]+rcs[0]));

    double p_split_pos = ausm::pressure_split_5_pos(Ml, alpha),
        p_split_neg = ausm::pressure_split_5_neg(Mr, alpha);
    
    double p12 = p_split_pos*pl + p_split_neg*pr -
        ausm::Ku*fa*p_split_pos*p_split_neg*(lcs[0]+rcs[0])*a12*(ur-ul);
    
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
 *
 * @note The expression for corrected eigen values is wrong in the reference stated above. Instead,
 * refer to ref. 13 of Druguet, Candler & Nompelis (2005).
 */
void NavierStokes::modified_sw_xflux(const State &lcs, const State &rcs, State &f) const
{
    // Calculate omega
    const double pl = get_p(lcs), pr = get_p(rcs);
    const double sigma2 = 0.5;
    const double omega = 0.5/( 1 + std::pow(sigma2*(pr-pl)/std::min(pl,pr), 2) );

    // pos and neg states which are inputs for calculating pos and neg jacobians
    // these are also used to calculate pos and neg eigenvalues
    dealii::Tensor<1,dim> vel_pos, vel_neg;
    for(int d=0; d<dim; d++){
        vel_pos[d] = (1-omega)*lcs[1+d]/lcs[0] + omega*rcs[1+d]/rcs[0];
        vel_neg[d] = (1-omega)*rcs[1+d]/rcs[0] + omega*lcs[1+d]/lcs[0];
    }

    const double rho_pos = (1-omega)*lcs[0] + omega*rcs[0],
        rho_neg = (1-omega)*rcs[0] + omega*lcs[0];
    const double p_pos = (1-omega)*pl + omega*pr,
        p_neg = (1-omega)*pr + omega*pl;
    const double a_pos = std::sqrt(gma_*p_pos/rho_pos),
        a_neg = std::sqrt(gma_*p_neg/rho_neg);
    const double H_pos = gma_*p_pos/((gma_-1)*rho_pos) + 0.5*dealii::scalar_product(vel_pos, vel_pos),
        H_neg = gma_*p_neg/((gma_-1)*rho_neg) + 0.5*dealii::scalar_product(vel_neg, vel_neg);

    // pos and neg eigenvector matrices
    dealii::Tensor<2,dim+2> K_pos, Kinv_pos, K_neg, Kinv_neg;
    get_xK(vel_pos, a_pos, H_pos, K_pos);
    get_xK(vel_neg, a_neg, H_neg, K_neg);
    get_xKinv(vel_pos, a_pos, H_pos, Kinv_pos);
    get_xKinv(vel_neg, a_neg, H_neg, Kinv_neg);

    // pos and neg eigenvalues, and their correction
    const double eps = 0.3;
    std::array<double, dim+2> eig_pos, eig_neg;
    eig_pos[0] = pos_smooth(vel_pos[0] - a_pos, eps*a_pos);
    for(int d=0; d<dim; d++) eig_pos[1+d] = pos_smooth(vel_pos[0], eps*a_pos);
    eig_pos[4] = pos_smooth(vel_pos[0] + a_pos, eps*a_pos);
    eig_neg[0] = neg_smooth(vel_neg[0] - a_neg, eps*a_neg);
    for(int d=0; d<dim; d++) eig_neg[1+d] = neg_smooth(vel_neg[0], eps*a_neg);
    eig_neg[4] = neg_smooth(vel_neg[0] + a_neg, eps*a_neg);

    dealii::Tensor<2,dim+2> A_pos, A_neg;
    A_pos = 0;
    A_neg = 0;
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            for(int k=0; k<dim+2; k++){
                A_pos[i][j] += K_pos[i][k]*eig_pos[k]*Kinv_pos[k][j];
                A_neg[i][j] += K_neg[i][k]*eig_neg[k]*Kinv_neg[k][j];
            }
        }
    }

    // calculate the flux
    for(int i=0; i<dim+2; i++){
        f[i] = 0;
        for(int j=0; j<dim+2; j++){
            f[i] += A_pos[i][j]*lcs[j] + A_neg[i][j]*rcs[j];
        }
    }
}



/**
 * Chandrashekhar x-directional surface (numerical) flux. See Gassner et al (2016) section 3.4,
 * particularly eqs (3.20) and (3.27-29) for the expression. This is essentially Chandrashekhar's
 * volume flux added with a stabilisation term similar to the 2nd term in a Rusanov flux.
 *
 * Much of the implementation follows directly from chandrashekhar_vol_flux(), just that here the
 * direction is fixed to point in x.
 *
 * For matrix based stabilisation, see Chandrashekhar (2013) and Berberich & Klingenberg (2021).
 *
 * @note See WJ-17-Feb-2022. This was added much later compared to other surface flux schemes.
 */
void NavierStokes::chandrashekhar_xflux(
    const State &lcs, const State &rcs, State &f
) const
{
    const double pl = get_p(lcs), pr = get_p(rcs);
    const double betal = 0.5*lcs[0]/pl, betar = 0.5*rcs[0]/pr, beta_ln = log_avg(betal, betar);
    
    const double p_hat = 0.5*(lcs[0]+rcs[0])/(betal+betar);
    const double rho_ln = log_avg(lcs[0], rcs[0]);

    dealii::Tensor<1,dim> vel_avg, vel_sq_avg; // sq for 'sq'uare
    dealii::Tensor<1,dim> vl, vr; // left and right velocities
    const double rho_invl = 1/lcs[0], rho_invr = 1/rcs[0];
    for(int d=0; d<dim; d++){
        vl[d] = lcs[1+d]*rho_invl;
        vr[d] = rcs[1+d]*rho_invr;
        vel_avg[d] = 0.5*(vl[d]+vr[d]);
        vel_sq_avg[d] = 0.5*(vl[d]*vl[d] + vr[d]*vr[d]);
    }

    double H_hat = 0.5/((gma_-1)*beta_ln) + p_hat/rho_ln; // initialise
    for(int d=0; d<dim; d++){
        H_hat += vel_avg[d]*vel_avg[d] - 0.5*vel_sq_avg[d];
    }

    // volume flux (symmetric part)
    f[0] = rho_ln*vel_avg[0];
    f[1] = rho_ln*vel_avg[0]*vel_avg[0] + p_hat;
    f[2] = rho_ln*vel_avg[0]*vel_avg[1];
    f[3] = rho_ln*vel_avg[0]*vel_avg[2];
    f[4] = rho_ln*vel_avg[0]*H_hat;

    // Hybrid matrix based stabilisation
    // see sections 6 and 8 of Chandrashekhar (2013), and section 2 of Berberich & Klingenberg (2021)
    // for robustness, Roe averaging is used here
    const double rho_sql = sqrt(lcs[0]), rho_sqr = sqrt(rcs[0]);
    const dealii::Tensor<1,dim> vi = 1/(rho_sql+rho_sqr)*(rho_sql*vl + rho_sqr*vr);
    const double Hi = ( (lcs[4]+pl)/rho_sql + (rcs[4]+pr)/rho_sqr )/(rho_sql + rho_sqr);
    const double ai = sqrt((gma_-1)*(Hi - 0.5*dealii::scalar_product(vi, vi)));
    dealii::Tensor<2,dim+2> K, Kinv, A; // eigenvector matrix, diffusion matrix
    get_xK(vi, ai, Hi, K);
    get_xKinv(vi, ai, Hi, Kinv);
    // eigen values: blending KES and Rusanov eigenvalues with flux blender value
    const double u_abs = fabs(vi[0]),
        lambda_max_rus = u_abs + ai,
        lambda2_blend = flux_blender_value*lambda_max_rus + (1-flux_blender_value)*u_abs;
    std::array<double, dim+2> eig = {
        lambda_max_rus,
        lambda2_blend,
        lambda2_blend,
        lambda2_blend,
        lambda_max_rus
    };
    A = 0;
    // A = K * Lambda * Kinv
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            for(int k=0; k<dim+2; k++){
                A[i][j] += K[i][k]*eig[k]*Kinv[k][j];
            }
        }
    }
    // stabilisation
    for(int i=0; i<dim+2; i++){
        for(int j=0; j<dim+2; j++){
            f[i] -= 0.5*A[i][j]*(rcs[j] - lcs[j]);
        }
    }
}



/**
 * See eq (3.10) and (3.28) of Gassner et al (2016).
 */
void NavierStokes::kennedy_gruber_xflux(
    const State &lcs, const State &rcs, State &f
) const
{
    dealii::Tensor<1,dim> vel_avg;
    dealii::Tensor<1,dim> vl, vr;
    const double rho_invl = 1/lcs[0], rho_invr = 1/rcs[0];
    for(int d=0; d<dim; d++){
        vl[d] = lcs[1+d]*rho_invl;
        vr[d] = rcs[1+d]*rho_invr;
        vel_avg[d] = 0.5*(vl[d]+vr[d]);
    }
    const double rho_avg = 0.5*(lcs[0] + rcs[0]);
    const double el = get_e(lcs), er = get_e(rcs);
    const double E_avg = 0.5*(
        el + er +
        0.5*(dealii::scalar_product(vl,vl) + dealii::scalar_product(vr,vr))
    );
    const double p_avg = (gma_-1)*0.5*(lcs[0]*el + rcs[0]*er);
    f[0] = vel_avg[0]*rho_avg;
    for(int d=0; d<dim; d++){
        f[1+d] = rho_avg*vel_avg[0]*vel_avg[d];
    }
    f[1] += p_avg;
    f[4] = vel_avg[0]*(rho_avg*E_avg + p_avg);

    // stabilisation term
    const double lambda_max = std::max(
        fabs(vl[0]) + std::sqrt(gma_*(gma_-1)*el),
        fabs(vr[0]) + std::sqrt(gma_*(gma_-1)*er)
    );
    for(cvar var: cvar_list) f[var] -= 0.5*lambda_max*(rcs[var] - lcs[var]);
}



/**
 * Blend of Rusanov and Kennedy-Gruber surface fluxes.
 */
void NavierStokes::rusanov_kennedy_gruber_blend_xflux(
    const State &lcs, const State &rcs, State &f
) const
{
    State f_rus, f_kg;
    rusanov_xflux(lcs, rcs, f_rus);
    kennedy_gruber_xflux(lcs, rcs, f_kg);
    for(cvar var: cvar_list) f[var] = flux_blender_value*f_rus[var] +
        (1-flux_blender_value)*f_kg[var];
}



/**
 * Ismail & Roe's surface flux. See ismail_roe_vol_flux() and eqs. (3.30-31) of Gassner et al
 * (2016).
 */
void NavierStokes::ismail_roe_xflux(
    const State &lcs, const State &rcs, State &f
) const
{
    const double gma_sqrt = sqrt(gma_); // square root of gamma
    const double al = get_a(lcs), ar = get_a(rcs);
    const double z1l = gma_sqrt/al, z1r = gma_sqrt/ar; // z1 of states 1 and 2
    const double z5l = lcs[0]*al/gma_sqrt, z5r = rcs[0]*ar/gma_sqrt; // z5 of the states

    const double z5_ln = log_avg(z5l, z5r);
    const double rho_hat = 0.5*(z1l + z1r)*z5_ln;
    const double p1_hat = (z5l+z5r)/(z1l+z1r),
        p2_hat = 0.5*( (gma_+1)*z5_ln/log_avg(z1l, z1r) + (gma_-1)*p1_hat )/(gma_);
    double h_hat = gma_*p2_hat/((gma_-1)*rho_hat); // initialised here, will be modified using velocities
    dealii::Tensor<1,dim> vl, vr, vel_hat;
    for(int d=0; d<dim; d++){
        vl[d] = lcs[1+d]/lcs[0], vr[d] = rcs[1+d]/rcs[0];
        vel_hat[d] = (z1l*vl[d] + z1r*vr[d])/(z1l + z1r);
        h_hat += 0.5*(vel_hat[d]*vel_hat[d]);
    }
    f[0] = rho_hat*vel_hat[0];
    for(int d=0; d<dim; d++) f[1+d] = rho_hat*vel_hat[0]*vel_hat[d];
    f[1] += p1_hat;
    f[4] = rho_hat*h_hat*vel_hat[0];

    // stabilisation term
    State cons_hat = {
        rho_hat,
        rho_hat*vel_hat[0],
        rho_hat*vel_hat[1],
        rho_hat*vel_hat[2],
        // p2_hat/(gma_-1) + rho_hat*dealii::scalar_product(vel_hat, vel_hat)
        rho_hat*h_hat - p2_hat
    };
    const double pl = get_p(lcs), pr = get_p(rcs);
    const double sl = log(pl) - gma_*log(lcs[0]), sr = log(pr) - gma_*log(rcs[0]);
    dealii::Tensor<1,dim+2> Vl, Vr; // entropy variables
    Vl[0] = (gma_-sl)/(gma_-1) - 0.5*lcs[0]*dealii::scalar_product(vl,vl)/pl;
    Vr[0] = (gma_-sr)/(gma_-1) - 0.5*rcs[0]*dealii::scalar_product(vr,vr)/pr;
    Vl[4] = -lcs[0]/pl;
    Vr[4] = -rcs[0]/pr;
    for(int d=0; d<dim; d++){
        Vl[1+d] = lcs[1+d]/pl;
        Vr[1+d] = rcs[1+d]/pr;
    }
    dealii::Tensor<2,dim+2> H; // using p2_hat for p
    // 1st col
    for(int i=0; i<dim+2; i++) H[i][0] = cons_hat[i];
    // 2-4 cols
    for(int d=0; d<dim; d++){
        State flux;
        dealii::Tensor<1,dim> dir; // initialised to 0
        dir[d] = 1.0;
        get_inv_flux(cons_hat, dir, flux);
        for(int i=0; i<dim+2; i++) H[i][1+d] = flux[i];
    }
    // last col
    for(int i=0; i<4; i++) H[i][4] = H[4][i];
    H[4][4] = rho_hat*h_hat*h_hat - gma_*p2_hat*p2_hat/(rho_hat*(gma_-1));

    const double lambda_max = fabs(vel_hat[0]) + sqrt(gma_*p2_hat/rho_hat);
    dealii::Tensor<1,dim+2> f_stab = 0.5*lambda_max*(H*(Vr-Vl));

    for(int i=0; i<dim+2; i++) f[i] -= f_stab[i];
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
void NavierStokes::chandrashekhar_vol_flux(
    const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
) const
{
    const double p1 = get_p(cs1), p2 = get_p(cs2);
    const double beta1 = 0.5*cs1[0]/p1, beta2 = 0.5*cs2[0]/p2, beta_ln = log_avg(beta1, beta2);
    
    const double p_hat = 0.5*(cs1[0]+cs2[0])/(beta1+beta2);
    const double rho_ln = log_avg(cs1[0], cs2[0]);
    
    dealii::Tensor<1,dim> vel_avg, vel_sq_avg; // sq for 'sq'uare
    double v1, v2; // temporary quantities
    const double rho_inv1 = 1/cs1[0], rho_inv2 = 1/cs2[0];
    for(int d=0; d<dim; d++){
        v1 = cs1[1+d]*rho_inv1;
        v2 = cs2[1+d]*rho_inv2;
        vel_avg[d] = 0.5*(v1+v2);
        vel_sq_avg[d] = 0.5*(v1*v1 + v2*v2);
    }
    
    double H_hat = 0.5/((gma_-1)*beta_ln) + p_hat/rho_ln; // initialise
    for(int d=0; d<dim; d++){
        H_hat += vel_avg[d]*vel_avg[d] - 0.5*vel_sq_avg[d];
    }
    
    const double vel_n = dealii::scalar_product(vel_avg, dir); // velocity in the direction 'dir'
    
    f[0] = vel_n*rho_ln;
    for(int d=0; d<dim; d++){
        f[1+d] = rho_ln*vel_n*vel_avg[d] + p_hat*dir[d];
    }
    f[4] = rho_ln*vel_n*H_hat;
}



/**
 * See eq (3.10) of Gassner et al (2016).
 */
void NavierStokes::kennedy_gruber_vol_flux(
    const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
) const
{
    dealii::Tensor<1,dim> vel_avg;
    dealii::Tensor<1,dim> vel1, vel2;
    const double rho_inv1 = 1/cs1[0], rho_inv2 = 1/cs2[0];
    for(int d=0; d<dim; d++){
        vel1[d] = cs1[1+d]*rho_inv1;
        vel2[d] = cs2[1+d]*rho_inv2;
        vel_avg[d] = 0.5*(vel1[d]+vel2[d]);
    }
    const double rho_avg = 0.5*(cs1[0] + cs2[0]);
    const double e1 = get_e(cs1), e2 = get_e(cs2);
    const double E_avg = 0.5*(
        e1 + e2 +
        0.5*(dealii::scalar_product(vel1, vel1) + dealii::scalar_product(vel2,vel2))
    );
    const double p_avg = (gma_-1)*0.5*(cs1[0]*e1 + cs2[0]*e2);
    const double vel_n = dealii::scalar_product(vel_avg, dir); // velocity in the direction 'dir'
    f[0] = vel_n*rho_avg;
    for(int d=0; d<dim; d++){
        f[1+d] = rho_avg*vel_n*vel_avg[d] + p_avg*dir[d];
    }
    f[4] = vel_n*(rho_avg*E_avg + p_avg);
}



/**
 * Ismail & Roe's volume flux. See eqs. (3.14) - (3.17) of Gassner et al (2016).
 *
 * To spare some calculations, the following relations are used
 * @f[
 * z_1 = \sqrt{\rho/p} = \sqrt{\gamma}/a \\
 * z_5 = \sqrt{p\rho} = \rho a/\sqrt{\gamma}
 * @f]
 *
 * @pre @p dir has to be a unit vector
 * @pre @p cs1 and @p cs2 must pass the test of NavierStokes::assert_positivity()
 */
void NavierStokes::ismail_roe_vol_flux(
    const State &cs1, const State &cs2, const dealii::Tensor<1,dim> &dir, State &f
) const
{
    const double gma_sqrt = sqrt(gma_); // square root of gamma
    const double a1 = get_a(cs1), a2 = get_a(cs2);
    const double z11 = gma_sqrt/a1, z12 = gma_sqrt/a2; // z1 of states 1 and 2
    const double z51 = cs1[0]*a1/gma_sqrt, z52 = cs2[0]*a2/gma_sqrt; // z5 of the states

    const double z5_ln = log_avg(z51, z52);
    const double rho_hat = 0.5*(z11 + z12)*z5_ln;
    const double p1_hat = (z51+z52)/(z11+z12),
        p2_hat = 0.5*( (gma_+1)*z5_ln/log_avg(z11, z12) + (gma_-1)*p1_hat )/(gma_);
    double h_hat = gma_*p2_hat/((gma_-1)*rho_hat); // initialised here, will be modified using velocities
    dealii::Tensor<1,dim> vel_hat;
    for(int d=0; d<dim; d++){
        const double v1 = cs1[1+d]/cs1[0], v2 = cs2[1+d]/cs2[0];
        vel_hat[d] = (z11*v1 + z12*v2)/(z11 + z12);
        h_hat += 0.5*(vel_hat[d]*vel_hat[d]);
    }

    const double vel_n = dealii::scalar_product(vel_hat, dir);
    f[0] = rho_hat*vel_n;
    for(int d=0; d<dim; d++){
        f[1+d] = rho_hat*vel_hat[d]*vel_n + p1_hat*dir[d];
    }
    f[4] = rho_hat*vel_n*h_hat;
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

