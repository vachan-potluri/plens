/**
 * @file varying_inflow.cc
 * @brief Spatially and temporally varying inflow boundary condition
 */

#include "varying_inflow.h"

using namespace BCs;

/**
 * Calculates the prescribed state at the location given through @p ldd and the current time (the
 * value in BCs::BC::time_). Uses BCs::BC::get_global_dof_id() internally along with
 * VaryingInflow::dof_locations.
 *
 * @pre Assumes that BCs::BC::time_ variable is set appropriately.
 */
void VaryingInflow::calculate_prescribed_state(const FaceLocalDoFData &ldd, State &s)
{
    // set the time
    p_pr_fn_.set_time(time_);
    T_pr_fn_.set_time(time_);
    velocity_pr_fn_.set_time(time_);

    // get dof location
    const Point<dim> loc = dof_locations.at(get_global_dof_id(ldd));

    // get primitive quantities
    const double p = p_pr_fn_.value(loc);
    const double T = T_pr_fn_.value(loc);
    Tensor<1,dim> velocity;
    for(int d=0; d<dim; d++) velocity[d] = velocity_pr_fn_.value(loc, d);

    // get conservative state
    ns_ptr_->prim_to_cons(p/(ns_ptr_->get_R()*T), velocity, p, s);
}



/**
 * Constructor. Calls the base constructor and sets the function parsers.
 *
 * @pre The velocity expression @p velocity_expr is assumed to have all the components, separated
 * by ";". This way, dealii's version of `FunctionParser::initialise()` will directly work even for
 * velocity.
 */
VaryingInflow::VaryingInflow(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav,
    const std::map<psize, Point<dim>> &dl,
    const std::string& p_expr,
    const std::string& T_expr,
    const std::string& velocity_expr,
    const NavierStokes* ns_ptr
):
BC("varying inflow", dh, gcv, gav),
dof_locations(dl),
ns_ptr_(ns_ptr),
p_pr_fn_(1),
T_pr_fn_(1),
velocity_pr_fn_(dim)
{
    // variables and constants
    std::string variables("x,y,z,t");
    std::map<std::string, double> constants; // empty

    // set the functions
    p_pr_fn_.initialize(
        variables, p_expr, constants, true
    );
    T_pr_fn_.initialize(
        variables, T_expr, constants, true
    );
    velocity_pr_fn_.initialize(
        variables, velocity_expr, constants, true
    ); // works if velocity_expr has components separated by ";"
}



/**
 * Exactly similar to BCs::UniformInflow::get_ghost_stage1(), except that the functions for
 * pressure temperature and velocity are evaluated to get the ghost state. The evaluation is left
 * for VaryingInflow::calculate_prescribed_state() to do.
 */
void VaryingInflow::get_ghost_stage1(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    State cons_pr;
    calculate_prescribed_state(ldd, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
}



/**
 * Sets the ghost state equal to prescribed state, like in BCs::UniformInflow::get_ghost_stage2().
 * Uses VaryingInflow::calculate_prescribed_state() internally to calculate the ghost state.
 *
 * @note @p normal and @p ldd are unused
 */
void VaryingInflow::get_ghost_stage2(
    const FaceLocalDoFData &ldd,
    const State &cons,
    const Tensor<1,dim> &normal,
    State &cons_gh
)
{
    State cons_pr;
    calculate_prescribed_state(ldd, cons_pr);
    cons_gh = cons_pr;
}



/**
 * Ghost getter for stage 3. Again, the implementation is exactly similar to that in
 * BCs::UniformInflow::get_ghost_stage3(). Uses
 * VaryingInflow::calculate_prescribed_state() internally to calculate the ghost state.
 */
void VaryingInflow::get_ghost_stage3(
    const FaceLocalDoFData &ldd,
    const CAvars &cav,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
)
{
    const State& cons = cav.get_state();
    State& cons_gh = cav_gh.get_state();
    State cons_pr;
    calculate_prescribed_state(ldd, cons_pr);
    for(cvar var: cvar_list) cons_gh[var] = 2*cons_pr[var] - cons[var];
    
    const Avars& av = cav.get_avars();
    Avars& av_gh = cav_gh.get_avars();
    av_gh = av;
}
