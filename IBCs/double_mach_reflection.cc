/**
 * @file double_mach_reflection.cc
 * @brief Sets initial condition for the double Mach reflection (DMR) case.
 */

#include "double_mach_reflection.h"

using namespace ICs;

/**
 * Constructor. Calls the base constructor and sets the shock position function.
 *
 * @note z-components of @p wedge_le_loc and @p shock_offset are ignored. @p wedge_angle is assumed
 * to be given in degrees.
 */
DoubleMachReflection::DoubleMachReflection(
    const DoFHandler<dim> &dh,
    const std::map<psize, Point<dim>> &dl,
    std::array<LA::MPI::Vector, 5> &gcv,
    const Point<dim>& wedge_le_loc,
    const double wedge_angle,
    const Tensor<1,dim>& shock_offset
):
IC(dh, dl, gcv),
shock_pos_function_(1)
{
    std::string variables("x,y");
    std::map<std::string, double> constants;
    constants["pi"] = numbers::PI;

    Point<dim> shock_point = wedge_le_loc + shock_offset;
    std::stringstream shock_pos_expr;
    shock_pos_expr << "(x - " << shock_point[0] << ") - "
        << "(y - " << shock_point[1] << ") * "
        << "cot((60+" << wedge_angle << ")*pi/180)";
    shock_pos_function_.initialize(variables, shock_pos_expr.str(), constants, false);
}



/**
 * Sets the IC. Algo is simple:
 * - For each owned dof
 *   - If DoubleMachReflection::shock_pos_function_ evaluates to positive value
 *     - Set pre shock state
 *   - Else
 *     - Set post-shock state
 *
 * The values used here are from Hennemenn et al (2021).
 */
void DoubleMachReflection::set()
{
    std::vector<psize> dof_ids(dof_handler.get_fe().dofs_per_cell); // dof ids of cell
    for(const auto &cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        cell->get_dof_indices(dof_ids);

        for(psize i: dof_ids){
            if(shock_pos_function_.value(dof_locations[i]) > 0){
                // pre shock
                g_cvars[0][i] = 1.4;
                g_cvars[1][i] = 0;
                g_cvars[2][i] = 0;
                g_cvars[3][i] = 0;
                g_cvars[4][i] = 1/0.4;
            }
            else{
                // post shock
                g_cvars[0][i] = 8;
                g_cvars[1][i] = 8*7.144709581221619;
                g_cvars[2][i] = -8*4.125;
                g_cvars[3][i] = 0;
                g_cvars[4][i] = 116.5/0.4 + 0.5*8*(
                    7.144709581221619*7.144709581221619 + 4.125*4.125
                );
            }
        }
    }
}
