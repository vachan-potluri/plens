/**
 * @file blender_calculator.h
 * @brief Calculates the value of blending parameter, cell-wise.
 */

#include "blender_calculator.h"



/**
 * Constructor. Parses all blender parameters. Calls `enter_subsection("blender parameters")` and
 * then exits that subsection.
 */
BlenderCalculator::BlenderCalculator(
    const usi d,
    const LA::MPI::Vector& var,
    ParameterHandler& prm
):
variable(var),
cbm(d)
{
    prm.enter_subsection("blender parameters");
    {
        threshold_factor = prm.get_double("threshold factor");
        threshold_exp_factor = prm.get_double("threshold exponent factor");
        sharpness_factor = prm.get_double("sharpness factor");
        alpha_min = prm.get_double("blender min value");
        alpha_max = prm.get_double("blender max value");
    }
    prm.leave_subsection();
}



/**
 * Calculates blender value for the given cell according to section 4 of Hennemann et al (2021).
 * Calls BlenderCalculator::get_troubled() to calculate the trouble and then returns the modified
 * blender value.
 */
double BlenderCalculator::get_blender(
    const DoFHandler<dim>::active_cell_iterator& cell
) const
{}



// * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * * //



/**
 * Calculates the troubled based on the modal coefficients provided in `modes`.
 *
 * @pre `modes` must have a size `(degree+1)^dim`
 */
 double BlenderCalculator::get_trouble(
    const std::vector<double>& modes
) const
{
    return 0;
}
