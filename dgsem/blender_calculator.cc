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
 * Calls BlenderCalculator::get_trouble() to calculate the trouble and then returns the modified
 * blender value.
 */
double BlenderCalculator::get_blender(
    const DoFHandler<dim>::active_cell_iterator& cell
) const
{
    const usi n_dofs_per_cell = cell->get_fe().dofs_per_cell;
    std::vector<psize> dof_ids(n_dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    // get the modal coefficients (or modes)
    std::vector<double> modes(n_dofs_per_cell);
    for(usi row=0; row<n_dofs_per_cell; row++){
        modes[row] = 0;
        for(usi col=0; col<n_dofs_per_cell; col++){
            modes[row] += cbm(row,col)*variable[dof_ids[col]];
        }
    } // loop over modes

    // calculate alpha
    const double trouble = get_trouble(modes);
    const double threshold = threshold_factor*std::pow(
        10,
        -threshold_exp_factor*(
            std::pow(cbm.degree+1,0.25)
        )
    );
    const double alpha_tilde = 1/( 1 + exp(-sharpness_factor*(trouble/threshold-1)) );

    if(alpha_tilde < alpha_min) return alpha_min;
    else if(alpha_tilde < alpha_max) return alpha_tilde;
    else return alpha_max;
}



// * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * * //



/**
 * Calculates the troubled based on the modal coefficients provided in `modes`. Currently, there
 * is ambiguity in the choice of algorithm for this function. The approach presented in Perrson &
 * Peraire (2006), I think, will not work for purely one-dimensional variations. An example for
 * this is given in WJ-25-May-2021.
 *
 * However, until a better algorithm is found, this will be implemented as is. The actual algorithm
 * implemented also uses the small modification proposed by Hennemann et al (2021) in section 4.
 * This function will need modification after some testing.
 *
 * @pre `modes` must have a size `(degree+1)^dim`
 */
 double BlenderCalculator::get_trouble(
    const std::vector<double>& modes
) const
{
    AssertThrow(
        modes.size() == std::pow(cbm.degree+1, dim),
        StandardExceptions::ExcMessage(
            "The size of modes provided doesn't match with the expected size: (degree+1)^dim."
        )
    );

    double e_tot = 0; // total energy of all modes
    for(double m: modes) e_tot += m*m;
    double e_Nm1 = 0; // energy of modes only upto (N-1)-th polynomial degree
    for(usi i=0; i<std::pow(cbm.degree,dim); i++) e_Nm1 += modes[i]*modes[i];
    double e_Nm2 = 0; // energy of modes only upto (N-2)-th polynomial degree
    for(usi i=0; i<std::pow(cbm.degree-1,dim); i++) e_Nm2 += modes[i]*modes[i];

    if(cbm.degree>1) return std::max(1-e_Nm1/e_tot, 1-e_Nm2/e_Nm1);
    else return 1-e_Nm1/e_tot;
}
