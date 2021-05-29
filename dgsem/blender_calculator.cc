/**
 * @file blender_calculator.h
 * @brief Calculates the value of blending parameter, cell-wise.
 */

#include "blender_calculator.h"



/**
 * Constructor. Populates BlenderCalculator::mode_indices_Nm1 and
 * BlenderCalculator::mode_indices_Nm2. These indices are set such that they give tensor product
 * polynomials with uni-directional polynomials of degree at most @f$N-1@f$ and @f$N-2@f$
 * respectively.
 *
 * @note BlenderCalculator::mode_indices_Nm2 is empty if `d=1`. And for `d=2`, even though this
 * array is populated, it is not used. See BlenderCalculator::get_trouble()
 *
 * @warning The construction of this object will be incomplete unless
 * BlenderCalculator::parse_parameters() is called. This is where the parameters are actually
 * read. If BlenderCalculator::get_blender() is called without a call to
 * BlenderCalculator::parse_parameters(), then an exception will be thrown.
 */
BlenderCalculator::BlenderCalculator(
    const usi d,
    const LA::MPI::Vector& var
):
parsed_params(false),
variable(var),
cbm(d)
{
    // for mode_indices_Nm1
    for(usi k=0; k<d; k++){
        for(usi j=0; j<d; j++){
            for(usi i=0; i<d; i++){
                usi index = i + j*(d+1) + k*(d+1)*(d+1);
                mode_indices_Nm1.emplace_back(index);
            }
        }
    }

    // for mode_indices_Nm2
    // although calculated only when d>1, these are actually used only if d>2
    if(d>1){
        for(usi k=0; k<d-1; k++){
            for(usi j=0; j<d-1; j++){
                for(usi i=0; i<d; i++){
                    usi index = i + j*(d+1) + k*(d+1)*(d+1);
                    mode_indices_Nm2.emplace_back(index);
                }
            }
        }
    }
}



/**
 * Parses the "blender parameters" section from the `prm` object provided. The calls made are
 * one `prm.enter_subsection("blender parameters")` and `prm.leave_subsection()` together with
 * the getters.
 */
void BlenderCalculator::parse_parameters(ParameterHandler& prm)
{
    parsed_params = true;

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
 *
 * @note For this function to be available, a call to BlenderCalculator::parse_parameters() is
 * necessary. Otherwise, the relevant parameters would have garbage value. An exception is thrown
 * is such a call is not made.
 */
double BlenderCalculator::get_blender(
    const DoFHandler<dim>::active_cell_iterator& cell
) const
{
    AssertThrow(
        parsed_params,
        StandardExceptions::ExcMessage(
            "You are trying to call get_blender() from a BlenderCalculator object that is not yet "
            "completely set. You must first parse the parameters through "
            "BlenderCalculator::parse_parameters() before calling this function."
        )
    );

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
 * Calculates the troubled based on the modal coefficients provided in `modes`. See the detailed
 * documentation for the formula used. Basically, it is a modification of Persson & Peraire
 * (2006)'s formula.
 *
 * @note The modification by Hennemann et al (2021) will be considered only if @f$N>2@f$. Otherwise
 * it doesn't make sense.
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
    for(usi i: mode_indices_Nm1) e_Nm1 += modes[i]*modes[i];
    double e_Nm2 = 0; // energy of modes only upto (N-2)-th polynomial degree
    for(usi i: mode_indices_Nm2) e_Nm2 += modes[i]*modes[i];

    if(cbm.degree>2) return std::max(1-e_Nm1/e_tot, 1-e_Nm2/e_Nm1);
    else return 1-e_Nm1/e_tot; // degree=1
}



#ifdef DEBUG
void BlenderCalculator::test()
{
    utilities::Testing t("BlenderCalculator", "class");

    ParameterHandler prm;
    prm.enter_subsection("blender parameters");
    {
        prm.declare_entry(
            "variable",
            "pxrho",
            Patterns::Selection("pxrho|p|rho"),
            "The variable used to calculate trouble (and hence the blender). Options: "
            "'pxrho|p|rho'. 'pxrho' means the product of pressure and density."
        );
        prm.declare_entry(
            "threshold factor",
            "0.5",
            Patterns::Double(1e-2),
            "The pre-factor for threshold calculation. See WJ-24-May-2021."
        );
        prm.declare_entry(
            "threshold exponent factor",
            "1.8",
            Patterns::Double(1e-2),
            "The pre-factor of exponent for threshold calculation. See WJ-24-May-2021. Range: "
            "[1e-2, infty)."
        );
        prm.declare_entry(
            "sharpness factor",
            "9.21024",
            Patterns::Double(1e-2),
            "The sharpness factor for calculating blender value from trouble and threshold. See "
            "WJ-24-May-2021. Range: [1e-2, infty)."
        );
        prm.declare_entry(
            "blender min value",
            "1e-3",
            Patterns::Double(1e-8,1),
            "Minimum value of blender (alpha min). See WJ-24-May-2021. Range: [1e-8, 1]."
        );
        prm.declare_entry(
            "blender max value",
            "0.5",
            Patterns::Double(0.5,1),
            "Maximum value of blender (alpha max). See WJ-24-May-2021. Range: [0.5, 1]."
        );
    }
    prm.leave_subsection(); // blender parameters
    prm.parse_input("input.prm", "", true); // ignores undefined entries in input.prm

    // parallel data structures being used just to be compatible
    // can comfortably test this function in serial as well
    parallel::distributed::Triangulation<dim> triang(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(triang); // generates [0,1]^dim with 1 cell
    DoFHandler<dim> dof_handler(triang);
    const usi fe_degree = 3;
    std::cout << "FE degree: " << fe_degree << "\n";
    FE_DGQ<dim> fe(fe_degree);
    dof_handler.distribute_dofs(fe);

    IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    LA::MPI::Vector var;
    var.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    std::map<unsigned int, Point<dim>> dof_locations;
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, dof_locations);

    t.new_block("testing construction");
    {
        BlenderCalculator bc(fe_degree, var);
        bc.parse_parameters(prm);
        std::cout << "Parameters read:\n" << bc.threshold_factor << "\n"
            << bc.threshold_exp_factor << "\n" << bc.sharpness_factor << "\n"
            << bc.alpha_min << "\n" << bc.alpha_max << "\n";
    }

    t.new_block("testing get_blender()");
    {
        // assuming dof handler has only 1 cell spanning [0,1]^3
        const auto& cell = dof_handler.begin_active();
        for(auto i: locally_owned_dofs){
            // step from 0 to 1 in x-direction
            var[i] = 0;
            if(dof_locations[i][0]>=0.5) var[i] = 1;
        }

        BlenderCalculator bc(fe_degree, var);
        bc.parse_parameters(prm);
        std::cout << "step from 0 to 1 in x-direction\n";
        std::cout << "Blender value: " << bc.get_blender(cell) << "\n";

        for(auto i: locally_owned_dofs){
            // linear variation in x, y and z
            var[i] = dof_locations[i][0] + 2*dof_locations[i][1] + 5*dof_locations[i][2];
        }
        std::cout << "linear variation in x, y and z\n";
        std::cout << "Blender value: " << bc.get_blender(cell) << "\n";
    }

    t.new_block("testing exception handling");
    {
        const auto& cell = dof_handler.begin_active();
        BlenderCalculator bc(fe_degree, var);
        bc.get_blender(cell); // exception: parse_parameters not called
    }
}
#endif
