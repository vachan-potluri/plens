/**
 * @file blender_calculator.h
 * @brief Calculates the value of blending parameter, cell-wise.
 */

#ifndef BLENDER_CALCULATOR_H
#define BLENDER_CALCULATOR_H

#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_handler.h>

#include "LA.h"
#include "change_of_basis_matrix.h"

using namespace dealii;

/**
 * Calculates the value of blending parameter, cell-wise. The inputs to this class are from the
 * 'blender parameters' section of the ParameterHandler provided for construction. The parameters
 * are stored only here, and not in PLENS because this is the only place where they are required.
 * However, the entry 'blender parameters/variable' is ignored and instead, a const reference to
 * the variable vector will also be taken in the constructor. All variables read from the
 * ParameterHandler are kept public in this class. Hence it is suggested to use a const version of
 * this class object.
 *
 * @note It will be assumed that the ParameterHandler is in the outer most scope from where, the
 * relevant section is accessible and is one level down.
 *
 * The conversion from nodal to modal basis happens using the class ChangeOfBasisMatrix. One of the
 * most important aspect is the calculation of trouble (@f$\mathbb{E}@f$). Hennemann et al (2021)
 * provide the formula only for 1D.
 */
class BlenderCalculator
{
    public:

    static constexpr usi dim = 3;

    /**
     * Const reference to the variable vector passed in the ctor
     */
    const LA::MPI::Vector& variable;

    /**
     * The threshold (pre-)factor
     */
    double threshold_factor;

    /**
     * The threshold exponent's factor
     */
    double threshold_exp_factor;

    /**
     * The sharpness factor
     */
    double sharpness_factor;

    /**
     * Min value of alpha (blender)
     */
    double alpha_min;

    /**
     * Max value of alpha (blender)
     */
    double alpha_max;

    /**
     * The change of basis matrix. The also stores the degree, so a separate variable for that is
     * not required.
     */
    const ChangeOfBasisMatrix<dim> cbm;

    BlenderCalculator(
        const usi d,
        const LA::MPI::Vector& var,
        ParameterHandler& prm
    );

    double get_blender(
        const DoFHandler<dim>::active_cell_iterator& cell
    ) const;
};

#endif
