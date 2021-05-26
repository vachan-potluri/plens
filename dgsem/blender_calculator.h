/**
 * @file blender_calculator.h
 * @brief Calculates the value of blending parameter, cell-wise.
 */

#ifndef BLENDER_CALCULATOR_H
#define BLENDER_CALCULATOR_H

#include <deal.II/base/parameter_handler.h>

#include "LA.h"

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
 * The conversion from nodal to modal basis happens using the class ChangeOfBasisMatrix. One of the
 * most important aspect is the calculation of trouble (@f$\mathbb{E}@f$). Hennemann et al (2021)
 * provide the formula only for 1D.
 */
class BlenderCalculator
{
    public:

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

    BlenderCalculator(
        const ParameterHandler& prm,
        const LA::MPI::Vector& var
    );
};

#endif
