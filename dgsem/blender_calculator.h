/**
 * @file blender_calculator.h
 * @brief Calculates the value of blending parameter, cell-wise.
 */

#ifndef BLENDER_CALCULATOR_H
#define BLENDER_CALCULATOR_H

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include "LA.h"
#include "change_of_basis_matrix.h"
#include "dtype_aliases.h"

#include <vector>
#include <cmath>
#include <algorithm>

#ifdef DEBUG
#include <deal.II/base/index_set.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <iostream>
#include <utilities/testing.h>
#endif

using namespace dealii;

/**
 * Calculates the value of blending parameter, cell-wise. The construction of this object is
 * complete only after BlenderCalculator::parse_parameters() is called. This function reads the
 * parameters under "blender parameters" section of the ParameterHandler provided. The parameters
 * are stored only here, and not in PLENS because this is the only place where they are required.
 * However, the entry 'blender parameters/variable' is ignored and instead, a const reference to
 * the variable vector will also be taken in the constructor. All variables read from the
 * ParameterHandler are kept private in this class, since these parameters are not generally
 * required anywhere else. If any such requirement comes up in the future, some getters will be
 * added.
 *
 * @note It will be assumed that the ParameterHandler is in the outer most scope from where, the
 * relevant section (viz. "blender parameters") is accessible and is one level down.
 *
 * The conversion from nodal to modal basis happens using the class ChangeOfBasisMatrix. One of the
 * most important aspect is the calculation of trouble (@f$\mathbb{E}@f$). Hennemann et al (2021)
 * provide the formula only for 1D. For multidimensional problems, the formula of Perrson &
 * Peraire (2006) is used with a modification. See WJ-26-May-2021, WJ-27-May-2021. Let @f$m@f$
 * represent the modal coefficients (or modes), then
 * @f[
 * \mathbb{E} = \max \left\{
 * 1-\frac{\sum_{i,j,k=0}^{N-1} m_{i,j,k}}{\sum_{i,j,k=0}^{N} m_{i,j,k}},
 * 1-\frac{\sum_{i,j,k=0}^{N-2} m_{i,j,k}}{\sum_{i,j,k=0}^{N-1} m_{i,j,k}},
 * \right\}
 * @f]
 * where the second term is the modification by Hennemann et al (2021). Note that the second term
 * makes sense only when @f$N>2@f$. Otherwise, it would give counter-intuitive results.
 */
class BlenderCalculator
{
    private:

    /**
     * A boolean variable depicting the state of the object. This will be set to true only when
     * BlenderCalculator::parse_parameters() has been called. BlenderCalculator::get_blender() can
     * be called only when this variable is `true`.
     */
    bool parsed_params;

    /**
     * The indices in a modes list that indicate modes of order upto @f$N-1@f$. This is required
     * because the modes list has to be in 3d index ordering (as against tensorial ordering) to be
     * usable with ChangeOfBasisMatrix. This vector is populated in the ctor.
     */
    std::vector<usi> mode_indices_Nm1;

    /**
     * Similar to BlenderCalculator::mode_indices_Nm1 but for modes of order upto @f$N-2@f$. This
     * is also populated in the ctor and is empty if the degree passed in the ctor is 1.
     */
    std::vector<usi> mode_indices_Nm2;

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

    double get_trouble(const std::vector<double>&) const;

    public:

    static constexpr usi dim = 3;

    /**
     * Const reference to the variable vector passed in the ctor
     */
    const LA::MPI::Vector& variable;

    /**
     * The change of basis matrix. The also stores the degree, so a separate variable for that is
     * not required.
     */
    const ChangeOfBasisMatrix<dim> cbm;

    BlenderCalculator(
        const usi d,
        const LA::MPI::Vector& var
    );

    void parse_parameters(ParameterHandler& prm);

    double get_blender(
        const DoFHandler<dim>::active_cell_iterator& cell
    ) const;

    /**
     * Returns the value of @f$\alpha_\text{max}@f$ (stored in BlenderCalculator::alpha_max).
     */
    inline double get_blender_max_value() const {return alpha_max;}

    #ifdef DEBUG
    static void test();
    #endif
};

#endif
