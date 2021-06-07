/**
 * @file rk_coeffs.h
 * @brief A class which gives coefficients for explicit RK updates.
 */

#ifndef RKCOEFFS_H
#define RKCOEFFS_H

#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>

#include <vector>
#include <array>

#include "dtype_aliases.h"

using namespace dealii;

/**
 * @class RKCoeffs
 * This class provides coefficients for explicit multi-step RK updates. There are two parameters
 * here: the order/degree of update algorithm and the number of steps/stages required in each
 * update. For order upto 3, the number of steps is equal to the order. For order greater than 3,
 * the number of steps are required to be greater than the order. To allow for this generalisation,
 * these two variables are kept separate.
 *
 * In each step, 3 coefficients are required:
 * 1. The coefficient of old time step solution
 * 1. The coefficient of previous RK stage solution
 * 1. The coefficient of the product of time step and residual
 *
 * These coefficients can be accessed using RKCoeffs::get().
 */
class RKCoeffs
{
    private:

    /**
     * Number of stages. Kept private, can be accessed using RKCoeffs::n_stages().
     */
    usi n_stages_;

    /**
     * The coefficients. This variable is kept private and can be accessed using RKCoeffs::get().
     * The size of this vector will be equal to the number of stages.
     */
    std::vector<
        std::array<double, 3>
    > coeffs_;

    public:

    /**
     * The order/degree of RK update.
     */
    const usi degree;



    /**
     * Constructor. Based on the order, sets the number of stages. Currently, only upto 3rd degree
     * are supported.
     */
    RKCoeffs(const usi d): degree(d)
    {
        AssertThrow(
            d >= 1 && d <= 3,
            StandardExceptions::ExcMessage(
                "Currently only upto 3rd order RK updates are supported."
            )
        );

        // Number of stages is equal to degree for upto 3rd order
        n_stages_ = d;

        // set the coefficients
        coeffs_.resize(d);
        if(d == 1){
            coeffs_ = {{1.0, 0.0 ,1.0}};
        }
        else if(d == 2){
            coeffs_ = {
                {1.0, 0.0, 1.0},
                {0.5, 0.5, 0.5}
            };
        }
        else{
            // d==3
            coeffs_ = {
                {1.0, 0.0, 1.0},
                {0.75, 0.25, 0.25},
                {1.0/3, 2.0/3, 2.0/3}
            };
        }
    }



    /**
     * Returns the number of stages
     */
    usi n_stages() const {return n_stages_;}



    /**
     * Accessor for coefficients.
     * @param[in] stage The stage of update
     * @param[in] id Must be 0, 1 or 2. See the class documentation for what these mean
     *
     * @warning No assertion on the ranges of `stage` and `id` are done. Thus if this function is
     * called with out of range values, then errors will show-up.
     */
    double get(const usi stage, const usi id) const {return coeffs_[stage][id];}
};

#endif
