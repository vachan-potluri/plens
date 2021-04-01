/**
 * @file IC.cc
 * @brief Base class for initial conditions
 */

#include "IC.h"

using namespace ICs;

/**
 * @brief Constructor. Sets the references IC::dof_handler and IC::g_cvars.
 * Sets IC::locally_owned_dofs_.
 */
IC::IC(const DoFHandler<dim> &dh, std::array<LA::MPI::Vector, 5> &gcv)
: dof_handler(dh), g_cvars(gcv)
{
    locally_owned_dofs_ = dof_handler.locally_owned_dofs();
}

/**
 * @brief Default virtual destructor.
 */
IC::~IC() = default;

/**
 * The main function of this class. Set to empty here. Will be overridden in implementations.
 */
void IC::set(){};

/**
 * Asserts the positivity of states set through IC::set() at all dofs. Internally uses
 * NavierStokes::assert_positivity().
 */
void IC::assert_positivity() const
{
    State cons;
    for(psize i: locally_owned_dofs_){
        for(cvar var: cvar_list) cons[var] = g_cvars[var][i];
        NavierStokes::assert_positivity(cons);
    }
}
