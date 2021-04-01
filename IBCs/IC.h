/**
 * @file IC.h
 * @brief Base class for initial conditions
 */

#ifndef IC_H
#define IC_H

#include <deal.II/dofs/dof_handler.h>

#include <dgsem/LA.h>
#include <modelling/state.h>
#include <modelling/navier_stokes.h>

#include <array>

using namespace dealii;
namespace ICs
{

/**
 * @class IC
 * @brief Base class for initial conditions
 *
 * This class takes a dof handler and solution vectors for construction. The solution vectors are
 * assumed to be conservative variables.
 *
 * The function IC::set() will be overridden in specific IC classes. Checks will be done on the
 * validity of IC being set using IC::assert_positivity().
 */
class IC
{
    public:
    static constexpr int dim = 3;

    /**
     * Dof handler of the problem. This will be subsequently used within set() functions of
     * implementations to get locations of dofs.
     */
    const DoFHandler<dim> &dof_handler;

    /**
     * The conservative variable vectors. The prefix `g_` is to indicate that the ordering of dofs
     * in these vectors is the global ordering.
     */
    std::array<LA::MPI::Vector, 5> &g_cvars;

    IC(const DoFHandler<dim> &dh, std::array<LA::MPI::Vector, 5> &gcv);
    virtual ~IC();
    virtual void set();
    void assert_positivity() const;
};

} // namespace ICs

#endif
