/**
 * @file bc_test_data.h
 * @brief Defines a struct with some variables useful for testing ICs, similar to BCTestData
 */

#ifndef ICTESTDATA_H
#define ICTESTDATA_H

#include <deal.II/base/index_set.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>

#include <array>

#include <dgsem/LA.h>

using namespace dealii;

namespace utilities
{

/**
 * @struct ICTestData
 * This struct is a simplification of utilities::BCTestData. It provides vectors for conservative
 * variables on a unit cube constructed using the divisions and fe degree provided in the
 * constructor.
 */
struct ICTestData
{
    static constexpr int dim = 3;
    
    parallel::distributed::Triangulation<dim> triang;
    
    DoFHandler<dim> dof_handler;
    FE_DGQ<dim> fe;
    
    std::array<LA::MPI::Vector, 5> g_cvars;

    ICTestData(const usi divisions, const usi degree)
    : triang(MPI_COMM_WORLD), fe(degree), dof_handler(triang)
    {
        GridGenerator::subdivided_hyper_cube(triang, divisions); // divisions cells in each dim
        dof_handler.distribute_dofs(fe);
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        for(cvar var: cvar_list){
            g_cvars[var].reinit(locally_owned_dofs, MPI_COMM_WORLD);
        }
    }
};

} // namespace utilities

#endif
