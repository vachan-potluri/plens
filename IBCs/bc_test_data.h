/**
 * @file bc_test_data.h
 * @brief Defines a struct with some variables useful for testing BCs
 */

#ifndef BCTEST_H
#define BCTEST_H

#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/index_set.h>

#include <array>

#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <modelling/var_enums.h>
#include <modelling/state.h>
#include <modelling/avars.h>

/**
 * @struct BCTestData
 * @brief Contains for variables useful for testing BCs
 *
 * A triangulation, dof handler and sample conservative and auxiliary variable vectors are stored
 * in this namespace. These can be used in testing BC classes by including this header
 */
using namespace dealii;

namespace utilities{

struct BCTestData
{
    static constexpr int dim = 3;
    
    parallel::distributed::Triangulation<dim> triang;
    
    DoFHandler<dim> dof_handler;
    FE_DGQ<dim> fe;
    
    std::array<LA::MPI::Vector, 5> g_cvars;
    std::array<LA::MPI::Vector, 9> g_avars;
    
    /**
     * @brief Constructor
     *
     * Takes the hyper cube refinement and fe degree as parameters
     */
    BCTestData(const usi refinement, const usi degree):
    triang(MPI_COMM_WORLD),
    fe(degree),
    dof_handler(triang)
    {
        GridGenerator::hyper_cube(triang);
        triang.refine_global(refinement); // 2**refinement cells in each direction
        
        dof_handler.distribute_dofs(fe);
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        
        // initialise petsc vectors
        for(cvar var: cvar_list){
            g_cvars[var].reinit(locally_owned_dofs, MPI_COMM_WORLD);
        }
        for(avar var: avar_list){
            g_avars[var].reinit(locally_owned_dofs, MPI_COMM_WORLD);
        }
        
        // give some arbitrary values to cvars and avars
        State cons({2,2,6,4,50});
        Avars av({10,11,12,13,14,15,16,17,18});
        for(psize i: locally_owned_dofs){
            for(cvar var: cvar_list) g_cvars[var] = cons[var]*(i+1);
            for(avar var: avar_list) g_avars[var] = av[var]*(i+1);
        }
    }
};

} // namespace utilities

#endif


