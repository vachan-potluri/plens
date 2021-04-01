/**
 * @file bc_test_data.h
 * @brief Defines a struct with some variables useful for testing BCs
 */

#ifndef BCTESTDATA_H
#define BCTESTDATA_H

#include <deal.II/base/index_set.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

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
 * in this namespace. These can be used in testing BC classes by including this header.
 *
 * Ghosted version of the vectors are also constructed. Further, periodicity is applied in
 * x-direction.
 */
using namespace dealii;

namespace utilities{

struct BCTestData
{
    static constexpr int dim = 3;
    
    parallel::distributed::Triangulation<dim> triang;
    
    DoFHandler<dim> dof_handler;
    FE_DGQ<dim> fe;
    
    std::array<LA::MPI::Vector, 5> g_cvars, gh_g_cvars;
    std::array<LA::MPI::Vector, 9> g_avars, gh_g_avars;

    std::vector<GridTools::PeriodicFacePair<
        parallel::distributed::Triangulation<dim>::cell_iterator>
    > matched_pairs;
    
    /**
     * @brief Constructor
     *
     * Takes the divisions and fe degree as parameters for generating subdivided hyper cube.
     * Generating a hyper cube and then refining doesn't work if periodicity is to be added. See
     * WJ-(26-27)-Mar-2021.
     */
    BCTestData(const usi divisions, const usi degree):
    triang(MPI_COMM_WORLD),
    fe(degree),
    dof_handler(triang)
    {
        GridGenerator::subdivided_hyper_cube(triang, divisions); // divisions cells in each dim

        // be default, hyper cube generates all boundaries with id 0
        GridTools::collect_periodic_faces(
            triang,
            0, // boundary id
            0, // direction,
            matched_pairs
        );
        triang.add_periodicity(matched_pairs);
        
        dof_handler.distribute_dofs(fe);
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        IndexSet locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
        
        // initialise petsc vectors
        for(cvar var: cvar_list){
            g_cvars[var].reinit(locally_owned_dofs, MPI_COMM_WORLD);
            gh_g_cvars[var].reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
        }
        for(avar var: avar_list){
            g_avars[var].reinit(locally_owned_dofs, MPI_COMM_WORLD);
            gh_g_avars[var].reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
        }
        
        // give some arbitrary values to cvars and avars
        State cons({2,2,6,4,50});
        Avars av({10,11,12,13,14,15,16,17,18});
        for(psize i: locally_owned_dofs){
            for(cvar var: cvar_list) g_cvars[var][i] = cons[var]*(i+1);
            for(avar var: avar_list) g_avars[var][i] = av[var]*(i+1);
        }

        // set ghosted vectors
        for(cvar var: cvar_list){
            g_cvars[var].compress(VectorOperation::insert);
            gh_g_cvars[var] = g_cvars[var];
        }
        for(avar var: avar_list){
            g_avars[var].compress(VectorOperation::insert);
            gh_g_avars[var] = g_avars[var];
        }
    }
};

} // namespace utilities

#endif


