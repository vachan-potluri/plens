#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/navier_stokes.h"
#include "dgsem/face_local_dof_data.h"
#include "dgsem/face_dof_info.h"
#include "IBCs/BC.h"
#include "IBCs/free.h"
#include "IBCs/outflow.h"
#include "IBCs/uniform_inflow.h"
#include "IBCs/uniform_temp_wall.h"
#include "IBCs/symmetry.h"
#include "IBCs/periodic.h"
#include "utilities/split_string.h"
#include "IBCs/piecewise_function.h"
#include "dgsem/plens.h"
#include "dgsem/plens_test.h"
#include "dgsem/metric_terms.h"
#include "dgsem/cell_dof_info.h"
#include "dgsem/change_of_basis_matrix.h"
#include "dgsem/blender_calculator.h"
#include "dgsem/rk_coeffs.h"
#include "dgsem/dtype_aliases.h"
#include "utilities/plens_git_revision.h"

#include <iostream>



using namespace dealii;
/**
 * The main function. Constructs the relevant PLENS object and runs it. Takes two mandatory command
 * line arguments:
 * - The high order mapping degree
 * - The fe degree
 *
 * Both these are required for the construction of PLENS object.
 */
int main(int argc, char** argv){    
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    #ifdef DEBUG
    // CAvars::test();
    // NavierStokes::test();
    // FaceLocalDoFData::test();
    // FaceDoFInfo::test();
    // BCs::BC::test();
    // BCs::Free::test();
    // BCs::Outflow::test();
    // BCs::UniformInflow::test();
    // BCs::UniformTempWall::test();
    // BCs::Symmetry::test();
    // BCs::Periodic::test();
    // utilities::split_string_test();
    // ICs::PiecewiseFunction::test();
    // PLENS::test();
    // plens_test pt(); // doesn't work
    // MetricTerms<3>::test();
    // CellDoFInfo::test();
    // ChangeOfBasisMatrix<3>::test(); // template argument doesn't matter for test() as it is static
    // BlenderCalculator::test();
    // RKCoeffs::test();
    // return 0;
    #endif

    // plens_test();
    // return 0;
    if(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
        std::cout << "Git revision information:\n"
            << "\tBranch: " << PLENS_GIT_BRANCH << "\n"
            << "\tGit hash (short): " << PLENS_GIT_SHORTREV << "\n"
            << "\tGit hash (full): " << PLENS_GIT_REVISION << "\n\n";
    }
    const usi n_args = argc-1;
    AssertThrow(
        n_args == 2,
        StandardExceptions::ExcMessage(
            "Exactly two command line (mandatory) arguments are expected for mapping high order "
            "degree and FE degree."
        )
    );

    const usi mhod = std::atoi(argv[1]);
    const usi fe_degree = std::atoi(argv[2]);

    PLENS problem(mhod, fe_degree);
    problem.read_mesh();
    problem.set_NS();
    problem.set_dof_handler();
    problem.set_sol_vecs();
    problem.set_IC();
    problem.set_BC();
    problem.read_time_settings();

    problem.run();
    return 0;
}
