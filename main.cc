#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/navier_stokes.h"
#include "dgsem/local_dof_data.h"
#include "dgsem/face_dof_info.h"
#include "IBCs/BC.h"
#include "IBCs/free.h"
#include <iostream>

int main(int argc, char** argv){
    std::cout << "Hello, World!\n";
    
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    #ifdef DEBUG
    CAvars::test();
    NavierStokes::test();
    LocalDoFData::test();
    FaceDoFInfo::test();
    BCs::BC::test();
    BCs::Free::test();
    #endif
    return 0;
}
