#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/navier_stokes.h"
#include "dgsem/local_dof_data.h"
#include "dgsem/face_dof_info.h"
#include <iostream>

int main(){
    std::cout << "Hello, World!\n";

    #ifdef DEBUG
    CAvars::test();
    NavierStokes::test();
    LocalDoFData::test();
    face_dof_info::test();
    #endif
    return 0;
}
