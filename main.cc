#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/navier_stokes.h"
#include "dgsem/ldof_data.h"
#include "dgsem/face_dof_info.h"
#include <iostream>

int main(){
    std::cout << "Hello, World!\n";

    #ifdef DEBUG
    cavars::test();
    NavierStokes::test();
    ldof_data::test();
    face_dof_info::test();
    #endif
    return 0;
}
