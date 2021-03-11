#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/NavierStokes.h"
#include "dgsem/ldof_data.h"
#include <iostream>

int main(){
    std::cout << "Hello, World!\n";

    #ifdef DEBUG
    cavars::test();
    NavierStokes::test();
    ldof_data::test();
    #endif
    return 0;
}
