#include "modelling/var_enums.h"
#include "modelling/cavars.h"
#include "modelling/NavierStokes.h"
#include <iostream>

int main(){
    std::cout << "Hello, World!\n";

    #ifdef DEBUG
    cavars::test();
    //NavierStokes::test();
    #endif
    return 0;
}
