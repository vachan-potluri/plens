#include "modelling/vars.h"
#include "modelling/NavierStokes.h"
#include <iostream>

int main(){
    std::cout << "Hello, World!\n";

    #ifdef DEBUG
    NavierStokes::test();
    #endif
    return 0;
}
