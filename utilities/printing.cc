/**
 * @file printing.h
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include "printing.h"



void utilities::print_state(const state &cons)
{
    print_state(cons, "\n", "\n", ", ");
}



void utilities::print_state(
    const state &cons, std::string &prefix, std::string &suffix, std::string &delim
)
{
    std::cout << prefix;
    for(double x: cons){
        printf("%e%s", x, delim);
    }
    std::cout << suffix;
}

