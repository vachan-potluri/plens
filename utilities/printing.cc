/**
 * @file printing.h
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include "printing.h"



void utilities::print_state(const state &cons)
{
    print_array1(cons, "\n", "\n", ", ");
}


/**
 * @brief Prints a 1D std::array
 */
void utilities::print_array1(
    const state &arr, const std::string &prefix, const std::string &suffix,
    const std::string &delim
)
{
    int N = arr.size();
    std::cout << prefix;
    for(int i=0; i<N; i++){
        printf("%e%s", arr[i], delim.c_str());
    }
    std::cout << suffix;
}

