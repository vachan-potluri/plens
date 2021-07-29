/**
 * @file printing.cc
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include "printing.h"

/**
 * @brief Calls utilities::print_array1<5>() with a specific prefix, delimiter and suffix.
 */
void utilities::print_state(const State &cons)
{
    print_array1<5>(cons, "", "\n", ", ");
}



/**
 * @brief Calls utilities::print_array1<9>() with a specific prefix, delimiter and suffix.
 */
void utilities::print_avars(const Avars &a)
{
    print_array1<9>(a, "", "\n", ", ");
}



/**
 * @brief Prints a 1D std::array. The size is taken as a template argument. The print format is
@verbatim
<prefix><element 0><delimiter>...<suffix>
@endverbatim
 */
template <int size>
void utilities::print_array1(
    const std::array<double, size> &arr, const std::string &prefix, const std::string &suffix,
    const std::string &delim
)
{
    std::cout << prefix;
    for(int i=0; i<size; i++){
        printf("%e%s", arr[i], delim.c_str());
    }
    std::cout << suffix;
}

