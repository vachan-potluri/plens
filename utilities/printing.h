/**
 * @file printing.h
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include <iostream>
#include <array>

#ifndef PRINTING_H
#define PRINTING_H

using state = std::array<double, 5>;

namespace utilities{
    void print_state(const state &cons);
    void print_array1(
        const state &cons, const std::string &prefix, const std::string &suffix,
        const std::string &delim
    );
}

#endif

