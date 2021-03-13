/**
 * @file printing.h
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include <iostream>
#include <array>

#include "modelling/state.h"
#include "modelling/avars.h"

#ifndef PRINTING_H
#define PRINTING_H

namespace utilities{
    void print_state(const State &cons);
    void print_avars(const avars &a);
    
    template <int size>
    void print_array1(
        const std::array<double, size> &arr, const std::string &prefix, const std::string &suffix,
        const std::string &delim
    );
}

#endif

