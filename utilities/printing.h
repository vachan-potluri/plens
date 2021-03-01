/**
 * @file printing.h
 * @brief Contains commonly used printing operations (generally used for debugging)
 */

#include <iostream>

#ifndef PRINTING_H
#define PRINTING_H

using state = NavierStokes::state;

namespace utilities{
    void print_state(const state &cons);
    void print_state(
        const state &cons, std::string &prefix, std::string &suffix, std::string &delim
    );
}

#endif

