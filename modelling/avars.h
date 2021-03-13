/**
 * @file avars.h
 * @brief Data structure (alias) for auxiliary variables
 */

#ifndef AVARS_H
#define AVARS_H

#include <array>

// the ordering will be as described in the enum avar
using Avars = std::array<double, 9>; // 6 stresses and 3 heat fluxes

#endif
 

