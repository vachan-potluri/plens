/**
 * @file state.h
 * @brief Defines a data structure for conservative/primitive state/flux
 */

#include <array>

#include <deal.II/base/tensor.h>

#include "var_enums.h"

#ifndef STATE_H
#define STATE_H

using State = std::array<double, 5>; // for conservative/primitive state array/vector

#endif

