/**
 * @file dtype_aliases.h
 * @brief Defines different aliases for data types that are used commonly.
 */

#ifndef DTYPE_ALIASES_H
#define DTYPE_ALIASES_H

#include <deal.II/base/types.h>

/// A size data type for 'p'lens. Generally used for cell index or global dof index. See
/// https://www.dealii.org/current/doxygen/deal.II/GlobalDoFIndex.html
using psize = dealii::types::global_dof_index;

/// Generally used for face id w.r.t, cell and dof id w.r.t. cell/face
using usi = unsigned short int;

#endif

