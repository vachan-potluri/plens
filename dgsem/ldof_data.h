/**
 * @file ldof_data.h
 * @brief A data structure to contain all the info required to describe a dof locally
 */

#ifndef LDOF_DATA_H
#define LDOF_DATA_H

#include "dtype_aliases.h"

#ifdef DEBUG
#include <iostream>
#include <utilities/testing.h>
#endif
/**
 * @struct ldof_data
 * @brief This struct is like a container to hold 'l'ocal dof data
 *
 * The struct contains
 * 1. Cell id
 * 2. Face id wrt cell @f$[0,6)@f$
 * 3. DoF id wrt face @f$[0, (N+1)^2)@f$
 *
 * This is mainly used in BC class for periodic BC setting or for setting spatially varying BC.
 */
struct ldof_data
{
    psize cell_id;
    usi face_id, dof_id;
    
    /**
     * @brief Constructor
     */
    ldof_data(const psize &c, const usi f, const usi d): cell_id(c), face_id(f), dof_id(d) {}
};

#endif

