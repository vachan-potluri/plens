/**
 * @file face_local_dof_data.h
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
 * @struct FaceLocalDoFData
 * @brief This struct is like a container to hold a face's local dof data
 *
 * The struct contains
 * 1. Cell id
 * 2. Face id wrt cell @f$[0,6)@f$
 * 3. DoF id wrt face @f$[0, (N+1)^2)@f$
 *
 * This is mainly used in BCs::Periodic or for setting spatially varying BC (not yet implemented).
 * Often, this is used in conjunction with FaceDoFInfo and DoFHandler<dim>::cell_iterator to get
 * the global dof from the stored local data. The cell id is generally `cell->index()` and hence
 * this class is of use only when
 * 1. You want to index cells that are strictly locally owned
 * 2. You want to use meshes that don't involve refinement
 */
struct FaceLocalDoFData
{
    psize cell_id;
    usi face_id, face_dof_id;
    
    /**
     * @brief Constructor taking all values
     *
     * @pre @p f @f$\in [0,6)@f$, @p dof_id @f$\in [0,(N+1)^2)@f$. No assertion is done.
     */
    FaceLocalDoFData(const psize &c, const usi f, const usi d):
        cell_id(c), face_id(f), face_dof_id(d)
    {}



    /**
     * @brief Default constructor
     */
    FaceLocalDoFData() = default;
    
    
    
    #ifdef DEBUG
    static void test()
    {
        utilities::Testing t("FaceLocalDoFData", "struct");
        
        {
            t.new_block("testing complete construction");
            FaceLocalDoFData ldd(1000, 5, 3);
            std::cout << ldd.cell_id << " " << ldd.face_id << " " << ldd.face_dof_id << "\n";
        }

        {
            t.new_block("testing default construction");
            FaceLocalDoFData ldd;
            std::cout << ldd.cell_id << " " << ldd.face_id << " " << ldd.face_dof_id << "\n";
        }
        
        {
            t.new_block("testing access");
            const FaceLocalDoFData ldd(1000, 5, 3);
            // ldd.cell_id = 100; // error
        }
    }
    #endif
};

#endif

