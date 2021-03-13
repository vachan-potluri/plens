/**
 * @file face_dof_info.h
 * @brief A class for giving local dofs on a face
 */

#ifndef FACE_DOF_INFO_H
#define FACE_DOF_INFO_H

#include "dtype_aliases.h"

/**
 * @struct face_dof_info
 * @brief This struct contains maps that give local dofs on a face
 *
 * This data comes in handy when looping over dofs on a face. This is mostly required for surface
 * flux assembly. This class is also used for periodic boundary condition. The usage of this data
 * is best illustrated by an example.
@code
std::vector<usi> dof_ids;
for(auto &cell: dof_handler.active_cell_iterators()){
    cell->get_dof_indices(dof_ids);
    for(usi face_id=0; face_id<6; face_id++){
        for(usi face_dof=0; face_dof<fe_face.dofs_per_face; face_dof++){
            usi ldof = face_dof_info.maps[face_id][face_dof]; // get cell local id of current dof
            usi gdof = dof_ids[ldof]; // global id of current dof
            state cons; // conservative state at current dof
            for(cvar var: cvar_list) cons[var] = g_cvars[var][gdof];
        } // loop over face dofs
    } // loop over faces
} // loop over cells
@endcode
 * The maps themselves are easy to construct. See the note <b>pens2D to plens</b> and also see
 * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html) class.
 */
struct face_dof_info
{
    static constexpr usi dim = 3; // dimension
    const usi degree_;
    
    /**
     * @brief Constructor. Sets the degree of FE interpolation.
     */
    face_dof_info(const usi degree): degree_(degree) {};
};

#endif

