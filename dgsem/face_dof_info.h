/**
 * @file face_dof_info.h
 * @brief A class for giving local dofs on a face
 */

#ifndef FACE_DOF_INFO_H
#define FACE_DOF_INFO_H

#include "dtype_aliases.h"

#ifdef DEBUG
#include <iostream>
#include <utilities/testing.h>
#endif

/**
 * @class FaceDoFInfo
 * @brief This class contains maps that map face-local dof id to cell-local dof id
 *
 * Simply put, for FE_DGQ<3> elements, this class provides 6 arrays. These arrays consist the
 * cell-locally ordered dofs lying on the faces of the element. And then it also provides inverse
 * maps so that given a face id in a cell, the face-local dof numbering can be converted back to
 * cell-local dof numbering.
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
            usi ldof = FaceDoFInfo.maps[face_id][face_dof]; // get cell local id of current dof
            usi gdof = dof_ids[ldof]; // global id of current dof
            state cons; // conservative state at current dof
            for(cvar var: cvar_list) cons[var] = g_cvars[var][gdof];
        } // loop over face dofs
    } // loop over faces
} // loop over cells
@endcode
 * The maps themselves are easy to construct. See the note <b>pens2D to plens</b> and also see
 * [FE_DGQ](https://www.dealii.org/current/doxygen/deal.II/classFE__DGQ.html) class.
 *
 * @todo Use `dealii::Table<N,T>` instead of nested arrays. CellDoFInfo does this. 22-May-2021.
 */
class FaceDoFInfo
{
    public:
    static constexpr usi dim = 3; // dimension
    const usi degree;
    /**
     * @brief Forward maps. Access: `maps[face id][face-local dof id]`. The result is cell-local dof
     * id.
     *
     * @waring This is a non-const public variable, do not accidentally modify its values.
     */
    std::array<std::vector<usi>, 2*dim> maps;
    
    /**
     * @brief Inverse maps. Access: `maps[face id][cell-local dof id]`. The result is face-local dof
     * id.
     *
     * @warning This is a non-const public variable. Hence it is recommended that entries of this
     * map are accessed using `std::map::at()` rather than `std::map::operator[]` because the
     * latter can accidentally insert elements is used carelessly
     */
    std::array<std::map<usi, usi>, 2*dim> inverse_maps;
    
    /**
     * @brief Constructor. Sets the degree of FE interpolation.
     */
    FaceDoFInfo(const usi deg): degree(deg)
    {
        form_maps();
    }
    
    private:
    /**
     * @brief Populates maps and inverse_maps.
     *
     * An example of maps for 3rd order (degree, @f$N=2@f$) element is given here
     * @f[
     * D_0 = \{ 0,3,6,9,12,15,18,21,24 \}\\
     * D_1 = \{ 2,5,8,11,14,17,20,23,26 \}\\
     * D_2 = \{ 0,1,2,9,10,11,18,19,20 \}\\
     * D_3 = \{ 6,7,8,15,16,17,24,25,26 \}\\
     * D_4 = \{ 0,1,2,3,4,5,6,7,8 \}\\
     * D_5 = \{ 18,19,20,21,22,23,24,25,26 \}\\
     * @f]
     */
    void form_maps()
    {
        const usi dofs_per_face = (degree+1)*(degree+1);
        
        // set size for maps; not reqd for inverse_maps
        for(usi face_id=0; face_id<2*dim; face_id++) maps[face_id].resize(dofs_per_face);
        
        usi face_dof, cell_dof; // ids of current dof w.r.t. face and cell
        
        // Face 0
        for(face_dof=0; face_dof<dofs_per_face; face_dof++){
            cell_dof = (degree+1)*face_dof;
            maps[0][face_dof] = cell_dof;
            inverse_maps[0][cell_dof] = face_dof;
        }
        
        // Face 1
        for(face_dof=0; face_dof<dofs_per_face; face_dof++){
            cell_dof = degree + (degree+1)*face_dof;
            maps[1][face_dof] = cell_dof;
            inverse_maps[1][cell_dof] = face_dof;
        }
        
        // Face 2
        face_dof = 0;
        for(usi i=0; i<=degree; i++){
            for(usi j=0; j<=degree; j++){
                cell_dof = i*dofs_per_face + j;
                maps[2][face_dof] = cell_dof;
                inverse_maps[2][cell_dof] = face_dof;
                face_dof++;
            }
        }
        
        // Face 3
        face_dof = 0;
        for(usi i=0; i<=degree; i++){
            for(usi j=0; j<=degree; j++){
                cell_dof = (degree+1)*degree + i*dofs_per_face + j;
                maps[3][face_dof] = cell_dof;
                inverse_maps[3][cell_dof] = face_dof;
                face_dof++;
            }
        }
        
        // Face 4
        for(face_dof=0; face_dof<dofs_per_face; face_dof++){
            cell_dof = face_dof;
            maps[4][face_dof] = cell_dof;
            inverse_maps[4][cell_dof] = face_dof;
        }
        
        // Face 5
        for(face_dof=0; face_dof<dofs_per_face; face_dof++){
            cell_dof = degree*dofs_per_face + face_dof;
            maps[5][face_dof] = cell_dof;
            inverse_maps[5][cell_dof] = face_dof;
        }
    }
    
    
    
    #ifdef DEBUG
    private:
    void print_maps()
    {
        const usi dofs_per_face = (degree+1)*(degree+1);
        
        std::cout << "Maps\n";
        for(usi face_id=0; face_id<6; face_id++){
            std::cout << "Face " << face_id << ": ";
            for(usi i=0; i<dofs_per_face; i++) std::cout << maps[face_id][i] << " ";
            std::cout << "\n";
        }
        
        std::cout << "Inverse maps\n";
        for(usi face_id=0; face_id<6; face_id++){
            std::cout << "Face " << face_id << ": ";
            for(const std::pair<usi, usi> &p: inverse_maps[face_id])
                std::cout << p.first << "-" << p.second << " ";
            std::cout << "\n";
        }
    }
    
    public:
    static void test()
    {
        utilities::Testing t("class", "FaceDoFInfo");
        
        {
            t.new_block("degree 1 maps");
            FaceDoFInfo fdi(1);
            fdi.print_maps();
        }
        
        {
            t.new_block("degree 2 maps");
            FaceDoFInfo fdi(2);
            fdi.print_maps();
        }
        
        {
            t.new_block("degree 3 maps");
            FaceDoFInfo fdi(3);
            fdi.print_maps();
        }
    }
    #endif
};

#endif

