/**
 * @file periodic.h
 * @brief Implements periodic BC
 */

#ifndef PERIODIC_H
#define PERIODIC_H

#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <modelling/state.h>
#include <modelling/avars.h>
#include <modelling/cavars.h>
#include <dgsem/LA.h>
#include <dgsem/dtype_aliases.h>
#include <dgsem/face_local_dof_data.h>
#include "BC.h"

#include <array>
#include <vector>
#include <map>
#include <exception>
#include <string>

#ifdef DEBUG
#include <iostream>
#include <memory>
#endif

namespace BCs
{

/**
 * @class Periodic
 * @brief Implements periodic BCs. See the documentation of BCs::BC.
 * See the note 'pens2D to plens'. The class additionally takes a vector
 * of PeriodicFacePair for construction. This can be generated using
 * GridTools::collect_periodic_faces(). Further, an integer argument is taken for construction
 * which defines which face in the face pair is the face in consideration for this BC object. This
 * class assumes the dof handler to be in
 * [standard orientation](https://www.dealii.org/current/doxygen/deal.II/DEALGlossary.html#GlossFaceOrientation).
 * Otherwise, the face dofs would require additional operations to figure out matching dofs. This
 * check is done in the constructor.
 *
 * Based on Periodic::fid and Periodic::per_paris, the constructor constructs a map
 * Periodic::cellid_to_pairid such that
 * `per_pairs[cellid_to_pairid[cell_id]].cell[fid]->index() = cell_id`.
 * This map can then be used along with Periodic::ofid_ in all the getters.
 *
 * For all the 3 stages, the ghost values are set to corresponding values on the periodic face
 * pair.
 *
 * It is mandatory that BC::g_cvars and BC::g_avars are ghosted vectors. This would enable access
 * to periodic cells/faces in case they are not owned by this process. Of course, it is also
 * mandatory that the relevant dofs of BC::g_cvars and BC::g_avars are set appropriately to include
 * these periodic ghosted cells. `add_periodicity()` function can be used for this purpose.
 */
class Periodic: public BC
{
    private:
    /**
     * 'O'ther 'f'ace id. As name suggests, ofid = 1 if fid = 0 and vice-versa.
     */
    usi ofid_;
    /**
     * Maps cell id to pair id. See the class documentation.
     */
    std::map<psize, psize> cellid_to_pairid_;

    public:
    /**
     * The vector containing periodic face pairs. See the class documentation.
     */
    const std::vector<
        GridTools::PeriodicFacePair<
            parallel::distributed::Triangulation<dim>::cell_iterator>
        >& per_pairs;
    /**
     * The 'f'ace id. Must be 0 or 1. This indicates which set of Periodic::per_pairs is in
     * consideration for this BC. See the class documentation.
     */
    const usi fid;

    Periodic(
        const DoFHandler<dim>& dh,
        const std::array<LA::MPI::Vector, 5>& gcv,
        const std::array<LA::MPI::Vector, 9>& gav,
        const std::vector<GridTools::PeriodicFacePair<
            parallel::distributed::Triangulation<dim>::cell_iterator>>& pairs,
        const usi id
    );

    private:
    void get_periodic_ldd(const FaceLocalDoFData& ldd, FaceLocalDoFData& pldd) const;

    public:
    virtual void get_ghost_stage1(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage2(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        State &cons_gh
    ) const override;
    
    virtual void get_ghost_stage3(
        const FaceLocalDoFData &ldd,
        const Tensor<1,dim> &normal,
        CAvars &ca_gh
    ) const override;

    #ifdef DEBUG
    static void test();
    #endif
};

} // namespace BCs

#endif
