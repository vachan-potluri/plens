/**
 * @file periodic.cc
 * @brief Implements periodic BC
 */

#include "periodic.h"

using namespace BCs;

/**
 * Constructor. Calls the base constructor and sets Periodic::per_pairs and Periodic::fid. If fid
 * is not 0 or 1, exception is raised. Populates Periodic::cellid_to_pairid_. Checks is the faces
 * linked through @p pairs have standard orientation. The orientation is a `std::bitset` object.
 * See https://en.cppreference.com/w/cpp/utility/bitset for details. Also see
 * https://www.dealii.org/current/doxygen/deal.II/namespaceGridTools.html#ac2a1903382c6cff07b33d456a641f6d9
 * for the meaning of `orientation` in periodic face pair.
 */
Periodic::Periodic(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav,
    const std::vector<GridTools::PeriodicFacePair<
        parallel::distributed::Triangulation<dim>::cell_iterator>>& pairs,
    const usi id
): BC(dh, gcv, gav), per_pairs(pairs), fid(id)
{
    AssertThrow(
        fid == 0 || fid == 1,
        StandardExceptions::ExcMessage(
            "The id used for constructing periodic BC must be 0 or 1"
        )
    );
    if(fid == 0) ofid_ = 1;
    else ofid_ = 0;

    psize pair_id(0);
    // check orientation and form the map between cell and pair ids
    for(const auto& pair: per_pairs){
        AssertThrow(
            pair.orientation[0] == true && // orientation
            pair.orientation[1] == false && // flip
            pair.orientation[2] == false, // rotation
            StandardExceptions::ExcMessage(
                "At least one periodic face pair doesn't have expected orientation. For periodic "
                "BC to work, the matched pairs must have same orientation, no rotation and no "
                "flip. This BC is designed to work mostly for cartesian meshes."
            )
        );
        cellid_to_pairid_[pair.cell[fid]->index()] = pair_id;
        pair_id++;
    }
}



/**
 * Given a LocalDoFData object (@p ldd), this function first checks if the data is valid. I.e.;
 * whether the cell id and face id given in @p ldd actually occur in one of Periodic::per_pairs.
 * Then, @p pldd is set such that it gives the dof linked to @p ldd through periodicity.
 *
 * In future, the assertions will probably be made only in debug mode.
 */
void Periodic::get_periodic_ldd(const LocalDoFData& ldd, LocalDoFData& pldd) const
{
    psize pair_id;
    // check if cell with given id is present in periodic face pairs
    try{
        pair_id = cellid_to_pairid_.at(ldd.cell_id);
    }
    catch(std::exception &e){
        std::string msg = "Exception caught: " + *(e.what());
        msg += ("\nThis was probably raised when trying to use a periodic face pair with invalid "
            "cell");
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(msg)
        );
    }
    // check if the face id of current pair matches with that given in ldd
    AssertThrow(
        ldd.face_id == per_pairs[pair_id].face_idx[fid],
        StandardExceptions::ExcMessage(
            "Face id given in LocalDofData doesn't match with the one stored in periodic face "
            "pair vector."
        )
    );

    // Both checks done, now set pldd
    pldd.cell_id = per_pairs[pair_id].cell[ofid_]->index();
    pldd.face_id = per_pairs[pair_id].face_idx[ofid_];
    // since both faces have exactly same orientation (asserted in ctor), the face-local dof ids
    // will be same
    pldd.face_dof_id = ldd.face_dof_id;
}



/**
 * Gets ghost conservative state using Periodic::get_periodic_ldd().
 */
void Periodic::get_ghost_stage1(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    LocalDoFData pldd;
    get_periodic_ldd(ldd, pldd);
    get_state(pldd, cons_gh);
}



/**
 * Gets ghost conservative state using Periodic::get_periodic_ldd().
 */
void Periodic::get_ghost_stage2(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    State &cons_gh
) const
{
    LocalDoFData pldd;
    get_periodic_ldd(ldd, pldd);
    get_state(pldd, cons_gh);
}



/**
 * Gets ghost conservative state and auxiliary variables using Periodic::get_periodic_ldd().
 */
void Periodic::get_ghost_stage3(
    const LocalDoFData &ldd,
    const Tensor<1,dim> &normal,
    CAvars &cav_gh
) const
{
    LocalDoFData pldd;
    get_periodic_ldd(ldd, pldd);
    get_cavars(pldd, cav_gh);
}



#ifdef DEBUG
void Periodic::test()
{
    utilities::Testing t("Periodic", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree

    std::unique_ptr<BC> bc_p = std::make_unique<Periodic>(
        bctd.dof_handler,
        bctd.gh_g_cvars,
        bctd.gh_g_avars,
        bctd.matched_pairs,
        0
    );
}
#endif
