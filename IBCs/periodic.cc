/**
 * @file periodic.cc
 * @brief Implements periodic BC
 */

#include "periodic.h"

using namespace BCs;

/**
 * Constructor. Calls the base constructor and sets Periodic::per_paris and Periodic::fid. If fid
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
