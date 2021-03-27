/**
 * @file periodic.cc
 * @brief Implements periodic BC
 */

#include "periodic.h"

using namespace BCs;

/**
 * Constructor. Calls the base constructor and sets Periodic::per_paris and Periodic::fid. If fid
 * is not 0 or 1, exception is raised. Also populates Periodic::cellid_to_pairid_
 */
Periodic::Periodic(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav,
    const std::vector<GridTools::PeriodicFacePair<DoFHandler<dim>::cell_iterator>>& pairs,
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
    // form the map between cell and pair ids
    for(const auto& pair: per_pairs){
        cellid_to_pairid_[pair.cell[fid]->index()] = pair_id;
        pair_id++;
    }
}



#ifdef DEBUG
void Periodic::test()
{
    utilities::Testing t("Periodic", "class");
    utilities::BCTestData bctd(2,2); // divisions and degree
}
#endif
