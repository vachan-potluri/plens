/**
 * @file nose_cylinder.cc
 * @brief A class to apply manifolds for geometries like cylinder with nose.
 */

#include "nose_cylinder.h"

using namespace ManifoldDescriptions;

/**
 * Constructor. Sets NoseCylinder::cyl_man_ using `axis_dir` and `bifurcation_point` and
 * NoseCylinder::sph_man_ using `nose_center`. See detailed documentation for description
 * of these three entities. Gives a warning if the nose center lies in spherical region.
 */
NoseCylinder::NoseCylinder(
    const Tensor<1,3>& axis_dir,
    const Point<3>& bifurcation_point,
    const Point<3>& nose_center
):
ManifoldDescription(),
axis_dir_(axis_dir),
bifurcation_point_(bifurcation_point),
nose_center_(nose_center),
cyl_man_(axis_dir, bifurcation_point),
sph_man_(nose_center)
{
    AssertThrow(
        scalar_product(axis_dir_, nose_center_-bifurcation_point_) < 0,
        StandardExceptions::ExcMessage(
            "The nose center lies in spherical region. This will produce a singularity in "
            "spherical manifold. It can be overcome using TransfiniteInterpolationManifold, but "
            "that is not currently implemented."
        )
    );
}



/**
 * Sets the manifold to `triang` according to the algorithm described in detailed documentation.
 */
void NoseCylinder::set(Triangulation<dim,dim> &triang)
{
    ConditionalOStream pcout(
        std::cout,
        (Utilities::MPI::this_mpi_process(triang.get_communicator())==0)
    );
    pcout << "Applying 'nose cylinder' manifold description\n";
    // manifold ids for cylindrical and spherical parts
    const int cyl_id(0), sph_id(1);

    double dotp; // dot product
    for(auto& cell: triang.active_cell_iterators()){
        dotp = scalar_product(axis_dir_, cell->center() - bifurcation_point_);
        if(dotp > 0){
            // sphere
            cell->set_all_manifold_ids(sph_id);
        }
        else if(dotp < 0){
            // cylinder
            cell->set_all_manifold_ids(cyl_id);
        }
        else{
            // cell center lies on bifucation plane, throw exception
            AssertThrow(
                false,
                StandardExceptions::ExcMessage(
                    "One of the cell centers in the triangulation provided lies on the "
                    "bifurcation plane. To avoid ambiguity, such cases are not allowed. See the "
                    "detailed class documentation for more information."
                )
            );
        }
    } // loop over active cells

    triang.set_manifold(sph_id, sph_man_);
    triang.set_manifold(cyl_id, cyl_man_);
}
