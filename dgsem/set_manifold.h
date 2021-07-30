/**
 * @file set_manifold.h
 * @brief Provides functions for setting different manifolds.
 * @deprecated These functions are not used anymore, classes in ManifoldDescriptions are currently
 * used.
 */

#ifndef SET_MANIFOLD_H
#define SET_MANIFOLD_H

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>

using namespace dealii;
namespace SetManifold
{

/**
 * Sets a cylindrical manifold for cylinder-flare objects, or any general rotated bodies. The
 * parameters `axis` and `axis_point` are used to generate the cylindrical manifold.
 */
inline void cylinder_flare(
    const Tensor<1,3>& axis,
    const Point<3>& axis_point,
    parallel::distributed::Triangulation<3>& triang
)
{
    CylindricalManifold<3> manifold(axis, axis_point);
    triang.set_all_manifold_ids(0);
    triang.set_manifold(0, manifold);
}



/**
 * Sets manifolds for a double cone geometry. The cone part is assigned cylindrical manifold and
 * the nose part is assigned spherical manifold. There are 3 objects required
 * @param axis[in] The cone axis direction. It is assumed that its direction points from the cone
 * to the nose.
 * @param separation_point[in] The point on the common axis which separates the conical and
 * cylindrical parts. It is assumed that a plane normal to the axis and passing through this point
 * acts as the bifurcator between the two manifolds. Consequently, it is also assumed that no cell
 * centers lie on this plane.
 * @param nose_center[in] The center for the spherical nose. This is required for generating the
 * spherical manifold.
 *
 * The algorithm used to identify cylindrical and spherical parts is as follows:
 * - Loop over all cells
 *  - If the dot product of (a) the line joining separation point to cell center, and (b) the axis
 *    is negative
 *    - Set manifold id for spherical manifold
 *  - Else
 *    - Set manifold id for cylindrical manifold
 */
inline void blunted_double_cone(
    const Tensor<1,3>& axis,
    const Point<3>& separation_point,
    const Point<3>& nose_center,
    parallel::distributed::Triangulation<3>& triang
)
{
    // manifold ids for cylindrical and spherical parts
    const int cyl_id(0), sph_id(1);

    double dotp; // dot product
    for(auto& cell: triang.active_cell_iterators()){
        dotp = scalar_product(axis, cell->center() - separation_point);
        if(dotp < 0) cell->set_all_manifold_ids(sph_id); // sphere
        else cell->set_all_manifold_ids(cyl_id); // cylinder
    } // loop over owned active cells

    CylindricalManifold<3> cyl_man(axis, separation_point);
    SphericalManifold<3> sph_man(nose_center);
    triang.set_manifold(sph_id, sph_man);
    triang.set_manifold(cyl_id, cyl_man);
}

} // namespace SetManifold

#endif
