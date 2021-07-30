/**
 * @file nose_cylinder.h
 * @brief A class to apply manifolds for geometries like cylinder with nose.
 */

#ifndef NOSE_CYLINDER_H
#define NOSE_CYLINDER_H

#include <iostream>

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>

#include "manifold_description.h"

using namespace dealii;

namespace ManifoldDescriptions{

/**
 * @class NoseCylinder
 * A class to apply manifolds for geometries like cylinder with nose. The common example
 * is for this is a blunted cone or blunted double cone. This class requires three information
 * pieces:
 * 1. The axis direction. For convenience, it will be assumed that the axis direction matches with
 *    the nose outward normal for a point on the axis.
 * 2. The separation or bifurcation point. Using this point, a bifurcating plane is constructed
 *    normal to the axis for demarcating between cylindrical and spherical manifold sections. It is
 *    assumed that this point lies on the axis. For cells with centers on the "outward" side of
 *    this plane, spherical manifold is applied, while the cells on "inward" side are assigned
 *    cylindrical manifold. It is asserted that no cell center lies exactly on the plane. This is
 *    not an unreasonable restriction as many meshing softwares make such meshes as separate parts
 *    and hence all the cells are generally completely on one side of the plane.
 * 3. The nose center. This is also assumed to lie on the axis. Generally, if the nose center lies
 *    to the left of bifurcation point, then the geometry following the nose is a diminishing
 *    truncated cone and if it is to the right, then the geometry following is a truncated cone
 *    with increasing radius (like in a blunted double cone geometry). If it exactly coincides with
 *    the bifurcation point, then the geometry is a combination of hemisphere and cylinder.
 *
 * The algorithm used to identify cylindrical and spherical parts is as follows:
 * - Loop over all cells
 *  - If the dot product of (a) the line joining separation point to cell center, and (b) the axis
 *    is positive
 *    - Set manifold id for spherical manifold
 *  - Else if it is negative
 *    - Set manifold id for cylindrical manifold
 *  - Else (it is 0)
 *    - Throw an exception
 *
 * @note If the nose portion contains nose center, then spherical manifold is not applied for
 * the cell containing nose center. This is because the spherical manifold functions have a
 * singularity at that point. See
 * https://www.dealii.org/current/doxygen/deal.II/classSphericalManifold.html
 * The alternative suggested is to use a TransfiniteInterpolationManifold. However, for the blunted
 * double cone geometry (which is where this class is mostly used), the nose center likes in the
 * cone part and hence this detail is avoided. An exception is thrown in case nose center is
 * observed in the spherical part. In future, this may be corrected.
 */
class NoseCylinder: public ManifoldDescription
{
    private:

    /**
     * Axis direction
     */
    const Tensor<1,3> axis_dir_;

    /**
     * Bifurcation point
     */
    const Point<3> bifurcation_point_;

    /**
     * Nose center
     */
    const Point<3> nose_center_;

    /**
     * The cylindrical manifold. Set in constructor.
     */
    CylindricalManifold<dim> cyl_man_;

    /**
     * The spherical manifold. Set in constructor.
     */
    SphericalManifold<dim> sph_man_;

    public:
    NoseCylinder(
        const Tensor<1,3>& axis_dir,
        const Point<3>& bifurcation_point,
        const Point<3>& nose_center
    );

    virtual void set(Triangulation<dim,dim> &triang) override;
};

} // namespace ManifoldDescriptions

#endif
