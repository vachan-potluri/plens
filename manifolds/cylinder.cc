/**
 * @file cylinder.h
 * @brief A class to apply cylindrical manifold
 */

#include "cylinder.h"

using namespace ManifoldDescriptions;

/**
 * Constructor. Forms Cylinder::cyl_man_ using the give parameters.
 * @param[in] axis_dir The axis direction
 * @param[in] axis_point A point on the axis
 */
Cylinder::Cylinder(const Tensor<1,dim>& axis_dir, const Point<dim>& axis_point)
: ManifoldDescription(), cyl_man_(axis_dir, axis_point)
{}



/**
 * Sets cylindrical manifold to all entities of `triang`.
 */
void Cylinder::set(Triangulation<dim,dim> &triang)
{
    triang.set_all_manifold_ids(0);
    triang.set_manifold(0, cyl_man_);
}
