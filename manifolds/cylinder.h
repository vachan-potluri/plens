/**
 * @file cylinder.h
 * @brief A class to apply cylindrical manifold
 */

#ifndef CYLINDER_H
#define CYLINDER_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/manifold_lib.h>

#include "manifold_description.h"

#include <iostream>

using namespace dealii;

namespace ManifoldDescriptions{

/**
 * @brief A class to apply cylindrical manifold. The inputs required are the axis direction and
 * a point on the axis. This can be used for all cylindrical geometries. E.g. truncated cone,
 * double cone.
 */
class Cylinder: public ManifoldDescription
{
    private:

    /**
     * The cylindrical manifold. Set in the constructor.
     */
    CylindricalManifold<dim> cyl_man_;

    public:

    Cylinder(const Tensor<1,dim>& axis_dir, const Point<dim>& axis_point);

    virtual void set(Triangulation<dim,dim> &triang) override;
};

} // namespace ManifoldDescriptions

#endif
