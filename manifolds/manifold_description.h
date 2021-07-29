/**
 * @file manifold_description.h
 * @brief A base class for classes that describe a generic manifold
 */

#ifndef MANIFOLD_DESCRIPTION_H
#define MANIFOLD_DESCRIPTION_H

#include <string>

#include <deal.II/grid/tria.h>

using namespace dealii;

/**
 * @class ManifoldDescription
 * @brief A base class for classes that describe a generic manifold and provide functions for
 * applying them on triangulations. This is an abstract class and the derivatives of this class do
 * the actual job. The need for this class is obvious. Often, it is not possible to describe an
 * entire triangulation using a single Manifold. In that case, this class allows for inheritance
 * in storing and using different information for different kinds of manifold descriptions.
 *
 * The function set() provided by this base class is overridden by derived classes. Three
 * dimensionality is assumed. The derived classes take the other required information and sotre
 * it within them.
 */
class ManifoldDescription
{
    public:

    /**
     * Dimensions.
     */
    static constexpr int dim = 3;

    /**
     * Blank constructor, does nothing.
     */
    ManifoldDescription() = default;

    /**
     * Pure virtual destructor.
     */
    virtual ~ManifoldDescription() = 0;

    /**
     * Pure virtual setter. This is the main function which applies the manifold(s) to the
     * triangulation provided. Here, the base triangulation class is used.
     */
    virtual void set(Triangulation<dim,dim> &triang) = 0;
};

#endif
