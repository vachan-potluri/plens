/**
 * @file metric_terms.h
 * @brief A class for storing DGSEM metric terms of a cell
 */

#ifndef METRIC_TERMS_H
#define METRIC_TERMS_H

#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/full_matrix.h>

#include <vector>
#include <array>

#include "dtype_aliases.h"
#include "cell_dof_info.h"

#ifdef DEBUG
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <iostream>
#include <cmath>
#include <functional>
#include <utilities/testing.h>
#endif

using namespace dealii;

/**
 * @class MetricTerms
 * This class computes and stores metric terms required for DGSEM algorithm for a cell. Using this
 * class is simple. Construct the class by passing an appropriate `FEValues` object. The
 * constructor internally calls MetricTerms::reinit(). When making this call, the `FEValues`
 * object must also have been re-initialised on the desired cell. A plain constructor is also
 * provided in case MetricTerms::reinit() is to be called separately.
 *
 * The most common way of using this with a `DoFHandler` is to use a `std::map` of `MetricTerms`
 * objects using the cell index as the key. This will however not work for adaptively refined
 * meshes where a cell index is not sufficient to uniquely identify a cell. However, we are not
 * interested in such cases and for our purposes, a simple map will suffice.
 *
 * The metric terms currently calculated and stored are (update the list when required)
 * - Contravariant vectors
 * - Subcell normals
 *
 * Access to stored quantities is provided directly through public variables. Hence it might be
 * dangerous to use non-const version of this object. However, a const version can be readily
 * constructed by passing an appropriate `FEValues` object for construction.
 *
 * @note Although this class is constructed from a template, only `dim=3` is currently supported.
 * The file `metric_terms.inst` instantiates this template class only for `dim=3`.
 */
template <int dim>
class MetricTerms
{
    public:

    /**
     * A variable holding the contravariant vectors. Access:
     * `JxContra_vecs[cell-local dof id][direction]`
     * gives contravariant vector in the direction `direction`.
     * Note here that cell-local dof ordering is used, rather than a tensor-product ordering.
     */
    std::vector<std::array<Tensor<1,dim>, dim>> JxContra_vecs;

    /**
     * An array holding the Jacobian determinants. The ordering is cell-local again, as against
     * tensorial ordering. Access:
     * `detJ[cell-local dof id]`
     */
    std::vector<double> detJ;

    /**
     * Subcell normal vectors. There is one important distinction between subcell normals and cell
     * normals as used in PLENS. The subcell normals are uni-directional: their doesn't change
     * depending on the subcell. This implies that the normal will point _into_ the subcell on the
     * "left" subcell face and will point _outward_ the subcell on the "right" subcell face.
     * There is one more point to note: these vectors are not unit vectors. These are actually the
     * area vectors, scaled by some weight factor. See TW1 notes 20-May-2021 or WJ-20-May-2021. So,
     * if a subcell shares the cell's face, then the normal given by this variable will point in
     * the direction of outward face normal if the face index is 1, 3 or 5. If the face index is 0,
     * 2 or 4, then the subcell normal and the face outward normal point in opposite directions.
     * Because this deals with the subcell, the ordering has to be tensorial. Access:
     * `subcell_normals[dir][i][j][k]`
     * gives @f$\vec{n}_{(i-1,i)jk}@f$ for `dir=0` and similarly for `dir=1` and `dir=2`. For
     * `dir=0`, `i=0,1,...,N+1` and `j,k=0,1,...,N`. With `dir=0`, when `i=0`, this gives
     * @f$\vec{n}_{(L,0)jk}@f$ and when `i=N+1`, this gives @f$\vec{n}_{(N,R)jk}@f$.
     */
    std::array< Table<3, Tensor<1,dim>>, dim > subcell_normals;

    /**
     * Plain constructor. Does nothing.
     */
    MetricTerms(){}
    MetricTerms(const FEValues<dim>&, const FullMatrix<double>&);
    void reinit(const FEValues<dim>&, const FullMatrix<double>&);

    #ifdef DEBUG
    static void test();
    static Point<dim> transform(const Point<dim>&);
    #endif
};

#endif
