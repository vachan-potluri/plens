/**
 * @file change_of_basis_matrix.h
 * @brief Provides a class for storing change of basis matrix
 */

#ifndef CHANGE_OF_BASIS_MATRIX_H
#define CHANGE_OF_BASIS_MATRIX_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/full_matrix.h>

#include "dtype_aliases.h"

#include <cmath>
#include <vector>

#ifdef DEBUG
#include <iostream>
#include <utilities/testing.h>
#endif

/**
 * @class ChangeOfBasisMatrix
 * @brief Provides a class for storing change of basis (from DG nodal to DG modal) matrix
 *
 * This class is used for the belnding parameter calculation as part of the subcell limiter. This
 * class computes a matrix @f$C@f$ such that
 * @f[
 * \sum_{i=0}^{N_p-1} u_i \phi_i = \sum_{j=0}^{N_p-1} m_j \psi_j \\
 * m_j = \sum_{k=0}^{N_p-1} C_{jk} u_k \\
 * C_{jk} = \int_{0}^{1} \int_{0}^{1} \int_{0}^{1} \psi_j \phi_k\, dx dy dz
 * @f]
 * where @f$\{\phi_i\}@f$ are the Lagrange interpolation basis functions, @f$\psi_j@f$ are the
 * Legendre interpolation basis functions and @f$N_p=(N+1)^3@f$ is the number of polynomials.
 *
 * The Lagrange interpolation basis is provided by FE_DGQ<3> and the Legendre basis by
 * FE_DGQLegendre<3>. The ordering of basis functions for the is the same as the nodal/support
 * points' ordering. The shape function ordering in the latter is lexicographic: see the
 * documentation of TensorProductPolynomials.
 *
 * Although this class is inteded for using only in 3d, the dimension parameter is taken as a
 * template. So this class can be used in any dimension.
 */
template <int dim>
class ChangeOfBasisMatrix
{
    public:

    /**
     * The degree of finite elements used.
     */
    const usi degree;

    /**
     * The number of shape function polynomials. Equals `pow(degree+1,3)`
     */
    const usi n_poly;

    private:

    /**
     * The actual matrix for change of basis. This is kept private and an operator overload is
     * provided by the class to access the elements of this variable.
     */
    FullMatrix<double> matrix_;

    public:

    /**
     * Constructs the matrix using FE_DGQ<1>(d) and FE_DGQLegendre<1>(d).
     */
    ChangeOfBasisMatrix(const usi d)
    :
    degree(d),
    n_poly(std::pow(d+1,dim)),
    matrix_(std::pow(d+1,dim)) // initialise to square matrix
    {
        const FE_DGQ<dim> fe_nodal(degree);
        const FE_DGQLegendre<dim> fe_modal(degree);
        const QGauss<dim> quad(degree+1); // gives exact quadrature
        const std::vector<Point<dim>> points = quad.get_points();
        const std::vector<double> weights = quad.get_weights();
        const usi n_qp = quad.size(); // number of quadrature points

        matrix_ = 0; // initialise
        for(usi i=0; i<n_poly; i++){
            for(usi j=0; j<n_poly; j++){
                // integrate the product of shape functions over [0,1]^dim
                for(usi k=0; k<n_qp; k++){
                    matrix_(i,j) += weights[k]*
                        fe_modal.shape_value(i, points[k])*
                        fe_nodal.shape_value(j, points[k]);
                } // loop over quad points
            } // loop over cols
        } // loop over rows
    }



    /**
     * Operator overload for accessing elements of ChangeOfBasisMatrix::matrix_. Returned by value.
     */
    double operator()(const usi i, const usi j) const
    {
        return matrix_(i,j);
    }



    /**
     * Returns a constant reference to the matrix stored.
     */
    const FullMatrix<double>& get() const
    {
        return matrix_;
    }



    #ifdef DEBUG
    static void test()
    {
        utilities::Testing t("ChangeOfBasisMatrix", "class");
        t.new_block("testing construction");

        const usi degree = 2;
        std::cout << "Degree: " << degree << "\n";
        const usi d = 2;
        ChangeOfBasisMatrix<d> C(degree);
        std::cout << "\n\nDimension: " << d << "\nMatrix:\n";
        for(usi row=0; row<C.n_poly; row++){
            for(usi col=0; col<C.n_poly; col++){
                std::cout << C(row,col) << "\t";
            }
            std::cout << "\n";
        }
    }
    #endif
};

#endif
