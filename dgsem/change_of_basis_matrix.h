/**
 * @file change_of_basis_matrix.h
 * @brief Provides a class for storing change of basis matrix
 */

#ifndef CHANGE_OF_BASIS_MATRIX_H
#define CHANGE_OF_BASIS_MATRIX_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/full_matrix.h>

/**
 * @class ChangeOfBasisMatrix
 * @brief Provides a class for storing change of basis (from DG nodal to DG modal) matrix
 *
 * This class is used for the belnding parameter calculation as part of the subcell limiter. This
 * class computes a matrix @f$C@f$ such that
 * @f[
 * \sum_{i=0}^{N_p-1} u_i \phi_i = \sum_{j=0}^{N_p-1} m_j \psi_j \\
 * m_j = \sum_{k=0}^{N_p-1} C_{jk} u_k \\
 * C_{jk} = \int_{0}^{1} \int_{0}^{1} \int_{0}^{1} \phi_j \psi_k\, dx dy dz
 * @f]
 * where @f$\{\phi_i\}@f$ are the Lagrange interpolation basis functions, @f$\psi_j@f$ are the
 * Legendre interpolation basis functions and @f$N_p=(N+1)^3@f$ is the number of polynomials.
 *
 * The Lagrange interpolation basis is provided by FE_DGQ<3> and the Legendre basis by
 * FE_DGQLegendre<3>. The ordering of basis functions for the is the same as the nodal/support
 * points' ordering. The shape function ordering in the latter is lexicographic: see the
 * documentation of TensorProductPolynomials.
 */
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
    n_poly((d+1)*(d+1)*(d+1)),
    matrix_((d+1)*(d+1)*(d+1)) // initialise to square matrix
    {}

    /**
     * Operator overload for accessing elements of ChangeOfBasisMatrix::matrix_. Returned by value.
     */
    double operator()(const usi i, const usi j) const
    {
        return matrix_(i,j);
    }
};

#endif
