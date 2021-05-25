/**
 * @file change_of_basis_matrix.h
 * @brief Provides a class for storing change of basis matrix
 */

#ifndef CHANGE_OF_BASIS_MATRIX_H
#define CHANGE_OF_BASIS_MATRIX_H

/**
 * @class ChangeOfBasisMatrix
 * @brief Provides a class for storing change of basis (from DG nodal to DG modal) matrix
 *
 * This class is used for the belnding parameter calculation as part of the subcell limiter. This
 * class computes a matrix @f$C@f$ such that
 * @f[
 * \sum_{i=0}^{N} u_i \phi_i = \sum_{j=0}^{N} m_j \psi_j \\
 * m_j = \sum_{k=0}^{N} C_{jk} u_k \\
 * C_{jk} = \int_{0}^{1} \phi_j \psi_k\, dx
 * @f]
 * where @f$\{\phi_i\}@f$ are the Lagrange interpolation basis functions and @f$\psi_j@f$ are the
 * Legendre interpolation basis functions.
 *
 * Although the matrix itself is calculated using 1D basis functions, since the 3D basis functions
 * of both FE_DGQ and FE_DGQLegendre are tensor product polynomials, this one-dimensional matrix
 * can also be used to convert basis from FE_DGQ<3> to FE_DGQLegendre<3>. For this, one would
 * require converting a 3D basis function index into 3 appropriate tensor indices. This can be
 * readily done using the CellDoFInfo class since the indexing of shape functions and the indexing
 * of support/nodal points is the same for FE_DGQ<3>.
 */
class ChangeOfBasisMatrix
{};

#endif
