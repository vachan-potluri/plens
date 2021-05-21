/**
 * @file metric_terms.h
 * @brief A class for storing DGSEM metric terms of a cell
 */

#include "metric_terms.h"

/**
 * Constructor using `FEValues`. Calls reinit() internally
 */
template <int dim>
MetricTerms<dim>::MetricTerms(const FEValues<dim>& fev)
{
    reinit(fev);
}



/**
 * The main function: computes and stores all the metric terms according to `fev` provided. This
 * function first resizes JxContra_vecs and sets its data freshly. For calculation of contravariant
 * vectors, `FEValuesBase::jacobian()` is used. See WJ-21-May-2021. This function call returns a
 * `DerivativeForm` object which has useful functions from which the vectors are calculated.
 *
 * First, the covariant vectors are obtained (see appendix B of Hennemann et al 2021):
 * @f[
 * \vec{a}_i = \frac{\partial \vec{x}}{\partial \xi^i}
 * @f]
 * where @f$\xi^i@f$ is the i-th coordinate in the reference space. Then, using these covariant
 * vectors, the contravariant vectors are obtained:
 * @f[
 * J\vec{a}^i = \vec{a}_j \times \vec{a}_k,\quad (i,j,k)\ \text{cyclic}
 * @f]
 *
 * The covariant vectors are the rows of the transpose of Jacobian matrix. These can be obtained
 * as `FEValuesBase::jacobian(qid).transpose().operator(i)`.
 *
 * @note There are different forms of contravariant vectors. Even specifically for curved meshes,
 * there are different formulations of them. E.g. see Kopriva, 2006, J. of Sci. Comp. However, all
 * of them can be shown to be equivalent to this form, see WJ-21-May-2021 or the notes in TW1 dated
 * 21-May-2021.
 *
 * @pre `fev` must be reinitialised on the appropriate cell on which metric terms are desired. It
 * is also assumed that the quadrature type and FE type are correctly set to `QGaussLobatto` and
 * `FE_DGQ` respectively.
 */
template <int dim>
void MetricTerms<dim>::reinit(const FEValues<dim>& fev)
{
    // resize the metric term containers
    JxContra_vecs.resize(fev.n_quadrature_points); // dofs_per_cell == n_quadrature_points

    // calculate the contravariant vectors
    std::array<Tensor<1,dim>, dim> co_vecs; // covariant vectors
    for(usi i=0; i<fev.n_quadrature_points; i++){
        const DerivativeForm<1,dim,dim>& J_mat = fev.jacobian(i); // jacobian matrix
        const DerivativeForm<1,dim,dim> J_mat_T = J_mat.transpose();

        // get covariant vecs
        for(usi dir=0; dir<dim; dir++) co_vecs[dir] = J_mat_T[dir];

        // calculate contravariant vecs
        // https://codereview.stackexchange.com/questions/57923/index-into-array-as-if-it-is-circular
        for(usi dir=0; dir<dim; dir++){
            // cyclic indices
            usi dir1 = (dir+1)%dim;
            usi dir2 = (dir+2)%dim;

            JxContra_vecs[i][dir] = cross_product_3d(co_vecs[dir1], co_vecs[dir2]);
        }
    } // loop over dofs == quad points
}
