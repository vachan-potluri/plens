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



#ifdef DEBUG
template <int dim>
void MetricTerms<dim>::test()
{
    utilities::Testing t("MetricTerms", "class");

    t.new_block("Testing reinit() function");

    Triangulation<dim> triang;
    Point<dim> p1,p2;
    p1[0] = 0; p1[1] = 0; p1[2] = 0;
    p2[0] = 4; p2[1] = 2; p2[2] = 1;
    GridGenerator::hyper_rectangle(triang, p1, p2);

    std::function<Point<dim>(const Point<dim>&)> transformation = [=](
        const Point<dim>& p
    ){
        MetricTerms<dim>::transform(p);
    };

    GridTools::transform(transformation, triang);

    DoFHandler<dim> dof_handler(triang);
    FE_DGQ<dim> fe(2); // 2nd order polynomial
    dof_handler.distribute_dofs(fe);

    FEValues<dim> fe_values(
        fe,
        QGaussLobatto<dim>(fe.degree+1),
        update_values | update_JxW_values | update_quadrature_points
    );

    const auto &cell = dof_handler.begin_active(); // first and only cell
    fe_values.reinit(cell);

    MetricTerms<dim> mt(fe_values);
    std::cout << "JxContra_vecs:\n";
    for(usi i=0; i<fe.dofs_per_cell; i++){
        std::cout << "\tDof " << i << "\n";
        for(usi dir=0; dir<dim; dir++){
            std::cout << "\t\tDirection " << dir << ": " << mt.JxContra_vecs[i][dir] << "\n";
        }
    }
}



/**
 * Transforms point `p`. Does a simple 30 degree rotation in xy plane. Note the points are
 * transformed, not the axes.
 */
template <int dim>
Point<dim> MetricTerms<dim>::transform(const Point<dim>& p)
{
    // if point lies at origin:
    if(p.norm() == 0) return Point<dim>();

    const double PI = 3.142;
    const double rot_angle = 30*PI/180;
    const double theta_xy;
    if(p[0] == 0) theta_xy = PI/2;
    else theta_xy = atan(p[1]/p[0]);
    const double l = p.norm();

    Point<dim> p_rot;
    p_rot[2] = p[2];
    p_rot[0] = l*cos(theta_xy + rot_angle);
    p_rot[1] = l*sin(theta_xy + rot_angle);

    return p_rot;
}
#endif
