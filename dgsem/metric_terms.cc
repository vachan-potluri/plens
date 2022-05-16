/**
 * @file metric_terms.h
 * @brief A class for storing DGSEM metric terms of a cell
 */

#include "metric_terms.h"

/**
 * Constructor using `FEValues`. Calls reinit() internally
 */
template <int dim>
MetricTerms<dim>::MetricTerms(const FEValues<dim>& fev, const FullMatrix<double>& Q)
{
    reinit(fev, Q);
}



/**
 * The main function: computes and stores all the metric terms according to `fev` provided. This
 * function first resizes MetricTerms::JxContra_vecs and sets its data freshly. For calculation of
 * contravariant vectors, `FEValuesBase::jacobian()` is used. See WJ-21-May-2021. This function
 * call returns a `DerivativeForm` object which has useful functions from which the vectors are
 * calculated.
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
 * The parameter `Q` is used for calculation of subcell normals. Its values must be consistent with
 * `fev`. No checks on this are done. The formula for normal calculation is (for `dir=0`)
 * @f[
 * \vec{n}_{(i-1,i)jk} = (J\vec{a}^1)_{0jk} +
 *      \sum_{l=0}^{i-1} \sum_{m=0}^{N} Q_{lm} (J\vec{a}^1)_{mjk}
 * @f]
 * When @f$i=0@f$ in above formula, we get @f$\vec{n}_{(L,0)jk}@f$ and when @f$i=N+1@f$, we get
 * @f$\vec{n}_{(N,R)jk}@f$. The formula is valid for both these end cases also.
 *
 * @pre `fev` must be reinitialised on the appropriate cell on which metric terms are desired. It
 * is also assumed that the quadrature type and FE type are correctly set to `QGaussLobatto` and
 * `FE_DGQ` respectively.
 * @pre Values in `Q` must be consistent with `fev`. It must be the product of 1D weights diagonal
 * matrix and the 1D differentiation matrix:
 * @f$Q_{ij} = w_i D_{ij} = w_i \frac{\partial l_j}{\partial \xi}(\xi_i)@f$
 * where @f$w@f$ and @f$l@f$ are the weights and shape functions corresponding to 1D LGL
 * quadrature.
 *
 * @note The current implementation is mathematically correct only for straight meshes. For curved
 * meshes, an interpolant of the coordinate transformation has to be used rather than the function
 * itself. See WJ-16-Mar-2022.
 */
template <int dim>
void MetricTerms<dim>::reinit(const FEValues<dim>& fev, const FullMatrix<double>& Q)
{
    // resize the metric term containers
    JxContra_vecs.resize(fev.n_quadrature_points); // dofs_per_cell == n_quadrature_points
    detJ.resize(fev.n_quadrature_points);
    Jinv.resize(fev.n_quadrature_points);
    const usi degree = fev.get_fe().degree;
    for(usi dir=0; dir<dim; dir++){
        TableIndices<dim> ti(degree+1,degree+1,degree+1);
        ti[dir] += 1; // one size extra in direction `dir`
        subcell_normals[dir].reinit(ti);
    }

    CellDoFInfo cdi(degree);

    // calculate the contravariant vectors
    std::array<Tensor<1,dim>, dim> co_vecs; // covariant vectors
    for(usi i=0; i<fev.n_quadrature_points; i++){
        const DerivativeForm<1,dim,dim>& J_mat = fev.jacobian(i); // jacobian matrix
        detJ[i] = J_mat.determinant();
        Jinv[i] = fev.inverse_jacobian(i);

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

    // subcell normals
    // here, `id*` is the tensor index corresponding to direction `dir*`
    for(usi dir=0; dir<dim; dir++){
        // other directions
        usi dir1 = (dir+1)%dim;
        usi dir2 = (dir+2)%dim;

        // loop over tensor indices
        for(usi id=0; id<=degree+1; id++){
            for(usi id1=0; id1<=degree; id1++){
                for(usi id2=0; id2<=degree; id2++){
                    TableIndices<dim> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;
                    // indices corrsponding to the first term in normal expression
                    TableIndices<dim> ti0;
                    ti0[dir] = 0;
                    ti0[dir1] = id1;
                    ti0[dir2] = id2;

                    // initialise normal
                    subcell_normals[dir](ti) = JxContra_vecs[cdi.tensorial_to_local(ti0)][dir];

                    for(usi l=0; l<=id-1; l++){
                        for(usi m=0; m<=degree; m++){
                            TableIndices<dim> tim;
                            tim[dir] = m;
                            tim[dir1] = id1;
                            tim[dir2] = id2;
                            subcell_normals[dir](ti) += Q(l,m)*
                                JxContra_vecs[cdi.tensorial_to_local(tim)][dir];
                        }
                    }
                } // loop over tensor index
            } // loop over tensor index
        } // loop over tensor index
    } // loop over direction
}



#ifdef DEBUG
template <int dim>
void MetricTerms<dim>::test()
{
    utilities::Testing t("MetricTerms", "class");

    t.new_block("Testing reinit() function");

    // create test triangulation
    Triangulation<dim> triang;
    Point<dim> p1,p2;
    p1[0] = 0; p1[1] = 0; p1[2] = 0;
    p2[0] = 4; p2[1] = 2; p2[2] = 1;
    GridGenerator::hyper_rectangle(triang, p1, p2);
    std::function<Point<dim>(const Point<dim>&)> transformation = MetricTerms<dim>::transform;
    GridTools::transform(transformation, triang);

    DoFHandler<dim> dof_handler(triang);
    FE_DGQ<dim> fe(2); // 2nd order polynomial
    dof_handler.distribute_dofs(fe);

    // create fevalues object
    FEValues<dim> fe_values(
        fe,
        QGaussLobatto<dim>(fe.degree+1),
        update_values | update_jacobians | update_quadrature_points
    );
    const auto &cell = dof_handler.begin_active(); // first and only cell
    fe_values.reinit(cell);

    // form Q matrix
    std::vector<double> w_1d(fe.degree+1);
    FullMatrix<double> D(fe.degree+1), Q(fe.degree+1);
    FE_DGQ<1> fe_1d(fe.degree);
    QGaussLobatto<1> quad_lgl_1d(fe.degree+1);
    for(usi i=0; i<=fe.degree; i++){
        w_1d[i] = quad_lgl_1d.weight(i);
    }
    const std::vector<Point<1>> &points_1d = quad_lgl_1d.get_points();
    for(usi row=0; row<=fe.degree; row++){
        for(usi col=0; col<=fe.degree; col++){
            D(row,col) = fe_1d.shape_grad(col, points_1d[row])[0];
            Q(row,col) = w_1d[row]*D(row,col);
        }
    }

    // invoke and print
    MetricTerms<dim> mt(fe_values, Q);
    std::cout << "JxContra_vecs:\n";
    for(usi i=0; i<fe.dofs_per_cell; i++){
        std::cout << "\tDof " << i << "\n";
        for(usi dir=0; dir<dim; dir++){
            std::cout << "\t\tDirection " << dir << ": " << mt.JxContra_vecs[i][dir] << "\n";
        }
    }

    std::cout << "Determinants:\n";
    for(usi i=0; i<fe.dofs_per_cell; i++){
        std::cout << "\tDof " << i << " : " << mt.detJ[i] << "\n";
    }

    std::cout << "Subcell normals:\n";
    for(usi dir=0; dir<dim; dir++){
        std::cout << "\tDirection " << dir << "\n";
        usi dir1 = (dir+1)%dim;
        usi dir2 = (dir+2)%dim;

        for(usi id=0; id<=fe.degree+1; id++){
            for(usi id1=0; id1<=fe.degree; id1++){
                for(usi id2=0; id2<=fe.degree; id2++){
                    TableIndices<dim> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;

                    std::cout << "\t\tTable indices " << ti << " : "
                        << mt.subcell_normals[dir](ti) << "\n";
                }
            }
        }
    }

    // Expected result (assuming the triang box has widths 4, 2 and 1)
    // Ja^1 = (2*cos r, 2*sin r, 0)
    // Ja^2 = (-4*sin r, 4*cos r, 0)
    // Ja^3 = (0, 0, 8)
    // where r is the rotation angle in transform function (assuming transform does 2d rotation
    // in xy plane)

    // For normals, all subcell normals must match with Ja^{1/2/3} depending on the direction
}



/**
 * Transforms point `p`. Does a simple 30 degree rotation in xy plane. Note the points are
 * transformed, not the axes.
 */
template <int dim>
Point<dim> MetricTerms<dim>::transform(const Point<dim>& p)
{
    const double PI = 3.141592653;
    const double rot_angle = 30*PI/180;

    Point<dim> p_rot;
    p_rot[2] = p[2];
    p_rot[0] = p[0]*std::cos(rot_angle) - p[1]*std::sin(rot_angle);
    p_rot[1] = p[0]*std::sin(rot_angle) + p[1]*std::cos(rot_angle);

    return p_rot;
}
#endif

#include "metric_terms.inst"
