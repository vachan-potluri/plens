/**
 * @file subcell_interpolator.cc
 * A class for performing subcell interpolation
 */

#include "subcell_interpolator.h"



/**
 * Constructor. Calculates the 1d weights and nodal locations
 */
SubcellInterpolator::SubcellInterpolator(
    const usi d,
    const std::array<LA::MPI::Vector, 5>& vecs,
    const std::array<LA::MPI::Vector, 5>& gh_vecs,
    const SlopeLimiter& s
):
degree(d),
dofs_per_cell((d+1)*(d+1)*(d+1)),
cdi(d),
w_1d(d+1),
node_loc_1d(d+1),
subcell_face_loc_1d(d+2),
gcrk_cvars(vecs),
gh_gcrk_cvars(gh_vecs),
slope_lim(s)
{
    QGaussLobatto<1> quad_lgl_1d(degree+1);
    for(usi i=0; i<=degree; i++){
        w_1d[i] = quad_lgl_1d.weight(i);
        node_loc_1d[i] = quad_lgl_1d.point(i)[0];
    }

    subcell_face_loc_1d[0] = 0;
    subcell_face_loc_1d[degree+2] = 1;
    for(usi i=1; i<=degree; i++){
        subcell_face_loc_1d[i] = subcell_face_loc_1d[i-1] + w_1d[i-1];
    }
}



/**
 * Calculates the conservative variable slopes within subcells. Currently, the slope for boundary
 * subcells is assumed to be zero. Can be modified later by using ghosted conservative variable
 * data available in SubcellInterpolator::gh_gcrk_cvars.
 */
void SubcellInterpolator::reinit(const DoFHandler<dim>::active_cell_iterator& cell)
{
    std::vector<psize> dof_ids(dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    std::vector<State> cell_states(dofs_per_cell);
    for(cvar var: cvar_list){
        for(usi i=0; i<dofs_per_cell; i++){
            cell_states[i][var] = gcrk_cvars[var][dof_ids[i]];
        }
    }

    for(usi dir=0; dir<dim; dir++){
        // complementary directions
        usi dir1 = (dir+1)%dim;
        usi dir2 = (dir+2)%dim;
        for(usi id1=0; id1<=degree; id1++){
            // loop over complementary direction 1
            for(usi id2=0; id2<=degree; id2++){
                // loop over complementary direction 2
                for(usi id=1; id<degree; id++){
                    // internal subcells in current direction
                    TableIndices<dim> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;
                    const usi ldof_this = cdi.tensorial_to_local(ti);
                    ti[dir] = id-1;
                    const usi ldof_left = cdi.tensorial_to_local(ti);
                    ti[dir] = id+1;
                    const usi ldof_right = cdi.tensorial_to_local(ti);

                    State cons_this = cell_states[ldof_this],
                        cons_left = cell_states[ldof_left],
                        cons_right = cell_states[ldof_right];

                    for(cvar var: cvar_list){
                        const double slope_left = (cons_this[var] - cons_left[var])/
                            (node_loc_1d[id] - node_loc_1d[id-1]);
                        const double slope_right = (cons_right[var] - cons_this[var])/
                            (node_loc_1d[id+1] - node_loc_1d[id]);
                        const double slope_avg = 0.5*(slope_left + slope_right),
                            slope_ratio = slope_left/slope_right;
                        cell_cvar_slopes[ldof_this][dir][var] =
                            slope_lim.value(slope_ratio)*slope_avg;
                    }
                }

                // now set slope for element interface subcells
                for(usi id=0; id<=degree; id+=degree){
                    TableIndices<dim> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;
                    const usi ldof = cdi.tensorial_to_local(ti);
                    for(cvar var: cvar_list){
                        cell_cvar_slopes[ldof][dir][var] = 0;
                    }
                }
            }
        }
    }
}



/**
 * Gives the left and right states for an internal subcell interface. The parameter @p ti is for
 * the subcell interface, rather than subcell itself. Its range is @f$[0,N]@f$ for complementary
 * directions and @f$[1,N]f@$ in the direction provided by the parameter @p dir.
 *
 * @note This function doesn't perform any checks on @p ti. This has to correspond to an internal
 * face and not any subcell surface that coincides with the element boundary.
 */
