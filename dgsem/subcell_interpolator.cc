/**
 * @file subcell_interpolator.cc
 * A class for performing subcell interpolation
 */

#include "subcell_interpolator.h"



/**
 * Constructor. Calculates the 1d weights, nodal and subcell face locations
 */
SubcellInterpolator::SubcellInterpolator(
    const usi d,
    const std::array<LA::MPI::Vector, 5>& cvar_vecs,
    const LA::MPI::Vector& alpha_vec,
    const SlopeLimiter& s,
    const NavierStokes* ns_p
):
degree(d),
dofs_per_cell((d+1)*(d+1)*(d+1)),
cdi(d),
w_1d(d+1),
node_loc_1d(d+1),
subcell_face_loc_1d(d+2),
gcrk_cvars(cvar_vecs),
gcrk_alpha(alpha_vec),
slope_lim(s),
ns_ptr(ns_p),
cell_states(dofs_per_cell),
cell_cvar_linslopes(dofs_per_cell)
{
    QGaussLobatto<1> quad_lgl_1d(degree+1);
    for(usi i=0; i<=degree; i++){
        w_1d[i] = quad_lgl_1d.weight(i);
        node_loc_1d[i] = quad_lgl_1d.point(i)[0];
    }

    subcell_face_loc_1d[0] = 0;
    subcell_face_loc_1d[degree+1] = 1;
    for(usi i=1; i<=degree; i++){
        subcell_face_loc_1d[i] = subcell_face_loc_1d[i-1] + w_1d[i-1];
    }
}



/**
 * Calculates the conservative variable linear slopes within subcells. First, linear interpolation
 * is done from subcell values to subcell face values. These are then used to get the slope within
 * subcell. Currently, for subcell faces on element boundary, ghost data is not used.
 */
void SubcellInterpolator::reinit(const DoFHandler<dim>::active_cell_iterator& cell)
{
    alpha = gcrk_alpha[cell->global_active_cell_index()];
    std::vector<psize> dof_ids(dofs_per_cell);
    cell->get_dof_indices(dof_ids);

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
                for(usi id=0; id<=degree; id++){
                    // loop over current dir
                    TableIndices<dim> ti;
                    ti[dir] = id;
                    ti[dir1] = id1;
                    ti[dir2] = id2;

                    // get states on the left and right subcell faces
                    const usi ldof_this = cdi.tensorial_to_local(ti);
                    State cons_this = cell_states[ldof_this], cons_left, cons_right;
                    if(id == 0) cons_left = cons_this; // subcell on left element face
                    else{
                        ti[dir] = id-1;
                        // initialise to left subcell value
                        cons_left = cell_states[cdi.tensorial_to_local(ti)];
                        for(cvar var: cvar_list){
                            cons_left[var] += (cons_this[var] - cons_left[var])*
                                (subcell_face_loc_1d[id] - node_loc_1d[id-1])/
                                (node_loc_1d[id] - node_loc_1d[id-1]);
                        }
                    }
                    if(id == degree) cons_right = cons_this; // subcell on right element face
                    else{
                        ti[dir] = id+1;
                        // initialise to right subcell value
                        cons_right = cell_states[cdi.tensorial_to_local(ti)];
                        for(cvar var: cvar_list){
                            cons_right[var] -= (cons_right[var] - cons_this[var])*
                                (node_loc_1d[id+1] - subcell_face_loc_1d[id+1])/
                                (node_loc_1d[id+1] - node_loc_1d[id]);
                        }
                    }

                    // set the slope
                    for(cvar var: cvar_list){
                        cell_cvar_linslopes[ldof_this][dir][var] =
                            (cons_right[var] - cons_left[var])/
                            (subcell_face_loc_1d[id+1] - subcell_face_loc_1d[id]);
                    }
                }
            }
        }
    }
}



/**
 * Gives the left and right states for an __internal__ subcell interface. The parameter @p ti is for
 * the subcell interface, rather than subcell itself. Its range is @f$[0,N]@f$ for complementary
 * directions and @f$[1,N]f@$ in the direction provided by the parameter @p dir.
 *
 * @note This function doesn't perform any checks on @p ti. This has to correspond to an internal
 * face and not any subcell surface that coincides with the element boundary.
 *
 * @pre This function has to be called after reinit(). Else, unexpected results may occur.
 */
void SubcellInterpolator::get_left_right_states(
    const TableIndices<dim>& ti,
    const usi dir,
    State& cl,
    State& cr
)
{
    // linear interpolation weight for the subcell face
    const double Lf = (node_loc_1d[ti[dir]] - subcell_face_loc_1d[ti[dir]])/
        (node_loc_1d[ti[dir]] - node_loc_1d[ti[dir]-1]);
    
    // state and linear slope of subcell lying right of this subcell face
    const usi ldof_right = cdi.tensorial_to_local(ti);
    const State& cons_right = cell_states[ldof_right];
    const State& linslopes_right = cell_cvar_linslopes[ldof_right][dir];

    TableIndices<dim> ti_left(ti);
    ti_left[dir] -= 1;
    const usi ldof_left = cdi.tensorial_to_local(ti_left);
    // left subcell state and linear slope
    const State& cons_left = cell_states[ldof_left];
    const State& linslopes_left = cell_cvar_linslopes[ldof_left][dir];

    const double factor = std::pow(alpha, 10);
    for(cvar var: cvar_list){
        // gradient parameter
        const double r_left = 2*linslopes_left[var]*
            (node_loc_1d[ti[dir]] - node_loc_1d[ti[dir]-1])/
            (cons_right[var] - cons_left[var]) - 1;
        const double r_right = 2*linslopes_right[var]*
            (node_loc_1d[ti[dir]] - node_loc_1d[ti[dir]-1])/
            (cons_right[var] - cons_left[var]) - 1;
        const double beta_left = factor*slope_lim.value(r_left),
            beta_right = factor*slope_lim.value(r_right);
        const double w_left = beta_left*Lf + (1-beta_left), w_right = beta_right*Lf;

        cl[var] = w_left*cons_left[var] + (1-w_left)*cons_right[var];
        cr[var] = w_right*cons_left[var] + (1-w_right)*cons_right[var];
    }
}
