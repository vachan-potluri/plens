/**
 * @file piecewise_function.h
 * @brief Sets initial condition based on piecewise defined function.
 */

#include "piecewise_function.h"

using namespace ICs;

/**
 * Returns the piece id based on the point @p p provided. The function is extended directly from
 * pens2D.
 * 
 * The algorithm used here uses the concept of unidirectional piece id. As the name suggests, it
 * is calculated independently for all directions. Its value is equal to the id of the right
 * interface in a direction of the cell identified. The resulting piece id is then obtained using
 * these unidirectional ids.
 * 
 * Further, to get this unidirectional id, the distance (with sign) between interface and point
 * coordinate is calculated for all interfaces. The unidirectional piece id is then determined
 * using the location where sign change occurs in this distance.
 *
 * Terms like left and right here will be used synonymously for all directions.
 *
 * There are no bounds for @p p. The interfaces (processed in constructor) alone are used to
 * determine piece id, there are no other bounds for @p p.
 *
 * It is expected that @p p is a cell center. This is because the way this class is designed to
 * operate.
 */
usi PiecewiseFunction::get_piece_id(const Point<dim> &p)
{
    // calculate the distances from interfaces
    std::array<std::vector<double>, dim> distances; // distances from interfaces
    int dir;
    for(dir=0; dir<dim; dir++){
        for(double loc: ilocs_[dir]) distances[dir].emplace_back(p[dir]-loc);
    }

    // determine unidirectional piece ids
    std::array<usi, 3> upid;
    for(dir=0; dir<dim; dir++){
        if(np_[dir] == 1){
            // single piece, no interfaces
            upid[dir] = 0;
        }
        else{
            // at least 2 pieces (1 interface)
            if(distances[dir][0] < 0){
                // lies left to the left most interface
                upid[dir] = 0;
            }
            else if(distances[dir][np_[dir]-2] > 0){
                // lies right to the right most interface
                upid[dir] = np_[dir] - 1;
            }
            else{
                // lies between left most and right most interfaces
                // loop 1: check if the point lies on any interface
                for(int iid=0; iid<np_[dir]-1; iid++){
                    AssertThrow(
                        distances[dir][iid] != 0,
                        StandardExceptions::ExcMessage(
                            "The provided point lies exactly on an interface."
                        )
                    );
                } // loop over all interfaces
                // loop 2: loop over interfaces to get unidirectional piece id
                for(int iid=1; iid<np_[dir]-1; iid++){
                    if(distances[dir][iid-1] > 0 && distances[dir][iid] < 0){
                        upid[dir] = iid;
                        break;
                    }
                } // loop over all except left most interface
            }
        }
    } // loop over dirs
    return upid[0] + upid[1]*np_[0] + upid[2]*np_[0]*np_[1];
}



/**
 * Constructor. Calls the base constructor and parses all the functions in the file @p filename.
 */
PiecewiseFunction::PiecewiseFunction(
    const DoFHandler<dim> &dh,
    const std::map<unsigned int, Point<dim>> &dl,
    std::array<LA::MPI::Vector, 5> &gcv,
    const std::string &filename,
    const NavierStokes *ns_ptr
)
: IC(dh, dl, gcv), ns_ptr_(ns_ptr)
{
    std::ifstream fn_file(filename);
    AssertThrow(
        fn_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file provided for reading IC functions."
        )
    );

    // Note the choice: primitive or conservative functions
    std::string line;
    try{
        std::getline(fn_file, line);
    } catch(...){
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "Unable to read first line of the IC file."
            )
        );
    }
    if(line == "c"){
        prim_fns_ = false;
    } else if(line == "p"){
        prim_fns_ = true;
    } else{
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "The first line in the file containing functions must be 'p' or 'c' to indicate "
                "whether the functions are for primitive or conservative variables."
            )
        );
    }

    // get the number of pieces
    for(int dir=0; dir<dim; dir++){
        try{
            std::getline(fn_file, line);
        } catch(...){
            AssertThrow(
                false,
                StandardExceptions::ExcMessage(
                    "Unable to read number of pieces in each direction from IC file."
                )
            );
        }
        np_[dir] = stod(line); // assuming line contains a number
        AssertThrow(
            np_[dir] >= 1,
            StandardExceptions::ExcMessage(
                "The number of pieces provided in each direction must be >= 1."
            )
        );
    }

    // get the interface locations
    std::vector<std::string> splits;
    for(int dir=0; dir<dim; dir++){
        try{
            std::getline(fn_file, line);
        } catch(...){
            AssertThrow(
                false,
                StandardExceptions::ExcMessage(
                    "Unable to read the interface locations in each direction from IC file."
                )
            );
        }
        utilities::split_string(line, " ", splits); // split at spaces
        AssertThrow(
            splits.size() == np_[dir]-1,
            StandardExceptions::ExcMessage(
                "The number of interface locations must be exactly one less than the number of "
                "pieces provided earlier for each direction."
            )
        );
        for(const auto &s: splits) ilocs_[dir].emplace_back(stod(s)); // assuming s is a number
    }

    // resize fpps_. This automatically initialises the pointed objects using default ctor of
    // FunctionParser which is a scalar function ctor
    fpps_.resize(np_[0]*np_[1]*np_[2]);

    // parse the functions
    std::string variables("x,y,z");
    std::map<std::string, double> constants; // empty

    for(usi pid=0; pid<fpps_.size(); pid++){
        for(cvar var: cvar_list){
            try{
                std::getline(fn_file, line);
            } catch(...){
                AssertThrow(
                    false,
                    StandardExceptions::ExcMessage(
                        "Unable to read functions from IC file. Probably insufficient number of "
                        "functions were provided in the file."
                    )
                );
            }
            fpps_[pid][var].reset(new FunctionParser<dim>());
            fpps_[pid][var]->initialize(variables, line, constants);
        }
    }
}



/**
 * Sets the IC. Algo:
 * - Loop over owned cells
 *   - Determine piece id for all dofs of that cell using cell center
 *   - If prim_fns_ is false
 *     - Directly evaluate the functions and set them in g_cvars
 *   - Else
 *     - Use NavierStokes::prim_to_cons to get the conservative state
 */
void PiecewiseFunction::set()
{
    Point<dim> center; // cell center
    usi pid; // piece id
    std::vector<unsigned int> dof_ids(dof_handler.get_fe().dofs_per_cell); // dof ids of cell
    State cons;
    double rho, p;
    Tensor<1,dim> vel;
    for(const auto &cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        center = cell->center();
        pid = get_piece_id(center);
        cell->get_dof_indices(dof_ids);

        for(psize i: dof_ids){
            if(!prim_fns_){
                // functions parsed in fpps_ are conservative functions
                for(cvar var: cvar_list){
                    g_cvars[var][i] = fpps_[pid][var]->value(dof_locations[i]);
                }
            }
            else{
                // functions parsed in fpps_ are primitive functions
                rho = fpps_[pid][0]->value(dof_locations[i]);
                for(int dir=0; dir<dim; dir++){
                    vel[dir] = fpps_[pid][dir+1]->value(dof_locations[i]);
                }
                p = fpps_[pid][4]->value(dof_locations[i]);
                ns_ptr_->prim_to_cons(rho, vel, p, cons);
                for(cvar var: cvar_list) g_cvars[var][i] = cons[var];
            }
        } // loop over cell (global) dofs
    } // loop over owned cells
    assert_positivity(); // calls IC::assert_positivity()
}



#ifdef DEBUG
void PiecewiseFunction::test()
{
    utilities::Testing t("PiecewiseFunction", "class");
    utilities::ICTestData ictd(5,2); // divisions, degree
    NavierStokes ns("air");

    // run these tests in serial in a folder where there is an IC file named ic_fns.txt
    {
        t.new_block("testing construction");
        std::unique_ptr<IC> icp = std::make_unique<PiecewiseFunction>(
            ictd.dof_handler,
            ictd.dof_locations,
            ictd.g_cvars,
            "ic_fns.txt",
            &ns
        );
    }
}
#endif
