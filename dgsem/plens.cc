/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * Constructor. Calls declare_parameters() and parses parameter file. The parameter file should be
 * named 'input.prm'. Asserts `mhod>0`. Sets ns_ptr and mapping_ptr to nullptr.
 */
PLENS::PLENS(
    const usi mhod,
    const usi fe_degree
)
:
mpi_comm(MPI_COMM_WORLD),
pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm)==0)),
triang(mpi_comm),
ns_ptr(nullptr),
mapping_ho_degree(mhod),
mapping_ptr(nullptr),
fe(fe_degree),
fe_face(fe_degree),
dof_handler(triang)
{
    declare_parameters();
    prm.parse_input("input.prm");

    AssertThrow(
        mhod > 0,
        StandardExceptions::ExcMessage("High order mapping degree must be >= 1.")
    );
}



/**
 * Destructor.
 */
PLENS::~PLENS()
{
    if(bid_list.size() > 0){
        for(auto cur_bc_pair: bc_list) delete cur_bc_pair.second;
    }
}



/**
 * Declares all parameters.
 *
 * A total of 12 BCs with ids 0 to 11 are declared in subsection BCs. If any of them are left unset
 * in 'input.prm', then the default type "none" is assumed. Data required for BCs is specified in
 * subsections BCs.bid<x> where x takes values 0 to 11. In the function set_BC(), the value of x
 * is set to all physical boundary ids and corresponding sections are read from prm file. Thus, any
 * physical boundary id cannot have "none" as its type.
 */
void PLENS::declare_parameters()
{
    pcout << "\nDeclaring parameters\n\n";

    prm.enter_subsection("mesh");
    {
        prm.declare_entry(
            "type",
            "straight",
            Patterns::Selection("curved|straight"),
            "Type of mesh. Options: curved|straight. If the mesh is curved, then it has to either "
            "be cylinder flare geometry or blunted double cone geometry. See the note "
            "'pens2D to plens' and entries around WJ-05-Apr-2021."
        );

        prm.declare_entry(
            "format",
            "msh",
            Patterns::Selection("msh|vtk"),
            "Format of mesh file. Only those that permit setting boundary id are of use. "
            "Options: msh|vtk"
        );

        prm.declare_entry(
            "file name",
            "",
            Patterns::FileName(),
            "Mesh file name. No checks are done on the format and existance of the file"
        );

        prm.declare_entry(
            "curved subtype",
            "cylinder flare",
            Patterns::Selection("cylinder flare|blunted double cone"),
            "If 'type' was set to 'curved', then this entry indicates which type of curved mesh "
            "is in question. Accordingly, the entries of that subsection are considered. Options "
            "'cylinder flare|blunted double cone'"
        );

        // Entries for curved mesh. See entries around WJ-05-Apr-2021 and the note
        // 'pens2D to plens'
        prm.enter_subsection("cylinder flare");
        {
            prm.declare_entry(
                "axis direction",
                "1 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a vector indicating axis direction. The vector may "
                "be scaled arbitrarily"
            );

            prm.declare_entry(
                "axis point",
                "0 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a coordinate indicating any point on the axis"
            );
        } // subsection cylinder flare
        prm.leave_subsection();

        prm.enter_subsection("blunted double cone");
        {
            prm.declare_entry(
                "axis direction",
                "1 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a vector indicating axis direction. The vector may "
                "be scaled arbitrarily. Here, it is assumed that axis points away from the nose"
            );

            prm.declare_entry(
                "separation point",
                "0 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a coordinate indicating the bifurcation between nose "
                "section and cone section. This point must also lie on the axis"
            );

            prm.declare_entry(
                "nose center",
                "0 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a coordinate indicating the nose center. This point, "
                "naturally, lies on the axis."
            );
        } // subsection blunted double cone
        prm.leave_subsection();
    }
    prm.leave_subsection(); // subsection mesh

    prm.enter_subsection("Navier-Stokes");
    {
        prm.declare_entry(
            "gas name",
            "air",
            Patterns::Selection("air|nitrogen|custom"),
            "Gas to be used in simulation. Options: 'air|nitrogen|custom'. For a 'custom' gas "
            "gamma, molecular weight, Prandtl number, mu0, T0 and S must be provided. For air and "
            "nitrogen, these entries are ignored."
        );

        prm.declare_entry(
            "inviscid",
            "true",
            Patterns::Bool(),
            "If gas name is either 'air' or 'nitrogen', then this value is passed to NavierStokes "
            "constructor. For a 'custom' gas, this setting has no effect and inviscid nature of "
            "such a gas can be set by assigning 0 value to mu0"
        );

        prm.declare_entry(
            "gamma",
            "1.4",
            Patterns::Double(1,3),
            "Specific heat ratio. Has to lie in [1,3]"
        );
        prm.declare_entry(
            "molecular weight",
            "0.028",
            Patterns::Double(1e-3),
            "Molecular weight in kg/mol. Has to lie in [1e-3, infty)"
        );
        prm.declare_entry(
            "Prandtl number",
            "0.71",
            Patterns::Double(1e-8),
            "Has to lie in [1e-8, infty)"
        );
        prm.declare_entry(
            "mu0",
            "1e-5",
            Patterns::Double(0),
            "Value of mu0 in SI units to be used in Sutherland's viscosity model. Has to lie in "
            "[0, infty). Zero value implies inviscid and non-conducting gas"
        );
        prm.declare_entry(
            "T0",
            "300",
            Patterns::Double(1e-4),
            "Value of T0 in K to be used in Sutherland's viscosity model. Has to lie in "
            "[1e-4, infty)"
        );
        prm.declare_entry(
            "S",
            "100",
            Patterns::Double(1e-4),
            "Value of S in K to be used in Sutherland's viscosity model. Has to lie in "
            "[1e-4, infty)"
        );

        prm.declare_entry(
            "auxiliary surface flux scheme",
            "BR1",
            Patterns::Selection("BR1"),
            "Options: 'BR1'"
        );
        prm.declare_entry(
            "auxiliary volume flux scheme",
            "BR1",
            Patterns::Selection("BR1"),
            "Options: 'BR1'"
        );
        prm.declare_entry(
            "inviscid surface flux scheme",
            "HLLC",
            Patterns::Selection("HLLC|Rusanov"),
            "Options: 'HLLC|Rusanov'"
        );
        prm.declare_entry(
            "inviscid volume flux scheme",
            "Chandrashekhar",
            Patterns::Selection("Chandrashekhar"),
            "Options: 'Chandrashekhar'"
        );
        prm.declare_entry(
            "diffusive surface flux scheme",
            "BR1",
            Patterns::Selection("BR1"),
            "Options: 'BR1'"
        );
        prm.declare_entry(
            "diffusive volume flux scheme",
            "BR1",
            Patterns::Selection("BR1"),
            "Options: 'BR1'"
        );
    }
    prm.leave_subsection(); // subsection Navier-Stokes

    prm.enter_subsection("BCs");
    {
        std::string base_name("bid"), cur_name;
        // for gmsh, physical ids must be strictly positive
        for(usi i=1; i<=n_bc_max; i++){
            cur_name = base_name + std::to_string(i);
            prm.enter_subsection(cur_name);
            {
                prm.declare_entry(
                    "type",
                    "none",
                    Patterns::Selection(
                        "none|free|outflow|uniform inflow|uniform temp wall|symmetry|periodic"
                    ),
                    "Type of BC. Options: 'none|free|outflow|uniform inflow|uniform temp wall"
                    "|symmetry|periodic'. 'none' type cannot be specified for a physical "
                    "boundary. For a periodic boundary, separate entries for both surfaces "
                    "have to be made"
                );

                prm.declare_entry(
                    "prescribed p",
                    "1e5",
                    Patterns::Double(1e-16),
                    "Presribed pressure. Must lie in [1e-16, infty). Relevant for: outflow, "
                    "uniform inflow"
                );

                prm.declare_entry(
                    "prescribed T",
                    "1",
                    Patterns::Double(1e-16),
                    "Prescribed temperature. Must lie in [1e-16, infty). Relevant for: "
                    "uniform inflow, uniform temp wall"
                );

                prm.declare_entry(
                    "prescribed velocity",
                    "0 0 0",
                    Patterns::List(Patterns::Double(), dim, dim, " "),
                    "Space-separated values for prescribed velocity. Relevant for: "
                    "uniform inflow, uniform temp wall. For uniform temp wall, this describes "
                    "the wall velocity."
                );

                prm.declare_entry(
                    "periodic direction",
                    "0",
                    Patterns::Integer(0,dim-1),
                    "Cartesian direction for periodic BC. See "
                    "https://www.dealii.org/current/doxygen/deal.II/namespaceGridTools.html#ab22eef800535f9e85a1723a6a36fd0f6. "
                    "Relevant for: periodic"
                );

                prm.declare_entry(
                    "periodic orientation",
                    "left",
                    Patterns::Selection("left|right"),
                    "Whether this face is the first (left) or second (right) surface of the "
                    "periodic pair. See "
                    "https://www.dealii.org/current/doxygen/deal.II/namespaceGridTools.html#aee88c4dce5066a41183b5dd70289b9df."
                    "Consequently, the orientation of surface with id 'other surface boundary id' "
                    "will be treated as complementary to this. Basically, left face to right face "
                    "must go along the periodic direction"
                );

                prm.declare_entry(
                    "other surface boundary id",
                    "1",
                    Patterns::Integer(1,n_bc_max),
                    "The boundary id of the other periodic surface. A BC entry for this face has "
                    "to be made separately with appropriate changes in periodic orientation. The "
                    "periodic direction value will remain same"
                );
            }
            prm.leave_subsection(); // bid<x> subsection
        } // loop over boundary ids
    }
    prm.leave_subsection(); // subsection BCs

    prm.enter_subsection("IC");
    {
        prm.declare_entry(
            "type",
            "piecewise function",
            Patterns::Selection("piecewise function"),
            "Options: 'piecewise function'"
        );

        prm.declare_entry(
            "file name",
            "ic.dat",
            Patterns::FileName(),
            "Relevant for: piecewise function. For piecewise function, this file must contain a "
            "list of functions of conservative variables in pieces of domain."
        );
    }
    prm.leave_subsection(); // subsection IC

    std::ofstream sample_file("sample_input_file.prm");
    AssertThrow(
        sample_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file for writing sample input parameters settings."
        )
    );
    prm.print_parameters(sample_file, ParameterHandler::Text);
    sample_file.close();
    pcout << "Sample input file named 'sample_input_file.prm' written\n";
}



/**
 * Reads mesh based on the settings in prm file.
 *
 * Supported mesh formats (msh|vtk) are shown in declare_parameters(). VTK format also supports
 * boundary ids, see
 * https://dealii.org/developer/doxygen/deal.II/group__simplex.html#ga058cd187cea704428ac1118410cd0fb8.
 * However, there seems to be a version mismatch between the vtk meshes that gmsh writes and what
 * dealii expects. See WJ-07-May-2021
 *
 * For straight edged meshes, the procedure is simple. mapping_ptr is set using
 * MappingQGeneric<dim>(1).
 *
 * For curved meshes, only "cylinder flare" and "blunted double cone" geometries are supported. For
 * cylinder flare, the entire mesh is assigned CylindricalManifold. For blunted double cone, the
 * nose part is assigned SphericalManifold and the rest is assigned CylindricalManifold. See note
 * __pens2D to plens__ and entries around WJ-05-Apr-2021. The algorithm for setting manifold on
 * blunted double cone is
 *
 * - Loop over all cells
 *  - If the dot product of line joining separation point to cell center and axis is negative
 *    - Set manifold id for spherical manifold
 *  - Else
 *    - Set manifold id for cylindrical manifold
 *
 * Where separation point is a point on the axis which bifurcates the cone section from sphere
 * section. The plane normal to axis and passing through the separation point is the bifurcator.
 *
 * For curved meshes, mapping_ptr is set using MappingQGeneric<dim>(mapping_ho_degree).
 *
 * @note It is assumed that the cells are also strictly bifurcated: no cells have bifurcator plane
 * passing through them in the middle. This can be ensured by meshing the two regions separately.
 *
 * If there are any periodic BCs to be set, then `triang` object must be modified. This will be
 * done in set_dof_handler().
 */
void PLENS::read_mesh()
{
    pcout << "Reading mesh ... ";
    std::string type, format, filename;

    // read all data
    prm.enter_subsection("mesh");
    {
        type = prm.get("type");
        format = prm.get("format");
        filename = prm.get("file name");
    } // subsection mesh

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triang);
    GridIn<dim>::Format fmt;

    // set format
    if(format == "msh") fmt = GridIn<dim>::Format::msh;
    else fmt = GridIn<dim>::Format::vtk;

    // common code for straight/curved meshes
    std::ifstream file(filename);
    AssertThrow(
        file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open mesh file"
        )
    );
    grid_in.read(file, fmt);
    file.close();

    // straight edge meshes
    if(type == "straight"){
        // additional settings for straight meshes
        mapping_ptr = std::make_unique<MappingQGeneric<dim>>(1);
    }
    else{
        // curved type
        mapping_ptr = std::make_unique<MappingQGeneric<dim>>(mapping_ho_degree);
        std::string curved_type = prm.get("curved subtype");
        if(curved_type == "cylinder flare"){
            std::string temp;
            std::vector<std::string> splits;

            Tensor<1,dim> axis;
            Point<dim> axis_p;
            prm.enter_subsection("cylinder flare");
            {
                // get axis and axis point
                temp = prm.get("axis direction"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) axis[d] = std::stod(splits[d]);

                temp = prm.get("axis point"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) axis_p[d] = std::stod(splits[d]);
            }
            prm.leave_subsection(); // cylinder flare

            CylindricalManifold<dim> manifold(axis, axis_p);
            triang.set_all_manifold_ids(0);
            triang.set_manifold(0, manifold);
        } // if cylinder flare
        else{
            // guaranteed to be blunted double cone
            std::string temp;
            std::vector<std::string> splits;

            Tensor<1,dim> axis;
            Point<dim> separation_p, nose_center;

            prm.enter_subsection("blunted double cone");
            {
                // get axis and axis point
                temp = prm.get("axis direction"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) axis[d] = std::stod(splits[d]);

                temp = prm.get("separation point"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) separation_p[d] = std::stod(splits[d]);

                temp = prm.get("nose center"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) nose_center[d] = std::stod(splits[d]);
            }
            prm.leave_subsection(); // blunted double cone

            CylindricalManifold<dim> cyl_man(axis, separation_p);
            SphericalManifold<dim> sph_man(nose_center);

            double dotp; // dot product
            for(auto cell: triang.active_cell_iterators()){
                if(!(cell->is_locally_owned())) continue;
                dotp = scalar_product(axis, cell->center() - separation_p);
                if(dotp < 0) cell->set_all_manifold_ids(0); // sphere
                else cell->set_all_manifold_ids(1); // cylinder
            } // loop over owned active cells

            triang.set_manifold(0, sph_man);
            triang.set_manifold(1, cyl_man);
        } // if blunted double cone
    }
    prm.leave_subsection(); // subsection mesh
    pcout << "Completed\n";
}



/**
 * Forms the NavierStokes object based on settings in prm file. See declare_parameters() for the
 * settings.
 */
void PLENS::set_NS()
{
    pcout << "Forming NavierStokes object ... ";
    prm.enter_subsection("Navier-Stokes");
    {
        // first read flux scheme settings since they are common
        std::string temp;
        NavierStokes::aux_surf_flux_scheme asfs;
        NavierStokes::aux_vol_flux_scheme avfs;
        NavierStokes::inv_surf_flux_scheme isfs;
        NavierStokes::inv_vol_flux_scheme ivfs;
        NavierStokes::dif_surf_flux_scheme dsfs;
        NavierStokes::dif_vol_flux_scheme dvfs;

        temp = prm.get("auxiliary surface flux scheme");
        if(temp == "BR1") asfs = NavierStokes::aux_surf_flux_scheme::BR1;

        temp = prm.get("auxiliary volume flux scheme");
        if(temp == "BR1") avfs = NavierStokes::aux_vol_flux_scheme::BR1;

        temp = prm.get("inviscid surface flux scheme");
        if(temp == "HLLC") isfs = NavierStokes::inv_surf_flux_scheme::hllc;
        else isfs = NavierStokes::inv_surf_flux_scheme::rusanov;

        temp = prm.get("inviscid volume flux scheme");
        if(temp == "Chandrashekhar") ivfs = NavierStokes::inv_vol_flux_scheme::chandrashekhar;

        temp = prm.get("diffusive surface flux scheme");
        if(temp == "BR1") dsfs = NavierStokes::dif_surf_flux_scheme::BR1;

        temp = prm.get("diffusive volume flux scheme");
        if(temp == "BR1") dvfs = NavierStokes::dif_vol_flux_scheme::BR1;

        std::string gas_name;
        bool inviscid;
        gas_name = prm.get("gas name");
        if(gas_name != "custom"){
            // 'air' or 'nitrogen'
            inviscid = prm.get_bool("inviscid");
            ns_ptr = std::make_unique<NavierStokes>(
                gas_name, inviscid, asfs, avfs, isfs, ivfs, dsfs, dvfs
            );
        }
        else{
            // custom gas, read modelling parameters
            double gma = prm.get_double("gamma");
            double MW = prm.get_double("molecular weight");
            double Pr = prm.get_double("Prandtl number");
            double mu0 = prm.get_double("mu0");
            double T0 = prm.get_double("T0");
            double S = prm.get_double("S");
            ns_ptr = std::make_unique<NavierStokes>(
                gma, MW, Pr, mu0, T0, S, asfs, avfs, isfs, ivfs, dsfs, dvfs
            );
        }
    }
    prm.leave_subsection(); // subsection Navier-Stokes
    pcout << "Completed\n";
}



/**
 * Sets the dof handler object. Here a loop over BC section of prm file is done to check if there
 * are any periodic BCs. If there are, then periodicity is added appropriately. Step-45 has detials
 * about adding periodicity to distributed meshes. See WJ-13-Apr-2021.
 *
 * Since there will be separate entries for each surface/boundary in a periodic BC, the
 * 'other surface boundary id' is stored in a set and any boundary id which lies in this set is
 * ignored. Of course, this assumes the two entries of periodic boundaries are consistent with
 * physical description in terms of their entries. Collecting and setting matched pairs only once
 * is sufficient.
 *
 * If the 'right' id and 'left' id of a pair match (i.e.; both faces have same id), then a more
 * restricted version of collect_periodic_faces() is called which takes only one boundary id.
 *
 * @note It is required for set_BC() that right id and left id be different. But that is not
 * required for this function. So this function doesn't assert that, set_BC() will. But this
 * function will print a warning stating so.
 *
 * @pre read_mesh() has to be called before this.
 */
void PLENS::set_dof_handler()
{
    pcout << "Setting DoFHandler object ...\n";
    prm.enter_subsection("BCs");
    {
        pcout << "Looking for periodic BCs\n";
        std::string base_name("bid"), cur_name, type;
        // boundary ids of periodic boundary which can be ignored
        std::set<usi> bid_periodic_ignore;
        for(usi i=1; i<=n_bc_max; i++){
            if(bid_periodic_ignore.count(i) == 1){
                // this id has already been accounted for while parsing its partner periodic
                // surface, hence skip
                continue;
            }

            cur_name = base_name + std::to_string(i);
            prm.enter_subsection(cur_name);
            {
                type = prm.get("type");
                if(type == "periodic"){
                    const usi periodic_direction = prm.get_integer("periodic direction");
                    const std::string orientation = prm.get("periodic orientation");
                    const usi other_id = prm.get_integer("other surface boundary id");
                    bid_periodic_ignore.emplace(other_id);

                    pcout << "Found bid " << i << " and " << other_id << " with periodic type\n";

                    std::vector<GridTools::PeriodicFacePair<
                        parallel::distributed::Triangulation<dim>::cell_iterator>
                    > matched_pairs;

                    if(i != other_id){
                        usi left_id, right_id;
                        if(orientation == "left"){
                            left_id = i;
                            right_id = other_id;
                        }
                        else{
                            // prm handler guarantees orientation is either left or right
                            left_id = other_id;
                            right_id = i;
                        }
                        GridTools::collect_periodic_faces(
                            triang,
                            left_id, // 'left' boundary id,
                            right_id, // 'right' boundary id
                            periodic_direction, // direction,
                            matched_pairs
                        );
                        pcout << "Formed matched pairs for left id " << left_id << " and "
                            << "right id " << right_id << "\n";
                    }
                    else{
                        pcout << "WARNING\n Your prm file entry for boundary id " << i << " which "
                            << "is of type 'periodic' has 'other surface boundary id' equal to "
                            << "the boundary id. This is ok for setting the dof handler, but will "
                            << "cause an exception to be thrown while setting BCs. "
                            << "See PLENS::set_BC() for detailed info.\n";
                        GridTools::collect_periodic_faces(
                            triang,
                            i, // 'left' and 'right' boundary id,
                            periodic_direction, // direction,
                            matched_pairs
                        );
                        pcout << "Formed matched pairs with id " << i << "\n";
                    }

                    triang.add_periodicity(matched_pairs);
                }
            }
            prm.leave_subsection(); // bid<x>
        }
    }
    prm.leave_subsection(); // subsection BCs

    dof_handler.distribute_dofs(fe);

    DoFTools::map_dofs_to_support_points(*mapping_ptr, dof_handler, dof_locations);
    pcout << "Completed\n";
}



/**
 * Initialises the solution vectors. An important difference from pens2D is that currently relevant
 * dofs are set directly to what dealii's functions return. So all dofs of neighboring and periodic
 * cells are added instead of just those lying on common face.
 *
 * @pre read_mesh() and set_dof_handler() must be called before this
 */
void PLENS::set_sol_vecs()
{
    pcout << "Initialising solution vectors ... ";
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    for(cvar var: cvar_list){
        g_cvars[var].reinit(locally_owned_dofs, mpi_comm);
        gold_cvars[var].reinit(locally_owned_dofs, mpi_comm);
        gcrk_cvars[var].reinit(locally_owned_dofs, mpi_comm);
        gh_gcrk_cvars[var].reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    }

    for(avar var: avar_list){
        gcrk_avars[var].reinit(locally_owned_dofs, mpi_comm);
        gh_gcrk_avars[var].reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    }
    pcout << "Completed\n";
}



/**
 * Sets the IC
 *
 * @pre read_mesh(), set_NS(), set_dof_handler() and set_sol_vecs() must be called before this
 */
void PLENS::set_IC()
{
    pcout << "Setting IC ... ";
    std::unique_ptr<ICs::IC> ic_ptr;
    prm.enter_subsection("IC");
    {
        std::string type = prm.get("type");
        if(type == "piecewise function"){
            std::string filename = prm.get("file name");
            ic_ptr = std::make_unique<ICs::PiecewiseFunction>(
                dof_handler,
                dof_locations,
                g_cvars,
                filename,
                ns_ptr.get()
            );
        }
    }
    prm.leave_subsection(); // subsection IC

    ic_ptr->set();
    pcout << "Completed\n";
}



/**
 * Sets the BC
 *
 * Every boundary with a unique boundary id is assigned a BC object. For periodic boundary
 * conditions, each boundary of the periodic pair is assigned a separate BC object.
 *
 * First, a loop over all faces owned by this process is used to determing the number of boundaries
 * this process is exposed to. The list PLENS::bid_list is thus populated. Then, based on this
 * list, the boundary condition objects are constructed by parsing the prm file. For periodic
 * boundary conditions, it might well be possible that the pairs are owned by different mpi
 * processes and hence separate objects are assigned to them.
 *
 * @note For the reasons described above, it is mandatory that the 'left' and 'right' boundary ids
 * for a periodic pair are different. Otherwise, it is impossible to set BCs::Periodic::fid
 * correctly. Moreover, consider a special case: this mpi process contains only the 'right' part
 * of a periodic boundary. In that case, it is inevitable that we have two entries for a periodic
 * pair.
 *
 * @warning Dynamic memory allocation will be done for BC objects
 */
void PLENS::set_BC()
{
    // get the boundary ids held by this process
    for(const auto &cell: dof_handler.active_cell_iterators()){
        if(!cell->is_locally_owned()) continue;

        for(usi fid=0; fid<n_faces_per_cell; fid++){
            const auto &face = cell->face(fid);
            if(face->at_boundary()){
                bid_list.emplace(face->boundary_id());
            }
        } // loop over faces
    } // loop over owned active cells

    // now construct the bc objects
    std::string subsec_name, type;
    for(auto cur_bid: bid_list){
        pcout << "Setting BC for bid " << cur_bid;

        subsec_name = "bid" + std::to_string(cur_bid);
        prm.enter_subsection("BCs");
        {
            prm.enter_subsection(subsec_name);
            {
                type = prm.get("type");
                AssertThrow(
                    type != "none",
                    StandardExceptions::ExcMessage(
                        "'none' type BC is not for ids specified in the mesh. It is to be used "
                        "only for ids not existing in the mesh, that too optionally. For ids "
                        "mentioned in mesh, a definite boundary type must be specified."
                    )
                );

                pcout << ": type " << type << "\n";

                if(type == "free"){
                    bc_list[cur_bid] = new BCs::Free(dof_handler, gcrk_cvars, gcrk_avars);
                }
                else if(type == "outflow"){
                    const double p = prm.get_double("prescribed p");
                    pcout << "\tPrescribed pressure: " << p << "\n";
                    bc_list[cur_bid] = new BCs::Outflow(
                        dof_handler,
                        gcrk_cvars,
                        gcrk_avars,
                        p,
                        ns_ptr.get()
                    );
                }
                else if(type == "uniform inflow"){
                    const double p = prm.get_double("prescribed p");
                    const double T = prm.get_double("prescribed T");
                    const std::string vel_str = prm.get("prescribed velocity");
                    std::vector<std::string> splits;
                    Tensor<1,dim> vel;
                    utilities::split_string(vel_str, " ", splits);
                    for(int d=0; d<dim; d++) vel[d] = std::stod(splits[d]);
                    pcout << "\t Prescribed p, T and U: " << p << " " << T << " " << vel << "\n";

                    State cons;
                    const double rho = p/(ns_ptr->get_R()*T);
                    ns_ptr->prim_to_cons(rho, vel, p, cons);

                    bc_list[cur_bid] = new BCs::UniformInflow(
                        dof_handler,
                        gcrk_cvars,
                        gcrk_avars,
                        cons
                    );
                }
                else if(type == "symmetry"){
                    bc_list[cur_bid] = new BCs::Symmetry(
                        dof_handler,
                        gcrk_cvars,
                        gcrk_avars,
                        ns_ptr.get()
                    );
                }
                else if(type == "uniform temp wall"){
                    const double T = prm.get_double("prescribed T");
                    const std::string vel_str = prm.get("prescribed velocity");
                    std::vector<std::string> splits;
                    Tensor<1,dim> vel;
                    utilities::split_string(vel_str, " ", splits);
                    for(int d=0; d<dim; d++) vel[d] = std::stod(splits[d]);
                    pcout << "\t Prescribed T and U: " << T << " " << vel << "\n";

                    bc_list[cur_bid] = new BCs::UniformTempWall(
                        dof_handler,
                        gcrk_cvars,
                        gcrk_avars,
                        T,
                        vel,
                        ns_ptr.get()
                    );
                }
                else if(type == "periodic"){
                    const usi other_id = prm.get_integer("other surface boundary id");
                    AssertThrow(
                        cur_bid != other_id,
                        StandardExceptions::ExcMessage(
                            "Two surfaces of a periodic pair of boundaries must have two "
                            "different boundary ids"
                        )
                    );
                    const usi per_dir = prm.get_integer("periodic direction");

                    // face orientation id
                    // foid == 0 ==> this face is left and 'other' face is right
                    // foid == 1 ==> opposite
                    usi foid;
                    usi left_id, right_id; // for forming matched pairs
                    const std::string per_orientation = prm.get("periodic orientation");

                    if(per_orientation == "left"){
                        foid = 0;
                        left_id = cur_bid;
                        right_id = other_id;
                    }
                    else{
                        foid = 1;
                        left_id = other_id;
                        right_id = cur_bid;
                    }

                    pcout << "\t'Left' id: " << left_id << " 'right' id: " << right_id
                        << " direction: " << per_dir << "\n";

                    std::vector<GridTools::PeriodicFacePair<
                        parallel::distributed::Triangulation<dim>::cell_iterator>
                    > matched_pairs;

                    GridTools::collect_periodic_faces(
                        triang,
                        left_id, // 'left' boundary id,
                        right_id, // 'right' boundary id
                        per_dir, // direction,
                        matched_pairs
                    );

                    bc_list[cur_bid] = new BCs::Periodic(
                        dof_handler,
                        gcrk_cvars,
                        gcrk_avars,
                        matched_pairs,
                        foid
                    );
                }
                else{
                    AssertThrow(
                        false,
                        StandardExceptions::ExcMessage("Invalid boundary type")
                    );
                }
            }
            prm.leave_subsection(); // bid<x>
        }
        prm.leave_subsection(); // BCs
    } // loop over bid_list
}



#ifdef DEBUG
void PLENS::test()
{
    utilities::Testing t("PLENS", "class");

    PLENS problem;
    problem.read_mesh();
}
#endif
