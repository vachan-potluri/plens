/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * Constructor. Calls declare_parameters() and parses parameter file. The parameter file should be
 * named 'input.prm'. Asserts `mhod>0`. Sets ns_ptr and mapping_ptr to nullptr. Completely
 * constructs PLENS::blender_calc by calling BlenderCalculator::parse_parameters().
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
mfld_desc_ptr(nullptr),
fe(fe_degree),
fe_face(fe_degree),
dof_handler(triang),
has_periodic_bc(false),
fdi(fe_degree),
w_1d(fe_degree+1),
ref_D_1d(fe_degree+1),
ref_Q_1d(fe_degree+1),
cdi(fe_degree),
blender_calc(fe_degree, gcrk_blender_var),
clk(mpi_comm),
timer(pcout, TimerOutput::never, TimerOutput::wall_times)
{
    declare_parameters();
    prm.parse_input("input.prm");
    blender_calc.parse_parameters(prm); // completes the construction of blender_calc

    AssertThrow(
        mhod > 0,
        StandardExceptions::ExcMessage("High order mapping degree must be >= 1.")
    );

    MPI_Barrier(mpi_comm);
}



/**
 * Destructor. Deletes the boundary condition object pointers held by PLENS::bc_list.
 */
PLENS::~PLENS()
{
    if(bid_list.size() > 0){
        for(auto cur_bc_pair: bc_list) delete cur_bc_pair.second;
    }
}



/**
 * Declares all parameters. Read the sample file generated after any simulation to understand what
 * this function does.
 *
 * A total of 12 BCs with ids 1 to 12 are declared in subsection BCs. If any of them are left unset
 * in 'input.prm', then the default type "none" is assumed. Data required for BCs is specified in
 * subsections BCs.bid<x> where x takes values 1 to 12. In the function set_BC(), the value of x
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
            "cylinder",
            Patterns::Selection("cylinder|nose cylinder"),
            "If 'type' was set to 'curved', then this entry indicates which type of curved mesh "
            "is in question. Accordingly, the entries of that subsection are considered. Options "
            "'cylinder|nose cylinder'. The subsequent settings are used to form a"
            " ManifoldDescription class."
        );

        // Entries for curved mesh. See entries around WJ-05-Apr-2021 and the note
        // 'pens2D to plens'. Also see WJ-30-Jul-2021.
        prm.enter_subsection("cylinder");
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

        prm.enter_subsection("nose cylinder");
        {
            prm.declare_entry(
                "axis direction",
                "1 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a vector indicating axis direction. The vector may "
                "be scaled arbitrarily. Here, it is assumed that axis points in the nose surface "
                "normal direction"
            );

            prm.declare_entry(
                "separation point",
                "0 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a coordinate indicating the bifurcation between nose "
                "section and cylinder section. This point must lie on the axis"
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
                        "none|free|outflow|uniform inflow|uniform temp wall|symmetry|empty|"
                        "periodic"
                    ),
                    "Type of BC. Options: 'none|free|outflow|uniform inflow|uniform temp wall"
                    "|symmetry|empty|periodic'. 'none' type cannot be specified for a physical "
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
            Patterns::Selection("piecewise function|from archive|from archive restart"),
            "Options: 'piecewise function|from archive|from archive restart'"
        );

        prm.declare_entry(
            "file name",
            "ic.dat",
            Patterns::FileName(),
            "Relevant for: 'piecewise function', 'from archive', 'from archive restart'. For "
            "'piecewise function', this file must contain a list of functions of conservative "
            "variables in pieces of domain. For 'from archive' and 'from archive restart', this "
            "file is the archive's filename."
        );

        prm.declare_entry(
            "archive mesh file name",
            "mesh.msh",
            Patterns::FileName(),
            "Relevant for: 'from archive'. The file name of archive's mesh."
        );

        prm.declare_entry(
            "archive fe degree",
            "1",
            Patterns::Integer(1,9),
            "Relevant for: 'from archive'. The fe degree of the solution stored in archive."
        );
    }
    prm.leave_subsection(); // subsection IC

    prm.enter_subsection("blender parameters");
    {
        prm.declare_entry(
            "variable",
            "pxrho",
            Patterns::Selection("pxrho|p|rho"),
            "The variable used to calculate trouble (and hence the blender). Options: "
            "'pxrho|p|rho'. 'pxrho' means the product of pressure and density."
        );
        prm.declare_entry(
            "threshold factor",
            "0.5",
            Patterns::Double(1e-2),
            "The pre-factor for threshold calculation. See WJ-24-May-2021."
        );
        prm.declare_entry(
            "threshold exponent factor",
            "1.8",
            Patterns::Double(1e-2),
            "The pre-factor of exponent for threshold calculation. See WJ-24-May-2021. Range: "
            "[1e-2, infty)."
        );
        prm.declare_entry(
            "sharpness factor",
            "9.21024",
            Patterns::Double(1e-2),
            "The sharpness factor for calculating blender value from trouble and threshold. See "
            "WJ-24-May-2021. Range: [1e-2, infty)."
        );
        prm.declare_entry(
            "blender min value",
            "1e-3",
            Patterns::Double(1e-8,1),
            "Minimum value of blender (alpha min). See WJ-24-May-2021. Range: [1e-8, 1]."
        );
        prm.declare_entry(
            "blender max value",
            "0.5",
            Patterns::Double(0.5,1),
            "Maximum value of blender (alpha max). See WJ-24-May-2021. Range: [0.5, 1]."
        );
    }
    prm.leave_subsection(); // blender parameters

    prm.enter_subsection("time integration");
    {
        prm.declare_entry(
            "RK order",
            "4",
            Patterns::Integer(3,4),
            "Order of RK time integration. Currently only 3rd and 4th order are supported."
        );
        prm.declare_entry(
            "Courant number",
            "0.1",
            Patterns::Double(1e-8,1),
            "Courant number. Range: [1e-8, 1]."
        );
        prm.declare_entry(
            "start time",
            "0",
            Patterns::Double(),
            "Simulation start time."
        );
        prm.declare_entry(
            "starting output counter",
            "0",
            Patterns::Integer(0),
            "The value of output counter at the beginning of simulation. Every time data output "
            "is done, this counter will be incremented by 1."
        );
        prm.declare_entry(
            "end time",
            "1",
            Patterns::Double(),
            "Simulation end time. If 'start time' is non-zero, then the end time provided here is "
            "treated as the absolute end time, and not relative to start time."
        );
    }
    prm.leave_subsection(); // time integration

    prm.enter_subsection("data output");
    {
        prm.declare_entry(
            "directory",
            "result",
            Patterns::DirectoryName(),
            "The directory into which all results will be put."
        );
        prm.declare_entry(
            "base file name",
            "output",
            Patterns::Anything(),
            "The base name to be used for output files."
        );
        prm.declare_entry(
            "write frequency",
            "100",
            Patterns::Integer(1),
            "The solution will be written after these many time steps."
        );
        prm.declare_entry(
            "calculate steady state error",
            "true",
            Patterns::Bool(),
            "Whether or not steady state error is to be calculated. rhoE is used for this purpose "
            "and if true, cell wise error vector will be added to data output and the error will "
            "be written to the file <base file name>.ss_error."
        );
    }
    prm.leave_subsection(); // data output

    prm.print_parameters("sample_input_file.prm", ParameterHandler::KeepDeclarationOrder);
    pcout << "Sample input file named 'sample_input_file.prm' written\n";

    MPI_Barrier(mpi_comm);
} // declare_parameters



/**
 * Reads mesh based on the settings in prm file.
 *
 * Supported mesh formats (msh|vtk) are shown in declare_parameters(). VTK format also supports
 * boundary ids, see
 * https://dealii.org/developer/doxygen/deal.II/group__simplex.html#ga058cd187cea704428ac1118410cd0fb8.
 * However, there seems to be a version mismatch between the vtk meshes that gmsh writes and what
 * dealii expects. See WJ-07-May-2021. As a result, the vtk format is never really used.
 *
 * For straight edged meshes, the procedure is simple. PLENS::mapping_ptr is set using
 * MappingQGeneric<dim>(1).
 *
 * For curved meshes, the classes in ManifoldDescriptions are used to set the curvatures in the
 * mesh. The parameter options in prm file are set according to the classes in this namespace. All
 * classes in this namespace provide a `set()` function which is used to assign manifold(s) to
 * PLENS::triang. PLENS::mapping_ptr is set using MappingQGeneric<dim>(mapping_ho_degree).
 *
 * If there are any periodic BCs to be set, then `triang` object must be modified so that the cells
 * connected through periodicity are marked as relevant. This will be done in set_dof_handler().
 *
 * When using ICs::FromArchive as the IC, PLENS::mfld_desc_ptr is passed for its constructor.
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
        if(curved_type == "cylinder"){
            std::string temp;
            std::vector<std::string> splits;

            Tensor<1,dim> axis;
            Point<dim> axis_p;
            prm.enter_subsection("cylinder");
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

            mfld_desc_ptr = std::make_unique<ManifoldDescriptions::Cylinder>(axis, axis_p);
            // SetManifold::cylinder_flare(axis, axis_p, triang);
        } // if cylinder flare
        else{
            // guaranteed to be nose cylinder
            std::string temp;
            std::vector<std::string> splits;

            Tensor<1,dim> axis;
            Point<dim> separation_p, nose_center;

            prm.enter_subsection("nose cylinder");
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

            mfld_desc_ptr = std::make_unique<ManifoldDescriptions::NoseCylinder>(
                axis,
                separation_p,
                nose_center
            );
            // SetManifold::blunted_double_cone(axis, separation_p, nose_center, triang);
        } // if blunted double cone

        // apply the manifold for curved mesh
        mfld_desc_ptr->set(triang);
    } // curved subtype
    prm.leave_subsection(); // subsection mesh
    pcout << "Completed\n";

    MPI_Barrier(mpi_comm);
} // read_mesh



/**
 * Forms the NavierStokes object PLENS::ns_ptr based on settings in prm file. See
 * declare_parameters() for the settings.
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

    MPI_Barrier(mpi_comm);
}



/**
 * Sets the dof handler object. Here a loop over BC section of prm file is done to check if there
 * are any periodic BCs. If there are, then periodicity is added appropriately. Step-45 has detials
 * about adding periodicity to distributed meshes. See WJ-13-Apr-2021.
 *
 * Since there will be separate entries for each surface/boundary in a periodic BC, the
 * 'other surface boundary id' is stored in a set and any boundary id which lies in this set is
 * ignored. Of course, this assumes the two entries of periodic boundaries are consistent with
 * physical description in terms of their entries for 'left' and 'right' ids. Collecting and
 * setting matched pairs only once is sufficient.
 *
 * @remark It was realised later the call to `Triangulation::add_periodicity()` must happen only
 * once. This means __all__ periodic boundary pairs must be collected in a single vector of
 * PerioficFacePair. To enable this, the variable PLENS::has_periodic_bc was introduced.
 *
 * If the 'right' id and 'left' id of a pair match (i.e.; both faces have same id), then a more
 * restricted version of collect_periodic_faces() is called which takes only one boundary id.
 *
 * @note It is required for set_BC() that right id and left id be different. But that is not
 * required for this function. So this function doesn't assert that, set_BC() will. But this
 * function will print a warning stating so.
 *
 * This function populates PLENS::dof_locations and calls form_neighbor_face_matchings(),
 * calc_metric_terms(). If PLENS::has_periodic_bc is `true`, then form_neighbor_face_matchings()
 * doesn't populate the PLENS::locally_relevant_dofs and instead, those will be obtained using
 * `DoFTools::extract_locally_relevant_dofs()` in set_sol_vecs().
 *
 * @pre read_mesh() has to be called before this.
 */
void PLENS::set_dof_handler()
{
    pcout << "Setting DoFHandler object ...\n";

    // the matched pairs will be appended to this list
    std::vector<GridTools::PeriodicFacePair<
        parallel::distributed::Triangulation<dim>::cell_iterator>
    > matched_pairs;

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
                    has_periodic_bc = true;
                    const usi periodic_direction = prm.get_integer("periodic direction");
                    const std::string orientation = prm.get("periodic orientation");
                    const usi other_id = prm.get_integer("other surface boundary id");
                    bid_periodic_ignore.emplace(other_id);

                    pcout << "Found bid " << i << " and " << other_id << " with periodic type\n";

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
                } // if type is periodic
            }
            prm.leave_subsection(); // bid<x>
        } // loop over BCs
    }
    prm.leave_subsection(); // subsection BCs

    // add all periodic face relations at once
    if(has_periodic_bc) triang.add_periodicity(matched_pairs);

    dof_handler.distribute_dofs(fe);

    DoFTools::map_dofs_to_support_points(*mapping_ptr, dof_handler, dof_locations);
    pcout << "Completed\n";

    MPI_Barrier(mpi_comm);

    form_neighbor_face_matchings(locally_relevant_dofs);
    calc_metric_terms();

    MPI_Barrier(mpi_comm);
} // set_dof_handler



/**
 * Initialises the solution vectors. If PLENS::has_periodic_bc is `false`, the locally relevant
 * dofs are already set when form_neighbor_face_matchings() was called from set_dof_handler().
 * Otherwise, they are set here using `DoFTools::extract_locally_relevant_dofs()`.
 *
 * For initialising PLENS::gcrk_alpha and PLENS::gh_gcrk_alpha, the functions of
 * `Utilities::MPI::Partitioner` are used. See the documentation of these variables and also
 * @ref cell_indices. The partitioner for this purpose is a cell index partitioner returned by
 * `Triangulation::global_active_cell_index_partitioner()`.
 *
 * @pre read_mesh() and set_dof_handler() must be called before this
 */
void PLENS::set_sol_vecs()
{
    pcout << "Initialising solution vectors ... ";
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    // if the problem has periodic BC, then face dof matching doesn't give relevant dofs
    if(has_periodic_bc){
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    }

    for(cvar var: cvar_list){
        g_cvars[var].reinit(locally_owned_dofs, mpi_comm);
        gcrk_cvars[var].reinit(locally_owned_dofs, mpi_comm);
        gcrk_rhs[var].reinit(locally_owned_dofs, mpi_comm);
        gprk_rhs[var].reinit(locally_owned_dofs, mpi_comm);
        gpprk_rhs[var].reinit(locally_owned_dofs, mpi_comm);
        gh_gcrk_cvars[var].reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    }

    for(avar var: avar_list){
        gcrk_avars[var].reinit(locally_owned_dofs, mpi_comm);
        gh_gcrk_avars[var].reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    }

    for(usi d=0; d<dim; d++){
        gcrk_vel[d].reinit(locally_owned_dofs, mpi_comm);
    }

    gcrk_mu.reinit(locally_owned_dofs, mpi_comm);
    gcrk_k.reinit(locally_owned_dofs, mpi_comm);
    gcrk_p.reinit(locally_owned_dofs, mpi_comm);
    gcrk_T.reinit(locally_owned_dofs, mpi_comm);
    gcrk_blender_var.reinit(locally_owned_dofs, mpi_comm);
    rhoE_old.reinit(locally_owned_dofs, mpi_comm);

    // the return type is a weak ptr, it must be converted to shared ptr for usage
    // see https://en.cppreference.com/w/cpp/memory/weak_ptr
    const std::shared_ptr<const Utilities::MPI::Partitioner> cell_partitioner =
        triang.global_active_cell_index_partitioner().lock();
    
    gcrk_alpha.reinit(cell_partitioner->locally_owned_range(), mpi_comm);
    gh_gcrk_alpha.reinit(
        cell_partitioner->locally_owned_range(),
        cell_partitioner->ghost_indices(),
        mpi_comm
    );

    pcout << "Completed\n";

    MPI_Barrier(mpi_comm);
}



/**
 * Sets the IC. If the IC type is ICs::FromArchive, then PLENS::mfld_desc_ptr and
 * PLENS::mapping_ptr are applied to archived triangulation. Generally this is always the case.
 * This may be changed in future according to requirement.
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
        else if(type == "from archive"){
            const std::string ar_mesh_filename = prm.get("archive mesh file name");
            const std::string ar_filename = prm.get("file name");
            const usi ar_fe_degree = prm.get_integer("archive fe degree");
            ic_ptr = std::make_unique<ICs::FromArchive>(
                dof_handler,
                dof_locations,
                g_cvars,
                mpi_comm,
                ar_mesh_filename,
                mfld_desc_ptr,
                mapping_ptr,
                ar_fe_degree,
                ar_filename
            );
        }
        else{
            // from archive restart
            const std::string ar_filename = prm.get("file name");
            ic_ptr = std::make_unique<ICs::FromArchiveRestart>(
                dof_handler,
                dof_locations,
                g_cvars,
                ar_filename
            );
        }
    }
    prm.leave_subsection(); // subsection IC

    ic_ptr->set();
    pcout << "Completed\n";

    MPI_Barrier(mpi_comm);
}



/**
 * Sets the BC
 *
 * Every boundary with a unique boundary id is assigned a BC object. For periodic boundary
 * conditions, each boundary of the periodic pair is assigned a separate BC object.
 *
 * First, a loop over all faces owned by this process is used to determine the number of boundaries
 * this process is exposed to. The list PLENS::bid_list is thus populated. Then, based on this
 * list, the boundary condition objects (PLENS::bc_list) are constructed by parsing the prm file.
 * For periodic boundary conditions, it might well be possible that the pairs are owned by
 * different mpi processes and hence separate objects are assigned to them.
 *
 * A physical boundary (i.e.; one which has an id assigned in mesh file) cannot be of type `none`.
 *
 * @note For the reasons described above, it is mandatory that the 'left' and 'right' boundary ids
 * for a periodic pair are different. Otherwise, it is impossible to set BCs::Periodic::fid
 * correctly. Moreover, consider a special case: this mpi process contains only the 'right' part
 * of a periodic boundary. In that case, it is inevitable that we have two entries for a periodic
 * pair.
 *
 * @warning Dynamic memory allocation will be done for BC objects. They are deleted in the
 * destructor.
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
                else if(type == "empty"){
                    bc_list[cur_bid] = new BCs::Empty(
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
                        gh_gcrk_cvars,
                        gh_gcrk_avars,
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

    pcout << "Completed setting BCs\n";
    MPI_Barrier(mpi_comm);
} // set_BC



/**
 * Reads all the time integration settings and stores them in appropriate variables.
 */
void PLENS::read_time_settings()
{
    prm.enter_subsection("time integration");
    {
        rk_order = prm.get_integer("RK order");
        cur_time = prm.get_double("start time");
        end_time = prm.get_double("end time");
        Co = prm.get_double("Courant number");
        output_counter = prm.get_integer("starting output counter");
    }
    prm.leave_subsection();

    prm.enter_subsection("data output");
    {
        write_freq = prm.get_integer("write frequency");
    }
    prm.leave_subsection();

    n_time_steps = 0;
}



/**
 * Runs the entire simulation. A copy of the prm file along with fe degree and mapping degree is
 * written to the output directory at the beginning. The fe and mapping degrees are added as
 * comments in the copy printed. The copy will be printed at
 * `<output directory>/simulation_parameters.prm`.
 */
void PLENS::run()
{
    std::string op_dir, filename;
    prm.enter_subsection("data output");
    {
        op_dir = prm.get("directory");
    }
    prm.leave_subsection();
    filename = op_dir + "/simulation_parameters.prm";
    pcout << "Writing read parameters to " << filename << "\n";
    prm.print_parameters(filename, ParameterHandler::KeepDeclarationOrder);
    if(Utilities::MPI::this_mpi_process(mpi_comm) == 0){
        std::ofstream file(filename, std::ios::app);
        file << "\n# FE degree " << fe.degree
            << "\n# Mapping degree " << mapping_ptr->get_degree() << "\n";
        file.close();
    }

    pcout << "\n\nStarting simulation\n";
    clk.start();
    while(cur_time < end_time) update();
    write();
    clk.stop();
}



// * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * * //



/**
 * Populates the data in PLENS::nei_face_matching_dofs. See also @ref face_assem and
 * plens_test::face_dof_matching_test(). Algorithm:
 * - Loop over all owned cells
 *   - Loop over all internal faces
 *     - Loop over all dofs on face
 *       - Loop over all neighbor dofs on the same face
 *         - If the dof locations match (upto tolerance level `tol` times minimum vertex distance),
 *           populate the map
 *
 * While doing this, also polulates the `loc_rel_dofs` parameter to store the ghost cell
 * dofs lying on the partition interface in case PLENS::has_periodic_bc is `false`.
 *
 * @pre Must be called after read_mesh() and set_dof_handler(). May be called at the end of
 * set_dof_handler() after PLENS::dof_locations are calculated.
 */
void PLENS::form_neighbor_face_matchings(
    IndexSet& loc_rel_dofs,
    const double tol
)
{
    pcout << "Matching neighbor side dof ids on internal faces ... ";
    std::vector<psize> dof_ids(fe.dofs_per_cell), dof_ids_nei(fe.dofs_per_cell);
    for(const auto &cell: dof_handler.active_cell_iterators()){
        if(!cell->is_locally_owned()) continue;
        cell->get_dof_indices(dof_ids);

        for(usi fid=0; fid<n_faces_per_cell; fid++){
            // set size of nei_face_matching_dofs
            nei_face_matching_dofs[cell->index()][fid].resize(fe_face.dofs_per_face, 0);
            if(cell->face(fid)->at_boundary()) continue;

            const auto &neighbor = cell->neighbor(fid);
            usi fid_nei = cell->neighbor_of_neighbor(fid);
            neighbor->get_dof_indices(dof_ids_nei);

            for(usi i=0; i<fe_face.dofs_per_face; i++){
                Point<dim> loc = dof_locations[dof_ids[fdi.maps[fid].at(i)]];

                for(usi j=0; j<fe_face.dofs_per_face; j++){
                    Point<dim> loc_nei = dof_locations[dof_ids_nei[fdi.maps[fid_nei].at(j)]];
                    Point<dim> diff(loc - loc_nei);
                    if(diff.norm() < tol*cell->minimum_vertex_distance()){
                        // match obtained
                        nei_face_matching_dofs[cell->index()][fid][i] = j;
                        // for some reason, add_index is only working when there are no periodic
                        // BCs
                        if(!has_periodic_bc){
                            loc_rel_dofs.add_index(dof_ids_nei[fdi.maps[fid_nei].at(j)]);
                        }

                        // print if i and j are not equal (abnormal match)
                        if(i != j){
                            std::cout << "\n\tCell id: " << cell->index()
                                << ", Cell-local face id: " << fid
                                << ", Face-local dof id: " << i
                                << "\n\tFound abnormal match with\n"
                                << "\tNeighbor cell id: " << neighbor->index()
                                << ", Neighbor-local face id: " << fid_nei
                                << ", Neighbor-face-local dof id: " << j << "\n";
                        }
                    }
                } // loop over neighbors dofs on same face
            } // loop over face dofs
        } // loop over internal faces
    } // loop over owned cells
    pcout << "Completed\n";
} // form_neighbor_face_matchings



/**
 * Calculates all metric terms. More specifically
 * 1. 1D quadrature weights @f$w_i@f$
 * 2. 1D Differentiation matrix @f$D@f$ (strong version) and also the matrix @f$Q@f$
 * 3. Other metric terms using MetricTerms and populate PLENS::metrics. Once these are obtained,
 *    the weights and matrices can be ignored.
 *
 * @f[
 * D_{ij} = \frac{\partial l_j}{\partial \xi}(\xi_i)\\
 * Q_{ij} = w_i D_{ij}
 * @f]
 *
 * @pre `dof_handler` must be ready to use. This function can be called at the end of
 * PLENS::set_dof_handler().
 */
void PLENS::calc_metric_terms()
{
    pcout << "Calculating all relevant metric terms ... ";
    // Set weights. Size of w_1d set in ctor
    QGaussLobatto<1> quad_lgl_1d(fe.degree+1);
    for(usi i=0; i<=fe.degree; i++){
        w_1d[i] = quad_lgl_1d.weight(i);
    }
    const std::vector<Point<1>> &points_1d = quad_lgl_1d.get_points();

    // Compute the 1D differentiation matrix in reference cell. Size of ref_D_1d set in ctor
    FE_DGQ<1> fe_1d(fe.degree);
    for(usi row=0; row<=fe.degree; row++){
        for(usi col=0; col<=fe.degree; col++){
            ref_D_1d(row,col) = fe_1d.shape_grad(col, points_1d[row])[0];
            ref_Q_1d(row,col) = w_1d[row]*ref_D_1d(row,col);
        }
    }

    // Calculate and store the metric terms
    FEValues<dim> fe_values(
        *mapping_ptr,
        fe,
        QGaussLobatto<dim>(fe.degree+1),
        update_values | update_jacobians | update_quadrature_points
    );
    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!cell->is_locally_owned()) continue;

        fe_values.reinit(cell);
        metrics.emplace(std::make_pair(
            cell->index(),
            MetricTerms<dim>(fe_values, ref_Q_1d)
        ));
    } // loop over owned cells
    pcout << "Completed\n";
} // calc_metric_terms



/**
 * Calculates the surface flux for a given `stage` into the variable `surf_flux_term`. The
 * information about stages is given in NavierStokes. See also @ref face_assem and BCs::BC. The
 * algorithm employed is as in pens2D.
 * - Loop over owned cells
 *   - Loop over faces
 *     - If face is at boundary
 *       - Calculate the surface flux using BC objects from PLENS::bc_list
 *     - Else (internal face)
 *       - If the neighbor is also owned by this process
 *         - Compute the flux using owner and neighbor solution.
 *         - Set the owner flux using outward (w.r.t. owner) normal
 *         - Set the neighbor flux using inward normal (reversed components)
 *       - Else (the neighbor is a ghost cell)
 *         - Use ghosted version of solution vectors
 *         - Compute the flux using outward normal (w.r.t. owner)
 *         - Neighbor flux will be computed by a different process, leave it as is
 *
 * If `stage == 1` or `stage == 2`, then only `gcrk_cvars` will be used to calculate the flux.
 * Auxiliary variables will be passed to the wrappers (since they require these) but they will not
 * be used inside the wrappers. If `stage == 3`, both `gcrk_cvars` and `gcrk_avars` will be used.
 * In all stages, the ghosted version of required vectors are also used.
 *
 * @pre This function assumes that `gh_gcrk_cvars` and `gh_gcrk_avars` are ready to use. Also,
 * assumes that the entire setup (including BCs) is complete.
 *
 * @note This function resizes `surf_flux_term`.
 */
void PLENS::calc_surf_flux(
    const usi stage,
    locly_ord_surf_flux_term_t<double> &surf_flux_term
) const
{
    AssertThrow(
        stage >=1 && stage <= 3,
        StandardExceptions::ExcMessage(
            "'stage' parameter must be 1, 2 or 3. Nothing else."
        )
    );

    // (re)set the size of surf_flux_term
    for(cvar var: cvar_list){
        for(const auto &cell: dof_handler.active_cell_iterators()){
            if(!cell->is_locally_owned()) continue;
            for(usi face_id=0; face_id<n_faces_per_cell; face_id++){
                surf_flux_term[var][cell->index()][face_id].resize(fe_face.dofs_per_face, 0);
            }
        }
    }

    // For internal faces not at processor boundary, flux is calculated only once from the owner
    // side. Owner is defined as the cell having lesser index. The neighbor flux is set using the
    // flux calculated from owner side by variable-wise multiplying with these signs. For stage 1,
    // the flux doesn't have sign reversal, while for stages 2 & 3, it does have a sign reversal
    // See https://stackoverflow.com/questions/12844475/why-cant-simple-initialize-with-braces-2d-stdarray
    const std::array<std::array<float, 5>, 3> reverse_flux_sign = {{
        {1,1,1,1,1},
        {-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1}
    }};

    const usi stage_id = stage-1; // to access wrappers from NavierStokes and BC
    FEFaceValues<dim> fe_face_values(
        *mapping_ptr,
        fe,
        QGaussLobatto<dim-1>(fe.degree+1),
        update_normal_vectors
    ); // to get normal vectors at dof locations on face
    std::vector<psize> dof_ids(fe.dofs_per_cell), dof_ids_nei(fe.dofs_per_cell);

    for(const auto &cell: dof_handler.active_cell_iterators()){
        if(!cell->is_locally_owned()) continue;

        cell->get_dof_indices(dof_ids);
        for(usi face_id=0; face_id<n_faces_per_cell; face_id++){
            const auto &face = cell->face(face_id);
            fe_face_values.reinit(cell, face_id);

            if(face->at_boundary()){
                // boundary face, use BC objects for flux
                for(usi face_dof=0; face_dof<fe_face.dofs_per_face; face_dof++){
                    FaceLocalDoFData ldd(cell->index(), face_id, face_dof);
                    usi bid = face->boundary_id();

                    // set inner cons and avars
                    State cons, cons_gh;
                    Avars av, av_gh;
                    psize global_dof_id = dof_ids[fdi.maps[face_id].at(face_dof)];
                    for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][global_dof_id];
                    for(avar var: avar_list) av[var] = gcrk_avars[var][global_dof_id];

                    // first get ghost state
                    CAvars cav(&cons, &av), cav_gh(&cons_gh, &av_gh);
                    Tensor<1,dim> normal = fe_face_values.normal_vector(face_dof);
                    bc_list.at(bid)->get_ghost_wrappers[stage_id](ldd, cav, normal,cav_gh);

                    // now get the flux
                    State flux;
                    ns_ptr->surf_flux_wrappers[stage_id](cav, cav_gh, normal, flux);

                    // set surf_flux_term object
                    for(cvar var: cvar_list){
                        surf_flux_term[var][cell->index()][face_id][face_dof] = flux[var];
                    }
                } // loop over face dofs
            } // boundary face

            else if(cell->neighbor(face_id)->is_ghost()){
                // internal face at processor boundary
                const auto &neighbor = cell->neighbor(face_id);
                usi face_id_nei = cell->neighbor_of_neighbor(face_id);
                neighbor->get_dof_indices(dof_ids_nei);

                for(usi face_dof=0; face_dof<fe_face.dofs_per_face; face_dof++){
                    // first get neighbor-side matching dof's global id
                    psize gdof_id = dof_ids[fdi.maps[face_id].at(face_dof)];
                    usi face_dof_nei = nei_face_matching_dofs.at(cell->index())[face_id][face_dof];
                    psize gdof_id_nei = dof_ids_nei[fdi.maps[face_id_nei].at(face_dof_nei)];

                    // use ghosted vectors to get neighbor state information
                    State cons, cons_nei;
                    Avars av, av_nei;
                    for(cvar var: cvar_list){
                        cons[var] = gcrk_cvars[var][gdof_id];
                        cons_nei[var] = gh_gcrk_cvars[var][gdof_id_nei];
                    }
                    for(avar var: avar_list){
                        av[var] = gcrk_avars[var][gdof_id];
                        av_nei[var] = gh_gcrk_avars[var][gdof_id_nei];
                    }

                    // get the flux
                    CAvars cav(&cons, &av), cav_nei(&cons_nei, &av_nei);
                    State flux;
                    Tensor<1,dim> normal = fe_face_values.normal_vector(face_dof);
                    ns_ptr->surf_flux_wrappers[stage_id](cav, cav_nei, normal, flux);

                    // set surf_flux_term entries for owner, neighbor's flux will be calculated
                    // by its own process
                    for(cvar var: cvar_list){
                        surf_flux_term[var][cell->index()][face_id][face_dof] = flux[var];
                    }
                } // loop over face dofs
            } // internal face at processor boundary

            else if(cell->index() < cell->neighbor_index(face_id)){
                // processor internal face
                const auto &neighbor = cell->neighbor(face_id);
                usi face_id_nei = cell->neighbor_of_neighbor(face_id);
                neighbor->get_dof_indices(dof_ids_nei);

                for(usi face_dof=0; face_dof<fe_face.dofs_per_face; face_dof++){
                    // first get neighbor-side matching dof's global id
                    psize gdof_id = dof_ids[fdi.maps[face_id].at(face_dof)];
                    usi face_dof_nei = nei_face_matching_dofs.at(cell->index())[face_id][face_dof];
                    psize gdof_id_nei = dof_ids_nei[fdi.maps[face_id_nei].at(face_dof_nei)];

                    // get neighbor state information
                    State cons, cons_nei;
                    Avars av, av_nei;
                    for(cvar var: cvar_list){
                        cons[var] = gcrk_cvars[var][gdof_id];
                        cons_nei[var] = gcrk_cvars[var][gdof_id_nei];
                    }
                    for(avar var: avar_list){
                        av[var] = gcrk_avars[var][gdof_id];
                        av_nei[var] = gcrk_avars[var][gdof_id_nei];
                    }

                    // get the flux
                    CAvars cav(&cons, &av), cav_nei(&cons_nei, &av_nei);
                    State flux;
                    Tensor<1,dim> normal = fe_face_values.normal_vector(face_dof);
                    ns_ptr->surf_flux_wrappers[stage_id](cav, cav_nei, normal, flux);

                    // set surf_flux_term entries for owner and neighbor
                    // reverse the flux for neighbor
                    for(cvar var: cvar_list){
                        surf_flux_term[var][cell->index()][face_id][face_dof] = flux[var];
                        surf_flux_term[var][neighbor->index()][face_id_nei][face_dof_nei]
                            = reverse_flux_sign[stage_id][var]*flux[var];
                    }
                } // loop over face dofs
            } // processor internal face

            else continue;
        } // loop over faces
    } // loop over owned cells
} // calc_surf_flux



/**
 * Calculates conservative variable gradients in a `cell`. The relevant formula is eq. (B.14) of
 * [1]. The volumetric terms are calculated using PLENS::gcrk_cvars and the surface flux is taken
 * from `s1_surf_flux` which is assumed to hold the conservative variable flux for stage 1 (viz.
 * the auxiliary flux).
 *
 * Suppose @f$\partial \rho/\partial x@f$ is to be calculated. Then, a differential equation is
 * constructed:
 * @f[
 * X + \frac{\partial}{\partial x}(-\rho) = 0
 * @f]
 * where @f$X=\partial \rho/\partial x@f$ is to be calculated. This then is in the form shown in
 * eq. (B.1) of [1]. Hence, eq. (B.14) is used to calculate the gradient. For every conservative
 * variable, 3 such equations corresponding to 3 gradients are used.
 *
 * The strategy for calculation is simple. For the two point volume fluxes, use
 * NavierStokes::get_aux_vol_flux(). For the special case of BR1 flux (which is the only option
 * available currently in the NavierStokes class), the direction passed in this function doesn't
 * matter. Also, unlike for stages 2 & 3 (viz. inviscid and diffusive fluxes), the NavierStokes
 * function for auxiliary flux doesn't give a "component", but rather, a value. Hence, this flux
 * must be multiplied by the component of contravariant vector in the gradient direction. See eq.
 * (B.15) of [1]. See also TW1 notes or WJ dated 20-May-2021.
 *
 * Note that for the surface contribution of those dofs lying on face, since this involves stage 1
 * flux, and hence has no "directional component" (as said above), the flux will have same sign in
 * the storage (i.e.; in `s1_surf_flux`) of both cells having the common interface. This is unlike
 * stages 2 and 3 where the flux for every boundary dof is stored (in, say `s2_surf_flux` or
 * `s3_surf_flux`) using normal flux with outward pointing face normal. This means that the
 * surface contribution expression in (B.14) of [1] can be directly used here (without worrying
 * about surface normals).
 *
 * Note that the mapping to reference cell has already been done and we are concerned only with
 * the reference cell-transformed equation and its derivative equations now. (B.14) is one such
 * equation. We need not worry about the physical cell at this point (of course, assuming the
 * metrics are properly calculated).
 *
 * As mentioned in @ref vol_contrib, a unified function for calculating residual is not possible.
 * Moreover, having this function separate also helps if in future, just the algorithm for
 * calculating auxiliary variables needs to be changed. The way I am currently doing is not what
 * is suggested in the DGSEM algorithm. Although refs [1-2] don't mention this explicitly, ref [3]
 * provides a different algo for second derivative calculation in section 2.3.3 with the following
 * comment:
 * _It is trivial to show that two applications of the first derivative operator satisfy the SBP
 * condition. In practice, this is not advisable, as the approximation using two first derivative
 * operations requires a much wider stencil (is less efficient), is less accurate, and leads to
 * only neutrally stable approximations_
 *
 * @param[in] cell The iterator corresponding to the cell in which gradients are to be calculated
 * @param[in] s1_surf_flux Stage 1 surface flux of conservative variables. Access:
 * `s1_surf_flux[var][cell id][cell-local face id][face-local dof id]`
 * @param[out] cons_grad Conservative variable gradients calculated for the cell. Access:
 * `cons_grad[cell-local dof id][dir][var]`
 *
 * @pre `cons_grad` must have the appropriate size. An assertion will be done.
 */
void PLENS::calc_cell_cons_grad(
    const DoFHandler<dim>::active_cell_iterator& cell,
    const locly_ord_surf_flux_term_t<double>& s1_surf_flux,
    std::vector<std::array<State, 3>>& cons_grad
) const
{
    AssertThrow(
        cons_grad.size() == fe.dofs_per_cell,
        StandardExceptions::ExcMessage(
            "The vector provided to store conservative gradients must have the proper size."
        )
    );

    // get cell dof indices
    std::vector<psize> dof_ids(fe.dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    // initialise cons_grad to 0
    for(usi dir=0; dir<dim; dir++){
        for(usi i=0; i<fe.dofs_per_cell; i++){
            for(cvar var: cvar_list) cons_grad[i][dir][var] = 0; // initialise to 0
        }
    }

    for(usi grad_dir=0; grad_dir<dim; grad_dir++){
        // first calculate volumetric contribution
        for(usi ldof_this=0; ldof_this<fe.dofs_per_cell; ldof_this++){
            State cons_this;
            for(cvar var: cvar_list) cons_this[var] = gcrk_cvars[var][dof_ids[ldof_this]];
            TableIndices<dim> ti_this;
            for(usi dir=0; dir<dim; dir++) ti_this[dir] = cdi.local_to_tensorial[ldof_this][dir];

            State flux; // flux between 'this' and 'other' states
            Tensor<1,dim> temp_dir; // temporary, doesn't matter for BR1 flux calculation

            for(usi m=0; m<=fe.degree; m++){
                for(usi m_dir=0; m_dir<dim; m_dir++){
                    State cons_other;
                    // strangely TableIndices doesn't have a copy ctor
                    TableIndices<dim> ti_other(ti_this[0], ti_this[1], ti_this[2]);
                    ti_other[m_dir] = m;
                    usi ldof_other = cdi.tensorial_to_local(ti_other);

                    for(cvar var: cvar_list){
                        cons_other[var] = gcrk_cvars[var][dof_ids[ldof_other]];
                    }
                    ns_ptr->get_aux_vol_flux(cons_this, cons_other, temp_dir, flux);
                    double JxContra_avg_comp = 0.5*(
                        metrics.at(cell->index()).JxContra_vecs[ldof_this][m_dir][grad_dir] +
                        metrics.at(cell->index()).JxContra_vecs[ldof_other][m_dir][grad_dir]
                    ); // component (of average contravariant vector) in gradient direction
                    for(cvar var: cvar_list){
                        cons_grad[ldof_this][grad_dir][var] +=
                            2*ref_D_1d(ti_this[m_dir],m)*flux[var]*JxContra_avg_comp;
                    }
                } // loop over three directions for m
            } // loop over m
        } // loop over cell dofs

        // now surface contribution for those dofs lying on face
        for(usi surf_dir=0; surf_dir<dim; surf_dir++){
            // loop over 'l'eft and 'r'ight faces
            for(usi lr_id=0; lr_id<=1; lr_id++){
                usi face_id = 2*surf_dir + lr_id; // the face id
                for(usi face_dof_id=0; face_dof_id<fe_face.dofs_per_face; face_dof_id++){
                    usi ldof = fdi.maps[face_id][face_dof_id];
                    State flux_in, flux_surf;
                    for(cvar var: cvar_list){
                        flux_in[var] = gcrk_cvars[var][dof_ids[ldof]]*
                            metrics.at(cell->index()).JxContra_vecs[ldof][surf_dir][grad_dir];
                        flux_surf[var] = s1_surf_flux[var].at(cell->index())[face_id][face_dof_id]*
                            metrics.at(cell->index()).JxContra_vecs[ldof][surf_dir][grad_dir];
                        // a hack for incorporating contributions from both left and right faces
                        cons_grad[ldof][grad_dir][var] -= std::pow(-1,lr_id)*
                            (flux_surf[var]-flux_in[var])/w_1d[lr_id*fe.degree];
                    }
                } // loop over face dofs
            } // loop over two surfaces with normal in 'surf_dir'
        } // loop over 3 directions for surface contributions

        // now divide by Jacobian determinant
        for(usi ldof=0; ldof<fe.dofs_per_cell; ldof++){
            for(cvar var: cvar_list){
                cons_grad[ldof][grad_dir][var] /= metrics.at(cell->index()).detJ[ldof];
            }
        } // loop over cell dofs
    } // loop over gradient directions
} // calc_cell_cons_grad



/**
 * Asserts the positivity of PLENS::gcrk_cvars. This is a requirement for using the DGSEM limiter
 * algo and NavierStokes functions. Uses NavierStokes::assert_positivity() for the task.
 */
void PLENS::assert_positivity() const
{
    State cons;
    for(auto i: locally_owned_dofs){
        for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
        ns_ptr->assert_positivity(cons);
    }
}



/**
 * Calculates the auxiliary variables and populates PLENS::gcrk_avars. Internally uses
 * calc_surf_flux() and calc_cell_cons_grad(). Obviously, PLENS::gcrk_cvars are used within these
 * functions to calculate the relevant quantities. The latter function gives the gradients of
 * conservative variables, from which the auxiliary variables are to be calculated. The relevant
 * formulae are
 * @f[
 * \tau_{ij} = \mu \left(
 * \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}
 * \right) - \frac{2\mu}{3} \delta_{ij} \frac{\partial u_k}{\partial x_k}\\
 * \vec{q} = -k\nabla T
 * @f]
 * The key is to obtain @f$\nabla \vec{u}@f$ and @f$\nabla T@f$ from the gradients of conservative
 * variables. For velocity
 * @f[
 * \nabla u_i = \frac{1}{\rho} \left[ \nabla (\rho u_i) - u_i \nabla \rho \right]
 * @f]
 * For energy
 * @f[
 * \rho \nabla e = \nabla(\rho E) - e \nabla \rho -
 * \frac{u_i u_i}{2} \nabla \rho - \rho u_i \nabla u_i
 * @f]
 * Since @f$e=c_v T@f$, and @f$c_v@f$ is constant, @f$\nabla T = \frac{\nabla e}{c_v}@f$.
 * All the scalar multiplication operations will be done using simple dof-wise multiplication.
 *
 * @pre Since division by density is required, positivity of density is required. Also, since
 * viscosity and thermal conductivity are calculated, positivity of energy is required. One way to
 * ensure this is to call assert_positivity()
 *
 * Algo
 * - Loop over owned dofs
 *   - Calculate viscosity and thermal conductivity
 * - Calculate surface flux for stage 1
 * - Loop over owned cells
 *   - Get conservative gradients in the cell
 *   - Loop over dofs
 *     - Calculate velocity gradients
 *     - Set stress auxiliary variables in PLENS::gcrk_avars (without factor @f$\mu@f$)
 *     - Calculate temperature gradient
 *     - Set heat flux auxiliary variables in PLENS::gcrk_avars (without factor @f$-k@f$)
 * - Multiply appropriate components of PLENS::gcrk_avars dof-wise with PLENS::gcrk_mu and
 * (-PLENS::gcrk_k)
 */
void PLENS::calc_aux_vars()
{
    // update viscosity and thermal conductivity
    const double cv = ns_ptr->get_cv();
    double e;
    State cons;
    for(psize i: locally_owned_dofs){
        for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
        e = ns_ptr->get_e(cons);
        gcrk_mu[i] = ns_ptr->get_mu(e/cv);
        gcrk_k[i] = ns_ptr->get_k(gcrk_mu[i]);
    } // loop over owned dofs
    gcrk_mu.compress(VectorOperation::insert);
    gcrk_k.compress(VectorOperation::insert);

    // calculate stage 1 flux
    locly_ord_surf_flux_term_t<double> s1_surf_flux;
    calc_surf_flux(1, s1_surf_flux);

    // now set (factored) avars, cell-by-cell
    std::vector<psize> dof_ids(fe.dofs_per_cell);
    std::vector<std::array<State, dim>> cons_grad(fe.dofs_per_cell);
    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        // calculate conservative gradients
        calc_cell_cons_grad(cell, s1_surf_flux, cons_grad);

        cell->get_dof_indices(dof_ids);

        for(usi i=0; i<fe.dofs_per_cell; i++){
            // cons declared above
            for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][dof_ids[i]];
            std::array<double, dim> vel{cons[1]/cons[0], cons[2]/cons[0], cons[3]/cons[0]};

            // calculate velocity gradient
            std::array<std::array<double, dim>, dim> vel_grad; // access: vel_grad[vel_dir][grad_dir]
            for(usi vel_dir=0; vel_dir<dim; vel_dir++){
                for(usi grad_dir=0; grad_dir<dim; grad_dir++){
                    vel_grad[vel_dir][grad_dir] = (
                        cons_grad[i][grad_dir][1+vel_dir] - // grad (rho u)
                        cons_grad[i][grad_dir][0]* // grad rho
                            vel[vel_dir] // u
                    )/cons[0];
                } // loop over gradient directions
            } // loop over velocity directions

            // set (factored) stress components
            // the loop is set such that it follows the ordering given in var_enums.h
            usi stress_id = 0;
            double vel_grad_trace = 0;
            for(usi d=0; d<dim; d++) vel_grad_trace += vel_grad[d][d];
            for(usi row=0; row<dim; row++){
                for(usi col=row; col<dim; col++){
                    if(col == row){
                        // somehow operator -= doesn't work together with operator=
                        // hence this split is required
                        gcrk_avars[stress_id][dof_ids[i]] = vel_grad[row][col] +
                            vel_grad[col][row] - 2.0/3*vel_grad_trace;
                    }
                    else{
                        gcrk_avars[stress_id][dof_ids[i]] = vel_grad[row][col] +
                            vel_grad[col][row];
                    }
                    stress_id++;
                }
            }

            // calculate temperature gradient
            std::array<double, dim> e_grad;
            e = ns_ptr->get_e(cons); // e declared above
            for(usi grad_dir=0; grad_dir<dim; grad_dir++){
                e_grad[grad_dir] = (cons_grad[i][grad_dir][4] - cons_grad[i][grad_dir][0]*e)/cons[0];
                for(usi vel_dir=0; vel_dir<dim; vel_dir++){
                    e_grad[grad_dir] -= (
                        0.5*cons_grad[i][grad_dir][0]*vel[vel_dir]*vel[vel_dir] +
                        cons[1+vel_dir]*vel_grad[vel_dir][grad_dir]
                    )/cons[0];
                }
            }

            // set (factored) heat flux components
            // the negative sign is added here so that the vectors can be scaled with gcrk_k
            // (instead of -gcrk_k)
            for(usi d=0; d<dim; d++) gcrk_avars[6+d][dof_ids[i]] = -e_grad[d]/cv;
        } // loop over cell dofs
    } // loop over owned cells

    // now compress and scale with mu and k, then set the ghosted vectors
    for(avar var: avar_list) gcrk_avars[var].compress(VectorOperation::insert);
    for(usi i=0; i<6; i++) gcrk_avars[i].scale(gcrk_mu);
    for(usi i=6; i<9; i++) gcrk_avars[i].scale(gcrk_k);
    for(avar var: avar_list) gh_gcrk_avars[var] = gcrk_avars[var];
} // calc_aux_vars



/**
 * Calculates the value of blender (@f$\alpha@f$). First PLENS::gcrk_blender_var is updated
 * according to "blender parameters/variable" entry of prm file. Then, the blender value is
 * calculated in each cell using BlenderCalculator::get_blender() and populated in
 * PLENS::gcrk_alpha. Then, PLENS::gh_gcrk_alpha is set for diffusing the value of alpha in each
 * cell.
 *
 * @note The indexing for PLENS::gcrk_alpha is based on `cell->global_active_cell_index()` which
 * was introduced in dealii-9.3.0, rather than `cell->index()` which is commonly used. This is
 * already explained in @ref cell_indices.
 *
 * The blender variable used in this function is simply based on PLENS::gcrk_cvars. No
 * fluxes/boundary conditions are required.
 */
void PLENS::calc_blender()
{
    // first update gcrk_blender_var
    prm.enter_subsection("blender parameters");
    {
        const std::string var_name = prm.get("variable");
        if(var_name == "rho") gcrk_blender_var = gcrk_cvars[0];
        else if(var_name == "p"){
            for(psize i: locally_owned_dofs){
                State cons;
                for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
                gcrk_blender_var[i] = ns_ptr->get_p(cons);
            }
        }
        else{
            // p times rho
            for(psize i: locally_owned_dofs){
                State cons;
                for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
                gcrk_blender_var[i] = ns_ptr->get_p(cons)*cons[0];
            }
        }
    }
    prm.leave_subsection();

    // now populate gcrk_alpha
    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        gcrk_alpha[cell->global_active_cell_index()] = blender_calc.get_blender(cell);
    } // loop over owned cells
    gcrk_alpha.compress(VectorOperation::insert);
    gh_gcrk_alpha = gcrk_alpha; // communicate

    // now diffuse
    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        for(usi face_id=0; face_id<n_faces_per_cell; face_id++){
            if(cell->face(face_id)->at_boundary()) continue;

            const double cur_alpha = gcrk_alpha[cell->global_active_cell_index()];
            const double nei_alpha =
                gh_gcrk_alpha[cell->neighbor(face_id)->global_active_cell_index()];
            gcrk_alpha[cell->global_active_cell_index()] = std::max(cur_alpha, 0.5*nei_alpha);
        }
    }
    gcrk_alpha.compress(VectorOperation::insert);
    gh_gcrk_alpha = gcrk_alpha; // for data output in write()
} // calc_blender



/**
 * Calculates high order cell residuals for stages 2 or 3 (inviscid and diffusive fluxes) in a
 * given `cell`. The algorithm is much similar to calc_cell_cons_grad() and uses PLENS::gcrk_cvars
 * and PLENS::gcrk_avars. Relevant equations are eqs. (B.14-15) of [1].
 *
 * The volumetric constribution is easily calculated by passing the two relevant states with the
 * average contravariant vector as the "direction" parameter to the volumetric flux functions of
 * NavierStokes class (see eq. (B.15) of [1] for what this means). Note that the average
 * contravariant vector need not be a unit vector. NavierStokes functions require a unit vector.
 * Hence, later, appropriately scale the flux obtained using the norm of contravariant vector.
 *
 * For the surface contribution, the numerical flux is already held by `s_surf_flux`. One important
 * note here. The algorithm described in appendix B of Hennemann et al (2021) solves the
 * transformed equations in a reference cell. The face normals to this cell are assumed there to be
 * in positive cartesian directions. However, `s_surf_flux` stores normal flux wrt "outward"
 * normals. This distinction is to be kept in mind. So for faces 1, 3 and 5, the numerical flux
 * stored in `s_surf_flux` can be used without sign change (scaling by norm of contravariant
 * vector would be required). For faces 0, 2 and 4, a sign change will also be required (because
 * outward normals and contravariant vectors face opposite to each other on these faces).
 *
 * Now for the internal flux, we simply need to get the volume flux in the direction of
 * contravariant vector followed by a scaling with its norm.
 *
 * @note Simply put, the contravariant vectors for dofs on faces give the face normal direction
 * such that it is outward for faces 1, 3 and 5, and inward for faces 0, 2 and 4. The sign changing
 * is required because contravariant and "outward" normal vectors don't match at the latter faces.
 *
 * The elements of `residual` will be reset here (not added to). The residual will be treated as
 * the "restriction" of global rhs to a cell, meaning its sign will be set such that it can be
 * treated as a rhs quantity. See for ref, chapter 2 of APS 1 report.
 *
 * @pre `stage` must be 2 or 3
 * @pre `residual` Must have the size `fe.dofs_per_cell`
 */
void PLENS::calc_cell_ho_residual(
    const usi stage,
    const DoFHandler<dim>::active_cell_iterator& cell,
    const locly_ord_surf_flux_term_t<double>& s_surf_flux,
    std::vector<State>& residual
) const
{
    AssertThrow(
        stage == 2 || stage == 3,
        StandardExceptions::ExcMessage(
            "Stage parameter passed to calc_cell_ho_residual() must be 2 or 3."
        )
    );
    AssertThrow(
        residual.size() == fe.dofs_per_cell,
        StandardExceptions::ExcMessage(
            "Residual vector passed to calc_cell_ho_residual() must have the size of dofs per cell"
        )
    );

    const usi stage_id = stage-1;
    // -1 for stage 2 (stage id 1) and +1 for stage 3 (stage id 2)
    const double stage_sign = std::pow(-1, stage_id);

    // get cell dof indices
    std::vector<psize> dof_ids(fe.dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    // reset values in residual to 0
    for(usi i=0; i<fe.dofs_per_cell; i++){
        for(cvar var: cvar_list) residual[i][var] = 0;
    }

    // set the residual
    // first volumetric contrib
    for(usi ldof_this=0; ldof_this<fe.dofs_per_cell; ldof_this++){
        State cons_this;
        for(cvar var: cvar_list) cons_this[var] = gcrk_cvars[var][dof_ids[ldof_this]];
        Avars av_this;
        for(avar var: avar_list) av_this[var] = gcrk_avars[var][dof_ids[ldof_this]];
        TableIndices<dim> ti_this;
        for(usi dir=0; dir<dim; dir++) ti_this[dir] = cdi.local_to_tensorial[ldof_this][dir];

        State flux; // flux between 'this' and 'other' states
        Tensor<1,dim> dir; // dir to be calculated based on avg of contravariant vectors
        double norm; // norm of the avg contravariant vector

        for(usi m=0; m<=fe.degree; m++){
            for(usi m_dir=0; m_dir<dim; m_dir++){
                State cons_other;
                Avars av_other;
                // strangely TableIndices doesn't have a copy ctor
                TableIndices<dim> ti_other(ti_this[0], ti_this[1], ti_this[2]);
                ti_other[m_dir] = m;
                usi ldof_other = cdi.tensorial_to_local(ti_other);

                for(cvar var: cvar_list) cons_other[var] = gcrk_cvars[var][dof_ids[ldof_other]];
                for(avar var: avar_list) av_other[var] = gcrk_avars[var][dof_ids[ldof_other]];
                dir = 0.5*(
                    metrics.at(cell->index()).JxContra_vecs[ldof_this][m_dir] +
                    metrics.at(cell->index()).JxContra_vecs[ldof_other][m_dir]
                ); // not unit vector yet
                norm = dir.norm();
                dir /= norm; // unit vector now
                CAvars cav_this(&cons_this, &av_this), cav_other(&cons_other, &av_other);
                ns_ptr->vol_flux_wrappers[stage_id](cav_this, cav_other, dir, flux);

                // for(cvar var: cvar_list) std::cout << "flux: " << flux[var] << "\n";

                // do an addition here, negative sign will be incorporated when scaling with
                // jacobian
                for(cvar var: cvar_list){
                    residual[ldof_this][var] += 2*ref_D_1d(ti_this[m_dir],m)*flux[var]*norm;
                }
            } // loop over m dir
        } // loop over m
    } // loop over cell dofs

    // now surface contribution
    for(usi surf_dir=0; surf_dir<dim; surf_dir++){
        // loop over 'l'eft and 'r'ight faces
        for(usi lr_id=0; lr_id<=1; lr_id++){
            usi face_id = 2*surf_dir + lr_id; // the face id
            for(usi face_dof_id=0; face_dof_id<fe_face.dofs_per_face; face_dof_id++){
                usi ldof = fdi.maps[face_id][face_dof_id];
                Tensor<1,dim> dir =
                    metrics.at(cell->index()).JxContra_vecs[ldof][surf_dir]; // not yet unit vector
                double norm = dir.norm();
                dir /= norm; // unit vector
                State flux_in, flux_surf;
                // surface flux
                for(cvar var: cvar_list){
                    // s_surf_flux assumes outward normal: true for right and false for left faces
                    // thus a sign changes is required for left faces since the algorithm assumes
                    // all face normals in positive cartesian directions (that's what the
                    // contravariant vectors give)
                    flux_surf[var] = -std::pow(-1, lr_id)* // -ve for left and +ve for right faces
                        s_surf_flux[var].at(cell->index())[face_id][face_dof_id];
                }
                // inner flux
                State cons;
                Avars av;
                for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][dof_ids[ldof]];
                for(avar var: avar_list) av[var] = gcrk_avars[var][dof_ids[ldof]];
                CAvars cav(&cons, &av);
                ns_ptr->flux_wrappers[stage_id](cav, dir, flux_in);

                for(cvar var: cvar_list){
                    residual[ldof][var] += -std::pow(-1,lr_id)*
                        (flux_surf[var]-flux_in[var])*norm/
                        w_1d[lr_id*fe.degree];
                }
            } // loop over face dofs
        } // loop over left/right faces
    } // loop over directions

    // multiply by sign and scale by jacobian
    for(usi i=0; i<fe.dofs_per_cell; i++){
        for(cvar var: cvar_list) residual[i][var] /= stage_sign*
            metrics.at(cell->index()).detJ[i];
    }
} // calc_cell_ho_residual



/**
 * Calculates the low order inviscid residual based on subcell update. This is the only place where
 * the subcell normals of metric terms are relevant. The algorithm used is slightly complicated.
 *
 * First, an outer loop of directions is followed by loops over dofs in complementary directions.
 * Then, a loop over the subcell interfaces of the outer direction (corresponding to outer loop) is
 * done where the required inter-subcell flux is calculated and its contributions are appropriately
 * added to the residuals. For each subcell interface, two dofs (subcells) get the flux
 * contribution. The number of such subcell normals (in each direction) will be @f$N+2@f$: one more
 * than the number of dofs. This is obvious because each dof is treated as a subcell and each
 * subcell will have two faces. The indexing of these subcell faces will follow the convention of
 * MetricTerms. See the member documentation of MetricTerms::subcell_normals for more details. Note
 * the boundary cases: when `i=0`, the subcell normal @f$\vec{n}_{(L,0)jk}@f$ is obtained, and when
 * `i=N+1` the subcell normal @f$\vec{n}_{(N,R)jk}@f$ is obtained when accessing subcell normals.
 *
 * For the boundary cases, the numerical flux will be obtained from `s2_surf_flux`. Again, remember
 * that `s2_surf_flux` stores fluxes assuming outward normals at all faces of a cell. These
 * directions match with the subcell normals (for those subcells lying on cell face) for faces 1, 3
 * and 5, while they will be opposite to subcell normals for the remaining faces (0, 2 and 4). See
 * for reference, the note in calc_cell_ho_residual().
 *
 * Algo:
 * - Loop over direction
 *   - Loop over ids in complementary direction 2 (id2=0 to id2=N) [start of internal contributions
 *     loop]
 *     - Loop over ids in complementary direction 1 (id1=0 to id1=N)
 *       - Loop over id in the direction (id=1 to id=N) (may be called direction 0)
 *         - Evaluate flux between states (id-1, id1, id2) and (id, id1, id2) [the indices have to
 *           be ordered properly, the case for direction=0 is taken as example here]
 *         - Add the contribution to dofs (id-1, id1, id2) and (id, id1, id2) [again, ordering must
 *           be modified appropriately]
 *   - For id==0 and id==N+1, use the surface flux from `s2_surf_flux` [surface contrib loop]
 *     - Loop over face dofs
 *       - Obtain the flux from `s2_surf_flux` and scale it appropriately with contravariant vector
 *         - Also, appropriate sign change is required because `s2_surf_flux` gives flux assuming
 *           outward pointing normal on cell faces, while contravariant vector points along
 *           reference cell axis
 *       - Add the contiribution to surface dof
 *         - The sign of contribution depends on whether the dof is on "left" face of "right" face
 * - Loop over all dofs
 *   - Scale residual by -jacobian (negative sign because above loops calculate divergence of flux
 *      but the RHS contains negative of the flux divergence)
 *
 * The elements of `residual` will be reset here (not added to). The residual will be treated as
 * the "restriction" of global rhs to a cell, meaning its sign will be set such that it can be
 * treated as a rhs quantity.
 *
 * @pre `s2_surf_flux` must be stage 2's surface flux
 * @pre `residual` Must have the size `fe.dofs_per_cell`
 */
void PLENS::calc_cell_lo_inv_residual(
    const DoFHandler<dim>::active_cell_iterator& cell,
    const locly_ord_surf_flux_term_t<double>& s2_surf_flux,
    std::vector<State>& residual
) const
{
    AssertThrow(
        residual.size() == fe.dofs_per_cell,
        StandardExceptions::ExcMessage(
            "Residual vector passed to calc_cell_lo_inv_residual() must have the size of dofs per "
            "cell."
        )
    );

    // reset values in residual to 0
    for(usi i=0; i<fe.dofs_per_cell; i++){
        for(cvar var: cvar_list) residual[i][var] = 0;
    }

    std::vector<psize> dof_ids(fe.dofs_per_cell);
    cell->get_dof_indices(dof_ids);

    // First the internal contributions
    for(usi dir=0; dir<dim; dir++){
        // complementary directions (required for internal contributions)
        usi dir1 = (dir+1)%dim;
        usi dir2 = (dir+2)%dim;
        for(usi id1=0; id1<=fe.degree; id1++){
            for(usi id2=0; id2<=fe.degree; id2++){
                for(usi id=1; id<=fe.degree; id++){
                    TableIndices<dim> ti_left, ti_right;
                    ti_left[dir] = id-1;
                    ti_left[dir1] = id1;
                    ti_left[dir2] = id2;
                    ti_right[dir] = ti_left[dir]+1; // equals "id"
                    ti_right[dir1] = ti_left[dir1];
                    ti_right[dir2] = ti_left[dir2];
                    usi ldof_left = cdi.tensorial_to_local(ti_left),
                        ldof_right = cdi.tensorial_to_local(ti_right);
                    
                    // set the "left" and "right" states
                    State cons_left, cons_right, flux;
                    for(cvar var: cvar_list){
                        cons_left[var] = gcrk_cvars[var][dof_ids[ldof_left]];
                        cons_right[var] = gcrk_cvars[var][dof_ids[ldof_right]];
                    }

                    // get normal between id-1 and id
                    Tensor<1,dim> normal_dir =
                        metrics.at(cell->index()).subcell_normals[dir](ti_right); // not unit vector yet
                    double normal_norm = normal_dir.norm();
                    normal_dir /= normal_norm; // unit vector now

                    // get the flux
                    ns_ptr->get_inv_surf_flux(cons_left, cons_right, normal_dir, flux);

                    // add the contribution to residual (eq. B.47 of Hennemann et al (2021))
                    // the negative sign for the surface divergence term will be assigned later,
                    // when dividing by jacobian
                    for(cvar var: cvar_list){
                        residual[ldof_left][var] += flux[var]*normal_norm/w_1d[id-1];
                        residual[ldof_right][var] -= flux[var]*normal_norm/w_1d[id];
                    }
                } // loop over internal ids in dir
            } // loop over ids in complementary dir 2
        } // loop over ids in complementary dir 1 (internal contributions loop)
    } // loop over directions

    // Now surface contributions
    for(usi dir=0; dir<dim; dir++){
        for(usi lr_id=0; lr_id<=1; lr_id++){
            usi face_id = 2*dir + lr_id;

            // sign of contribution (for subcell flux divergence)
            // negative for "left" faces (corresponding to (L,0) flux)
            // positive for "right" faces (corresponding to (N,R) flux)
            float contrib_sign = -std::pow(-1, lr_id);

            // flux sign: to align the flux provided by s2_surf_flux in the contravariant dir
            // s2_surf_flux provides flux on cell faces assuming outward normals
            // for "left" faces, contravariant vector points inwards ==> sign change required
            // for "right" faces, no sign change required
            float flux_sign = -std::pow(-1, lr_id);

            for(usi face_dof_id=0; face_dof_id<fe_face.dofs_per_face; face_dof_id++){
                usi ldof = fdi.maps[face_id].at(face_dof_id);

                // norm of contravariant vector
                // at end points, subcell normals equal contravariant vectors
                double normal_norm = metrics.at(cell->index()).JxContra_vecs[ldof][dir].norm();

                State flux; // flux wrt contravariant vector direction
                for(cvar var: cvar_list){
                    flux[var] = flux_sign*
                        s2_surf_flux[var].at(cell->index())[face_id][face_dof_id];
                    residual[ldof][var] += contrib_sign*flux[var]*normal_norm/
                        w_1d[lr_id*fe.degree];
                }
            } // loop over face dofs
        } // loop over "left" and "right" cell surfaces in dir
    } // loop over directions

    // Scale by negative jacobian (-ve because divergence of flux has -ve sign on RHS)
    for(usi i=0; i<fe.dofs_per_cell; i++){
        for(cvar var: cvar_list){
            residual[i][var] /= -metrics.at(cell->index()).detJ[i];
        }
    }
} // calc_cell_lo_inv_residual



/**
 * Calculates the rhs for all dofs and populates PLENS::gcrk_rhs. This function takes ownership of
 * the entire update residual calculation process. In other words, this function calculates the
 * residual using gcrk_cvars right from asserting positivity to calculating low order inviscid
 * contribution. The residual for conservative variables is set here using this simple algorithm
 * - Assert positivity
 * - Calculate auxiliary variables
 * - Calculate blender
 * - Calculate surface fluxes for stages 2 and 3
 * - Loop over owned cells
 *   - Get high order inviscid and viscous residual
 *   - Get low order inviscid residual
 *   - Calculate the final residual (using blender value)
 *   - Set the residual in gcrk_rhs
 */
void PLENS::calc_rhs()
{
    assert_positivity();
    {
        TimerOutput::Scope timer_section(timer, "Calculate auxiliary variables");
        calc_aux_vars();
    }
    {
        TimerOutput::Scope timer_section(timer, "Calculate blender");
        calc_blender();
    }
    
    // calculate stages 2 & 3 flux
    locly_ord_surf_flux_term_t<double> s2_surf_flux, s3_surf_flux;
    {
        TimerOutput::Scope timer_section(timer, "Calculate stages 2&3 surface fluxes");
        calc_surf_flux(2, s2_surf_flux);
        calc_surf_flux(3, s3_surf_flux);
    }

    // initialise variables
    std::vector<State> ho_inv_residual(fe.dofs_per_cell),
        ho_dif_residual(fe.dofs_per_cell),
        lo_inv_residual(fe.dofs_per_cell);
    std::vector<psize> dof_ids(fe.dofs_per_cell);

    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        {
            TimerOutput::Scope timer_section(timer, "Calculate stages 2&3 ho residual");
            calc_cell_ho_residual(
                2,
                cell,
                s2_surf_flux,
                ho_inv_residual
            ); // high order inviscid

            calc_cell_ho_residual(
                3,
                cell,
                s3_surf_flux,
                ho_dif_residual
            ); // high order viscous
        }

        {
            TimerOutput::Scope timer_section(timer, "Calculate lo inviscid residual");
            calc_cell_lo_inv_residual(
                cell,
                s2_surf_flux,
                lo_inv_residual
            );
        }

        cell->get_dof_indices(dof_ids);
        const double alpha = gcrk_alpha[cell->global_active_cell_index()];

        for(usi i=0; i<fe.dofs_per_cell; i++){
            for(cvar var: cvar_list){
                gcrk_rhs[var][dof_ids[i]] = ho_dif_residual[i][var] +
                    alpha*lo_inv_residual[i][var] +
                    (1-alpha)*ho_inv_residual[i][var];
            }
        } // loop over dofs
    } // loop over owned cells

    for(cvar var: cvar_list) gcrk_rhs[var].compress(VectorOperation::insert);
} // calc_rhs()



/**
 * Calculates time step based on PLENS::gcrk_cvars. The minimum of stable time step over all cells (across
 * all processes) is taken as the time step. The stable time step formula is taken from Hesthaven's
 * book.
 * @f[
 * \Delta t = \frac{h}{N^2} \frac{C}{ \left( |u|+|a| \right) + N^2\nu/h }
 * @f]
 * where @f$C@f$ is a constant and @f$h@f$ is the cell size.
 *
 * In this function, @f$C@f$ is taken 1 and @f$h@f$ is taken as the minimum vertex distance.
 *
 * @remark It was noted during pens2D project that the actual expression of Hethaven which uses
 * @f$\mu@f$ in place of @f$\nu@f$ in the above formula is dimensionally incorrect.
 *
 * @note Although time step calculation is required only once during an update, this function uses
 * PLENS::gcrk_cvars and not PLENS::g_cvars because PLENS::gcrk_mu will also be required. Thus, it
 * is expected that this function is called during the first RK step, before the first update and
 * after calling calc_aux_vars() because that is where PLENS::gcrk_mu will be set.
 *
 * @pre This function assumes PLENS::gcrk_mu is already set. Also assumes positivity of density.
 */
void PLENS::calc_time_step()
{
    TimerOutput::Scope timer_section(timer, "Calculate time step");
    std::vector<psize> dof_ids(fe.dofs_per_cell);

    double this_proc_step = 1e6; // stable time step (factored) for this process
    for(const auto& cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;

        cell->get_dof_indices(dof_ids);
        const double length = cell->minimum_vertex_distance();

        for(psize i: dof_ids){
            State cons;
            for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
            double a = ns_ptr->get_a(cons);
            double max_vel = std::max({
                fabs(cons[1]/cons[0]),
                fabs(cons[2]/cons[0]),
                fabs(cons[3]/cons[0])
            });
            // factoring out 1/N^2
            double cur_step = length/(max_vel + a +
                fe.degree*fe.degree*gcrk_mu[i]/(length*cons[0])
            );
            if(cur_step < this_proc_step) this_proc_step = cur_step;
        }
    } // loop over owned cells
    MPI_Barrier(mpi_comm);

    // first perform reduction (into "time_step" of 0-th process)
    MPI_Reduce(&this_proc_step, &time_step, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);

    // now multiply by Co/N^2 and broadcast
    if(Utilities::MPI::this_mpi_process(mpi_comm) == 0){
        // use smaller Co for first few time steps
        if(n_time_steps < 1) time_step *= 0.1*Co/(fe.degree*fe.degree);
        else time_step *= Co/(fe.degree*fe.degree);
    }
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, mpi_comm);
} // calc_time_step



/**
 * Calculates pressure, velocity and temperature from gcrk_cvars. Populates PLENS::gcrk_p,
 * PLENS::gcrk_T and PLENS::gcrk_vel.
 *
 * @pre PLENS::assert_positivity() must be successful before calling this function.
 */
void PLENS::post_process()
{
    for(psize i: locally_owned_dofs){
        State cons;
        for(cvar var: cvar_list) cons[var] = gcrk_cvars[var][i];
        gcrk_p[i] = ns_ptr->get_p(cons);
        gcrk_T[i] = ns_ptr->get_e(cons)/ns_ptr->get_cv();
        for(usi d=0; d<dim; d++) gcrk_vel[d][i] = gcrk_cvars[1+d][i]/gcrk_cvars[0][i];
    }

    gcrk_p.compress(VectorOperation::insert);
    gcrk_T.compress(VectorOperation::insert);
    for(usi d=0; d<dim; d++) gcrk_vel[d].compress(VectorOperation::insert);
}



/**
 * Calculates the steady state error. The parameter `cell_ss_error` is populated using
 * VectorTools::integrate_difference() and the global error is calculated using
 * VectorTools::compute_global_error() (which is returned). The global steady state error is
 * calculated as follows
 * @f[
 * e = \frac{\lVert (\rho E)^{n+1} - (\rho E)^n \rVert_{\Omega}}
 * {\lVert (\rho E)^{n+1} \rVert_{\Omega}}
 * @f]
 * Here, the norm is taken in L2 sense. The @f$n@f$ solution is taken as PLENS::rhoE_old and
 * @f$n+1@f$ is taken from PLENS::g_cvars. Thus, this function must be called after update() to get
 * the expected result (so that PLENS::g_cvars have the desired value).
 *
 * For both numerator and demoniator, VectorTools::integrate_difference() is called with
 * Function::ZeroFunction.
 *
 * @pre `cell_ss_error` must have a size `triang.n_active_cells()`. No checks on this are done.
 */
double PLENS::calc_ss_error(Vector<double>& cell_ss_error) const
{
    // first get the denominator
    VectorTools::integrate_difference(
        *mapping_ptr,
        dof_handler,
        g_cvars[4],
        Functions::ZeroFunction<dim>(),
        cell_ss_error,
        QGauss<dim>(fe.degree+1),
        VectorTools::NormType::L2_norm
    );

    const double denom = VectorTools::compute_global_error(
        triang,
        cell_ss_error,
        VectorTools::NormType::L2_norm
    );

    // now for numerator, compute a temporary dof vector for the difference
    LA::MPI::Vector rhoE_diff(locally_owned_dofs, mpi_comm);
    for(psize i: locally_owned_dofs) rhoE_diff[i] = g_cvars[4][i] - rhoE_old[i];
    rhoE_diff.compress(VectorOperation::insert);

    VectorTools::integrate_difference(
        *mapping_ptr,
        dof_handler,
        rhoE_diff,
        Functions::ZeroFunction<dim>(),
        cell_ss_error,
        QGauss<dim>(fe.degree+1),
        VectorTools::NormType::L2_norm
    );

    const double numer = VectorTools::compute_global_error(
        triang,
        cell_ss_error,
        VectorTools::NormType::L2_norm
    );
    
    return numer/denom; // denom guaranteed to be positive after passing assert_positivity()
}



/**
 * Performs a serialisation using solution transfer. If I understand correctly, serialisation a
 * process that converts dealii vectors to writable format. The purpose for saving solution vectors
 * is apparent. The file saved by this function can be read again. Before reading again, the
 * triangulation has to be read separately, and also its manifolds have to be set. Then, to read
 * the solution vectors, a solution transfer object has to be constructed again and initialised
 * with the dof handler. The dof handler used for reading must be attached to the triangulation
 * constructed like said above. Doing this procedure for parallel computations is slightly tricky.
 * Currently, this is not implemented, but is in plan. See ICs::FromArchive.
 *
 * According to dealii, writing a solution transfer requires ghosted vectors and reading requires
 * non-ghosted vectors. So, PLENS::gh_gcrk_cvars will be used here.
 *
 * @pre PLENS::gh_gcrk_cvars must be ready to use here.
 */
void PLENS::do_solution_transfer(const std::string& filename)
{
    std::vector<const LA::MPI::Vector*> gh_cvar_ptrs; // pointers to ghosted cvar vectors
    for(cvar var: cvar_list){
        gh_cvar_ptrs.emplace_back(&gh_gcrk_cvars[var]);
    }

    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> sol_trans(dof_handler);
    sol_trans.prepare_for_serialization(gh_cvar_ptrs);
    triang.save(filename);
}



/**
 * Writes the data. The following variables are written:
 * - PLENS::gh_gcrk_cvars
 * - PLENS::gh_gcrk_avars
 * - PLENS::gcrk_mu and PLENS::gcrk_k
 * - PLENS::gcrk_alpha
 * - Subdomain id (processor id)
 * - All other vectors calculated in PLENS::post_process()
 * - Cell wise steady state error vector (if requested)
 *
 * Since current RK solution is being written, it is expected that this, like
 * PLENS::calc_time_step(), will be called during the first RK step, after calculation of blender,
 * and before the first update.
 *
 * The location and naming of files is described in PLENS::declare_parameters(). Additionally, to
 * keep a log of data output time, a file `<base file name>.times` is written. This file contains
 * the output counter and the current time so that one can know at what times did the data output
 * happen. Additionally, a file `<base file name>.ss_error` is written to output the steady state
 * error at the output instant if the error calculation is requested. This calculation of error is
 * done using calc_ss_error()
 *
 * @note Although the documentation says `DataOut::add_data_vector()` requires a ghosted vector, a
 * non-ghosted but compressed vector also seems to do the job.
 *
 * This function writes the following files
 * - vtu files for individual processor data
 * - pvtu files for compiling processor data
 * - pvd files for compiling pvtu files across all outputs
 * - ".ar" archive files by solution transfer
 *
 * To be frank, the pvd file needs to be written only once after the entire simulation. However,
 * the cost incurred in writing it everytime is very low. Moreover, this is helpful in case the
 * simulation needs to be killed.
 *
 * This function also prints the timer output everytime it is called. The output will happen on
 * pcout. If the timer output data is of importance and needs to be stored, then the pcout output
 * can be logged into a file.
 */
void PLENS::write()
{
    timer.print_wall_time_statistics(mpi_comm);
    post_process();
    std::string op_dir, base_filename;
    bool calculate_ss_error;
    prm.enter_subsection("data output");
    {
        op_dir = prm.get("directory");
        base_filename = prm.get("base file name");
        calculate_ss_error = prm.get_bool("calculate steady state error");
    }
    prm.leave_subsection();

    DataOut<dim> data_out;
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    data_out.attach_dof_handler(dof_handler);
    for(cvar var: cvar_list) data_out.add_data_vector(gh_gcrk_cvars[var], cvar_names[var]);
    for(avar var: avar_list) data_out.add_data_vector(gh_gcrk_avars[var], avar_names[var]);

    data_out.add_data_vector(gcrk_mu, "mu");
    data_out.add_data_vector(gcrk_k, "k");
    data_out.add_data_vector(gcrk_p, "p");
    data_out.add_data_vector(gcrk_T, "T");
    data_out.add_data_vector(gcrk_vel[0], "u");
    data_out.add_data_vector(gcrk_vel[1], "v");
    data_out.add_data_vector(gcrk_vel[2], "w");

    // subdomain id
    Vector<float> subdom(triang.n_active_cells());
    for(float &x: subdom) x = triang.locally_owned_subdomain();
    data_out.add_data_vector(subdom, "Subdomain");

    // for alpha, direct addition not possible, see WJ-15-Jun-2021 and
    // https://groups.google.com/g/dealii/c/_lmP3VCLBsw
    // also see WJ-10-Jul-2021 and
    // https://groups.google.com/g/dealii/c/hF4AfBqnTdk
    Vector<float> alpha(triang.n_active_cells());
    for(auto &cell: dof_handler.active_cell_iterators()){
        if(!(cell->is_locally_owned())) continue;
        alpha[cell->active_cell_index()] = gh_gcrk_alpha[cell->global_active_cell_index()];
    }
    data_out.add_data_vector(alpha, "alpha");

    double ss_error;
    Vector<double> cell_ss_error(triang.n_active_cells());
    if(calculate_ss_error){
        ss_error = calc_ss_error(cell_ss_error);
        data_out.add_data_vector(cell_ss_error, "steady_state_error");
    }

    data_out.build_patches(
        *mapping_ptr,
        fe.degree,
        DataOut<dim>::CurvedCellRegion::curved_inner_cells
    );

    std::string master_filename = data_out.write_vtu_with_pvtu_record(
        op_dir + "/",
        base_filename,
        output_counter,
        mpi_comm,
        6
    ); // n_groups set to default value 0 (one file per processor)

    // pvd file
    if(Utilities::MPI::this_mpi_process(mpi_comm) == 0){
        times_and_names.emplace_back(
            cur_time,
            master_filename // name relative to pvd file path
        );
        std::ofstream pvd_file(op_dir + "/" + base_filename + ".pvd");
        AssertThrow(
            pvd_file.good(),
            StandardExceptions::ExcMessage(
                "Unable to open pvd file. Make sure the specified output directory exists."
            )
        );
        DataOutBase::write_pvd_record(pvd_file, times_and_names);
        pvd_file.close();
    } // root process

    // solution transfer
    // the file name for this is same as master filename, with a different extension
    // ".pvtu" substring is the last 5 characters
    std::string archive_filename = op_dir + "/" +
        master_filename.substr(0, master_filename.size()-5) + ".ar";
    do_solution_transfer(archive_filename);

    // append current time and output counter in <base file name>.times
    // also write steady state error if requested
    if(Utilities::MPI::this_mpi_process(mpi_comm) == 0){
        std::ofstream time_file(op_dir + "/" + base_filename + ".times", std::ios::app);
        AssertThrow(
            time_file.good(),
            StandardExceptions::ExcMessage("Unable to open times file.")
        );
        time_file << output_counter << " " << cur_time << " " << clk.wall_time() << "\n";
        time_file.close();

        if(calculate_ss_error){
            std::ofstream error_file(op_dir + "/" + base_filename + ".ss_error", std::ios::app);
            AssertThrow(
                error_file.good(),
                StandardExceptions::ExcMessage("Unable to open steady state error file.")
            );
            error_file << output_counter << " " << cur_time << " " << ss_error << "\n";
            error_file.close();
        }
    }
    output_counter++;
} // write()



/**
 * The update function for 5 stage 3 register RK4 algorithm. This is not intended to be used
 * separately, but through the wrapper function update(). See RK4Stage5Register3 class
 * documentation for details of the exact algorithm. Although the algorithm is a 3 register one,
 * since the solution data variable (viz. PLENS::gcrk_cvars) is not being overwritten here, the
 * implementation will be done like a 3N scheme: one register for solution and 3 registers for
 * residuals. See also WJ-09-Jun-2021.
 */
void PLENS::update_rk4()
{
    // stage 1
    calc_rhs();
    // time step must be calculated after calling calc_rhs(), otherwise gcrk_mu would be unset and
    // positivity would be unasserted
    calc_time_step();
    if(time_step > (end_time - cur_time)) time_step = end_time - cur_time;
    pcout << "Current time: " << cur_time
        << ", time step: " << time_step
        << ", elapsed wall time: " << clk.wall_time() << "\n";
    if(n_time_steps%write_freq == 0){
        write();
        pcout << "Writing solution\n";
    }
    pcout << "\tRK stage 1\n";
    for(cvar var: cvar_list){
        for(psize i: locally_owned_dofs){
            gcrk_cvars[var][i] = gcrk_cvars[var][i] +
                time_step*rk4_coeffs.a_outer[0]*gcrk_rhs[var][i];
            gprk_rhs[var][i] = gcrk_rhs[var][i];
        }
        gcrk_cvars[var].compress(VectorOperation::insert);
        gh_gcrk_cvars[var] = gcrk_cvars[var];
    }
    MPI_Barrier(mpi_comm);

    // stage 2
    calc_rhs();
    pcout << "\tRK stage 2\n";
    for(cvar var: cvar_list){
        for(psize i: locally_owned_dofs){
            gcrk_cvars[var][i] = gcrk_cvars[var][i] + time_step*(
                (rk4_coeffs.a_inner[0]-rk4_coeffs.a_outer[0])*gprk_rhs[var][i] +
                rk4_coeffs.a_outer[1]*gcrk_rhs[var][i]
            );
            gpprk_rhs[var][i] = gprk_rhs[var][i];
            gprk_rhs[var][i] = gcrk_rhs[var][i];
        }
        gcrk_cvars[var].compress(VectorOperation::insert);
        gh_gcrk_cvars[var] = gcrk_cvars[var];
    }
    MPI_Barrier(mpi_comm);

    // stages 3-5
    for(usi rk_stage=3; rk_stage<=5; rk_stage++){
        calc_rhs();
        pcout << "\tRK stage " << rk_stage << "\n";
        for(cvar var: cvar_list){
            for(psize i: locally_owned_dofs){
                gcrk_cvars[var][i] = gcrk_cvars[var][i] + time_step*(
                    (rk4_coeffs.b[rk_stage-3]-rk4_coeffs.a_inner[rk_stage-3])*gpprk_rhs[var][i] +
                    (rk4_coeffs.a_inner[rk_stage-2]-rk4_coeffs.a_outer[rk_stage-2])*
                        gprk_rhs[var][i] +
                    rk4_coeffs.a_outer[rk_stage-1]*gcrk_rhs[var][i]
                );
                gpprk_rhs[var][i] = gprk_rhs[var][i];
                gprk_rhs[var][i] = gcrk_rhs[var][i];
            }
            gcrk_cvars[var].compress(VectorOperation::insert);
            gh_gcrk_cvars[var] = gcrk_cvars[var];
        }
        MPI_Barrier(mpi_comm);
    }
} // update_rk4



/**
 * Updates the solution for one time step. This is a wrapper class around update_rk4() and
 * update_rk3(). Either of them is called depending on the RK degree provided in input file.
 * PLENS::g_cvars is used like a buffer here before and after the update. 
 */
void PLENS::update()
{
    // initialise
    for(cvar var: cvar_list){
        for(psize i: locally_owned_dofs){
            gcrk_cvars[var][i] = g_cvars[var][i];
        }
        gcrk_cvars[var].compress(VectorOperation::insert);
        gh_gcrk_cvars[var] = gcrk_cvars[var];
    }

    update_rk4();

    for(psize i: locally_owned_dofs) rhoE_old[i] = g_cvars[4][i];
    rhoE_old.compress(VectorOperation::insert);
    
    for(cvar var: cvar_list){
        for(psize i: locally_owned_dofs){
            g_cvars[var][i] = gcrk_cvars[var][i];
        }
    }

    cur_time += time_step;
    n_time_steps++;
} // update()



#ifdef DEBUG
void PLENS::test()
{
    utilities::Testing t("PLENS", "class");

    PLENS problem;
    problem.read_mesh();
}
#endif
