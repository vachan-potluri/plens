/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * Constructor. Calls declare_parameters() and parses parameter file. The parameter file should be
 * named 'input.prm'.
 */
PLENS::PLENS()
:
mpi_comm(MPI_COMM_WORLD),
pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm)==0)),
triang(mpi_comm),
ns_ptr(nullptr)
{
    declare_parameters();
    prm.parse_input("input.prm");
}



/**
 * Destructor.
 */
PLENS::~PLENS()
{}



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
        usi n_bc_max(12);
        for(usi i=0; i<n_bc_max; i++){
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
                    "|symmetry|periodic'. 'none' type cannot be specified for a physical boundary"
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
                    "right periodic boundary id",
                    "1",
                    Patterns::Integer(0,n_bc_max-1),
                    "The 'right' boundary id of the other periodic boundary of the mesh; the "
                    "boundary id of present section is treated as 'left'. See "
                    "https://www.dealii.org/current/doxygen/deal.II/namespaceGridTools.html#ab22eef800535f9e85a1723a6a36fd0f6. "
                    "Relevant for: periodic. It is important that ONLY ONE entry per periodic "
                    "boundary pair be present in the prm file."
                );
            }
            prm.leave_subsection(); // bid<x> subsection
        } // loop over boundary ids
    }
    prm.leave_subsection(); // subsection BCs

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
 * https://dealii.org/developer/doxygen/deal.II/group__simplex.html#ga058cd187cea704428ac1118410cd0fb8
 *
 * For straight edged meshes, the procedure is simple.
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
 * @note It is assumed that the cells are also strictly bifurcated: no cells have bifurcator plane
 * passing through them in the middle. This can be ensured by meshing the two regions separately.
 *
 * If there are any periodic BCs to be set, then `triang` object must be modified. This will be
 * done when BCs are set.
 */
void PLENS::read_mesh()
{
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
    }
    else{
        // curved type
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
                for(usi d=0; d<dim; d++) axis[d] = stod(splits[d]);

                temp = prm.get("axis point"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) axis_p[d] = stod(splits[d]);
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
                for(usi d=0; d<dim; d++) axis[d] = stod(splits[d]);

                temp = prm.get("separation point"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) separation_p[d] = stod(splits[d]);

                temp = prm.get("nose center"); // prm file guarantees this has exactly 3 doubles
                utilities::split_string(temp, " ", splits);
                for(usi d=0; d<dim; d++) nose_center[d] = stod(splits[d]);
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
}



/**
 * Forms the NavierStokes object based on settings in prm file. See declare_parameters() for the
 * settings.
 */
void PLENS::set_NS()
{
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
}



/**
 * Sets the dof handler object. Here a loop over BC section of prm file is done to check if there
 * are any periodic BCs. If there are, then periodicity is added appropriately.
 *
 * @note It is assumed that there is only one entry in the prm file for a pair of periodic
 * boundaries.
 *
 * @pre read_mesh() has to be called before this
 */
void PLENS::set_dof_handler()
{}



#ifdef DEBUG
void PLENS::test()
{
    utilities::Testing t("PLENS", "class");

    PLENS problem;
    problem.read_mesh();
}
#endif
