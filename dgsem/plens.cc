/**
 * @file plens.cc
 * @brief The main class of this project.
 */

#include "plens.h"

/**
 * Constructor.
 */
PLENS::PLENS()
:
mpi_comm(MPI_COMM_WORLD),
pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm)==0)),
triang(mpi_comm)
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
            "Type of mesh. Options: curved|straight"
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
            Patterns::Anything(),
            "Mesh file name. No checks are done on the format"
        );
    } // subsection mesh
    prm.leave_subsection();

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

    // straight edge meshes
    if(type == "straight"){
        std::ifstream file(filename);
        AssertThrow(
            file.good(),
            StandardExceptions::ExcMessage(
                "Unable to open mesh file"
            )
        );
        grid_in.read(file, fmt);
        file.close();
    }
    else{
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "Curved meshes currently unsupported"
            )
        );
    }
    prm.leave_subsection();
}



#ifdef DEBUG
void PLENS::test()
{
    utilities::Testing t("PLENS", "class");

    PLENS problem;
    problem.read_mesh();
}
#endif
