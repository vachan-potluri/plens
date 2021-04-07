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
            Patterns::Anything(),
            "Mesh file name. No checks are done on the format"
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

        prm.enter_subsection("blunted double cone");
        {
            prm.declare_entry(
                "axis direction",
                "1 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a vector indicating axis direction. The vector may "
                "be scaled arbitrarily"
            );

            prm.declare_entry(
                "nose center",
                "0 0 0",
                Patterns::List(Patterns::Double(), dim, dim, " "),
                "Space-separated elements of a coordinate indicating the center of blunted nose"
            );

            prm.declare_entry(
                "blunt angle",
                "0",
                Patterns::Double(),
                "Angle from axis upto which blunting is done. This angle generally matches with "
                "the tangent angle to the cone that follows blunted section"
            );
        } // subsection blunted double cone
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
 * For curved meshes, only "cylinder flare" and "blunted double cone" geometries are supported. For
 * cylinder flare, the entire mesh is assigned CylindricalManifold. For blunted double cone, the
 * nose part is assigned SphericalManifold and the rest is assigned CylindricalManifold. See note
 * __pens2D to plens__ and entries around WJ-05-Apr-2021. The algorithm for setting manifold on
 * blunted double cone is
 *
 * - Loop over all cells
 *  - If the angle between line joining cell center to nose center and axis is less than blunt
 *    angle
 *    - Set manifold id 0 (for spherical manifold)
 *  - Else
 *    - Set manifold id 1 (for cylindrical manifold)
 *
 * Where blunt angle is the angle from axis upto which blunting is done. This angle generally
 * matches with the tangent angle to the cone that follows blunted section. Note that it is
 * implicitly assumed that the cone originating from nose center with cone angle equalling blunt
 * angle (this cone divides the mesh into two manifolds) doesn't pass through any cell. Meaning,
 * only cell interfaces lie on this dividing cone. This can be ensured by meshing the two regions
 * independently.
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
