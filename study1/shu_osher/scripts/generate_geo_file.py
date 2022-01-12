# This script generates a gmsh script of a 1d mesh for Shu-Osher problem.

import argparse

parser = argparse.ArgumentParser(
    description = "A script to generate gmsh file to generate mesh for Shu-Osher problem."
)
parser.add_argument(
    "n_cells",
    type=int,
    nargs=3,
    help="Number of cells in the 3 directions in the mesh"
)
parser.add_argument("dest", help="Location where the gmsh file should be written")
parser.add_argument(
    "-f",
    "--filename",
    help="Name to be given to the file ('.geo' can be skipped). Default: 'mesh.geo'.",
    default="mesh.geo",
    action="store"
)
args = parser.parse_args()

content = """
w = 1;
Point(1) = {{-5, 0, 0}};
Point(2) = {{5, 0, 0}};
Point(3) = {{5, w, 0}};
Point(4) = {{-5, w, 0}};

Line(1) = {{1,2}};
Line(2) = {{2,3}};
Line(3) = {{3,4}};
Line(4) = {{4,1}};

Transfinite Line{{1,3}} = {};
Transfinite Line{{2,4}} = {};

Line Loop(1) = {{1,2,3,4}};
Plane Surface(1) = {{1}};
Transfinite Surface{{1}};
Recombine Surface{{1}};

out[] = Extrude {{0, 0, w}}{{
    Surface{{1}}; Layers{{ {} }}; Recombine;
}};

Physical Volume("volume", 1) = {{out[1]}};
Physical Surface("left", 1) = {{out[5]}}; // left
Physical Surface("right", 2) = {{out[3]}}; // right
Physical Surface("symmetry", 3) = {{out[2], out[4], 1, out[0]}}; // top and bottom, front and back

""".format(args.n_cells[0]+1, args.n_cells[1]+1, args.n_cells[2])

filename = args.filename
if filename[-4:] != ".geo":
    filename += ".geo"
dest = args.dest
if dest[-1] != "/":
    dest += "/"
full_filename = dest + filename
f = open(full_filename, "w") # overwrite content
f.write(content)
f.close()
print("Written file {}".format(full_filename))

