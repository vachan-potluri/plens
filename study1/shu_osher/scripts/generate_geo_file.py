# This script generates a gmsh script of a 1d mesh for Shu-Osher problem.

## 16-Jun-2022: Modification
# The number of cells in x direction are divided among the regions [-5,-4] and [-4,5] so that
# an interface always exists at x=-4. So the integral part of "n_cells/10" will go to [-5,-4], and
# remaining to [-4,5].

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

n_x1 = int(args.n_cells[0]/10)
n_x2 = args.n_cells[0] - n_x1
content = """
w = 1;
Point(1) = {{-5, 0, 0}};
Point(2) = {{-4, 0, 0}};
Point(3) = {{5, 0, 0}};
Point(4) = {{5, w, 0}};
Point(5) = {{-4, w, 0}};
Point(6) = {{-5, w, 0}};

Line(1) = {{1,2}};
Line(2) = {{2,3}};
Line(3) = {{3,4}};
Line(4) = {{4,5}};
Line(5) = {{5,6}};
Line(6) = {{6,1}};
Line(7) = {{2,5}};

Transfinite Line{{1,5}} = {};
Transfinite Line{{2,4}} = {};
Transfinite Line{{3,6,7}} = {};

Line Loop(1) = {{1,7,5,6}};
Line Loop(2) = {{2,3,4,-7}};
For i In {{1:2}}
    Plane Surface(i) = {{i}};
    Transfinite Surface{{i}};
    Recombine Surface{{i}};
    out~{{i}}[] = Extrude {{0, 0, w}}{{
        Surface{{i}}; Layers{{ {} }}; Recombine;
    }};
EndFor

Physical Volume("volume", 1) = {{out_1[1], out_2[1]}};
Physical Surface("left", 1) = {{out_1[5]}}; // left
Physical Surface("right", 2) = {{out_2[3]}}; // right
Physical Surface("symmetry", 3) = {{
    out_1[2], out_1[4], 1, out_1[0],
    out_2[2], out_2[4], 2, out_2[0]
}}; // top and bottom, front and back

""".format(
    n_x1+1,
    n_x2+1,
    args.n_cells[1]+1,
    args.n_cells[2]
)

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

