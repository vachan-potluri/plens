import argparse
parser = argparse.ArgumentParser(
    description = "A script to generate gmsh file to generate mesh for forward facing step case."
)
parser.add_argument(
    "hinv",
    type=float,
    help="Inverse of dx and dy in the mesh"
)
parser.add_argument(
    "nz",
    type=int,
    help="Number of cells in the z direction"
)
args = parser.parse_args()

# check if 1/h is divisible by 6 and print a warning
if args.hinv % 6 != 0:
    # prints in yellow
    print(
        "\033[93m" +
        "Warning: the value of 1/h provided ({}) is not divisible by 6.".format(
            args.hinv
        ) +
        "\033[0m"
    )

content = """
h = 1/{}; // element size (1/h must ideally be divisible by 6)
nz = {}; // number of elements in z direction

// bottom
Point(1) = {{0, 0, 0}};
Point(2) = {{1/6, 0, 0}};
Point(3) = {{4, 0, 0}};

// top
Point(4) = {{0, 2, 0}};
Point(5) = {{1/6, 2, 0}};
Point(6) = {{4, 2, 0}};



// x-dir
Line(1) = {{1,2}};
Line(2) = {{2,3}};
Line(3) = {{4,5}};
Line(4) = {{5,6}};

// y-dir
Line(5) = {{1,4}};
Line(6) = {{2,5}};
Line(7) = {{3,6}};



Line Loop(1) = {{1,6,-3,-5}};
Line Loop(2) = {{2,7,-4,-6}};



Transfinite Line{{1,3}} = Ceil((1/6)/h) + 1;
Transfinite Line{{2,4}} = Ceil((4-1/6)/h) + 1;
Transfinite Line{{5,6,7}} = Ceil(2/h) + 1;



// 2d Mesh
For i In {{1:2}}
    Plane Surface(i) = {{i}};
    Transfinite Surface{{i}};
    Recombine Surface{{i}};
EndFor



// 3d Mesh
out1[] = Extrude {{0,0,0.5}} {{
    Surface{{1}}; Layers{{nz}}; Recombine;
}};
out2[] = Extrude {{0,0,0.5}} {{
    Surface{{2}}; Layers{{nz}}; Recombine;
}};

Physical Volume("vol",1) = {{out1[1], out2[1]}};
Physical Surface("post_shock",1) = {{out1[2], out1[4], out1[5]}};
Physical Surface("moving_shock",2) = {{out2[4]}};
Physical Surface("pre_shock",3) = {{out2[3]}};
Physical Surface("wall",4) = {{out2[2]}};
Physical Surface("back",5) = {{1,2}};
Physical Surface("front",6) = {{out1[0],out2[0]}};
""".format(
    args.hinv,
    args.nz
)

mesh_filename = "mesh_generated.geo"
f = open(mesh_filename, "w") # overwrite content
f.write(content)
f.close()
print("Written file {}. Copy/move it to your desired location.".format(mesh_filename))