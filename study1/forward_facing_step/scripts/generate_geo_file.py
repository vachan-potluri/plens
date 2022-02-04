import numpy as np
from numpy.polynomial import Polynomial
def get_prog_ratio(length, start_size, n_cells):
    # S = sum (length)
    # a = start size
    # S/a = \frac{r^n - 1}{r-1}
    # S, a and n are known, r is required. A nonlinear polynomial equation.
    coeffs = np.zeros(n_cells+1)
    coeffs[-1] = 1 # r^n
    coeffs[0] = -1+length/start_size # constant
    coeffs[1] = -length/start_size # r^1
    # https://numpy.org/doc/stable/reference/routines.polynomials.html
    p = Polynomial(coeffs)
    r = p.roots()
    real_roots = r.real[abs(r.imag)<1e-4]
    valid_roots = real_roots[real_roots>1]
    assert len(valid_roots) == 1, "Found more than one valid roots"
    return valid_roots[0]

import argparse
parser = argparse.ArgumentParser(
    description = "A script to generate gmsh file to generate mesh for forward facing step case."
)
parser.add_argument(
    "h",
    type=float,
    help="dx and dy in the mesh"
)
parser.add_argument(
    "nz",
    type=int,
    help="Number of cells in the z direction"
)
parser.add_argument(
    "cr_ratio",
    type=float,
    help="Corner refinement ratio. The size at corner is calculated as 'h' divided by this ratio."
)
args = parser.parse_args()



content = """
h = {}; // element size



// Points
Point(1) =  {{0,    0,  0}};
Point(2) =  {{0.4,  0,  0}};
Point(3) =  {{0.6,  0,  0}};

Point(4) =  {{0,    0.2,  0}};
Point(5) =  {{0.4,  0.2,  0}};
Point(6) =  {{0.6,  0.2,  0}};
Point(7) =  {{0.8,  0.2,  0}};
Point(8) =  {{3,    0.2,  0}};

Point(9) =  {{0,    0.4,  0}};
Point(10) = {{0.4,  0.4,  0}};
Point(11) = {{0.6,  0.4,  0}};
Point(12) = {{0.8,  0.4,  0}};
Point(13) = {{3,    0.4,  0}};

Point(14) = {{0,    1,  0}};
Point(15) = {{0.4,  1,  0}};
Point(16) = {{0.6,  1,  0}};
Point(17) = {{0.8,  1,  0}};
Point(18) = {{3,    1,  0}};



// Lines
// x-dir
Line(1) =  {{1,  2}};
Line(2) =  {{2,  3}};

Line(3) =  {{4,  5}};
Line(4) =  {{5,  6}};
Line(5) =  {{6,  7}};
Line(6) =  {{7,  8}};

Line(7) =  {{9,  10}};
Line(8) =  {{10, 11}};
Line(9) =  {{11, 12}};
Line(10) = {{12, 13}};

Line(11) = {{14, 15}};
Line(12) = {{15, 16}};
Line(13) = {{16, 17}};
Line(14) = {{17, 18}};

// y-dir
Line(15) = {{1,  4}};
Line(16) = {{4,  9}};
Line(17) = {{9,  14}};

Line(18) = {{2,  5}};
Line(19) = {{5,  10}};
Line(20) = {{10, 15}};

Line(21) = {{3,  6}};
Line(22) = {{6,  11}};
Line(23) = {{11, 16}};

Line(24) = {{7,  12}};
Line(25) = {{12, 17}};

Line(26) = {{8,  13}};
Line(27) = {{13, 18}};




// Line loops
Curve Loop(1) =  {{1,   18,  -3,  -15}};
Curve Loop(2) =  {{2,   21,  -4,  -18}};
Curve Loop(3) =  {{3,   19,  -7,  -16}};
Curve Loop(4) =  {{4,   22,  -8,  -19}};
Curve Loop(5) =  {{5,   24,  -9,  -22}};
Curve Loop(6) =  {{6,   26, -10,  -24}};
Curve Loop(7) =  {{7,   20, -11,  -17}};
Curve Loop(8) =  {{8,   23, -12,  -20}};
Curve Loop(9) =  {{9,   25, -13,  -23}};
Curve Loop(10) = {{10,  27, -14,  -25}};



// Transfinite specs
// x-dir
Transfinite Curve{{1,3,7,11}} = 0.4/h + 1;
Transfinite Curve{{2,8,12,  9,13}} = 0.2/h + 1;
Transfinite Curve{{6,10,14}} = 2.2/h + 1;

// y-dir
Transfinite Curve{{15,18,  16,19,24,26}} = 0.2/h + 1;
Transfinite Curve{{17,20,23,25,27}} = 0.6/h + 1;

// originating at corner
Transfinite Curve{{-4,5,-21,22}} = 0.2/h + 1 Using Progression {};



// 2d Mesh
For i In {{1:10}}
    Plane Surface(i) = {{i}};
    Transfinite Surface{{i}};
    Recombine Surface{{i}};
EndFor

// extrude
// and boundaries
// 1-wall, 2-outflow, 3-inflow, 4-symmetry
Physical Surface(4) = {{}};
Physical Volume(1) = {{}};
For i In {{1:10}}
    out[] = Extrude {{0,0,1}} {{
        Surface{{i}}; Layers{{{}}}; Recombine;
    }};
    
    Physical Surface(4) += {{i,out[0]}};
    Physical Volume(1) += {{out[1]}};
    If (i == 1)
        Physical Surface(1) = {{out[2]}};
        Physical Surface(3) = {{out[5]}};
    EndIf
    If (i == 2)
        Physical Surface(1) += {{out[2], out[3]}};
    EndIf
    If (i == 3)
        Physical Surface(3) += {{out[5]}};
    EndIf
    If (i == 5)
        Physical Surface(1) += {{out[2]}};
    EndIf
    If (i == 6)
        Physical Surface(1) += {{out[2]}};
        Physical Surface(2) = {{out[3]}};
    EndIf
    If (i == 7)
        Physical Surface(1) += {{out[4]}};
        Physical Surface(3) += {{out[5]}};
    EndIf
    If (i == 8 || i == 9)
        Physical Surface(1) += {{out[4]}};
    EndIf
    If (i == 10)
        Physical Surface(1) += {{out[4]}};
        Physical Surface(2) += {{out[3]}};
    EndIf
EndFor
""".format(
    args.h,
    get_prog_ratio(0.2, args.h/args.cr_ratio, int(0.2/args.h)),
    args.nz
)

mesh_filename = "mesh_generated.geo"
f = open(mesh_filename, "w") # overwrite content
f.write(content)
f.close()
print("Written file {}. Copy/move it to your desired location.".format(mesh_filename))