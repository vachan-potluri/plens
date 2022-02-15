h = 1/12; // element size (1/h must be divisible by 6)
nz = 1; // number of elements in z direction

// bottom
Point(1) = {0, 0, 0};
Point(2) = {1/6, 0, 0};
Point(3) = {4, 0, 0};

// top
Point(4) = {0, 2, 0};
Point(5) = {1/6, 2, 0};
Point(6) = {4, 2, 0};



// x-dir
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {4,5};
Line(4) = {5,6};

// y-dir
Line(5) = {1,4};
Line(6) = {2,5};
Line(7) = {3,6};



Line Loop(1) = {1,6,-3,-5};
Line Loop(2) = {2,7,-4,-6};



Transfinite Line{1,3} = (1/6)/h + 1;
Transfinite Line{2,4} = (4-1/6)/h + 1;
Transfinite Line{5,6,7} = 2/h + 1;



// 2d Mesh
For i In {1:2}
    Plane Surface(i) = {i};
    Transfinite Surface{i};
    Recombine Surface{i};
EndFor



// 3d Mesh
out1[] = Extrude {0,0,0.5} {
    Surface{1}; Layers{nz}; Recombine;
};
out2[] = Extrude {0,0,0.5} {
    Surface{2}; Layers{nz}; Recombine;
};

Physical Volume("vol",1) = {out1[1], out2[1]};
Physical Surface("post_shock",1) = {out1[2], out1[4], out1[5]};
Physical Surface("moving_shock",2) = {out2[4]};
Physical Surface("pre_shock",3) = {out2[3]};
Physical Surface("wall",4) = {out2[2]};
Physical Surface("back",5) = {1,2};
Physical Surface("front",6) = {out1[0],out2[0]};
