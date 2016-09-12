ChL=0.05;
// Default 2d algo (1) snd MeshAdapt (4) produce nonsymetric mesh
// Delaunay (5) produce symmetric mesh
Mesh.Algorithm=5;

Point(1) = {0,0,0,ChL};
Point(2) = {0,1,0,ChL};
Point(3) = {1,1,0,ChL};
Point(4) = {1,0,0,ChL};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = {9};

Physical Surface("plane") = {10};

Physical Line(".bc_inflow") = {7};
Physical Line(".bc_outflow") = {5};
