ChL=0.05;
// Default 2d algo (1) snd MeshAdapt (4) produce nonsymetric mesh
// Delaunay (5) produce symmetric mesh
Mesh.Algorithm=5;

Point(1) = {0,0,0,ChL};
Point(2) = {0,0,1,ChL};
Point(3) = {1,0,1,ChL};
Point(4) = {1,0,0,ChL};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

Line Loop(5) = {1, 2, 3, 4};

Plane Surface(6) = {5};
Physical Surface("plane") = {6};

Physical Line(".bottom") = {1};
Physical Line(".right") = {2};
Physical Line(".top") = {3};
Physical Line(".left") = {4};
