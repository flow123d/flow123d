Ly=0.2;
Lx=0.4;
ChL=0.02;
// Default 2d algo (1) and MeshAdapt (4) produce nonsymetric mesh
// Delaunay (5) produce symmetric mesh
//Mesh.Algorithm=1;

Point(1) = {0,0,0,ChL};
Point(2) = {0,Ly,0,ChL};
Point(3) = {Lx,Ly,0,ChL};
Point(4) = {Lx,0,0,ChL};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = {9};


Physical Surface("bulk") = {10};
Physical Line(".left") = {5};
Physical Line(".right") = {7};


