lx = 16;
ly = 1;
lc = 0.2;

// Default Frontal 2d algo (1) produce nearly structured mesh
// MeshAdapt (4) produce ugly mesh
// Delaunay (5) produce nearly symmetric random mesh
Mesh.Algorithm=1;

Point(1) = {0,0,0,lc};
Point(2) = {0,ly,0,lc};
Point(3) = {lx,ly,0,lc};
Point(4) = {lx,0,0,lc};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = {9};

Physical Surface("domain") = {10};

Physical Line(".right") = {7};
Physical Line(".left") = {5};
