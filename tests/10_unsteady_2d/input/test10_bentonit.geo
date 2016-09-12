ChL=0.005;
// Default 2d algo (1) snd MeshAdapt (4) produce nonsymetric mesh
// Delaunay (5) produce symmetric mesh
Mesh.Algorithm=5;

Point(1) = {0,0, 0.0,ChL};
Point(4) = {0,0,-0.25,ChL};

Line(5) = {1,4};


Physical Line("plane") = {5};

Physical Point(".bc_top") = {1};
Physical Point(".bc_bottom") = {4};
