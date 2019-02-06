ChL=0.05;
// Default 2d algo (1) snd MeshAdapt (4) produce nonsymetric mesh
// Delaunay (5) produce symmetric mesh
Mesh.Algorithm=5;

Point(1) = {0,0,0,ChL};
Point(2) = {0,1,0,ChL};
Point(3) = {1,1,0,ChL};
Point(4) = {1,0,0,ChL};
Point(5) = {0,0.5,0,ChL};
Point(6) = {0.5,1,0,ChL};
Point(7) = {1,0.5,0,ChL};
Point(8) = {0.5,0,0,ChL};
Point(9) = {0.5,0.5,0,ChL};
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 6};
Line(4) = {6, 3};
Line(5) = {3, 7};
Line(6) = {7, 4};
Line(7) = {4, 8};
Line(8) = {8, 1};
Line(9) = {5, 9};
Line(10) = {6, 9};
Line(11) = {7, 9};
Line(12) = {8, 9};

Line Loop(20) = {1, 9, -12, 8};
Line Loop(21) = {2, 3, 10, -9};
Line Loop(22) = {-10, 4, 5, 11};
Line Loop(23) = {12, -11, 6, 7};

Plane Surface(30) = {20};
Physical Surface("square_sw") = {30};
Plane Surface(31) = {21};
Physical Surface("square_nw") = {31};
Plane Surface(32) = {22};
Physical Surface("square_ne") = {32};
Plane Surface(33) = {23};
Physical Surface("square_se") = {33};

Physical Line(".bc_west1") = {1};
Physical Line(".bc_west2") = {2};
Physical Line(".bc_north1") = {3};
Physical Line(".bc_north2") = {4};
Physical Line(".bc_east1") = {5};
Physical Line(".bc_east2") = {6};
Physical Line(".bc_south1") = {7};
Physical Line(".bc_south2") = {8};

Mesh 2;

Save "square4x_1x1.msh";
Exit;
