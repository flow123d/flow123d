/*
  Cube mesh with incompatible fractures.
  Coherence - removes duplicate nodes in mesh.
*/
meshC = 0.35;    //cube

//cube
Point(1) = {0, 0, 0, meshC};
Point(2) = {1, 0, 0, meshC};
Point(3) = {1, 1, 0, meshC};
Point(4) = {0, 1, 0, meshC};
Point(5) = {0, 1, 1, meshC};
Point(6) = {1, 1, 1, meshC};
Point(7) = {0, 0, 1, meshC};
Point(8) = {1, 0, 1, meshC};
Point(9) = {0.5, 0.5, 1, meshC};
Point(10) = {0.5, 0.5, 0, meshC};

Line(1) = {5, 4};
Line(2) = {7, 1};
Line(3) = {1, 4};
Line(4) = {5, 7};
Line(5) = {7, 8};
Line(6) = {8, 6};
Line(7) = {6, 5};
Line(8) = {6, 3};
Line(9) = {3, 4};
Line(10) = {1, 2};
Line(11) = {2, 8};
Line(12) = {2, 3};
Line(13) = {5, 9};
Line(14) = {9, 6};
Line(15) = {8, 9};
Line(16) = {9, 7};
Line(17) = {3, 10};
Line(18) = {10, 4};
Line(19) = {1, 10};
Line(20) = {10, 2};
Line Loop(21) = {17, 18, -9};
Plane Surface(22) = {21};
Line Loop(23) = {17, 20, 12};
Plane Surface(24) = {23};
Line Loop(25) = {19, 20, -10};
Plane Surface(26) = {25};
Line Loop(27) = {3, -18, -19};
Plane Surface(28) = {27};
Line Loop(29) = {13, 16, -4};
Plane Surface(30) = {29};
Line Loop(31) = {15, 14, -6};
Plane Surface(32) = {31};
Line Loop(33) = {7, 13, 14};
Plane Surface(34) = {33};
Line Loop(35) = {16, 5, 15};
Plane Surface(36) = {35};
Line Loop(37) = {1, -3, -2, -4};
Plane Surface(38) = {37};
Line Loop(39) = {5, -11, -10, -2};
Plane Surface(40) = {39};
Line Loop(41) = {1, -9, -8, 7};
Plane Surface(42) = {41};
Line Loop(43) = {6, 8, -12, 11};
Plane Surface(44) = {43};
//cube
Surface Loop(49) = {38, 28, 30, 44, 32, 24, 42, 22, 34, 40, 36, 26};
Volume(50) = {49};

// fracts
Line Loop(45) = {20, 11, 15, -13, 1, -18};
Plane Surface(46) = {45};
Line Loop(47) = {17, -19, -2, -16, 14, 8};
Plane Surface(48) = {47};



//bulk
Physical Surface("2d_fracture_1") = {46};
Physical Surface("2d_fracture_2") = {48};
Physical Volume("3d_cube") = {50};

//boundary
Physical Line(".2d_fracture_1") = {13, 15, 11, 20, 18, 1};
Physical Line(".2d_fracture_2") = {16, 14, 8, 17, 19, 2};
Physical Surface(".3d_cube") = {38, 30, 32, 34, 40, 44, 42, 22, 24, 26, 28, 36};

Mesh 2;
Mesh 3;

Save "cube_2f_incomp_SurfaceComp.msh";
Exit;
