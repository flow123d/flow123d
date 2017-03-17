/*
  Cube mesh with incompatible fractures.
  Coherence - removes duplicate nodes in mesh.
*/
meshC = 0.25;    //cube
meshF1 = 0.25;   //fract 1
meshF2 = 0.15;   //fract 2
//cube
Point(1) = {0, 0, 0, meshC};
Point(2) = {1, 0, 0, meshC};
Point(3) = {1, 1, 0, meshC};
Point(4) = {0, 1, 0, meshC};
Point(5) = {0, 1, 1, meshC};
Point(6) = {1, 1, 1, meshC};
Point(7) = {0, 0, 1, meshC};
Point(8) = {1, 0, 1, meshC};

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

Line Loop(22) = {1, -3, -2, -4};
Plane Surface(23) = {22};
Line Loop(24) = {5, 6, 7, 4};
Plane Surface(25) = {24};
Line Loop(26) = {10, 11, -5, 2};
Plane Surface(27) = {26};
Line Loop(28) = {6, 8, -12, 11};
Plane Surface(29) = {28};
Line Loop(30) = {9, -1, -7, 8};
Plane Surface(31) = {30};
Line Loop(32) = {9, -3, 10, 12};
Plane Surface(33) = {32};
Surface Loop(34) = {23, 33, 31, 25, 27, 29};
Volume(35) = {34};


//fract 1
Point(100) = {0, 0, 0, meshF1};
Point(101) = {0, 0, 1, meshF1};
Point(102) = {1, 1, 1, meshF1};
Point(103) = {1, 1, 0, meshF1};

Line(104) = {100, 101};
Line(105) = {101, 102};
Line(106) = {102, 103};
Line(107) = {103, 100};

Line Loop(108) = {104, 105, 106, 107};
Plane Surface(109) = {108};

//fract 2
Point(200) = {0, 1, 0, meshF2};
Point(201) = {0, 1, 1, meshF2};
Point(202) = {1, 0, 1, meshF2};
Point(203) = {1, 0, 0, meshF2};

Line(204) = {200, 201};
Line(205) = {201, 202};
Line(206) = {202, 203};
Line(207) = {203, 200};

Line Loop(208) = {204, 205, 206, 207};
Plane Surface(209) = {208};


//bulk
Physical Surface("2d_fracture_1") = {109};
Physical Surface("2d_fracture_2") = {209};
Physical Volume("3d_cube") = {35};

//boundary
Physical Line(".2d_fracture_1") = {104, 105, 106, 107};
Physical Line(".2d_fracture_2") = {204, 205, 206, 207};
Physical Surface(".3d_cube") = { 23, 25, 27, 29, 31, 33};

Mesh 2;
Mesh 3;

Save "cube_2f_incomp.msh";
Exit;