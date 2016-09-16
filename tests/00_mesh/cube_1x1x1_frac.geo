mesh = 0.2;
Point(1) = {0, 0, -1, mesh};
Point(2) = {1, 0, -1, mesh};
Point(3) = {1, 1, -1, mesh};
Point(4) = {0, 1, -1, mesh};
Point(5) = {0, 1, 0, mesh};
Point(6) = {1, 1, 0, mesh};
Point(7) = {0, 0, 0, mesh};
Point(8) = {1, 0, 0, mesh};

Line(9) = {7, 8};
Line(10) = {8, 6};
Line(11) = {6, 5};
Line(12) = {5, 7};
Line(13) = {1, 2};
Line(14) = {2, 3};
Line(15) = {3, 4};
Line(16) = {4, 1};
Line(17) = {7, 1};
Line(18) = {8, 2};
Line(19) = {6, 3};
Line(20) = {5, 4};
Line Loop(21) = {12, 9, 10, 11};
Plane Surface(22) = {21};
Line Loop(23) = {19, -14, -18, 10};
Plane Surface(24) = {23};
Line Loop(25) = {11, 20, -15, -19};
Plane Surface(26) = {25};
Line Loop(27) = {16, -17, -12, 20};
Plane Surface(28) = {27};
Line Loop(29) = {15, 16, 13, 14};
Plane Surface(30) = {29};
Line Loop(31) = {9, 18, -13, -17};
Plane Surface(32) = {31};
Surface Loop(33) = {32, 22, 28, 30, 26, 24};
Volume(34) = {33};

Physical Surface("2d_fraction") = {30};
Physical Volume("3d_cube") = {34};

Physical Line(".2d_fraction") = {13, 16, 15, 14};
Physical Surface(".3d_cube") = {32, 28, 26, 24, 22};


