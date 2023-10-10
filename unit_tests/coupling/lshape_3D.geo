// fine_step = 2.75e-4;
// mesh = 0.95;
Point(1) = {0, 0, 0, mesh};
Point(2) = {2, 0, 0, mesh};
Point(3) = {2, 1, 0, mesh};
Point(4) = {1, 1, 0, fine_step};
Point(5) = {1, 2, 0, mesh};
Point(6) = {0, 2, 0, mesh};

Point(7) = {0, 0, 2, mesh};
Point(8) = {2, 0, 2, mesh};
Point(9) = {2, 1, 2, mesh};
Point(10) = {1, 1, 2, fine_step};
Point(11) = {1, 2, 2, mesh};
Point(12) = {0, 2, 2, mesh};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 7};

Line(13) = {1, 7};
Line(14) = {2, 8};
Line(15) = {3, 9};
Line(16) = {4, 10};
Line(17) = {5, 11};
Line(18) = {6, 12};

Line Loop(19) = {1, 2, 3, 4, 5, 6};

Line Loop(20) = {12, 11, 10, 9, 8, 7};

Line Loop(21) = {13, 7, -14, -1};
Line Loop(22) = {14, 8, -15, -2};
Line Loop(23) = {15, 9, -16, -3};
Line Loop(24) = {16, 10, -17, -4};
Line Loop(25) = {17, 11, -18, -5};
Line Loop(26) = {18, 12, -13, -6};

Plane Surface(19) = {19};
Plane Surface(20) = {20};
Plane Surface(21) = {21};
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};
Plane Surface(25) = {25};
Plane Surface(26) = {26};

Surface Loop(27) = {19, 21, 22, 23, 24, 25, 26, 20};

Volume(28) = {27};

Physical Volume("main_volume") = {28};
