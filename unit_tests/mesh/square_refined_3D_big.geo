fine_step = 6.2699116e-5;
mesh = 0.21659695;

Point(1) = {0, 0, 0, fine_step};
Point(2) = {1, 0, 0, mesh};
Point(3) = {1, 1, 0, fine_step};
Point(4) = {0, 1, 0, mesh};
Point(5) = {0, 0, 1, fine_step};
Point(6) = {1, 0, 1, mesh};
Point(7) = {1, 1, 1, fine_step};
Point(8) = {0, 1, 1, mesh};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(13) = {4, 1, 2, 3};
Line Loop(14) = {8, 5, 6, 7};
Line Loop(15) = {3, 12, -7, -11};
Line Loop(16) = {1, 10, -5, -9};
Line Loop(17) = {10, 6, -11, -2};
Line Loop(18) = {9, -8, -12, 4};

Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};
Plane Surface(17) = {17};
Plane Surface(18) = {18};

Surface Loop(19) = {18, 13, 15, 14, 17, 16};

Volume(20) = {19};

Physical Volume("main_volume") = {20};
