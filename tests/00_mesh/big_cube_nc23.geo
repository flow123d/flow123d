// Big cube with diagonal fractures.
// Noncompatible 2d-3d and compatible 2d-2d



size=10000; // 10km

//step = 0.0450; // big mesh 412248 el
step = 0.09; // middle 59496 el
//step = 0.4; // small 1465 el

step_3d = step * size;
step_2d = step * size;


// left, front, bottom corner
X0 = -size;
X1 = -size;
X2 = -size;

// right, rear, top corner
Y0 = size;
Y1 = size;
Y2 = size;

/*
Point(1) = {X0, Y1, X2, step_3d};
Point(2) = {X0, X1, X2, step_3d};

Point(3) = {Y0, Y1, X2, step_3d};
Point(4) = {X0, Y1, X2, step_3d};
Point(5) = {X0, Y1, Y2, step_3d};
Point(6) = {Y0, Y1, Y2, step_3d};
Point(7) = {X0, X1, Y2, step_3d};
Point(8) = {Y0, X1, Y2, step_3d};

Point(1) = {X0, Y1, X2, step_2d};
Point(2) = {X0, X1, X2, step_2d};

Point(3) = {Y0, Y1, X2, step_2d};
Point(4) = {X0, Y1, X2, step_2d};
Point(5) = {X0, Y1, Y2, step_2d};
Point(6) = {Y0, Y1, Y2, step_2d};
Point(7) = {X0, X1, Y2, step_2d};
Point(8) = {Y0, X1, Y2, step_2d};
*/



Point(1) = {X0, Y1, X2, step_3d};
Point(2) = {X0, X1, X2, step_3d};
Point(3) = {Y0, X1, X2, step_3d};
Point(4) = {Y0, Y1, X2, step_3d};
Point(5) = {X0, Y1, Y2, step_3d};
Point(6) = {X0, X1, Y2, step_3d};
Point(7) = {Y0, X1, Y2, step_3d};
Point(8) = {Y0, Y1, Y2, step_3d};

Point(11) = {X0, X1, X2, step_2d};
Point(12) = {Y0, X1, X2, step_2d};
Point(13) = {Y0, Y1, X2, step_2d};
Point(14) = {X0, Y1, X2, step_2d};
Point(15) = {X0, Y1, Y2, step_2d};
Point(16) = {Y0, Y1, Y2, step_2d};
Point(17) = {X0, X1, Y2, step_2d};
Point(18) = {Y0, X1, Y2, step_2d};
Point(19) = {0.5*(X0 + Y0), 0.5*(X1 + Y1), X2, step_2d};
Point(20) = {0.5*(X0 + Y0), 0.5*(X1 + Y1), Y2, step_2d};
Line(1) = {1, 2};
Line(2) = {6, 2};
Line(3) = {1, 5};
Line(4) = {4, 1};
Line(5) = {5, 8};
Line(6) = {4, 8};
Line(7) = {5, 6};
Line(8) = {2, 3};
Line(9) = {4, 3};
Line(10) = {3, 7};
Line(11) = {7, 8};
Line(12) = {7, 6};
Line(24) = {15, 20};
Line(25) = {20, 18};
Line(26) = {18, 12};
Line(27) = {12, 19};
Line(28) = {19, 20};
Line(29) = {19, 14};
Line(30) = {14, 15};
Line(31) = {17, 11};
Line(32) = {17, 20};
Line(33) = {11, 19};
Line(34) = {20, 16};
Line(35) = {13, 16};
Line(36) = {13, 19};
Line Loop(14) = {5, -11, 12, -7};
Plane Surface(14) = {14};
Line Loop(16) = {3, 7, 2, -1};
Plane Surface(16) = {16};
Line Loop(18) = {8, 10, 12, 2};
Plane Surface(18) = {18};
Line Loop(20) = {4, 1, 8, -9};
Plane Surface(20) = {20};
Line Loop(22) = {6, -11, -10, -9};
Plane Surface(22) = {22};
Line Loop(24) = {5, -6, 4, 3};
Plane Surface(24) = {24};
Line Loop(38) = {24, -28, 29, 30};
Plane Surface(38) = {38};
Line Loop(40) = {28, 34, -35, 36};
Plane Surface(40) = {40};
Line Loop(42) = {28, 25, 26, 27};
Plane Surface(42) = {42};
Line Loop(44) = {28, -32, 31, 33};
Plane Surface(44) = {44};
Surface Loop(26) = {14, 16, 18, 20, 22, 24};
Volume(26) = {26};
Physical Surface(".3d_sides") = {16, 18, 22, 24};
Physical Surface(".3d_top") = {14};
Physical Surface(".3d_bottom") = {20};
Physical Surface("2d") = {38, 40, 42, 44};
Physical Volume("3d") = {26};
