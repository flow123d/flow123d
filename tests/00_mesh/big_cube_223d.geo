// Big cube with diagonal fractures.
// Noncompatible both 2d-3d and 2d-2d

size=10000; // 10km

//step = 0.0450; // big mesh 412248 el
step = 0.09; // middle 59496 el
step = 0.11; 
//step = 0.4; // small 1465 el

step_3d = step * size;
step_2d = step * size;



// left, front, bottom corner
X0 = -size;
X1 = -size;
X2 = -size;

// right, rear, top corner
Y0 = 1*size;
Y1 = 1*size;
Y2 = 1*size;


Point(2) = {X0, X1, X2, step_3d};
Point(3) = {Y0, X1, X2, step_3d};
Point(4) = {Y0, Y1, X2, step_3d};
Point(1) = {X0, Y1, X2, step_3d};
Point(5) = {X0, Y1, Y2, step_3d};
Point(8) = {Y0, Y1, Y2, step_3d};
Point(6) = {X0, X1, Y2, step_3d};
Point(7) = {Y0, X1, Y2, step_3d};


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
Surface Loop(26) = {24, 14, 22, 18, 20, 16};
Volume(26) = {26};

Physical Volume("3d") = {26};


Point(9)  = {X0, Y1, X2, step_2d};
Point(10) = {X0, X1, X2, step_2d};
Point(11) = {Y0, X1, X2, step_2d};
Point(12) = {Y0, Y1, X2, step_2d};

Point(13) = {X0, Y1, Y2, step_2d};
Point(14) = {X0, X1, Y2, step_2d};
Point(15) = {Y0, X1, Y2, step_2d};
Point(16) = {Y0, Y1, Y2, step_2d};


Line(51) = {13, 15};
Line(52) = {15, 11};
Line(53) = {11, 9};
Line(54) = {9, 13};
Line(55) = {14, 16};
Line(56) = {16, 12};
Line(57) = {12, 10};
Line(58) = {10, 14};
Line Loop(59) = {51, 52, 53, 54};
Plane Surface(40) = {59};
Line Loop(511) = {55, 56, 57, 58};
Plane Surface(41) = {511};


Physical Surface("2d") = {40, 41};

Physical Surface(".3d_sides") = {24, 22, 18, 16};
Physical Line(".2d") = {58, 54, 52, 56};
Physical Surface(".3d_top") = {14};
Physical Surface(".3d_bottom") = {20};
