size=10000; // 10km

h1=0.4*size; // step on 1D fracture
h2=0.4*size; // step on 2D domain


// left, front, bottom corner
X0 = -size;
X1 = -size;
X2 = -size;

// right, rear, top corner
Y0 = 1*size;
Y1 = 1*size;
Y2 = 1*size;


Point(2) = {X0, X1, X2, h2};
Point(3) = {Y0, X1, X2, h2};
Point(4) = {Y0, Y1, X2, h2};
Point(1) = {X0, Y1, X2, h2};
Point(5) = {X0, Y1, Y2, h2};
Point(8) = {Y0, Y1, Y2, h2};
Point(6) = {X0, X1, Y2, h2};
Point(7) = {Y0, X1, Y2, h2};


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
Physical Surface(".3d_sides") = {24, 22, 18, 16};
Physical Surface(".3d_top") = {14};
Physical Surface(".3d_bottom") = {20};

Point(11) = {X0, X1, X2, h1};
Point(12) = {Y0, X1, X2, h1};
Point(13) = {Y0, Y1, X2, h1};
Point(14) = {X0, Y1, X2, h1};
Point(15) = {X0, Y1, Y2, h1};
Point(16) = {Y0, Y1, Y2, h1};
Point(17) = {X0, X1, Y2, h1};
Point(18) = {Y0, X1, Y2, h1};
Point(19) = {0.5*(X0 + Y0), 0.5*(X1 + Y1), X2, h1};
Point(20) = {0.5*(X0 + Y0), 0.5*(X1 + Y1), Y2, h1};

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
Line Loop(37) = {24, -28, 29, 30};
Plane Surface(38) = {37};
Line Loop(39) = {28, 34, -35, 36};
Plane Surface(40) = {39};
Line Loop(41) = {28, 25, 26, 27};
Plane Surface(42) = {41};
Line Loop(43) = {28, -32, 31, 33};
Plane Surface(44) = {43};
Physical Surface("2d") = {38, 40, 42, 44};
