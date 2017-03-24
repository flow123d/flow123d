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


Physical Surface(".3d_sides") = {24, 22, 18, 16};
Physical Surface(".3d_top") = {14};
Physical Surface(".3d_bottom") = {20};
Line(27) = {5, 20};
Line(28) = {20, 7};
Line(29) = {20, 8};
