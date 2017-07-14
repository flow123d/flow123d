h1=0.1; // step on 1D fracture
h2=0.1; // step on 2D domain

// Number of extrusion layers in Z direction
z_layers = 2;

z=2*h2;

Point(1) = {-1,1,0,h2};
Point(2) = {-1,-1,0,h2};
Point(3) = {1,-1,0,h2};
Point(4) = {1,1,0,h2};
Point(5) = {-1,0,0,h1};
Point(6) = {1,0,0,h1};

Point(11) = {-1,1,z,h2};
Point(12) = {-1,-1,z,h2};
Point(13) = {1,-1,z,h2};
Point(14) = {1,1,z,h2};
Point(15) = {-1,0,z,h1};
Point(16) = {1,0,z,h1};

Line(1) = {1, 11};
Line(2) = {5, 15};
Line(3) = {2, 12};
Line(4) = {3, 13};
Line(5) = {6, 16};
Line(6) = {4, 14};
Line(7) = {1, 4};
Line(8) = {4, 6};
Line(9) = {6, 3};
Line(10) = {3, 2};
Line(11) = {2, 5};
Line(12) = {5, 1};
Line(13) = {11, 14};
Line(14) = {14, 16};
Line(15) = {16, 13};
Line(16) = {13, 12};
Line(17) = {12, 15};
Line(18) = {15, 11};
Line(19) = {5, 6};
Line(20) = {16, 15};
Line Loop(21) = {18, 13, 14, 20};
Plane Surface(22) = {21};
Line Loop(23) = {20, -17, -16, -15};
Plane Surface(24) = {23};
Line Loop(25) = {7, 8, -19, 12};
Plane Surface(26) = {25};
Line Loop(27) = {19, 9, 10, 11};
Plane Surface(28) = {27};
Line Loop(29) = {19, 5, 20, -2};
Plane Surface(30) = {29};
Line Loop(31) = {17, -2, -11, 3};
Plane Surface(32) = {31};
Line Loop(33) = {18, -1, -12, 2};
Plane Surface(34) = {33};
Line Loop(35) = {13, -6, -7, 1};
Plane Surface(36) = {35};
Line Loop(37) = {14, -5, -8, 6};
Plane Surface(38) = {37};
Line Loop(39) = {15, -4, -9, 5};
Plane Surface(40) = {39};
Line Loop(41) = {16, -3, -10, 4};
Plane Surface(42) = {41};
Surface Loop(43) = {24, 32, 28, 40, 42, 30};
Volume(44) = {43};
Surface Loop(45) = {22, 34, 36, 38, 26, 30};
Volume(46) = {45};

Physical Surface("2d") = {30};
Physical Surface(".3d_bottom_top") = {42, 36};
Physical Volume("3d") = {44, 46};
Physical Line(".2d") = {2, 5};
