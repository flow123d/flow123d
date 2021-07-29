// simple 3D cube

// following option is necessary for generation of very fine meshes
// random factor times mesh step has to be greater then machine precision 
Mesh.RandomFactor=1e-5;

// left, front, bottom corner
X0 = 0;
X1 = 0;
X2 = 0;

// right, rear, top corner
Y0 = 1;
Y1 = 1;
Y2 = 1;

// mesh step
mesh = 0.2;

Point(2) = {X0, X1, X2, mesh};
Point(3) = {Y0, X1, X2, mesh};
Point(4) = {Y0, Y1, X2, mesh};
Point(1) = {X0, Y1, X2, mesh};
Point(5) = {X0, Y1, Y2, mesh};
Point(8) = {Y0, Y1, Y2, mesh};
Point(6) = {X0, X1, Y2, mesh};
Point(7) = {Y0, X1, Y2, mesh};

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

Physical Volume("bulk") = {26};
Physical Surface(".top_z") = {14};
Physical Surface(".bottom_z") = {20};
Physical Surface(".front_y") = {18};
Physical Surface(".back_y") = {24};
Physical Surface(".left_x") = {16};
Physical Surface(".right_x") = {22};
