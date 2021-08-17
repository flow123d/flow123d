d = 0.01;
cl = 0.1;
Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Point(5) = {0.5-d/2,   1, 0, cl/10};
Point(6) = {0.5+d/2,   1, 0, cl/10};
Point(7) = {0.5-d/2, 0.5, 0, cl/10};
Point(8) = {0.5+d/2, 0.5, 0, cl/10};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {4, 1};
Line(5) = {6, 8};
Line(6) = {8, 7};
Line(7) = {7, 5};
Line(8) = {5, 4};
Line(9) = {5, 6};

Line Loop(7) = {1, 2, 3, 5, 6, 7, 8, 4};
Line Loop(8) = {-7, -6, -5, -9};
Plane Surface(1) = { 7 };
Plane Surface(2) = { 8 };


Physical Surface("bulk") = { 1 };
Physical Surface("fracture") = { 2 };
Physical Line(".left_x") = { 4 };
Physical Line(".right_x") = { 2 };
Physical Line(".top_y") = { 3,8 };
Physical Line(".bottom_y") = { 1 };
Physical Line(".fracture_out") = { 9 };
//Physical Line(".fracture_in") = { 6 };
