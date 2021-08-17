cl = 0.3;
Point(1) = {-0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0.5, 1, 0, cl/10};
Point(5) = {0.5, 0.5, 0, cl/10};
Point(7) = {0, 1, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(6) = {4, 7};
Line(7) = {7, 1};

Line Loop(7) = {1, 2, 3, 4, -4, 6, 7};
Plane Surface(1) = { 7 };

Physical Surface("bulk") = { 1 };
Physical Line(".left_x") = { 7 };
Physical Line(".right_x") = { 2 };
Physical Line(".top_y") = { 3,6 };
Physical Line(".bottom_y") = { 1 };
Physical Line("fracture") = { 4,5 };
Physical Point(".fracture_out") = { 4 };
Physical Point(".fracture_in") = { 5 };
