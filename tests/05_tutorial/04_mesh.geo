cl = 4.0; // characteristic length

Point(1) = {0, 0, 0, cl};
Point(2) = {100, 0, 0, cl};
Point(3) = {100, 0, 100, cl};
Point(4) = {0, 0, 100, cl};
Point(5) = {0, 0, 40, cl};
Point(6) = {0, 0, 60, cl};
Point(7) = {100, 0, 70, cl};
Point(8) = {100, 0, 20, cl};
Point(9) = {40, 0, 32, cl};
Point(10) = {70, 0, 32, cl};
Point(11) = {80, 0, 50, cl};
Point(12) = {50, 0, 65, cl};
Point(13) = {70, 0, 85, cl};
Point(14) = {90, 0, 85, cl};

Line(1) = {2, 1};
Line(2) = {5, 1};
Line(3) = {6, 5};
Line(4) = {4, 6};
Line(5) = {3, 4};
Line(6) = {7, 3};
Line(7) = {8, 7};
Line(8) = {2, 8};
Line(9) = {8, 9};
Line(10) = {5, 9};
Line(11) = {9, 10};
Line(12) = {11, 10};
Line(13) = {6, 12};
Line(14) = {7, 12};
Line(15) = {14, 13};
Line(16) = {12, 13};

Line Loop(19) = { 1, -2, -3, -4, -5, -6, -7, -8 };
Plane Surface(20) = { 19 };

Line { 9:16 } In Surface { 20 };


Physical Surface("rock") = {20};
Physical Line("flow_fracture1") = {13, 14};
Physical Line("flow_fracture2") = {10, 9};
Physical Line("deadend_fracture1") = {16, 15};
Physical Line("deadend_fracture2") = {11, 12};
Physical Line(".left") = {4, 3, 2};
Physical Line(".right") = {6, 7, 8};
Physical Line(".top") = {5};
Physical Line(".bottom") = {1};
Physical Point(".left_points") = {6, 5};
Physical Point(".right_points") = {7, 8};
