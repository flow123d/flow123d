h1=0.1; // step on 1D fracture
h2=0.1; // step on 2D domain

Point(1) = {-1,1,0,h2};
Point(2) = {-1,-1,0,h2};
Point(3) = {1,-1,0,h2};
Point(4) = {1,1,0,h2};

Point(5) = {-1,0,0,h1};
Point(6) = {1,0,0,h1};

Line(1) = {5, 6};
Line(2) = {1, 4};
Line(7) = {2, 3};

Line(8) = {1, 5};
Line(9) = {5, 2};
Line(10) = {3, 6};
Line(11) = {6, 4};
Line Loop(12) = {8, 1, 11, -2};
Plane Surface(13) = {12};
Line Loop(14) = {9, 7, 10, -1};
Plane Surface(15) = {14};

Physical Line("1d") = {1};
Physical Surface("2d") = {13, 15};

Physical Line(".2d_bottom_top") = {7, 2};
Physical Point(".1d") = {5, 6};
