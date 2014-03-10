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

Line(8) = {1, 2};
Line(9) = {3, 4};
Line Loop(10) = {8, 7, 9, -2};
Plane Surface(11) = {10};

Physical Line("1d") = {1};
Physical Surface("2d") = {11};

Physical Line(".2d_bottom_top") = {7, 2};
Physical Point(".1d") = {5, 6};

