mesh = 0.2;
Point(1) = {0, 0, 0, mesh};
Point(2) = {1, 0, 0, mesh};
Point(3) = {1, 1, 0, mesh};
Point(4) = {0, 1, 0, mesh};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

Physical Surface("plane") = {6};
Physical Line(".plane") = {3, 2, 1};

Physical Line("channel") = {4};
Physical Point(".channel") = {4, 1};
