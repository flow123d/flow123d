lc = 0.2;
r  = 0.5;
R  = 2;
Point (1) = {R-r, 0, 0, lc};
Point (2) = {R, 0, 0, lc};
Point (3) = {R-r, r, 0, lc};
Point (4) = {R-2*r, 0, 0, lc};
Point (5) = {R-r, -r, 0, lc};
Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line Loop (3) = {-4, -3, -2, -1};
Plane Surface (3) = {3};

Point (6) = {0, 0, R-r, lc};
Point (7) = {0, 0, R, lc};
Point (8) = {0, r, R-r, lc};
Point (9) = {0, 0, R-2*r, lc};
Point (10) = {0, -r, R-r, lc};
Circle (5) = {7, 6, 8};
Circle (6) = {8, 6, 9};
Circle (7) = {9, 6, 10};
Circle (8) = {10, 6, 7};
Line Loop (4) = {5, 6, 7, 8};
Plane Surface (4) = {4};

Point (11) = {0, 0, 0};
Point (12) = {0, r, 0};
Point (13) = {0, -r, 0};

Circle (9)  = {2, 11, 7};
Circle (10) = {3, 12, 8};
Line Loop (5) = {1, 10, -5, -9};
Ruled Surface (5) = {5};

Circle (11)  = {4, 11, 9};
Line Loop (6) = {11, -6, -10, 2};
Ruled Surface (6) = {6};

Circle (12)  = {5, 13, 10};
Line Loop (7) = {4, 9, -8, -12};
Ruled Surface (7) = {7};

Line Loop (8) = {12, -7, -11, 3};
Ruled Surface (8) = {8};

Surface Loop (1) = {3, 4, 5, 6, 7, 8};
Volume (1) = {1};

Physical Surface("inflow")  = { 3 };
Physical Surface("wall")    = { 5, 6, 7, 8 };
Physical Surface("outflow") = { 4 };
Physical Volume(0)          = { 1 };
