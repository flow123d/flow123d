x = 0.1;
a = 1;
b = 2;
Point(1) = {a, 0, 0, x};
Point(2) = {a, b, 0, x};
Point(3) = {0, b, 0, x};
Point(4) = {0, 0, 0, x};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Surface(7) = {6};
