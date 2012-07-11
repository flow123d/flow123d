
cl1 = 9.85;

Point(1) = {0, 0, 0, cl1};
Point(2) = {50, 0, 0, cl1};
Point(3) = {100, 0, 0, cl1};
Point(4) = {100, 100, 0, cl1};
Point(5) = {50, 100, 0, cl1};
Point(6) = {0, 100, 0, cl1};

Point(11) = {0, 0, 100, cl1};
Point(12) = {50, 0, 100, cl1};
Point(13) = {100, 0, 100, cl1};
Point(14) = {100, 100, 100, cl1};
Point(15) = {50, 100, 100, cl1};
Point(16) = {0, 100, 100, cl1};


Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};
Line(8) = {1, 11};
Line(9) = {3, 13};
Line(10) = {4, 14};
Line(11) = {6, 16};
Line(12) = {11, 12};
Line(13) = {12, 15};
Line(14) = {15, 16};
Line(15) = {16, 11};
Line(16) = {12, 13};
Line(17) = {13, 14};
Line(18) = {14, 15};

Line Loop(19) = {1, 2, 3, 4};
Plane Surface(20) = {19};
Line Loop(21) = {2, -7, -6, -5};
Plane Surface(22) = {21};
Line Loop(23) = {4, 8, -15, -11};
Plane Surface(24) = {23};
Line Loop(25) = {8, 12, 16, -9, -5, -1};
Plane Surface(26) = {25};
Line Loop(27) = {6, 10, -17, -9};
Plane Surface(28) = {27};
Line Loop(29) = {10, 18, 14, -11, -3, -7};
Plane Surface(30) = {29};
Line Loop(31) = {16, 17, 18, -13};
Plane Surface(32) = {31};
Line Loop(33) = {13, 14, 15, 12};
Plane Surface(34) = {33};

Surface Loop(35) = {26, 24, 20, 22, 30, 28, 32, 34};
Volume(36) = {35};


Physical Surface(1) = {34};
Physical Surface(2) = {22};

Physical Volume(1) = {36};

