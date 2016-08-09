cl_1 = 30;
Point(1) = {0, 250, 0, cl_1};
Point(2) = {0, 250, -300, cl_1};
Point(3) = {0, -250, 0, cl_1};
Point(4) = {0, -250, -300, cl_1};
Line(2) = {3, 4};
Line(4) = {2, 1};
Line(101) = {1, 3};
Line(102) = {4, 2};

cl_2 = 2;
Point(50) = {0, 0, -39, cl_2};
Point(51) = {0, 0, -37.2, cl_2};
Point(52) = {0, -1.8, -39, cl_2};
Point(54) = {0, 0, -40.8, cl_2};
Point(53) = {0, 1.8, -39, cl_2};


Circle(334) = {51, 50, 53};
Circle(335) = {53, 50, 54};
Circle(337) = {54, 50, 52};
Circle(338) = {52, 50, 51};

Line Loop(339) = {101, 2, 102, 4};
Line Loop(340) = {335, 337, 338, 334};
Plane Surface(1) = {339, 340};
Physical Line(".surface", 101) = {101};
Physical Line(".base", 102) = {102};
Physical Surface("rock", 1) = {1};
Physical Line(".tunnel", 103) = {334, 338, 337, 335};
