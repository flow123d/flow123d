d = 0.02;
cl = 0.1;
clf = 0.008;

Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 1, 0, cl };
Point(4) = { 0, 1, 0, cl };
Point(5) = { 0.5, 0.3-d/2, 0, clf };
Point(6) = { 1, 0.3-d/2, 0, clf };
Point(7) = { 1, 0.3+d/2, 0, clf };
Point(8) = { 0.5, 0.3+d/2, 0, clf };
Point(9) = { 0.5, 0.7-d/2, 0, clf };
Point(10) = { 1, 0.7-d/2, 0, clf };
Point(11) = { 1, 0.7+d/2, 0, clf };
Point(12) = { 0.5, 0.7+d/2, 0, clf };

Point(71) = { 1, 0.3+d/2+cl, 0, cl };
Point(101)= { 1, 0.7-d/2-cl, 0, cl };
Point(102) = { 0.5, 0.3+d/2+cl, 0, cl };
Point(103) = { 0.5, 0.7-d/2-cl, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 6 };
Line(3) = { 6, 7 };
Line(41) = { 7, 71 };
Line(42) = { 71, 101 };
Line(43) = { 101, 10 };
Line(5) = { 10, 11 };
Line(6) = { 11, 3 };
Line(7) = { 3, 4 };
Line(8) = { 4, 1 };
Line(9) = { 5, 6 };
Line(10) = { 7, 8 };
Line(11) = { 8, 5 };
Line(12) = { 9, 10 };
Line(13) = { 12, 11 };
Line(14) = { 12, 9 };
Line(101) = { 71, 102 };
Line(102) = { 102, 103 };
Line(103) = { 103, 101 };

Line Loop(1) = { 1, 2, -9, -11, -10, 41, 101, 102, 103, 43, -12, -14, 13, 6, 7, 8 };
Line Loop(11) = { 101:103, -42 };
Line Loop(2) = { 9, 3, 10, 11 };
Line Loop(3) = { 12, 5, -13, 14 };

Plane Surface(1) = { 1 };
Plane Surface(11) = { 11 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };

Physical Surface("rock") = { 1, 11 };
Physical Surface("fracture_lower") = { 2 };
Physical Surface("fracture_upper") = { 3 };
Physical Line(".bottom") = { 1 };
Physical Line(".top") = { 7 };
Physical Line(".left") = { 8 };
Physical Line(".right") = { 2, 41:43, 6 };
Physical Line(".right_fl") = { 3 };
Physical Line(".right_fu") = { 5 };



