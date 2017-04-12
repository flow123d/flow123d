cl = 0.05;

Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 1, 0, cl };
Point(4) = { 0, 1, 0, cl };
Point(5) = { 0.5, 0.3, 0, cl };
Point(6) = { 1, 0.3, 0, cl };
Point(7) = { 0.5, 0.7, 0, cl };
Point(8) = { 1, 0.7, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 6 };
Line(3) = { 6, 8 };
Line(4) = { 8, 3 };
Line(5) = { 3, 4 };
Line(6) = { 4, 1 };
Line(7) = { 5, 6 };
Line(8) = { 7, 8 };

Line Loop(1) = { 1:6 };

Plane Surface(1) = { 1 };

Line { 7,8 } In Surface { 1 };

Physical Surface("rock") = { 1 };
Physical Line("fracture_lower") = { 7 };
Physical Line("fracture_upper") = { 8 };
Physical Line(".bottom") = { 1 };
Physical Line(".top") = { 5 };
Physical Line(".left") = { 6 };
Physical Line(".right") = { 2:4 };
Physical Point(".right_fl") = { 6 };
Physical Point(".right_fu") = { 8 };
