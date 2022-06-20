cl = 0.1;

Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 0.5, 0, cl };
Point(4) = { 1, 1, 0, cl };
Point(5) = { 0, 1, 0, cl };
Point(6) = { 0.25, 0.5, 0, cl };
Point(7) = { 0.25, 0.25, 0, cl };
Point(8) = { 0.5, 0.5, 0, cl };
Point(9) = { 0.75, 0.75, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 5 };
Line(5) = { 5, 1 };
Line(7) = { 6, 8 };
Line(8) = { 8, 3 };
Line(9) = { 7, 8 };
Line(10) = { 8, 9 };

Line Loop(1) = { 1, 2, 3, 4, 5 };
Plane Surface(1) = { 1 };

Line { 7:10 } In Surface { 1 };

Physical Point(".right_frac") = { 3 };
Physical Line(".top") = { 4 };
Physical Line(".bottom") = { 1 };
Physical Line(".left") = { 5 };
Physical Line(".right") = { 2, 3 };
Physical Line("fracture_x") = { 7, 8 };
Physical Line("fracture_xy") = { 9, 10 };
Physical Surface("rock") = { 1 };

