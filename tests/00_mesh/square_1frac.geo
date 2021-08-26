cl = 0.05;

Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 1, 0, cl };
Point(4) = { 0, 1, 0, cl };
Point(5) = { 0.5, 0.5, 0, cl };
Point(6) = { 1, 0.5, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 6 };
Line(3) = { 6, 3 };
Line(4) = { 3, 4 };
Line(5) = { 4, 1 };

Line(6) = {5, 6};

Line Loop(1) = { 1:5 };

Plane Surface(1) = { 1 };

Line { 6 } In Surface { 1 };

Physical Surface("rock") = { 1 };
Physical Line("fracture") = { 6 };
Physical Line(".bottom") = { 1 };
Physical Line(".top") = { 4 };
Physical Line(".left") = { 5 };
Physical Line(".right") = { 2, 3 };
Physical Point(".right_fl") = { 6 };


