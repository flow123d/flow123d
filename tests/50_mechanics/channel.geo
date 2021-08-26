l = 25;
h = 1;
d = 1e-5;
cl = 0.1;

Point(1) = { 0, 0, 0, cl };
Point(2) = { l, 0, 0, cl };
Point(3) = { l, h/2, 0, cl };
Point(4) = { 0, h/2, 0, cl };
Point(5) = { l, h, 0, cl };
Point(6) = { 0, h, 0, cl };

Line(1) = { 1,2 };
Line(2) = { 2,3 };
Line(3) = { 3,4 };
Line(4) = { 4,1 };
Line(5) = { 3,5 };
Line(6) = { 5,6 };
Line(7) = { 6,4 };

Line Loop(1) = { 1:4 };
Line Loop(2) = { -3, 5, 6, 7 };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };

Physical Surface("rock") = { 1,2 };
Physical Line("fracture") = { 3 };
Physical Line(".top") = { 6 };
Physical Line(".bottom") = { 1 };
Physical Line(".left") = { 4,7 };
Physical Line(".right") = { 2,5 };
Physical Point(".fracture_left") = { 4 };
Physical Point(".fracture_right") = { 3 };
