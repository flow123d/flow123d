l1 = 0.1;
l2 = 0.1;

lc = 0.02;

Point(1) = { 0, 0, 0, lc };
Point(2) = { l1, 0, 0, lc };
Point(3) = { l1, 1, 0, lc };
Point(4) = { 0, 1, 0, lc };
Point(5) = { l1+l2, 0, 0, lc };
Point(6) = { l1+l2, 1, 0, lc };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };
Line(5) = { 2, 5 };
Line(6) = { 5, 6 };
Line(7) = { 6, 3 };

Line Loop (10) = { 1:4 };
Line Loop (20) = { 5, 6, 7, -2 };

Plane Surface(100) = { 10 };
Plane Surface(200) = { 20 };

Physical Point(".fracture_bottom") = { 2 };
Physical Point(".fracture_top") = { 3 };
Physical Line("fracture") = { 2 };
Physical Line(".left") = { 4 };
Physical Line(".right") = { 6 };
Physical Line(".rock_bottom") = { 1, 5 };
Physical Line(".rock_top") = { 3, 7 };
Physical Surface("rock") = { 100, 200 };
