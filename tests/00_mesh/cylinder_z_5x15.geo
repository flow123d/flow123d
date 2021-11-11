r = 5;    // radius
h = 15;   // height
cl = 3; // characteristic length of mesh elements

Point(1) = { 0, 0, 0, cl };
Point(2) = { r, 0, 0, cl };
Point(3) = { r*Cos(Pi*2/3), r*Sin(Pi*2/3), 0, cl };
Point(4) = { r*Cos(Pi*4/3), r*Sin(Pi*4/3), 0, cl };

Circle(1) = { 2, 1, 3 };
Circle(2) = { 3, 1, 4 };
Circle(3) = { 4, 1, 2 };
Line Loop(1) = { 1:3 };
Plane Surface(1) = { 1 };

Extrude{ 0, 0, -h }{ Surface{ 1 }; Layers{ h/cl}; }

Physical Volume("cylinder") = { 1 };
Physical Surface(".top") = { 1 };
Physical Surface(".bottom") = { 20 };
Physical Surface(".side") = { 11, 15, 19 };
