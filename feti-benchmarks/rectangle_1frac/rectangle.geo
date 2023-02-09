L = 6;
H = 2;
nx = 2;
ny = 1;

Point(0) = { 0, 0, 0, H/ny };
Point(1) = { 0, H, 0, H/ny };

Line(1) = { 0,1 };
Extrude{ L/2, 0, 0 }{ Line{ 1 }; Layers{ nx/2 }; } // left half domain
Extrude{ L/2, 0, 0 }{ Line{ 2 }; Layers{ nx/2 }; } // right half domain

Physical Surface("rock") = { 5, 9 };
Physical Line("fracture") = { 2 };
Physical Line(".left") = { 1 };
Physical Line(".right") = { 6 };
Physical Line(".bottom") = { 3, 7 };
Physical Line(".top") = { 4, 8 };
Physical Point(".f_bottom") = { 2 };
Physical Point(".f_top") = { 3 };
