// refinement parameter
layers = 20;

// dimensions of domain
lx = 1;
ly = 1;


Point(0) = { lx, 0, 0, 1/layers };
Point(1) = { lx, ly, 0, 1/layers };

Line(1) = { 0, 1 };
Extrude{ -lx, 0, 0 }{ Line{1}; Layers{ layers }; }

Physical Surface("domain") = { 5 };
Physical Line(".left") = { 2 };
Physical Line(".right") = { 1 };
Physical Line(".top") = { 4 };
Physical Line(".bottom") = { 3 };

