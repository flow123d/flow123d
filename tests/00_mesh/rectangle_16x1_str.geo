lx = 16;
ly = 1;
lc = 0.2;

Point(1) = { 0,0,0,lc };
Point(2) = { 0,ly,0,lc };

Line(1) = { 1,2 };

Extrude{ lx,0,0 }{ Line{1}; Layers{lx/lc}; }

Physical Surface("domain") = { 5 };
Physical Line(".left") = { 1 };
