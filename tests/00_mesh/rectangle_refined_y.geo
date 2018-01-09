lx = 1;
lyf = 0.00001;
ly  = 0.1;

clx = 0.1;
cly = 0.0001;


Point(1) = { 0, 0, 0, cly };
Point(2) = { 0, lyf, 0, cly };
Point(3) = { 0, ly, 0, clx };
Line(1) = { 1,2 };
Line(2) = { 2,3 };

Extrude{ lx, 0, 0 }{ Line{ 1,2 }; Layers{ lx/clx }; }

Physical Surface("rock") = { 6,10 };
Physical Line("fracture") = { 4 };
Physical Point(".left") = { 1 };
Physical Point(".right") = { 4 };
