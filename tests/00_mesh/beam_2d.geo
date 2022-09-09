cl1 = 0.06;
lx=1;
ly=0.06;

Point(1) = {0, 0, 0, cl1};
Point(2) = {0, ly, 0, cl1};
Line(10) = {1,2};
Extrude{ lx, 0, 0 }{ Line{ 10 }; Layers{ lx/cl1 }; }
Physical Line(".left") = {10};
Physical Line(".right") = {11};
Physical Line(".top") = {13};
Physical Line(".bottom") = {12};
Physical Surface("plane") = {14};
