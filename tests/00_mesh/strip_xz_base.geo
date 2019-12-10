cl1 = 0.02;
lx=3*cl1;
ztop=0;
zbot=-0.25;
d = 0.025;

Point(1) = {0, 0, zbot, cl1};
Point(2) = {lx, 0, zbot, cl1};
Point(3) = {lx, 0, ztop, cl1};
Point(4) = {0, 0, ztop, cl1};

Point(5) = {lx, 0, zbot-d, cl1};
Point(6) = {0, 0, zbot-d, cl1};

Line(20) = {1,2};
Line(21) = {2,3};
Line(22) = {3,4};
Line(23) = {4,1};

Line(24) = {2,5};
Line(25) = {5,6};
Line(26) = {6,1};

Line Loop(30) = {20, 21, 22, 23};
Plane Surface(30) = {30};

Line Loop(31) = {20, 24, 25, 26};
Plane Surface(31) = {31};

Physical Line(".left") = {23};
Physical Line(".right") = {21};
Physical Line(".top") = {22};
Physical Line(".bottom") = {20};
Physical Surface("plane") = {30};
Physical Surface("base") = {31};

//Mesh 1;
//Mesh 2;

