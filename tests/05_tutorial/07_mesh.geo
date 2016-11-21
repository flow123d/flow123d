l = 5000;
l1 = 1000;
l2 = 4000;

h = -5000;
h1 = -4100;
h2 = -4500;
hv = -4200;

// nastaveni zjemneni site
lc = 100;
lc2 = 0.1;
Field[1] = MathEval;
//Field[1].F = "10*Max(1,Min(Sqrt((x-1000)^2+(y-2500)^2),Sqrt((x-4000)^2+(y-2500)^2)))";
Field[1].F = "Max(0.1,Min(1000,1e4*(Exp((x-1000)^2+Max(0,Abs(z+4150)-50)^2)-1),1e4*(Exp((x-4000)^2+Max(0,Abs(z+4150)-50)^2)-1)))";

Field[2] = BoundaryLayer;
Field[2].hwall_n = 1;
Field[2].hwall_t = 10;
Field[2].thickness = 10;
Field[2].ratio = 10;
Field[2].hfar = 200;
Field[2].EdgesList = { 16, 18 };
Background Field = 2;
Mesh.CharacteristicLengthFromPoints = 0;

Point(1)  = { 0, 0, 0 };
Point(2)  = { l1, 0, 0 };
Point(25) = { l1, 0, -100 };
Point(3)  = { l2, 0, 0 };
Point(35) = { l2, 0, -100 };
Point(4)  = { l, 0, 0 };
Point(5)  = { 0, 0, h1 };
Point(6)  = { l1, 0, h1 };
Point(7)  = { l2, 0, h1 };
Point(8)  = { l, 0, h1 };
Point(9)  = { l1, 0, hv };
Point(10) = { l2, 0, hv };
Point(11) = { 0, 0, h2 };
Point(12) = { l1, 0, h2 };
Point(13) = { l2, 0, h2 };
Point(14) = { l, 0, h2 };
Point(15) = { 0, 0, h };
Point(16) = { l, 0, h };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 5, 6 };
Line(5) = { 6, 7 };
Line(6) = { 7, 8 };
Line(7) = { 11, 12 };
Line(8) = { 12, 13 };
Line(9) = { 13, 14 };
Line(10) = { 15, 16 };
Line(11) = { 1, 5 };
Line(12) = { 2, 25 };
Line(125)= { 25, 6 };
Line(13) = { 3, 35 };
Line(135)= { 35, 7 };
Line(14) = { 4, 8 };
Line(15) = { 5, 11 };
Line(16) = { 6, 9 };
Line(17) = { 9, 12 };
Line(18) = { 7, 10 };
Line(19) = { 10, 13 };
Line(20) = { 8, 14 };
Line(21) = { 11, 15 };
Line(22) = { 14, 16 };

Line Loop(10) = { 1, 12, 125, -4, -11 };
Line Loop(20) = { 2, 13, 135, -5, -125, -12 };
Line Loop(30) = { 3, 14, -6, -135, -13 };
Line Loop(40) = { 4, 16, 17, -7, -15 };
Line Loop(50) = { 5, 18, 19, -8, -17, -16 };
Line Loop(60) = { 6, 20, -9, -19, -18 };
Line Loop(70) = { 7, 8, 9, 22, -10, -21 };

Plane Surface(10) = { 10 };
Plane Surface(20) = { 20 };
Plane Surface(30) = { 30 };
Plane Surface(40) = { 40 };
Plane Surface(50) = { 50 };
Plane Surface(60) = { 60 };
Plane Surface(70) = { 70 };

Physical Surface("near_surface") = { 10, 20, 30 };
Physical Surface("exchanger") = { 40, 50, 60 };
Physical Surface("deep") = { 70 };

Physical Line("well1_surface") = { 12 };
Physical Line("well1_middle")  = { 125 };
Physical Line("well1_deep")   = { 16 };
Physical Line("well2_surface") = { 13 };
Physical Line("well2_middle")  = { 135 };
Physical Line("well2_deep")   = { 18 };


Physical Line(".surface") = { 1, 2, 3 };
Physical Line(".deep") = { 10 };
Physical Line(".sides") = { 11, 15, 21, 14, 20, 22 };
Physical Point(".well1_surface") = { 2 };
Physical Point(".well2_surface") = { 3 };
Physical Point(".well1_exchanger") = { 9 };
Physical Point(".well2_exchanger") = { 10 };




