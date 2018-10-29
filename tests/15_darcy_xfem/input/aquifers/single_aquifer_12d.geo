/* triangle inside cube
*/
mesh1 = 1.4;
mesh2 = 0.8;
nref = 1;

h1 = 4;
h2 = 8;

// well 1
Point(1) = {2.8, 2.5, h1, mesh1};
Point(2) = {2.8, 2.5, h2, mesh1};

// well 2
Point(3) = {4.9, 5.4, h1, mesh1};
Point(4) = {4.9, 5.4, h2, mesh1};

// well 3
Point(5) = {2.9, 7.4, h1, mesh1};
Point(6) = {2.9, 7.4, h2, mesh1};

// well 4
Point(7) = {7.3, 7.8, h1, mesh1};
Point(8) = {7.3, 7.8, h2, mesh1};

// well 5
Point(9) =  {7.4, 2.8, h1, mesh1};
Point(10) = {7.4, 2.8, h2, mesh1};

Line(101) = {1, 2};
Line(102) = {3, 4};
Line(103) = {5, 6};
Line(104) = {7, 8};
Line(105) = {9, 10};


Physical Line("well_1") = {101};
Physical Line("well_2") = {102};
Physical Line("well_3") = {103};
Physical Line("well_4") = {104};
Physical Line("well_5") = {105};

Physical Point(".well_1_top") = {2};
Physical Point(".well_2_top") = {4};
Physical Point(".well_3_top") = {6};
Physical Point(".well_4_top") = {8};
Physical Point(".well_5_top") = {10};

Physical Point(".well_1_bot") = {1};
Physical Point(".well_2_bot") = {3};
Physical Point(".well_3_bot") = {5};
Physical Point(".well_4_bot") = {7};
Physical Point(".well_5_bot") = {9};



zz = 5;
a = 200;

// aquifer 1
Point(a+1) = {0, 0,   zz, mesh2};
Point(a+2) = {10, 0,  zz, mesh2};
Point(a+3) = {0, 10,  zz, mesh2};
Point(a+4) = {10, 10, zz, mesh2};

Line(a+1) = {a+1, a+2};
Line(a+2) = {a+2, a+4};
Line(a+3) = {a+4, a+3};
Line(a+4) = {a+3, a+1};

Line Loop(a+17) = {a+1, a+2, a+3, a+4};
Plane Surface(a+17) = {a+17};

Physical Surface("aquifer") = {a+17};
Physical Line(".aquifer") = {a+1,a+2,a+3,a+4};


Mesh 1;
Mesh 2;


Save "single_aquifer_12d.msh";
Exit;
