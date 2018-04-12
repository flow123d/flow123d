/* triangle inside cube
*/
mesh1 = 1.4;
mesh2 = 0.2;    // serie: 0.8, 0.6, 0.4, 0.3, 0.2
nref = 1;

h1 = 4;
h2 = 8;

// well 1
Point(1) = {4.1, 4.3, h1, mesh1};
Point(2) = {4.1, 4.3, h2, mesh1};

// well 4
Point(7) = {5.7, 5.9, h1, mesh1};
Point(8) = {5.7, 5.9, h2, mesh1};

Line(101) = {1, 2};
Line(102) = {7, 8};


Physical Line("well_1") = {101};
Physical Line("well_2") = {102};

Physical Point(".well_1_top") = {2};
Physical Point(".well_2_top") = {8};

Physical Point(".well_1_bot") = {1};
Physical Point(".well_2_bot") = {7};



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


Save "single_aquifer_12d_2w.msh";
Exit;
