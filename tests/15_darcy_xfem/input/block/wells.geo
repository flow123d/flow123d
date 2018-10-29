/* triangle inside cube
*/
mesh1 = 5;

h1 = -1;
h2 = 4;

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

Mesh 1;

Save "wells.msh";
Exit;
