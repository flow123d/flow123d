/* triangle inside cube
*/
mesh1 = 2.7;
mesh2 = 3.8;
// aquifer 1
Point(1) = {0, 0, 0, mesh2};
Point(2) = {10, 0, 0, mesh2};
Point(3) = {0, 10, 0, mesh2};
Point(4) = {10, 10, 0, mesh2};

// aquifer 2
Point(5) = {0, 0, 10, mesh2};
Point(6) = {10, 0, 10, mesh2};
Point(7) = {0, 10, 10, mesh2};
Point(8) = {10, 10, 10, mesh2};

// aquifer 3
Point(9) = {0, 0, 4, mesh2};
Point(10) = {10, 0, 6, mesh2};
Point(11) = {0, 10, 4, mesh2};
Point(12) = {10, 10, 6, mesh2};

// well 1
Point(13) = {2.4, 4.4, 12, mesh1};
Point(14) = {2.4, 4.4, -2, mesh1};

// well 2
Point(15) = {2.2, 8.4, 12, mesh1};
Point(16) = {8.2, 8.4, -2, mesh1};

// well 3
Point(17) = {8.4, 4.4, 12, mesh1};
Point(18) = {8.4, 4.4, -2, mesh1};


Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 5};

Line(9) = {9, 10};
Line(10) = {10, 12};
Line(11) = {12, 11};
Line(12) = {11, 9};

// wells
Line(13) = {13, 14};
Line(14) = {15, 16};
Line(15) = {17, 18};

Line Loop(17) = {1, 2, 3, 4};
Plane Surface(17) = {17};
Line Loop(19) = {5, 6, 7, 8};
Plane Surface(19) = {19};
Line Loop(21) = {9, 10, 11, 12};
Plane Surface(21) = {21};

Physical Surface("aquifer1") = {17};
Physical Surface("aquifer2") = {19};
Physical Surface("aquifer3") = {21};

Physical Line("well1") = {13};
Physical Line("well2") = {14};
Physical Line("well3") = {15};


Mesh 1;
Mesh 2;

Geometry.Tolerance = 1e-9;
Coherence Mesh;

Save "prolongation_12d_01.msh";
Exit;