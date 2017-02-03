/*
  Multiple fractures vs cube.
*/
mesh = 0.16;
mesh_t = 0.12;
Point(1) = {0, 0, 0, mesh};
Point(2) = {1, 0, 0, mesh};
Point(3) = {1, 1, 0, mesh};
Point(4) = {0, 1, 0, mesh};
Point(5) = {0, 1, 1, mesh};
Point(6) = {1, 1, 1, mesh};
Point(7) = {0, 0, 1, mesh};
Point(8) = {1, 0, 1, mesh};
Point(9) = {0.5, 0.5, 0, mesh};
Point(10) = {0.5, 0.5, 1, mesh};
Line(1) = {5, 4};
Line(2) = {7, 1};
Line(3) = {1, 4};
Line(4) = {5, 7};
Line(5) = {7, 8};
Line(6) = {8, 6};
Line(7) = {6, 5};
Line(8) = {6, 3};
Line(9) = {3, 4};
Line(10) = {1, 2};
Line(11) = {2, 8};
Line(12) = {2, 3};
Line(13) = {5, 8};
Line(14) = {4, 2};
Line(15) = {6, 7};
Line(16) = {3, 1};
Line(17) = {9, 10};
Line Loop(18) = {13, -11, -14, -1};
Plane Surface(19) = {18};
Line Loop(20) = {15, 2, -16, -8};
Plane Surface(21) = {20};
Line Loop(22) = {1, -3, -2, -4};
Plane Surface(23) = {22};
Line Loop(24) = {5, 6, 7, 4};
Plane Surface(25) = {24};
Line Loop(26) = {10, 11, -5, 2};
Plane Surface(27) = {26};
Line Loop(28) = {6, 8, -12, 11};
Plane Surface(29) = {28};
Line Loop(30) = {9, -1, -7, 8};
Plane Surface(31) = {30};
Line Loop(32) = {9, -3, 10, 12};
Plane Surface(33) = {32};
Surface Loop(34) = {23, 33, 31, 25, 27, 29};
Volume(35) = {34};


// duplicate parallel down
p0[] = Translate {-0.25, -0.25, 0} {
  Duplicata { Surface{19}; }
};

// duplicate parallel up
p1[] = Translate {0.25, 0.25, 0} {
  Duplicata { Surface{19}; }
};

// triangle fracture
Point(32) = {1.5, 1.4, 0, mesh_t};
Point(33) = {-0.5, 0.3, -0.1, mesh_t};
Point(34) = {-0.1, -0.8, 2.3, mesh_t};
Line(53) = {33, 32};
Line(54) = {32, 34};
Line(55) = {34, 33};
Line Loop(56) = {54, 55, 53};
Plane Surface(57) = {56};


// duplicate parallel up
p2[] = Translate {0.4, 0.3, 1.3} {Rotate {{0, 1, 0}, {0, 0, 0}, 3*Pi/3} {
  Duplicata { Surface{57}; }
}};


Physical Surface("2Dfracture1") = {19, p0[0],p1[0]};
Physical Surface("2Dfracture2") = {21};
Physical Surface("triangles") = {57, p2[0]};
Physical Volume("volume") = {35};
Physical Line("channel") = {17};


Mesh 2;
Mesh 3;

Geometry.Tolerance = 1e-9;
Coherence Mesh;

Save "test1_incomp_mult.msh";
Exit



