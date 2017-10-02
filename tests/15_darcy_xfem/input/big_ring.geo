/* circle aquifer
*/
// mesh2 = 0.8;
mesh2 = 0.25;
R = 5;
ctrX = 3.33;
ctrY = 3.33;
rho = 0.03;
mesh_rho = (6.28*rho)/50.0;
// aquifer
Point(1) = {ctrX, ctrY, 0, mesh2};
Point(2) = {ctrX+R, ctrY, 0, mesh2};
Point(3) = {ctrX-R, ctrY, 0, mesh2};

Point(4) = {ctrX, ctrY, 0, mesh_rho};
Point(5) = {ctrX+rho, ctrY, 0, mesh_rho};
Point(6) = {ctrX-rho, ctrY, 0, mesh_rho};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

Circle(3) = {5,4,6};
Circle(4) = {6,4,5};

Line(5) = {2,5};
Line(6) = {3,6};

Line Loop(8) = {3, -6, -1, 5};
Plane Surface(8) = {8};
Line Loop(10) = {5, -4, -6, 2};
Plane Surface(10) = {10};


Physical Surface("aquifer") = {8,10};
Physical Line(".aquifer") = {1,2};
Physical Line(".well") = {3,4};

Mesh 2;

Geometry.Tolerance = 1e-9;
Coherence Mesh;

Save "big_ring.msh";
Exit;
