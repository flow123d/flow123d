/* circle aquifer
*/
// mesh2 = 0.8;
mesh2 = 1;
R = 5;
ctrX = 3.33;
ctrY = 3.33;
// aquifer
Point(1) = {ctrX, ctrY, 0, mesh2};
Point(2) = {ctrX+R, ctrY, 0, mesh2};
Point(3) = {ctrX-R, ctrY, 0, mesh2};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

Line Loop(20) = {1,2} ;
Plane Surface(20) = {20};

Physical Surface("aquifer") = {20};
Physical Line(".aquifer") = {1,2};

Mesh 2;

Geometry.Tolerance = 1e-9;
Coherence Mesh;

Save "big_circle.msh";
Exit;
