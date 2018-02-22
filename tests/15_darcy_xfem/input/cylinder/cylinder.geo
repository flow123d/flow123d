/* circle aquifer
*/
mesh2 = 1.5;
meshC = 0.05;
R = 5;
ctrX = 3.33;
ctrY = 3.33;
cylZ = 2;
// aquifer
Point(1) = {ctrX, ctrY, 0, meshC};
Point(2) = {ctrX-R, ctrY, 0, mesh2};
Point(3) = {ctrX, ctrY+R, 0, mesh2};
Point(4) = {ctrX+R, ctrY, 0, mesh2};
Point(5) = {ctrX, ctrY-R, 0, mesh2};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(20) = {1,2,3,4} ;
Plane Surface(21) = {20};


surfaceVector[] = Extrude {0, 0, cylZ} {
  Surface{21};
//   Layers{3};
  Recombine;
};


Point(10) = {ctrX, ctrY, 0, meshC};
Point(11) = {ctrX, ctrY, cylZ, meshC};

Field[1] = Attractor;
Field[1].NodesList = {10, 11};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = meshC;
Field[2].LcMax = 16*meshC;
Field[2].DistMin = mesC;
Field[2].DistMax = R;

Background Field = 2;

/* surfaceVector contains in the following order:
[0] - front surface (opposed to source surface)
[1] - extruded volume
[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
[3] - right surface (belonging to 2nd line in "Line Loop (6)")
[4] - top surface (belonging to 3rd line in "Line Loop (6)")
[5] - left surface (belonging to 4th line in "Line Loop (6)") */
Physical Surface(".top") = surfaceVector[0];
Physical Volume("aquifer") = surfaceVector[1];
Physical Surface(".aquifer") = {surfaceVector[2],surfaceVector[3],surfaceVector[4],surfaceVector[5]};
Physical Surface(".bottom") = {21}; // from Plane Surface (21) ...
    
//Physical Surface("aquifer") = {21};
//Physical Line(".aquifer") = {1,2};
  
Mesh 3;

// Geometry.Tolerance = 1e-9;
// Coherence Mesh;
// 
// Save "cylinder.msh";
// Exit;
