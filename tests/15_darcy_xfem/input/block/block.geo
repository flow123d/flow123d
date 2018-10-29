/* circle aquifer
*/
mesh2 = 1;

zz = 0;
height = 3;
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


surfaceVector[] = Extrude {0, 0, height} {
  Surface{a+17};
//   Layers{3};
  Recombine;
};


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
Physical Surface(".bottom") = {a+17}; // from Plane Surface (21) ...
    
//Physical Surface("aquifer") = {21};
//Physical Line(".aquifer") = {1,2};
  
Mesh 3;

// Geometry.Tolerance = 1e-9;
// Coherence Mesh;
// 
Save "block.msh";
Exit;
