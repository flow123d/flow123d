h1=0.1; // step on 1D fracture
h2=0.1; // step on 2D domain

// Number of extrusion layers in Z direction
z_layers = 2;

Point(1) = {-1,1,0,h2};
Point(2) = {-1,-1,0,h2};
Point(3) = {1,-1,0,h2};
Point(4) = {1,1,0,h2};

Point(5) = {-1,0,0,h1};
Point(6) = {1,0,0,h1};

Line(1) = {5, 6};
Line(2) = {1, 4};
Line(7) = {2, 3};

Line(8) = {1, 2};
Line(9) = {3, 4};
Line Loop(10) = {8, 7, 9, -2};

Plane Surface(11) = {10};
// , Line{1}, Line{7}, Line{2}, Point{5}, Point{6}
SurfaceList[] = Extrude{ 0, 0, z_layers*h2 } { Surface{11}; Layers{ z_layers }; };
Printf("%f %f %f   %f %f %f", 
       SurfaceList[0], SurfaceList[1], SurfaceList[2], 
       SurfaceList[3], SurfaceList[4], SurfaceList[5]); 

/* SurfaceList contains new geometric entities:
 * SurfaceList[0] = 33 - shifted surface 11 )
 * SurfaceList[1] = 1  - volume
 * SurfaceList[2] = 20 - extruded line 8, 1. in Line Loop
 * SurfaceList[3] = 24 - extruded line 7, 2. in Line Loop
 * SurfaceList[4] = 28 - extruded line 9, 3. in Line Loop
 * SurfaceList[5] = 32 - extruded line 2, 4. in Line Loop
 */   

//LineList[] = {1, 2, 7};
LineList[] = Extrude{ 0, 0, z_layers*h2 } { Line{ 1 }; Layers{ z_layers }; };
/*
 * LineList[0] = 34 - shifted line 1 )
 * LineList[1] = 37  - surface
 * LineList[2] = 35 - extruded point 5
 * LineList[3] = 36 - extruded point 6
 */   


Physical Surface("2d") = {LineList[1]};
Physical Volume("3d") = {SurfaceList[1]};

Physical Surface(".3d_bottom_top") = {SurfaceList[3], SurfaceList[5]};
Physical Line(".2d") = {LineList[2], LineList[3]};

