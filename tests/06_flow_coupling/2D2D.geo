// domain dimensions
domain_x    = 1.50;
domain_y    = 1.0;

// We consider fracture placed in the middle of domain,
// intersecting only horizontal boundary. The fracture central
// part may hav different material properties, effectively model
// fracture with ends inside

// fracture width
delta=0.02;

// fracture length relative to domain_y
main_rel_length=1.0;

// fracture angle  with y-axis (in degrees)
// we assume that fracture intersects only horizontal boundary
angle=45;

// mesh step on fracture
h_fracture=0.005;

// mesh step elsewhere
h_domain=0.005;

////////////////////////////////////////////////////////////////////////////////////////////////////
pi = 3.141592653589793;
deg_to_rad = (2 * pi) / 360.0;
angle_rad = angle * deg_to_rad;
tan = Sin(angle_rad) / Cos(angle_rad);

// length of intersection with x-axis
delta_x = delta/ Cos(angle_rad);

// coordinates of left border of the fracture
// fracture is placed in the middle
fracture_bottom_x=(domain_x - tan * domain_y - delta_x)/2;
fracture_top_x=(domain_x + tan * domain_y - delta_x)/2;


// corner points
Point(1) = {0,          0,              0,h_domain};
Point(2) = {domain_x,   0,              0,h_domain};
Point(3) = {domain_x,   domain_y,       0,h_domain};
Point(4) = {0,          domain_y,       0,h_domain};


// fracture corners (intersect with boundary)
Point(5) = {fracture_bottom_x, 0, 0, h_fracture};
Point(8) = {fracture_top_x, domain_y, 0, h_fracture};

Point(9) = {fracture_bottom_x+delta_x, 0, 0, h_fracture};
Point(12) = {fracture_top_x+delta_x, domain_y, 0, h_fracture};

//left
Line(1) = {1, 5};
Line(2) = {5, 8};
Line(3) = {8, 4};
Line(4) = {4, 1};
Line Loop(17) = { 1, 2, 3, 4};
Plane Surface(18) = {17};
Physical Surface("matrix_left") = {18};

//right
Line(5) = {3, 12};
Line(6) = {12, 9};
Line(7) = {9, 2};
Line(8) = {2, 3};
Line Loop(19) = { 5, 6, 7, 8};
Plane Surface(20) = {19};
Physical Surface("matrix_right") = {20};

// fracture
Line(10) = {8, 12};
Line(12) = {9, 5};
Line Loop(21) = { 2, 10, 6, 12};
Plane Surface(22) = {21};
Physical Surface("fracture") = {22};

// boundary
Physical Line(".top") = {3, 5};
Physical Line(".bottom") = {1, 7};
Physical Line(".left") = {4};
Physical Line(".right") = {8};

Physical Line(".top_fracture") = {10};
Physical Line(".bottom_fracture") = {12};

