pi = 3.141592653589793;
prevod = (2 * pi) / 360.0;

rozsah_x    = 1.50;
rozsah_y    = 1.0;
rozevreni_x = 0.01;

uhel        = 45.0;
delka_pukliny_x = Cos(uhel * prevod)/Sin(uhel * prevod) * rozsah_y;
odstup_y    = 0.1;
odstup_x    = (odstup_y/rozsah_y) * delka_pukliny_x;

diskretizace = 0.02;

Point(1) = {0,0,0,diskretizace};
Point(2) = {rozsah_x - rozevreni_x,0,0,diskretizace};
Point(3) = {rozsah_x - rozevreni_x,rozsah_y,0,diskretizace};
Point(4) = {0,rozsah_y,0,diskretizace};

Point(5) = {0+(rozsah_x - delka_pukliny_x)/2.0,	0,	0,	diskretizace};
Point(7) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0,	rozsah_y,	0,	diskretizace};
Point(9) = {0+(rozsah_x - delka_pukliny_x)/2.0 + odstup_x,	odstup_y,	0,	diskretizace};
Point(11) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0 - odstup_x,	rozsah_y - odstup_y,	0,	diskretizace};

Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 7};
Line(5) = {7, 4};
Line(6) = {4, 1};
Line(7) = {5, 9};
Line(8) = {9, 11};
Line(9) = {11, 7};

Line Loop(10) = {2, 3, 4, -9, -8, -7};
Line Loop(11) = {1, 7, 8, 9, 5, 6};
Plane Surface(12) = {11};
Plane Surface(13) = {10};

// 2D blok
Physical Surface("2d_porous_block_left") = {12};
Physical Surface("2d_porous_block_right") = {13};
// 1D puklina
Physical Line("1d_fracture_ends") = {7, 9};
Physical Line("1d_fracture") = {8};


// dole vlevo puklina vpravo
Physical Line(".2d_bottom_left") = {1};
Physical Point(".1d_bottom_fracture") = {5};
Physical Line(".2d_bottom_right") = {2};
// nahore vlevo puklina vpravo
Physical Line(".2d_top_left") = {5};
Physical Point(".1d_top_fracture") = {7};
Physical Line(".2d_top_right") = {4};



