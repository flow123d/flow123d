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
Point(2) = {rozsah_x,0,0,diskretizace};
Point(3) = {rozsah_x,rozsah_y,0,diskretizace};
Point(4) = {0,rozsah_y,0,diskretizace};


Point(5) = {0+(rozsah_x - delka_pukliny_x)/2.0,	0,	0,	diskretizace};
Point(6) = {0+(rozsah_x - delka_pukliny_x)/2.0 + rozevreni_x,	0,	0,	diskretizace};

Point(7) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0,	rozsah_y,	0,	diskretizace};
Point(8) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0 + rozevreni_x,	rozsah_y,	0,	diskretizace};

Point(9) = {0+(rozsah_x - delka_pukliny_x)/2.0 + odstup_x,	odstup_y,	0,	diskretizace};
Point(10) = {0+(rozsah_x - delka_pukliny_x)/2.0 + rozevreni_x +odstup_x,	odstup_y,	0,	diskretizace};

Point(11) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0 - odstup_x,	rozsah_y - odstup_y,	0,	diskretizace};
Point(12) = {rozsah_x-(rozsah_x - delka_pukliny_x)/2.0 + rozevreni_x- odstup_x,	rozsah_y - odstup_y,	0,	diskretizace};


Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 3};
Line(5) = {3, 8};
Line(6) = {8, 7};
Line(7) = {7, 4};
Line(8) = {4, 1};
Line(9) = {9, 10};
Line(10) = {12, 11};
Line(11) = {5, 9};
Line(12) = {9, 11};
Line(13) = {11, 7};
Line(14) = {6, 10};
Line(15) = {10, 12};
Line(16) = {12, 8};


Line Loop(17) = {1, 11, 12, 13, 7, 8};
Plane Surface(18) = {17};
Line Loop(19) = {3, 4, 5, -16, -15, -14};
Plane Surface(20) = {19};
Line Loop(21) = {2, 14, -9, -11};
Plane Surface(22) = {21};
Line Loop(23) = {12, -10, -15, -9};
Plane Surface(24) = {23};
Line Loop(25) = {10, 13, -6, -16};
Plane Surface(26) = {25};

// 2D blok
Physical Surface("2d_porous_block_left") = {18};
Physical Surface("2d_porous_block_right") = {20};
// 1D puklina
Physical Surface("2d_fracture_ends") = {22, 26};
Physical Surface("2d_fracture") = {24};

// dolni okraj vlevo, puklina, vpravo
Physical Line(".2d_bottom_left") = {1};
Physical Line(".2d_bottom_fracture") = {2};
Physical Line(".2d_bottom_right") = {3};
// horni okraje vlevo, puklina, vpravo
Physical Line(".2d_top_left") = {7};
Physical Line(".2d_top_fracture") = {6};
Physical Line(".2d_top_right") = {5};
