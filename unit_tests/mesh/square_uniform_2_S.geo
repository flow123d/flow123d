mesh = 0.0294;   # 3k elem
# mesh = 0.00895;  # 30k elem
# mesh = 0.00281;  # 300k elem

Point(1) = {0, 0, 0, mesh};
Point(2) = {1, 0, 0, mesh};
Point(3) = {1, 1, 0, mesh};
Point(4) = {0, 1, 0, mesh};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 3};

Line Loop(7) = {4, 1, 2, 3};
Plane Surface(7) = {7};

Line { 5 } In Surface { 7 };

Physical Line(".boundary") = {1, 2, 3, 4};
Physical Line("fracture") = {5};
Physical Surface("plane") = {7};
