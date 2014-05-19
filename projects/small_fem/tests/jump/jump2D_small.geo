l   = 1;
d   = 2;
L   = l / 2;
cl  = L / d;

Point(1) = { 0,  0, 0, cl};
Point(2) = { 0, +L, 0, cl};
Point(3) = {-L, +L, 0, cl};
Point(4) = {-L,  0, 0, cl};

Point(5) = {+L, +L, 0, cl};
Point(6) = {+L,  0, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 1};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {1, 5, 6, 7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line(5)    = {3};
Physical Line(6)    = {6};
Physical Surface(7) = {1};
Physical Surface(8) = {2};
