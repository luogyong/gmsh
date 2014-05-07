r  = 0.33;
R1 = 0.66;
R2 = 0.99;
lc = 0.1;

Point(0) = {0, 0, 0, lc};

Point(11) = {+r,  0, 0, lc};
Point(12) = { 0, +r, 0, lc};
Point(13) = {-r,  0, 0, lc};
Point(14) = { 0, -r, 0, lc};

Point(21) = {+R1,  0, 0, lc};
Point(22) = { 0, +R1, 0, lc};
Point(23) = {-R1,  0, 0, lc};
Point(24) = { 0, -R1, 0, lc};

Point(31) = {+R2,  0, 0, lc};
Point(32) = { 0, +R2, 0, lc};
Point(33) = {-R2,  0, 0, lc};
Point(34) = { 0, -R2, 0, lc};

Circle(11) = {11, 0, 12};
Circle(12) = {12, 0, 13};
Circle(13) = {13, 0, 14};
Circle(14) = {14, 0, 11};

Circle(21) = {21, 0, 22};
Circle(22) = {22, 0, 23};
Circle(23) = {23, 0, 24};
Circle(24) = {24, 0, 21};

Circle(31) = {31, 0, 32};
Circle(32) = {32, 0, 33};
Circle(33) = {33, 0, 34};
Circle(34) = {34, 0, 31};

Line Loop(41) = {11, 12, 13, 14};
Line Loop(42) = {21, 22, 23, 24};
Line Loop(43) = {31, 32, 33, 34};

Plane Surface(51) = {42, 41};
Plane Surface(52) = {43, 42};

// --- //
Physical Line(4) = {11, 12, 13, 14}; // Source
Physical Line(5) = {21, 22, 23, 24}; // DDM border
Physical Line(6) = {31, 32, 33, 34}; // Infinity

Physical Surface(7) = {51}; // Domain 1
Physical Surface(8) = {52}; // Domain 2
