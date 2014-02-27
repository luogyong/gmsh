r = 1;
R = 10;
D = 8;

d  = 5;
cl = r / d;

Point(0) = { 0,  0, 0, cl};

Point(1) = {+r,  0, 0, cl};
Point(2) = { 0, +r, 0, cl};
Point(3) = {-r,  0, 0, cl};
Point(4) = { 0, -r, 0, cl};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Point(11) = {+R, -R, 0, cl};
Point(12) = {+R, +R, 0, cl};
Point(13) = {-R, +R, 0, cl};
Point(14) = {-R, -R, 0, cl};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};

Point(21) = {+R + D, -R - D, 0, cl};
Point(22) = {+R + D, +R + D, 0, cl};
Point(23) = {-R - D, +R + D, 0, cl};
Point(24) = {-R - D, -R - D, 0, cl};

Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 21};

Line Loop(1) = { 1,  2,  3,  4};
Line Loop(2) = {11, 12, 13, 14};
Line Loop(3) = {22, 23, 24, 21};

Plane Surface(1) = {2, 1};
Plane Surface(2) = {2, 3};

Physical Line(4)    = {21, 22, 23, 24};
Physical Line(5)    = { 1,  2,  3,  4};
Physical Line(6)    = {11, 12, 13, 14};
Physical Surface(7) = {1};
Physical Surface(8) = {2};
