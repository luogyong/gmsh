cl = 1;
l  = 1;

Point(1) = {+l, -l, -l, cl};
Point(2) = {+l, +l, -l, cl};
Point(3) = {-l, +l, -l, cl};
Point(4) = {-l, -l, -l, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1)     = {1, 2, 3, 4};
Plane Surface(1) = {1};

//Recombine Surface {1};

Extrude {0, 0, 2 * l} {
  Surface{1};
}

Physical Surface(5) = {21};
Physical Surface(6) = {13};
Physical Volume(7)  = {1};
