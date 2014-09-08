Lc = 0.1;
x = 1.;
y = 1.;
z = 1.;
n = 5;

Point(1) = { 0. , 0. , 0. , Lc};
Point(2) = { x  , 0. , 0. , Lc};
Point(3) = { x  , y  , 0. , Lc};
Point(4) = { 0. , y  , 0. , Lc};

Point(5) = { 0. , 0. , z , Lc};
Point(6) = { x  , 0. , z , Lc};
Point(7) = { x  , y  , z , Lc};
Point(8) = { 0. , y  , z , Lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {1,5};
Line(10)= {2,6};
Line(11) ={3,7};
Line(12) ={4,8};

Line Loop(1) = {-1,-4,-3,-2};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {1,10,-5,-9};
Line Loop(4) = {2,11,-6,-10};
Line Loop(5) = {3,12,-7,-11};
Line Loop(6) = {4,9,-8,-12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = n+1;
Transfinite Surface{1,2,3,4,5,6};
Transfinite Volume{1};
Recombine Surface{1,2,3,4,5,6};
Recombine Volume{1};

Physical Volume(101) = {1};
Physical Surface(1584) = {6};
Physical Surface(2376) = {4};
