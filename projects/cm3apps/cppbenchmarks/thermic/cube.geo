// Test case a SCB with a vertical load at its free extremity
// Size

//definition of unit
mm = 1e-03;

// volum fraction

x=1*mm;
y=1*mm;
z=1*mm;

// Characteristic length
Lc1=z/2.5;

// definition of points
Point(1) = { 0.0 , 0.0 , 0.0 , Lc1};
Point(2) = {  x  , 0.0 , 0.0 , Lc1};
Point(3) = {  x  ,  y  , 0.0 , Lc1};
Point(4) = { 0.0 ,  y  , 0.0 , Lc1};
Point(5) = { 0.0 , 0.0 , z , Lc1};
Point(6) = {  x  , 0.0 , z , Lc1};
Point(7) = {  x  ,  y  , z , Lc1};
Point(8) = { 0.0 ,  y  , z , Lc1};

// Line between points
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
Line(11)= {3,7};
Line(12)= {4,8};

// Surface definition
Line Loop(1) = {1,2,3,4};
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

//VOlume

Surface Loop(7) = {1,2,3,4,5,6};
Volume(1) = {7};

// Physical objects to applied BC and material
Physical Surface(1234) = {1};
Physical Surface(5678) = {2};
Physical Surface(1265) = {3};
Physical Surface(2376) = {4};
Physical Surface(3487) = {5};
Physical Surface(4158) = {6};
Physical Line(12) = {1};
Physical Line(23) = {2};
Physical Line(34) = {3};
Physical Line(41) = {4};
Physical Line(56) = {5};
Physical Line(67) = {6};
Physical Line(78) = {7};
Physical Line(85) = {8};
Physical Line(15) = {9};
Physical Line(26) = {10};
Physical Line(37) = {11};
Physical Line(48) = {12};

Physical Point(1) ={1};
Physical Point(2) ={2};
Physical Point(3) ={3};
Physical Point(4) ={4};
Physical Point(5) ={5};
Physical Point(6) ={6};
Physical Point(7) ={7};
Physical Point(8) ={8};

Physical Volume(10) ={1};

// define transfinite mesh
Transfinite Line {1,2,3,4,5,6,7,8} = 3;
Transfinite Line {9,10,11,12} = 3;
Transfinite Surface {1,2,3,4,5,6} ;
Recombine Surface {1,2,3,4,5,6} ;
Transfinite Volume {1};

