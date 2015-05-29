// User constant //
Include "waveguide.dat";

// Geometry (2D) //
// Points
Point(1) = { 0,  0, 0, LC};
Point(2) = {LX,  0, 0, LC};
Point(3) = {LX, LY, 0, LC};
Point(4) = { 0, LY, 0, LC};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Faces
Line  Loop(1)    = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Extrusion //
ext[] = Extrude{0, 0, LZ}{Surface{1};};

// Physicals //
Physical  Volume(1) = {ext[1]};                    // Volume
Physical Surface(2) = {1, ext[0], ext[2], ext[4]}; // Walls
Physical Surface(3) = {ext[5]};                    // Source
Physical Surface(4) = {ext[3]};                    // Infinity
