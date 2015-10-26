// Data //
LC = 0.1;
LX = 1;
LY = 1;

// Geometry (2D) //
// Points
Point(1) = {LX * 0.0, LY * 0.0, 0, LC};
Point(2) = {LX * 0.5, LY * 0.0, 0, LC};
Point(3) = {LX * 1.0, LY * 0.0, 0, LC};

Point(4) = {LX * 1.0, LY * 0.5, 0, LC};

Point(5) = {LX * 1.0, LY * 1.0, 0, LC};
Point(6) = {LX * 0.5, LY * 1.0, 0, LC};
Point(7) = {LX * 0.0, LY * 1.0, 0, LC};

Point(8) = {LX * 0.0, LY * 0.5, 0, LC};
Point(9) = {LX * 0.5, LY * 0.5, 0, LC};

// Lines
Line(1)  = {1, 2};
Line(2)  = {2, 3};

Line(3)  = {8, 9};
Line(4)  = {9, 4};

Line(5)  = {7, 6};
Line(6)  = {6, 5};

Line(7)  = {1, 8};
Line(8)  = {8, 7};
Line(9)  = {2, 9};
Line(10) = {9, 6};
Line(11) = {3, 4};
Line(12) = {4, 5};

// Faces
Line Loop(1) = { 1,  9,  -3, -7};
Line Loop(2) = { 2, 11,  -4, -9};
Line Loop(3) = {12, -6, -10,  4};
Line Loop(4) = {10, -5,  -8,  3};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

// Physical
Physical Surface(11) = {1};     // Volume 1
Physical    Line(21) = {7};     // Src 1
Physical    Line(31) = {1};     // Zero 1
Physical    Line(41) = {9, 3};  // DDM 1

Physical Surface(12) = {2};     // Volume 2
Physical    Line(22) = {11};    // Infinity 2
Physical    Line(32) = {2};     // Zero 2
Physical    Line(42) = {9, 4};  // DDM 2

Physical Surface(13) = {3};     // Volume 3
Physical    Line(23) = {12};    // Infinity 3
Physical    Line(33) = {6};     // Zero 3
Physical    Line(43) = {4, 10}; // DDM 3

Physical Surface(14) = {4};     // Volume 4
Physical    Line(24) = {8};     // Src 4
Physical    Line(34) = {5};     // Zero 4
Physical    Line(44) = {10, 3}; // DDM 4
