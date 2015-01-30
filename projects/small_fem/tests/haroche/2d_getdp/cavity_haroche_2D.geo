// Data //
Include "cavity_haroche_2D.dat";

// Geo //
// Points
Point(1)  = {0,             L_cav / 2 - R,                   0, paramaille_air};
Point(2)  = {0,             0,                               0, paramaille_air};

Point(3)  = {0,             L_cav / 2,                       0, paramaille_mir};
Point(4)  = {radius_mirror, apert,                           0, paramaille_mir};

Point(5)  = {radius_mirror, L_cav / 2 + thick_mirror_center, 0, paramaille_mir};
Point(6)  = {0,             L_cav / 2 + thick_mirror_center, 0, paramaille_mir};
Point(7)  = {0,             box_y,                           0, paramaille_air};

Point(8)  = {box_x,         box_y,                           0, paramaille_air};
Point(9)  = {box_x + pml_x, box_y,                           0, paramaille_pml};
Point(10) = {0,             box_y + pml_y,                   0, paramaille_pml};
Point(11) = {box_x,         box_y + pml_y,                   0, paramaille_pml};
Point(12) = {box_x + pml_x, box_y + pml_y,                   0, paramaille_pml};
Point(13) = {box_x,         0,                               0, paramaille_air};
Point(14) = {box_x + pml_x, 0,                               0, paramaille_pml};

// Lines
Line(1)  = { 2,  3};
Line(2)  = { 6,  7};
Line(3)  = { 7, 10};
Line(4)  = { 4,  5};
Line(5)  = {13,  8};
Line(6)  = { 8, 11};
Line(7)  = {14,  9};
Line(8)  = { 9, 12};
Line(9)  = { 2, 13};
Line(10) = {13, 14};
Line(11) = { 6,  5};
Line(12) = { 7,  8};
Line(13) = { 8,  9};
Line(14) = {10, 11};
Line(15) = {11, 12};

Circle(16) = {3, 1, 4};

// Surfaces
Line  Loop(17)    = {1, 16, 4, -11, 2, 12, -5, -9};
Plane Surface(18) = {17};
Line  Loop(19)    = {5, 13, -7, -10};
Plane Surface(20) = {19};
Line  Loop(21)    = {6, 15, -8, -13};
Plane Surface(22) = {21};
Line  Loop(23)    = {3, 14, -6, -12};
Plane Surface(24) = {23};

// Physicals
Physical Line(101)     = {1, 2, 3};       // Neumann OY
Physical Line(102)     = {9, 10};         // Neumann OX
Physical Line(103)     = {16};//, 4, 11}; // Dirichlet
Physical Line(104)     = {7, 8, 14, 15};  // Ext PML
Physical Surface(1000) = {20};            // PML_X
Physical Surface(2000) = {22};            // PML_XY
Physical Surface(3000) = {24};            // PML_Y
Physical Surface(4000) = {18};            // Air
Physical Point(10000)  = {2};             // Print point

// Display
BoundingBox {0, box_x + pml_x, 0, box_y + pml_y, 0, 0};

// Options
Mesh.Optimize = 1;
