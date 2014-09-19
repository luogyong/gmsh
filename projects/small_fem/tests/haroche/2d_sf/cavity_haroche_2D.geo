// Data //
// User
DefineConstant[ F_HAR = {51.099e9, Name "Input/00Haroche/00Frequency"}      ];
DefineConstant[ S_PML = {1,        Name "Input/01Geometry/00PML size"}      ];
DefineConstant[ D_PML = {1,        Name "Input/01Geometry/01PML distance"}  ];
DefineConstant[ MSH_A = {10,       Name "Input/02Mesh/00Size Air"}          ];
DefineConstant[ MSH_P = {10,       Name "Input/02Mesh/00Size PML"}          ];
DefineConstant[ MSH_M = {10,       Name "Input/02Mesh/00Size Mirror"}       ];
DefineConstant[ ORDER = {1,        Name "Input/02Mesh/01Order"}             ];
DefineConstant[ TYPE  = {3,        Name "Input/02Mesh/02Type",
                                   Choices {3="Triangle", 4="Quadrangle"} } ];
// Scaling
nm       = 1e-9;
mm       = 1e-3;
epsilon0 = 8.854187817e-3 * nm;
mu0      = 400* Pi * nm;
cel      = 1 / (Sqrt[epsilon0 * mu0]);

// Haroche Wavelength //
DefineConstant[ lambda_haroche = {cel / F_HAR,
                                  Name "Input/00Haroche/01Wavelength",
                                  ReadOnly 1} ];

// Geomtrical Parameters //
// Mirror
R_small             = 39.4   * mm;
R_big               = 40.6   * mm;
R                   = R_big;
L_cav               = 27.57  * mm;
thick_mirror_center =  1.415 * mm;
radius_mirror       = 25     * mm;

// PML
dist2PML_x = D_PML * lambda_haroche;
dist2PML_y = D_PML * lambda_haroche;
pml_x      = S_PML * lambda_haroche;
pml_y      = S_PML * lambda_haroche;

// Rest
box_x = radius_mirror + dist2PML_x;
box_y = L_cav / 2 + thick_mirror_center + dist2PML_y;
apert = Sqrt[R^2 - radius_mirror^2] + L_cav / 2 - R;

// Mesh
paramaille_air  = lambda_haroche / MSH_A;
paramaille_pml  = lambda_haroche / MSH_P;
paramaille_mir  = lambda_haroche / MSH_M;

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
Physical Line(101)     = {1, 2, 3};      // Neumann OY
Physical Line(102)     = {9, 10};        // Neumann OX
Physical Line(103)     = {16};//, 4, 11};    // Dirichel
Physical Line(104)     = {7, 8, 14, 15}; // Ext PML
Physical Surface(1000) = {20};           // PML_X
Physical Surface(2000) = {22};           // PML_XY
Physical Surface(3000) = {24};           // PML_Y
Physical Surface(4000) = {18};           // Air

// Display
BoundingBox {0, box_x + pml_x, 0, box_y + pml_y, 0, 0};

// Options
If(TYPE == 4)
  Recombine Surface "*";
EndIf

If(ORDER > 1)
  Mesh.HighOrderOptimize = 1;
EndIf

Mesh.Optimize     = 1;
Mesh.ElementOrder = ORDER;
