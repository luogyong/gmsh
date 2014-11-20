// Include //
Include "corner3d_data.geo";

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

// Extrusion //
Extrude {0, 0, LZ} {
  Surface{1:4};
}

// Mesh //

Transfinite Surface "*";
If(HEX == 1)
  Recombine Surface "*";
EndIf
Transfinite Volume  "*";

// Physicals //
If(StrCmp(OnelabAction, "check")) // only mesh if not in onelab check mode
  Mesh 3;
  CreateDir Str(DIR);

  // Vol 1
  Delete Physicals;

  Physical  Volume( 100) = {1};         // Omega 1
  Physical Surface(1000) = {33};        // Gamma Source
  Physical Surface( 200) = {1, 21, 34}; // Gamma Wall
  Physical Surface(3000) = {29};        // Sigma 14
  Physical Surface(4000) = {25};        // Sigma 12
  Physical    Line(  10) = {24};        // Kappa

  Printf("Meshing corner3d subdomain %g...", 0);
  Save StrCat(MSH_NAME, Sprintf("%g.msh", 0));

  // Vol 2
  Delete Physicals;

  Physical  Volume( 101) = {2};         // Omega 2
  Physical Surface(2001) = {47};        // Gamma Neuman
  Physical Surface( 201) = {56, 43, 2}; // Gamma Wall
  Physical Surface(3001) = {25};        // Sigma 21
  Physical Surface(4001) = {51};        // Sigma 23
  Physical    Line(  11) = {24};        // Kappa

  Printf("Meshing corner3d subdomain %g...", 1);
  Save StrCat(MSH_NAME, Sprintf("%g.msh", 1));

  // Vol 3
  Delete Physicals;

  Physical  Volume( 102) = {3};         // Omega 3
  Physical Surface(2002) = {65};        // Gamma Neuman
  Physical Surface( 202) = {78, 69, 3}; // Gamma Wall
  Physical Surface(3002) = {51};        // Sigma 32
  Physical Surface(4002) = {73};        // Sigma 34
  Physical    Line(  12) = {24};        // Kappa

  Printf("Meshing corner3d subdomain %g...", 2);
  Save StrCat(MSH_NAME, Sprintf("%g.msh", 2));

  // Vol 4
  Delete Physicals;

  Physical  Volume( 103) = {4};          // Omega 4
  Physical Surface(1003) = {95};         // Gamma Source
  Physical Surface( 203) = {100, 91, 4}; // Gamma Wall
  Physical Surface(3003) = {73};         // Sigma 43
  Physical Surface(4003) = {29};         // Sigma 41
  Physical    Line(  13) = {24};         // Kappa

  Printf("Meshing corner3d subdomain %g...", 3);
  Save StrCat(MSH_NAME, Sprintf("%g.msh", 3));
EndIf
