// User constant //
DIM = 2;
Include "guide_data.geo";

// Geometry //
// Points
For n In {1:NDOM+1}
  Point(n)            = {+LX / NDOM * (n - 1),  0, 0.5, LC};
  Point(n + NDOM + 1) = {+LX / NDOM * (n - 1), LY, 0.5, LC};
EndFor

// Lines
For n In {1:NDOM+1}
  Line(n) = {n, n + NDOM + 1};
EndFor

For n In {1:NDOM}
  Line(n + (NDOM + 1)       ) = {n           , n        + 1};
  Line(n + (NDOM + 1) + NDOM) = {n + NDOM + 1, n + NDOM + 2};
EndFor

// Faces
For n In {1:NDOM}
  Line Loop(n) = {n, n + (NDOM + 1) + NDOM, -(n + 1), -(n + (NDOM + 1))};
  Plane Surface(n) = {n};
EndFor

// Mesh //
If(STRUCT == 1)
  Transfinite Surface "*";
EndIf

// Physicals //
For n In {1:NDOM}
  Physical Surface(n) = {n};
EndFor

For n In {1:NDOM+1}
  Physical Line(n + NDOM) = {n};
EndFor

Physical Line(2 * NDOM + 2) = {1 + (NDOM + 1):2 * (NDOM + (NDOM + 1))};

If(StrCmp(OnelabAction, "check")) // only mesh if not in onelab check mode
  Printf("Meshing waveguide full...");
  Mesh 3 ;
  CreateDir Str(DIR);
  Save StrCat(MSH_NAME, "all.msh");
  Save "guide2d.msh";
  Printf("Done.");
EndIf
