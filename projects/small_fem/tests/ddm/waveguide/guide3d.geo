// User constant //
DIM = 3;
Include "guide_data.geo";

// Geometry (2D) //
// Points
For n In {1:NDOM+1}
  Point(n)            = {+LX / NDOM * (n - 1),  0, 0, LC};
  Point(n + NDOM + 1) = {+LX / NDOM * (n - 1), LY, 0, LC};
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

// Extrusion //
ext[] = Extrude {0, 0, LZ} {
  Surface{1:NDOM};
};

// Mesh //
If(STRUCT == 1)
  Transfinite Surface "*";
  Transfinite Volume  "*";
EndIf

// Physicals //
zero[] = {};
For n In {1:NDOM}
  zero[n-1] = n;
EndFor

For n In {0:NDOM-1}
  zero[0 + n * 3 + NDOM] = ext[0 + n * 6];
  zero[1 + n * 3 + NDOM] = ext[3 + n * 6];
  zero[2 + n * 3 + NDOM] = ext[5 + n * 6];

  Physical Volume (n + 1)        = {ext[1 + n * 6]};
  Physical Surface(n + 1 + NDOM) = {ext[2 + n * 6]};
EndFor

Physical Surface(2 * NDOM + 1) = {ext[4 + (NDOM - 1) * 6]};
Physical Surface(2 * NDOM + 2) = {zero[]};

If(StrCmp(OnelabAction, "check")) // only mesh if not in onelab check mode
  Printf("Meshing waveguide full...");
  Mesh 3 ;
  CreateDir Str(DIR);
  Save StrCat(MSH_NAME, "all.msh");
  Save "guide3d.msh";
  Printf("Done.");
EndIf
