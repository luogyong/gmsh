// User constant //
DefineConstant[
  X    = {2,         Name "Geometry/00X dimension"},
  Y    = {1,         Name "Geometry/01Y dimension"},
  Z    = {1,         Name "Geometry/02Z dimension"},
  NDOM = {2,         Name "Geometry/03Number of subdomain"},

  K    = {10,        Name "Geometry/04Wavenumber"},
  L    = {2* Pi / K, Name "Geometry/05Wavelength", ReadOnly 1},

  NpL  = {10,        Name "Geometry/06Mesh points per wavelength"}
];

// Number of wavelength
NLX  = Round(X/L);
NLY  = Round(Y/L);
NLZ  = Round(Z/L);

// Number of mesh element
NpLX = NpL * NLX;
NpLY = NpL * NLY;
NpLZ = NpL * NLZ;

// Number of mesh element per subdomain
NY   = NpLY;
NZ   = NpLZ;
NX   = Round(NpLX / NDOM) + 2;

// Geometry (2D) //
// Points
For n In {1:NDOM+1}
  Point(n)            = {+X / NDOM * (n - 1), 0, 0, 1};
  Point(n + NDOM + 1) = {+X / NDOM * (n - 1), Y, 0, 1};
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
For n In {1:NDOM+1}
  Transfinite Line {n} = NY Using Progression 1;
EndFor

For n In {1:NDOM}
  Transfinite Line(n + (NDOM + 1)       ) = NX Using Progression 1;
  Transfinite Line(n + (NDOM + 1) + NDOM) = NX Using Progression 1;
EndFor

For n In {1:NDOM}
  Transfinite Surface {n};
EndFor

// Extrusion //
ext[] = Extrude {0, 0, Z} {
  Surface{1:NDOM};
  Layers{NZ};
};

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
