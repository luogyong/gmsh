// User constant //
DefineConstant[
  X  = {2,           Name "Geometry/00X dimension"},
  Y  = {1,           Name "Geometry/01Y dimension"},
  Z  = {1,           Name "Geometry/02Z dimension"},
  N  = {4,           Name "Geometry/03Number of subdomain"},

  K  = {10,          Name "Geometry/04Wavenumber"},
  LB = {2* Pi / K,   Name "Geometry/05Wavelength", ReadOnly 1},

  NL = {5,          Name "Geometry/06Points per wavelength"},
  NX = {X / LB * NL, Name "Geometry/07Mesh division", ReadOnly 1}
];

NY = NX / 2;
NZ = NX / 2;
CL = LB / NL;

// Geometry (2D) //
// Points
For n In {1:N+1}
  Point(n)         = {+X / N * (n - 1), 0, 0, CL};
  Point(n + N + 1) = {+X / N * (n - 1), Y, 0, CL};
EndFor

// Lines
For n In {1:N+1}
  Line(n) = {n, n + N + 1};
EndFor

For n In {1:N}
  Line(n + (N + 1)    ) = {n        , n     + 1};
  Line(n + (N + 1) + N) = {n + N + 1, n + N + 2};
EndFor

// Faces
For n In {1:N}
Line Loop(n) = {n, n + (N + 1) + N, -(n + 1), -(n + (N + 1))};
  Plane Surface(n) = {n};
EndFor

// Mesh //
For n In {1:N+1}
  Transfinite Line {n} = NY Using Progression 1;
EndFor

For n In {1:N}
  Transfinite Line(n + (N + 1)    ) = NX Using Progression 1;
  Transfinite Line(n + (N + 1) + N) = NX Using Progression 1;
EndFor

For n In {1:N}
  Transfinite Surface {n};
EndFor

// Extrusion //
ext[] = Extrude {0, 0, Z} {
  Surface{1:N};
  Layers{NZ};
};

// Physicals //
zero[] = {};
For n In {1:N}
  zero[n-1] = n;
EndFor

For n In {0:N-1}
  zero[0 + n * 3 + N] = ext[0 + n * 6];
  zero[1 + n * 3 + N] = ext[3 + n * 6];
  zero[2 + n * 3 + N] = ext[5 + n * 6];

  Physical Volume (n + 1)     = {ext[1 + n * 6]};
  Physical Surface(n + 1 + N) = {ext[2 + n * 6]};
EndFor

Physical Surface(2 * N + 1) = {ext[4 + (N - 1) * 6]};
Physical Surface(2 * N + 2) = {zero[]};
