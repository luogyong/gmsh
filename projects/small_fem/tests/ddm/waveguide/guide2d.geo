// User constant //
DefineConstant[
  X  = {2,         Name "Geometry/00X dimension"},
  Y  = {1,         Name "Geometry/01Y dimension"},
  N  = {4,         Name "Geometry/02Number of subdomain"},

  K  = {10,        Name "Geometry/02Wavenumber"},
  LB = {2* Pi / K, Name "Geometry/03Wavelength", ReadOnly 1},

  NL = {10,        Name "Geometry/04Points per wavelength"},
  NX = {X / LB * NL, Name "Geometry/05Mesh division", ReadOnly 1}
];

NY = NX;
CL = LB / NL;

// Geometry //
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

// Physicals
For n In {1:N}
  Physical Surface(n) = {n};
EndFor

For n In {1:N+1}
  Physical Line(n + N) = {n};
EndFor

Physical Line(2 * N + 2) = {1 + (N + 1):2 * (N + (N + 1))};

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
