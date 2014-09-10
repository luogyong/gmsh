DefineConstant[ N_DOM = {2,   Name "Geometry/Number of domain"} ];
DefineConstant[ r     = {1/3, Name "Geometry/Small radius"}     ];
DefineConstant[ R     = {3/3, Name "Geometry/Big radius"}       ];
DefineConstant[ lc    = {0.1, Name "Geometry/Mesh density"}     ];

For i In {1:N_DOM}
  R~{i} = r + i * (R - r) / (N_DOM);
EndFor

// Points //
Point(0) = {0, 0, 0, lc};

Point(1001) = {+r,  0, 0, lc};
Point(1002) = { 0, +r, 0, lc};
Point(1003) = {-r,  0, 0, lc};
Point(1004) = { 0, -r, 0, lc};

For i In {1:N_DOM}
  Point((i + 1) * 1000 + 1) = {+R~{i},  0, 0, lc};
  Point((i + 1) * 1000 + 2) = { 0, +R~{i}, 0, lc};
  Point((i + 1) * 1000 + 3) = {-R~{i},  0, 0, lc};
  Point((i + 1) * 1000 + 4) = { 0, -R~{i}, 0, lc};
EndFor

// Circles //
For i In {0:N_DOM}
  Circle((i + 1) * 1000 + 1) = {(i + 1) * 1000 + 1,  0, (i + 1) * 1000 + 2};
  Circle((i + 1) * 1000 + 2) = {(i + 1) * 1000 + 2,  0, (i + 1) * 1000 + 3};
  Circle((i + 1) * 1000 + 3) = {(i + 1) * 1000 + 3,  0, (i + 1) * 1000 + 4};
  Circle((i + 1) * 1000 + 4) = {(i + 1) * 1000 + 4,  0, (i + 1) * 1000 + 1};
EndFor

// Line Loops & Surfaces //
For i In {0:N_DOM}
  Line Loop((i + 1) * 1000) = {(i + 1) * 1000 + 1, (i + 1) * 1000 + 2,
                               (i + 1) * 1000 + 3, (i + 1) * 1000 + 4};
EndFor

For i In {0:N_DOM - 1}
  Plane Surface((i + 1) * 1000) = {(i + 2) * 1000, (i + 1) * 1000};
EndFor

// Physicals //
Physical Line(1001) = {1001, 1002, 1003, 1004};             // Src

For i In {1:N_DOM - 1}
  Physical Line(1000 + i + 1) = {(i + 1) * 1000 + 1,
                                 (i + 1) * 1000 + 2,
                                 (i + 1) * 1000 + 3,
                                 (i + 1) * 1000 + 4};       // DDM
EndFor

Physical Line(1000 + N_DOM + 1) = {(N_DOM + 1) * 1000 + 1,
                                   (N_DOM + 1) * 1000 + 2,
                                   (N_DOM + 1) * 1000 + 3,
                                   (N_DOM + 1) * 1000 + 4}; // Inf

For i In {1:N_DOM}
  Physical Surface(3000 + i) = {i * 1000};                  // Dom
EndFor
