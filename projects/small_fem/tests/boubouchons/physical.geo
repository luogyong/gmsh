// Physical when no DDM is used //
//////////////////////////////////

// Grab Air & Rods
For i In {0:RodN - 1}
  Air[i] = Vol~{i + 1}[1];
  Rod[i] = Vol~{i + 1}[0];
EndFor

For i In {0:EndN - 1}
  Air[RodN + i + 0 * EndN] = End~{i}[0];
  Air[RodN + i + 1 * EndN] = End~{i}[1];
EndFor

// Physicals
Physical  Volume(1000) = Pml~{6}[]; // Pml XYZ
Physical  Volume(1001) = Pml~{5}[]; // Pml XZ
Physical  Volume(1002) = Pml~{2}[]; // Pml YZ
Physical  Volume(1003) = Pml~{4}[]; // Pml XY
Physical  Volume(1004) = Pml~{0}[]; // Pml Z
Physical  Volume(1005) = Pml~{1}[]; // Pml Y
Physical  Volume(1006) = Pml~{3}[]; // Pml X

Physical  Volume(1007) = Air[];     // Air
Physical  Volume(1008) = Rod[];     // Rods

Physical  Volume(1009) = Src~{0};   // Source
Physical    Line(1011) = Src~{1};   // Source line

// Boundaries for 'Rod' (Right / Left)
For i In {0:RodN - 1}
  Physical Surface(20000 + i + 1) = DdmSortRod~{i}~{0}[];
  Physical Surface(30000 + i + 1) = DdmSortRod~{i}~{1}[];
EndFor

// Boundaries for 'End' (Right / Left)
For i In {0:(2 * EndN - 1)}
  Physical Surface(20000 + i + 1 + RodN) = DdmSortEnd~{i}~{0}[];
  Physical Surface(30000 + i + 1 + RodN) = DdmSortEnd~{i}~{1}[];
EndFor

// Boundaries for 'PmlX' (Right / Left)
For i In {0:(2 * PmlN - 1)}
  Physical Surface(20000 + i + 1 + RodN + 2 * EndN) = DdmSortPml~{i}~{0}[];
  Physical Surface(30000 + i + 1 + RodN + 2 * EndN) = DdmSortPml~{i}~{1}[];
EndFor

// Clear
Air = {};
Rod = {};
