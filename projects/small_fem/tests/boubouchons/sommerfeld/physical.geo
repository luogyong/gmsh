// Physical when no DDM is used //
//////////////////////////////////

// Grab Air & Rods //
For i In {0:RodN - 1}
  Air[i] = Vol~{i + 1}[1];
  Rod[i] = Vol~{i + 1}[0];
EndFor
For i In {0:EndN - 1}
  Air[RodN + i + 0 * EndN] = End~{i}[0];
  Air[RodN + i + 1 * EndN] = End~{i}[1];
EndFor

// Grab Infinity //
Inf = {};
For i In {0:CellN - 1}
  Inf += Infinity~{i}[];
EndFor

// Physicals //
Physical   Volume(1007) = Air[];     // Air
Physical   Volume(1008) = Rod[];     // Rods
Physical  Surface(1009) = Src~{1}[]; // Source
Physical  Surface(1010) = Inf[];     // Infinity

// Boundaries for DDM
For i In {0:(DomN - 1)}
  Physical Surface(40000 + i + 1) = { DdmBoundary~{Range~{i}[0]}~{1}[],
                                      DdmBoundary~{Range~{i}[1]}~{0}[] };
EndFor

// Clear //
Air = {};
Rod = {};
Inf = {};