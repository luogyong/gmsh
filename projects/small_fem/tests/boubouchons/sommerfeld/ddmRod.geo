// Getting cells boundaries for DDM //
//////////////////////////////////////

// Rod //
For i In {0:RodN - 1}
  VolSurf~{i} = Boundary{ Volume{Vol~{i + 1}[1]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{0}~{0}[i] = VolSurf~{i}[3]; // Right
  DdmRod~{0}~{1}[i] = VolSurf~{i}[5]; // Left
EndFor

// Clear //
For i In {0:RodN - 1}
  VolSurf~{i} = {};
EndFor
