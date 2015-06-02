// Getting ends boundaries for DDM //
/////////////////////////////////////

// Ends //
For i In {0:EndN - 1}
  EndSurf~{i}~{0} = Boundary{ Volume{End~{i}[0]}; };
  EndSurf~{i}~{1} = Boundary{ Volume{End~{i}[1]}; };
EndFor

For i In {0:EndN - 1}
  DdmEnd~{0}~{0}[i + 0 * EndN] = EndSurf~{i}~{0}[0]; // Right
  DdmEnd~{0}~{0}[i + 1 * EndN] = EndSurf~{i}~{1}[1]; // Right

  DdmEnd~{0}~{1}[i + 0 * EndN] = EndSurf~{i}~{0}[1]; // Left
  DdmEnd~{0}~{1}[i + 1 * EndN] = EndSurf~{i}~{1}[0]; // Left
EndFor

// Remove last boundaries //
DdmEnd~{0}~{1}[1 * EndN - 1] = 0; // Left
DdmEnd~{0}~{0}[2 * EndN - 1] = 0; // Right

// Clear //
For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor
