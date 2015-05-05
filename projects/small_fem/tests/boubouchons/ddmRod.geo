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

// Pml Z //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = Boundary{ Volume{Pml~{0}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{1}~{0}[i] = PmlZSurf~{i}[3]; // Right
  DdmRod~{1}~{1}[i] = PmlZSurf~{i}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{2}~{0}[i] = PmlZSurf~{i + RodN + 2 * EndN}[3]; // Right
  DdmRod~{2}~{1}[i] = PmlZSurf~{i + RodN + 2 * EndN}[5]; // Left
EndFor

// Pml Y //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlYSurf~{i} = Boundary{ Volume{Pml~{1}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{3}~{0}[i] = PmlYSurf~{i}[3]; // Right
  DdmRod~{3}~{1}[i] = PmlYSurf~{i}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{4}~{0}[i] = PmlYSurf~{i + RodN + 2 * EndN}[5]; // Right
  DdmRod~{4}~{1}[i] = PmlYSurf~{i + RodN + 2 * EndN}[3]; // Left
EndFor

// Pml YZ //
For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = Boundary{ Volume{Pml~{2}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{5}~{0}[i] = PmlYZSurf~{i}[3]; // Right
  DdmRod~{5}~{1}[i] = PmlYZSurf~{i}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{6}~{0}[i] = PmlYZSurf~{i + 1 * RodN + 2 * EndN}[3]; // Right
  DdmRod~{6}~{1}[i] = PmlYZSurf~{i + 1 * RodN + 2 * EndN}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{7}~{0}[i] = PmlYZSurf~{i + 2 * RodN + 4 * EndN}[5]; // Right
  DdmRod~{7}~{1}[i] = PmlYZSurf~{i + 2 * RodN + 4 * EndN}[3]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{8}~{0}[i] = PmlYZSurf~{i + 3 * RodN + 6 * EndN}[5]; // Right
  DdmRod~{8}~{1}[i] = PmlYZSurf~{i + 3 * RodN + 6 * EndN}[3]; // Left
EndFor

// Clear //
For i In {0:RodN - 1}
  VolSurf~{i} = {};
EndFor

For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = {};
  PmlYSurf~{i} = {};
EndFor

For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = {};
EndFor
