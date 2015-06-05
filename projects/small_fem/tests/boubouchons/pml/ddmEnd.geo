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

// Pml Z //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = Boundary{ Volume{Pml~{0}[i]}; };
EndFor

For i In {0:EndN - 1}
  DdmEnd~{1}~{0}[i + 0 * EndN] = PmlZSurf~{i + 1 * RodN + 0 * EndN}[2]; // Right
  DdmEnd~{1}~{0}[i + 1 * EndN] = PmlZSurf~{i + 1 * RodN + 1 * EndN}[4]; // Right

  DdmEnd~{1}~{1}[i + 0 * EndN] = PmlZSurf~{i + 1 * RodN + 0 * EndN}[4]; // Left
  DdmEnd~{1}~{1}[i + 1 * EndN] = PmlZSurf~{i + 1 * RodN + 1 * EndN}[2]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{2}~{0}[i + 0 * EndN] = PmlZSurf~{i + 2 * RodN + 2 * EndN}[2]; // Right
  DdmEnd~{2}~{0}[i + 1 * EndN] = PmlZSurf~{i + 2 * RodN + 3 * EndN}[4]; // Right

  DdmEnd~{2}~{1}[i + 0 * EndN] = PmlZSurf~{i + 2 * RodN + 2 * EndN}[4]; // Left
  DdmEnd~{2}~{1}[i + 1 * EndN] = PmlZSurf~{i + 2 * RodN + 3 * EndN}[2]; // Left
EndFor

// Pml Y //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlYSurf~{i} = Boundary{ Volume{Pml~{1}[i]}; };
EndFor

For i In {0:EndN - 1}
  DdmEnd~{3}~{0}[i + 0 * EndN] = PmlYSurf~{i + 1 * RodN + 0 * EndN}[2]; // Right
  DdmEnd~{3}~{0}[i + 1 * EndN] = PmlYSurf~{i + 1 * RodN + 1 * EndN}[4]; // Right

  DdmEnd~{3}~{1}[i + 0 * EndN] = PmlYSurf~{i + 1 * RodN + 0 * EndN}[4]; // Left
  DdmEnd~{3}~{1}[i + 1 * EndN] = PmlYSurf~{i + 1 * RodN + 1 * EndN}[2]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{4}~{0}[i + 0 * EndN] = PmlYSurf~{i + 2 * RodN + 2 * EndN}[2]; // Right
  DdmEnd~{4}~{0}[i + 1 * EndN] = PmlYSurf~{i + 2 * RodN + 3 * EndN}[4]; // Right

  DdmEnd~{4}~{1}[i + 0 * EndN] = PmlYSurf~{i + 2 * RodN + 2 * EndN}[4]; // Left
  DdmEnd~{4}~{1}[i + 1 * EndN] = PmlYSurf~{i + 2 * RodN + 3 * EndN}[2]; // Left
EndFor

// Pml YZ //
For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = Boundary{ Volume{Pml~{2}[i]}; };
EndFor

For i In {0:EndN - 1}
  DdmEnd~{5}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 1 * RodN + 0 * EndN}[5]; // Right
  DdmEnd~{5}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 1 * RodN + 1 * EndN}[5]; // Right

  DdmEnd~{5}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 1 * RodN + 0 * EndN}[3]; // Left
  DdmEnd~{5}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 1 * RodN + 1 * EndN}[3]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{6}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 2 * RodN + 2 * EndN}[3]; // Right
  DdmEnd~{6}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 2 * RodN + 3 * EndN}[3]; // Right

  DdmEnd~{6}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 2 * RodN + 2 * EndN}[5]; // Left
  DdmEnd~{6}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 2 * RodN + 3 * EndN}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{7}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 3 * RodN + 4 * EndN}[3]; // Right
  DdmEnd~{7}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 3 * RodN + 5 * EndN}[3]; // Right

  DdmEnd~{7}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 3 * RodN + 4 * EndN}[5]; // Left
  DdmEnd~{7}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 3 * RodN + 5 * EndN}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{8}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 4 * RodN + 6 * EndN}[5]; // Right
  DdmEnd~{8}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 4 * RodN + 7 * EndN}[5]; // Right

  DdmEnd~{8}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 4 * RodN + 6 * EndN}[3]; // Left
  DdmEnd~{8}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 4 * RodN + 7 * EndN}[3]; // Left
EndFor

// Clear //
For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor

For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = {};
  PmlYSurf~{i} = {};
EndFor

For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = {};
EndFor
