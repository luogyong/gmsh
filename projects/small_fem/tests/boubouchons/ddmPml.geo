// Getting PML X boundaries for DDM //
////////////////+/////////////////////

// Pml X //
For i In {0:(2 * PmlN) - 1}
  PmlXSurf~{i} = Boundary{ Volume{Pml~{3}[i]}; };
EndFor

For i In {0:PmlN - 1}
  DdmPml~{0}~{0}[i + 0 * PmlN] = PmlXSurf~{i + 0 * PmlN}[0]; // Right
  DdmPml~{0}~{0}[i + 1 * PmlN] = PmlXSurf~{i + 1 * PmlN}[1]; // Right

  DdmPml~{0}~{1}[i + 0 * PmlN] = PmlXSurf~{i + 0 * PmlN}[1]; // Left
  DdmPml~{0}~{1}[i + 1 * PmlN] = PmlXSurf~{i + 1 * PmlN}[0]; // Left
EndFor

// Pml XZ //
For i In {0:(4 * PmlN) - 1}
  PmlXZSurf~{i} = Boundary{ Volume{Pml~{5}[i]}; };
EndFor

For i In {0:PmlN - 1}
  DdmPml~{1}~{0}[i + 0 * PmlN] = PmlXZSurf~{i + 0 * PmlN}[0]; // Right
  DdmPml~{1}~{0}[i + 1 * PmlN] = PmlXZSurf~{i + 1 * PmlN}[1]; // Right

  DdmPml~{1}~{1}[i + 0 * PmlN] = PmlXZSurf~{i + 0 * PmlN}[1]; // Left
  DdmPml~{1}~{1}[i + 1 * PmlN] = PmlXZSurf~{i + 1 * PmlN}[0]; // Left
EndFor

For i In {0:5}
  Printf("%f: %f", i, PmlXZSurf~{0 + 1 * PmlN}[i]);
EndFor

For i In {0:(2 * PmlN) - 1}
  //Printf("%f: %f", i, DdmPml~{0}~{0}[i]);
EndFor
/*
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
*/
// Clear //
For i In {0:(2 * PmlN) - 1}
  PmlXSurf~{i} = {};
EndFor
For i In {0:(4 * PmlN) - 1}
  PmlXZSurf~{i} = {};
EndFor
