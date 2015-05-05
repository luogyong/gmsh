// Getting PML X boundaries for DDM //
//////////////////////////////////////

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
  // Down
  DdmPml~{1}~{0}[i + 0 * PmlN] = PmlXZSurf~{i * 2 + 0 + 0 * PmlN}[0]; // Right
  DdmPml~{1}~{0}[i + 1 * PmlN] = PmlXZSurf~{i * 2 + 0 + 2 * PmlN}[1]; // Right
  DdmPml~{1}~{1}[i + 0 * PmlN] = PmlXZSurf~{i * 2 + 0 + 0 * PmlN}[1]; // Left
  DdmPml~{1}~{1}[i + 1 * PmlN] = PmlXZSurf~{i * 2 + 0 + 2 * PmlN}[0]; // Left

  // Up
  DdmPml~{2}~{0}[i + 0 * PmlN] = PmlXZSurf~{i * 2 + 1 + 0 * PmlN}[0]; // Right
  DdmPml~{2}~{0}[i + 1 * PmlN] = PmlXZSurf~{i * 2 + 1 + 2 * PmlN}[1]; // Right
  DdmPml~{2}~{1}[i + 0 * PmlN] = PmlXZSurf~{i * 2 + 1 + 0 * PmlN}[1]; // Left
  DdmPml~{2}~{1}[i + 1 * PmlN] = PmlXZSurf~{i * 2 + 1 + 2 * PmlN}[0]; // Left
EndFor

// Pml XY //
For i In {0:(4 * PmlN) - 1}
  PmlXYSurf~{i} = Boundary{ Volume{Pml~{4}[i]}; };
EndFor

For i In {0:PmlN - 1}
  // Rear
  DdmPml~{3}~{0}[i + 0 * PmlN] = PmlXYSurf~{i * 2 + 0 + 0 * PmlN}[0]; // Right
  DdmPml~{3}~{0}[i + 1 * PmlN] = PmlXYSurf~{i * 2 + 0 + 2 * PmlN}[1]; // Right
  DdmPml~{3}~{1}[i + 0 * PmlN] = PmlXYSurf~{i * 2 + 0 + 0 * PmlN}[1]; // Left
  DdmPml~{3}~{1}[i + 1 * PmlN] = PmlXYSurf~{i * 2 + 0 + 2 * PmlN}[0]; // Left

  // Front
  DdmPml~{4}~{0}[i + 0 * PmlN] = PmlXYSurf~{i * 2 + 1 + 0 * PmlN}[0]; // Right
  DdmPml~{4}~{0}[i + 1 * PmlN] = PmlXYSurf~{i * 2 + 1 + 2 * PmlN}[1]; // Right
  DdmPml~{4}~{1}[i + 0 * PmlN] = PmlXYSurf~{i * 2 + 1 + 0 * PmlN}[1]; // Left
  DdmPml~{4}~{1}[i + 1 * PmlN] = PmlXYSurf~{i * 2 + 1 + 2 * PmlN}[0]; // Left
EndFor

// Pml XYZ //
For i In {0:(8 * PmlN) - 1}
  PmlXYZSurf~{i} = Boundary{ Volume{Pml~{6}[i]}; };
EndFor

For i In {0:PmlN - 1}
  // Rear Down
  DdmPml~{5}~{0}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 0 + 0 * PmlN}[0]; // Right
  DdmPml~{5}~{0}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 0 + 4 * PmlN}[1]; // Right
  DdmPml~{5}~{1}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 0 + 0 * PmlN}[1]; // Left
  DdmPml~{5}~{1}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 0 + 4 * PmlN}[0]; // Left

  // Rear Up
  DdmPml~{6}~{0}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 1 + 0 * PmlN}[0]; // Right
  DdmPml~{6}~{0}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 1 + 4 * PmlN}[1]; // Right
  DdmPml~{6}~{1}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 1 + 0 * PmlN}[1]; // Left
  DdmPml~{6}~{1}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 1 + 4 * PmlN}[0]; // Left

  // Front Down
  DdmPml~{7}~{0}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 2 + 0 * PmlN}[0]; // Right
  DdmPml~{7}~{0}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 2 + 4 * PmlN}[1]; // Right
  DdmPml~{7}~{1}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 2 + 0 * PmlN}[1]; // Left
  DdmPml~{7}~{1}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 2 + 4 * PmlN}[0]; // Left

  // Front Up
  DdmPml~{8}~{0}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 3 + 0 * PmlN}[0]; // Right
  DdmPml~{8}~{0}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 3 + 4 * PmlN}[1]; // Right
  DdmPml~{8}~{1}[i + 0 * PmlN] = PmlXYZSurf~{i * 4 + 3 + 0 * PmlN}[1]; // Left
  DdmPml~{8}~{1}[i + 1 * PmlN] = PmlXYZSurf~{i * 4 + 3 + 4 * PmlN}[0]; // Left
EndFor

// Clear //
For i In {0:(2 * PmlN) - 1}
  PmlXSurf~{i} = {};
EndFor

For i In {0:(4 * PmlN) - 1}
  PmlXZSurf~{i} = {};
  PmlXYSurf~{i} = {};
EndFor

For i In {0:(8 * PmlN) - 1}
  PmlXYZSurf~{i} = {};
EndFor
