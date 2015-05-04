// Physical when DDM is used //
///////////////////////////////

// Grab Air & Rods
For i In {0:RodN - 1}
  Air[i] = Vol~{i + 1}[1];
  Rod[i] = Vol~{i + 1}[0];
EndFor

For i In {0:EndN - 1}
  eEnd[i + 0 * EndN] = End~{i}[0];
  eEnd[i + 1 * EndN] = End~{i}[1];
EndFor

// Physicals
/*
Physical Volume(1000) = Pml~{6}[]; // Pml XYZ
Physical Volume(1001) = Pml~{5}[]; // Pml XZ
Physical Volume(1002) = Pml~{2}[]; // Pml YZ
Physical Volume(1003) = Pml~{4}[]; // Pml XY
Physical Volume(1004) = Pml~{0}[]; // Pml Z
Physical Volume(1005) = Pml~{1}[]; // Pml Y
Physical Volume(1006) = Pml~{3}[]; // Pml X
*/

Physical Volume(10) = Src~{0};   // Source
Physical   Line(11) = Src~{1};   // Source line

For i In {0:RodN - 1}
  Physical  Volume(10000 + i) = Air[i];                // Air
  Physical  Volume(20000 + i) = Rod[i];                // Rod

  Physical Surface(30000 + i) =   DdmRod~{0}~{0}[i];   // Right Air
  Physical Surface(40000 + i) = { DdmRod~{1}~{0}[i],
                                  DdmRod~{2}~{0}[i] }; // Right Pml Z
  Physical Surface(50000 + i) = { DdmRod~{3}~{0}[i],
                                  DdmRod~{4}~{0}[i] }; // Right Pml Y
  Physical Surface(60000 + i) = { DdmRod~{5}~{0}[i],
                                  DdmRod~{6}~{0}[i],
                                  DdmRod~{7}~{0}[i],
                                  DdmRod~{8}~{0}[i] }; // Right Pml YZ
EndFor

For i In {0:(2 * EndN) - 1}
  Physical  Volume( 70000 + i) = eEnd[i];               // End

  Physical Surface( 80000 + i) =   DdmEnd~{0}~{0}[i];   // Right End
  Physical Surface( 90000 + i) = { DdmEnd~{1}~{0}[i],
                                   DdmEnd~{2}~{0}[i] }; // Right End Pml Z
  Physical Surface(100000 + i) = { DdmEnd~{3}~{0}[i],
                                   DdmEnd~{4}~{0}[i] }; // Right End Pml Y
  Physical Surface(110000 + i) = { DdmEnd~{5}~{0}[i],
                                   DdmEnd~{6}~{0}[i],
                                   DdmEnd~{7}~{0}[i],
                                   DdmEnd~{8}~{0}[i] }; // Right End Pml YZ
EndFor

For i In {0:(2 * PmlN) - 1}
  Physical  Volume(120000 + i) = Pml~{3}[i]; // Pml X
EndFor

// Clear
Air  = {};
Rod  = {};
eEnd = {};