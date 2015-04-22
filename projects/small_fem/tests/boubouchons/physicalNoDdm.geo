// Physical when no DDM is used //
//////////////////////////////////

// Grab Air & Rods
For i In {0:RodN - 1}
  Air[i] = Vol~{i + 1}[1];
  Rod[i] = Vol~{i + 1}[0];
EndFor

Air[RodN + 0] = End~{0}[0];
Air[RodN + 1] = End~{0}[1];

// Physicals
Physical Volume(1000) = Pml~{6}[]; // Pml XYZ
Physical Volume(1001) = Pml~{5}[]; // Pml XZ
Physical Volume(1002) = Pml~{2}[]; // Pml YZ
Physical Volume(1003) = Pml~{4}[]; // Pml XY
Physical Volume(1004) = Pml~{0}[]; // Pml Z
Physical Volume(1005) = Pml~{1}[]; // Pml Y
Physical Volume(1006) = Pml~{3}[]; // Pml X

Physical Volume(1007) = Air[];     // Air
Physical Volume(1008) = Rod[];     // Rods

Physical Volume(1009) = Src~{0};   // Source
Physical   Line(1011) = Src~{1};   // Source line

// Clear
Air = {};
Rod = {};
