// Physical when DDM is used //
///////////////////////////////

// Grab Air & Rods //
For i In {0:RodN - 1}
  Air[i] = Vol~{i + 1}[1];
  Rod[i] = Vol~{i + 1}[0];
EndFor

// Grab Ends //
For i In {0:EndN - 1}
  eEnd[i + 0 * EndN] = End~{i}[0];
  eEnd[i + 1 * EndN] = End~{i}[1];
EndFor

// Source //
Physical Volume(10) = Src~{0};   // Source volume
Physical   Line(11) = Src~{1};   // Source line

// Air & Rods //
Off = 0;
For i In {0:RodN - 1}
  // Volumes
  Physical  Volume(110000 + Off + i) =   Air[i];            // Air
  Physical  Volume(120000 + Off + i) =   Rod[i];            // Rod
  Physical  Volume(130000 + Off + i) = {
    Pml~{0}[i + 0 * RodN + 0 * EndN],
    Pml~{0}[i + 1 * RodN + 2 * EndN]
  };                                                        // Pml Z
  Physical  Volume(140000 + Off + i) = {
    Pml~{1}[i + 0 * RodN + 0 * EndN],
    Pml~{1}[i + 1 * RodN + 2 * EndN]
  };                                                        // Pml Y
  Physical  Volume(150000 + Off + i) = {
    Pml~{2}[i + 0 * RodN + 0 * EndN],
    Pml~{2}[i + 1 * RodN + 2 * EndN],
    Pml~{2}[i + 2 * RodN + 4 * EndN],
    Pml~{2}[i + 3 * RodN + 6 * EndN]
  };                                                        // Pml YZ

  // Right boundaries
  Physical Surface(210000 + Off + i) =   DdmRod~{0}~{0}[i]; // Air
  Physical Surface(230000 + Off + i) = {
    DdmRod~{1}~{0}[i],
    DdmRod~{2}~{0}[i]
  };                                                        // Pml Z
  Physical Surface(240000 + Off + i) = {
    DdmRod~{3}~{0}[i],
    DdmRod~{4}~{0}[i]
  };                                                        // Pml Y
  Physical Surface(250000 + Off + i) = {
    DdmRod~{5}~{0}[i],
    DdmRod~{6}~{0}[i],
    DdmRod~{7}~{0}[i],
    DdmRod~{8}~{0}[i]
  };                                                        // Pml YZ

  // Left boundaries
  Physical Surface(310000 + Off + i) =   DdmRod~{0}~{1}[i]; // Air
  Physical Surface(330000 + Off + i) = {
    DdmRod~{1}~{1}[i],
    DdmRod~{2}~{1}[i]
  };                                                        // Pml Z
  Physical Surface(340000 + Off + i) = {
    DdmRod~{3}~{1}[i],
    DdmRod~{4}~{1}[i]
  };                                                        // Pml Y
  Physical Surface(350000 + Off + i) = {
    DdmRod~{5}~{1}[i],
    DdmRod~{6}~{1}[i],
    DdmRod~{7}~{1}[i],
    DdmRod~{8}~{1}[i]
  };                                                        // Pml YZ
EndFor

// Ends //
Off = RodN;
For i In {0:(2 * EndN) - 1}
  // Volumes
  Physical  Volume(110000 + Off + i) =   eEnd[i];           // End
  Physical  Volume(130000 + Off + i) = {
    Pml~{0}[i + 1 * RodN + 0 * EndN],
    Pml~{0}[i + 2 * RodN + 2 * EndN]
  };                                                        // Pml Z
  Physical  Volume(140000 + Off + i) = {
    Pml~{1}[i + 1 * RodN + 0 * EndN],
    Pml~{1}[i + 2 * RodN + 2 * EndN]
  };                                                        // Pml Y
  Physical  Volume(150000 + Off + i) = {
    Pml~{2}[i + 1 * RodN + 0 * EndN],
    Pml~{2}[i + 2 * RodN + 2 * EndN],
    Pml~{2}[i + 3 * RodN + 4 * EndN],
    Pml~{2}[i + 4 * RodN + 6 * EndN]
  };                                                        // Pml YZ

  // Right boundaries
  Physical Surface(210000 + Off + i) =   DdmEnd~{0}~{0}[i]; // End
  Physical Surface(230000 + Off + i) = {
    DdmEnd~{1}~{0}[i],
    DdmEnd~{2}~{0}[i] };                                    // Pml Z
  Physical Surface(240000 + Off + i) = {
    DdmEnd~{3}~{0}[i],
    DdmEnd~{4}~{0}[i] };                                    // Pml Y
  Physical Surface(250000 + Off + i) = {
    DdmEnd~{5}~{0}[i],
    DdmEnd~{6}~{0}[i],
    DdmEnd~{7}~{0}[i],
    DdmEnd~{8}~{0}[i]
  };                                                        // Pml YZ

  // Left boundaries
  Physical Surface(310000 + Off + i) =   DdmEnd~{0}~{1}[i]; // End
  Physical Surface(330000 + Off + i) = {
    DdmEnd~{1}~{1}[i],
    DdmEnd~{2}~{1}[i] };                                    // Pml Z
  Physical Surface(340000 + Off + i) = {
    DdmEnd~{3}~{1}[i],
    DdmEnd~{4}~{1}[i] };                                    // Pml Y
  Physical Surface(350000 + Off + i) = {
    DdmEnd~{5}~{1}[i],
    DdmEnd~{6}~{1}[i],
    DdmEnd~{7}~{1}[i],
    DdmEnd~{8}~{1}[i]
  };                                                        // Pml YZ
EndFor

// Pml X //
Off = RodN + 2 * EndN;
For i In {0:(2 * PmlN) - 1}
  // Volumes
  Physical  Volume(110000 + Off + i) =   Pml~{3}[i];        // Pml X
  Physical  Volume(130000 + Off + i) = {
    Pml~{5}[i * 2 + 0],
    Pml~{5}[i * 2 + 1]
  };                                                        // Pml XZ
  Physical  Volume(140000 + Off + i) = {
    Pml~{4}[i * 2 + 0],
    Pml~{4}[i * 2 + 1]
  };                                                        // Pml XY
  Physical  Volume(150000 + Off + i) = {
    Pml~{6}[i * 4 + 0],
    Pml~{6}[i * 4 + 1],
    Pml~{6}[i * 4 + 2],
    Pml~{6}[i * 4 + 3]
  };                                                        // Pml XYZ

  // Right boundaries
  Physical Surface(210000 + Off + i) =   DdmPml~{0}~{0}[i]; // Pml X
  Physical Surface(230000 + Off + i) = {
    DdmPml~{1}~{0}[i],
    DdmPml~{2}~{0}[i]
  };                                                        // Pml XZ
  Physical Surface(240000 + Off + i) = {
    DdmPml~{3}~{0}[i],
    DdmPml~{4}~{0}[i]
  };                                                        // Pml XY
  Physical Surface(250000 + Off + i) = {
    DdmPml~{5}~{0}[i],
    DdmPml~{6}~{0}[i],
    DdmPml~{7}~{0}[i],
    DdmPml~{8}~{0}[i]
  };                                                        // Pml XYZ

  // Left boundaries
  Physical Surface(310000 + Off + i) =   DdmPml~{0}~{1}[i]; // Pml X
  Physical Surface(330000 + Off + i) = {
    DdmPml~{1}~{1}[i],
    DdmPml~{2}~{1}[i]
  };                                                        // Pml XZ
  Physical Surface(340000 + Off + i) = {
    DdmPml~{3}~{1}[i],
    DdmPml~{4}~{1}[i]
  };                                                        // Pml XY
  Physical Surface(350000 + Off + i) = {
    DdmPml~{5}~{1}[i],
    DdmPml~{6}~{1}[i],
    DdmPml~{7}~{1}[i],
    DdmPml~{8}~{1}[i]
  };                                                        // Pml XYZ
EndFor

// Clear
Air  = {};
Rod  = {};
eEnd = {};
Off  = {};