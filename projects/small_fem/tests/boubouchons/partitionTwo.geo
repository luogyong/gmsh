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

// Air & Rods //
For i In {0:RodN - 1}
  If(i < (RodN / 2))
    part = 1;
  EndIf
  If(i >= (RodN / 2))
    part = 2;
  EndIf

  SetPartition (part) {
    Volume{
      Air[i],                           // Air
      Rod[i],                           // Rod
      Pml~{0}[i + 0 * RodN + 0 * EndN], // Pml Z
      Pml~{0}[i + 1 * RodN + 2 * EndN], // -----
      Pml~{1}[i + 0 * RodN + 0 * EndN], // Pml Y
      Pml~{1}[i + 1 * RodN + 2 * EndN], // -----
      Pml~{2}[i + 0 * RodN + 0 * EndN], // Pml YZ
      Pml~{2}[i + 1 * RodN + 2 * EndN], // ------
      Pml~{2}[i + 2 * RodN + 4 * EndN], // ------
      Pml~{2}[i + 3 * RodN + 6 * EndN]  // ------
    };
  }
EndFor

// Source //
SetPartition (1) { Surface{ Src~{1}[] }; } // Source

// Ends //
For i In {0:(2 * EndN) - 1}
  If(i < (EndN))
    part = 1;
  EndIf
  If(i >= (EndN))
    part = 2;
  EndIf

  SetPartition (part) {
    Volume{
      eEnd[i],                          // End
      Pml~{0}[i + 1 * RodN + 0 * EndN], // Pml Z
      Pml~{0}[i + 2 * RodN + 2 * EndN], // -----
      Pml~{1}[i + 1 * RodN + 0 * EndN], // Pml Y
      Pml~{1}[i + 2 * RodN + 2 * EndN], // -----
      Pml~{2}[i + 1 * RodN + 0 * EndN], // Pml YZ
      Pml~{2}[i + 2 * RodN + 2 * EndN], // ------
      Pml~{2}[i + 3 * RodN + 4 * EndN], // ------
      Pml~{2}[i + 4 * RodN + 6 * EndN]  // ------
    };
  }
EndFor

// Pml X //
For i In {0:(2 * PmlN) - 1}
  If(i < (PmlN))
    part = 1;
  EndIf
  If(i >= (PmlN))
    part = 2;
  EndIf

  SetPartition (part) {
    Volume{
      Pml~{3}[i],         // Pml X
      Pml~{5}[i * 2 + 0], // Pml XZ
      Pml~{5}[i * 2 + 1], // ------
      Pml~{4}[i * 2 + 0], // Pml XY
      Pml~{4}[i * 2 + 1], // ------
      Pml~{6}[i * 4 + 0], // Pml XYZ
      Pml~{6}[i * 4 + 1], // -------
      Pml~{6}[i * 4 + 2], // -------
      Pml~{6}[i * 4 + 3]  // -------
    };
  }
EndFor

// Clear
Air  = {};
Rod  = {};
eEnd = {};
