// Reordering volumes //
////////////////////////////////////////////////////////////////////////
// Before : |    Rods    |Ends (both sides)|    Pmls (both sides)   | //
// After  : |Pmls (right)|Ends (right)|Rods|Ends (left)|Pmls (right)| //
////////////////////////////////////////////////////////////////////////

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

// Iterate on Cells //
Off = 0;
For i In {0:(PmlN - 1)}
  Cell~{i + Off} = {
    Pml~{3}[(PmlN - 1 - i)],         // Pml X
    Pml~{5}[(PmlN - 1 - i) * 2 + 0], // Pml XZ
    Pml~{5}[(PmlN - 1 - i) * 2 + 1], // ------
    Pml~{4}[(PmlN - 1 - i) * 2 + 0], // Pml XY
    Pml~{4}[(PmlN - 1 - i) * 2 + 1], // ------
    Pml~{6}[(PmlN - 1 - i) * 4 + 0], // Pml XYZ
    Pml~{6}[(PmlN - 1 - i) * 4 + 1], // -------
    Pml~{6}[(PmlN - 1 - i) * 4 + 2], // -------
    Pml~{6}[(PmlN - 1 - i) * 4 + 3]  // -------
  };
EndFor

Off = PmlN;
For i In {0:(EndN - 1)}
  Cell~{i + Off} = {
    eEnd[(EndN - 1 - i)],                          // End
    Pml~{0}[(EndN - 1 - i) + 1 * RodN + 0 * EndN], // Pml Z
    Pml~{0}[(EndN - 1 - i) + 2 * RodN + 2 * EndN], // -----
    Pml~{1}[(EndN - 1 - i) + 1 * RodN + 0 * EndN], // Pml Y
    Pml~{1}[(EndN - 1 - i) + 2 * RodN + 2 * EndN], // -----
    Pml~{2}[(EndN - 1 - i) + 1 * RodN + 0 * EndN], // Pml YZ
    Pml~{2}[(EndN - 1 - i) + 2 * RodN + 2 * EndN], // ------
    Pml~{2}[(EndN - 1 - i) + 3 * RodN + 4 * EndN], // ------
    Pml~{2}[(EndN - 1 - i) + 4 * RodN + 6 * EndN]  // ------
  };
EndFor

Off = PmlN + EndN;
For i In {0:(RodN - 1)}
  Cell~{i + Off} = {
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
EndFor

Off = PmlN + EndN + RodN;
For i In {0:(EndN - 1)}
  Cell~{i + Off} = {
    eEnd[(i + EndN)],                          // End
    Pml~{0}[(i + EndN) + 1 * RodN + 0 * EndN], // Pml Z
    Pml~{0}[(i + EndN) + 2 * RodN + 2 * EndN], // -----
    Pml~{1}[(i + EndN) + 1 * RodN + 0 * EndN], // Pml Y
    Pml~{1}[(i + EndN) + 2 * RodN + 2 * EndN], // -----
    Pml~{2}[(i + EndN) + 1 * RodN + 0 * EndN], // Pml YZ
    Pml~{2}[(i + EndN) + 2 * RodN + 2 * EndN], // ------
    Pml~{2}[(i + EndN) + 3 * RodN + 4 * EndN], // ------
    Pml~{2}[(i + EndN) + 4 * RodN + 6 * EndN]  // ------
  };
EndFor

Off = PmlN + EndN + RodN + EndN;
For i In {0:(PmlN - 1)}
  Cell~{i + Off} = {
    Pml~{3}[(i + PmlN)],         // Pml X
    Pml~{5}[(i + PmlN) * 2 + 0], // Pml XZ
    Pml~{5}[(i + PmlN) * 2 + 1], // ------
    Pml~{4}[(i + PmlN) * 2 + 0], // Pml XY
    Pml~{4}[(i + PmlN) * 2 + 1], // ------
    Pml~{6}[(i + PmlN) * 4 + 0], // Pml XYZ
    Pml~{6}[(i + PmlN) * 4 + 1], // -------
    Pml~{6}[(i + PmlN) * 4 + 2], // -------
    Pml~{6}[(i + PmlN) * 4 + 3]  // -------
  };
EndFor

// Clear //
Air  = {};
Rod  = {};
eEnd = {};
Off  = {};
