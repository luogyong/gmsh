// Reordering volumes //
/////////////////////////////////////////////////////////////////
// Before : |     Rods     | Ends (both sides) |               //
// After  : | Ends (right) |        Rods       | Ends (left) | //
/////////////////////////////////////////////////////////////////

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
For i In {0:(EndN - 1)}
  Cell~{i + Off} = eEnd[(EndN - 1 - i)];
EndFor

Off = EndN;
For i In {0:(RodN - 1)}
  Cell~{i + Off} = { Air[i], Rod[i] };
EndFor

Off = EndN + RodN;
For i In {0:(EndN - 1)}
  Cell~{i + Off} = eEnd[(i + EndN)];
EndFor

// Infinity //
Off = 0;
For i In {0:(EndN - 1)}
  Infinity~{i + Off} = InfinityEnd~{(EndN - i - 1) * 2}[];
EndFor

Off = EndN;
For i In {0:(RodN - 1)}
  Infinity~{i + Off} = InfinityRod~{i}[];
EndFor

Off = EndN + RodN;
For i In {0:(EndN - 1)}
  Infinity~{i + Off} = InfinityEnd~{i * 2 + 1}[];
EndFor

// Clear //
Air  = {};
Rod  = {};
eEnd = {};
Off  = {};
