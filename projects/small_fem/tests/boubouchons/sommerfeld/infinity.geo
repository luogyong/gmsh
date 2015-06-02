// Infinity Boundaries //
/////////////////////////

// Iterate on Rods //
For i In {1:RodN}
  Bnd = {};
  Bnd = Boundary{ Volume{Vol~{i}[1]}; };

  InfinityRod~{i - 1} = {Bnd[0], Bnd[1], Bnd[2], Bnd[4]};
EndFor

// Iterate on Ends //
For i In {0:EndN - 1}
  BndR = {};
  BndL = {};
  BndR = Boundary{ Volume{End~{i}[0]}; };
  BndL = Boundary{ Volume{End~{i}[1]}; };

  InfinityEnd~{i * 2 + 0} = {BndR[2], BndR[3], BndR[4], BndR[5]};
  InfinityEnd~{i * 2 + 1} = {BndL[2], BndL[3], BndL[4], BndL[5]};

  If(i == (EndN - 1)) // Last Ends
    InfinityEnd~{i * 2 + 0} += BndR[1];
    InfinityEnd~{i * 2 + 1} += BndL[1];
  EndIf
EndFor

// Clear //
Bnd  = {};
BndR = {};
BndL = {};
