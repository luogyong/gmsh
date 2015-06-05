// PML Z //
///////////
// Pml~{0}[] = [Pml(Z) ...]

// Grab Surfaces //
// VolSurf~{i}[0:5] = [Air(Z-), Air(Z+), Air(Y-), Air(X+), Air(Y+), Air(X-)]
For i In {1:RodN}
  VolSurf~{i} = Boundary{ Volume{Vol~{i}[1]}; };
EndFor

// EndSurf~{i}~{0}[0:5] = [End(X+), End(X-), End(Z-), End(Y-), End(Z+), End(Y+)]
// EndSurf~{i}~{1}[0:5] = [End(X-), End(X+), End(Z-), End(Y+), End(Z+), End(Y-)]
For i In {0:EndN - 1}
  EndSurf~{i}~{0} = Boundary{ Volume{End~{i}[0]}; };
  EndSurf~{i}~{1} = Boundary{ Volume{End~{i}[1]}; };
EndFor


// Z- //
For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[0];
EndFor

For i In {0:EndN - 1}
  Tmp~{0}[RodN + i + 0 * EndN] = EndSurf~{i}~{0}[2];
  Tmp~{0}[RodN + i + 1 * EndN] = EndSurf~{i}~{1}[2];
EndFor

Tmp~{1} = Extrude{0, 0, -PmlSize}{ Surface{Tmp~{0}[]}; };
For i In {0:RodN + 2 * EndN - 1}
  Pml~{0}[i + 0 * (RodN + 2 * EndN)] = Tmp~{1}[i * 2 + 1];
EndFor

// Z+ //
For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[1];
EndFor

For i In {0:EndN - 1}
  Tmp~{0}[RodN + i + 0 * EndN] = EndSurf~{i}~{0}[4];
  Tmp~{0}[RodN + i + 1 * EndN] = EndSurf~{i}~{1}[4];
EndFor

Tmp~{1} = Extrude{0, 0, +PmlSize}{ Surface{Tmp~{0}[]}; };
For i In {0:RodN + 2 * EndN - 1}
  Pml~{0}[i + 1 * (RodN + 2 * EndN)] = Tmp~{1}[i * 2 + 1];
EndFor

// Clear useless variables //
For i In {1:RodN}
  VolSurf~{i} = {};
EndFor

For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor

Tmp~{0} = {};
Tmp~{1} = {};
