// PML Y //
///////////
// Pml~{1:2}[] = [... Pml(Y) Pml(YZ) ...]

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

// PML Z surfaces //
For i In {0:2 * (RodN + 2 * EndN) - 1}
  Tmp~{i} = Boundary{ Volume{Pml~{0}[i]}; };
EndFor

For i In {0:(RodN - 1)}
  Bnd~{0}[i] = Tmp~{i}[2];
  Bnd~{1}[i] = Tmp~{i}[4];
EndFor

For i In {RodN:(RodN + EndN - 1)}
  Bnd~{0}[i + 0 * EndN] = Tmp~{i + 0 * EndN}[3];
  Bnd~{0}[i + 1 * EndN] = Tmp~{i + 1 * EndN}[5];
  Bnd~{1}[i + 0 * EndN] = Tmp~{i + 0 * EndN}[5];
  Bnd~{1}[i + 1 * EndN] = Tmp~{i + 1 * EndN}[3];
EndFor

For i In {(RodN + 2 * EndN):(2 * RodN + 2 * EndN - 1)}
  Bnd~{0}[i] = Tmp~{i}[2];
  Bnd~{1}[i] = Tmp~{i}[4];
EndFor

For i In {(2 * RodN + 2 * EndN):(2 * RodN + 3 * EndN - 1)}
  Bnd~{0}[i + 0 * EndN] = Tmp~{i + 0 * EndN}[5];
  Bnd~{0}[i + 1 * EndN] = Tmp~{i + 1 * EndN}[3];
  Bnd~{1}[i + 0 * EndN] = Tmp~{i + 0 * EndN}[3];
  Bnd~{1}[i + 1 * EndN] = Tmp~{i + 1 * EndN}[5];
EndFor

// Y- //
Tmp~{0} = {};
Tmp~{1} = {};

For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[2];
EndFor

For i In {0:EndN - 1}
  Tmp~{0}[RodN + i + 0 * EndN] = EndSurf~{i}~{0}[3];
  Tmp~{0}[RodN + i + 1 * EndN] = EndSurf~{i}~{1}[5];
EndFor

Tmp~{1} = Extrude{0, -PmlSize, 0}{ Surface{Tmp~{0}[], Bnd~{0}[]}; };

For i In {0:RodN + 2 * EndN - 1}
  Pml~{1}[i + 0*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 0*2*(RodN + 2*EndN)];
  Pml~{2}[i + 0*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 1*2*(RodN + 2*EndN)];
  Pml~{2}[i + 1*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 2*2*(RodN + 2*EndN)];
EndFor

// Y+ //
Tmp~{0} = {};
Tmp~{1} = {};

For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[4];
EndFor

For i In {0:EndN - 1}
  Tmp~{0}[RodN + i + 0 * EndN] = EndSurf~{i}~{0}[5];
  Tmp~{0}[RodN + i + 1 * EndN] = EndSurf~{i}~{1}[3];
EndFor

Tmp~{1} = Extrude{0, +PmlSize, 0}{ Surface{Tmp~{0}[], Bnd~{1}[]}; };

For i In {0:RodN + 2 * EndN - 1}
  Pml~{1}[i + 1*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 0*2*(RodN + 2*EndN)];
  Pml~{2}[i + 2*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 1*2*(RodN + 2*EndN)];
  Pml~{2}[i + 3*(RodN + 2*EndN)] = Tmp~{1}[(2*i + 1) + 2*2*(RodN + 2*EndN)];
EndFor

// Clear useless variables //
For i In {1:RodN}
  VolSurf~{i} = {};
EndFor

For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor

For i In {0:2 * (RodN + 2) - 1}
  Tmp~{i} = {};
EndFor

Bnd~{0} = {};
Bnd~{1} = {};
