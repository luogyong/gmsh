// Grab Surfaces //
///////////////////
// VolSurf~{i}[0:5] = [Air(Z-), Air(Z+), Air(Y-), Air(X+), Air(Y+), Air(X-)]
For i In {1:RodN}
  VolSurf~{i} = Boundary{ Volume{Vol~{i}[1]}; };
EndFor

// EndSurf~{0}[0:5] = [End(X+), End(X-), End(Z-), End(Y-), End(Z+), End(Y+)]
// EndSurf~{1}[0:5] = [End(X-), End(X+), End(Z-), End(Y+), End(Z+), End(Y-)]
EndSurf~{0} = Boundary{ Volume{End~{0}}; };
EndSurf~{1} = Boundary{ Volume{End~{1}}; };

// PML Z surfaces
For i In {0:2 * (RodN + 2) - 1}
  Tmp~{i} = Boundary{ Volume{Pml~{0}[i]}; };
EndFor

For i In {0:RodN - 1}
  Bnd~{0}[i + 0 * RodN] = Tmp~{i + 0 * (RodN + 2)}[2];
  Bnd~{0}[i + 1 * RodN] = Tmp~{i + 1 * (RodN + 2)}[2];
  Bnd~{1}[i + 0 * RodN] = Tmp~{i + 0 * (RodN + 2)}[4];
  Bnd~{1}[i + 1 * RodN] = Tmp~{i + 1 * (RodN + 2)}[4];
EndFor

Bnd~{0}[2 * RodN + 0] = Tmp~{RodN + 0 + 0 * (RodN + 2)}[3];
Bnd~{0}[2 * RodN + 1] = Tmp~{RodN + 1 + 0 * (RodN + 2)}[5];
Bnd~{0}[2 * RodN + 2] = Tmp~{RodN + 0 + 1 * (RodN + 2)}[5];
Bnd~{0}[2 * RodN + 3] = Tmp~{RodN + 1 + 1 * (RodN + 2)}[3];

Bnd~{1}[2 * RodN + 0] = Tmp~{RodN + 0 + 0 * (RodN + 2)}[5];
Bnd~{1}[2 * RodN + 1] = Tmp~{RodN + 1 + 0 * (RodN + 2)}[3];
Bnd~{1}[2 * RodN + 2] = Tmp~{RodN + 0 + 1 * (RodN + 2)}[3];
Bnd~{1}[2 * RodN + 3] = Tmp~{RodN + 1 + 1 * (RodN + 2)}[5];

// PML Y //
///////////
// Pml~{1:2}[] = [... Pml(Y) Pml(YZ) ...]

// Y- //
Tmp~{0} = {};
Tmp~{1} = {};

For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[2];
EndFor

Tmp~{0}[RodN + 0] = EndSurf~{0}[3];
Tmp~{0}[RodN + 1] = EndSurf~{1}[5];

Tmp~{1} = Extrude{0, -PmlSize, 0}{ Surface{Tmp~{0}[], Bnd~{0}[]}; };

For i In {0:RodN + 1}
  Pml~{1}[i + 0 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 0 * 2 * (RodN + 2)];
  Pml~{2}[i + 0 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 1 * 2 * (RodN + 2)];
  Pml~{2}[i + 1 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 2 * 2 * (RodN + 2)];
EndFor

// Y+ //
Tmp~{0} = {};
Tmp~{1} = {};

For i In {0:RodN - 1}
  Tmp~{0}[i] = VolSurf~{i + 1}[4];
EndFor

Tmp~{0}[RodN + 0] = EndSurf~{0}[5];
Tmp~{0}[RodN + 1] = EndSurf~{1}[3];

Tmp~{1} = Extrude{0, +PmlSize, 0}{ Surface{Tmp~{0}[], Bnd~{1}[]}; };

For i In {0:RodN + 1}
  Pml~{1}[i + 1 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 0 * 2 * (RodN + 2)];
  Pml~{2}[i + 2 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 1 * 2 * (RodN + 2)];
  Pml~{2}[i + 3 * (RodN + 2)] = Tmp~{1}[(i * 2 + 1) + 2 * 2 * (RodN + 2)];
EndFor

// Clear useless variables //
/////////////////////////////
For i In {1:RodN}
  VolSurf~{i} = {};
EndFor

EndSurf~{0} = {};
EndSurf~{1} = {};

For i In {0:2 * (RodN + 2) - 1}
  Tmp~{i} = {};
EndFor

Bnd~{0} = {};
Bnd~{1} = {};
