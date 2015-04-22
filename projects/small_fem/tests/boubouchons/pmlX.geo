// Grab Surfaces //
///////////////////
// VolSurf~{i}[0:5] = [Air(Z-), Air(Z+), Air(Y-), Air(X+), Air(Y+), Air(X-)]
For i In {1:RodN}
  VolSurf~{i} = Boundary{ Volume{Vol~{i}[1]}; };
EndFor

// EndSurf~{0}[0:5] = [End(X+), End(X-), End(Z-), End(Y-), End(Z+), End(Y+)]
// EndSurf~{1}[0:5] = [End(X-), End(X+), End(Z-), End(Y+), End(Z+), End(Y-)]
EndSurf~{0} = Boundary{ Volume{End~{0}[0]}; };
EndSurf~{1} = Boundary{ Volume{End~{0}[1]}; };

// PML ZY surfaces
For i In {0:4 * (RodN + 2) - 1}
  Tmp~{0}~{i} = Boundary{ Volume{Pml~{2}[i]}; };
EndFor

Bnd~{0} = {};
Bnd~{1} = {};

Bnd~{0}[0] = Tmp~{0}~{2 * RodN + 0 + 0}[3];
Bnd~{0}[1] = Tmp~{0}~{2 * RodN + 0 + 2}[5];
Bnd~{0}[2] = Tmp~{0}~{4 * RodN + 4 + 0}[5];
Bnd~{0}[3] = Tmp~{0}~{4 * RodN + 4 + 2}[3];

Bnd~{1}[0] = Tmp~{0}~{2 * RodN + 0 + 1}[5];
Bnd~{1}[1] = Tmp~{0}~{2 * RodN + 0 + 3}[3];
Bnd~{1}[2] = Tmp~{0}~{4 * RodN + 4 + 1}[3];
Bnd~{1}[3] = Tmp~{0}~{4 * RodN + 4 + 3}[5];

// PML Y surfaces
Tmp~{0}~{0} = Boundary{ Volume{Pml~{1}[0 * (RodN + 2) + RodN + 0]}; };
Tmp~{0}~{1} = Boundary{ Volume{Pml~{1}[1 * (RodN + 2) + RodN + 0]}; };
Tmp~{1}~{0} = Boundary{ Volume{Pml~{1}[0 * (RodN + 2) + RodN + 1]}; };
Tmp~{1}~{1} = Boundary{ Volume{Pml~{1}[1 * (RodN + 2) + RodN + 1]}; };

Bnd~{2} = {};
Bnd~{3} = {};

Bnd~{2}[0] = Tmp~{0}~{0}[4];
Bnd~{2}[1] = Tmp~{0}~{1}[4];

Bnd~{3}[0] = Tmp~{1}~{0}[4];
Bnd~{3}[1] = Tmp~{1}~{1}[4];

// PML Z surfaces
Tmp~{0}~{0} = Boundary{ Volume{Pml~{0}[0 * (RodN + 2) + RodN + 0]}; };
Tmp~{0}~{1} = Boundary{ Volume{Pml~{0}[1 * (RodN + 2) + RodN + 0]}; };
Tmp~{1}~{0} = Boundary{ Volume{Pml~{0}[0 * (RodN + 2) + RodN + 1]}; };
Tmp~{1}~{1} = Boundary{ Volume{Pml~{0}[1 * (RodN + 2) + RodN + 1]}; };

Bnd~{4} = {};
Bnd~{5} = {};

Bnd~{4}[0] = Tmp~{0}~{0}[4];
Bnd~{4}[1] = Tmp~{0}~{1}[4];
Bnd~{5}[0] = Tmp~{1}~{0}[4];
Bnd~{5}[1] = Tmp~{1}~{1}[4];

// PML X //
///////////
// Pml~{3:6}[] = [... Pml(X) Pml(XY) Pml(XZ) Pml(XYZ)]

// X- //
Tmp~{0} = {};
Tmp~{0} = Extrude{-PmlSize, 0, 0}{ Surface{EndSurf~{0}[1],
                                           Bnd~{2}[], Bnd~{4}[], Bnd~{0}[]}; };

Pml~{3}[0] = Tmp~{0}[0 * 2 + 1]; // X
Pml~{4}[0] = Tmp~{0}[1 * 2 + 1]; // XY
Pml~{4}[1] = Tmp~{0}[2 * 2 + 1]; // XY
Pml~{5}[0] = Tmp~{0}[3 * 2 + 1]; // XZ
Pml~{5}[1] = Tmp~{0}[4 * 2 + 1]; // XZ
Pml~{6}[0] = Tmp~{0}[5 * 2 + 1]; // XYZ
Pml~{6}[1] = Tmp~{0}[6 * 2 + 1]; // XYZ
Pml~{6}[2] = Tmp~{0}[7 * 2 + 1]; // XYZ
Pml~{6}[3] = Tmp~{0}[8 * 2 + 1]; // XYZ

// X+ //
Tmp~{0} = {};
Tmp~{0} = Extrude{+PmlSize, 0, 0}{ Surface{EndSurf~{1}[1],
                                           Bnd~{3}[], Bnd~{5}[], Bnd~{1}[]}; };
Pml~{3}[1] = Tmp~{0}[0 * 2 + 1]; // X
Pml~{4}[2] = Tmp~{0}[1 * 2 + 1]; // XY
Pml~{4}[3] = Tmp~{0}[2 * 2 + 1]; // XY
Pml~{5}[2] = Tmp~{0}[3 * 2 + 1]; // XZ
Pml~{5}[3] = Tmp~{0}[4 * 2 + 1]; // XZ
Pml~{6}[4] = Tmp~{0}[5 * 2 + 1]; // XYZ
Pml~{6}[5] = Tmp~{0}[6 * 2 + 1]; // XYZ
Pml~{6}[6] = Tmp~{0}[7 * 2 + 1]; // XYZ
Pml~{6}[7] = Tmp~{0}[8 * 2 + 1]; // XYZ

// Clear useless variables //
/////////////////////////////
For i In {1:RodN}
  VolSurf~{i} = {};
EndFor

EndSurf~{0} = {};
EndSurf~{1} = {};

For i In {0:4 * (RodN + 2) - 1}
  Tmp~{0}~{i} = {};
EndFor
Tmp~{1}~{0} = {};
Tmp~{1}~{1} = {};

Bnd~{0} = {};
Bnd~{1} = {};
Bnd~{2} = {};
Bnd~{3} = {};
Bnd~{4} = {};
Bnd~{5} = {};

Tmp~{0} = {};
