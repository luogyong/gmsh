// PML X //
///////////
// Pml~{3:6}[] = [... Pml(X) Pml(XY) Pml(XZ) Pml(XYZ)]

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

// PML ZY surfaces //
For i In {0:(4 * (RodN + 2 * EndN) - 1)}
  Tmp~{i} = Boundary{ Volume{Pml~{2}[i]}; };
EndFor

Bnd~{0} = {};
Bnd~{1} = {};

Bnd~{0}[0] = Tmp~{1 * RodN + 1 * EndN - 1}[3];
Bnd~{0}[1] = Tmp~{2 * RodN + 3 * EndN - 1}[5];
Bnd~{0}[2] = Tmp~{3 * RodN + 5 * EndN - 1}[5];
Bnd~{0}[3] = Tmp~{4 * RodN + 7 * EndN - 1}[3];

Bnd~{1}[0] = Tmp~{1 * RodN + 2 * EndN - 1}[5];
Bnd~{1}[1] = Tmp~{2 * RodN + 4 * EndN - 1}[3];
Bnd~{1}[2] = Tmp~{3 * RodN + 6 * EndN - 1}[3];
Bnd~{1}[3] = Tmp~{4 * RodN + 8 * EndN - 1}[5];

// PML Y surfaces //
Tmp~{0} = Boundary{ Volume{Pml~{1}[1 * RodN + 1 * EndN - 1]}; };
Tmp~{1} = Boundary{ Volume{Pml~{1}[1 * RodN + 2 * EndN - 1]}; };
Tmp~{2} = Boundary{ Volume{Pml~{1}[2 * RodN + 3 * EndN - 1]}; };
Tmp~{3} = Boundary{ Volume{Pml~{1}[2 * RodN + 4 * EndN - 1]}; };

Bnd~{2} = {};
Bnd~{3} = {};

Bnd~{2}[0] = Tmp~{0}[4];
Bnd~{2}[1] = Tmp~{2}[4];

Bnd~{3}[0] = Tmp~{1}[4];
Bnd~{3}[1] = Tmp~{3}[4];

// PML Z surfaces //
Tmp~{0} = Boundary{ Volume{Pml~{0}[1 * RodN + 1 * EndN - 1]}; };
Tmp~{1} = Boundary{ Volume{Pml~{0}[1 * RodN + 2 * EndN - 1]}; };
Tmp~{2} = Boundary{ Volume{Pml~{0}[2 * RodN + 3 * EndN - 1]}; };
Tmp~{3} = Boundary{ Volume{Pml~{0}[2 * RodN + 4 * EndN - 1]}; };

Bnd~{4} = {};
Bnd~{5} = {};

Bnd~{4}[0] = Tmp~{0}[4];
Bnd~{4}[1] = Tmp~{2}[4];

Bnd~{5}[0] = Tmp~{1}[4];
Bnd~{5}[1] = Tmp~{3}[4];


// Number of extrusion //
PmlN = Ceil(PmlSize / RodP);

// There is at least on extrusion //
If(PmlN == 0)
  PmlN = 1;
EndIf


// X- //
// Extrude once
Tmp~{0} = {};
Tmp~{0} = Extrude{-RodP, 0, 0}{
  Surface{EndSurf~{EndN - 1}~{0}[1], Bnd~{2}[], Bnd~{4}[], Bnd~{0}[]};
};

Pml~{3}[0] = Tmp~{0}[0 * 2 + 1]; // X
Pml~{4}[0] = Tmp~{0}[1 * 2 + 1]; // XY
Pml~{4}[1] = Tmp~{0}[2 * 2 + 1]; // XY
Pml~{5}[0] = Tmp~{0}[3 * 2 + 1]; // XZ
Pml~{5}[1] = Tmp~{0}[4 * 2 + 1]; // XZ
Pml~{6}[0] = Tmp~{0}[5 * 2 + 1]; // XYZ
Pml~{6}[1] = Tmp~{0}[6 * 2 + 1]; // XYZ
Pml~{6}[2] = Tmp~{0}[7 * 2 + 1]; // XYZ
Pml~{6}[3] = Tmp~{0}[8 * 2 + 1]; // XYZ

// Continue extruding
For i In {1:PmlN - 1}
  Tmp~{0} = Boundary{ Volume{Pml~{3}[1 * (i - 1) + 0]}; };
  Tmp~{1} = Boundary{ Volume{Pml~{4}[2 * (i - 1) + 0]}; };
  Tmp~{2} = Boundary{ Volume{Pml~{4}[2 * (i - 1) + 1]}; };
  Tmp~{3} = Boundary{ Volume{Pml~{5}[2 * (i - 1) + 0]}; };
  Tmp~{4} = Boundary{ Volume{Pml~{5}[2 * (i - 1) + 1]}; };
  Tmp~{5} = Boundary{ Volume{Pml~{6}[4 * (i - 1) + 0]}; };
  Tmp~{6} = Boundary{ Volume{Pml~{6}[4 * (i - 1) + 1]}; };
  Tmp~{7} = Boundary{ Volume{Pml~{6}[4 * (i - 1) + 2]}; };
  Tmp~{8} = Boundary{ Volume{Pml~{6}[4 * (i - 1) + 3]}; };

  Tmp~{0} = Extrude{-RodP, 0, 0}{
    Surface{Tmp~{0}[1],
            Tmp~{1}[1], Tmp~{2}[1],
            Tmp~{3}[1], Tmp~{4}[1],
            Tmp~{5}[1], Tmp~{6}[1], Tmp~{7}[1], Tmp~{8}[1]};
  };

  Pml~{3}[1 * i + 0] = Tmp~{0}[0 * 2 + 1]; // X
  Pml~{4}[2 * i + 0] = Tmp~{0}[1 * 2 + 1]; // XY
  Pml~{4}[2 * i + 1] = Tmp~{0}[2 * 2 + 1]; // XY
  Pml~{5}[2 * i + 0] = Tmp~{0}[3 * 2 + 1]; // XZ
  Pml~{5}[2 * i + 1] = Tmp~{0}[4 * 2 + 1]; // XZ
  Pml~{6}[4 * i + 0] = Tmp~{0}[5 * 2 + 1]; // XYZ
  Pml~{6}[4 * i + 1] = Tmp~{0}[6 * 2 + 1]; // XYZ
  Pml~{6}[4 * i + 2] = Tmp~{0}[7 * 2 + 1]; // XYZ
  Pml~{6}[4 * i + 3] = Tmp~{0}[8 * 2 + 1]; // XYZ
EndFor


// X+ //
// Extrude once
Tmp~{0} = {};
Tmp~{0} = Extrude{+RodP, 0, 0}{
  Surface{EndSurf~{EndN - 1}~{1}[1], Bnd~{3}[], Bnd~{5}[], Bnd~{1}[]};
};

Pml~{3}[1 * PmlN + 0] = Tmp~{0}[0 * 2 + 1]; // X
Pml~{4}[2 * PmlN + 0] = Tmp~{0}[1 * 2 + 1]; // XY
Pml~{4}[2 * PmlN + 1] = Tmp~{0}[2 * 2 + 1]; // XY
Pml~{5}[2 * PmlN + 0] = Tmp~{0}[3 * 2 + 1]; // XZ
Pml~{5}[2 * PmlN + 1] = Tmp~{0}[4 * 2 + 1]; // XZ
Pml~{6}[4 * PmlN + 0] = Tmp~{0}[5 * 2 + 1]; // XYZ
Pml~{6}[4 * PmlN + 1] = Tmp~{0}[6 * 2 + 1]; // XYZ
Pml~{6}[4 * PmlN + 2] = Tmp~{0}[7 * 2 + 1]; // XYZ
Pml~{6}[4 * PmlN + 3] = Tmp~{0}[8 * 2 + 1]; // XYZ

// Continue extruding
For i In {1:PmlN - 1}
  Tmp~{0} = Boundary{ Volume{Pml~{3}[1 * PmlN + 1 * (i - 1) + 0]}; };
  Tmp~{1} = Boundary{ Volume{Pml~{4}[2 * PmlN + 2 * (i - 1) + 0]}; };
  Tmp~{2} = Boundary{ Volume{Pml~{4}[2 * PmlN + 2 * (i - 1) + 1]}; };
  Tmp~{3} = Boundary{ Volume{Pml~{5}[2 * PmlN + 2 * (i - 1) + 0]}; };
  Tmp~{4} = Boundary{ Volume{Pml~{5}[2 * PmlN + 2 * (i - 1) + 1]}; };
  Tmp~{5} = Boundary{ Volume{Pml~{6}[4 * PmlN + 4 * (i - 1) + 0]}; };
  Tmp~{6} = Boundary{ Volume{Pml~{6}[4 * PmlN + 4 * (i - 1) + 1]}; };
  Tmp~{7} = Boundary{ Volume{Pml~{6}[4 * PmlN + 4 * (i - 1) + 2]}; };
  Tmp~{8} = Boundary{ Volume{Pml~{6}[4 * PmlN + 4 * (i - 1) + 3]}; };

  Tmp~{0} = Extrude{+RodP, 0, 0}{
    Surface{Tmp~{0}[1],
            Tmp~{1}[1], Tmp~{2}[1],
            Tmp~{3}[1], Tmp~{4}[1],
            Tmp~{5}[1], Tmp~{6}[1], Tmp~{7}[1], Tmp~{8}[1]};
  };

  Pml~{3}[1 * PmlN + 1 * i + 0] = Tmp~{0}[0 * 2 + 1]; // X
  Pml~{4}[2 * PmlN + 2 * i + 0] = Tmp~{0}[1 * 2 + 1]; // XY
  Pml~{4}[2 * PmlN + 2 * i + 1] = Tmp~{0}[2 * 2 + 1]; // XY
  Pml~{5}[2 * PmlN + 2 * i + 0] = Tmp~{0}[3 * 2 + 1]; // XZ
  Pml~{5}[2 * PmlN + 2 * i + 1] = Tmp~{0}[4 * 2 + 1]; // XZ
  Pml~{6}[4 * PmlN + 4 * i + 0] = Tmp~{0}[5 * 2 + 1]; // XYZ
  Pml~{6}[4 * PmlN + 4 * i + 1] = Tmp~{0}[6 * 2 + 1]; // XYZ
  Pml~{6}[4 * PmlN + 4 * i + 2] = Tmp~{0}[7 * 2 + 1]; // XYZ
  Pml~{6}[4 * PmlN + 4 * i + 3] = Tmp~{0}[8 * 2 + 1]; // XYZ
EndFor

// Clear useless variables //
For i In {1:RodN}
  VolSurf~{i} = {};
EndFor

For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor

For i In {0:7}
  Tmp~{i} = {};
EndFor

For i In {0:(4 * (RodN + 2 * EndN) - 1)}
  Tmp~{i} = {};
EndFor

Bnd~{0} = {};
Bnd~{1} = {};
Bnd~{2} = {};
Bnd~{3} = {};
Bnd~{4} = {};
Bnd~{5} = {};
