// Getting cells and PML boundaries for DDM //
//////////////////////////////////////////////

// DdmRod~{i}~{j}[k]
//  * k denotes the kth rod volume
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume
//  * i = 0: 'Air'              part of kth volume
//      = 1: 'PmlZ'  down       part of kth volume
//      = 2: 'PmlZ'    up       part of kth volume
//      = 3: 'PmlY'  rear       part of kth volume
//      = 4: 'PmlY'  front      part of kth volume
//      = 5: 'PmlYZ' rear  down part of kth volume
//      = 6: 'PmlYZ' rear    up part of kth volume
//      = 7: 'PmlYZ' front down part of kth volume
//      = 8: 'PmlYZ' front   up part of kth volume

// DdmEnd~{i}~{j}[k]
//  * k denotes the kth end volume
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume
//  * i = 0: 'Air'              part of kth volume
//      = 1: 'PmlZ'  down       part of kth volume
//      = 2: 'PmlZ'    up       part of kth volume
//      = 3: 'PmlY'  rear       part of kth volume
//      = 4: 'PmlY'  front      part of kth volume
//      = 5: 'PmlYZ' rear  down part of kth volume
//      = 6: 'PmlYZ' rear    up part of kth volume
//      = 7: 'PmlYZ' front down part of kth volume
//      = 8: 'PmlYZ' front   up part of kth volume


// Rod //
For i In {0:RodN - 1}
  VolSurf~{i} = Boundary{ Volume{Vol~{i + 1}[1]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{0}~{0}[i] = VolSurf~{i}[3]; // Right
  DdmRod~{0}~{1}[i] = VolSurf~{i}[5]; // Left
EndFor

// Ends //
For i In {0:EndN - 1}
  EndSurf~{i}~{0} = Boundary{ Volume{End~{i}[0]}; };
  EndSurf~{i}~{1} = Boundary{ Volume{End~{i}[1]}; };
EndFor

For i In {0:EndN - 1}
  DdmEnd~{0}~{0}[i + 0 * EndN] = EndSurf~{i}~{0}[0]; // Right
  DdmEnd~{0}~{0}[i + 1 * EndN] = EndSurf~{i}~{1}[1]; // Right

  DdmEnd~{0}~{1}[i + 0 * EndN] = EndSurf~{i}~{0}[1]; // Left
  DdmEnd~{0}~{1}[i + 1 * EndN] = EndSurf~{i}~{1}[0]; // Left
EndFor

// Pml Z //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = Boundary{ Volume{Pml~{0}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{1}~{0}[i] = PmlZSurf~{i}[3]; // Right
  DdmRod~{1}~{1}[i] = PmlZSurf~{i}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{1}~{0}[i + 0 * EndN] = PmlZSurf~{i + 1 * RodN + 0 * EndN}[2]; // Right
  DdmEnd~{1}~{0}[i + 1 * EndN] = PmlZSurf~{i + 1 * RodN + 1 * EndN}[4]; // Right

  DdmEnd~{1}~{1}[i + 0 * EndN] = PmlZSurf~{i + 1 * RodN + 0 * EndN}[4]; // Left
  DdmEnd~{1}~{1}[i + 1 * EndN] = PmlZSurf~{i + 1 * RodN + 1 * EndN}[2]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{2}~{0}[i] = PmlZSurf~{i + RodN + 2 * EndN}[3]; // Right
  DdmRod~{2}~{1}[i] = PmlZSurf~{i + RodN + 2 * EndN}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{2}~{0}[i + 0 * EndN] = PmlZSurf~{i + 2 * RodN + 2 * EndN}[2]; // Right
  DdmEnd~{2}~{0}[i + 1 * EndN] = PmlZSurf~{i + 2 * RodN + 3 * EndN}[4]; // Right

  DdmEnd~{2}~{1}[i + 0 * EndN] = PmlZSurf~{i + 2 * RodN + 2 * EndN}[4]; // Left
  DdmEnd~{2}~{1}[i + 1 * EndN] = PmlZSurf~{i + 2 * RodN + 3 * EndN}[2]; // Left
EndFor

// Pml Y //
For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlYSurf~{i} = Boundary{ Volume{Pml~{1}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{3}~{0}[i] = PmlYSurf~{i}[3]; // Right
  DdmRod~{3}~{1}[i] = PmlYSurf~{i}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{3}~{0}[i + 0 * EndN] = PmlYSurf~{i + 1 * RodN + 0 * EndN}[2]; // Right
  DdmEnd~{3}~{0}[i + 1 * EndN] = PmlYSurf~{i + 1 * RodN + 1 * EndN}[4]; // Right

  DdmEnd~{3}~{1}[i + 0 * EndN] = PmlYSurf~{i + 1 * RodN + 0 * EndN}[4]; // Left
  DdmEnd~{3}~{1}[i + 1 * EndN] = PmlYSurf~{i + 1 * RodN + 1 * EndN}[2]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{4}~{0}[i] = PmlYSurf~{i + RodN + 2 * EndN}[5]; // Right
  DdmRod~{4}~{1}[i] = PmlYSurf~{i + RodN + 2 * EndN}[3]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{4}~{0}[i + 0 * EndN] = PmlYSurf~{i + 2 * RodN + 2 * EndN}[2]; // Right
  DdmEnd~{4}~{0}[i + 1 * EndN] = PmlYSurf~{i + 2 * RodN + 3 * EndN}[4]; // Right

  DdmEnd~{4}~{1}[i + 0 * EndN] = PmlYSurf~{i + 2 * RodN + 2 * EndN}[4]; // Left
  DdmEnd~{4}~{1}[i + 1 * EndN] = PmlYSurf~{i + 2 * RodN + 3 * EndN}[2]; // Left
EndFor

// Pml YZ //
For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = Boundary{ Volume{Pml~{2}[i]}; };
EndFor

For i In {0:RodN - 1}
  DdmRod~{5}~{0}[i] = PmlYZSurf~{i}[3]; // Right
  DdmRod~{5}~{1}[i] = PmlYZSurf~{i}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{5}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 1 * RodN + 0 * EndN}[5]; // Right
  DdmEnd~{5}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 1 * RodN + 1 * EndN}[5]; // Right

  DdmEnd~{5}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 1 * RodN + 0 * EndN}[3]; // Left
  DdmEnd~{5}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 1 * RodN + 1 * EndN}[3]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{6}~{0}[i] = PmlYZSurf~{i + 1 * RodN + 2 * EndN}[3]; // Right
  DdmRod~{6}~{1}[i] = PmlYZSurf~{i + 1 * RodN + 2 * EndN}[5]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{6}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 2 * RodN + 2 * EndN}[3]; // Right
  DdmEnd~{6}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 2 * RodN + 3 * EndN}[3]; // Right

  DdmEnd~{6}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 2 * RodN + 2 * EndN}[5]; // Left
  DdmEnd~{6}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 2 * RodN + 3 * EndN}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{7}~{0}[i] = PmlYZSurf~{i + 2 * RodN + 4 * EndN}[5]; // Right
  DdmRod~{7}~{1}[i] = PmlYZSurf~{i + 2 * RodN + 4 * EndN}[3]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{7}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 3 * RodN + 4 * EndN}[3]; // Right
  DdmEnd~{7}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 3 * RodN + 5 * EndN}[3]; // Right

  DdmEnd~{7}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 3 * RodN + 4 * EndN}[5]; // Left
  DdmEnd~{7}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 3 * RodN + 5 * EndN}[5]; // Left
EndFor

For i In {0:RodN - 1}
  DdmRod~{8}~{0}[i] = PmlYZSurf~{i + 3 * RodN + 6 * EndN}[5]; // Right
  DdmRod~{8}~{1}[i] = PmlYZSurf~{i + 3 * RodN + 6 * EndN}[3]; // Left
EndFor

For i In {0:EndN - 1}
  DdmEnd~{8}~{0}[i + 0 * EndN] = PmlYZSurf~{i + 4 * RodN + 6 * EndN}[5]; // Right
  DdmEnd~{8}~{0}[i + 1 * EndN] = PmlYZSurf~{i + 4 * RodN + 7 * EndN}[5]; // Right

  DdmEnd~{8}~{1}[i + 0 * EndN] = PmlYZSurf~{i + 4 * RodN + 6 * EndN}[3]; // Left
  DdmEnd~{8}~{1}[i + 1 * EndN] = PmlYZSurf~{i + 4 * RodN + 7 * EndN}[3]; // Left
EndFor

For i In {0:5}
  //Printf("%f: %f", i, PmlYZSurf~{0 + 4 * RodN + 7 * EndN}[i]);
EndFor

For i In {0:2*EndN - 1}
  //Printf("%f: %f", i, DdmEnd~{8}~{0}[i]);
EndFor


// Clear //
For i In {0:RodN - 1}
  VolSurf~{i} = {};
EndFor

For i In {0:EndN - 1}
  EndSurf~{i}~{0} = {};
  EndSurf~{i}~{1} = {};
EndFor

For i In {0:(2 * RodN + 4 * EndN) - 1}
  PmlZSurf~{i} = {};
  PmlYSurf~{i} = {};
EndFor

For i In {0:(4 * RodN + 8 * EndN) - 1}
  PmlYZSurf~{i} = {};
EndFor
