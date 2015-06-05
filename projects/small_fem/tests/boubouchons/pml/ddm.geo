// Getting boundaries for DDM //
////////////////////////////////

// DdmRod~{i}~{j}[k]
//  * k denotes the kth Rod volume
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
//  * k denotes the kth End volume
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

// DdmPml~{i}~{j}[k]
//  * k denotes the kth Pml X volume
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume
//  * i = 0: 'PmlX'              part of kth volume
//      = 1: 'PmlXZ'  down       part of kth volume
//      = 2: 'PmlXZ'    up       part of kth volume
//      = 3: 'PmlXY'  rear       part of kth volume
//      = 4: 'PmlXY'  front      part of kth volume
//      = 5: 'PmlXYZ' rear  down part of kth volume
//      = 6: 'PmlXYZ' rear    up part of kth volume
//      = 7: 'PmlXYZ' front down part of kth volume
//      = 8: 'PmlXYZ' front   up part of kth volume

// DdmBoundary~{k}~{j}[]
//  * k denotes the kth Cell volume
//  * j = 0: for right boundary of kth volume
//      = 1: for  left boundary of kth volume
//  * Surface indices are in list stored at ~{k}~{j}

// Ddm boundaries //
Include "ddmRod.geo";
Include "ddmEnd.geo";
Include "ddmPml.geo";

// Agregate //
For i In {0:1}
  For j In {0:(PmlN - 1)}
    DdmBoundary~{j}~{i} = {
      DdmPml~{0}~{i}[PmlN - 1 - j],
      DdmPml~{1}~{i}[PmlN - 1 - j],
      DdmPml~{2}~{i}[PmlN - 1 - j],
      DdmPml~{3}~{i}[PmlN - 1 - j],
      DdmPml~{4}~{i}[PmlN - 1 - j],
      DdmPml~{5}~{i}[PmlN - 1 - j],
      DdmPml~{6}~{i}[PmlN - 1 - j],
      DdmPml~{7}~{i}[PmlN - 1 - j],
      DdmPml~{8}~{i}[PmlN - 1 - j]
    };
  EndFor

  For j In {0:(EndN - 1)}
    DdmBoundary~{j + PmlN}~{i} = {
      DdmEnd~{0}~{i}[EndN - 1 - j],
      DdmEnd~{1}~{i}[EndN - 1 - j],
      DdmEnd~{2}~{i}[EndN - 1 - j],
      DdmEnd~{3}~{i}[EndN - 1 - j],
      DdmEnd~{4}~{i}[EndN - 1 - j],
      DdmEnd~{5}~{i}[EndN - 1 - j],
      DdmEnd~{6}~{i}[EndN - 1 - j],
      DdmEnd~{7}~{i}[EndN - 1 - j],
      DdmEnd~{8}~{i}[EndN - 1 - j]
    };
  EndFor

  For j In {0:(RodN - 1)}
    DdmBoundary~{j + PmlN + EndN}~{i} = {
      DdmRod~{0}~{i}[j],
      DdmRod~{1}~{i}[j],
      DdmRod~{2}~{i}[j],
      DdmRod~{3}~{i}[j],
      DdmRod~{4}~{i}[j],
      DdmRod~{5}~{i}[j],
      DdmRod~{6}~{i}[j],
      DdmRod~{7}~{i}[j],
      DdmRod~{8}~{i}[j]
    };
  EndFor

  For j In {0:(EndN - 1)}
    DdmBoundary~{j + PmlN + EndN + RodN}~{i} = {
      DdmEnd~{0}~{i}[j + EndN],
      DdmEnd~{1}~{i}[j + EndN],
      DdmEnd~{2}~{i}[j + EndN],
      DdmEnd~{3}~{i}[j + EndN],
      DdmEnd~{4}~{i}[j + EndN],
      DdmEnd~{5}~{i}[j + EndN],
      DdmEnd~{6}~{i}[j + EndN],
      DdmEnd~{7}~{i}[j + EndN],
      DdmEnd~{8}~{i}[j + EndN]
    };
  EndFor

  For j In {0:(PmlN - 1)}
    DdmBoundary~{j + PmlN + EndN + RodN + EndN}~{i} = {
      DdmPml~{0}~{i}[j + PmlN],
      DdmPml~{1}~{i}[j + PmlN],
      DdmPml~{2}~{i}[j + PmlN],
      DdmPml~{3}~{i}[j + PmlN],
      DdmPml~{4}~{i}[j + PmlN],
      DdmPml~{5}~{i}[j + PmlN],
      DdmPml~{6}~{i}[j + PmlN],
      DdmPml~{7}~{i}[j + PmlN],
      DdmPml~{8}~{i}[j + PmlN]
    };
  EndFor
EndFor
