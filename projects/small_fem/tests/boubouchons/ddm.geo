// Getting boundaries for DDM //
////////////////////////////////

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

// DdmSurf~{j}[k]
//  * k denotes a unique volume in the geometry
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume

Include "ddmRod.geo";
Include "ddmEnd.geo";
Include "ddmPml.geo";

For i In {0:1}
  For j In {0:RodN - 1}
    DdmSortRod~{j}~{i} = {
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

  For j In {0:(2 * EndN - 1)}
    DdmSortEnd~{j}~{i} = {
      DdmEnd~{0}~{i}[j],
      DdmEnd~{1}~{i}[j],
      DdmEnd~{2}~{i}[j],
      DdmEnd~{3}~{i}[j],
      DdmEnd~{4}~{i}[j],
      DdmEnd~{5}~{i}[j],
      DdmEnd~{6}~{i}[j],
      DdmEnd~{7}~{i}[j],
      DdmEnd~{8}~{i}[j]
    };
  EndFor

  For j In {0:(2 * PmlN - 1)}
    DdmSortPml~{j}~{i} = {
      DdmPml~{0}~{i}[j],
      DdmPml~{1}~{i}[j],
      DdmPml~{2}~{i}[j],
      DdmPml~{3}~{i}[j],
      DdmPml~{4}~{i}[j],
      DdmPml~{5}~{i}[j],
      DdmPml~{6}~{i}[j],
      DdmPml~{7}~{i}[j],
      DdmPml~{8}~{i}[j]
    };
  EndFor
EndFor
