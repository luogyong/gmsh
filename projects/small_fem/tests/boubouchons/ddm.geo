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
  DdmSurf~{i} = {
    DdmRod~{0}~{i}[],
    DdmRod~{1}~{i}[],
    DdmRod~{2}~{i}[],
    DdmRod~{3}~{i}[],
    DdmRod~{4}~{i}[],
    DdmRod~{5}~{i}[],
    DdmRod~{6}~{i}[],
    DdmRod~{7}~{i}[],
    DdmRod~{8}~{i}[],

    DdmEnd~{0}~{i}[],
    DdmEnd~{1}~{i}[],
    DdmEnd~{2}~{i}[],
    DdmEnd~{3}~{i}[],
    DdmEnd~{4}~{i}[],
    DdmEnd~{5}~{i}[],
    DdmEnd~{6}~{i}[],
    DdmEnd~{7}~{i}[],
    DdmEnd~{8}~{i}[],

    DdmPml~{0}~{i}[],
    DdmPml~{1}~{i}[],
    DdmPml~{2}~{i}[],
    DdmPml~{3}~{i}[],
    DdmPml~{4}~{i}[],
    DdmPml~{5}~{i}[],
    DdmPml~{6}~{i}[],
    DdmPml~{7}~{i}[],
    DdmPml~{8}~{i}[]
  };
EndFor
