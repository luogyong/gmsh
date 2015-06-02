// Getting boundaries for DDM //
////////////////////////////////

// DdmRod~{i}~{j}[k]
//  * k denotes the kth Rod volume
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume
//  * i = 0: 'Air'              part of kth volume

// DdmEnd~{i}~{j}[k]
//  * k denotes the kth End volume
//  * j = 0: for right boundary (along x axis) of kth volume
//      = 1: for  left boundary (along x axis) of kth volume
//  * i = 0: 'Air'              part of kth volume

// DdmBoundary~{k}~{j}[]
//  * k denotes the kth Cell volume
//  * j = 0: for right boundary of kth volume
//      = 1: for  left boundary of kth volume
//  * Surface indices are in list stored at ~{k}~{j}

// Ddm boundaries //
Include "ddmRod.geo";
Include "ddmEnd.geo";

// Agregate //
For i In {0:1}
  For j In {0:(EndN - 1)}
    DdmBoundary~{j}~{i} = DdmEnd~{0}~{i}[EndN - 1 - j];
  EndFor

  For j In {0:(RodN - 1)}
    DdmBoundary~{j + EndN}~{i} = DdmRod~{0}~{i}[j];
  EndFor

  For j In {0:(EndN - 1)}
    DdmBoundary~{j + EndN + RodN}~{i} = DdmEnd~{0}~{i}[j + EndN];
  EndFor
EndFor
