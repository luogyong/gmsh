// User Data //
///////////////
Include "boubouchons.dat";

// Create geomtetry //
//////////////////////
Include "cell.geo";
Include "grid.geo";
Include "ends.geo";
Include  "src.geo";
Include  "pml.geo";

// Center geomtry //
////////////////////
Include "center.geo";

// Surfaces for DDM //
//////////////////////
If(UseDDM == 1)
  Include "ddm.geo";
EndIf

// Physicals //
///////////////
//If(UseDDM == 0)
Include "physicalNoDdm.geo";
//EndIf

If(UseDDM == 1)
  //Include "physicalDdm.geo";
EndIf

Mesh 3;

Include "partition.geo";

Save "boubouchons.msh";
