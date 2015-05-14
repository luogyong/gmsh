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
Include "ddm.geo";

// Physicals //
///////////////
Include "physical.geo";

// Mesh (with partitions) //
////////////////////////////
If(DoMesh == 1)
  Mesh    3;
  Include "partition.geo";
  Save    "boubouchons.msh";
EndIf
