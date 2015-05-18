// User Data //
///////////////
Include "boubouchons.dat";

// Create geomtetry //
//////////////////////
Include "cell.geo";
Include "grid.geo";
Include "ends.geo";
//Include  "srcLine.geo";
Include  "srcSphere.geo";
Include  "pml.geo";

// Center geomtry //
////////////////////
Include "center.geo";

// Surfaces for DDM //
//////////////////////
//Include "ddm.geo";
Include "ddmTwo.geo";

// Physicals //
///////////////
//Include "physical.geo";
Include "physicalTwo.geo";

// Mesh (with partitions) //
////////////////////////////
If(DoMesh == 1)
  Mesh    3;
  //Include "partition.geo";
  Include "partitionTwo.geo";
  Save    "boubouchons.msh";
EndIf
