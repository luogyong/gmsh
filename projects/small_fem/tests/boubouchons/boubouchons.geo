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

// Number of cells //
/////////////////////
CellN = (2*PmlN + 2*EndN + RodN);

// Sub-domains size for DDM (in cells) //
/////////////////////////////////////////
DdmN  = Ceil(CellN / DomN);

If(DdmN != Floor(CellN / DomN))
  Error("Size of sub-domains is not integer (%g): abort!", CellN / DomN);
  Abort;
EndIf

// Surfaces for DDM //
//////////////////////
Include "ddm.geo";

// Range //
///////////
For i In {0:(DomN - 1)}
  Range~{i} = {i * DdmN, ((i + 1) * DdmN - 1)};
EndFor

// Reorder stuffs //
////////////////////
Include "reorder.geo";

// Physicals //
///////////////
Include "physical.geo";

// Stats //
///////////
Include "stats.geo";

// Mesh (with partitions) //
////////////////////////////
If(DoMesh == 1)
  Mesh.Algorithm   = 6;
  Mesh.Algorithm3D = 1;
  Mesh.Optimize    = 1;

  If(MeshOrder > 1)
    Mesh.HighOrderOptimize = 1;
    Mesh.ElementOrder      = MeshOrder;
  EndIf

  Mesh    3;
  Include "partition.geo";
  Save    "boubouchons.msh";
EndIf
