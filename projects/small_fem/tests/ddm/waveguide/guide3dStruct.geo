// User constant //
DIM = 3;
Include "guide_data.geo";

// Works from NDOM == 2 //
If(NDOM != 2)
  Error("NDOM must be equal to 2");
  Abort;
EndIf

// Geometry (1D) //
// Points
Point(1) = {0,  0, 0, LC};
Point(2) = {0, LY, 0, LC};

// Line
Line(1) = {1, 2};
Transfinite Line {1} = (Ceil(LY / LAMBDA * N_LAMBDA)) Using Progression 1;

// Extrusion 2D //
Extrude {LX / NDOM, 0, 0} {
  Line{1};
  Layers{Ceil((LX / NDOM) / LAMBDA * N_LAMBDA)};
}

Extrude {LX / NDOM, 0, 0} {
  Line{2};
  Layers{Ceil((LX / NDOM) / LAMBDA * N_LAMBDA)};
}

// Extrusion 3D //
Extrude {0, 0, LZ} {
  Surface{5, 9};
  Layers{Ceil(LZ / LAMBDA * N_LAMBDA)};
}

// Physicals //
Physical Volume (1) = {1};
Physical Volume (2) = {2};
Physical Surface(3) = {18};
Physical Surface(4) = {26};
Physical Surface(5) = {48};
Physical Surface(6) = {5, 22, 30, 31, 9, 44, 52, 53};

If(StrCmp(OnelabAction, "check")) // only mesh if not in onelab check mode
  Printf("Meshing waveguide full...");
  Mesh 3 ;
  CreateDir Str(DIR);
  Save StrCat(MSH_NAME, "all.msh");
  Save "guide3d.msh";
  Printf("Done.");
EndIf
