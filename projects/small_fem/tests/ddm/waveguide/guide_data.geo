// User constant //
DefineConstant[
  LX       = {2,               Name "Input/00Geometry/00X dimension"},
  LY       = {1,               Name "Input/00Geometry/01Y dimension"},

  K        = {10,              Name "Input/01Mesh/00Wavenumber"},
  LAMBDA   = {2* Pi / K,       Name "Input/01Mesh/01Wavelength",   ReadOnly 1},
  N_LAMBDA = {10,              Name "Input/01Mesh/02Points per wavelength"},
  LC       = {LAMBDA/N_LAMBDA, Name "Input/01Mesh/03Mesh density", ReadOnly 1},
  STRUCT   = {1,               Name "Input/01Mesh/04Structured", Choices {0,1}},

  NDOM     = {2,               Name "Input/02DDM/00Number of subdomain"}
];

If(DIM == 2)
  DefineConstant[
    LZ     = {1, Name "Input/00Geometry/02Z dimension", Visible 0}
  ];
EndIf
If(DIM == 3)
  DefineConstant[
    LZ     = {1, Name "Input/00Geometry/02Z dimension"}
  ];
EndIf

// base msh filename
MSH_BASE_NAME = "mesh_";
// directory for output files
DIR = "out/";
// prefix for (split) mesh files (one for each partition)
MSH_NAME = StrCat(DIR, MSH_BASE_NAME);
