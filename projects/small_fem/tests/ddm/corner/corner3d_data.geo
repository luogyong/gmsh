DefineConstant[
  // Analysis type
  ANALYSIS = {1, Name "Input/00Type of analysis", ReadOnly 1,
    Choices {0="Helmholtz", 1="Maxwell"}},

  // Dimensions
  LX = {2, Name "Input/01Geometry/00X length", ReadOnly 1},
  LY = {1, Name "Input/01Geometry/01Y length", ReadOnly 1},
  LZ = {1, Name "Input/01Geometry/02Z length", ReadOnly 1},

  // Frequency
  WAVENUMBER = {5, Name "Input/01Physics/00Wavenumber"},
  LAMBDA = {2*Pi/WAVENUMBER, Name "Input/01Physics/01Wavelength", ReadOnly 1},
  MODE_M = {1, Name "Input/01Physics/02Y Wavenumber", ReadOnly 1},
  MODE_N = {1, Name "Input/01Physics/03Z Wavenumber", ReadOnly 1},

  // Mesh
  N_LAMBDA = {10, Name "Input/02Mesh/00Points per wavelength"},
  LC = {LAMBDA/N_LAMBDA, Name "Input/02Mesh/01Mesh density", ReadOnly 1},
  HEX = {1, Name "Input/02Mesh/02Hexahedra", Choices {0, 1}},

  // Number of subdomains
  N_DOM = {4, Name "Input/03DDM/00Number of subdomains", ReadOnly 1}
];

// Base msh filename
MSH_BASE_NAME = "mesh_";

// Directory for output files
DIR = "out/";

// prefix for (split) mesh files (one for each partition)
MSH_NAME = StrCat(DIR, MSH_BASE_NAME);
