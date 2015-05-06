DefineConstant[
  // Analysis type
  ANALYSIS = {1, Name "Input/00Type of analysis",
    Choices {0="Helmholtz", 1="Maxwell"}},
  // frequency
  WAVENUMBER = {5, Min 0.1, Max 31.5, Step 0.1, Name "Input/0Wavenumber"},
  LAMBDA = {2*Pi/WAVENUMBER, Name "Input/1Wavelength", ReadOnly 1},
  // number of points per wavelength
  N_LAMBDA = {10, Name "Input/2Points per wavelength"},
  // geomtrical element order
  ELEMENT_ORDER = {1, Name "Input/Geometrical order"},
  // incident angle
  THETA_INC = {0, Min 0., Max 2*Pi, Step 0.1, Name "Input/3Incident angle"},
  // geometry
  R_INT = {1/3, Name "Input/4Internal radius"},
  R_EXT = {3/3, Name "Input/5External radius"},
  // number of subdomains
  N_DOM = {4, Min 1, Max 20, Step 1, Name "Input/04Number of subdomains"},
  // base msh filename
  MSH_BASE_NAME = "mesh_",
  // directory for output files
  DIR = "out/"
];

LC = LAMBDA/N_LAMBDA;
Printf("LC: %f", LC);

// prefix for (split) mesh files (one for each partition)
MSH_NAME = StrCat(DIR, MSH_BASE_NAME) ;
