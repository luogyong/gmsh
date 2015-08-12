// Step File
Merge "tesla.stp";

// Mesh
C = 299792458 * 1000; // ! [mm/s] !
F = 1.3e9;
L = C / F;

DefineConstant[ N = {10, Name "Input/00Mesh/00Density [per wavelength (mm)]"} ];
DefineConstant[ O = { 2, Name "Input/00Mesh/01Order [-]"} ];

Characteristic Length {1, 2, 3, 4, 5, 6, 7, 8} = L / N;
Mesh.Optimize            = 1;
Mesh.ElementOrder        = O;
If(O > 1)
  Mesh.HighOrderOptimize = 1;
EndIf

// Volume
Physical Volume(7)  = {1};

// Borders Metallic
Physical Surface(5) = {3, 4, 5, 6};

// Borders In-Out
Physical Surface(4) = {1, 2};

// Print point
Physical Point(100) = {1};
