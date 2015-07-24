// Step File
Merge "tesla_1cell_UQparam.stp";

// Mesh
C = 299792458 * 1000; // ! [mm/s] !
F = 1.3e9;
L = C / F;
N = 10;

Characteristic Length {1, 2, 3, 4, 5, 6, 7, 8} = L / N;
Mesh.Optimize          = 1;
Mesh.HighOrderOptimize = 1;

// Volume
Physical Volume(7)  = {1};

// Borders Metallic
Physical Surface(5) = {3, 4, 5, 6};

// Borders In-Out
Physical Surface(4) = {1, 2};

// Print point
Physical Point(100) = {1};
