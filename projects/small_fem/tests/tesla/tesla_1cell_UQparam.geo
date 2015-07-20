// Step File
Merge "tesla_1cell_UQparam.stp";

// Mesh
Characteristic Length {1, 2, 3, 4, 5, 6, 7, 8} = 10;
Mesh.ElementOrder      = 2;
Mesh.Optimize          = 1;
Mesh.HighOrderOptimize = 1;

// Volume
Physical Volume(7)  = {1};

// Borders
Physical Surface(5) = {1, 2, 3, 4, 5, 6};

// Print point
Physical Point(100) = {1};
