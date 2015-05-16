// Splitting geometry in two for DDM //
///////////////////////////////////////

// Rod //
Surf = {};
Surf = Boundary{ Volume{Vol~{RodN / 2}[1]}; };
DdmRod~{0} = Surf[3];

// Pml Z down //
Surf = {};
Surf = Boundary{ Volume{Pml~{0}[RodN / 2]}; };

DdmRod~{1} = Surf[5];

// Pml Z up //
Surf = {};
Surf = Boundary{ Volume{Pml~{0}[RodN / 2 + RodN + 2 * EndN]}; };

DdmRod~{2} = Surf[5];

// Pml Y rear //
Surf = {};
Surf = Boundary{ Volume{Pml~{1}[RodN / 2]}; };

DdmRod~{3} = Surf[5];

// Pml Y front //
Surf = {};
Surf = Boundary{ Volume{Pml~{1}[RodN / 2 + RodN + 2 * EndN]}; };

DdmRod~{4} = Surf[3];

// Pml YZ rear down //
Surf = {};
Surf = Boundary{ Volume{Pml~{2}[RodN / 2 + 0 * (RodN + 2 * EndN)]}; };

DdmRod~{5} = Surf[5];

// Pml YZ rear up //
Surf = {};
Surf = Boundary{ Volume{Pml~{2}[RodN / 2 + 1 * (RodN + 2 * EndN)]}; };

DdmRod~{6} = Surf[5];

// Pml YZ front down //
Surf = {};
Surf = Boundary{ Volume{Pml~{2}[RodN / 2 + 2 * (RodN + 2 * EndN)]}; };

DdmRod~{7} = Surf[3];

// Pml YZ front up //
Surf = {};
Surf = Boundary{ Volume{Pml~{2}[RodN / 2 + 3 * (RodN + 2 * EndN)]}; };

DdmRod~{8} = Surf[3];

// Clear //
Surf = {};
