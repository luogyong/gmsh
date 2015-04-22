// User Data //
///////////////

// Physical constants //
////////////////////////
C0   = 299792458;            // [m/s]
Mu0  = 4 * Pi * 1e-7;        // [H/m]
Eps0 = 1. / (C0 * C0 * Mu0); // [F/m]


// Rods //
//////////
RodN = 10;    // [-]
RodR = 0.005; // [m]
RodP = 0.015; // [m]
RodL = 0.03;  // [m]


// Source //
////////////
Freq   = 6e9;        // [Hz]
Lambda = C0 / Freq;  // [m]

SrcR = 0.0005;       // [m]
SrcL = 0.005;        // [m]
SrcX = -RodP + RodR; // [m]

IsSrcParallel = 0;   // [bool]


// Air //
/////////
AirX = RodP;            // [m]
AirY = RodR * 9;        // [m]
AirZ = RodR * 6 + RodL; // [m]


// Pml //
/////////
PmlSize = 0.8 * Lambda; // [m]


// Mesh //
//////////
MeshAir = Lambda / 7;
MeshRod = Lambda / 10;
MeshSrc = Lambda / 50;
