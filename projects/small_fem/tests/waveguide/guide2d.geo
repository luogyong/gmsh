// User constant //
DIM = 2;
Include "guide.dat";

// Geometry (1D) //
Point(1) = {0,  0, 0, LC};
Point(2) = {LX, 0, 0, LC};

Line(1) = {1, 2};
Transfinite Line {1} = (Ceil(LX / LAMBDA * N_LAMBDA)) Using Progression 1;

// Geometry (2D) //
// Extrusion
Extrude {0, LY, 0} {
  Line{1};
  Layers{Ceil(LY / LAMBDA * N_LAMBDA)};
}

// Physicals //
Physical    Line(4) = {4};    // Infinity
Physical    Line(5) = {3};    // Source
Physical    Line(6) = {1, 2}; // Walls
Physical Surface(7) = {5};    // Volume
