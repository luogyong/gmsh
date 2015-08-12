// User constant //
DIM = 3;
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

// Geometry (3D) //
// Extrusion
Extrude {0, 0, LZ} {
  Surface{5};
  Layers{Ceil(LZ / LAMBDA * N_LAMBDA)};
}

// Physicals //
Physical Surface(4) = {18};            // Infinity
Physical Surface(5) = {26};            // Source
Physical Surface(6) = {5, 22, 27, 14}; // Walls
Physical  Volume(7) = {1};             // Volume
