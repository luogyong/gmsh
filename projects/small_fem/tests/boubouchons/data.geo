// User Data //
///////////////
lambda = 0.1;

MeshAir = lambda / 5;
MeshRod = lambda / 10;
MeshSrc = lambda / 20;

RodR = 0.005;
RodL = 0.03;
RodN = 3;
RodP = 0.015;

PmlSize = 0.1 * lambda;
PmlDist = 0.1 * lambda;

AirX = RodP;
AirY = RodP + 0.1 * lambda;
AirZ = RodL + 0.2 * lambda;

SrcR = 0.0005;
SrcX = -RodP + RodR;
SrcL = 0.005;
IsSrcParallel = 0;
