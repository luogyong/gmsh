// User Data //
///////////////
lambda = 0.1;

MeshAir = lambda / 5;
MeshRod = lambda / 10;

RodR = 0.005;
RodL = 0.03;
RodN = 3;
RodP = 0.015;

PmlSize = 0.1 * lambda;
PmlDist = 0.1 * lambda;

AirX = RodP;
AirY = RodP + 0.1 * lambda;
AirZ = RodL + 0.2 * lambda;

// Create Geomtetry //
//////////////////////
Include "cell.geo";
Include "grid.geo";
Include "ends.geo";
Include  "pml.geo";

// Center //
////////////
For i In {0:RodN - 1}
  allVol[i * 2 + 0] = Vol~{i + 1}[0];
  allVol[i * 2 + 1] = Vol~{i + 1}[1];
EndFor

Translate{-RodP * (RodN / 2 - 0.5), 0, 0}{
  Volume{
    allVol[],
    End~{0}[],
    End~{1}[],
    Pml~{0}[],
    Pml~{1}[],
    Pml~{2}[],
    Pml~{3}[],
    Pml~{4}[],
    Pml~{5}[],
    Pml~{6}[]
  };
}

allVol = {};
