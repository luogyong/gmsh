// Center //
////////////
// Grab all volumes
For i In {0:RodN - 1}
  allVol[i * 2 + 0] = Vol~{i + 1}[0];
  allVol[i * 2 + 1] = Vol~{i + 1}[1];
EndFor

// Translate
Translate{-RodP * (RodN / 2 - 0.5), 0, 0}{
  Volume{
    allVol[],
    End~{0}[],
    Src~{0},
    Pml~{0}[],
    Pml~{1}[],
    Pml~{2}[],
    Pml~{3}[],
    Pml~{4}[],
    Pml~{5}[],
    Pml~{6}[]
  };
}

// Clear
allVol = {};
