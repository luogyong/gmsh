// Center //
////////////
// Grab all volumes
For i In {0:RodN - 1}
  allVol[i * 2 + 0] = Vol~{i + 1}[0];
  allVol[i * 2 + 1] = Vol~{i + 1}[1];
EndFor

// Grab all ends
For i In {0:EndN - 1}
  allEnd[i * 2 + 0] = End~{i}[0];
  allEnd[i * 2 + 1] = End~{i}[1];
EndFor

// Translate
Translate{-RodP * (RodN / 2 - 0.5), 0, 0}{
  Volume{allVol[], allEnd[], Src~{0}[]};
}

// Clear
allVol = {};
allEnd = {};
