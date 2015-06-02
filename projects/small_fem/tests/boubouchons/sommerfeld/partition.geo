// Parititions //
/////////////////

// Iterate //
For i In {0:(DomN - 1)}
  TmpV = {};
  TmpS = {};
  Off  = {};
  InsertSource = 0;

  For j In {0:(DdmN - 1)}
    Off   = j + Range~{i}[0];
    TmpV += Cell~{Off}[];
    TmpS += Infinity~{Off}[];

    If(Off == (EndN - 1))
      InsertSource = 1;
    EndIf
  EndFor

  SetPartition(i + 1){ Volume{ TmpV[] }; Surface{ TmpS[] }; }

  If(InsertSource == 1)
    SetPartition(i + 1){ Surface{ Src~{1}[] }; }
  EndIf
EndFor

// Clear //
TmpV         = {};
TmpS         = {};
Off          = {};
InsertSource = {};