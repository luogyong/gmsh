// Parititions //
/////////////////

// Iterate //
For i In {0:(DomN - 1)}
  Tmp = {};
  Off = {};
  InsertSource = 0;

  For j In {0:(DdmN - 1)}
    Off  = j + Range~{i}[0];
    Tmp += Cell~{Off}[];

    If(Off == (PmlN + EndN - 1))
      InsertSource = 1;
    EndIf
  EndFor

  SetPartition(i + 1){ Volume{ Tmp[] }; }

  If(InsertSource == 1)
    SetPartition(i + 1){ Surface{ Src~{1}[] }; }
  EndIf
EndFor

// Clear //
Tmp          = {};
Off          = {};
InsertSource = {};