// Translate Cell //
////////////////////
// Vol~{i}[0, 1] = [Rod(i), Air(i)]
For i In {2:RodN}
  Vol~{i} =
    Translate{RodP, 0, 0}{
      Duplicata{
        Volume{Vol~{i-1}[0], Vol~{i-1}[1]};
      }
    };
EndFor
