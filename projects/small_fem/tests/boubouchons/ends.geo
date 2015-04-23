// Ends //
//////////
// End~{i} = {0:EndN - 1}}[0:1] = [End(i)-, End(i)+]

// Grab Surfaces //
VolSurf~{1}    = Boundary{ Volume{Vol~{1}[1]};    };
VolSurf~{RodN} = Boundary{ Volume{Vol~{RodN}[1]}; };

// Number of extrusion //
EndN = Ceil(AirX / RodP);

// There is at least on extrusion //
If(EndN == 0)
  EndN = 1;
EndIf

// Extrude once //
Tmp~{0} = Extrude{-RodP, 0, 0}{ Surface{VolSurf~{1}[5]};    };
Tmp~{1} = Extrude{+RodP, 0, 0}{ Surface{VolSurf~{RodN}[3]}; };

End~{0}[0] = Tmp~{0}[1];
End~{0}[1] = Tmp~{1}[1];

// Continue extruding //
For i In {1:EndN - 1}
  Tmp~{0} = Boundary{ Volume{End~{i - 1}[0]}; };
  Tmp~{1} = Boundary{ Volume{End~{i - 1}[1]}; };

  Tmp~{0} = Extrude{-RodP, 0, 0}{ Surface{Tmp~{0}[1]}; };
  Tmp~{1} = Extrude{+RodP, 0, 0}{ Surface{Tmp~{1}[1]}; };

  End~{i}[0] = Tmp~{0}[1];
  End~{i}[1] = Tmp~{1}[1];
EndFor

// Clear useless variables //
VolSurf~{1}    = {};
VolSurf~{RodN} = {};
Tmp~{0}        = {};
Tmp~{1}        = {};
