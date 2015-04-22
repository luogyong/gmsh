// Grab Surfaces //
///////////////////
VolSurf~{1}    = Boundary{ Volume{Vol~{1}[1]};    };
VolSurf~{RodN} = Boundary{ Volume{Vol~{RodN}[1]}; };

// Ends //
//////////
// End~{0}[0:1] = [End-, End+]
Tmp~{0} = Extrude{-AirX, 0, 0}{ Surface{VolSurf~{1}[5]};    };
Tmp~{1} = Extrude{+AirX, 0, 0}{ Surface{VolSurf~{RodN}[3]}; };

End~{0}[0] = Tmp~{0}[1];
End~{0}[1] = Tmp~{1}[1];

// Clear useless variables //
/////////////////////////////
VolSurf~{1}    = {};
VolSurf~{RodN} = {};
Tmp~{0}        = {};
Tmp~{1}        = {};
