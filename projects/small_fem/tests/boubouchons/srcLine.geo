// 2D Source //
///////////////
SrcPt~{0} = newp;  Point(SrcPt~{0}) = {SrcX,        0,     -SrcL / 2, MeshSrc};
SrcPt~{1} = newp;  Point(SrcPt~{1}) = {SrcX + SrcR, 0,     -SrcL / 2, MeshSrc};
SrcPt~{2} = newp;  Point(SrcPt~{2}) = {SrcX,        +SrcR, -SrcL / 2, MeshSrc};
SrcPt~{3} = newp;  Point(SrcPt~{3}) = {SrcX - SrcR, 0,     -SrcL / 2, MeshSrc};
SrcPt~{4} = newp;  Point(SrcPt~{4}) = {SrcX,        -SrcR, -SrcL / 2, MeshSrc};

SrcCr~{1} = newl; Circle(SrcCr~{1}) = {SrcPt~{1}, SrcPt~{0}, SrcPt~{2}};
SrcCr~{2} = newl; Circle(SrcCr~{2}) = {SrcPt~{2}, SrcPt~{0}, SrcPt~{3}};
SrcCr~{3} = newl; Circle(SrcCr~{3}) = {SrcPt~{3}, SrcPt~{0}, SrcPt~{4}};
SrcCr~{4} = newl; Circle(SrcCr~{4}) = {SrcPt~{4}, SrcPt~{0}, SrcPt~{1}};

SrcLl = newll; Line     Loop(SrcLl) = {SrcCr~{1}, SrcCr~{2},
                                       SrcCr~{3}, SrcCr~{4}};
SrcSf = news;  Plane Surface(SrcSf) = {SrcLl};

// 3D Source //
///////////////
Tmp     = Extrude{0, 0, SrcL}{ Surface{SrcSf}; };
Src~{0} = Tmp[1];

// Get Source Boundary //
/////////////////////////
Src~{1} = Boundary{ Volume{Src~{0}}; };

// Rotate //
////////////
If(IsSrcParallel == 0)
  Rotate{{1, 0, 0}, {0, 0, 0}, Pi / 2}{ Volume{Src~{0}}; }
EndIf

// Insert in End~{0}[0] //
//////////////////////////
BndEnd = Boundary{ Volume{End~{0}[0]}; };
BndSrc = Boundary{ Volume{Src~{0}};    };

Delete{ Volume{End~{0}[0]}; }

EndSl = newsl; Surface Loop(EndSl) = {BndEnd[]};
SrcSl = newsl; Surface Loop(SrcSl) = {BndSrc[]};

Volume(End~{0}[0]) = {EndSl, SrcSl};

// Remove Volume //
///////////////////
Delete{ Volume{Src~{0}}; }

// Clear //
///////////
SrcPt~{0} = {};
SrcPt~{1} = {};
SrcPt~{2} = {};
SrcPt~{3} = {};
SrcPt~{4} = {};
SrcCr~{1} = {};
SrcCr~{2} = {};
SrcCr~{3} = {};
SrcCr~{4} = {};
SrcLl     = {};
SrcSf     = {};
Tmp       = {};
BndEnd    = {};
BndSrc    = {};
EndSl     = {};
SrcSl     = {};
