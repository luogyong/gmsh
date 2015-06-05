// Spherical source //
//////////////////////

// 2D half circle //
Pt~{0} = newp;    Point(Pt~{0})  = {SrcX,        0,     0, MeshSrc};
Pt~{1} = newp;    Point(Pt~{1})  = {SrcX + SrcR, 0,     0, MeshSrc};
Pt~{2} = newp;    Point(Pt~{2})  = {SrcX,        +SrcR, 0, MeshSrc};
Pt~{3} = newp;    Point(Pt~{3})  = {SrcX - SrcR, 0,     0, MeshSrc};

Ln~{1} = newl;    Circle(Ln~{1}) = {Pt~{1}, Pt~{0}, Pt~{2}};
Ln~{2} = newl;    Circle(Ln~{2}) = {Pt~{2}, Pt~{0}, Pt~{3}};
Ln~{3} = newl;      Line(Ln~{3}) = {Pt~{3},            Pt~{1}};

Ll = newll;    Line     Loop(Ll) = {Ln~{1}, Ln~{2}, Ln~{3}};
Pl = news;     Plane Surface(Pl) = {Ll};

// Extrude 2 quarters of sphere //
Ext~{1} = Extrude {{1, 0, 0}, {0, 0, 0}, -Pi / 2} { Surface{Pl}; };
Ext~{2} = Extrude {{1, 0, 0}, {0, 0, 0}, +Pi / 2} { Surface{Pl}; };

// Symmetry on surfaces for 2 missing quarters //
Bnd~{1} = Boundary{ Volume{Ext~{1}[1]}; };
Bnd~{2} = Boundary{ Volume{Ext~{2}[1]}; };
Bnd~{3} = Symmetry{0, 1, 0, 0}{ Duplicata{ Surface{Bnd~{1}[]}; } };
Bnd~{4} = Symmetry{0, 1, 0, 0}{ Duplicata{ Surface{Bnd~{2}[]}; } };

// Missing volumes //
Sl~{3} = newsl; Surface Loop(Sl~{3}) = {Ext~{1}[0],
                                        Bnd~{3}[{0, 2, 3}]};
Sl~{4} = newsl; Surface Loop(Sl~{4}) = {Ext~{2}[0],
                                        Bnd~{3}[0],
                                        Bnd~{4}[{2, 3}]};

Vl~{3} = newv;        Volume(Vl~{3}) = {Sl~{3}};
Vl~{4} = newv;        Volume(Vl~{4}) = {Sl~{4}};

// Sphere & boundary //
Src~{0} = {Ext~{1}[1], Ext~{2}[1], Vl~{3}, Vl~{4}};
Src~{1} = CombinedBoundary{ Volume{Src~{0}[]}; };

// Insert in End~{0}[0] //
//////////////////////////
BndEnd = Boundary{ Volume{End~{0}[0]}; };
Delete{ Volume{End~{0}[0]}; }

EndSl = newsl; Surface Loop(EndSl) = {BndEnd[]};
SphSl = newsl; Surface Loop(SphSl) = {Src~{1}[]};
Volume(End~{0}[0]) = {EndSl, SphSl};

// Clear //
Pt~{0}  = {};
Pt~{1}  = {};
Pt~{2}  = {};
Pt~{3}  = {};
Ln~{1}  = {};
Ln~{2}  = {};
Ln~{3}  = {};
Ll      = {};
Pl      = {};
Ext~{1} = {};
Ext~{2} = {};
Bnd~{1} = {};
Bnd~{2} = {};
Bnd~{3} = {};
Bnd~{4} = {};
Sl~{3}  = {};
Sl~{4}  = {};
Vl~{3}  = {};
Vl~{4}  = {};
BndEnd  = {};
EndSl   = {};
EndSl   = {};
