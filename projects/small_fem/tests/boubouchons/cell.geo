// Gmsh Data //
///////////////
Geometry.ExtrudeReturnLateralEntities = 0; // Easier to count

// Zero //
//////////
Point(0) = {0, 0, 0};

// 2D Cell //
/////////////
Point(1) = {-AirX / 2, -AirY / 2, -AirZ / 2, MeshAir};
Point(2) = {+AirX / 2, -AirY / 2, -AirZ / 2, MeshAir};
Point(3) = {+AirX / 2, +AirY / 2, -AirZ / 2, MeshAir};
Point(4) = {-AirX / 2, +AirY / 2, -AirZ / 2, MeshAir};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};

Line Loop(1)     = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 3D Shell //
//////////////
Extrude{0, 0, AirZ}{ Surface{1}; }

BndAir = Boundary{ Volume{1}; };
Delete{ Volume{1}; }

// Rod //
/////////
RodPt~{0} = newp;  Point(RodPt~{0}) = {0,     0,     -RodL / 2, MeshRod};
RodPt~{1} = newp;  Point(RodPt~{1}) = {+RodR, 0,     -RodL / 2, MeshRod};
RodPt~{2} = newp;  Point(RodPt~{2}) = {0,     +RodR, -RodL / 2, MeshRod};
RodPt~{3} = newp;  Point(RodPt~{3}) = {-RodR, 0,     -RodL / 2, MeshRod};
RodPt~{4} = newp;  Point(RodPt~{4}) = {0,     -RodR, -RodL / 2, MeshRod};

RodCr~{1} = newl; Circle(RodCr~{1}) = {RodPt~{1}, RodPt~{0}, RodPt~{2}};
RodCr~{2} = newl; Circle(RodCr~{2}) = {RodPt~{2}, RodPt~{0}, RodPt~{3}};
RodCr~{3} = newl; Circle(RodCr~{3}) = {RodPt~{3}, RodPt~{0}, RodPt~{4}};
RodCr~{4} = newl; Circle(RodCr~{4}) = {RodPt~{4}, RodPt~{0}, RodPt~{1}};

RodLl = newll; Line     Loop(RodLl) = {RodCr~{1}, RodCr~{2},
                                       RodCr~{3}, RodCr~{4}};
RodSf = news;  Plane Surface(RodSf) = {RodLl};

Extrude{0, 0, RodL}{ Surface{RodSf}; }

// 3D Cell //
/////////////
BndRod = Boundary{ Volume{1}; };

AirSl  = newsl; Surface Loop(AirSl) = {BndAir[0], BndAir[1], BndAir[2],
                                       BndAir[3], BndAir[4], BndAir[5]};
RodSl  = newsl; Surface Loop(RodSl) = {BndRod[0], BndRod[1], BndRod[2],
                                       BndRod[3], BndRod[4], BndRod[5]};
Volume(2) = {AirSl, RodSl};

// Save Volumes //
//////////////////
Vol~{1} = {1, 2}; // Vol~{1}[0, 1] = [Rod(1), Air(1)]

// Clear useless variables //
/////////////////////////////
BndAir    = {};
RodPt~{0} = {};
RodPt~{1} = {};
RodPt~{2} = {};
RodPt~{3} = {};
RodPt~{4} = {};
RodCr~{1} = {};
RodCr~{2} = {};
RodCr~{3} = {};
RodCr~{4} = {};
RodLl     = {};
RodSf     = {};
BndRod    = {};
AirSl     = {};
RodSl     = {};
