Mesh.MshFileVersion = 2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// Parametros de geometria
D_inlet = 0.2e-3;
L_inlet = 1e-3; 
H_inlet = 0;
H_outlet = 0.2e-3;
W_outlet = 20e-3;
L_outlet = 1e-3;
R1 = 2e-3;
R2 = 1e-3;

// Parametro de malha
MeshFactor = 1.7e-4;

// Inlet cylinder 
// Cylinder 

pInletCenter = newp; Point(pInletCenter) = {0,0,0,MeshFactor}; //center
pInletRight = newp; Point(pInletRight) = {0,0,D_inlet/2,MeshFactor};
pInletUp = newp; Point(pInletUp) = {0,D_inlet/2,0,MeshFactor};
pInletDown = newp; Point(pInletDown) = {0,-D_inlet/2,0,MeshFactor};

cInletDown = newl; Circle(cInletDown) = {pInletDown,pInletCenter,pInletRight};
cInletUp = newl; Circle(cInletUp) = {pInletRight,pInletCenter,pInletUp};
lInletUp = newl; Line(lInletUp) = {pInletUp,pInletCenter};
lInletDown = newl; Line(lInletDown) = {pInletCenter,pInletDown};

llInlet = newll; Line Loop(llInlet) = {cInletDown,cInletUp,lInletUp,lInletDown};

sInlet = news; Plane Surface(sInlet) = {llInlet};

Extrude {-L_inlet,0,0} {
  Surface{sInlet};
}


// R1 Wall
pR1Center = newp; Point(pR1Center) ={R1,-D_inlet/2,0,MeshFactor}; //center
pR1Down = newp; Point(pR1Down) ={R1,-D_inlet/2-R1,0,MeshFactor}; //down
pR1Right = newp; Point(pR1Right) ={2*R1,-D_inlet/2,0,MeshFactor}; //right

cR1Left = newl; Circle(cR1Left) = {pInletDown,pR1Center,pR1Down};
cR1Right = newl; Circle(cR1Right) = {pR1Down,pR1Center,pR1Right};
lR1Left = newl; Line(lR1Left) = {pInletDown,pR1Center};
lR1Right = newl; Line(lR1Right) = {pR1Center,pR1Right};

llR1 = newll; Line Loop(llR1) = {lR1Left,lR1Right,-cR1Right,-cR1Left};
sR1 = news; Plane Surface(sR1) = {llR1};

//R2 Wall
pR2Right = newp; Point(pR2Right) = {2*R1,-D_inlet/2,W_outlet/2,MeshFactor}; //right
pR2Center = newp; Point(pR2Center) = {2*R1-R2,-D_inlet/2,W_outlet/2,MeshFactor}; //center
pR2Down = newp; Point(pR2Down) = {2*R1-R2,-D_inlet/2-R2,W_outlet/2,MeshFactor}; //down
pR2Left = newp; Point(pR2Left) = {2*R1-2*R2,-D_inlet/2,W_outlet/2,MeshFactor}; //left

cR2Left = newl; Circle(cR2Left) = {pR2Left,pR2Center,pR2Down};
cR2Right = newl; Circle(cR2Right) = {pR2Down,pR2Center,pR2Right};
lR2Left = newl; Line(lR2Left) = {pR2Left,pR2Center};
lR2Right = newl; Line(lR2Right) = {pR2Center,pR2Right};

llR2 = newll; Line Loop(llR2) = {lR2Left,lR2Right,-cR2Right,-cR2Left};
sR2 = news; Plane Surface(sR2) = {llR2};

//Combine R1 & R2
lR12Left = newl; Line(lR12Left) = {pInletDown,pR2Left}; // Left
lR12Down = newl; Line(lR12Down) = {pR1Down,pR2Down};// Down
lR12Right = newl; Line(lR12Right) = {pR1Right,pR2Right};// Right

// llR12Up = newll; Line Loop(llR12Up) = {lR1Left,lR1Right,lR12Right,lR2Right,lR2Left,lR12Left};
// sR12Up = news; Plane Surface(sR12Up) = {llR12Up};

llR12Right = newll; Line Loop(llR12Right) = {cR1Right,lR12Down,cR2Right,lR12Right};
sR12Right = news; Surface(sR12Right) = {llR12Right};

llR12Left = newll; Line Loop(llR12Left) = {cR1Left,lR12Down,cR2Left,lR12Left};
sR12Left = news; Surface(sR12Left) = {llR12Left};

// slR12 =newsl; Surface Loop(slR12) = {sR1,sR2,sR12Left,sR12Right};
// vR12 = newv; Volume(vR12) = {slR12};

// Outlet 
pOutletCenterUp = newp; Point(pOutletCenterUp) ={2*R1+L_outlet,D_inlet/2,0,MeshFactor};
pOutletWallUp = newp; Point(pOutletWallUp) = {2*R1+L_outlet,+D_inlet/2,W_outlet/2,MeshFactor};
pOutletCenterDown = newp; Point(pOutletCenterDown) ={2*R1+L_outlet,D_inlet/2-H_outlet,0,MeshFactor};
pOutletWallDown = newp; Point(pOutletWallDown) = {2*R1+L_outlet,+D_inlet/2-H_outlet,W_outlet/2,MeshFactor};

lOutletDown = newl; Line(lOutletDown) ={pOutletCenterDown,pOutletWallDown};
lOutletWall = newl; Line(lOutletWall) ={pOutletWallDown,pOutletWallUp};
lOutletUp = newl; Line(lOutletUp) ={pOutletWallUp,pOutletCenterUp};
lOutletCenter = newl; Line(lOutletCenter) ={pOutletCenterUp,pOutletCenterDown};

llOutlet = newll; Line Loop(llOutlet)={lOutletDown,lOutletWall,lOutletUp,lOutletCenter};
sOutlet = news; Surface(sOutlet) = {llOutlet};

//Slot Center
lSlotCenterUp = newl; Line(lSlotCenterUp) = {pInletUp,pOutletCenterUp};
lSlotCenterDown = newl; Line(lSlotCenterDown) = {pR1Right,pOutletCenterDown};

llSlotCenter= newll; Line Loop(llSlotCenter)={lInletDown,lInletUp,lSlotCenterUp,lOutletCenter,lSlotCenterDown,lR1Right,lR1Left};
sSlotCenter = news; Surface(sSlotCenter) = {llSlotCenter};

//Slot Back
pR2Slot = newp; Point(pR2Slot) = {2*R1-2*R2,-D_inlet/2+H_outlet,W_outlet/2,MeshFactor};

lSlotBackUp = newl; Line(lSlotBackUp) = {pInletUp,pR2Slot};
lSlotBackWall = newl; Line(lSlotBackWall) = {pR2Left,pR2Slot};

llSlotBack = newll; Line Loop(llSlotBack) = {cInletDown,cInletUp,lSlotBackUp,lSlotBackWall,lR12Left};
sSlotBack = news; Surface(sSlotBack) = {llSlotBack};

//Slot Wall
lSlotWallUp = newl; Line(lSlotWallUp) = {pR2Slot,pOutletWallUp};
lSlotWallDown = newl; Line(lSlotWallDown) = {pOutletWallDown,pR2Right};

llSlotWall = newll; Line Loop(llSlotWall) = {lSlotBackWall,lSlotWallUp,lOutletWall,lSlotWallDown,lR2Right,lR2Left};
sSlotWall = news; Surface(sSlotWall) = {llSlotWall};
//Slot Up
llSlotUp = newll; Line Loop(llSlotUp) = {lSlotCenterUp,lOutletUp,lSlotWallUp,lSlotBackUp};
sSlotUp = news; Surface(sSlotUp) = {llSlotUp};

//Slot down
llSlotDown = newll; Line Loop(llSlotDown) = {lR12Right,lSlotCenterDown,lOutletDown,lSlotWallDown};
sSlotDown = news; Surface(sSlotDown) = {llSlotDown};

//Create Volume
slCoating = newsl; Surface Loop(slCoating)={sR1,sR2,sR12Right,sR12Left,sOutlet,sSlotCenter,sSlotBack,sSlotWall,sSlotUp,sSlotDown,sInlet};
vCoating = newv; Volume(vCoating)={slCoating};

v() = BooleanFragments{ Volume{1}; Delete;}{ Volume{vCoating}; Delete; };


Physical Surface("Inlet") = {11};
Physical Surface("Outlet") = {sOutlet};
Physical Surface("Wall") = {8,7,sR2,sSlotWall,sSlotBack,sSlotUp,sSlotDown,sR12Left,sR12Right};
Physical Surface("Symmetry") = {9,10,sR1,sSlotCenter};
Physical Volume("Fluid") = {1,vCoating};

// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances

// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax

Field[1] = Distance;
Field[1].CurvesList = {lOutletUp,lOutletDown};//outlet lines

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor/4.5;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = L_outlet;
Field[2].DistMax =L_outlet+R1;

Field[3] = Distance;
Field[3].CurvesList = {cInletDown,cInletUp,lInletDown,lInletUp};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = MeshFactor / 6;
Field[4].SizeMax = MeshFactor;
Field[4].DistMin = L_inlet;
Field[4].DistMax =L_inlet*5;

Field[5] = Min;
Field[5].FieldsList = {2,4};
Background Field = 5;

