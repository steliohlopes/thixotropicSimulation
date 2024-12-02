// Gmsh project created on Wed May 29 13:27:32 2024
Mesh.MshFileVersion = 2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// Parametros de geometria
R_inlet = 10e-3;
L_inlet = 80e-3; 
H_outlet = 1.5e-3;
W_outlet = 168e-3; //half slot
alfa = 30*(Pi/180); //radianos
L_outlet = 50e-3;
R1 = 10e-3;
R2 = 5e-3;
slotCenter = W_outlet*Tan(alfa);
// Parametro de malha
MeshFactor = 4e-3;

// R1
pR1Back = newp; Point(pR1Back) = {0,0,0,MeshFactor};
pR1Center = newp; Point(pR1Center) ={R1,0,0,MeshFactor}; 
pR1Up = newp; Point(pR1Up) ={R1,R1,0,MeshFactor}; 
pR1Forth = newp; Point(pR1Forth) = {Sqrt((R1)^2-(H_outlet/2)^2)+R1,H_outlet/2,0,MeshFactor};

cR1Back = newc; Circle(cR1Back) = {pR1Back,pR1Center,pR1Up};
cR1Forth = newc; Circle(cR1Forth) = {pR1Up,pR1Center,pR1Forth};

// R2
pR2Center = newp; Point(pR2Center) = {2*R1+slotCenter-R2,0,W_outlet,MeshFactor}; 
pR2Back = newp; Point(pR2Back) = {2*R1+slotCenter-2*R2,0,W_outlet,MeshFactor};
pR2Up = newp; Point(pR2Up) = {2*R1+slotCenter-R2,R2,W_outlet,MeshFactor}; 
pR2Forth = newp; Point(pR2Forth) = {2*R1+slotCenter-R2 +Sqrt((R2)^2-(H_outlet/2)^2),H_outlet/2,W_outlet,MeshFactor}; 

cR2Back = newc; Circle(cR2Back) = {pR2Back,pR2Center,pR2Up};
cR2Forth = newc; Circle(cR2Forth) = {pR2Up,pR2Center,pR2Forth};
//R1 & R2

lR1R2Back = newl; Line(lR1R2Back) = {pR1Back,pR2Back};

lR1R2Forth = newl; Line(lR1R2Forth) = {pR1Forth,pR2Forth};
lR1R2Up = newl; Line(lR1R2Up) = {pR1Up,pR2Up};

llR1R2Back = newll; Line Loop(llR1R2Back) = {cR1Back,lR1R2Up,cR2Back,lR1R2Back};
sR1R2Back = news; Surface(sR1R2Back) = {llR1R2Back};

llR1R2Forth = newll; Line Loop(llR1R2Forth) = {cR1Forth,lR1R2Up,cR2Forth,lR1R2Forth};
sR1R2Forth = news; Surface(sR1R2Forth) = {llR1R2Forth};

//Outlet

pOutletCenterUp = newp; Point(pOutletCenterUp) = {2*R1+slotCenter+L_outlet,H_outlet/2,0,MeshFactor};
pOutletWallUp = newp; Point(pOutletWallUp) = {2*R1+slotCenter+L_outlet,H_outlet/2,W_outlet,MeshFactor};

pOutletCenterDown = newp; Point(pOutletCenterDown) = {2*R1+slotCenter+L_outlet,0,0,MeshFactor};
pOutletWallDown = newp; Point(pOutletWallDown) = {2*R1+slotCenter+L_outlet,0,W_outlet,MeshFactor};


lOutletUp = newl; Line(lOutletUp) = {pOutletCenterUp,pOutletWallUp};
lOutletDown = newl; Line(lOutletDown) = {pOutletCenterDown,pOutletWallDown};
lOutletWall = newl; Line(lOutletWall) = {pOutletWallUp,pOutletWallDown};
lOutletCenter = newl; Line(lOutletCenter) = {pOutletCenterUp,pOutletCenterDown};
llOutlet = newll; Line Loop(llOutlet) = {lOutletUp,lOutletWall,lOutletDown,lOutletCenter};
sOutlet = news; Surface(sOutlet) = {llOutlet};

//Wall

lWallDown = newl;  Line(lWallDown) = {pR2Center,pOutletWallDown};
lWallDown2 = newl;  Line(lWallDown2) = {pR2Center,pR2Back};
lWallUp = newl;  Line(lWallUp) = {pR2Forth,pOutletWallUp};

llWall = newll; Line Loop(llWall) = {lWallDown,lWallDown2,cR2Back,cR2Forth,lWallUp,lOutletWall};
sWall = news; Surface(sWall) = {llWall};

//Center

lCenterDown = newl;  Line(lCenterDown) = {pR1Center,pOutletCenterDown};
lCenterDown2 = newl;  Line(lCenterDown2) = {pR1Center,pR1Back};
lCenterUp = newl;  Line(lCenterUp) = {pR1Forth,pOutletCenterUp};

llCenter = newll; Curve Loop(llCenter) = {lCenterDown,lCenterDown2,cR1Back,cR1Forth,lCenterUp,lOutletCenter};
sCenter = news; Surface(sCenter) = {llCenter};

//Up Slot

llUp = newll; Line Loop(llUp) = {lR1R2Forth,lWallUp,lOutletUp,lCenterUp};
sUp = news; Plane Surface(sUp) = {llUp};

//Down Slot

llDown = newll; Curve Loop(llDown) = {lR1R2Back,lWallDown2,lWallDown,lOutletDown,lCenterDown,lCenterDown2};
sDown = news;Plane Surface(sDown) = {llDown};

slSlot = newsl; Surface Loop(slSlot) = {sR1R2Forth,sR1R2Back,sDown,sOutlet,sUp,sCenter,sWall};
vSlot = newv; Volume(vSlot) = {slSlot};

// INLET CYLINDER

pInletRight = newp; Point(pInletRight) = {R1,0,R1,MeshFactor};

cInlet = newc; Circle(cInlet) = {pR1Up,pR1Center,pInletRight};

lInletLeft = newl; Line(lInletLeft) = {pR1Up,pR1Center};
lInletDown = newl; Line(lInletDown) = {pR1Center,pInletRight};

llInlet = newll; Line Loop(llInlet) = {cInlet,lInletLeft,lInletDown};
sInlet = news; Plane Surface(sInlet) = {llInlet};

Extrude {-L_inlet,0,0} {
    Surface{sInlet};
}

v() = BooleanDifference{ Volume{34};Delete;}{ Volume{vSlot};};
Coherence;

Physical Surface("Inlet") = {36};
Physical Surface("Outlet") = {41};
Physical Surface("Wall") = {32,33,39,38,42,44};
Physical Surface("SymmetryZ") = {37,43};
Physical Surface("SymmetryY") = {35,40};
Physical Volume("Fluid") = {33,34};


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
Field[1].CurvesList = {44,49};//outlet lines

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor/3;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = W_outlet;
Field[2].DistMax =W_outlet+R2;

Field[3] = Distance;
Field[3].CurvesList = {29,33,35};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = MeshFactor;
Field[4].SizeMax = MeshFactor;
Field[4].DistMin = L_inlet;
Field[4].DistMax =L_inlet+W_outlet;

Field[5] = Min;
Field[5].FieldsList = {2,4};
Background Field = 5;//+