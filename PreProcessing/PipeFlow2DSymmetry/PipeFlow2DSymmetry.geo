Mesh.MshFileVersion = 2; // Version of the MSH file format to use

// Parametros de geometria
L = 24e-3;
l = 4e-3;
R = 100e-6;
r = 50e-6;
r0 = 0;

// refinamento da garganta
Div=100;
dx=l/Div;

// Parametro de malha
MeshFactor = 1.7e-5;


Point(1) ={L/2, R, 0,MeshFactor};
Point(2) = {L/2, r0, 0,MeshFactor};


p = newp;
Point(p) = {-L/2,r0 , 0,MeshFactor};
p = newp;
Point(p) = {-L/2,R , 0,MeshFactor};

// Criando pontos da garganta
For x In {-l/2:l/2:dx}
    p = newp;
    f =r + (R-r) *Sin(Pi*Fabs(x)/l); 
    Point(p) = {x,f , 0,MeshFactor};
EndFor

// Linhas
l = newl;
Line(l) = {newp-1, 1};

For point In {1:newp-2}
    l = newl;
    Line(l) = {point, point+1} ; 
EndFor

Line Loop(1) = {1:l};
Plane Surface(1) = (1);

// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax

Field[1] = Distance;
Field[1].CurvesList = {6:104};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor / 2;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = r*3;
Field[2].DistMax = r*40;

// Field[3] = Distance;
// Field[3].CurvesList = {4};

// Field[4] = Threshold;
// Field[4].InField = 3;
// Field[4].SizeMin = MeshFactor / 3;
// Field[4].SizeMax = MeshFactor;
// Field[4].DistMin = r/10;
// Field[4].DistMax = r*50;


Field[5] = Min;
Field[5].FieldsList = {2};
Background Field = 5;

Physical Curve("Inlet") = {4};

Physical Curve("Outlet") = {2};

Physical Curve("Wall") = {1,5:104};

Physical Curve("Symmetry") = {3};

Physical Surface("Fluid") = {1};
