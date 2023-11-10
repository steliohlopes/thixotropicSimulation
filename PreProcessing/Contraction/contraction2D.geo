Mesh.MshFileVersion = 2; // Version of the MSH file format to use

// Parametros de geometria
L = 12e-3;
l = 2e-3;
R = 100e-6;
r = 50e-6;
r0 = -R;

// refinamento da garganta
Div=50;
dx=l/Div;

// Parametro de malha
MeshFactor = 2e-5;


Point(1) ={L/2, R, 0,MeshFactor};
Point(2) = {L/2, r0, 0,MeshFactor};

// Criando pontos da garganta
For x In {-l/2:0:dx}
    p = newp;
    f =r + (R-r) *Sin(Pi*Fabs(x)/l); 
    Point(p) = {x*(-1),f*(-1) , 0,MeshFactor};
EndFor

p = newp;
Point(p) = {-L/2,-r , 0,MeshFactor};
p = newp;
Point(p) = {-L/2,r , 0,MeshFactor};

// Criando pontos da garganta
For x In {0:l/2:dx}
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
Field[1].CurvesList = {3:28,31:56};
Field[1].Sampling = 20;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor / 2;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = r/5;
Field[2].DistMax = r*50;

// Field[3] = Distance;
// Field[3].CurvesList = {4};
// Field[3].Sampling = 100;

// Field[4] = Threshold;
// Field[4].InField = 3;
// Field[4].SizeMin = MeshFactor / 3;
// Field[4].SizeMax = MeshFactor;
// Field[4].DistMin = r/10;
// Field[4].DistMax = r*50;


Field[5] = Min;
Field[5].FieldsList = {2};
Background Field = 5;

Physical Curve("Inlet") = {30};

Physical Curve("Outlet") = {2};

Physical Curve("Wall") = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 1, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 29, 31};

Physical Surface("Fluid") = {1};

