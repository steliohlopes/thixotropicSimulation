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
MeshFactor = 3.2e-5;


Point(1) ={L/2, R, 0,MeshFactor};
Point(2) = {L/2, r0, 0,MeshFactor};

// Criando pontos da garganta
For x In {-l/2:l/2:dx}
    p = newp;
    f =r + (R-r) *Sin(Pi*Fabs(x)/l); 
    Point(p) = {x*(-1),f*(-1) , 0,MeshFactor};
EndFor

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
Field[1].CurvesList = {4:53,57:106};
// Field[1].CurvesList = {131:327:4, 407:603:4, 683:879:4,
//     959:1155:4,117,393,669,945};
Field[1].Sampling = 20;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor / 3;
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

Physical Curve("Inlet") = {55};

Physical Curve("Outlet") = {2};

Physical Curve("Wall") = {56, 54, 53, 57, 52, 58, 51, 59, 50, 60, 49, 61, 48, 62, 47, 63, 46, 64, 45, 65, 44, 66, 43, 67, 42, 68, 41, 69, 40, 70, 39, 71, 38, 72, 37, 73, 36, 74, 35, 75, 34, 76, 33, 77, 32, 78, 31, 79, 30, 80, 29, 81, 28, 82, 27, 83, 26, 84, 25, 85, 24, 86, 23, 87, 22, 88, 21, 89, 20, 90, 19, 91, 18, 92, 17, 93, 16, 94, 15, 95, 14, 96, 13, 97, 12, 98, 11, 99, 10, 100, 9, 101, 8, 102, 7, 103, 6, 104, 5, 105, 4, 106, 3, 1};

Physical Surface("Fluid") = {1};

