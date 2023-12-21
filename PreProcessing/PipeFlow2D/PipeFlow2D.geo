Mesh.MshFileVersion = 2; // Version of the MSH file format to use

// Parametros de geometria
L = 24e-3;
l = 4e-3;
R = 100e-6;
r = 50e-6;
r0 = -R;

// refinamento da garganta
Div=100;
dx=l/Div;

// Parametro de malha
MeshFactor = 2e-5;


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
Field[1].CurvesList = {4:102,106:204};
Field[1].Sampling = 20;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = MeshFactor / 2;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = r*3;
Field[2].DistMax = r*40;

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

Physical Curve("Inlet") = {104};

Physical Curve("Outlet") = {2};

Physical Curve("Wall") = {66, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 79, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 91, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 53, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 15, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 28, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 40, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 105, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 167, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 168, 180, 166, 165, 164, 163, 162, 161, 160, 159, 158, 157, 156, 192, 204, 203, 202, 201, 200, 199, 198, 197, 196, 195, 194, 193, 155, 191, 190, 189, 188, 187, 186, 185, 184, 183, 182, 181, 117, 129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 130, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 142, 154, 153, 152, 151, 150, 149, 148, 147, 146, 145, 144, 143, 3, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 1};

Physical Surface("Fluid") = {1};
