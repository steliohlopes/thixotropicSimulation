Mesh.MshFileVersion = 2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// Parametros de geometria
D_inlet=2e-3;
L_inlet=2e-3; //Verificar medida
H_inlet=0.2e-3;
H_outlet=H_inlet;
W_outlet=10e-3;
L_outlet=2e-3;
R=2e-3;
alpha=45*(Pi/180); //radianos

l1 = 2*R/(Tan(Pi/2 - alpha));
l2 = R*Cos(alpha);
l3 = R*Sin(alpha);

// Parametro de malha
MeshFactor = 2e-4;

// Cylinder 
Point(1) = {0,0,0,MeshFactor}; //center
Point(2) = {0,0,D_inlet/2,MeshFactor};
Point(3) = {0,D_inlet/2,0,MeshFactor};
Point(4) = {0,-D_inlet/2,0,MeshFactor};

Circle(1) = {4,1,2};
Circle(2) = {2,1,3};
Line(3) = {3,1};
Line(4) = {1,4};

Line Loop(5) = {1,2,3,4};

Plane Surface(6) = {5};

Extrude {L_inlet,0,0} {
  Surface{6};
}

// //Side Wall outlet
Point(20) = {L_inlet+R+l2+l1,(D_inlet/2),W_outlet/2,MeshFactor};
Point(21) = {L_inlet,(D_inlet/2)+H_inlet,W_outlet/2,MeshFactor};
Point(22) = {L_inlet,(D_inlet/2)-R,W_outlet/2,MeshFactor};
Point(23) = {L_inlet+R+l2+l1+L_outlet,(D_inlet/2)+H_inlet,W_outlet/2,MeshFactor};
Point(25) = {L_inlet+R+l2+l1+L_outlet,(D_inlet/2),W_outlet/2,MeshFactor};
Point(26) = {L_inlet+R+l2,(D_inlet/2)-R-l3,W_outlet/2,MeshFactor};
Point(27) = {L_inlet+R,(D_inlet/2)-R,W_outlet/2,MeshFactor}; //center

Line(25) = {22,21};
Line(26) = {21,23};
Line(27) = {23,25};
Line(28) = {25,20};
Line(29) = {20,26};
Circle(30)= {26,27,22};

Line Loop(35) = {25:30};
Plane Surface(35) = (35);

Extrude {0,0,-W_outlet/2} { 
  Surface{35}; 
}



v() = BooleanFragments{ Volume{1}; Delete;}{ Volume{2}; Delete; };

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
Field[1].CurvesList = {17,20};//outlet lines

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 5e-5;
Field[2].SizeMax = MeshFactor;
Field[2].DistMin = L_outlet;
Field[2].DistMax =L_outlet+l1;

Field[3] = Distance;
Field[3].CurvesList = {11,3,6,9};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = MeshFactor / 4.1;
Field[4].SizeMax = MeshFactor;
Field[4].DistMin = L_inlet;
Field[4].DistMax =L_inlet*3;



Field[5] = Min;
Field[5].FieldsList = {2,4};
Background Field = 5;

Physical Surface("Inlet") = {5};
Physical Surface("Outlet") = {9};
Physical Surface("Wall") = {1,2,7,12,11,8,10,13};
Physical Surface("Symmetry") = {14,3,4};
Physical Volume("Fluid") = {1:2};





