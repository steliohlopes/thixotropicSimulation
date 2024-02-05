Mesh.MshFileVersion = 2; // Version of the MSH file format to use

// Parametros de geometria
L = 24e-3;
l = 4e-3;
R = 100e-6;
r = 50e-6;
r0 = -R;
h = 8*R;
// refinamento da garganta
Div=100;
dx=l/Div;

// Parametro de malha
MeshFactor = 3e-5;


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


Extrude {0,0,h}{
    Surface{1};
}







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

// Field[1] = Distance;
// Field[1].CurvesList = {4:102,106:204};
// Field[1].Sampling = 20;

// Field[2] = Threshold;
// Field[2].InField = 1;
// Field[2].SizeMin = MeshFactor / 2;
// Field[2].SizeMax = MeshFactor;
// Field[2].DistMin = r*3;
// Field[2].DistMax = r*40;

// Field[3] = Distance;
// Field[3].CurvesList = {4};
// Field[3].Sampling = 100;

// Field[4] = Threshold;
// Field[4].InField = 3;
// Field[4].SizeMin = MeshFactor / 3;
// Field[4].SizeMax = MeshFactor;
// Field[4].DistMin = r/10;
// Field[4].DistMax = r*50;


// Field[5] = Min;
// Field[5].FieldsList = {2};
// Background Field = 5;


Physical Volume("Fluid") = {1};
//+
Physical Surface("Inlet") = {825};
//+
Physical Surface("Outlet") = {417};
//+
Physical Surface("Wall") = {1226, 1, 413, 977, 929, 933, 937, 941, 945, 949, 953, 957, 961, 965, 969, 973, 925, 981, 985, 989, 993, 997, 1001, 1005, 1009, 1013, 1017, 1021, 1025, 877, 829, 833, 837, 841, 845, 849, 853, 857, 861, 865, 869, 873, 1029, 881, 885, 889, 893, 897, 901, 905, 909, 913, 917, 921, 1181, 1133, 1137, 1141, 1145, 1149, 1153, 1157, 1161, 1165, 1169, 1173, 1177, 1129, 1185, 1189, 1193, 1197, 1201, 1205, 1209, 1213, 1217, 1221, 1225, 1081, 1033, 1037, 1041, 1045, 1049, 1053, 1057, 1061, 1065, 1069, 1073, 1077, 821, 1085, 1089, 1093, 1097, 1101, 1105, 1109, 1113, 1117, 1121, 1125, 573, 525, 529, 533, 537, 541, 545, 549, 553, 557, 561, 565, 569, 521, 577, 581, 585, 589, 593, 597, 601, 605, 609, 613, 617, 473, 425, 429, 433, 437, 441, 445, 449, 453, 457, 461, 465, 469, 621, 477, 481, 485, 489, 493, 497, 501, 505, 509, 513, 517, 773, 725, 729, 733, 737, 741, 745, 749, 753, 757, 761, 765, 769, 721, 777, 781, 785, 789, 793, 797, 801, 805, 809, 813, 817, 673, 625, 629, 633, 637, 641, 645, 649, 653, 657, 661, 665, 669, 421, 677, 681, 685, 689, 693, 697, 701, 705, 709, 713, 717};
