// Gmsh project created on Mon Jul 15 18:34:03 2019

Mesh.MshFileVersion = 2; // Version of the MSH file format to use
Lx0 = 0;
Lxf = 2;
Ly0 = 0;
Lyf = 0.2;
Lz0 = 0;
Lzf = 2;
nW = (Lyf-Ly0)*5;
nL = (Lxf-Lx0)*5;
MeshFactor = 0.05;

Point(1) = {Lx0, Ly0, Lz0,MeshFactor};
Point(2) = {Lxf, Ly0, Lz0,MeshFactor};
Point(3) = {Lxf, Lyf, Lz0,MeshFactor};
Point(4) = {Lx0, Lyf, Lz0,MeshFactor};
Point(5) = {Lx0, Ly0, Lzf,MeshFactor};
Point(6) = {Lxf, Ly0, Lzf,MeshFactor};
Point(7) = {Lxf, Lyf, Lzf,MeshFactor};
Point(8) = {Lx0, Lyf, Lzf,MeshFactor};

Line(1) = {1, 2}; 
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6}; 
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5}; 
Line(10) = {4, 8};
Line(11) = {6, 2};
Line(12) = {3, 7};

Line Loop(1) = {-1, -4, -3, -2};
Plane Surface(1) = (1);
Line Loop(2) = {1, -11, -5, -9};
Plane Surface(2) = (2);
Line Loop(3) = {4, 9, -8, -10};
Plane Surface(3) = (3);
Line Loop(4) = {2, 12, -6, 11};
Plane Surface(4) = (4);
Line Loop(5) = {3, 10, -7, -12};
Plane Surface(5) = (5);
Line Loop(6) = {5, 6, 7, 8};
Plane Surface(6) = (6);

Surface Loop(1) = {1,6,2,5,4,3};
Volume(1) = {1};


Physical Surface("Inlet") = {3};
Physical Surface("Outlet") = {4};
Physical Surface("Wall") = {2,5,1,6};
Physical Volume("Fluid") = {1};