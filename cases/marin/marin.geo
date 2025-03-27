// Gmsh project created on Tue Mar  4 16:11:14 2025
SetFactory("OpenCASCADE");

//
Box(1) = {0, 0, 0, 3.22, 1, 1};
Box(2) = {2.3955, 0.2985, 0.0, 0.161, 0.403, 0.161};

//
Delete {
  Volume{1}; Volume{2};
  Surface{5}; Surface{11};
}

// Bottom
MeshSize{2, 4, 6 ,8} = 1.0;

// Top
MeshSize{1, 3, 5 ,7} = 2.0;

// Box
MeshSize{9,10,11,12,13,14,15,16} = 0.5;


// Bottom surface
Curve Loop(13) = {4, 11, -8, -9};
Curve Loop(14) = {16, 23, -20, -21};

Plane Surface(11) = {13,14};

// Volume
Surface Loop(100) = {1,2,3,4,6,7,8,9,10,11,12};
Volume(100) = {100};

Physical Surface(1) = {1,2,3,4,6,7,8,9,10,11,12};
Physical Volume(100) = {100};
