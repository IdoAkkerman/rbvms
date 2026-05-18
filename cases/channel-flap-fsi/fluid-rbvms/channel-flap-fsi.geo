SetFactory("OpenCASCADE");

// --------------------------------------------------
// Vertices from blockMeshDict
// --------------------------------------------------
X1 = -0.2; // pre rigid body
X2 =  0.0; // begin rigid body
X3 =  0.2; // end rigid body
X4 =  0.6; // wake

Y1 = -0.13; // distance to bottom
Y2 = -0.01; // half body height
Y3 =  0.01; // half body height
Y4 =  0.12; // distance to top

// --------------------------------------------------
// Geometry parameters
// --------------------------------------------------
L  = X4 - X1;    // channel length
H  = Y4 - Y1;    // channel height

Lf = X3 - X2;    // flap length
Hf = Y3 - Y2;    // flap thickness

lc = 0.01;    // mesh size

// --------------------------------------------------
// Channel rectangle
// --------------------------------------------------

Rectangle(1) = {X1, Y1, 0, L, H};

// --------------------------------------------------
// Flap rectangle
// --------------------------------------------------

Rectangle(2) = {X2, Y2, 0, Lf, Hf};

// --------------------------------------------------
// Fluid domain = channel minus flap
// --------------------------------------------------

BooleanDifference(3) = { Surface{1}; Delete; }
                       { Surface{2}; Delete; };

// --------------------------------------------------
// Physical groups
// --------------------------------------------------

// Fluid volume
Physical Surface("Fluid") = {3};

// --------------------------------------------------
// Boundary identification
// --------------------------------------------------

// After BooleanDifference, Gmsh creates new curves.
// Use the GUI once ("Tools -> Visibility") or
// `gmsh quickstart.geo` to inspect IDs if needed.
//
// The IDs below are the typical ones generated
// by OpenCASCADE BooleanDifference.
// --------------------------------------------------

// inlet
Physical Curve("Inlet") = {2};

// outlet
Physical Curve("Outlet") = {3};

// top wall
Physical Curve("TopWall") = {4};

// bottom wall
Physical Curve("BottomWall") = {1};

// flap interface
Physical Curve("FlapInterface") = {5, 6, 7, 8};

// --------------------------------------------------
// Mesh settings
// --------------------------------------------------

Mesh.Algorithm = 6;
Mesh.RecombineAll = 0;

// Optional refinement near flap
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc/5;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.02;
Field[2].DistMax = 0.15;

Background Field = 2;
