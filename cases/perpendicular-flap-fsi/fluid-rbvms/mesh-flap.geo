// Gmsh project created on Tue Jun 16 15:31:47 2026
SetFactory("OpenCASCADE");

// =====================================================
// Geometry parameters
// =====================================================

x0 = -3.0;
x1 = -0.05;
x2 =  0.05;
x3 =  3.0;

y0 = 0.0;
y1 = 1.0;
y2 = 4.0;

zlen = 1.0;

// =====================================================
// Mesh parameters (from blockMesh)
// =====================================================

h1 = 30;
h2 = 3;

v1 = 15;
v2 = 30;

refine = 1;

// =====================================================
// Points
// =====================================================

Point(1)  = {x0,y0,0};
Point(2)  = {x1,y0,0};
Point(3)  = {x2,y0,0};
Point(4)  = {x3,y0,0};

Point(5)  = {x0,y1,0};
Point(6)  = {x1,y1,0};
Point(7)  = {x2,y1,0};
Point(8)  = {x3,y1,0};

Point(9)  = {x0,y2,0};
Point(10) = {x1,y2,0};
Point(11) = {x2,y2,0};
Point(12) = {x3,y2,0};

// =====================================================
// Curves
// =====================================================

// bottom
Line(1) = {1,2};
Line(2) = {3,4};

// y = 1
Line(3) = {5,6};
Line(4) = {6,7};
Line(5) = {7,8};

// top
Line(6) = {9,10};
Line(7) = {10,11};
Line(8) = {11,12};

// verticals
Line(9)  = {1,5};
Line(10) = {5,9};

Line(11) = {2,6};
Line(12) = {6,10};

Line(13) = {3,7};
Line(14) = {7,11};

Line(15) = {4,8};
Line(16) = {8,12};

// flap edges
Line(18) = {2,6}; // left flap side (same geometry as 11)
Line(19) = {6,7}; // flap top (same geometry as 4)
Line(20) = {7,3}; // right flap side

// =====================================================
// Surfaces (5 OpenFOAM blocks)
// =====================================================

// lower left
Curve Loop(1) = {1,11,-3,-9};
Plane Surface(1) = {1};

// lower right
Curve Loop(2) = {2,15,-5,-13};
Plane Surface(2) = {2};

// upper left
Curve Loop(3) = {3,12,-6,-10};
Plane Surface(3) = {3};

// upper center
Curve Loop(4) = {4,14,-7,-12};
Plane Surface(4) = {4};

// upper right
Curve Loop(5) = {5,16,-8,-14};
Plane Surface(5) = {5};

// =====================================================
// Structured mesh
// =====================================================

// horizontal grading

Transfinite Curve {1,3,6} = refine*(h1 + 1) Using Progression 0.95;
Transfinite Curve {2,5,8} = refine*(h1 + 1) Using Progression 1.05;
Transfinite Curve {4,7}   = refine*(h2 + 1) Using Progression 1.0;

// vertical grading

Transfinite Curve {9,11,13,15} = refine*(v1 + 1) Using Bump 0.8;

Transfinite Curve {10,12,14,16} = refine*(v2 + 1) Using Progression 1.02;

// transfinite surfaces

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};

Recombine Surface {1,2,3,4,5};

// =====================================================
// Physical groups
// =====================================================

//+
Physical Curve("inlet", 1) = {10, 9};
//+
Physical Curve("outlet0", 20) = {15};
//+
Physical Curve("outlet1", 21) = {16};
//+
Physical Curve("flap", 3) = {11, 4, 13};
//+
Physical Curve("lowerWall", 4) = {1, 2};
//+
Physical Curve("upperWall", 5) = {6, 7, 8};
//+
Physical Surface("fluid", 1) = {3, 4, 5, 2, 1};
