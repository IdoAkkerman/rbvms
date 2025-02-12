// Ventofoils in box

Geometry.OldNewReg=0; //makes numbering continious

//************************Things to adapt*******************************************************
//Import the Ventofoils
Include "ventofoil-2D.geo";






foilamount=1; //insert number of foils (must be matching the amount in Ventofoil.txt)

//accuracies domain
lcb=30; //mesh size for box

//domain extensions
xb=-6; //x coordinate of 1 corner point of box
yb =-6; //y coordinate of 1 corner point of box
zb=0; //z coordinate of 1 corner point of box
dxb=32; //(dx) length of box
dyb=12; //(dy) width of box
angle=0*Pi; //rotation angle of domain in rad (positive anticlockwise)
//************************End Adapting***********************************************************

//Domain definition -------------------------------------------------
//corner points of domain
p901=newp; Point(p901)={xb,yb,zb,lcb};
p902=newp; Point(p902)={xb+dxb,yb,zb,lcb};
p903=newp; Point(p903)={xb+dxb,yb+dyb,zb,lcb};
p904=newp; Point(p904)={xb,yb+dyb,zb,lcb};

//edges of domain
l901=newl; Line(l901)={p901,p902};
l902=newl; Line(l902)={p902,p903};
l903=newl; Line(l903)={p903,p904};
l904=newl; Line(l904)={p904,p901};

//edge loops of planes of domain
cl901=newcl; Curve Loop(cl901)={l901,l902,l903,l904};

//Surface with foil cut out
s900=news; Plane Surface(s900)={cl901,2};

//count amount of points, lines, surfaces
end_point_list[]={Point{:}};
end_point_amount= #end_point_list[];
end_line_list[]={Line{:}};
end_line_amount= #end_line_list[];
end_surface_list[]={Surface{:}};
end_surface_amount= #end_surface_list[];
Printf("FOIL WITH DOMAIN end amount point, line, surface",end_point_amount, end_line_amount, end_surface_amount);

//Physical elements domain
Physical Line(1) = {l904}; //inflow
Physical Line(2) = {l901, l903}; //sides
Physical Line(3) = {l902}; //outflow

//Physical elements VentoFoil
Physical Line(10)={8,9,11,12}; //no slip at foil
Physical Line(11)={7};  // suction
Physical Line(12)={10}; // blowing

Physical Surface(100) = {s900}; //domain with foil cut out

//Mesh.MeshSizeFactor = 0.1;
Mesh.Smoothing = 100;


// Mesh
//Field[1] = MathEval;
//Field[1].F = "(0.1 + ((7-0.2*x)/64)*y*y) + (1.3-tanh(0.1*(x+8))+0.02*x)";




Field[1] = Distance;
Field[1].CurvesList = {7,8,9,10,11,12};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcb / 30;
Field[2].SizeMax = lcb / 3;
Field[2].DistMin = 0.1;
Field[2].DistMax = 1.0;


//rp1 = newp; Point(rp1)={2.5,1.4,0,lcb};
rp1 = newp; Point(rp1)={0,0,0,lcb};
rp2 = newp; Point(rp2)={11,6,0,lcb};
rl=newl; Line(rl)={rp1,rp2};

Field[3] = Distance;
Field[3].CurvesList = {rl};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lcb / 10;
Field[4].SizeMax = lcb ;
Field[4].DistMin = 2;
Field[4].DistMax = 3.0;

//edges of domain

Field[5] = Threshold;
Field[5].InField = 3;
Field[5].SizeMin = lcb / 15;
Field[5].SizeMax = lcb / 3;
Field[5].DistMin = 0.5;
Field[5].DistMax = 1.5;

Field[6] = Min;
Field[6].FieldsList = {2, 4, 5};
Background Field = 6;

