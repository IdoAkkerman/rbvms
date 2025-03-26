SetFactory("OpenCASCADE");

// Parameter
x0 =-1.0;
x1 = 7.0;
y0 = -1.0;
y1 = 1.0;

H=0.2;
L=2;

// Pointss
p1=newp; Point(p1) = {x0,  y0, 0.0, 0.2};
ps=newp; Point(ps)  = {0.0, y0, 0.0, 0.1};
pm1=newp; Point(pm1) = {L/3,   y0+H/L, 0.0, 0.1};
pm2=newp; Point(pm2) = {2*L/3, y0+H/(2*L), 0.0, 0.1};
pe=newp; Point(pe)  = {L, y0, 0.0, 0.1};
p2=newp; Point(p2) = {x1,  y0, 0.0, 0.3};
p3=newp; Point(p3) = {x1,  y1, 0.0, 0.5};
p4=newp; Point(p4) = {x0,  y1, 0.0, 0.4};

// Lines
l1s=newl; Line(l1s) = {p1, ps};
lse=newl; Spline(lse) = {ps, pm1, pm2, pe};
le2=newl; Line(le2) = {pe, p2};
l23=newl; Line(l23) = {p2, p3};
l34=newl; Line(l34) = {p3, p4};
l41=newl; Line(l41) = {p4, p1};

cl1=newcl; Curve Loop(cl1) = {l1s, lse, le2, l23, l34, l41};
s1=news; Plane Surface(s1)={cl1};

// Physical
Physical Line(1) = {l1s, lse, le2};
Physical Line(2) = {l23};
Physical Line(3) = {l34};
Physical Line(4) = {l41};
Physical Surface(100) = {s1};
