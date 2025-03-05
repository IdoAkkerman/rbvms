// Parameter
x = 3.22;
y = 1.0;
z = 1.0;

// Pointss
p1=newp; Point(p1) = {0.0, 0.0, 0.0, 1.0};
p2=newp; Point(p2) = {  x, 0.0, 0.0, 1.0};
p3=newp; Point(p3) = {  x,   y, 0.0, 1.0};
p4=newp; Point(p4) = {0.0,   y, 0.0, 1.0};

p5=newp; Point(p5) = {0.0, 0.0, z, 1.0};
p6=newp; Point(p6) = {  x, 0.0, z, 1.0};
p7=newp; Point(p7) = {  x,   y, z, 1.0};
p8=newp; Point(p8) = {0.0,   y, z, 1.0};

// Lines
l11=newl; Line(l11) = {p1, p2};
l12=newl; Line(l12) = {p2, p3};
l13=newl; Line(l13) = {p3, p4};
l14=newl; Line(l14) = {p4, p1};

l21=newl; Line(l21) = {p1, p5};
l22=newl; Line(l22) = {p2, p6};
l23=newl; Line(l23) = {p3, p7};
l24=newl; Line(l24) = {p4, p8};

l31=newl; Line(l31) = {p5, p6};
l32=newl; Line(l32) = {p6, p7};
l33=newl; Line(l33) = {p7, p8};
l34=newl; Line(l34) = {p8, p5};

// Curve Loop
cl1=newcl; Curve Loop(cl1) = {l11, l12, l13, l14};
cl2=newcl; Curve Loop(cl2) = {l11, l22, -l31, -l21};
cl3=newcl; Curve Loop(cl3) = {l12, l23, -l32, -l22};
cl4=newcl; Curve Loop(cl4) = {l13, l24, -l33, -l23};
cl5=newcl; Curve Loop(cl5) = {l14, l21, -l34, -l24};
cl6=newcl; Curve Loop(cl6) = {l31, l32, l33, l34};

// Surface
s1=news; Plane Surface(s1)={cl1};
s2=news; Plane Surface(s2)={cl2};
s3=news; Plane Surface(s3)={cl3};
s4=news; Plane Surface(s4)={cl4};
s5=news; Plane Surface(s5)={cl5};
s6=news; Plane Surface(s6)={cl6};

// Volume
Surface Loop(1) = {s1, s2, s3, s4, s5, s6};
Volume(1) = {1};


Physical Surface(1) = {s1, s2, s3, s4, s5, s6};
Physical Volume(100) = {1};
