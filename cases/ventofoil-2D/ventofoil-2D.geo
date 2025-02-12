//Ventofoil


//Notes-----------------------------------------------------
//pay attention to the origin of the geometry of stepfile, depending on it it appears in gmsh
//a step file or a iges file can be included
//just export the surfaces (no curves, volumes, points) of the geometry from Rhino with origin 


//Begin -------------------------
//counts all points, lines, surfaces present at start
start_point_list[]={Point{:}};
start_line_list[]={Line{:}};
start_surface_list[]={Surface{:}};
start_point_amount= #start_point_list[];
start_line_amount= #start_line_list[];
start_surface_amount= #start_surface_list[];
Printf("FOIL sizes start point, line, surface",start_point_amount,start_line_amount, start_surface_amount);
SetFactory("OpenCASCADE"); //makes translation and rotation possible
Geometry.OCCTargetUnit='M';
// Import file ----------------------------------------------
//***************Things to adapt ***************
Merge "VentoFoil-2D.stp"; //change here the filename
foilamount =1; //input here the number of Foils (right now max. 6, include more translate-rotate-lines and coordinates if more foils)

//Translation and Rotation of 1st VentoFoil
x1=00; //x coordinate of foil
y1=0; //y coordinate of foil
z1=0; //z coordinate of foil from waterplane area
angle1= 1/6*Pi; //rotation angle in rad (positive anticlockwise) 
//***************** End adapting ***************


Geometry.OCCSewFaces=1; //Heals model that super close lines are identical
//Counting of points, lines, surfaces after import
import_point_list[]={Point{:}};
import_point_amount= #import_point_list[];
import_line_list[]={Line{:}};
import_line_amount= #import_line_list[];
import_surface_list[]={Surface{:}};
import_surface_amount= #import_surface_list[];
Printf("FOIL import point, lines, surfaces", import_point_amount, import_line_amount, import_surface_amount);

//surface loop
cl10=newcl; Curve Loop(cl10)={start_line_amount+1:import_line_amount};
s0=news; Surface(s0)={cl10};
Rotate {{0,0,1}, {x1,y1,z1}, angle1} { Surface{s0};} 


Delete{Line{start_line_amount+2:import_line_amount};} 
Delete{Point{:};} 
Delete{Surface{s0};}

end_point_list[]={Point{:}};
end_point_amount= #end_point_list[];
end_line_list[]={Line{:}};
end_line_amount= #end_line_list[];
end_surface_list[]={Surface{:}};
end_surface_amount= #end_surface_list[];
Printf("FOIL end amount point, line, surface",end_point_amount, end_line_amount, end_surface_amount);

