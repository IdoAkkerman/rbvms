//--------------------------------------------------------------
// Solution function for the flow around a naca foil
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC naca-foil.c
// gcc   -shared -o libfun.so -fPIC naca-foil.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0] - 0.5;
   double y = coord[1];
   const double eps = 1e-6;

   double r = sqrt(x*x + 10*y*y);
   double r0 = 0.65;

   // Base flow + BC
   u[0] = fmin(1.0, 5*fmax(r - r0, 0.0));
   u[1] = 0.0;
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
}

double mu(double *coord, int dim, double time)
{
   return 0.01;
}

