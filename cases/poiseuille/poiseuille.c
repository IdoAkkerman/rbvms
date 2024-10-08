//--------------------------------------------------------------
// Solution function for the von Karman vortex street
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC lid-driven-cavity.c
// gcc   -shared -o libfun.so -fPIC lid-driven-cavity.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0];
   double y = coord[1];

   u[0] = 4*y*(1.0-y);
   u[1] = 0.0;
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 16.0;
   force[1] = 0.0;
}

double mu(double *coord, int dim, double time)
{
   return 1.0;
}

