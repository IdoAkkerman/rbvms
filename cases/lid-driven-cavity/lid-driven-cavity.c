//--------------------------------------------------------------
// Solution function for the lid driven cavity
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
   const double eps = 1e-6;

   u[0] = pow(4.0*x*(1.0-x),0.01)*y*y*y;
   u[1] = 0.0;

   if (y     < eps) u[0] = 0.0;
   if (x     < eps) u[0] = 0.0;
   if (1.0-x < eps) u[0] = 0.0;
   if (1.0-y < eps)
   {
      u[0] = 1.0;
      if (x     < eps) u[0] = 0.5;
      if (1.0-x < eps) u[0] = 0.5;
   }
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

