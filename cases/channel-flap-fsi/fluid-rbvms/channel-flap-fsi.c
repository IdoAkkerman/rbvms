//--------------------------------------------------------------
// To compile run, for example:
// mpicc -shared -o channel-flap-fsi.so -fPIC channel-flap-fsi.c
// gcc   -shared -o channel-flap-fsi.so -fPIC channel-flap-fsi.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   // Base flow
   u[0] = 1.0;
   u[1] = 0.0;

   // BC
   double x = coord[0];
   double y = coord[1];
   const double eps = 0.01;
   if ((x > 0.0 - eps) && (x < 0.2 +eps) &&
       (y > -0.01 - eps) && (y < 0.01 +eps))
   {
      u[0] = 0.0;
   }
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
}

double mu(double *coord, int dim, double time)
{
   return 0.001;
}

