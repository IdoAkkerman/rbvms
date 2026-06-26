//--------------------------------------------------------------
// Solution function for the von Karman vortex street
//
// To compile run, for example:
// mpicc -shared -o libvkvs.so -fPIC vkvs.c
// gcc   -shared -o libvkvs.so -fPIC vkvs.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0];
   double y = coord[1];
   double r = sqrt(x*x + y*y);
   double r0 = 0.5;

   // Base flow + BC
   u[0] = fmin(1.0, r - r0);
   u[1] = 0.0;

   // Disturbance 1
   x = coord[0];
   y = coord[1] - 3.0;
   r = sqrt(x*x + y*y);
   r0 = 1.0;
   u[0] += 2.0*fmax(0.0, r0 - r);
   u[1] -= 0.5*fmax(0.0, r0 - r);

   // Disturbance 2
   x = coord[0] + 3.0;
   y = coord[1];
   r = sqrt(x*x + y*y);
   r0 = 1.0;
   u[0] -= 0.5*fmax(0.0, r0 - r);
   u[1] += 2.0*fmax(0.0, r0 - r);
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
}

double mu(double *coord, int dim, double time)
{
   double mu_f = 0.01;
   double mu_0 = 0.01;
   double s   = -log(mu_f/mu_0)/1.0;

   return fmax(mu_f, mu_0 *exp(-s*time));
}

