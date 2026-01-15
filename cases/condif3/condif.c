//--------------------------------------------------------------
// Solution function for convection-diffusion
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC condif.c
// gcc   -shared -o libfun.so -fPIC condif.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

#ifndef MU
#define MU 0.0
#endif

double sol_phi(double *coord, int dim, double time)
{
   double r = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);

   return 0.5*(1.0 - cos(2*M_PI*fmin(1.0, r)));
}

void advection(double *coord, int dim, double time, double *adv, int vdim)
{
   adv[0] = coord[1];
   adv[1] = -coord[0];
}

double force(double *coord, int dim, double time)
{
   return 0.0;
}

double mu(double *coord, int dim, double time)
{
   return MU;
}

