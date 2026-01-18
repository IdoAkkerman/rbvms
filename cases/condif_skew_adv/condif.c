//--------------------------------------------------------------
// Solution function for convection-diffusion
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC condif.c
// gcc   -shared -o libfun.so -fPIC condif.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

#ifndef ADVEC_ANGLE
#define ADVEC_ANGLE 30
#endif

#ifndef MU
#define MU 0.0001
#endif

double sol_phi(double *coord, int dim, double time)
{
   if (coord[0]+coord[1]<0.5)
      return 1.0;
   return 0.0;
}

void advection(double *coord, int dim, double time, double *adv, int vdim)
{
   adv[0] = cos(ADVEC_ANGLE*(M_PI/180));
   adv[1] = sin(ADVEC_ANGLE*(M_PI/180));
}

double force(double *coord, int dim, double time)
{
   return 0.0;
}

double mu(double *coord, int dim, double time)
{
   return MU;
}

