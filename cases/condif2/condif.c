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
#define MU 0.00001
#endif

double sol_phi(double *coord, int dim, double time)
{
   return 0.0;
}

void advection(double *coord, int dim, double time, double *adv, int vdim)
{
   adv[0] = 0;
   adv[1] = -1;
}

double force(double *coord, int dim, double time)
{
   return 1.0;
}

double mu(double *coord, int dim, double time)
{
   return MU;
}

