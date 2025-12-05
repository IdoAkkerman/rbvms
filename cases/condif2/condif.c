//--------------------------------------------------------------
// Solution function for convection-diffusion
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC condif.c
// gcc   -shared -o libfun.so -fPIC condif.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>


double sol_phi(double *coord, int dim, double time)
{
   return 0.0;
}

void advection(double *coord, int dim, double time, double *adv, int vdim)
{
   double r = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);

   adv[0] = -coord[0]/r;//-sqrt(2.0)/2.0;
   adv[1] = -coord[1]/r;//-sqrt(2.0)/2.0;
}

double force(double *coord, int dim, double time)
{
   return 1.0;
}

double mu(double *coord, int dim, double time)
{
   return 0.00001;
}

