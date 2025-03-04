//--------------------------------------------------------------
// Solution function for the MARIN dambreak problem
//
// To compile run, for example:
// mpicc -shared -o rbvms-ls.so -fPIC marin.c
// gcc   -shared -o rbvms-ls.so -fPIC marin.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>


void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   u[0] = 0.0;
   u[1] = 0.0;
   u[2] = 0.0;
}

double sol_phi(double *coord, int dim, double time)
{
   double x = coord[0] - 1.0;
   double z = coord[2] - 0.55;
   double r = sqrt(x*x + z*z);
   //return fmin(-x,-z);
   return ((x > 0) && (z >0)) ? -r : fmin(-x,-z);
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
   force[2] = 0.0;
}

double mu(double *coord, int dim, double time)
{
   double mu_f = 0.01;
   double mu_0 = 0.01;
   double s   = -log(mu_f/mu_0)/1.0;

   return fmax(mu_f, mu_0 *exp(-s*time));
}

