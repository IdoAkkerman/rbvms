//--------------------------------------------------------------
// Solution function for convection-diffusion -- manufactured solution
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC condif_ms.c
// gcc   -shared -o libfun.so -fPIC condif_ms.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void advection(double *x, int dim, double t, double *adv, int vdim)
{
   adv[0] = sqrt(2.0)/2.0;
   adv[1] = sqrt(2.0)/2.0;
}

double mu(double *x, int dim, double t)
{
   return 0.01;
}

double sol_phi(double *x, int dim, double t)
{
   return sin(M_PI*x[0])*sin(M_PI*x[1]);
}

void grad_phi(double *x, int dim, double t, double *grad, int vdim)
{
   grad[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
   grad[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}

double lap_phi(double *x, int dim, double t)
{
   return -2*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
}

double force(double *x, int dim, double t)
{
   double a[dim];
   double grad[dim];
   advection(x, dim, t, a, dim);
   grad_phi(x, dim, t, grad, dim);

   double res = - mu(x, dim, t)*lap_phi(x, dim, t);
   for (int i = 0; i < dim; i++)
   {
      res += a[i]*grad[i];
   }

   return res;
}
