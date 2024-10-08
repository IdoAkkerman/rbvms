//--------------------------------------------------------------
// Solution function for the von Karman vortex street
//
// To compile run, for example:
// mpicc -shared -o libfun.so -fPIC channel.c
// gcc   -shared -o libfun.so -fPIC channel.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0];
   double y = coord[1];
   double z = coord[2];

   u[0] = 20.0*(2.0 + sin(2.0*M_PI*x/10.0) + + sin(2.0*M_PI*y/2.0))*z*(1.0-z);
   u[1] = sin(2.0*M_PI*2.0*z);
   u[2] = sin(2.0*M_PI*2.0*z);
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   double x = coord[0];
   double y = coord[1];
   double z = coord[2];

   force[0] = 1.0;
   force[1] = 0.1*sin(2.0*M_PI*2.0*z)*sin(time)*exp(-0.1*time);
   force[2] = 0.5*sin(2.0*M_PI*2.0*z)*cos(time)*exp(-0.1*time);
}

double mu(double *coord, int dim, double time)
{
   return 1.0/180.0;
}

