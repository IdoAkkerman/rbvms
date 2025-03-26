//--------------------------------------------------------------
// Solution function for the Cahouet problem
//
// To compile run, for example:
// mpicc -shared -o rbvms-ls.so -fPIC cahouet.c
// gcc   -shared -o rbvms-ls.so -fPIC cahouet.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>


void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0];
   double y = coord[1];

   u[0] = 0.43*sqrt(9.81);
   u[1] = 0.0;

   if (time < 0.0)
   {
      return;
   }
   else if (y < -0.999)
   {
      u[0] = 0.0;
   }
   else if (((x > -0.1) && (x < 2.1)) && (y < -0.205))
   {
      u[0] = 0.0;
   }
}

double sol_phi(double *coord, int dim, double time)
{
   return -coord[1];
}


double sol_p(double *coord, int dim, double time)
{
   double p;

   if (coord[1] > 0.0)
   {
      p = -9.81*1*coord[1];
   }
   else
   {
      p = -9.81*1000*coord[1];
   }

   return p;
}






void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = -9.81;
}

double suction(double *coord, int dim, double time)
{
   double h;

   if (coord[1] > 0.0)
   {
      h = -9.81*1*coord[1];
   }
   else
   {
      h = -9.81*1000*coord[1];
   }

   return h;
}


double mu(double *coord, int dim, double time)
{
   return 0.01;
}
