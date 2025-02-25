//--------------------------------------------------------------
// Function for foil simulation
//
// To compile run, for example:
// mpicc -shared -o foil.so -fPIC foil.c
// gcc   -shared -o foil.so -fPIC foil.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

// Solution: IC & BC
// t = -1.0 at IC
void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0]-1.15;
   double y = coord[1]-0.65;
   double r = sqrt(x*x + y*y);
   double uinf = 10.0;

   // Base flow + BC
   u[0] = (r < 1.5 && time > 0.0) ?  0.0 : uinf;
   u[1] = 0.0;
}

// Force: body force in N/m^3
void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
}

// Suction: pressure in Pa
double suction(double *coord, int dim, double time)
{
   return -250.0;
}

// Blowing: normal velocity in m/s
double blowing(double *coord, int dim, double time)
{
   return 1.978*fmin(1.0, 10*exp(10.0*time/1.0));
}

// Density of air in kg/m^3 at 273 K (20C) & 1.01325 bar (1 atm)
double rho(double *coord, int dim, double time)
{
   return 1.204;
}

// Dynamic viscosity of air kg/ms at 273 K (20C) & 1.01325 bar (1 atm)
double mu(double *coord, int dim, double time)
{
   return 1.825e-5*fmax(1.0, 1e4*exp(-100.0*time/1.0));
}
