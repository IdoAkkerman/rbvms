// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "tau.hpp"

using namespace mfem;
using namespace RBVMS;

/* Evaluate the stabilisation parameter at @a ip.
    - tau[0] approximates the momentum operator and
      is used for velocity reconstruction.
    - tau[1] approximates the pressure poisson operator and
      is used for pressure reconstruction.
    - tau[2] reports the CFL number.*/
void Tau::Eval(Vector &tau,
               ElementTransformation &T,
               const IntegrationPoint &ip)
{
   real_t mu = c_mu.Eval(T, ip);
   c_adv.Eval(u, T, ip);

   MultAtB(T.InverseJacobian(),T.InverseJacobian(),Gij);

   real_t tau_t = Ct/(dt*dt);
   real_t tau_c = 0.0;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_c += Gij(i,j)*u[i]*u[j];
      }
   }
   real_t tau_d = 0.0;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_d  += Cd*Gij(i,j)*Gij(i,j)*mu*mu;
      }
   }

   tau[0] = 1.0/sqrt(tau_t + tau_c + tau_d);
   tau[1] = 1.0/(tau[0]*Gij.Trace());
   tau[2] = sqrt(tau_c/tau_t);
}
