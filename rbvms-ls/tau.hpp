// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_TAU_HPP
#define RBVMS_TAU_HPP

#include "mfem.hpp"

using namespace mfem;

namespace RBVMS
{

/** This Class defines the stabilisation parameter for reconstructing the
    small-scale velocity.
     - tau[0] approximates the momentum operator and
       is used for velocity reconstruction.
     - tau[1] approximates the pressure poisson operator and
       is used for pressure reconstruction.
*/
class Tau: public VectorCoefficient
{
protected:
   /// The advection field
   VectorCoefficient &c_adv;

   /// The diffusion parameter field
   Coefficient &c_mu;
   real_t dt = -1.0;

   /// Constants used in the stabilisation parameter
   real_t Cd, Ct;

   /// Dimension of the problem
   int dim;

   /// Velocity vector
   Vector u;

   // Metric tensor
   DenseMatrix Gij;

public:
   /** Construct a stabilized confection-diffusion integrator with:
       - @a adv the convection velocity.
       - @a mu the diffusion coefficient.
       - @a Cd the constant used to scale the diffusive part
       - @a Ct the constant used to scale the temporal part
   */
   Tau(VectorCoefficient &adv, Coefficient &mu,
       real_t Cd = 12.0, real_t Ct = 1.0)
      : VectorCoefficient(3), c_adv(adv), c_mu(mu),
        Cd(Cd), Ct(Ct)
   {
      dim = c_adv.GetVDim();
      u.SetSize(dim);
      Gij.SetSize(dim);
   };

   /// Set the timestep @a dt
   void SetTimeStep(const double dt_) { dt = dt_;};

   /** Evaluate the stabilisation parameter at @a ip.
     - tau[0] approximates the momentum operator and
       is used for velocity reconstruction.
     - tau[1] approximates the pressure poisson operator and
       is used for pressure reconstruction.*/
   virtual void Eval(Vector &tau,
                     ElementTransformation &T,
                     const IntegrationPoint &ip) override;

};

} // namespace RBVMS

#endif
