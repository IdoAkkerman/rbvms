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

/// This Class defines a generic stabilisation parameter.
class Tau: public Coefficient
{
protected:
   /// The advection field
   VectorCoefficient *adv;
   /// The diffusion parameter field
   Coefficient *kappa;
   /// Dimension of the problem
   int dim;
   /// Velocity vector
   Vector a;

public:
   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient*/
   Tau (VectorCoefficient *a_, Coefficient *k) : adv(a_), kappa(k)
   {
      dim = adv->GetVDim();
      a.SetSize(dim);
   };
   /// Simple Constructor
   Tau () : adv(nullptr), kappa(nullptr) {};

   /// Set the convection coefficient
   virtual void SetConvection(VectorCoefficient *a_)
   {
      adv = a_;
      dim = adv->GetVDim();
      a.SetSize(dim);
   };

   /// Set the convection coefficient
   virtual void SetDiffusion(Coefficient *k_) { kappa = k_; };

   /// Flag for printing
   bool print = true;
};

/** This Class defines the stabilisation parameter for the multi-dimensional
    convection-diffusion problem.
    When @a k_fac =2 the parameter is defined as given in:

    Franca, L.P., Frey, S.L., & Hughes, T.J.R.
    Stabilized finite element methods:
    I. Application to the advective-diffusive model.
    Computer Methods in Applied Mechanics and Engineering, 95(2), 253-276.

    This also works for the convection-diffusion part of the navier-Stokes problem.
    When @a k_fac = 4 the parameter is defined as given in:

    Franca, L.P., Frey, S.L.,
    Stabilized finite element methods:
    II. The incompressible Navier-Stokes equations.
    Computer Methods in Applied Mechanics and Engineering, 99(2-3), 209-233.
*/
class FFH92Tau: public Tau
{
protected:
   /// User provided inverse estimate of the elements
   real_t Ci = -1.0;

   /// If @a Ci is negative the is inverse estimate computed
   Coefficient *invEst_cf = nullptr;
   bool own_ie = false;

   // Routine to get the inverse estimate at each point
   real_t GetInverseEstimate(ElementTransformation &T,
                             const IntegrationPoint &ip, real_t scale = 1.0);

   /// The norm used for the velocity vector
   real_t p = 2.0;

   /// The facor used for computing the element Peclet/Reynolds number
   real_t k_fac = 2.0;

   /// Temp variable
   Vector row;

   /// Element size in different directions
   Vector h;

   /** Returns element size according to:

       Harari, I, & Hughes, T.J.R.
       What are C and h?: Inequalities for the analysis and design of
       finite element methods.
       Computer methods in applied mechanics and engineering 97(2), 157-192.
   */
   real_t GetElementSize(ElementTransformation &T);

public:
   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a ie_cf for computing the inverse estimates
       - @a f factor for computing the element Pe/Re number (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2) */
   FFH92Tau (VectorCoefficient *a, Coefficient *k, Coefficient *ie_cf,
             real_t f = 2.0, real_t norm_p = 2.0)
      :  Tau(a,k), invEst_cf(ie_cf), Ci(-1.0), k_fac(f), p(norm_p) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a fes to provide to coefficient for computing the inverse estimates
       - @a f factor for computing the element Pe/Re number (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2) */
   FFH92Tau (VectorCoefficient *a, Coefficient *k,
             FiniteElementSpace *fes,
             real_t f = 2.0, real_t norm_p = 2.0)
      :  Tau(a,k), Ci(-1.0), k_fac(f), p(norm_p)
   {
      invEst_cf = NULL;//new InverseEstimateCoefficient(fes);
      own_ie = true;
   };

   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a c_explicity provided inverse estimate (default = 1.0/12.0)
       - @a f factor for computing the element Pe/Re number (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)*/
   FFH92Tau (VectorCoefficient *a, Coefficient *k,
             real_t c_ = 1.0/12.0, real_t f = 2.0, real_t norm_p = 2.0)
      :  Tau(a,k), Ci(c_), k_fac(f), p(norm_p)
   {
      invEst_cf = NULL;
   };

   /** Construct a stabilized confection-diffusion integrator with:
       - @a ie_cf for computing the inverse estimatestes
       - @a f factor for Pe/Re definition (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FFH92Tau (Coefficient *ie_cf,
             real_t f = 2.0, real_t norm_p = 2.0)
      : invEst_cf(ie_cf), Ci(-1.0), k_fac(f), p(norm_p) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a fes to provide to coefficient for computing the inverse estimates
       - @a f factor for Pe/Re definition (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FFH92Tau (FiniteElementSpace *fes,
             real_t f = 2.0, real_t norm_p = 2.0)
      : Ci(-1.0), k_fac(f), p(norm_p)
   {
      invEst_cf = NULL;//new InverseEstimateCoefficient(fes);
      own_ie = true;
   };

   /** Construct a stabilized confection-diffusion integrator with:
       - @a c_explicity provided inverse estimate (default = 1.0/12.0)
       - @a f factor for computing the element Pe/Re number (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FFH92Tau (real_t c_ = 1.0/12.0, real_t f = 2.0, real_t norm_p = 2.0)
      : Ci(c_), k_fac(f), p(norm_p)
   {
      invEst_cf = NULL;
   };

   /// Evaluate the coefficient at @a ip.
   virtual real_t Eval(ElementTransformation &T,
                       const IntegrationPoint &ip) override;

   // Destructor
   ~FFH92Tau()
   { if (own_ie) { delete invEst_cf; } }
};

/** This Class defines the stabilisation parameter for the multi-dimensional
    convection-diffusion problem.

    This also works for the convection-diffusion part of the navier-Stokes problem.

    Franca, L.P., Frey, S.L.,
    Stabilized finite element methods:
    II. The incompressible Navier-Stokes equations.
    Computer Methods in Applied Mechanics and Engineering, 99(2-3), 209-233.
*/
class FF91Delta: public FFH92Tau
{
protected:
   /// Overall scalling parameter
   real_t lambda = 1.0;

public:
   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a ie_cf for computing the inverse estimates
       - @a f factor for computing the element Pe/Re number (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2) */
   FF91Delta (VectorCoefficient *a, Coefficient *k, Coefficient *ie_cf,
              real_t l = 1.0, real_t f = 2.0, real_t norm_p = 2.0)
      :  FFH92Tau(a,k,ie_cf,f,norm_p), lambda(l) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a fes for computing the inverse estimates
       - @a l overall scalling factor for delta (default = 1)
       - @a f factor for computing the element Pe/Re number (default = 4)
       - @a p which norm to use for the velocity magnitude  (default = 2) */
   FF91Delta (VectorCoefficient *a, Coefficient *k,
              FiniteElementSpace *fes,
              real_t l = 1.0,
              real_t f = 4.0, real_t norm_p = 2.0)
      :  FFH92Tau(a,k,fes,f,norm_p), lambda(l) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a a the convection velocity
       - @a d the diffusion coefficient
       - @a c_explicity provided inverse estimate (default = 1.0/12.0)
       - @a l overall scalling factor for delta (default = 1)
       - @a f factor for computing the element Pe/Re number (default = 4)
       - @a p which norm to use for the velocity magnitude (default = 2)*/
   FF91Delta (VectorCoefficient *a, Coefficient *k,
              real_t c_ = 1.0/12.0, real_t l = 1.0,
              real_t f = 4.0, real_t norm_p = 2.0)
      : FFH92Tau(a,k,c_,f,norm_p), lambda(l) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a ie_cf for computing the inverse estimatestes
       - @a f factor for Pe/Re definition (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FF91Delta (Coefficient *ie_cf,
              real_t l = 1.0, real_t f = 2.0, real_t norm_p = 2.0)
      : FFH92Tau(ie_cf,f,norm_p), lambda(l) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a fes for computing the inverse estimates
       - @a f factor for Pe/Re definition (default = 2)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FF91Delta (FiniteElementSpace *fes,
              real_t l = 1.0, real_t f = 4.0, real_t norm_p = 2.0)
      : FFH92Tau(fes,f,norm_p), lambda(l) {};

   /** Construct a stabilized confection-diffusion integrator with:
       - @a c_explicity provided inverse estimate (default = 1.0/12.0)
       - @a l overall scalling factor for delta (default = 1)
       - @a f factor for computing the element Pe/Re number (default = 4)
       - @a p which norm to use for the velocity magnitude (default = 2)
       Convection and diffusion need to be specified later using
       SetConvection and SetDiffusion, respectivly*/
   FF91Delta (real_t c_ = 1.0/12.0, real_t l = 1.0,
              real_t f = 2.0, real_t norm_p = 2.0)
      : FFH92Tau(c_,f,norm_p), lambda(l) {};

   /// Evaluate the coefficient at @a ip.
   virtual real_t Eval(ElementTransformation &T,
                       const IntegrationPoint &ip) override;

   // Destructor
   ~FF91Delta() {};
};

} // namespace RBVMS

#endif
