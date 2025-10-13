// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef SUPG_INTEGRATOR_HPP
#define SUPG_INTEGRATOR_HPP

#include "mfem.hpp"

using namespace mfem;

/** This Class defines an integrator for stabilized multi-dimensional
    convection-reaction equation.

     $(a \cdot \nabla u, v) + (k \nabla u,\nabla v)
    + \sum (a \cdot \nabla u + k u, \tau (a \cdot \nabla v - k \Delta v - f))_e - (f,  v)$
*/
class StabConvDifIntegrator : public NonlinearFormIntegrator
{
private:
   // Physical parameters
   VectorCoefficient *adv_cf;
   Coefficient *mu_cf;
   Coefficient *force_cf;

   /// Discontinuity capturing parameters
   real_t kdc0;  // inconsistent part
   real_t kdc1;  // consistent part

   /// The stabilization parameter
   void SetDim(int dim);

   /// The stabilization parameter
   real_t GetTau(real_t &k, Vector &a, DenseMatrix &Gij);

   /// The disconituity capturing parameter
   real_t GetKdc(Vector &a, real_t &res, Vector &dphidx, DenseMatrix &Gij);

   /// Temporary variables
   int dim;
   Vector a, dphidx, shape, trail, test;
   DenseMatrix dshape, Gij;

public:
   /// Constructor
   StabConvDifIntegrator(VectorCoefficient &a,
                         Coefficient &m,
                         Coefficient &f,
                         real_t k0 = 0.01,
                         real_t k1 = 0.25) : adv_cf(&a), mu_cf(&m), force_cf(&f)
   {
      kdc0 = k0;
      kdc1 = k1;
      dim = -1;
   };

   /// Destructor
   ~StabConvDifIntegrator() {};

   /// Set the penalty parameter for pinning the zero level-set
   void SetInconsistentDC(real_t k0) { kdc0 = k0; };

   /// Set the penalty parameter for pinning the zero level-set
   void SetConsistentDC(real_t k1) { kdc1 = k1; };

   /// Compute the element nonlinear residual
   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun,
                                      Vector &elvect) override;

   /// Compute the element jacobian
   virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun,
                                    DenseMatrix &elmat) override;

   /// Give element integration rule
   static const IntegrationRule &GetRule(const FiniteElement &trial_fe,
                                         const FiniteElement &test_fe,
                                         ElementTransformation &Trans);
};

#endif
