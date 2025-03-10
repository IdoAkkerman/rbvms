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

#ifndef RBVMS_DIST_SOLVER_HPP
#define RBVMS_DIST_SOLVER_HPP

#include "mfem.hpp"
#include "miniapps/common/mfem-common.hpp"

using namespace mfem;
using namespace common;

/// This class manages the interface
class Heaviside
{
public:
   static real_t epsilon;
   static real_t eps;

   static real_t h(Vector &grad_phi, ElementTransformation &Tr);
   static real_t rphi(real_t &phi, real_t &h);

   static real_t step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
   static real_t sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
   static real_t dirac(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
};

/** This Class defines an integrator for stabilized multi-dimensional
    convection-reaction equation.

     $(a \cdot \nabla u, v) + (\nabla u, s v)
    + \sum (a \cdot \nabla u + s u, \tau (a \cdot \nabla v + s v))_e$

     $(f, \nabla v)
    + \sum (f, \tau (a \cdot \nabla v + s v))_e$
*/
class StabConvReactIntegrator : public BilinearFormIntegrator,
   public LinearFormIntegrator
{
protected:
   /// The advection field
   VectorCoefficient *adv;

   /// The reaction parameter and force fields
   Coefficient *react, *force;

   /// The stabilization parameter
   real_t GetTau(real_t &k, Vector &a, ElementTransformation &T);

private:
   Vector shape, adshape, trail, test;
   DenseMatrix dshape;

public:
   StabConvReactIntegrator(VectorCoefficient *a,
                           Coefficient *r,
                           Coefficient *f);

   ~StabConvReactIntegrator();

   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect);

   virtual void AssembleElementMatrix(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      DenseMatrix &elmat);

   static const IntegrationRule &GetRule(const FiniteElement &trial_fe,
                                         const FiniteElement &test_fe,
                                         ElementTransformation &Trans);


   using LinearFormIntegrator::AssembleRHSElementVect;
};










//
//
class ConvectionDistanceSolver : public common::DistanceSolver
{
public:
   ConvectionDistanceSolver() {}

   void ComputeScalarDistance(Coefficient &zero_level_set,
                              ParGridFunction &distance);

};


#endif
