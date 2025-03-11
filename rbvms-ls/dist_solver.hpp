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

/** This class manages the interface it implements the smooth Heaviside function
    and derived functions such as the smooth sign and smooth dirac functions.
*/
class Heaviside
{
public:
   static real_t epsilon;
   static real_t eps;

   static real_t h(Vector &grad_phi, ElementTransformation &Tr);
   static real_t rphi(real_t &phi, real_t &h);

   static real_t step(real_t &rphi);
   static real_t sign(real_t &rphi);
   static real_t dirac(real_t &rphi, real_t &h);

   static real_t step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
   static real_t sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
   static real_t dirac(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
};

/** This Class defines an integrator for stabilized multi-dimensional
    convection-reaction equation.

     $(a \cdot \nabla u, v) + (k u,v)
    + \sum (a \cdot \nabla u + k u, \tau (a \cdot \nabla v + k v))_e - (f,  v)$
*/
class StabConvReactIntegrator : public NonlinearFormIntegrator
{
protected:
   ///
   real_t lambda = 100.0;

   /// The function that defines the interface location
   GridFunction *ls_gf;

   /// The stabilization parameter
   real_t GetTau(real_t &k, Vector &a, ElementTransformation &T);

private:
   Vector shape, trail, test;
   DenseMatrix dshape;

public:
   StabConvReactIntegrator(){};

   ~StabConvReactIntegrator(){};

   void SetZeroLevelSet(ParGridFunction &zero_level_set)
   {
      ls_gf = &zero_level_set;
   }

   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun,
                                      Vector &elvect) override;

   virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun,
                                    DenseMatrix &elmat) override;

   static const IntegrationRule &GetRule(const FiniteElement &trial_fe,
                                         const FiniteElement &test_fe,
                                         ElementTransformation &Trans);
};


/**
*/
class ConvectionDistanceSolver : public common::DistanceSolver
{
private:
   ParNonlinearForm form;
   FGMRESSolver gmres;
   NewtonSolver newton_solver;
  // ConvectionCoefficient a_cf;
 //  ReactionCoefficient k_cf;
 //  ForceCoefficient f_cf;
   StabConvReactIntegrator integrator;
   Solver *prec;
   ParGridFunction phi0_gf, phi_gf;

public:
   ConvectionDistanceSolver(ParFiniteElementSpace &space,
                            real_t lambda);
   ~ConvectionDistanceSolver()
   {
      delete prec;
   }

   // Set linear solver parameters
   void SetLinearRelTol(real_t rtol) { gmres.SetRelTol(rtol); }
   void SetLinearAbsTol(real_t atol) { gmres.SetAbsTol(atol); }
   void SetLinearMaxIter(int maxiter) { gmres.SetMaxIter(maxiter); }

   // Set nonlinear solver parameters
   void SetNonlinearRelTol(real_t rtol) { newton_solver.SetRelTol(rtol); }
   void SetNonlinearAbsTol(real_t atol) { newton_solver.SetAbsTol(atol); }
   void SetNonlinearMaxIter(int maxiter) { newton_solver.SetMaxIter(maxiter); }

   void ComputeScalarDistance(Coefficient &zero_level_set,
                              ParGridFunction &distance);

};


#endif
