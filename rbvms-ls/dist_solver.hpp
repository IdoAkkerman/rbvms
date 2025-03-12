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

/** This class manages the smooth interface functions,
    such as the smooth heaviside, smooth sign and smooth dirac functions.*/
class Heaviside
{
public:
   /// Small number to prevent division by zero
   static real_t epsilon;

   /// Half width of the interface in terms of elements
   static real_t eps;

   /// Compute element size
   static real_t h(Vector &grad_phi, DenseMatrix &Gij);
   static real_t h(Vector &grad_phi, ElementTransformation &Tr);

   /// Compute relative distance wrt zero level
   static real_t rphi(real_t &phi, real_t &h);

   /// Smooth step function, going from 0 to 1 in the interval [-1, 1]
   static real_t step(real_t &rphi);
   static real_t step(real_t &phi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);

   /// Smooth sign function, going from -1 to 1 in the interval [-1, 1]
   static real_t sign(real_t &rphi);
   static real_t sign(real_t &rphi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);

   /// Smooth dirac function, going from 0 to 0 in the interval [-1, 1]
   static real_t dirac(real_t &rphi, real_t &h);

   /// Integral of this function over phi is 1
   static real_t dirac(real_t &rphi, Vector &grad_phi, DenseMatrix &Gij);
   static real_t dirac(real_t &phi, Vector &grad_phi, ElementTransformation &Tr);
};

/** This Class defines an integrator for stabilized multi-dimensional
    convection-reaction equation.

     $(a \cdot \nabla u, v) + (k u,v)
    + \sum (a \cdot \nabla u + k u, \tau (a \cdot \nabla v + k v))_e - (f,  v)$
*/
class StabConvReactIntegrator : public NonlinearFormIntegrator
{
private:
   /// Penalty parameter for pinning the zero level-set
   real_t lambda;

   /// The function that defines the interface location
   const GridFunction *ls_gf;

   /// The stabilization parameter
   real_t GetTau(real_t &k, Vector &a, DenseMatrix &Gij);

   /// The disconituity capturing parameter
   real_t GetKdc(real_t &res, Vector &dphidx, DenseMatrix &Gij);

   /// Temporary variables
   Vector shape, trail, test;
   DenseMatrix dshape, Gij;

public:
   /// Constructor
   StabConvReactIntegrator(real_t l = 1.0){lambda = l;};

   /// Destructor
   ~StabConvReactIntegrator(){};

   /// Set the penalty parameter for pinning the zero level-set
   void SetPenalty(real_t l){lambda = l;};

   /// Provide the Gridfunction that specifies the zero level-set
   void SetZeroLevelSet(const GridFunction *zero_level_set)
   {
      ls_gf = zero_level_set;
   }

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


/** This Class defines a redistancing algorithm based on convection.
    The convective equation
     $S_{\eps} (\phi_0)\frac{\nabla \phi}{\|\nabla \phi\|} \cdot \nabla \phi
      + \lambda \delta_{eps}(\phi_0) (\phi-phi_0) 
      = S_{\eps} (\phi_0)$
    is solved using SUPG and Discontinuity capturing.
*/
class ConvectionDistanceSolver : public common::DistanceSolver
{
private:
   /// Abstract nonlinear formulation
   ParNonlinearForm form;

   /// Linear solver
   FGMRESSolver gmres;

   /// Nonlinear solver
   NewtonSolver newton_solver;

   /// Defines the element level integrals for the convection problem
   /// using SUPG and DC
   StabConvReactIntegrator integrator;

   /// Solver
   Solver *prec = nullptr;
   Vector zero, sol;

public:
   /// Constructor
   ConvectionDistanceSolver(ParFiniteElementSpace &space,
                            real_t lambda);

   /// Denstructor
   ~ConvectionDistanceSolver()
   {
       if (prec) delete prec;
   };

   // Set linear solver parameters
   void SetPenalty(real_t lambda){ integrator.SetPenalty(lambda); };

   // Set linear solver parameters
   void SetLinearRelTol(real_t rtol) { gmres.SetRelTol(rtol); }
   void SetLinearAbsTol(real_t atol) { gmres.SetAbsTol(atol); }
   void SetLinearMaxIter(int maxiter) { gmres.SetMaxIter(maxiter); }

   void SetLinearPreconditioner(Solver &pc)
   {
       if (prec) delete prec;
       prec = nullptr;
       gmres.SetPreconditioner(pc);
   }

   // Set nonlinear solver parameters
   void SetNonlinearRelTol(real_t rtol) { newton_solver.SetRelTol(rtol); }
   void SetNonlinearAbsTol(real_t atol) { newton_solver.SetAbsTol(atol); }
   void SetNonlinearMaxIter(int maxiter) { newton_solver.SetMaxIter(maxiter); }

   // Compute distance field for given level-set
   void ComputeScalarDistance(Coefficient &zero_level_set,
                              ParGridFunction &distance);
};


#endif
