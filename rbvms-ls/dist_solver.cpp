// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "dist_solver.hpp"

using namespace mfem;

// Default values
real_t Heaviside::epsilon = 1e-10;
real_t Heaviside::eps = 2.0;

// Compute element size
real_t Heaviside::h(Vector &grad_phi, DenseMatrix &Gij)
{
   int dim = grad_phi.Size();
   real_t gGg = 0.0;
   for (int j = 0; j < dim; j++)
   {
      real_t gj = grad_phi[j];
      for (int i = 0; i < dim; i++)
      {
         gGg += Gij(i,j)*grad_phi[i]*gj;
      }
   }
   return grad_phi.Norml2()/sqrt(fmax(gGg, epsilon));
}

real_t Heaviside::h(Vector &grad_phi, ElementTransformation &Tr)
{
   // Metric tensor
   DenseMatrix Gij(grad_phi.Size());
   MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);

   return h(grad_phi, Gij);
}

// Compute relative distance wrt zero level
real_t Heaviside::rphi(real_t &phi, real_t &h)
{
   return phi/(eps*h);
}

// Smooth step function
real_t Heaviside::step(real_t &rphi)
{
   if (rphi < -1.0)
   {
      return 0.0;
   }
   else if (rphi > 1.0)
   {
      return 1.0;
   }
   else
   {
      return (1.0 + sin(M_PI*rphi/2))/2;
   }
}

real_t Heaviside::step(real_t &phi, Vector &grad_phi, DenseMatrix &Gij)
{
   real_t h = Heaviside::h(grad_phi, Gij);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::step(rphi);
}

real_t Heaviside::step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::step(rphi);
}

// Smooth sign function
real_t Heaviside::sign(real_t &rphi)
{
   return 2*Heaviside::step(rphi) - 1.0;
}

real_t Heaviside::sign(real_t &phi, Vector &grad_phi, DenseMatrix &Gij)
{
   return 2*Heaviside::step(phi, grad_phi, Gij) - 1.0;
}

real_t Heaviside::sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   return 2*Heaviside::step(phi, grad_phi, Tr) - 1.0;
}

// Smooth dirac distribution
real_t Heaviside::dirac(real_t &rphi, real_t &h)
{
   if (rphi < -1.0)
   {
      return 0.0;
   }
   else if (rphi > 1.0)
   {
      return 0.0;
   }
   else
   {
      return cos(M_PI*rphi/2)*M_PI/(4*eps*h);
   }
}

real_t Heaviside::dirac(real_t &phi, Vector &grad_phi,
                        DenseMatrix &Gij)
{
   real_t h = Heaviside::h(grad_phi, Gij);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::dirac(rphi, h);
}

real_t Heaviside::dirac(real_t &phi, Vector &grad_phi,
                        ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);
   return Heaviside::dirac(rphi, h);
}

// Define the integration rule based on FE order
const IntegrationRule &StabConvReactIntegrator::GetRule(
   const FiniteElement &trial_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans)
{
   int order = trial_fe.GetOrder() + test_fe.GetOrder();
   return IntRules.Get(trial_fe.GetGeomType(), order);
}

// Define convective tau
real_t StabConvReactIntegrator::GetTau(real_t &k, Vector &a, DenseMatrix &Gij)
{
   // Reaction part
   real_t  tau_i2 = k*k;

   // Convective part
   int dim = Gij.Width();
   for (int j = 0; j < dim; j++)
   {
      real_t aj = a[j];
      for (int i = 0; i < dim; i++)
      {
         tau_i2 += Gij(i,j)*a[i]*aj;
      }
   }

   // Momentum stabilisation parameter
   return 1.0/sqrt(fmax(tau_i2, 1e-20));
}

// Define convective kdc
real_t StabConvReactIntegrator::GetKdc(real_t &res,
                                       Vector &dphidx,
                                       DenseMatrix &Gij)
{
   real_t h = 1.0/sqrt(Gij.Trace()/dphidx.Size());
   return 0.01*h + 0.25*h*fabs(res)/fmax(dphidx.Norml2(), Heaviside::epsilon);
}

// Compute the element nonlinear residual
void StabConvReactIntegrator::AssembleElementVector(const FiniteElement &el,
                                                    ElementTransformation &Trans,
                                                    const Vector &elfun,
                                                    Vector &elvect)
{
   int nd = el.GetDof();
   int dim = el.GetDim();
   real_t w,k,tau,f,phi,phi0,res;
   Vector a(dim), dphidx(dim), dphidx0(dim);

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   test.SetSize(nd);
   Gij.SetSize(dim);

   const IntegrationRule *ir = NonlinearFormIntegrator::IntRule ?
                               NonlinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   elvect = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      MultAtB(Trans.InverseJacobian(),Trans.InverseJacobian(),Gij);

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      phi = shape*elfun;

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      // Interface params
      phi0 = ls_gf->GetValue(Trans, ip);
      ls_gf-> GetGradient(Trans, dphidx0);

      real_t h = Heaviside::h(dphidx0, Gij);
      real_t rphi = Heaviside::rphi(phi0,h);
      real_t Se = Heaviside::sign(rphi);
      real_t de = Heaviside::dirac(rphi,h);

      f = Se + lambda*de*phi0;
      k = lambda*de;
      a.Set(Se/fmax(dphidx.Norml2(),Heaviside::epsilon), dphidx);

      tau = GetTau(k, a, Gij);

      // Strong residual
      res = a*dphidx + k*phi - f;

      // Compute test function
      dshape.Mult(a, test); // Add Convection stabilization term
      test.Add(k, shape);   // Add Reaction stabilization term
      test *= tau;          // Scale stabilization term with parameter
      test += shape;        // Galerkin term

      // Weak residual
      elvect.Add(w*res, test);

      // Artificial diffusion
      //h = 1.0/fmax(Trans.InverseJacobian().FNorm(), 10e-10);
      //real_t kdc = 0.25*h*fabs(res)/fmax(dphidx.Norml2(), 10e-10);
      //kdc += 0.01*h;
      real_t kdc= GetKdc(res, dphidx, Gij);
      dshape.Mult(dphidx, test);
      elvect.Add(w*kdc, test);
   }
}

// Compute the element jacobian
void StabConvReactIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                  ElementTransformation &Trans,
                                                  const Vector &elfun,
                                                  DenseMatrix &elmat)
{
   int nd = el.GetDof();
   int dim = el.GetDim();
   real_t w,k,f,phi,phi0,res,tau;
   Vector a(dim), dphidx(dim), dphidx0(dim);

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   trail.SetSize(nd);
   test.SetSize(nd);
   Gij.SetSize(dim);

   const IntegrationRule *ir = NonlinearFormIntegrator::IntRule ?
                               NonlinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   elmat = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;
      MultAtB(Trans.InverseJacobian(),Trans.InverseJacobian(),Gij);

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      phi = shape*elfun;

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      // Interface params
      phi0 = ls_gf->GetValue(Trans, ip);
      ls_gf-> GetGradient(Trans, dphidx0);

      real_t h = Heaviside::h(dphidx0, Gij);
      real_t rphi = Heaviside::rphi(phi0,h);
      real_t Se = Heaviside::sign(rphi);
      real_t de = Heaviside::dirac(rphi,h);

      f = Se + lambda*de*phi0;
      k = lambda*de;
      a.Set(Se/fmax(dphidx.Norml2(),Heaviside::epsilon), dphidx);

      tau = GetTau(k, a, Gij);

      // Compute trail function
      dshape.Mult(a, trail); // Add Convection term
      trail.Add(k, shape);   // Add Reaction term

      // Compute test function
      test.Set(tau, trail);   // Add stabilization term
      test += shape;          // Add Galerkin term

      // Stablization term
      AddMult_a_VWt(w, test, trail, elmat);

      // Strong residual
      phi = shape*elfun;
      dshape.MultTranspose(elfun, dphidx);
      res = a*dphidx + k*phi - f;

      // Artificial diffusion
     // h = 1.0/fmax(Trans.InverseJacobian().FNorm(), Heaviside::epsilon);
     // real_t kdc = 0.25*h*fabs(res)/fmax(dphidx.Norml2(), Heaviside::epsilon);
      //kdc += 0.01*h;
      real_t kdc= GetKdc(res, dphidx, Gij);
      AddMult_a_AAt(w*kdc, dshape, elmat);
   }
}

ConvectionDistanceSolver::ConvectionDistanceSolver(ParFiniteElementSpace &space,
                                                   real_t lambda)
   : form(&space), gmres(space.GetComm()),newton_solver(space.GetComm())
{
   form.AddDomainIntegrator(&integrator);

   // Set up the Jacobian solver
   gmres.iterative_mode = false;
   gmres.SetPrintLevel(-1);

   // Default values
   gmres.SetRelTol(1e-4);
   gmres.SetAbsTol(1e-12);
   gmres.SetMaxIter(100);
   prec = new HypreSmoother();
   gmres.SetPreconditioner(*prec);

   // Set up the Newton solver
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(1);
   newton_solver.SetSolver(gmres);
   newton_solver.SetOperator(form);

   // Default values
   newton_solver.SetRelTol(1e-4);
   newton_solver.SetAbsTol(1e-12);
   newton_solver.SetMaxIter(10);

}

// Solver for the distance field given a zero level-set coefficient
void ConvectionDistanceSolver::ComputeScalarDistance(Coefficient
                                                     &zero_level_set,
                                                     ParGridFunction &distance)
{
   GridFunctionCoefficient* ls_gfcf = dynamic_cast<GridFunctionCoefficient*>
                                      (&zero_level_set);
   integrator.SetZeroLevelSet(ls_gfcf->GetGridFunction());

   distance.GetTrueDofs(sol);
   newton_solver.Mult(zero, sol);
   distance.Distribute(sol);
}
