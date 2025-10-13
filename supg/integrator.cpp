// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "integrator.hpp"

using namespace mfem;

// Define the integration rule based on FE order
const IntegrationRule &StabConvDifIntegrator::GetRule(
   const FiniteElement &trial_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans)
{
   int order = trial_fe.GetOrder() + test_fe.GetOrder();
   return IntRules.Get(trial_fe.GetGeomType(), order);
}

// Define convective tau
real_t StabConvDifIntegrator::GetTau(real_t &k, Vector &a, DenseMatrix &Gij)
{
   double CI = 1.0/12.0;
   double tau = 0.0;
   int dim = Gij.Width();
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau += Gij(i,j)*a[i]*a[j] + CI*k*k*Gij(i,j)*Gij(i,j);
      }
   }
   return 1.0/sqrt(tau);
}

// Define convective kdc
real_t StabConvDifIntegrator::GetKdc(Vector &a,
                                       real_t &res,
                                       Vector &dphidx,
                                       DenseMatrix &Gij)
{
   real_t h = 1.0/sqrt(Gij.Trace()/Gij.Width());
   return kdc0*h*a.Norml2() + kdc1*h*fabs(res)/dphidx.Norml2();
}

// Set the dimension of temporary variables
void StabConvDifIntegrator::SetDim(int d)
{
   if (dim != d )
   {
      dim = d;
      Gij.SetSize(dim);
      a.SetSize(dim);
      dphidx.SetSize(dim);
   }
}

// Compute the element nonlinear residual
void StabConvDifIntegrator::AssembleElementVector(const FiniteElement &el,
                                                    ElementTransformation &Trans,
                                                    const Vector &elfun,
                                                    Vector &elvect)
{
   real_t w,mu,tau,f,phi,res,kdc;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   test.SetSize(nd);

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

      // Force, reaction and convection parameter
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f =  force_cf->Eval(Trans, ip);

      // Strong residual
      res = a*dphidx + mu*phi - f;
      res = phi - f;
     // res *= -1.0;
      // Numerical parameters
      tau = GetTau(mu, a, Gij);
      kdc = GetKdc(a, res, dphidx, Gij);

      // Compute test function
      dshape.Mult(a, test);  // Add Convection stabilization term
      test.Add(mu, shape);   // Add Reaction stabilization term
      test *= tau;           // Scale Stabilization term with parameter
      test *= 0.0;
      test += shape;         // Galerkin term

      // Weak residual
      elvect.Add(w*res, test);

      // Artificial diffusion
      dshape.Mult(dphidx, test);
      //elvect.Add(w*kdc, test);
   }
 // elvect.Print(std::cout, 5555);
}

// Compute the element jacobian
void StabConvDifIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                  ElementTransformation &Trans,
                                                  const Vector &elfun,
                                                  DenseMatrix &elmat)
{
   real_t w,mu,tau,f,phi,res,kdc;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   trail.SetSize(nd);
   test.SetSize(nd);

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

      // Force, reaction and convection parameter
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f =  force_cf->Eval(Trans, ip);

      // Strong residual
      res = a*dphidx + mu*phi - f;

      // Numerical parameters
      tau = GetTau(mu, a, Gij);
      kdc = GetKdc(a, res, dphidx, Gij);

      // Compute trail function
      dshape.Mult(a, trail);  // Add Convection term
      trail.Add(mu, shape);   // Add Reaction term

      // Compute test function
      test.Set(tau, trail);   // Add Stabilization term
      test += shape;          // Add Galerkin term

      // Stablization term
  //    AddMult_a_VWt(w, test, trail, elmat);
      AddMult_a_VVt(w, shape, elmat);

      // Artificial diffusion
     // AddMult_a_AAt(w*kdc, dshape, elmat);
   }
}
