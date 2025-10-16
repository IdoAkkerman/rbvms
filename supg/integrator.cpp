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
   double CI = 1./12.0;
   double tau = 1e-10;
   int dim = Gij.Width();
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau += Gij(i,j)*a[i]*a[j] + CI*CI*k*k*Gij(i,j)*Gij(i,j);
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
   return kdc0*h*a.Norml2() + kdc1*h*fabs(res)/(dphidx.Norml2() + 1e-10);
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
   real_t w,mu,f,phi,dphidx2,res,res_red;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   lshape.SetSize(nd);
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

      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      el.CalcPhysLinLaplacian(Trans, lshape);
      dphidx2 = lshape*elfun;

      // Force, reaction and convection parameter
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f =  force_cf->Eval(Trans, ip);

      // Galerkin terms
      res = a*dphidx - f;
      elvect.Add(w*res, shape);

      // Diffusion (Galerkin + artificial)
      res = a*dphidx - mu*dphidx2 - f;
      dshape.Mult(dphidx, test);
      elvect.Add(w*(GetKdc(a, res, dphidx, Gij)+mu), test);

      // Compute stabilized test function
      dshape.Mult(a, test);             // Add Convection stabilization term
      test.Add((int) type*mu, lshape);  // Add Reaction stabilization term

      // Stabilized terms
      elvect.Add(w*GetTau(mu, a, Gij)*res, test);
   }
}

// Compute the element jacobian
void StabConvDifIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                  ElementTransformation &Trans,
                                                  const Vector &elfun,
                                                  DenseMatrix &elmat)
{
   real_t w,mu,f,phi,dphidx2,res;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   lshape.SetSize(nd);
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

      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      el.CalcPhysLinLaplacian(Trans, lshape);
      dphidx2 = lshape*elfun;

      // Force, reaction and convection parameter
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f =  force_cf->Eval(Trans, ip);

      // Convection Galerkin terms
      dshape.Mult(a, trail);
      AddMult_a_VWt(w, shape, trail, elmat);

      // Diffusion (Galerkin + artificial)
      res = a*dphidx - mu*dphidx2 - f;
      AddMult_a_AAt(w*(GetKdc(a, res, dphidx, Gij)+mu), dshape, elmat);

      // Compute stabilized test functions
      dshape.Mult(a, test);            // Add Convection term
      test.Add((int) type*mu, lshape); // Add Diffusion term

      // Compute stabilized trail functions
      dshape.Mult(a, trail);          // Add Convection term
      trail.Add(-mu, lshape);         // Add Diffusion term

      AddMult_a_VWt(w*GetTau(mu, a, Gij), test, trail, elmat);
   }
}
