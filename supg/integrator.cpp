// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
// -----------------------------------------------------------------------------
// SUPG integrator with r-switch (smooth-min) MATRIX-BASED stabilization parameter
// Steady-state convection–diffusion with discontinuity-capturing diffusion.
//
// τ is metric-based via G = J^{-T} J^{-1}:
//
//   tau_adv  = 1 / sqrt( a^T G a )
//   tau_diff = 1 / ( C_d * k * sqrt(G:G) )
//
//   tau = ( tau_adv^{-r} + tau_diff^{-r} )^{-1/r}    (r=2 typical)
//
// -----------------------------------------------------------------------------

#include "integrator.hpp"

#include <algorithm>
#include <cmath>

using namespace mfem;

// -----------------------------------------------------------------------------
// Integration rule based on FE order
// -----------------------------------------------------------------------------
const IntegrationRule &StabConvDifIntegrator::GetRule(
   const FiniteElement &trial_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans)
{
   int order = trial_fe.GetOrder() + test_fe.GetOrder();
   return IntRules.Get(trial_fe.GetGeomType(), order);
}

// -----------------------------------------------------------------------------
// r-switch MATRIX-BASED SUPG stabilization parameter (steady-state)
// -----------------------------------------------------------------------------
real_t StabConvDifIntegrator::GetTau(real_t &k, Vector &a, DenseMatrix &Gij)
{
   // Diffusion scaling constant (often same order as CI=1/12 used in baseline tau)
   const double Cd  = 1.0 / 12.0;

   // r-switch exponent (r=2 is a common smooth-min choice)
   const double r   = 2.0;

   // Numerical safety
   const double eps = 1e-14;

   const int dim = Gij.Width();

   // Compute:
   //   aGa   = a^T G a
   //   Gfro2 = G:G = sum_ij G_ij^2  (squared Frobenius norm)
   double aGa   = 0.0;
   double Gfro2 = 0.0;

   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         const double gij = Gij(i, j);
         aGa   += gij * a[i] * a[j];
         Gfro2 += gij * gij;
      }
   }

   aGa   = std::max(aGa, eps);
   Gfro2 = std::max(Gfro2, eps);

   // Ensure nonnegative diffusion coefficient
   const double kk = std::max((double)k, 0.0);

   // Characteristic timescales
   const double tau_adv  = 1.0 / std::sqrt(aGa);
   const double tau_diff = 1.0 / (Cd * kk * std::sqrt(Gfro2) + eps);

   // r-switch blend: tau = (tau_adv^{-r} + tau_diff^{-r})^{-1/r}
   const double inv_tau_adv_r  = std::pow(1.0 / std::max(tau_adv,  eps), r);
   const double inv_tau_diff_r = std::pow(1.0 / std::max(tau_diff, eps), r);

   double inv_tau_r = inv_tau_adv_r + inv_tau_diff_r;
   inv_tau_r = std::max(inv_tau_r, eps);

   return (real_t) std::pow(inv_tau_r, -1.0 / r);
}

// -----------------------------------------------------------------------------
// Discontinuity-capturing diffusion (unchanged from your baseline)
// -----------------------------------------------------------------------------
real_t StabConvDifIntegrator::GetKdc(Vector &a,
                                    real_t &res,
                                    Vector &dphidx,
                                    DenseMatrix &Gij)
{
   real_t h = 1.0 / sqrt(Gij.Trace() / Gij.Width());
   return kdc0 * h * a.Norml2()
        + kdc1 * h * fabs(res) / (dphidx.Norml2() + 1e-10);
}

// -----------------------------------------------------------------------------
// Set dimension of temporary variables
// -----------------------------------------------------------------------------
void StabConvDifIntegrator::SetDim(int d)
{
   if (dim != d)
   {
      dim = d;
      Gij.SetSize(dim);
      a.SetSize(dim);
      dphidx.SetSize(dim);
   }
}

// -----------------------------------------------------------------------------
// Assemble element nonlinear residual
// -----------------------------------------------------------------------------
void StabConvDifIntegrator::AssembleElementVector(const FiniteElement &el,
                                                  ElementTransformation &Trans,
                                                  const Vector &elfun,
                                                  Vector &elvect)
{
   real_t w, mu, f, phi, dphidx2, res;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd, dim);
   lshape.SetSize(nd);
   test.SetSize(nd);

   const IntegrationRule *ir =
      NonlinearFormIntegrator::IntRule ?
      NonlinearFormIntegrator::IntRule :
      &GetRule(el, el, Trans);

   elvect = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint(&ip);
      w = Trans.Weight() * ip.weight;

      // Metric tensor G = J^{-T} J^{-1}
      MultAtB(Trans.InverseJacobian(),
              Trans.InverseJacobian(),
              Gij);

      // Shape functions
      el.CalcPhysShape(Trans, shape);
      phi = shape * elfun;

      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      el.CalcPhysLinLaplacian(Trans, lshape);
      dphidx2 = lshape * elfun;

      // Coefficients
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f  = force_cf->Eval(Trans, ip);

      // Galerkin convection term
      res = a * dphidx - f;
      elvect.Add(w * res, shape);

      // Diffusion (Galerkin + artificial/DC)
      res = a * dphidx - mu * dphidx2 - f;
      dshape.Mult(dphidx, test);
      elvect.Add(w * (GetKdc(a, res, dphidx, Gij) + mu), test);

      // Stabilized test function: (a·∇v + type*mu*Δv)
      dshape.Mult(a, test);
      test.Add((int)type * mu, lshape);

      // SUPG stabilization term
      elvect.Add(w * GetTau(mu, a, Gij) * res, test);
   }
}

// -----------------------------------------------------------------------------
// Assemble element Jacobian (Newton consistency)
// -----------------------------------------------------------------------------
void StabConvDifIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                ElementTransformation &Trans,
                                                const Vector &elfun,
                                                DenseMatrix &elmat)
{
   real_t w, mu, f, phi, dphidx2, res;

   SetDim(el.GetDim());
   int nd = el.GetDof();

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd, dim);
   lshape.SetSize(nd);
   trail.SetSize(nd);
   test.SetSize(nd);

   const IntegrationRule *ir =
      NonlinearFormIntegrator::IntRule ?
      NonlinearFormIntegrator::IntRule :
      &GetRule(el, el, Trans);

   elmat = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint(&ip);
      w = Trans.Weight() * ip.weight;

      // Metric tensor G = J^{-T} J^{-1}
      MultAtB(Trans.InverseJacobian(),
              Trans.InverseJacobian(),
              Gij);

      // Shape functions
      el.CalcPhysShape(Trans, shape);
      phi = shape * elfun;

      el.CalcPhysDShape(Trans, dshape);
      dshape.MultTranspose(elfun, dphidx);

      el.CalcPhysLinLaplacian(Trans, lshape);
      dphidx2 = lshape * elfun;

      // Coefficients
      adv_cf->Eval(a, Trans, ip);
      mu = mu_cf->Eval(Trans, ip);
      f  = force_cf->Eval(Trans, ip);

      // Galerkin convection contribution
      dshape.Mult(a, trail);
      AddMult_a_VWt(w, shape, trail, elmat);

      // Diffusion (Galerkin + artificial/DC)
      res = a * dphidx - mu * dphidx2 - f;
      AddMult_a_AAt(w * (GetKdc(a, res, dphidx, Gij) + mu),
                    dshape, elmat);

      // Stabilized test function
      dshape.Mult(a, test);
      test.Add((int)type * mu, lshape);

      // Stabilized trial function
      dshape.Mult(a, trail);
      trail.Add(-mu, lshape);

      // SUPG Jacobian contribution
      AddMult_a_VWt(w * GetTau(mu, a, Gij), test, trail, elmat);
   }
}
