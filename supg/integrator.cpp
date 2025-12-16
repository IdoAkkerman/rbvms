// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "integrator.hpp"

using namespace mfem;

// -----------------------------------------------------------------------------
// Integration rule
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
// MATRIX–TENSOR SUPG STABILIZATION PARAMETER
//
//   tau^{-2} = a^T G a + C_I^2 * k^2 * (G : G)
//   G = J^{-T} J^{-1}
//
// This formulation:
//   - is dimension-independent (2D / 3D)
//   - captures mesh anisotropy naturally
//   - avoids ad-hoc scalar mesh sizes
// -----------------------------------------------------------------------------
real_t StabConvDifIntegrator::GetTau(real_t &k, Vector &a, DenseMatrix &Gij)
{
   const double CI = 1.0 / 12.0;
   const int dim = Gij.Width();

   double tau_conv = 0.0; // a^T G a
   double tau_diff = 0.0; // G : G

   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         const double gij = Gij(i, j);
         tau_conv += gij * a[i] * a[j];
         tau_diff += gij * gij;
      }
   }

   double denom = tau_conv + (CI * CI) * k * k * tau_diff;

   // Numerical safety
   if (denom < 1e-14) { denom = 1e-14; }

   return 1.0 / sqrt(denom);
}

// -----------------------------------------------------------------------------
// Discontinuity-capturing diffusion (unchanged)
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

      // Galerkin convection
      res = a * dphidx - f;
      elvect.Add(w * res, shape);

      // Diffusion (Galerkin + artificial)
      res = a * dphidx - mu * dphidx2 - f;
      dshape.Mult(dphidx, test);
      elvect.Add(w * (GetKdc(a, res, dphidx, Gij) + mu), test);

      // Stabilized test function
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

      // Metric tensor
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

      // Galerkin convection
      dshape.Mult(a, trail);
      AddMult_a_VWt(w, shape, trail, elmat);

      // Diffusion
      res = a * dphidx - mu * dphidx2 - f;
      AddMult_a_AAt(w * (GetKdc(a, res, dphidx, Gij) + mu),
                    dshape, elmat);

      // Stabilized test & trial
      dshape.Mult(a, test);
      test.Add((int)type * mu, lshape);

      dshape.Mult(a, trail);
      trail.Add(-mu, lshape);

      // SUPG Jacobian contribution
      AddMult_a_VWt(w * GetTau(mu, a, Gij),
                    test, trail, elmat);
   }
}
