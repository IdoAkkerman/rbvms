// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "dist_solver.hpp"

using namespace mfem;


real_t Heaviside::epsilon = 1e-10;
real_t Heaviside::eps = 4.0;
//real_t Heaviside::rho0 = 1.0;
//real_t Heaviside::rho1 = 1000.0;

real_t Heaviside::h(Vector &grad_phi, ElementTransformation &Tr)
{
   // Metric tensor
   DenseMatrix Gij;
   MultAtB(Tr.InverseJacobian(),Tr.InverseJacobian(),Gij);

   ///
   int dim = Gij.Width();
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

real_t Heaviside::rphi(real_t &phi, real_t &h)
{
   return phi/(eps*h);
}

//
real_t Heaviside::step(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);

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

//
real_t Heaviside::sign(real_t &phi, Vector &grad_phi, ElementTransformation &Tr)
{
   return 2*Heaviside::step(phi, grad_phi, Tr) - 1.0;
}

//
real_t Heaviside::dirac(real_t &phi, Vector &grad_phi,
                        ElementTransformation &Tr)
{
   real_t h = Heaviside::h(grad_phi, Tr);
   real_t rphi = Heaviside::rphi(phi, h);

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

StabConvReactIntegrator::StabConvReactIntegrator(VectorCoefficient *a,
                                                 Coefficient *k,
                                                 Coefficient *f)
   : adv(a), react(k), force(f)
{
}

StabConvReactIntegrator::~StabConvReactIntegrator()
{
}

const IntegrationRule &StabConvReactIntegrator::GetRule(
   const FiniteElement &trial_fe,
   const FiniteElement &test_fe,
   ElementTransformation &Trans)
{
   int order = trial_fe.GetOrder() + test_fe.GetOrder();
   return IntRules.Get(trial_fe.GetGeomType(), order);
}

real_t StabConvReactIntegrator::GetTau(real_t &k, Vector &a,
                                 ElementTransformation &T)
{
   real_t Cd = 6.0;
   real_t Ct = 1.0;

   // Metric tensor
   DenseMatrix Gij;
   MultAtB(T.InverseJacobian(),T.InverseJacobian(),Gij);

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


void StabConvReactIntegrator::AssembleRHSElementVect(const FiniteElement &el,
                                                     ElementTransformation &Trans,
                                                     Vector &elvect)
{
   int nd = el.GetDof();
   int dim = el.GetDim();
   real_t w,k,tau,f;
   Vector a(dim);

   elvect.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   adshape.SetSize(nd);
   test.SetSize(nd);

   const IntegrationRule *ir = LinearFormIntegrator::IntRule ?
                               LinearFormIntegrator::IntRule : &GetRule(el, el, Trans);

   elvect = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint (&ip);
      w = Trans.Weight() * ip.weight;

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);

      // Evaluate coefficients
      f = force->Eval(Trans, ip);

      // Galerkin term
      elvect.Add(w*f, shape);

      // Calculate shapes
      el.CalcPhysDShape(Trans, dshape);

      // Evaluate coefficients
      k = react->Eval(Trans, ip);
      adv->Eval(a, Trans, ip);
      tau = GetTau(k, a, Trans);

      // Advective derivative
      dshape.Mult(a, adshape);

      // Stablization term
      elvect.Add(w*f*tau, test);
   }
}

void StabConvReactIntegrator::AssembleElementMatrix(const FiniteElement &el,
                                                    ElementTransformation &Trans,
                                                    DenseMatrix &elmat)
{
   int nd = el.GetDof();
   int dim = el.GetDim();
   real_t w,k,tau;
   Vector a(dim);

   elmat.SetSize(nd);
   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   adshape.SetSize(nd);
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

      // Calculate shapes
      el.CalcPhysShape(Trans, shape);
      el.CalcPhysDShape(Trans, dshape);

//      real_t phi = sh_p*(*elsol[1]);


      // Evaluate coefficients
      k = react->Eval(Trans, ip);
      adv->Eval(a, Trans, ip);

      // Galerkin convection term
      dshape.Mult(a, adshape);
      AddMult_a_VWt(w, shape, adshape, elmat);

      // Galerkin diffusion term
      AddMult_a_AAt(w*k, dshape, elmat);

      // Evaluate coefficients
      tau = GetTau(k, a, Trans);

      // Stablization term
      AddMult_a_VWt(w*tau, test, trail, elmat);
   }
}




//
void ConvectionDistanceSolver::ComputeScalarDistance(Coefficient
                                                     &zero_level_set,
                                                     ParGridFunction &distance)
{
   ParFiniteElementSpace &pfes = *distance.ParFESpace();

   auto check_h1 = dynamic_cast<const H1_FECollection *>(pfes.FEColl());
   MFEM_VERIFY(check_h1 && pfes.GetVDim() == 1,
               "This solver supports only scalar H1 spaces.");

   // Compute average mesh size (assumes similar cells).
   ParMesh &pmesh = *pfes.GetParMesh();

   // Step 0 - transform the input level set into a source-type bump.
   ParGridFunction source(&pfes);
   source.ProjectCoefficient(zero_level_set);

   /*
      int amg_print_level = 0;

      // Solver.
      CGSolver cg(MPI_COMM_WORLD);
      cg.SetRelTol(1e-12);
      cg.SetMaxIter(100);
      cg.SetPrintLevel(print_level);
      OperatorPtr A;
      Vector B, X;

      // Step 1 - diffuse.
      ParGridFunction diffused_source(&pfes);
      for (int i = 0; i < diffuse_iter; i++)
      {
         // Set up RHS.
         ParLinearForm b(&pfes);
         GridFunctionCoefficient src_coeff(&source);
         b.AddDomainIntegrator(new DomainLFIntegrator(src_coeff));
         b.Assemble();

         // Diffusion and mass terms in the LHS.
         ParBilinearForm a_d(&pfes);
         a_d.AddDomainIntegrator(new MassIntegrator);
         ConstantCoefficient t_coeff(parameter_t);
         a_d.AddDomainIntegrator(new DiffusionIntegrator(t_coeff));
         a_d.Assemble();

         // Solve with Dirichlet BC.
         Array<int> ess_tdof_list;
         if (pmesh.bdr_attributes.Size())
         {
            Array<int> ess_bdr(pmesh.bdr_attributes.Max());
            ess_bdr = 1;
            pfes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
         }
         ParGridFunction u_dirichlet(&pfes);
         u_dirichlet = 0.0;
         a_d.FormLinearSystem(ess_tdof_list, u_dirichlet, b, A, X, B);
         auto *prec = new HypreBoomerAMG;
         prec->SetPrintLevel(amg_print_level);
         cg.SetPreconditioner(*prec);
         cg.SetOperator(*A);
         cg.Mult(B, X);
         a_d.RecoverFEMSolution(X, b, u_dirichlet);
         delete prec;

         // Diffusion and mass terms in the LHS.
         ParBilinearForm a_n(&pfes);
         a_n.AddDomainIntegrator(new MassIntegrator);
         a_n.AddDomainIntegrator(new DiffusionIntegrator(t_coeff));
         a_n.Assemble();

         // Solve with Neumann BC.
         ParGridFunction u_neumann(&pfes);
         ess_tdof_list.DeleteAll();
         a_n.FormLinearSystem(ess_tdof_list, u_neumann, b, A, X, B);
         auto *prec2 = new HypreBoomerAMG;
         prec2->SetPrintLevel(amg_print_level);
         cg.SetPreconditioner(*prec2);
         cg.SetOperator(*A);
         cg.Mult(B, X);
         a_n.RecoverFEMSolution(X, b, u_neumann);
         delete prec2;

         for (int ii = 0; ii < diffused_source.Size(); ii++)
         {
            // This assumes that the magnitudes of the two solutions are somewhat
            // similar; otherwise one of the solutions would dominate and the BC
            // won't look correct. To avoid this, it's good to have the source
            // away from the boundary (i.e. have more resolution).
            diffused_source(ii) = 0.5 * (u_neumann(ii) + u_dirichlet(ii));
         }
         source = diffused_source;
      }

      // Step 2 - solve for the distance using the normalized gradient.
      {
         // RHS - normalized gradient.
         ParLinearForm b2(&pfes);
         NormalizedGradCoefficient grad_u(diffused_source, pmesh.Dimension());
         b2.AddDomainIntegrator(new DomainLFGradIntegrator(grad_u));
         b2.Assemble();

         // LHS - diffusion.
         ParBilinearForm a2(&pfes);
         a2.AddDomainIntegrator(new DiffusionIntegrator);
         a2.Assemble();

         // No BC.
         Array<int> no_ess_tdofs;

         a2.FormLinearSystem(no_ess_tdofs, distance, b2, A, X, B);

         auto *prec = new HypreBoomerAMG;
         prec->SetPrintLevel(amg_print_level);
         OrthoSolver ortho(pfes.GetComm());
         ortho.SetSolver(*prec);
         cg.SetPreconditioner(ortho);
         cg.SetOperator(*A);
         cg.Mult(B, X);
         a2.RecoverFEMSolution(X, b2, distance);
         delete prec;
      }

      // Shift the distance values to have minimum at zero.
      // Shift the distance to conserve global volume
      real_t d_min_loc = 0.0;// distance.Min();
      real_t d_min_glob;
      MPI_Allreduce(&d_min_loc, &d_min_glob, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_MIN, pfes.GetComm());
      distance -= d_min_glob;
   */
}

