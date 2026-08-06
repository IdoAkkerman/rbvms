
// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "mesh-motion.hpp"

using namespace mfem;
using namespace RBVMS;

MeshMotion::MeshMotion(precice::Participant &part,
                       ParMesh &pm,
                       const char *mName,
                       Array<int> bdr_is_fsi)
   : participant(part), pmesh(pm), meshName(mName)
{
   MFEM_VERIFY(participant.getMeshDimensions(meshName) == dim,
               "MFEM and Precice dimension don't match!");

   // Get nodes gridfunction
   // pmesh.EnsureNodes();
   nodes = pmesh.GetNodes();

   // Get associated parallel fespace
   pfes = dynamic_cast<ParFiniteElementSpace *>(nodes->FESpace());
   cout << "Mesh motion ordering" << pfes->GetOrdering() << endl;
   MFEM_VERIFY(pfes, "FESpace should be a parallel");

   // Get boundary dofs
   pfes->GetEssentialTrueDofs(bdr_is_fsi,fsi_dofs);
   // pfes->GetEssentialTrueDofs(bdr_is_fsi,fsi_dofs1);
   vertexSize = fsi_dofs.Size()/dim;

   // Get boundary coordinates
   std::vector<double>  vertices(vertexSize * dim);
   vertexIDs.resize(vertexSize);

   if (pfes->GetOrdering() == Ordering::byNODES)
   {
      for (int j = 0; j < dim; j++)
      {
         for (int i = 0; i < vertexSize; i++)
         {
            vertices.at(j + i*dim) = nodes->Elem(fsi_dofs[i + j*vertexSize]);
         }
      }
   }
   else if (pfes->GetOrdering() == Ordering::byVDIM)
   {
      for (int i = 0; i < vertexSize; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            vertices.at(j + i*dim) = nodes->Elem(fsi_dofs[i*dim + j]);
         }
      }
   }
   else
   {
      mfem_error("Used FESpace Ordering not handled by precice init.");
   }

   // Communicate boundary coordinates with precice
   participant.setMeshVertices(meshName, vertices, vertexIDs);

   // Check, resize and initialize data members to communicate
   MFEM_VERIFY(participant.getDataDimensions(meshName,"Force") == dim,
               "MFEM and Precice dimension don't match!");
   MFEM_VERIFY(participant.getDataDimensions(meshName,"Displacement") == dim,
               "MFEM and Precice dimension don't match!");

   forces.resize(vertexSize*dim);
   disp.resize(vertexSize*dim);
   disp0.resize(vertexSize*dim);

   for (int i = 0; i < disp0.size(); i++)
   {
      disp0[i] = 0.0;
   }

   // Setup mesh motion problem

   pgf_d.SetSpace(pfes); pgf_d = 0.0;
   //   pgf_d0.SetSpace(pfes); // tbd
   pgf_um.SetSpace(pfes); pgf_um = 0.0;

   Array<int> mm_bdr_dof;
   pfes->GetBoundaryTrueDofs(mm_bdr_dof);

   a = new ParBilinearForm(pfes);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

   f = new ParLinearForm(pfes);
   //  f->AddBdrFaceIntegrator(new EvalTraction(), bdr_is_fsi);
}

void MeshMotion::Solve(real_t dt)
{
   // Get boundary displacement
   participant.readData(meshName,
                        "Displacement",
                        vertexIDs,
                        dt,
                        disp);

   // Set boundary condition for incremental displacement
   pgf_d = 0.0;
   if (pfes->GetOrdering() == Ordering::byNODES)
   {
      for (int j = 0; j < dim; j++)
      {
         for (int i = 0; i < vertexSize; i++)
         {
            pgf_d[fsi_dofs[i + j*vertexSize]] = disp[j + i*dim] - disp0[j + i*dim] ;

         }
      }
   }
   else if (pfes->GetOrdering() == Ordering::byVDIM)
   {
      for (int i = 0; i < vertexSize; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            pgf_d[fsi_dofs[i*dim + j]] = disp[j + i*dim] - disp0[j + i*dim] ;
         }
      }
   }

   // Compute linear system
   a->Update();
   a->Assemble();

   Array<int> bdr_dof;
   pfes->GetBoundaryTrueDofs(bdr_dof);

   Vector b(pgf_d.Size()); b = 0.0;
   HypreParMatrix A;
   Vector B, X;
   a->FormLinearSystem(bdr_dof, pgf_d, b, A, X, B);

   // Solve linear system
   HypreBoomerAMG *amg = new HypreBoomerAMG(A);
   amg->SetPrintLevel(0);
   HyprePCG *pcg = new HyprePCG(A);
   pcg->SetTol(1e-3);
   pcg->SetMaxIter(500);
   pcg->SetPrintLevel(3);
   pcg->SetPreconditioner(*amg);

   pcg->Mult(B, X);
   a->RecoverFEMSolution(X, b, pgf_d);

   // Extract mesh velocity
   pgf_um = pgf_d;
   pgf_um /= dt;

   // Set the new mesh
   nodes0 = (*nodes);
   SetTimeLevel(1.0);

   // Remember the boundary displacement for incremental solves
   disp0 = disp;
}

// ## LOR
// =============================================================
// void MeshMotion::TransferForcesToHO(ParGridFunction &forces_ho)
// {
   
//    if (pfes->GetOrdering() == Ordering::byVDIM)
//    {
//       cout << "Transfering forces to forces_ho byVDIM" << endl;
//       for (int j = 0; j < dim; j++)
//       {
//          for (int i = 0; i < vertexSize; i++)
//          {
//             forces_ho[fsi_dofs[i + j*vertexSize]] = forces[j + i*dim];

//          }
//       }
//    }
//    else if (pfes->GetOrdering() == Ordering::byNODES)
//    {
//       cout << "Transfering forces to forces_ho byNODES" << endl;
//       for (int i = 0; i < vertexSize; i++)
//       {
//          for (int j = 0; j < dim; j++)
//          {
//             forces_ho[fsi_dofs[i*dim + j]] = forces[j + i*dim];
//          }
//       }
//    }
// }

// void MeshMotion::TransferForcesFromLOR(ParGridFunction &forces_lor)
// {
//    if (pfes->GetOrdering() == Ordering::byVDIM)
//    {
//       cout << "Transfering forces from LOR byVDIM" << endl;
//       for (int j = 0; j < dim; j++)
//       {
//          for (int i = 0; i < vertexSize; i++)
//          {
//             forces[j + i*dim] = forces_lor[fsi_dofs[i + j*vertexSize]];

//          }
//       }
//    }
//    else if (pfes->GetOrdering() == Ordering::byNODES)
//    {
//       cout << "Transfering forces from LOR byNODES" << endl;
//       for (int i = 0; i < vertexSize; i++)
//       {
//          for (int j = 0; j < dim; j++)
//          {
//             forces[j + i*dim] = forces_lor[fsi_dofs[i*dim + j]];
//          }
//       }
//    }
// }
// =============================================================

void MeshMotion::SetTimeLevel(double alpha)
{
   add(nodes0, alpha, pgf_d, *nodes);
}

void MeshMotion::SetVelocityBCs(Vector &x)
{
   for (int i = 0; i < fsi_dofs.Size(); i++)
   {
      x[fsi_dofs[i]] = pgf_um[fsi_dofs[i]];
   }
}

void MeshMotion::SetForce(Vector &f)
{
   if (pfes->GetOrdering() == Ordering::byNODES)
   {
      for (int j = 0; j < dim; j++)
      {
         for (int i = 0; i < vertexSize; i++)
         {
            forces[j + i*dim] = f[fsi_dofs[i + j*vertexSize]];
         }
      }
   }
   else if (pfes->GetOrdering() == Ordering::byVDIM)
   {
      for (int i = 0; i < vertexSize; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            forces[j + i*dim] = f[fsi_dofs[i*dim + j]];
         }
      }
   }
}

ForceExtraction::ForceExtraction (Array<ParFiniteElementSpace *> spaces_,
                                  Array<int> bdr_is_fsi_,
                                  Coefficient &mu_)
   :
   ParGridFunction(spaces_[0]), spaces(spaces_),
   bdr_is_fsi(bdr_is_fsi_),c_mu(mu_)
{

}

void ForceExtraction::ComputeBoundaryForce(BlockVector &xp)
{


   ParMesh &pmesh = *pfes->GetParMesh();
   int dim  = pmesh.Dimension();

   Vector &tmp = (*this);
   tmp = 0.0;

   // Compute force
   //(*this) = 0.0;
   //this->operator=(0.0);
   Array<int> vdofs0, vdofs1;
   Vector el_u, el_p;
   Vector nor(dim);
   for (int i = 0; i < pmesh.GetNBE(); i++)
   {
      const int bdr_attr = pmesh.GetBdrAttribute(i);
      // if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }
      //std::cout<<bdr_attr<<"  --> "<<bdr_is_fsi[bdr_attr]<<std::endl;
      if (bdr_is_fsi[bdr_attr-1] == 0) { continue; }

      // ElementTransformation *T;
      FaceElementTransformations *Tr = pmesh.GetBdrFaceTransformations(i);
      if (Tr == NULL) { continue; }
      const FiniteElement &el0 = *spaces[0] ->GetFE(Tr->Elem1No);
      const FiniteElement &el1 = *spaces[1] ->GetFE(Tr->Elem1No);
      spaces[0] -> GetElementVDofs (Tr -> Elem1No, vdofs0);
      spaces[1] -> GetElementVDofs (Tr -> Elem1No, vdofs1);

      xp.GetBlock(0).GetSubVector(vdofs0, el_u);
      xp.GetBlock(1).GetSubVector(vdofs1, el_p);
      const IntegrationRule *ir = &IntRules.Get(Tr->GetGeometryType(),
                                                2*el0.GetOrder());

      int dof_u = el0.GetDof();
      int dof_p = el1.GetDof();

      Vector elemvect(dof_u*dim); elemvect = 0.0;

      /// Solution & Residual vector
      DenseMatrix elf_u, elv_u;
      elf_u.UseExternalData(el_u.GetData(), dof_u, dim);
      elv_u.UseExternalData(elemvect.GetData(), dof_u, dim);

      /// Shape function data
      Vector  sh_u(dof_u),sh_p(dof_p), traction(dim);
      DenseMatrix shg_u(dof_u, dim), grad_u(dim,dim);


      for (int p = 0; p < ir->GetNPoints(); p++)
      {
         const IntegrationPoint &ip = ir->IntPoint(p);
         Tr->SetAllIntPoints(&ip);

         // Access the neighboring element's integration point
         const IntegrationPoint &eip = Tr->GetElement1IntPoint();

         real_t mu_val = c_mu.Eval(*Tr->Elem1, eip);
         real_t w = ip.weight * Tr->Weight();
         CalcOrtho(Tr->Jacobian(), nor);
         nor /= nor.Norml2();
         el0.CalcPhysShape(*Tr->Elem1, sh_u);
         el0.CalcPhysDShape(*Tr->Elem1, shg_u);
         MultAtB(elf_u, shg_u, grad_u);
         grad_u.Symmetrize();           // Grad to strain

         el1.CalcPhysShape(*Tr->Elem1, sh_p);
         real_t pressure = sh_p*el_p;

         // Traction
         grad_u.Mult(nor, traction);
         traction *= -2*mu_val;                   // Consistency
         traction.Add(pressure, nor);             // Pressure
         AddMult_a_VWt(w, sh_u, traction, elv_u); //
      }
      AddElementVector (vdofs0, elemvect);
   }
}
