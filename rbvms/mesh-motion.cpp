
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
   pmesh.EnsureNodes();
   nodes = pmesh.GetNodes();

   // Get associated parallel fespace
   pfes = dynamic_cast<ParFiniteElementSpace *>(nodes->FESpace());
   MFEM_VERIFY(pfes, "FESpace should be a parallel");

   // Get boundary dofs
   pfes->GetEssentialTrueDofs(bdr_is_fsi,fsi_dofs);
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

void MeshMotion::SetTimeLevel(double alpha)
{
   add(nodes0, -alpha, pgf_d, *nodes);
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
