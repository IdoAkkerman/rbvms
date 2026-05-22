
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

   //b = new ParLinearForm(pfes);
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
   //b->Update();
   a->Update();

   //b->Assemble();
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
   pcg->SetTol(1e-8);
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



//358
/*
  // Configure precice
  precice::Participant precice(std::string(precice_solverName),
                               std::string(precice_configFile), 0, 1); //??,rank,size);
  const std::string meshName(precice_meshName);
  std::vector<int>     vertexIDs;
  int vertexSize;
  Array<int> fsi_dofs;
  //{
  MFEM_VERIFY(precice.getMeshDimensions(meshName) == dim,
              "MFEM and Precice dimension don't match!");

  // Get nodes gridfunction
  pmesh.EnsureNodes();
  GridFunction *nodes = pmesh.GetNodes();

  // Get boundary dofs
  ParFiniteElementSpace *pfes = dynamic_cast<ParFiniteElementSpace *>
                                (nodes->FESpace());

  MFEM_VERIFY(pfes, "FESpace should be a parallel");
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

  if (vertexSize == 4)
  {
     std::cout<<"----------------------------------------------\n";
     for (int i = 0; i < vertexSize; i++)
     {
        std::cout<<vertices[dim * i]<<" "<<vertices[dim * i + 1]<<std::endl;
     }
  }
     // Set boundary coordinates
  precice.setMeshVertices(meshName, vertices, vertexIDs);

  // Set vectors
//  const int forceDim = precice.getDataDimensions(meshName,"Force");
//  const int dispDim = precice.getDataDimensions(meshName,"Displacement");

  MFEM_VERIFY(precice.getDataDimensions(meshName,"Force") == dim,
              "MFEM and Precice dimension don't match!");
  MFEM_VERIFY(precice.getDataDimensions(meshName,"Displacement") == dim,
              "MFEM and Precice dimension don't match!");

  std::vector<double> forces(vertexSize*dim);
  std::vector<double> disp(vertexSize*dim);
  std::vector<double> disp0(vertexSize*dim);

  for (int i = 0; i < disp0.size(); i++)
  {
     disp0[i] = 0.0;
  }
  mfem::out<<"forceDim = "<<forces.size()<<std::endl;
  mfem::out<<"dispDim  = "<<disp .size()<<std::endl;

  // Add fsi boundary to weak or strong
  if (fsi_strong)
  {
     strong_bdr.Append(fsi_bdr);
  }
  else
  {
     weak_bdr.Append(fsi_bdr);
  }
*/
//-------- 453



//-------- 533
/*
   ParGridFunction x_d(pfes); x_d = 0.0;
   ParGridFunction x_d0(pfes); x_d0 = 0.0;
   ParGridFunction x_um(pfes);
   form.SetMeshVelocity(&x_um);

   Array<int> mm_bdr_dof;
   pfes->GetBoundaryTrueDofs(mm_bdr_dof);

   PowerDetCoefficient lambda_func(1.0);
   PowerDetCoefficient mu_func(1.0);

   ParBilinearForm *a_mm = new ParBilinearForm(pfes);
   a_mm->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

   ParLinearForm *b_mm = new ParLinearForm(pfes);
*/
// 555


//   HypreParMatrix A_mm; // 665
//   Vector B_mm, X_mm;  // 666



// 696
/*
      precice.readData(meshName,
                      "Displacement",
                      vertexIDs,
                      dt_used,
                      disp);

     // Set FSI boundary displacement
     x_d = 0.0;
     if (pfes->GetOrdering() == Ordering::byNODES)
     {
        for (int j = 0; j < dim; j++)
        {
           for (int i = 0; i < vertexSize; i++)
           {
              x_d[fsi_dofs[i + j*vertexSize]] = disp[j + i*dim] - disp0[j + i*dim] ;

           }
        }
     }
     else if (pfes->GetOrdering() == Ordering::byVDIM)
     {
        for (int i = 0; i < vertexSize; i++)
        {
           for (int j = 0; j < dim; j++)
           {
              x_d[fsi_dofs[i*dim + j]] = disp[j + i*dim] - disp0[j + i*dim] ;
           }
        }
     }

     b_mm->Update();
     a_mm->Update();

     b_mm->Assemble();
     a_mm->Assemble();
     a_mm->FormLinearSystem(mm_bdr_dof, x_d, *b_mm, A_mm, X_mm, B_mm);
     HypreBoomerAMG *amg = new HypreBoomerAMG(A_mm);
     amg->SetPrintLevel(0);
     HyprePCG *pcg = new HyprePCG(A_mm);
     pcg->SetTol(1e-8);
     pcg->SetMaxIter(500);
     pcg->SetPrintLevel(3);
     pcg->SetPreconditioner(*amg);

     pcg->Mult(B_mm, X_mm);
     a_mm->RecoverFEMSolution(X_mm, *b_mm, x_d);


     disp0=disp;
     (*nodes) += x_d;

     // Extract mesh velocity
     x_um = x_d;
     x_um /= dt_used;
*/
//761
