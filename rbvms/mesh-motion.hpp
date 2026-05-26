// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_MESH_MOTION_HPP
#define RBVMS_MESH_MOTION_HPP

#include "mfem.hpp"
#include "precice/precice.hpp"

using namespace std;
using namespace mfem;

namespace RBVMS
{


/// A coefficient that is scales with the mesh-size
class PowerDetCoefficient : public Coefficient
{
public:
   real_t power;

   /// c is value of constant function
   explicit PowerDetCoefficient(real_t p = 1.0) { power=p; }

   /// Evaluate the coefficient at @a ip.
   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      return std::pow(T.Weight(), -power);
   }
};

/** This class is a

.
*/
class MeshMotion
{
public:
   precice::Participant &participant;
   ParMesh &pmesh;

   int dim = pmesh.Dimension();
   std::vector<int>     vertexIDs;
   int vertexSize;
   Array<int> fsi_dofs;
   std::string meshName;

   GridFunction *nodes;
   ParFiniteElementSpace *pfes;

   std::vector<double> forces, disp, disp0;

   ParGridFunction pgf_d, pgf_d0, pgf_um;

   PowerDetCoefficient lambda_func;
   PowerDetCoefficient mu_func;

   ParBilinearForm *a;
   Vector nodes0;

   ParLinearForm *f;

public:
   MeshMotion(precice::Participant &part,
              ParMesh &pmesh,
              const char *meshName,
              Array<int> bdr_is_fsi);


   // ParGridFunction *GetMotion(){ return *pgf_d;};
   //ParGridFunction *GetVelocity(){ return *pgf_um;};

   ParGridFunction &GetMotion() { return pgf_d;};
   ParGridFunction &GetVelocity() { return pgf_um;};

   void Solve(double dt);

   void SetTimeLevel(double alpha);

   void SetVelocityBCs(Vector &x);

   void SetForce(Vector &x);
};






} // namespace RBVMS

#endif
