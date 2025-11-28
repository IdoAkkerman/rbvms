// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
//   _____  ______      ____  __  _____
//   |  __ \|  _ \ \    / /  \/  |/ ____|
//   | |__) | |_) \ \  / /| \  / | (___
//   |  _  /|  _ < \ \/ / | |\/| |\___ \
//   | | \ \| |_) | \  /  | |  | |____) |
//   |_|  \_\____/   \/   |_|  |_|_____/
//
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------
#include <sys/stat.h>
#include "mfem.hpp"
#include "../util/coefficients.hpp"
#include "../util/solver.hpp"

#include "integrator.hpp"

using namespace std;
using namespace mfem;

extern void printInfo();
extern void line(int len);

// Routine for checking duplicity of boundary conditions
void CheckBoundaries(Array<bool> &bnd_flags,
                     Array<int> &bc_bnds)
{
   int amax = bnd_flags.Size();
   // Check strong boundaries
   for (int b = 0; b < bc_bnds.Size(); b++)
   {
      int bnd = bc_bnds[b];
      if ( bnd < 0 || bnd > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flags[bnd])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flags[bnd] = true;
   }
}

int main(int argc, char *argv[])
{
   // Initialize MPI and HYPRE and print info
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   printInfo();

   // Parse command-line options.
   OptionsParser args(argc, argv);

   // Mesh and discretization parameters
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   int order = 1;
   int ref_levels = 0;
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");

   // Problem parameters
   Array<int> strong_bdr;
   Array<int> weak_bdr;

   Array<int> master_bdr;
   Array<int> slave_bdr;

   const char *lib_file = "libfun.so";
   real_t mu_param = 0.1;

   args.AddOption(&strong_bdr, "-sbc", "--strong-bdr",
                  "Boundaries where Dirichelet BCs are enforced strongly.");
   args.AddOption(&weak_bdr, "-wbc", "--weak-bdr",
                  "Boundaries where Dirichelet BCs are enforced weakly.");
   args.AddOption(&master_bdr, "-mbc", "--master-bdr",
                  "Periodic master boundaries.");
   args.AddOption(&slave_bdr, "-sbc", "--slave-bdr",
                  "Periodic slave boundaries.");
   args.AddOption(&lib_file, "-l", "--lib",
                  "Library file for case specific function definitions:\n\t"
                  " - Initial condition\n\t"
                  " - Boundary condition\n\t"
                  " - Forcing\n\t"
                  " - Diffusion\n\t");
   args.AddOption(&mu_param, "-mu", "--dyn-visc",
                  "Sets the dynamic diffusion parameters, should be positive.");

   // Artificial diffusion parameters
   real_t kdc0 = 0.0;
   real_t kdc1 = 0.1;

   args.AddOption(&kdc0, "-k0", "--dc0",
                  "Inconsistent diffusion parameter of the NS formulation.");
   args.AddOption(&kdc1, "-k1", "--kdc1",
                  "Consistent diffusion parameter of tha NS formulation.");

   // Solver parameters
   double GMRES_RelTol = 1e-5;
   int    GMRES_MaxIter = 500;
   double Newton_RelTol = 1e-5;
   int    Newton_MaxIter = 10;

   args.AddOption(&GMRES_RelTol, "-lt", "--linear-tolerance",
                  "Relative tolerance for the GMRES solver.");
   args.AddOption(&GMRES_MaxIter, "-li", "--linear-itermax",
                  "Maximum iteration count for the GMRES solver.");
   args.AddOption(&Newton_RelTol, "-nt", "--newton-tolerance",
                  "Relative tolerance for the Newton solver.");
   args.AddOption(&Newton_MaxIter, "-ni", "--newton-itermax",
                  "Maximum iteration count for the Newton solver.");

   // Solution input/output params
   const char *vis_dir = "solution";
   args.AddOption(&vis_dir, "-vd", "--vis-dir",
                  "Directory for visualization files.\n\t");

   // Parse parameters
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   // Read the mesh from the given mesh file.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // Refine mesh
   {
      if (mesh.NURBSext && (strlen(ref_file) != 0))
      {
         mesh.RefineNURBSFromFile(ref_file);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
      if (Mpi::Root()) { mesh.PrintInfo(); }
   }

   // Partition mesh
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();

   // Boundary conditions
   if (Mpi::Root())
   {
      if (strong_bdr.Size()>0) {cout<<"Strong  = "; strong_bdr.Print();}
      if (weak_bdr.Size()>0) {cout<<"Weak    = "; weak_bdr.Print();}
      if (master_bdr.Size()>0) {cout<<"Periodic (master) = "; master_bdr.Print();}
      if (slave_bdr.Size()>0) {cout<<"Periodic (slave)  = "; slave_bdr.Print();}
   }

   Array<bool> bnd_flag(pmesh.bdr_attributes.Max()+1);
   bnd_flag = true;
   for (int b = 0; b < pmesh.bdr_attributes.Size(); b++)
   {
      bnd_flag[pmesh.bdr_attributes[b]] = false;
   }
   CheckBoundaries(bnd_flag, strong_bdr);
   CheckBoundaries(bnd_flag, weak_bdr);
   CheckBoundaries(bnd_flag, master_bdr);
   CheckBoundaries(bnd_flag, slave_bdr);

   MFEM_VERIFY(master_bdr.Size() == master_bdr.Size(),
               "Master-slave count do not match.");
   for (int b = 0; b < bnd_flag.Size(); b++)
   {
      MFEM_VERIFY(bnd_flag[b],
                  "Not all boundaries have a boundary condition set.");
   }

   // Define a finite element space on the mesh.
   FiniteElementCollection* fec = nullptr;
   NURBSExtension* ext = NULL;
   if (pmesh.NURBSext && order > 1)
   {
      fec = new NURBSFECollection(order);
      ext = new NURBSExtension(pmesh.NURBSext,order);
   }
   else
   {
      fec = new H1_FECollection(abs(order), dim);
   }

   ParFiniteElementSpace* space;
   space = new ParFiniteElementSpace(&pmesh, ext, fec);

   // Report the degree of freedoms used
   {
      Array<int> tdof(num_procs),udof(num_procs);
      tdof = 0;
      tdof[myid] = space->TrueVSize();
      MPI_Reduce(tdof.GetData(), udof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);


      int dof_t = space->GlobalTrueVSize();
      if (Mpi::Root())
      {
         mfem::out << "Number of finite element unknowns:"<< dof_t << endl;
         if (num_procs > 1)
         {
            mfem::out << "Number of finite element unknowns per partition:\n";
            udof.Print(mfem::out, num_procs);
         }
      }
   }

   // Define the gridfunction and solution vector
   ParGridFunction phi_gf(space);
   LibCoefficient sol_phi(lib_file, "sol_phi");
   phi_gf.ProjectCoefficient(sol_phi);
   Vector xp;
   phi_gf.GetTrueDofs(xp);

   // Define the visualisation output
   VisItDataCollection vdc("step", &pmesh);
   vdc.SetPrefixPath(vis_dir);
   vdc.RegisterField("phi", &phi_gf);
   vdc.SetCycle(0);
   vdc.Save();

   // Define the physical parameters
   LibVectorCoefficient adv(dim, lib_file, "advection");
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibCoefficient force(lib_file, "force");

   // Define the inverse estimate
   InverseEstimateCoefficient inv_est(space);

   // Define weak form and evolution
   StabConvDifIntegrator integrator(adv, mu, force, inv_est);
   ParNonlinearForm form(space);
   form.AddDomainIntegrator(&integrator);
   form.UseExternalIntegrators();

   Array<int> ess_bdr(space->GetMesh()->bdr_attributes.Max());
   ess_bdr = 0;
   for (int b = 0; b < strong_bdr.Size(); ++b)
   {
      ess_bdr[strong_bdr[b]-1] = 1;
   }
   ess_bdr.Print();
   form.SetEssentialBC(ess_bdr);
   //form.SetWeakBC   (weak_bdr);

   Solver* pc_mom  = nullptr;
   HypreILU* ilu_mom = new HypreILU();
   pc_mom  = ilu_mom;

   // Set up the Jacobian solver
   FGMRESSolver gmres(MPI_COMM_WORLD);
   gmres.iterative_mode = false;
   gmres.SetRelTol(GMRES_RelTol);
   gmres.SetMaxIter(GMRES_MaxIter);
   gmres.SetKDim(GMRES_MaxIter+1);
   gmres.SetPrintLevel(3);
   gmres.SetPreconditioner(*pc_mom);

   // Set up the Newton solver
   NewtonSolver newton_solver(MPI_COMM_WORLD);
   newton_solver.SetOperator(form);
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(1);
   newton_solver.SetRelTol(Newton_RelTol);
   newton_solver.SetMaxIter(Newton_MaxIter);
   newton_solver.SetSolver(gmres);

   // Solver nonlinear system
   Vector zero(space->TrueVSize());
   zero = 0.0;
   newton_solver.Mult(zero, xp);
   phi_gf.Distribute(xp);

   // Compute errors
   LibVectorCoefficient sol_grad(dim, lib_file, "grad_phi", false);
   if (sol_grad.Foundfunction())
   {
      int order_quad = max(2, 2*order+1);
      const IntegrationRule *irs[Geometry::NumGeom];
      for (int i=0; i < Geometry::NumGeom; ++i)
      {
         irs[i] = &(IntRules.Get(i, order_quad));
      }

      double err_phi  = phi_gf.ComputeL2Error(sol_phi, irs);
      double norm_phi = ComputeGlobalLpNorm(2., sol_phi, pmesh, irs);
      std::cout << "|| phi_h - phi_ex || / || phi_ex || = " << err_phi / norm_phi << "\n";

      err_phi  = phi_gf.ComputeGradError(&sol_grad, irs);
      norm_phi =  ComputeGlobalLpNorm(2., sol_grad, pmesh, irs);
      std::cout << "||grad phi_h - grad phi_ex || / || grad phi_ex || = " << err_phi / norm_phi << "\n";
   }

   // Write solution
   vdc.SetCycle(1);
   vdc.Save();

   {
      char vishost[] = "localhost";
      int visport = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << pmesh << phi_gf << flush;
   }

   // Free the used memory.
   delete fec;
   delete space;

   return 0;
}
