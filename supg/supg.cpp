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

#include <fenv.h>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

using namespace std;
using namespace mfem;

extern void printInfo();
extern void line(int len);

namespace fs = std::filesystem;
/**
 * Coefficient class for mesh-agnostic interpolation of a GridFunction.
 *
 * This class allows evaluating a GridFunction at arbitrary physical coordinates
 * by using an mfem::KDTree for spatial acceleration. It first finds the
 * element with the closest center and then checks its neighborhood.
 */
class GridFunctionInterpCoefficient : public Coefficient
{
private:
   const GridFunction *gf;
   KDTreeBase<int, real_t> *kdtree;
   Table *vtoel;
   mutable InverseElementTransformation inv_tr;

public:
   GridFunctionInterpCoefficient(const GridFunction *gf_)
      : gf(gf_), kdtree(nullptr), vtoel(nullptr)
   {
      inv_tr.SetInitialGuessType(InverseElementTransformation::ClosestRefNode);
      Mesh *mesh = gf->FESpace()->GetMesh();
      int sdim = mesh->SpaceDimension();
      if (sdim == 1) kdtree = new KDTree1D();
      else if (sdim == 2) kdtree = new KDTree2D();
      else if (sdim == 3) kdtree = new KDTree3D();

      if (kdtree)
      {
         Vector center(sdim);
         for (int i = 0; i < mesh->GetNE(); i++)
         {
            mesh->GetElementTransformation(i)->Transform(
               Geometries.GetCenter(mesh->GetElementBaseGeometry(i)), center);
            kdtree->AddPoint(center.GetData(), i);
         }
         kdtree->Sort();
      }
      vtoel = mesh->GetVertexToElementTable();
   }

   virtual ~GridFunctionInterpCoefficient()
   {
      delete kdtree;
      delete vtoel;
   }

   virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip)
   {
      Mesh *mesh = gf->FESpace()->GetMesh();
      int sdim = mesh->SpaceDimension();
      Vector x(sdim);
      T.Transform(ip, x);

      if (x.Norml2() < 1e-12) return 0.0; // Robustness for singular origin in NURBS

      IntegrationPoint ip_ref;
      if (kdtree && mesh->GetNE() > 0)
      {
         int closest_el = kdtree->FindClosestPoint(x.GetData());
         inv_tr.SetTransformation(*mesh->GetElementTransformation(closest_el));
         if (inv_tr.Transform(x, ip_ref) == InverseElementTransformation::Inside)
         {
            return gf->GetValue(closest_el, ip_ref);
         }

         if (vtoel)
         {
            Array<int> vertices;
            mesh->GetElementVertices(closest_el, vertices);
            for (int i = 0; i < vertices.Size(); i++)
            {
               int v = vertices[i];
               int ne = vtoel->RowSize(v);
               const int *els = vtoel->GetRow(v);
               for (int j = 0; j < ne; j++)
               {
                  int el = els[j];
                  if (el == closest_el) continue;
                  inv_tr.SetTransformation(*mesh->GetElementTransformation(el));
                  if (inv_tr.Transform(x, ip_ref) == InverseElementTransformation::Inside)
                  {
                     return gf->GetValue(el, ip_ref);
                  }
               }
            }
         }
      }

      // Fallback for robustness
      Array<int> elem_ids(1);
      Array<IntegrationPoint> ips(1);
      DenseMatrix point_mat(sdim, 1);
      for (int i=0; i<sdim; i++) { point_mat(i,0) = x(i); }
      mesh->FindPoints(point_mat, elem_ids, ips, false);

      if (elem_ids[0] >= 0)
      {
         return gf->GetValue(elem_ids[0], ips[0]);
      }
      return 0.0;
   }
};

bool fileExists(const std::string& filename) {
    return fs::exists(filename);
}

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
#ifdef __linux__
   fedisableexcept(FE_DIVBYZERO);
#endif
   printInfo();

   // Parse command-line options.
   OptionsParser args(argc, argv);

   // Mesh and discretization parameters
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   const char *ref_mesh_file = "";
   const char *ref_sol_file = "";
   const char *save_mesh_file = "";
   const char *save_sol_file = "";
   int order = 1;
   int ref_levels = 0;
   int plotter = 0;

   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");
   args.AddOption(&plotter, "-pl","--plot",
                  "Plot the dataset by using 1.");
   args.AddOption(&ref_mesh_file, "-rm", "--ref-mesh", "Reference mesh file.");
   args.AddOption(&ref_sol_file, "-rs", "--ref-solution", "Reference solution file.");
   args.AddOption(&save_mesh_file, "-sm", "--save-mesh", "File to save the current mesh.");
   args.AddOption(&save_sol_file, "-ss", "--save-solution", "File to save the current solution.");

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

   // Save current mesh for future reference (before partitioning and clearing)
   if (strlen(save_mesh_file) > 0 && Mpi::Root())
   {
      mesh.Save(save_mesh_file);
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
   FiniteElementCollection* ifec;
   ifec = new H1_FECollection(abs(order), dim);
   ParFiniteElementSpace* ispace;
   ispace = new ParFiniteElementSpace(&pmesh, ifec);
   ParGridFunction phi_igf(ispace);
   ParGridFunction err_igf(ispace);
   GridFunctionCoefficient phi_gf_cf(&phi_gf);
   phi_igf.ProjectCoefficient(phi_gf_cf);
   err_igf = 0.0;
   VisItDataCollection vdc("step", &pmesh);
   vdc.SetPrefixPath(vis_dir);
   vdc.RegisterField("phi", &phi_igf);
   vdc.RegisterField("error", &err_igf);
   vdc.SetCycle(0);
   vdc.Save();

   // Define the physical parameters
   LibVectorCoefficient adv(dim, lib_file, "advection");
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibCoefficient force(lib_file, "force");

   // Define weak form and evolution
   StabConvDifIntegrator integrator(adv, mu, force);
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

   // Load reference solution (if provided)
   Mesh *ref_mesh = nullptr;
   FiniteElementCollection *ref_fec = nullptr;
   FiniteElementSpace *ref_fes = nullptr;
   GridFunction *ref_gf = nullptr;
   Coefficient *ref_coeff = &sol_phi;

   if (strlen(ref_mesh_file) > 0 && strlen(ref_sol_file) > 0)
   {
      if (Mpi::Root()) { cout << "Loading reference mesh: " << ref_mesh_file << endl; }
      ref_mesh = new Mesh(ref_mesh_file, 1, 1);

      if (Mpi::Root()) { cout << "Loading reference solution: " << ref_sol_file << endl; }
      std::ifstream in(ref_sol_file);
      ref_gf = new GridFunction(ref_mesh, in);

      ref_coeff = new GridFunctionInterpCoefficient(ref_gf);
   }

   if (plotter==1)
   {
      int order_quad = max(2, 2*order+1);
      const IntegrationRule *irs[Geometry::NumGeom];
      for (int i=0; i < Geometry::NumGeom; ++i)
      {
         irs[i] = &(IntRules.Get(i, order_quad));
      }


      double l2_err_phi  = phi_gf.ComputeL2Error(*ref_coeff, irs);
      double norm_phi = ComputeGlobalLpNorm(2., *ref_coeff, pmesh, irs);
      double l2_err_norm = l2_err_phi/norm_phi;

      GradientGridFunctionCoefficient exgrad(&phi_gf);
      double h1_err_phi = phi_gf.ComputeH1Error(ref_coeff, &exgrad, irs);

      double h_local = integrator.GetMinH();
      double h;
      MPI_Reduce(&h_local, &h, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

      if (Mpi::Root())
      {
         std::cout << "L2 error: || phi_h - phi_ex || / || phi_ex || = " << l2_err_norm << "\n";
         std::cout << "H1 error: sqrt(norm_u^2+norm_du^2) = " << h1_err_phi << "\n";
         std::string filename = "plot.csv";
         if (fileExists(filename)) {
            std::ofstream outfile("plot.csv", std::ios::app);
            if (outfile.is_open()) {
               outfile << h << ',' << l2_err_norm << ',' << h1_err_phi << "\n";
               outfile.close();
            }
            else {
               std::cerr << "Unable to open file for appending.\n";
            }
         }
         else {
            std::ofstream file("plot.csv");
            if (file.is_open()) {
               file << 'h' << ',' << "L2_error" << ',' << "H1_error" << "\n";
               file << h << ',' << l2_err_norm << ',' << h1_err_phi << "\n";
               file.close();
            }
         }
      }
   }

   // Write solution
   phi_igf.ProjectCoefficient(phi_gf_cf);
   err_igf.ProjectCoefficient(*ref_coeff);
   err_igf -= phi_igf;
   for (int i = 0; i < err_igf.Size(); i++)
   {
      err_igf(i) = std::abs(err_igf(i));
   }
   vdc.SetCycle(1);
   vdc.Save();

   // Save solution for future reference
   if (strlen(save_sol_file) > 0)
   {
      phi_gf.SaveAsOne(save_sol_file);
   }

   {
      char vishost[] = "localhost";
      int visport = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << pmesh << phi_gf << flush;
   }

   // Free the used memory.
   if (ref_coeff != &sol_phi) { delete ref_coeff; }
   delete ref_gf;
   delete ref_fes;
   delete ref_fec;
   delete ref_mesh;

   delete fec;
   delete space;
   delete ifec;
   delete ispace;
   delete ilu_mom;

   return 0;
}
