//                                MFEM Example 5
//
// Compile with: make ex5
//
// Sample runs:  ex5 -m ../data/square-disc.mesh
//               ex5 -m ../data/star.mesh
//               ex5 -m ../data/star.mesh -pa
//               ex5 -m ../data/beam-tet.mesh
//               ex5 -m ../data/beam-hex.mesh
//               ex5 -m ../data/beam-hex.mesh -pa
//               ex5 -m ../data/escher.mesh
//               ex5 -m ../data/fichera.mesh
//
// Device sample runs:
//               ex5 -m ../data/star.mesh -pa -d cuda
//               ex5 -m ../data/star.mesh -pa -d raja-cuda
//               ex5 -m ../data/star.mesh -pa -d raja-omp
//               ex5 -m ../data/beam-hex.mesh -pa -d cuda
//
// Description:  This example code solves a 3D mixed Stokes problem
//               corresponding to the saddle point system
//                                  ω  - curl u        = 0
//                                curl ω   +    grad p = f
//                                      -div u         = 0
//
//               with natural boundary condition u_t = <given tangential velocity>.
//               and essentail boundary condition u_n = <given normal velocity>.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace mfem;

// Define the forcing terms / boundary conditions
void fFun(const Vector & x, Vector & f);
double divFunc(const Vector & x);
void zero_velocity(const Vector & x, Vector & v);
void lid_velocity(const Vector & x, Vector & v);

class VectorAsColumnOperator : public Operator
{
private:
    const Vector &c;

public:
    VectorAsColumnOperator(const Vector &c_)
        : Operator(c_.Size(), 1), c(c_) {}

    void Mult(const Vector &x, Vector &y) const override
    {
        y.SetSize(c.Size());
        y = c;
        y *= x(0);
    }

    void MultTranspose(const Vector &x, Vector &y) const override
    {
        y.SetSize(1);
        y(0) = InnerProduct(c, x);
    }
};

class Block3x3LUPreconditioner : public Solver
{
public:
   Block3x3LUPreconditioner(const Array<int> & offsets_)
      : Solver(offsets_.Last()),
        nBlocks(offsets_.Size() - 1),
        offsets(0),
        ops(nBlocks, nBlocks)
   {
      MFEM_VERIFY(offsets_.Size() == 4, "Need 3 blocks");

      ops = static_cast<Operator *>(NULL);
      offsets.MakeRef(offsets_);

      tmp0.SetSize(offsets[1] - offsets[0]); // block 0
      tmp1.SetSize(offsets[2] - offsets[1]); // block 1
      tmp2.SetSize(offsets[3] - offsets[2]); // block 2

      coeffs.SetSize(nBlocks, nBlocks);
      coeffs = 1.0;
   }

   void SetBlock(int iRow, int iCol, Operator *op, double c = 1.0)
   {
      MFEM_VERIFY(offsets[iRow+1] - offsets[iRow] == op->NumRows() &&
                  offsets[iCol+1] - offsets[iCol] == op->NumCols(),
                  "incompatible Operator dimensions");

      ops(iRow, iCol) = op;
      coeffs(iRow, iCol) = c;
   }

   virtual void SetOperator(const Operator &op) { }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      sol.Update(y.GetData(), offsets);
      rhs.Update(x.GetData(), offsets);

      y = 0.0;

      // --- Step 1 ---
      ops(0, 0)->Mult(rhs.GetBlock(0), sol.GetBlock(0));
      sol.GetBlock(0) *= coeffs(0, 0);
      // std::cout << "updated omega\n";

      // --- Step 2 ---
      // ops(1, 0)->Mult(sol.GetBlock(0), tmp1);
      // tmp1.Add(-1.0, rhs.GetBlock(1));
      // tmp1 *= -1.0;
      // ops(1, 1)->Mult(tmp1, sol.GetBlock(1));
      tmp1.SetSize(ops(1, 0)->NumRows());
      // std::cout << "updated U0\n";
      ops(1, 0)->Mult(sol.GetBlock(0), tmp1);
      tmp1 *= coeffs(1, 0);
      // std::cout << "updated U1\n";

      Vector r1(rhs.GetBlock(1));  // deep copy
      // std::cout << "updated U2\n";
      r1 -= tmp1;
      // std::cout << "updated U3\n";

      ops(1, 1)->Mult(r1, sol.GetBlock(1));
      sol.GetBlock(1) *= coeffs(1, 1);
      // std::cout << "updated U4\n";

      // --- Step 3 ---
      // ops(2, 1)->Mult(sol.GetBlock(1), tmp2);
      // tmp2.Add(-1.0, rhs.GetBlock(2));
      // tmp2 *= -1.0;
      // ops(2, 2)->Mult(tmp2, sol.GetBlock(2));
      tmp2.SetSize(ops(2, 1)->NumRows());
      ops(2, 1)->Mult(sol.GetBlock(1), tmp2);
      tmp2 *= coeffs(2, 1);

      Vector r2(rhs.GetBlock(2));  // deep copy
      r2 -= tmp2;

      ops(2, 2)->Mult(r2, sol.GetBlock(2));
      sol.GetBlock(2) *= coeffs(2, 2);
      // std::cout << "updated p\n";

      // --- Back substitution ---

      // x1 = x1 - A12 x2
      // ops(1, 2)->Mult(sol.GetBlock(2), tmp1);
      // ops(1, 1)->Mult(tmp1, tmp1);
      // sol.GetBlock(1).Add(1.0, tmp1);
      tmp1.SetSize(ops(1, 2)->NumRows());
      ops(1, 2)->Mult(sol.GetBlock(2), tmp1);
      tmp1 *= coeffs(1, 2);

      // apply M2^{-1} safely (no in-place!)
      Vector tmp1b(tmp1);
      ops(1, 1)->Mult(tmp1b, tmp1);
      tmp1 *= coeffs(1, 1);

      sol.GetBlock(1) += tmp1;
      // std::cout << "back updated U\n";

      // x0 = x0 - A01 x1
      // ops(0, 1)->Mult(sol.GetBlock(1), tmp0);
      // ops(0, 0)->Mult(tmp0, tmp0);
      // sol.GetBlock(0).Add(1.0, tmp0);
      tmp0.SetSize(ops(0, 1)->NumRows());
      ops(0, 1)->Mult(sol.GetBlock(1), tmp0);
      tmp0 *= coeffs(0, 1);

      // apply M1^{-1} safely
      Vector tmp0b(tmp0);
      ops(0, 0)->Mult(tmp0b, tmp0);
      tmp0 *= coeffs(0, 0);

      sol.GetBlock(0) += tmp0;
      // std::cout << "back updated omega\n";
   }

private:
   int nBlocks;
   Array<int> offsets;
   Array2D<Operator *> ops;
   Array2D<double> coeffs;

   mutable BlockVector sol;
   mutable BlockVector rhs;

   mutable Vector tmp0, tmp1, tmp2;
};

class RotationalConvectionLFIntegrator
    : public LinearFormIntegrator
{
private:
    VectorCoefficient &omega;
    VectorCoefficient &u;

public:
    RotationalConvectionLFIntegrator(
        VectorCoefficient &omega_,
        VectorCoefficient &u_)
        : omega(omega_), u(u_) {}

    virtual void AssembleRHSElementVect(
        const FiniteElement &el,
        ElementTransformation &Tr,
        Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();

        elvect.SetSize(dof);
        elvect = 0.0;

        DenseMatrix vshape(dof, dim);

        Vector omega_val(dim);
        Vector u_val(dim);
        Vector cross(dim);

        // high quadrature order
        const IntegrationRule *ir =
            IntRule ? IntRule :
            &IntRules.Get(
                el.GetGeomType(),
                4 * el.GetOrder() + 4);

        for (int q = 0; q < ir->GetNPoints(); q++)
        {
            const IntegrationPoint &ip = ir->IntPoint(q);

            Tr.SetIntPoint(&ip);

            const double w = ip.weight * Tr.Weight();

            // RT basis in physical coordinates
            el.CalcVShape(Tr, vshape);

            omega.Eval(omega_val, Tr, ip);
            u.Eval(u_val, Tr, ip);

            // omega x u
            cross(0) =
                omega_val(1)*u_val(2)
              - omega_val(2)*u_val(1);

            cross(1) =
                omega_val(2)*u_val(0)
              - omega_val(0)*u_val(2);

            cross(2) =
                omega_val(0)*u_val(1)
              - omega_val(1)*u_val(0);

            for (int i = 0; i < dof; i++)
            {
                double val = 0.0;

                for (int d = 0; d < dim; d++)
                {
                    val += vshape(i,d) * cross(d);
                }

                elvect(i) += w * val;
            }
        }
    }
};

int main(int argc, char *argv[])
{
   StopWatch chrono;

   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   bool verbose = (myid == 0);

   // 2. Parse command-line options.
   const char *mesh_file = "/home/suyash/Uni/postdoc/rbvms/mfem/data/inline-hex.mesh";
   int order = 1;
   bool par_format = false;
   const char *device_config = "cpu";
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&par_format, "-pf", "--parallel-format", "-sf",
                  "--serial-format",
                  "Format to use when saving the results for VisIt.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (verbose)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (verbose)
   {
      args.PrintOptions(cout);
   }

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 5. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 10,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   {
      int par_ref_levels = 1;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }

   double dt = 0.01;
   int Ndt = 1000;
   int saveInterval = 10;
   double Re = 1000;
   double picard_tol = 1e-2;

   // 7. Define a finite element space on the mesh. Full de Rham complex
   FiniteElementCollection *hcurl_coll(new ND_FECollection(order, dim));
   FiniteElementCollection *hdiv_coll(new RT_FECollection(order - 1, dim));
   FiniteElementCollection *l2_coll(new L2_FECollection(order - 1, dim));

   ParFiniteElementSpace *N_space = new ParFiniteElementSpace(pmesh, hcurl_coll);
   ParFiniteElementSpace *R_space = new ParFiniteElementSpace(pmesh, hdiv_coll);
   ParFiniteElementSpace *W_space = new ParFiniteElementSpace(pmesh, l2_coll);

   HYPRE_BigInt dimN = N_space->GlobalTrueVSize();
   HYPRE_BigInt dimR = R_space->GlobalTrueVSize();
   HYPRE_BigInt dimW = W_space->GlobalTrueVSize();

   // 8. Define the BlockStructure of the problem, i.e. define the array of
   //    offsets for each variable. The last component of the Array is the sum
   //    of the dimensions of each block.
   Array<int> block_offsets(4); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = N_space->GetVSize();
   block_offsets[2] = R_space->GetVSize();
   block_offsets[3] = W_space->GetVSize();
   // block_offsets[4] = 1;
   block_offsets.PartialSum();

   Array<int> block_trueOffsets(4); // number of variables + 1
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = N_space->TrueVSize();
   block_trueOffsets[2] = R_space->TrueVSize();
   block_trueOffsets[3] = W_space->TrueVSize();
   // block_trueOffsets[4] = 1;
   block_trueOffsets.PartialSum();

   if (verbose)
   {
   std::cout << "***********************************************************\n";
   std::cout << "dim(N) = " << dimN << "\n";
   std::cout << "dim(R) = " << dimR << "\n";
   std::cout << "dim(W) = " << dimW << "\n";
   std::cout << "dim(N+R+W) = " << dimN + dimR + dimW << "\n";
   std::cout << "dim mesh = " << dim << "\n";
   std::cout << "***********************************************************\n";
   // 9. Define the rhs of the PDE.
   // cout << "Boundary attributes: ";
   // mesh->bdr_attributes.Print(cout);
   }

   ParGridFunction omega(N_space);
   ParGridFunction u(R_space);
   ParGridFunction p(W_space);

   omega = 0.0;
   u = 0.0;
   p = 0.0;

   VectorFunctionCoefficient zero_vel(dim, zero_velocity);
   VectorFunctionCoefficient fcoeff(dim, fFun);
   FunctionCoefficient div_u(divFunc);
   VectorFunctionCoefficient lid_vel(dim, lid_velocity);

   // Weak bcs
   Array<int> wall_marker(pmesh->bdr_attributes.Max());
   wall_marker = 1;   // mark all boundaries
   // wall_marker[0] = 0;   // example: attribute 5 is not a wall

   // Strong bcs
   Array<int> ess_tdof_list;
   R_space->GetEssentialTrueDofs(wall_marker, ess_tdof_list);
   u.ProjectBdrCoefficient(zero_vel, wall_marker);

   Array<int> lid_marker(pmesh->bdr_attributes.Max());
   lid_marker = 0;
   lid_marker[5] = 1;   // example: attribute 5

   // 10. Allocate memory (x, rhs)
   MemoryType mt = device.GetMemoryType();
   BlockVector x(block_offsets, mt), rhs(block_offsets, mt);
   BlockVector trueX(block_trueOffsets, mt), trueRhs(block_trueOffsets, mt);

   ParLinearForm *fform(new ParLinearForm);
   fform->Update(R_space, rhs.GetBlock(1), 0);
   fform->AddDomainIntegrator(new VectorFEDomainLFIntegrator(fcoeff));
   fform->Assemble();
   fform->SyncAliasMemory(rhs);
   fform->ParallelAssemble(trueRhs.GetBlock(1));
   trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);
   
   ParLinearForm *bform(new ParLinearForm);
   bform->Update(N_space, rhs.GetBlock(0), 0);
   bform->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(lid_vel), lid_marker);
   bform->Assemble();
   bform->SyncAliasMemory(rhs);
   bform->ParallelAssemble(trueRhs.GetBlock(0));
   trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);

   ParLinearForm *dform(new ParLinearForm);
   dform->Update(W_space, rhs.GetBlock(2), 0);
   dform->AddDomainIntegrator(new DomainLFIntegrator(div_u));
   dform->Assemble();
   dform->SyncAliasMemory(rhs);
   dform->ParallelAssemble(trueRhs.GetBlock(2));
   trueRhs.GetBlock(2).SyncAliasMemory(trueRhs);

   // ParLinearForm c_form(W_space);
   // ConstantCoefficient one(1.0);
   // c_form.AddDomainIntegrator(new DomainLFIntegrator(one));
   // c_form.Assemble();
   // Vector c = *c_form.ParallelAssemble();
   // trueRhs.GetBlock(3) = 0.0;

   // Operator *C_lambda = new VectorAsColumnOperator(c);
   // Operator *C_lambda_T = new TransposeOperator(C_lambda);

   if (verbose){
      std::cout << "assembled right-hand side" << "\n";
   }

   // 11. Assemble the finite element matrices for the Stokes operator

   ParBilinearForm *M1(new ParBilinearForm(N_space));
   ParBilinearForm *M2(new ParBilinearForm(R_space));
   ParBilinearForm *M3(new ParBilinearForm(W_space));
   ParMixedBilinearForm *M2_E21(new ParMixedBilinearForm(N_space, R_space));
   ParMixedBilinearForm *M3_E32(new ParMixedBilinearForm(R_space, W_space));

   M1->AddDomainIntegrator(new VectorFEMassIntegrator());
   M1->Assemble();
   M1->Finalize();

   if (verbose){
      std::cout << "assembled M1" << "\n";
   }

   M2->AddDomainIntegrator(new VectorFEMassIntegrator());
   M2->Assemble();
   M2->SpMat() *= 1.0/dt;
   M2->EliminateEssentialBC(ess_tdof_list, u, *fform);
   M2->Finalize();

   if (verbose){
      std::cout << "assembled M2" << "\n";
   }

   M3->AddDomainIntegrator(new MassIntegrator());
   M3->Assemble();
   M3->Finalize();

   if (verbose){
      std::cout << "assembled M3" << "\n";
   }

   M2_E21->AddDomainIntegrator(new VectorFECurlIntegrator());
   M2_E21->Assemble();
   // M2_E21->SpMat() *= 0.5;
   M2_E21->EliminateTestDofs(ess_tdof_list);
   M2_E21->Finalize();

   if (verbose){
   std::cout << "assembled M2_E21" << "\n";
   }

   M3_E32->AddDomainIntegrator(new VectorFEDivergenceIntegrator());
   M3_E32->Assemble();
   M3_E32->EliminateTrialDofs(ess_tdof_list, u, *dform);
   M3_E32->Finalize();

   if (verbose){
   std::cout << "assembled M3_E32" << "\n";
   }

   trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);
   trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);
   trueRhs.GetBlock(2).SyncAliasMemory(trueRhs);

   BlockOperator stokesOp(block_trueOffsets);

   HypreParMatrix *A = NULL;
   HypreParMatrix *B = NULL;
   HypreParMatrix *C = NULL;
   HypreParMatrix *D = NULL;
   HypreParMatrix *Mp = NULL;

   TransposeOperator *Ct = NULL;
   TransposeOperator *Dt = NULL;

   A = M1->ParallelAssemble();
   B = M2->ParallelAssemble();
   Mp = M3->ParallelAssemble();
   C = M2_E21->ParallelAssemble();
   D = M3_E32->ParallelAssemble();

   Ct = new TransposeOperator(C);
   Dt = new TransposeOperator(D);

   stokesOp.SetBlock(0, 0, A);
   stokesOp.SetBlock(1, 1, B);
   stokesOp.SetBlock(0, 1, Ct, -1.0);
   stokesOp.SetBlock(1, 0, C, 1/(2*Re));
   stokesOp.SetBlock(1, 2, Dt, -1.0);
   stokesOp.SetBlock(2, 1, D);
   // stokesOp.SetBlock(2, 3, C_lambda);
   // stokesOp.SetBlock(3, 2, C_lambda_T);

   if (verbose){
      std::cout << "assembled stokes operator" << "\n";
   }

   // ------------------------------------------------------------------
   // 10. Build block preconditioner (parallel)
   // ------------------------------------------------------------------

   CGSolver *M1_inv = new CGSolver(MPI_COMM_WORLD);
   // HypreBoomerAMG *M1_inv = new HypreBoomerAMG;
   M1_inv->SetOperator(*A);          // A = vorticity mass matrix
   M1_inv->SetRelTol(1e-8);
   M1_inv->SetMaxIter(50);
   M1_inv->SetPrintLevel(0);
   M1_inv->iterative_mode = true;

   CGSolver *M2_inv = new CGSolver(MPI_COMM_WORLD);
   M2_inv->SetOperator(*B);          // B = velocity mass matrix
   M2_inv->SetRelTol(1e-8);
   M2_inv->SetMaxIter(50);
   M2_inv->SetPrintLevel(0);
   M2_inv->iterative_mode = true;

   CGSolver *M3_inv = new CGSolver(MPI_COMM_WORLD);
   M3_inv->SetOperator(*Mp);          // Mp = pressure mass matrix
   M3_inv->SetRelTol(1e-8);
   M3_inv->SetMaxIter(50);
   M3_inv->SetPrintLevel(0);
   M3_inv->iterative_mode = true;

   ProductOperator *MuInvD = new ProductOperator(M2_inv, Dt, false, false);
   ProductOperator *S = new ProductOperator(D, MuInvD, false, false);

   CGSolver *S_inv = new CGSolver(MPI_COMM_WORLD);
   S_inv->SetOperator(*S);          // Shur complement
   S_inv->SetPreconditioner(*M3_inv);
   S_inv->SetRelTol(1e-5);
   S_inv->SetAbsTol(1e-6);
   S_inv->SetMaxIter(50);
   S_inv->SetPrintLevel(0);
   S_inv->iterative_mode = true;

   Block3x3LUPreconditioner *prec = new Block3x3LUPreconditioner(block_trueOffsets);

   prec->SetBlock(0, 0, M1_inv);
   prec->SetBlock(1, 1, M2_inv);
   prec->SetBlock(2, 2, S_inv);
   prec->SetBlock(0, 1, Ct);
   prec->SetBlock(1, 0, C, 1/(2*Re));
   prec->SetBlock(1, 2, Dt);
   prec->SetBlock(2, 1, D);

   if (verbose){
      std::cout << "constructed block preconditioner\n";
   }

   // 14. Save data in the ParaView format
   ParaViewDataCollection paraview_dc("Example5-Parallel_NavierStokes", pmesh);
   paraview_dc.SetPrefixPath("ParaView_NavierStokes_par");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("vorticity", &omega);
   paraview_dc.RegisterField("velocity", &u);
   paraview_dc.RegisterField("pressure", &p);
   paraview_dc.Save();

   // 11. Solve the linear system A x = rhs using a iterative solver

   chrono.Clear();
   chrono.Start();
   GMRESSolver solver(MPI_COMM_WORLD);
   solver.SetOperator(stokesOp);
   solver.SetPreconditioner(*prec);
   solver.SetAbsTol(1e-10);
   // solver.SetRelTol(1e-10);
   solver.SetMaxIter(100);
   solver.SetPrintLevel(verbose);
   solver.SetKDim(100);
   solver.iterative_mode = true;

   Vector u_current(R_space->GetTrueVSize());
   Vector omega_current(N_space->GetTrueVSize());
   Vector M2_u_current(R_space->GetTrueVSize());
   Vector M2_E21_omega_current(R_space->GetTrueVSize());
   BlockVector trueRhs0(trueRhs);

   ParGridFunction omega_old(N_space);
   ParGridFunction u_old(R_space);

   // Midpoint fields
   ParGridFunction omega_mid(N_space);
   ParGridFunction u_mid(R_space);

   // Picard iteration storage
   ParGridFunction omega_prev_iter(N_space);
   ParGridFunction u_prev_iter(R_space);
   ParGridFunction u_diff(R_space);

   // Convection linear form
   ParLinearForm nform(R_space);

   // coefficients referencing midpoint fields
   VectorGridFunctionCoefficient omega_coeff(&omega_mid);
   VectorGridFunctionCoefficient u_coeff(&u_mid);

   nform.AddDomainIntegrator(new RotationalConvectionLFIntegrator(omega_coeff, u_coeff));

   // Temporary vectors
   Vector true_conv(R_space->TrueVSize());

   for (int ti = 0; ti <= Ndt; ti++)
   {
      if (verbose)
      {
         std::cout << "########  Time step " << ti + 1 << "/" << Ndt << " ########\n";
      }
      
      // --------------------------------------------------------------
      // Store previous timestep
      // --------------------------------------------------------------
      omega_old = omega;
      u_old = u;

      double resd = 1.0;

      while (resd > picard_tol)
      {
         trueRhs = trueRhs0;

         u.GetTrueDofs(u_current);
         omega.GetTrueDofs(omega_current);
         C->Mult(omega_current, M2_E21_omega_current);
         B->Mult(u_current, M2_u_current);

         trueRhs.GetBlock(1).Add(1.0, M2_u_current);
         trueRhs.GetBlock(1).Add(-1.0/(2*Re), M2_E21_omega_current);

         // Update midpoint fields
         omega_mid = omega_old;
         omega_mid += omega;
         omega_mid *= 0.5;

         u_mid = u_old;
         u_mid += u;
         u_mid *= 0.5;

         // Update nonlinear convection term
         nform = 0.0;
         nform.Assemble();
         nform.ParallelAssemble(true_conv);
         for (int i = 0; i < ess_tdof_list.Size(); i++)
         {
            true_conv(ess_tdof_list[i]) = 0.0;
         }
         trueRhs.GetBlock(1).Add(-1.0, true_conv);

         // Update solution
         omega_prev_iter = omega;
         u_prev_iter = u;

         // Solve system
         solver.Mult(trueRhs, trueX);

         if (verbose)
         {
         if (solver.GetConverged())
         {
            std::cout << "solver converged in " << solver.GetNumIterations()
                     << " iterations with a residual norm of "
                     << solver.GetFinalNorm() << ".\n";
         }
         else
         {
            std::cout << "solver did not converge in " << solver.GetNumIterations()
                     << " iterations. Residual norm is " << solver.GetFinalNorm()
                     << ".\n";
         }
         }

         omega.Distribute(&trueX.GetBlock(0));
         u.Distribute(&trueX.GetBlock(1));
         p.Distribute(&trueX.GetBlock(2));

         u_diff = u;
         u_diff -= u_prev_iter;
         resd = u_diff.Norml2();

         if (verbose)
         {
            std::cout << "######## Picard residual: " << resd << " ########\n";
         }
      }

      if ((ti) % saveInterval == 0)
      {
         paraview_dc.SetCycle(ti + 1);
         paraview_dc.SetTime(ti*dt);
         paraview_dc.Save();
      }

   }
   chrono.Stop();

   if (verbose)
   {
   std::cout << "solver took " << chrono.RealTime() << "s.\n";
   }


   // 15. Free the used memory.
   delete fform;
   delete bform;
   delete dform;
   delete A;
   delete B;
   delete C;
   delete D;
   delete Ct;
   delete Dt;
   delete M1;
   delete M2;
   delete M2_E21;
   delete M3_E32;
   delete N_space;
   delete R_space;
   delete W_space;
   delete hdiv_coll;
   delete hcurl_coll;
   delete l2_coll;
   delete pmesh;
   delete S;
   delete S_inv;
   delete M1_inv;
   delete M2_inv;
   delete M3_inv;
   delete prec;
   // delete C_lambda;
   // delete C_lambda_T;
   // // delete MuInvD;
   // delete DtMuInvD;

   return 0;
}

void fFun(const Vector & x, Vector & f)
{
   f.SetSize(3);

   f(0) = 0*sin(2*M_PI*x[0])*sin(2*M_PI*x[1])*sin(2*M_PI*x[2]);
   f(1) = 0*x[1];
   f(2) = 0*x[2];
}

void zero_velocity(const Vector &x, Vector &v)
{
    v.SetSize(3);
    v(0) = 0.0;
    v(1) = 0.0;
    v(2) = 0.0;
}

void lid_velocity(const Vector &x, Vector &v)
{
    v.SetSize(3);
    v(0) = 1.0;
    v(1) = 0.0;
    v(2) = 0.0;
}

double divFunc(const Vector &x)
{
    return 0.0;
}
