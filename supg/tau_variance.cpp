// This file implements a small test to measure the variance of the stabilization
// parameter tau under permutations of element vertex ordering. The current
// method in supg/integrator.cpp uses GeomToPerfGeomJac, which should make tau
// invariant to vertex ordering for a single linear element because it maps to
// a perfect reference (equilateral triangle / regular tetrahedron, etc.).

#include "mfem.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <sstream>

using namespace mfem;
using std::cout;
using std::endl;

static std::string BuildSingleTriMeshText(const std::array<std::array<double,2>,3> &v)
{
   // Reproduce mfem/data/ref-triangle.mesh layout but with custom vertex order
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n2\n\n";
   os << "elements\n1\n";
   os << "1 2 0 1 2\n\n";
   os << "boundary\n3\n";
   os << "1 1 0 1\n";
   os << "2 1 1 2\n";
   os << "3 1 2 0\n\n";
   os << "vertices\n3\n2\n";
   os.setf(std::ios::fixed); os.precision(15);
   for (int i = 0; i < 3; i++)
   {
      os << v[i][0] << ' ' << v[i][1] << "\n";
   }
   os << "\n";
   return os.str();
}

static std::string BuildSingleTetMeshText(const std::array<std::array<double,3>,4> &v)
{
   // Reproduce mfem/data/ref-tetrahedron.mesh layout but with custom vertex order
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n3\n\n";
   os << "elements\n1\n";
   // Element: attribute=1, geometry=4 (TETRAHEDRON), vertices 0 1 2 3
   os << "1 4 0 1 2 3\n\n";
   // Boundary: 4 triangular faces with attributes 1..4
   os << "boundary\n4\n";
   os << "1 2 1 2 3\n";
   os << "2 2 0 3 2\n";
   os << "3 2 0 1 3\n";
   os << "4 2 0 2 1\n\n";
   // Vertices (4) with dimension flag 3
   os << "vertices\n4\n3\n";
   os.setf(std::ios::fixed); os.precision(15);
   for (int i = 0; i < 4; i++)
   {
      os << v[i][0] << ' ' << v[i][1] << ' ' << v[i][2] << "\n";
   }
   os << "\n";
   return os.str();
}

static double TauFromG(const DenseMatrix &G, const Vector &a, double k, double Ch)
{
   // Mirror supg/integrator.cpp: StabConvDifIntegrator::GetTau
   // Production code does: Ch = 1/Ch; and then adds (Ch*Ch*k*k) inside the i,j loop
   // i.e. acc = 1e-10 + sum_{i,j} ( G(i,j) a[i] a[j] + (1/Ch)^2 k^2 )
   double invCh = 1.0 / Ch;
   double tau_acc = 1e-10;
   const int dim = G.Width();
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_acc += G(i,j) * a[i] * a[j] + invCh * invCh * k * k;
      }
   }
   return 1.0 / std::sqrt(tau_acc);
}

static void ComputeTausForMeshText(const std::string &mesh_text,
                                   double &tau_old, double &tau_new,
                                   const Vector &a, double k, double Ch)
{
   std::istringstream is(mesh_text);
   Mesh mesh(is, 1, 0, false);

   // Uncomment below to check for negative determinant
   // int ninv = mesh.CheckElementOrientation(false);
   // std::cout << "Inverted elements (without fixing): " << ninv << "\n";

   // One element mesh, pick centroid of reference element (triangle/tetra)
   ElementTransformation *T = mesh.GetElementTransformation(0);
   IntegrationPoint ip;
   if (mesh.Dimension() == 2)
   {
      ip.Set2(1.0/3.0, 1.0/3.0);
   }
   else if (mesh.Dimension() == 3)
   {
      ip.Set3(1.0/4.0, 1.0/4.0, 1.0/4.0);
   }
   else
   {
      MFEM_ABORT("Unsupported mesh dimension in tau_variance test.");
   }
   T->SetIntPoint(&ip);

   const DenseMatrix &invJ = T->InverseJacobian();

   // Old method: G = invJ^T invJ
   DenseMatrix G_old(invJ.Width());
   MultAtB(invJ, invJ, G_old);
   tau_old = TauFromG(G_old, a, k, Ch);

   // New method using GeomToPerfGeomJac
   const DenseMatrix &A = Geometries.GetGeomToPerfGeomJac(T->GetGeometryType());
   DenseMatrix invJ_perf(invJ.Height(), invJ.Width());
   Mult(A, invJ, invJ_perf);
   DenseMatrix G_new(invJ.Width());
   MultAtB(invJ_perf, invJ_perf, G_new);
   tau_new = TauFromG(G_new, a, k, Ch);
}

// -------------------- Small utilities and geometry-specific runners --------------------

struct Stats
{
   double mean = 0.0;
   double var  = 0.0;
   double mn   = 0.0;
   double mx   = 0.0;
};

static Stats ComputeStats(const std::vector<double> &x)
{
   MFEM_VERIFY(!x.empty(), "ComputeStats requires non-empty input");
   Stats s{};
   for (double v : x) s.mean += v;
   s.mean /= x.size();
   for (double v : x) { const double d = v - s.mean; s.var += d*d; }
   s.var /= x.size();
   s.mn = x[0]; s.mx = x[0];
   for (double v : x) { if (v < s.mn) s.mn = v; if (v > s.mx) s.mx = v; }
   return s;
}

template <size_t N>
static void PrintPermutationValues(const std::vector<std::array<int,N>> &perms,
                                   const std::vector<double> &taus_old,
                                   const std::vector<double> &taus_new)
{
   using std::cout; using std::endl;
   cout << "\nValues per permutation (";
   for (size_t i = 0; i < N; i++) { cout << 'p' << i << (i+1<N ? ' ' : ')'); }
   cout << ":" << endl;
   for (size_t i = 0; i < perms.size(); i++)
   {
      cout << "perm (";
      for (size_t j = 0; j < N; j++) { cout << perms[i][j] << (j+1<N ? ' ' : ')'); }
      cout << " : old=" << taus_old[i] << ", new=" << taus_new[i] << endl;
   }
}

static int RunTriangleTest(const Vector &a2, double k, double Ch)
{
   using std::cout; using std::endl;

   const std::array<std::array<double,2>,3> base = {{{0.0,0.0},{1.0,0.0},{0.0,1.0}}};
   std::array<int,3> perm = {0,1,2};

   std::vector<double> taus_old, taus_new;
   std::vector<std::array<int,3>> perms_record;

   do {
      std::array<std::array<double,2>,3> v = { base[perm[0]], base[perm[1]], base[perm[2]] };
      std::string mesh_text = BuildSingleTriMeshText(v);
      double to = 0.0, tn = 0.0;
      ComputeTausForMeshText(mesh_text, to, tn, a2, k, Ch);
      taus_old.push_back(to); taus_new.push_back(tn); perms_record.push_back(perm);
   } while (std::next_permutation(perm.begin(), perm.end()));

   const Stats so = ComputeStats(taus_old);
   const Stats sn = ComputeStats(taus_new);

   cout << "Tau variance test over all 6 permutations of a single triangle" << endl;
   cout << "Advective velocity a = (" << a2[0] << ", " << a2[1] << ")";
   cout << ", k = " << k << ", Ch = " << Ch << endl;

   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << so.mean << ", var=" << so.var
        << ", min=" << so.mn << ", max=" << so.mx << endl;
   cout << "New method (perf geom):   mean=" << sn.mean << ", var=" << sn.var
        << ", min=" << sn.mn << ", max=" << sn.mx << endl;

   PrintPermutationValues(perms_record, taus_old, taus_new);

   // Same failure criterion used previously
   return (sn.var > 1e-28 && sn.var > 1e-6 * so.var) ? 1 : 0;
}

static int RunTetrahedronTest(const Vector &a3, double k, double Ch)
{
   using std::cout; using std::endl;

   const std::array<std::array<double,3>,4> base = {{{{0.0,0.0,0.0}},{{1.0,0.0,0.0}},{{0.0,1.0,0.0}},{{0.0,0.0,1.0}}}};
   std::array<int,4> perm = {0,1,2,3};

   std::vector<double> taus_old, taus_new;
   std::vector<std::array<int,4>> perms_record;

   do {
      std::array<std::array<double,3>,4> v = { base[perm[0]], base[perm[1]], base[perm[2]], base[perm[3]] };
      std::string mesh_text = BuildSingleTetMeshText(v);
      double to = 0.0, tn = 0.0;
      ComputeTausForMeshText(mesh_text, to, tn, a3, k, Ch);
      taus_old.push_back(to); taus_new.push_back(tn); perms_record.push_back(perm);
   } while (std::next_permutation(perm.begin(), perm.end()));

   const Stats so = ComputeStats(taus_old);
   const Stats sn = ComputeStats(taus_new);

   cout << "\nTau variance test over all 24 permutations of a single tetrahedron" << endl;
   cout << "Advective velocity a = (" << a3[0] << ", " << a3[1] << ", " << a3[2] << ")";
   cout << ", k = " << k << ", Ch = " << Ch << endl;

   cout << "Old method (raw invJ):    mean=" << so.mean << ", var=" << so.var
        << ", min=" << so.mn << ", max=" << so.mx << endl;
   cout << "New method (perf geom):   mean=" << sn.mean << ", var=" << sn.var
        << ", min=" << sn.mn << ", max=" << sn.mx << endl;

   PrintPermutationValues(perms_record, taus_old, taus_new);

   return (sn.var > 1e-28 && sn.var > 1e-6 * so.var) ? 1 : 0;
}

int main(int argc, char** argv)
{
   // Advection and diffusion parameters to play around
   double ax = 1.0, ay = 0.5, k = 0.0, Ch = 12.0;
   if (argc >= 3) { ax = atof(argv[1]); ay = atof(argv[2]); }
   if (argc >= 4) { k = atof(argv[3]); }
   if (argc >= 5) { Ch = atof(argv[4]); }

   // Triangle
   Vector a2(2); a2[0] = ax; a2[1] = ay;
   int fail = RunTriangleTest(a2, k, Ch);

   // Tetrahedron
   Vector a3(3); a3[0] = ax; a3[1] = ay; a3[2] = 0.2;
   fail |= RunTetrahedronTest(a3, k, Ch);

   return fail;
}
