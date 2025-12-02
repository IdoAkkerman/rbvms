// This file implements a small test to measure the variance of the stabilization
// parameter tau under permutations of triangle vertex ordering. The current
// method in supg/integrator.cpp uses GeomToPerfGeomJac, which should make tau
// invariant to vertex ordering for a single linear triangle because it maps to
// an equilateral (perfect) reference.

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

static double TauFromG(const DenseMatrix &G, const Vector &a, double k)
{
   // Tau formula copied from supg/integrator.cpp
   // Standalone, so we don't have to instantiate
   const double CI = 1.0/12.0;
   double tau_acc = 1e-10;
   const int dim = G.Width();
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_acc += G(i,j) * a[i] * a[j] + CI*CI * k * k * G(i,j) * G(i,j);
      }
   }
   return 1.0 / std::sqrt(tau_acc);
}

static void ComputeTausForMeshText(const std::string &mesh_text,
                                   double &tau_old, double &tau_new,
                                   const Vector &a, double k)
{
   std::istringstream is(mesh_text);
   Mesh mesh(is, 1, 0, false); // generate_nodes=1 to get Nodes if needed, do not fix orientation

   // Uncomment below to check for negative determinant
   // int ninv = mesh.CheckElementOrientation(false);
   // std::cout << "Inverted elements (without fixing): " << ninv << "\n";

   // One element mesh, just pick centroid of reference triangle
   ElementTransformation *T = mesh.GetElementTransformation(0);
   IntegrationPoint ip; ip.Set2(1.0/3.0, 1.0/3.0);
   T->SetIntPoint(&ip);

   const DenseMatrix &invJ = T->InverseJacobian();

   // Old method: G = invJ^T invJ
   DenseMatrix G_old(invJ.Width());
   MultAtB(invJ, invJ, G_old);
   tau_old = TauFromG(G_old, a, k);

   // New method using GeomToPerfGeomJac
   const DenseMatrix &A = Geometries.GetGeomToPerfGeomJac(T->GetGeometryType());
   DenseMatrix invJ_perf(invJ.Height(), invJ.Width());
   Mult(A, invJ, invJ_perf);
   DenseMatrix G_new(invJ.Width());
   MultAtB(invJ_perf, invJ_perf, G_new);
   tau_new = TauFromG(G_new, a, k);
}

int main(int argc, char** argv)
{
   // Advection and diffusion parameters to play around
   double ax = 1.0, ay = 0.5, k = 0.0;
   if (argc >= 3) { ax = atof(argv[1]); ay = atof(argv[2]); }
   if (argc >= 4) { k = atof(argv[3]); }

   Vector a(2); a[0] = ax; a[1] = ay;

   // Base vertex coordinates identical to mfem/data/ref-triangle.mesh
   const std::array<std::array<double,2>,3> base = {{{0.0,0.0},{1.0,0.0},{0.0,1.0}}};
   std::array<int,3> perm = {0,1,2};

   std::vector<double> taus_old, taus_new;
   std::vector<std::array<int,3>> perms_record;
   do {
      std::array<std::array<double,2>,3> v = { base[perm[0]], base[perm[1]], base[perm[2]] };
      std::string mesh_text = BuildSingleTriMeshText(v);
      double to=0.0, tn=0.0;
      ComputeTausForMeshText(mesh_text, to, tn, a, k);
      taus_old.push_back(to);
      taus_new.push_back(tn);
      perms_record.push_back(perm);
   } while (std::next_permutation(perm.begin(), perm.end()));

   auto stats = [](const std::vector<double> &x){
      double mean = 0.0; for (double v : x) mean += v; mean /= x.size();
      double var = 0.0; for (double v : x) { double d=v-mean; var += d*d; } var /= x.size();
      double mn = x[0], mx = x[0];
      for (double v : x) { if (v<mn) mn=v; if (v>mx) mx=v; }
      return std::tuple<double,double,double,double>(mean,var,mn,mx);
   };

   auto [m_old,v_old,min_old,max_old] = stats(taus_old);
   auto [m_new,v_new,min_new,max_new] = stats(taus_new);

   cout << "Tau variance test over all 6 permutations of a single triangle" << endl;
   cout << "Advective velocity a = (" << a[0] << ", " << a[1] << ")";
   cout << ", k = " << k << endl;

   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << m_old << ", var=" << v_old
        << ", min=" << min_old << ", max=" << max_old << endl;
   cout << "New method (perf geom):   mean=" << m_new << ", var=" << v_new
        << ", min=" << min_new << ", max=" << max_new << endl;

   cout << "\nValues per permutation (p0 p1 p2):" << endl;
   for (size_t i = 0; i < perms_record.size(); i++)
   {
      const auto &p = perms_record[i];
      cout << "perm (" << p[0] << ' ' << p[1] << ' ' << p[2] << ") : old="
           << taus_old[i] << ", new=" << taus_new[i] << endl;
   }

   // Return non-zero if new variance is unexpectedly large compared to old.
   // Threshold is arbitrary but ensures the test can be used in CI if desired.
   return (v_new > 1e-28 && v_new > 1e-6 * v_old) ? 1 : 0;
}
