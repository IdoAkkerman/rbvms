// This file implements a small test to measure the variance of the stabilization
// parameter tau under permutations of element vertex ordering. The current
// method in supg/integrator.cpp uses GeomToPerfGeomJac, which should make tau
// invariant to vertex ordering for a single linear element because it maps to
// a perfect reference (equilateral triangle / regular tetrahedron, etc.).

#include "mfem.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace mfem;
using std::cout;
using std::endl;


/// Create a mesh with one triangle, with different element vertices
static Mesh *BuildSingleTriMesh(const int shift)
{

   // Reproduce mfem/data/ref-triangle.mesh layout but with custom vertex order
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n2\n\n";
   os << "elements\n1\n";
   os << "1 2";
   for (int i = 0; i < 3; i++)
   {
      os << " "<<(i + shift)%3;
   }
   os << "\n\n";
   os << "boundary\n3\n";
   os << "1 1 0 1\n";
   os << "2 1 1 2\n";
   os << "3 1 2 0\n\n";
   os << "vertices\n3\n2\n";
   os << "0 0\n";
   os << "1 0\n";
   os << "0 1\n";
   os << "\n";
   std::istringstream is(os.str());
   return new Mesh(is, 1, 0, false);
}

/// Create a mesh with one tetrahadron, with different element vertices
static Mesh *BuildSingleTetMesh(const int shift)
{
   // Reproduce mfem/data/ref-tetrahedron.mesh layout but with custom vertex order
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n3\n\n";
   os << "elements\n1\n";
   // Element: attribute=1, geometry=4 (TETRAHEDRON), vertices 0 1 2 3
   os << "1 4";
   for (int i = 0; i < 4; i++)
   {
      os << " "<<(i + shift)%4;
   }
   os << "\n\n";
   // Boundary: 4 triangular faces with attributes 1..4
   os << "boundary\n4\n";
   os << "1 2 1 2 3\n";
   os << "2 2 0 3 2\n";
   os << "3 2 0 1 3\n";
   os << "4 2 0 2 1\n\n";
   // Vertices (4) with dimension flag 3
   os << "vertices\n4\n3\n";
   os << "0 0 0\n";
   os << "1 0 0\n";
   os << "0 1 0\n";
   os << "0 0 1\n";
   os << "\n";
   std::istringstream is(os.str());
   return new Mesh(is, 1, 0, false);
}

static double TauFromG(const DenseMatrix &G, const Vector &a)
{
   double tau_acc = 1e-10;
   const int dim = G.Width();
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_acc += G(i,j) * a[i] * a[j];
      }
   }
   return 1.0 / std::sqrt(tau_acc);
}

static void ComputeConvectiveTau(Mesh *mesh,
                                 double &tau_old, double &tau_new,
                                 const Vector &a)
{
   // Uncomment below to check for negative determinant
   // int ninv = mesh.CheckElementOrientation(false);
   // std::cout << "Inverted elements (without fixing): " << ninv << "\n";

   // One element mesh, pick a canonical interior point for the reference element
   ElementTransformation *T = mesh->GetElementTransformation(0);
   IntegrationPoint ip;
   switch (T->GetGeometryType())
   {
      case Geometry::TRIANGLE: ip.Set2(1.0/3.0, 1.0/3.0); break;                // centroid
      case Geometry::TETRAHEDRON: ip.Set3(1.0/4.0, 1.0/4.0, 1.0/4.0); break;    // centroid
      case Geometry::PRISM: ip.Set3(1.0/3.0, 1.0/3.0, 0.5); break;              // tri-centroid x mid-height
      case Geometry::PYRAMID: ip.Set3(0.375, 0.375, 0.25); break;               // ref pyramid centroid?
      default: MFEM_ABORT("Unsupported geometry in tau_variance_alt test.");
   }
   T->SetIntPoint(&ip);

   const DenseMatrix &invJ = T->InverseJacobian();

   // Old method: G = invJ^T invJ
   DenseMatrix G_old(invJ.Width());
   MultAtB(invJ, invJ, G_old);
   tau_old = TauFromG(G_old, a);

   // New method using GeomToPerfGeomJac
   const DenseMatrix &A = Geometries.GetGeomToPerfGeomJac(T->GetGeometryType());
   DenseMatrix invJ_perf(invJ.Height(), invJ.Width());
   Mult(A, invJ, invJ_perf);
   DenseMatrix G_new(invJ.Width());
   MultAtB(invJ_perf, invJ_perf, G_new);
   tau_new = TauFromG(G_new, a);
}

static void ComputeDiffusive(Mesh *mesh,
                             double &tr_old, double &tr_new,
                             double &nf_old, double &nf_new)
{
   // One element mesh, pick a canonical interior point for the reference element
   ElementTransformation *T = mesh->GetElementTransformation(0);
   IntegrationPoint ip;
   switch (T->GetGeometryType())
   {
      case Geometry::TRIANGLE: ip.Set2(1.0/3.0, 1.0/3.0); break;
      case Geometry::TETRAHEDRON: ip.Set3(1.0/4.0, 1.0/4.0, 1.0/4.0); break;
      case Geometry::PRISM: ip.Set3(1.0/3.0, 1.0/3.0, 0.5); break;
      case Geometry::PYRAMID: ip.Set3(0.375, 0.375, 0.25); break;
      default: MFEM_ABORT("Unsupported geometry in tau_variance_alt test.");
   }
   T->SetIntPoint(&ip);

   const DenseMatrix &invJ = T->InverseJacobian();

   // Old method: G = invJ^T invJ
   DenseMatrix G_old(invJ.Width());
   MultAtB(invJ, invJ, G_old);
   tr_old = G_old.Trace();
   nf_old = G_old.FNorm();

   // New method using GeomToPerfGeomJac
   const DenseMatrix &A = Geometries.GetGeomToPerfGeomJac(T->GetGeometryType());
   DenseMatrix invJ_perf(invJ.Height(), invJ.Width());
   Mult(A, invJ, invJ_perf);
   DenseMatrix G_new(invJ.Width());
   MultAtB(invJ_perf, invJ_perf, G_new);
   tr_new = G_new.Trace();
   nf_new = G_new.FNorm();
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


static void PrintPermutationValues(const std::vector<double> &taus_old,
                                   const std::vector<double> &taus_new)
{
   for (size_t i = 0; i < taus_old.size(); i++)
   {
      cout << "perm (";
      for (size_t j = 0; j < taus_old.size(); j++)
      {
         cout <<" "<< (i+j)%taus_old.size();
      }
      cout << " )";
      cout << " : old=" << taus_old[i] << ", new=" << taus_new[i] << endl;
   }
}

static int RunTriangleTest()
{
   Mesh *mesh[3];
   for (size_t i = 0; i < 3; i++)
   {
      mesh[i] = BuildSingleTriMesh(i);
   }

   int num_samples = 101;
   Vector a(2);
   std::ofstream out ("conv_tau_2D.dat");
   for (size_t i = 0; i < num_samples; i++)
   {
      double theta = (2*M_PI*i)/(num_samples-1);
      a[0] = sin(theta);
      a[1] = cos(theta);
      out<< theta ;
      double to, tn;
      for (size_t i = 0; i < 3; i++)
      {
         ComputeConvectiveTau(mesh[i], to, tn, a);
         out<< " "<<to<<" "<<tn;
      }
      out<<"\n";
   }
   out.close();

   std::vector<double> trs_old, trs_new;
   std::vector<double> nfs_old, nfs_new;
   double to, tn, no, nn;
   for (size_t i = 0; i < 3; i++)
   {
      ComputeDiffusive(mesh[i], to, tn, no, nn);
      trs_old.push_back(to);
      trs_new.push_back(tn);
      nfs_old.push_back(no);
      nfs_new.push_back(nn);
   }

   Stats sto = ComputeStats(trs_old);
   Stats stn = ComputeStats(trs_new);
   Stats sno = ComputeStats(nfs_old);
   Stats snn = ComputeStats(nfs_new);

   cout << "Trace variance test over all 3 permutations of a single triangle" << endl;

   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sto.mean << ", var=" << sto.var
        << ", min=" << sto.mn << ", max=" << sto.mx << endl;
   cout << "New method (perf geom):   mean=" << stn.mean << ", var=" << stn.var
        << ", min=" << stn.mn << ", max=" << stn.mx << endl;

   PrintPermutationValues(trs_old, trs_new);

   cout << "FNorm variance test over all 3 permutations of a single triangle" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sno.mean << ", var=" << sno.var
        << ", min=" << sno.mn << ", max=" << sno.mx << endl;
   cout << "New method (perf geom):   mean=" << snn.mean << ", var=" << snn.var
        << ", min=" << snn.mn << ", max=" << snn.mx << endl;

   PrintPermutationValues(nfs_old, nfs_new);

   for (size_t i = 0; i < 3; i++)
   {
      delete mesh[i];
   }
   // Same failure criterion used previously
   return 0;//(sn.var > 1e-28 && sn.var > 1e-6 * so.var) ? 1 : 0;
}

static int RunTetrahedronTest()
{
   Mesh *mesh[4];
   for (size_t i = 0; i < 4; i++)
   {
      mesh[i] = BuildSingleTetMesh(i);
   }

   int num_samples = 101;
   Vector a(3);

   std::ofstream out ("conv_tau_3D.dat");
   for (size_t ii = 0; ii < num_samples; ii++)
   {
      double phi = (2*M_PI*ii)/(num_samples-1);
      for (size_t i = 0; i < num_samples; i++)
      {
         double theta = (2*M_PI*i)/(num_samples-1);

         a[0] = sin(theta);
         a[1] = cos(theta)*sin(phi);
         a[2] = cos(theta)*cos(phi);
         out<< theta <<" " << phi;
         double to, tn;
         for (size_t i = 0; i < 4; i++)
         {
            ComputeConvectiveTau(mesh[i], to, tn, a);
            out<< " "<<to<<" "<<tn;
         }
         out<<"\n";
      }
      out<<"\n";
   }

   std::vector<double> trs_old, trs_new;
   std::vector<double> nfs_old, nfs_new;

   double to, tn, no, nn;
   for (size_t i = 0; i < 4; i++)
   {
      ComputeDiffusive(mesh[i], to, tn, no, nn);
      trs_old.push_back(to);
      trs_new.push_back(tn);
      nfs_old.push_back(no);
      nfs_new.push_back(nn);
   }

   Stats sto = ComputeStats(trs_old);
   Stats stn = ComputeStats(trs_new);
   Stats sno = ComputeStats(nfs_old);
   Stats snn = ComputeStats(nfs_new);


   cout << "Trace variance test over all 4 permutations of a single tetrahedron" << endl;

   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sto.mean << ", var=" << sto.var
        << ", min=" << sto.mn << ", max=" << sto.mx << endl;
   cout << "New method (perf geom):   mean=" << stn.mean << ", var=" << stn.var
        << ", min=" << stn.mn << ", max=" << stn.mx << endl;

   PrintPermutationValues(trs_old, trs_new);

   cout << "FNorm variance test over all 4 permutations of a single tetrahedron" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sno.mean << ", var=" << sno.var
        << ", min=" << sno.mn << ", max=" << sno.mx << endl;
   cout << "New method (perf geom):   mean=" << snn.mean << ", var=" << snn.var
        << ", min=" << snn.mn << ", max=" << snn.mx << endl;

   PrintPermutationValues(nfs_old, nfs_new);

   for (size_t i = 0; i < 4; i++)
   {
      delete mesh[i];
   }
   
   return 0;
}

/// Create a mesh with one right prism, with different element vertex permutations
static Mesh *BuildSinglePrismMesh(const int rot)
{
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n3\n\n";
   os << "elements\n1\n";
   // Prism local vertex order: bottom tri (0,1,2), top tri (3,4,5)
   // We only rotate cyclically the (0,1,2) and (3,4,5) blocks by 'rot' in [0,2]
   int bot[3] = {0,1,2};
   int top[3] = {3,4,5};
   int b0 = bot[(0+rot)%3], b1 = bot[(1+rot)%3], b2 = bot[(2+rot)%3];
   int t0 = top[(0+rot)%3], t1 = top[(1+rot)%3], t2 = top[(2+rot)%3];
   os << "1 6 " << b0 << " " << b1 << " " << b2 << " " << t0 << " " << t1 << " " << t2 << "\n\n";
   os << "boundary\n5\n";
   // bottom triangle
   os << "1 2 " << b0 << " " << b2 << " " << b1 << "\n";
   // top triangle
   os << "1 2 " << t0 << " " << t1 << " " << t2 << "\n";
   // side quads
   os << "1 3 " << b0 << " " << b1 << " " << t1 << " " << t0 << "\n";
   os << "1 3 " << b1 << " " << b2 << " " << t2 << " " << t1 << "\n";
   os << "1 3 " << b2 << " " << b0 << " " << t0 << " " << t2 << "\n\n";
   os << "vertices\n6\n3\n";
   os << "0 0 0\n";
   os << "1 0 0\n";
   os << "0 1 0\n";
   os << "0 0 1\n";
   os << "1 0 1\n";
   os << "0 1 1\n\n";
   std::istringstream is(os.str());
   return new Mesh(is, 1, 0, false);
}

/// Create a mesh with one right pyramid, with different element vertex permutations
static Mesh *BuildSinglePyramidMesh(const int rot)
{
   std::ostringstream os;
   os << "MFEM mesh v1.0\n\n";
   os << "dimension\n3\n\n";
   os << "elements\n1\n";
   // Pyramid local order: square base (0,1,2,3) + apex (4)
   // Permute by cyclic rotation of base only by 'rot' in [0,3]
   int base[4] = {0,1,2,3};
   int b0 = base[(0+rot)%4], b1 = base[(1+rot)%4], b2 = base[(2+rot)%4], b3 = base[(3+rot)%4];
   int apex = 4;
   os << "1 7 " << b0 << " " << b1 << " " << b2 << " " << b3 << " " << apex << "\n\n";
   os << "boundary\n5\n";
   // base quad
   os << "1 3 " << b3 << " " << b2 << " " << b1 << " " << b0 << "\n";
   // side triangles
   os << "1 2 " << b0 << " " << b1 << " " << apex << "\n";
   os << "1 2 " << b1 << " " << b2 << " " << apex << "\n";
   os << "1 2 " << b2 << " " << b3 << " " << apex << "\n";
   os << "1 2 " << b3 << " " << b0 << " " << apex << "\n\n";
   os << "vertices\n5\n3\n";
   os << "0 0 0\n";
   os << "1 0 0\n";
   os << "1 1 0\n";
   os << "0 1 0\n";
   os << "0 0 1\n\n";
   std::istringstream is(os.str());
   return new Mesh(is, 1, 0, false);
}

static int RunPrismTest()
{
   Mesh *mesh[3];
   for (size_t i = 0; i < 3; i++) { mesh[i] = BuildSinglePrismMesh((int)i); }

   int num_samples = 101;
   Vector a(3);

   std::ofstream out("conv_tau_prism.dat");
   for (size_t ii = 0; ii < num_samples; ii++)
   {
      double phi = (2*M_PI*ii)/(num_samples-1);
      for (size_t i = 0; i < num_samples; i++)
      {
         double theta = (2*M_PI*i)/(num_samples-1);
         a[0] = sin(theta);
         a[1] = cos(theta)*sin(phi);
         a[2] = cos(theta)*cos(phi);
         out << theta << " " << phi;
         double to, tn;
         for (size_t p = 0; p < 3; p++)
         {
            ComputeConvectiveTau(mesh[p], to, tn, a);
            out << " " << to << " " << tn;
         }
         out << "\n";
      }
      out << "\n";
   }

   std::vector<double> trs_old, trs_new;
   std::vector<double> nfs_old, nfs_new;
   double to, tn, no, nn;
   for (size_t i = 0; i < 3; i++)
   {
      ComputeDiffusive(mesh[i], to, tn, no, nn);
      trs_old.push_back(to); trs_new.push_back(tn);
      nfs_old.push_back(no); nfs_new.push_back(nn);
   }

   Stats sto = ComputeStats(trs_old);
   Stats stn = ComputeStats(trs_new);
   Stats sno = ComputeStats(nfs_old);
   Stats snn = ComputeStats(nfs_new);

   cout << "Trace variance test over 3 cyclic base-rotations of a single prism" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sto.mean << ", var=" << sto.var
        << ", min=" << sto.mn << ", max=" << sto.mx << endl;
   cout << "New method (perf geom):   mean=" << stn.mean << ", var=" << stn.var
        << ", min=" << stn.mn << ", max=" << stn.mx << endl;
   PrintPermutationValues(trs_old, trs_new);

   cout << "FNorm variance test over 3 cyclic base-rotations of a single prism" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sno.mean << ", var=" << sno.var
        << ", min=" << sno.mn << ", max=" << sno.mx << endl;
   cout << "New method (perf geom):   mean=" << snn.mean << ", var=" << snn.var
        << ", min=" << snn.mn << ", max=" << snn.mx << endl;
   PrintPermutationValues(nfs_old, nfs_new);

   for (size_t i = 0; i < 3; i++) { delete mesh[i]; }
   return 0;
}

static int RunPyramidTest()
{
   Mesh *mesh[4];
   for (size_t i = 0; i < 4; i++) { mesh[i] = BuildSinglePyramidMesh((int)i); }

   int num_samples = 101;
   Vector a(3);

   std::ofstream out("conv_tau_pyramid.dat");
   for (size_t ii = 0; ii < num_samples; ii++)
   {
      double phi = (2*M_PI*ii)/(num_samples-1);
      for (size_t i = 0; i < num_samples; i++)
      {
         double theta = (2*M_PI*i)/(num_samples-1);
         a[0] = sin(theta);
         a[1] = cos(theta)*sin(phi);
         a[2] = cos(theta)*cos(phi);
         out << theta << " " << phi;
         double to, tn;
         for (size_t p = 0; p < 4; p++)
         {
            ComputeConvectiveTau(mesh[p], to, tn, a);
            out << " " << to << " " << tn;
         }
         out << "\n";
      }
      out << "\n";
   }

   std::vector<double> trs_old, trs_new;
   std::vector<double> nfs_old, nfs_new;
   double to, tn, no, nn;
   for (size_t i = 0; i < 4; i++)
   {
      ComputeDiffusive(mesh[i], to, tn, no, nn);
      trs_old.push_back(to); trs_new.push_back(tn);
      nfs_old.push_back(no); nfs_new.push_back(nn);
   }

   Stats sto = ComputeStats(trs_old);
   Stats stn = ComputeStats(trs_new);
   Stats sno = ComputeStats(nfs_old);
   Stats snn = ComputeStats(nfs_new);

   cout << "Trace variance test over 4 cyclic base-rotations of a single pyramid" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sto.mean << ", var=" << sto.var
        << ", min=" << sto.mn << ", max=" << sto.mx << endl;
   cout << "New method (perf geom):   mean=" << stn.mean << ", var=" << stn.var
        << ", min=" << stn.mn << ", max=" << stn.mx << endl;
   PrintPermutationValues(trs_old, trs_new);

   cout << "FNorm variance test over 4 cyclic base-rotations of a single pyramid" << endl;
   cout.setf(std::ios::scientific); cout.precision(6);
   cout << "Old method (raw invJ):    mean=" << sno.mean << ", var=" << sno.var
        << ", min=" << sno.mn << ", max=" << sno.mx << endl;
   cout << "New method (perf geom):   mean=" << snn.mean << ", var=" << snn.var
        << ", min=" << snn.mn << ", max=" << snn.mx << endl;
   PrintPermutationValues(nfs_old, nfs_new);

   for (size_t i = 0; i < 4; i++) { delete mesh[i]; }
   return 0;
}

int main(int argc, char** argv)
{
   int fail = RunTriangleTest();
   fail |= RunTetrahedronTest();
   fail |= RunPrismTest();
   fail |= RunPyramidTest();

   return fail;
}
