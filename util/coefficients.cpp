// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "coefficients.hpp"
#include <dlfcn.h>

using namespace std;
using namespace mfem;

extern "C" void
dsygvx_(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO,
                           int *N, double *A, int *LDA, double *B, int *LDB,
                           double *VL, double *VU, int *IL, int *IU,
                           double *ABSTOL, int *M, double *W, double *Z,
                           int *LDZ, double *WORK, int *LWORK,int *IWORK,
                           int *IFAIL, int *INFO);

void dsygvx_Eigensystem(DenseMatrix &a, DenseMatrix &b,
                        Vector &ev, DenseMatrix *evect,
                        char RANGE, real_t VL, real_t VU, int IL, int IU)
{
#ifdef MFEM_USE_LAPACK
   ev.SetSize(a.Width());
   int       ITYPE    = 1;
   char      JOBZ     = 'N';
   char      UPLO     = 'U';
   int       N        = a.Width();
   real_t   *A        = new real_t[N*N];
   int       LDA      = N;
   real_t   *B        = new real_t[N*N];
   int       LDB      = N;
   real_t    ABSTOL   = 0.0;
   int       M;
   real_t   *W        = ev.GetData();
   real_t   *Z        = NULL;
   int       LDZ      = 1;
   real_t    QWORK;
   real_t   *WORK     = NULL;
   int       LWORK    = -1; // query optimal (double) workspace size
   int      *IWORK    = new int[5*N];
   int      *IFAIL    = new int[N];
   int       INFO;

   if (evect) // Compute eigenvectors too
   {
      evect->SetSize(N);

      JOBZ     = 'V';
      Z        = evect->Data();
      LDZ      = N;
   }

   int hw = a.Height() * a.Width();
   real_t *data_a = a.Data();
   real_t *data_b = b.Data();

   for (int i = 0; i < hw; i++)
   {
      A[i] = data_a[i];
      B[i] = data_b[i];
   }

//   MFEM_LAPACK_PREFIX()
   dsygvx_( &ITYPE, &JOBZ, &RANGE, &UPLO,
            &N, A, &LDA, B, &LDB, &VL, &VU, &IL, &IU,
            &ABSTOL, &M, W, Z, &LDZ, &QWORK, &LWORK,
            IWORK, IFAIL, &INFO );

   LWORK  = (int) QWORK;
   WORK  = new real_t[LWORK];

   //MFEM_LAPACK_PREFIX()
   dsygvx_( &ITYPE, &JOBZ, &RANGE, &UPLO,
            &N, A, &LDA, B, &LDB, &VL, &VU, &IL, &IU,
            &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK,
            IWORK, IFAIL, &INFO );

   if (INFO != 0)
   {
      mfem::err << "dsyevr_Eigensystem(...): DSYEVR error code: "
                << INFO << endl;
      mfem_error();
   }

   if (evect) // Compute eigenvectors too
   {
      evect->SetSize(N,M);
   }
   delete [] IFAIL;
   delete [] IWORK;
   delete [] WORK;
   delete [] B;
   if (evect == NULL) { delete [] A; }

#else
   MFEM_CONTRACT_VAR(a);
   MFEM_CONTRACT_VAR(ev);
   MFEM_CONTRACT_VAR(evect);
#endif
}

real_t Eigenvalue(DenseMatrix &a, DenseMatrix &b, int i = -1)
{
#ifdef MFEM_USE_LAPACK
   if (i < 0) { i = a.Width(); }
   Vector ev;
   dsygvx_Eigensystem(a, b, ev, NULL, 'I', 0.0, 0.0, i, i);
   return ev[0];

#else
  // MFEM_CONTRACT_VAR(ns);
  // MFEM_CONTRACT_VAR(tol);
   mfem_error("DenseMatrix::Eigenvalue: Compiled without LAPACK");
   return 0.0;
#endif
}

InverseEstimateCoefficient::InverseEstimateCoefficient(FiniteElementSpace *f)
   : fes(f), ir(NULL), Q(NULL)
{
   ComputeInverseEstimates();
}

InverseEstimateCoefficient::InverseEstimateCoefficient(FiniteElementSpace *f,
                                                       Coefficient &q)
   : fes(f), ir(NULL), Q(NULL)
{
   ComputeInverseEstimates();
}

GridFunction *InverseEstimateCoefficient::GetGridFunction()
{
   FiniteElementCollection* fec_ec = new L2_FECollection(0,
                                                         fes ->GetMesh()->Dimension());
   FiniteElementSpace *fes_ec = new FiniteElementSpace(fes ->GetMesh(), fec_ec);
   GridFunction *gf = new GridFunction(fes_ec, elemInvEst.GetData());
   gf->MakeOwner(fec_ec);
   return gf;
}

void InverseEstimateCoefficient::ComputeInverseEstimates()
{
   elemInvEst.SetSize(fes -> GetNE());
   SetIntRule(*fes->GetFE(0));
   for (int i = 0; i < fes -> GetNE(); i++)
   {
      elemInvEst[i] = ElementInverseEstimate(*fes->GetFE(i),
                                             *fes->GetElementTransformation(i));
   }
}

void InverseEstimateCoefficient::SetIntRule(const FiniteElement &el)
{
   ir = &IntRules.Get(el.GetGeomType(), 2*el.GetOrder());
}

real_t InverseEstimateCoefficient::ElementInverseEstimate(
   const FiniteElement &el,
   ElementTransformation &Trans)
{
   if (el.GetOrder() < 2)
   {
      return std::numeric_limits<real_t>::min();
   }

   int nd = el.GetDof();
   int dim = el.GetDim();

   shape.SetSize(nd);
   dshape.SetSize(nd,dim);
   laplace.SetSize(nd);

   lapmat.SetSize(nd,nd);
   bimat.SetSize(nd,nd);
   ovec.SetSize(nd);

   real_t w,q = 1.0;

   bimat = 0.0;
   lapmat = 0.0;
   ovec = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint(&ip);
      w = Trans.Weight()*ip.weight;
      if (Q)
      {
         q = Q->Eval(Trans, ip);
      }

      el.CalcPhysDShape(Trans, dshape);
      AddMult_a_AAt(w*q, dshape, lapmat);

      el.CalcPhysLaplacian(Trans, laplace);
      AddMult_a_VVt(w*q*q, laplace, bimat);

      el.CalcPhysShape(Trans, shape);
      ovec.Add(w, shape);
   }
   ovec *= 1.0/ovec.Norml2();

   // Correct nullspace
   AddMultVVt(ovec, lapmat);

   // Return largest eigenvalue
   return Eigenvalue(bimat, lapmat);
}

// Get the function from the library
void LibCoefficient::GetLibFunction(string libName,
                                    vector<string> funNames,
                                    bool required)
{
   // Open library
   libHandle = dlopen (libName.c_str(), RTLD_LAZY);
   if (!libHandle)
   {
      cout<<"Library "<<libName<<" not found."<<endl;
      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout <<"Functions:";
      for (int i = 0; i < funNames.size(); i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" will be set to val = "<<val<< endl;
      return;
   }

   // Search Function
   TDFunction = nullptr;
   for (int i = 0; i < funNames.size(); i++)
   {
      TDFunction = (TDFunPtr)dlsym(libHandle, funNames[i].c_str());
      if (TDFunction) { break; }
   }

   // Check if function is found
   if (!TDFunction)
   {
      cout <<"Functions:";
      for (int i = 0; i < funNames.size(); i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" can not be found in "<<libName<<endl;

      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout<<"Function will be set to val = "<<val<< endl;
   }
}

// Evaluate
real_t LibCoefficient::Eval(ElementTransformation &T,
                            const IntegrationPoint &ip)
{
   // Homegenous if not defined
   if (!TDFunction) { return val; }

   // Evaluate library function
   T.Transform(ip, x);
   return ((*TDFunction)(x.GetData(),x.Size(),GetTime()));
}

// Destructor
LibCoefficient::~LibCoefficient()
{
   if (libHandle) { dlclose(libHandle); }
}


// Get the function from the library
void LibVectorCoefficient::GetLibFunction(string libName,
                                          vector<string> funNames,
                                          bool required)
{
   // Open library
   libHandle = dlopen (libName.c_str(), RTLD_LAZY);
   if (!libHandle)
   {
      cout<<"Library "<<libName<<" not found."<<endl;
      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout <<"Functions:";
      for (int i = 0; i < funNames.size(); i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" will be set to homogenous."<<endl;
      return;
   }

   // Search Function
   TDFunction = nullptr;
   for (int i = 0; i < funNames.size(); i++)
   {
      TDFunction = (TDFunPtr)dlsym(libHandle, funNames[i].c_str());
      if (TDFunction) { break; }
   }

   // Check if function is found
   if (!TDFunction)
   {
      cout <<"Functions:";
      for (int i = 0; i < funNames.size(); i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" can not be found in "<<libName<<endl;

      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout<<"Function will be set to homogenous."<<endl;
   }
}

// Evaluate coefficient
void LibVectorCoefficient::Eval(Vector &V,
                                ElementTransformation &T,
                                const IntegrationPoint &ip)
{
   V.SetSize(vdim);
   // Homegenous if not defined
   if (!TDFunction)
   {
      V = 0.0;
      return;
   }

   // Evaluate library function
   T.Transform(ip, x);
   (*TDFunction)(x.GetData(), x.Size(), GetTime(), V.GetData(), V.Size());
}

//  Destructor
LibVectorCoefficient::~LibVectorCoefficient()
{
   if (libHandle) { dlclose(libHandle); }
}
