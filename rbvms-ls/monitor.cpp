// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "monitor.hpp"

using namespace mfem;
using namespace RBVMS;

// Constructor
GeneralResidualMonitor::GeneralResidualMonitor(MPI_Comm comm,
                                               const std::string& prefix_,
                                               int print_lvl)
   : prefix(prefix_)
{
   if (Mpi::Root())
   {
      print_level = print_lvl;
   }
   else
   {
      print_level = -1;
   }
}

// Print residual
void GeneralResidualMonitor::MonitorResidual(int it,
                                             real_t norm,
                                             const Vector &r,
                                             bool final)
{
   if (it == 0)
   {
      norm0 = norm;
   }

   if ((print_level > 0) && (( it%print_level == 0) || final))
   {
      mfem::out<<prefix<<" iteration "<<std::setw(3)<<it
               <<": ||r|| = "
               <<std::setw(8)<<std::defaultfloat<<std::setprecision(3)
               <<norm
               <<", ||r||/||r_0|| = "
               <<std::setw(6)<<std::fixed<<std::setprecision(2)
               <<100*norm/norm0<<" %\n";
   }
}

// Constructor
SystemResidualMonitor::SystemResidualMonitor(MPI_Comm comm,
                                             const std::string& prefix_,
                                             int print_lvl,
                                             Array<int> &offsets)
   : prefix(prefix_), bOffsets(offsets)
{
   if (Mpi::Root())
   {
      print_level = print_lvl;
   }
   else
   {
      print_level = -1;
   }
   nvar = bOffsets.Size()-1;
   norm0.SetSize(nvar);
}

// Print residual
void SystemResidualMonitor::MonitorResidual(int it,
                                            real_t norm,
                                            const Vector &r,
                                            bool final)
{
   Vector vnorm(nvar);

   for (int i = 0; i < nvar; ++i)
   {
      Vector r_i(r.GetData() + bOffsets[i], bOffsets[i+1] - bOffsets[i]);
      vnorm[i] = sqrt(InnerProduct(MPI_COMM_WORLD, r_i, r_i));
      if (it == 0) { norm0[i] = vnorm[i]; }
   }

   bool print = (print_level > 0) &&  ((it%print_level == 0) || final);
   if (print)
   {
      mfem::out << prefix << " iteration " << std::setw(3) << it <<"\n"
                << " ||r||  \t"<< "||r||/||r_0||\n";
      for (int i = 0; i < nvar; ++i)
      {
         mfem::out<<std::setw(8)<<std::defaultfloat<<std::setprecision(4)
                  <<vnorm[i]<<"\t"
                  <<std::setw(8)<<std::fixed<<std::setprecision(2)
                  <<100*vnorm[i]/norm0[i]<<" %\n";
      }
   }
   mfem::out<<std::flush;
}
