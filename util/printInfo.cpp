// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#if __has_include("buildInfo.hpp")
#include "buildInfo.hpp"
#else
#include <string>
#include <fstream>
#include <sstream>
std::istringstream buildInfo(R"~(
------------------------------------
No build information
------------------------------------
)~");
#endif

#if defined(_WIN32)
   #include <winsock.h>
#else
   #include <unistd.h>
#endif

#include "mfem.hpp"

using namespace std;
using namespace mfem;

void line(int len)
{
   for (int b=0; b<len; ++b) { cout<<"-"; }
   cout<<"\n";
}

// Function for printing compile time and runtime information.
void printInfo()
{
   if (Mpi::Root())
   {
      // Print header
      line(80);
      cout<<R"(    _____  ______      ____  __  _____  )"<<endl;
      cout<<R"(   |  __ \|  _ \ \    / /  \/  |/ ____| )"<<endl;
      cout<<R"(   | |__) | |_) \ \  / /| \  / | (___   )"<<endl;
      cout<<R"(   |  _  /|  _ < \ \/ / | |\/| |\___ \  )"<<endl;
      cout<<R"(   | | \ \| |_) | \  /  | |  | |____) | )"<<endl;
      cout<<R"(   |_|  \_\____/   \/   |_|  |_|_____/  )"<<endl;
      cout<<R"(                                        )"<<endl;

      // Build info
      line(80);
      cout<<"Compile time info\n";
      line(80);
      cout<< buildInfo.str() << endl;

      // Run info
      line(80);
      cout<<"Run time info"<<endl;
      line(80); cout<<endl;
      time_t     now = time(0);
      struct tm  tstruct = *localtime(&now);
      char       time[80], host[80];
      strftime(time, sizeof(time), "%Y-%m-%d.%X", &tstruct);
      gethostname(host,sizeof(host));

      cout<<"Time: "<<time<<endl;
      cout<<"Numer of MPI ranks "<<Mpi::WorldSize()<<endl;

      cout<<"List  of hosts\n0: "<<host<<endl;

      // Receive hostnames from non-root nodes
      for (int i = 1; i < Mpi::WorldSize(); i++)
      {
         MPI_Status status;
         MPI_Recv (&host, sizeof(host), MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
         cout<<i<<": "<<host<<endl;
      }

      cout<<endl; line(80); cout<<endl;
   }
   else
   {
      // Send hostnames from non-root nodes to root
      char host[80];
      gethostname(host,sizeof(host));
      MPI_Send (&host, sizeof(host), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
   }
}

/// This class help monitor the convergence of the linear Krylov solve.
class GeneralResidualMonitor : public IterativeSolverMonitor
{
private:
   const std::string prefix;
   int interval;
   mutable real_t norm0;

public:
   /// Constructor
   GeneralResidualMonitor(const std::string& prefix_,
                          int print_iv)
      : prefix(prefix_), interval(print_iv)
   {
      if (Mpi::Root()) { interval = -1; }
   }

   /// Print residual
   virtual void MonitorResidual(int it,
                                real_t norm,
                                const Vector &r,
                                bool final)
   {
      if (interval < 0) { return; }
      if (it == 0) { norm0 = norm; }

      if ( ( it%interval == 0) || final )
      {
         mfem::out<<prefix<<" iteration "<<std::setw(3)<<it
                  <<std::setw(8)<<std::defaultfloat<<std::setprecision(3)
                  <<": ||r|| = "<<norm
                  <<std::setw(6)<<std::fixed<<std::setprecision(2)
                  <<", ||r||/||r_0|| = "<<100*norm/norm0<<" %\n";
      }
   }
};
