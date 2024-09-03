// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#if __has_include("buildInfo.hpp")
#include "buildInfo.hpp"
#else
#include "noInfo.hpp"
#endif

#if defined(_WIN32)
#include <winsock.h>
#else
#include <unistd.h>
#endif

#include "mfem.hpp"

using namespace std;
using namespace mfem;

// Function for printing compile time and runtime information.
void printInfo()
{
   if (Mpi::Root())
   {
      // Print header
      cout<<"----------------------------------------\n";
      cout<<R"(    _____  ______      ____  __  _____  )"<<endl;
      cout<<R"(   |  __ \|  _ \ \    / /  \/  |/ ____| )"<<endl;
      cout<<R"(   | |__) | |_) \ \  / /| \  / | (___   )"<<endl;
      cout<<R"(   |  _  /|  _ < \ \/ / | |\/| |\___ \  )"<<endl;
      cout<<R"(   | | \ \| |_) | \  /  | |  | |____) | )"<<endl;
      cout<<R"(   |_|  \_\____/   \/   |_|  |_|_____/  )"<<endl;
      cout<<R"(                                        )"<<endl;

      // Build info
      cout<<"----------------------------------------\n";
      cout<<"Compile time info\n";
      cout<<"----------------------------------------\n";
      cout<< buildInfo.str() << endl;

      // Run info
      cout<<"----------------------------------------\n";
      cout<<"Run time info"<<endl;
      cout<<"----------------------------------------\n\n";
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

      cout<<"\n----------------------------------------\n\n";
   }
   else
   {
      // Send hostnames from non-root nodes to root
      char host[80];
      gethostname(host,sizeof(host));
      MPI_Send (&host, sizeof(host), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
   }
}
