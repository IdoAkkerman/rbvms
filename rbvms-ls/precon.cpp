// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "precon.hpp"

using namespace mfem;
using namespace RBVMS;

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{
   BlockOperator *jacobian = (BlockOperator *) &op;

   for (int i = 0; i < prec.Size(); ++i)
   {
      if (prec[i])
      {
         prec[i]->SetOperator(jacobian->GetBlock(i,i));
         SetDiagonalBlock(i, prec[i]);
      }
      for (int j = i+1; j < prec.Size(); ++j)
      {
         SetBlock(j,i, &jacobian->GetBlock(j,i));
      }
   }
}

// Destructor
JacobianPreconditioner::~JacobianPreconditioner()
{
   for (int i = 0; i < prec.Size(); ++i)
   {
      if (prec[i]) { delete prec[i]; }
   }
}

