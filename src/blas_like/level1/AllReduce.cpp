/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void AllReduce( AbstractDistMatrix<T>& A, mpi::Comm comm, mpi::Op op )
{
    DEBUG_ONLY(CallStackEntry cse("AllReduce"))
    if( !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int localSize = localHeight*localWidth;
    if( localHeight == A.LDim() )
    {
        mpi::AllReduce( A.Buffer(), localSize, op, comm );
    }
    else
    {
        std::vector<T> buf( localSize );

        // Pack
        copy::util::InterleaveMatrix
        ( localHeight, localWidth,
          A.LockedBuffer(), 1, A.LDim(),
          buf.data(),       1, localHeight );

        mpi::AllReduce( buf.data(), localSize, op, comm );

        // Unpack
        copy::util::InterleaveMatrix
        ( localHeight, localWidth,
          buf.data(), 1, localHeight,
          A.Buffer(), 1, A.LDim() );
    }
}

template<typename T>
void AllReduce( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, mpi::Op op )
{
    DEBUG_ONLY(CallStackEntry cse("AllReduce"))
    if( !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int localSize = localHeight*localWidth;
    if( localHeight == A.LDim() )
    {
        mpi::AllReduce( A.Buffer(), localSize, op, comm );
    }
    else
    {
        std::vector<T> buf( localSize );

        // Pack
        copy::util::InterleaveMatrix
        ( localHeight, localWidth,
          A.LockedBuffer(), 1, A.LDim(),
          buf.data(),       1, localHeight );

        mpi::AllReduce( buf.data(), localSize, op, comm );

        // Unpack
        copy::util::InterleaveMatrix
        ( localHeight, localWidth,
          buf.data(), 1, localHeight,
          A.Buffer(), 1, A.LDim() );
    }
}

#define PROTO(T) \
  template void AllReduce \
  ( AbstractDistMatrix<T>& A, mpi::Comm comm, mpi::Op op ); \
  template void AllReduce \
  ( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, mpi::Op op );

#include "El/macros/Instantiate.h"

} // namespace El
