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
void Broadcast( AbstractDistMatrix<T>& A, mpi::Comm comm, Int rank )
{
    DEBUG_ONLY(CallStackEntry cse("Broadcast"))
    if( !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int localSize = localHeight*localWidth;
    if( localHeight == A.LDim() )
    {
        mpi::Broadcast( A.Buffer(), localSize, rank, comm );
    }
    else
    {
        std::vector<T> buf( localSize );

        // Pack
        if( mpi::Rank(comm) == rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              A.LockedBuffer(), 1, A.LDim(),
              buf.data(),       1, localHeight );

        mpi::Broadcast( buf.data(), localSize, rank, comm );

        // Unpack
        if( mpi::Rank(comm) != rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              buf.data(), 1, localHeight,
              A.Buffer(), 1, A.LDim() );
    }
}

template<typename T>
void Broadcast( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, Int rank )
{
    DEBUG_ONLY(CallStackEntry cse("Broadcast"))
    if( !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int localSize = localHeight*localWidth;
    if( localHeight == A.LDim() )
    {
        mpi::Broadcast( A.Buffer(), localSize, rank, comm );
    }
    else
    {
        std::vector<T> buf( localSize );

        // Pack
        if( mpi::Rank(comm) == rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              A.LockedBuffer(), 1, A.LDim(),
              buf.data(),       1, localHeight );

        mpi::Broadcast( buf.data(), localSize, rank, comm );

        // Unpack
        if( mpi::Rank(comm) != rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              buf.data(), 1, localHeight,
              A.Buffer(), 1, A.LDim() );
    }
}


#define PROTO(T) \
  template void Broadcast \
  ( AbstractDistMatrix<T>& A, mpi::Comm comm, Int rank ); \
  template void Broadcast \
  ( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, Int rank );

#include "El/macros/Instantiate.h"

} // namespace El
