/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T>
void Gather
( const AbstractDistMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::Gather"))
    AssertSameGrids( A, B );

    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    if( A.CrossSize() != 1 )
        LogicError("Non-trivial cross teams not yet supported");
    if( !A.Grid().InGrid() )
        return;

    // Translate the root into our DistComm (if possible)
    const Int root = B.Root();
    const Int target = mpi::Translate( B.CrossComm(), root, A.DistComm() );
    if( target == mpi::UNDEFINED )
        return;

    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    const Int mLocalMax = MaxLength(m,colStride);
    const Int nLocalMax = MaxLength(n,rowStride);
    const Int pkgSize = mpi::Pad( mLocalMax*nLocalMax );
    const Int numDist = A.DistSize();

    std::vector<T> buffer;
    T *sendBuf, *recvBuf;
    if( B.CrossRank() == root )
    {
        buffer.resize( (numDist+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        buffer.resize( pkgSize );
        sendBuf = &buffer[0];
        recvBuf = 0;
    }

    // Pack
    copy::util::InterleaveMatrix
    ( A.LocalHeight(), A.LocalWidth(),
      A.LockedBuffer(), 1, A.LDim(),
      sendBuf,          1, A.LocalHeight() );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, target, A.DistComm() );

    // Unpack
    if( B.CrossRank() == root )
        copy::util::StridedUnpack
        ( m, n,
          A.ColAlign(), colStride,
          A.RowAlign(), rowStride,
          recvBuf, pkgSize,
          B.Buffer(), B.LDim() );
}

#define PROTO(T) \
  template void Gather \
  ( const AbstractDistMatrix<T>& A, \
          DistMatrix<T,CIRC,CIRC>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
