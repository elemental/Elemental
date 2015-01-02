/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void Filter
( const DistMatrix<T,Collect<U>(),Collect<V>()>& A,
        DistMatrix<T,        U,           V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::Filter"))
    AssertSameGrids( A, B );

    B.Resize( A.Height(), A.Width() );
    if( !B.Participating() )
        return;

    const Int colShift = B.ColShift();
    const Int rowShift = B.RowShift();
    util::InterleaveMatrix
    ( B.LocalHeight(), B.LocalWidth(),
      A.LockedBuffer(colShift,rowShift), B.ColStride(), B.RowStride()*A.LDim(),
      B.Buffer(),                        1,             B.LDim() );
}

template<typename T,Dist U,Dist V>
void Filter
( const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A,
        BlockDistMatrix<T,        U,           V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::Filter"))
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void Filter \
  ( const DistMatrix<T,Collect<U>(),Collect<V>()>& A, \
          DistMatrix<T,        U,           V   >& B ); \
  template void Filter \
  ( const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A, \
          BlockDistMatrix<T,        U,           V   >& B );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR) 

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
